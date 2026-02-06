#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import logging
import re
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import GroupKFold

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from ml.calibration import apply_calibration, fit_oof_calibrator
from ml.config import FeaturesConfig, config_to_dict, load_atlas_cfg, load_config
from ml.data.bigbind import load_train_holdout_from_bigbind
from ml.evaluate import DEFAULT_FRACTIONS, evaluate_holdout_metrics
from ml.featurize import BigBindFeaturizer
from ml.fpocket_bigbind import merge_fpocket_metrics_on_pocket, precompute_fpocket_for_bigbind_df
from ml.hard_negatives import merge_hard_negatives
from ml.labels import apply_label_smoothing, derive_sample_weight
from ml.models import build_model
from ml.registry import (
    append_registry_index,
    collect_env_versions,
    compute_dataset_hash,
    get_git_sha,
    write_registry_record,
)


_FEATURES_ALL_VARIANTS_NAME = "features_all_variants.csv"
_ONBITS_ALL_VARIANTS_NAME = "morgan_onbits_all_variants.csv"
_FEATURE_AUDIT_META_COLUMNS = (
    "row_number",
    "source_index",
    "active",
    "lig_smiles",
    "ex_rec_pdb",
    "pocket",
    "murcko_scaffold",
)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run ML feature ablation experiments with outer holdout-by-train_pdb "
            "and inner scaffold-grouped validation."
        )
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to ml config file (JSON or YAML).",
    )
    return parser


def _setup_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    return logging.getLogger("ml.experiment_runner")


def _resolve_run_id(config_run_id: str | None) -> str:
    if config_run_id:
        return config_run_id
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def _normalize_groups(murcko_scaffolds: list[str]) -> np.ndarray:
    groups: list[str] = []
    for idx, scaffold in enumerate(murcko_scaffolds):
        text = str(scaffold or "").strip()
        groups.append(text if text else f"__empty_scaffold_{idx}")
    return np.asarray(groups, dtype=object)


def _inner_grouped_metrics(
    *,
    X: sparse.csr_matrix,
    y: np.ndarray,
    murcko_scaffolds: list[str],
    n_splits_requested: int,
    model_family: str,
    model_params: dict[str, Any],
    calibration_enabled: bool,
    calibration_method: str,
    calibration_cv_folds: int,
    calibration_seed: int,
    random_seed: int,
    ranking_query_groups: np.ndarray | None = None,
    sample_weight: np.ndarray | None = None,
) -> dict[str, float | int]:
    groups = _normalize_groups(murcko_scaffolds)
    unique_groups = np.unique(groups)
    n_splits = min(max(2, int(n_splits_requested)), int(unique_groups.shape[0]), int(y.shape[0]))
    if n_splits < 2:
        return {"inner_folds": 0}

    splitter = GroupKFold(n_splits=n_splits)
    fold_reports: list[dict[str, float | int]] = []
    for fold_idx, (train_idx, val_idx) in enumerate(splitter.split(X, y, groups)):
        y_train_fold = y[train_idx]
        if np.unique(y_train_fold).size < 2:
            continue
        y_val_fold = y[val_idx]
        if np.unique(y_val_fold).size < 2:
            continue

        model = build_model(
            model_family=model_family,
            model_params=model_params,
            random_seed=random_seed + fold_idx,
        )
        fit_kwargs: dict[str, Any] = {}
        if ranking_query_groups is not None:
            fit_kwargs["query_groups"] = ranking_query_groups[train_idx]
        if sample_weight is not None:
            fit_kwargs["sample_weight"] = sample_weight[train_idx]
        model.fit(X[train_idx], y_train_fold, **fit_kwargs)
        val_prob_uncal = model.predict_proba(X[val_idx])[:, 1]
        val_prob_active = np.asarray(val_prob_uncal, dtype=float)
        if calibration_enabled:
            train_groups = [murcko_scaffolds[idx] for idx in train_idx.tolist()]

            def _builder(seed_value: int):
                return build_model(
                    model_family=model_family,
                    model_params=model_params,
                    random_seed=seed_value,
                )

            calibrator, _ = fit_oof_calibrator(
                X=X[train_idx],
                y=y_train_fold,
                groups=np.asarray(train_groups, dtype=object),
                base_model_builder=_builder,
                method=calibration_method,
                cv_folds=calibration_cv_folds,
                seed=calibration_seed + fold_idx,
                fit_query_groups=(
                    ranking_query_groups[train_idx] if ranking_query_groups is not None else None
                ),
                fit_sample_weights=sample_weight[train_idx] if sample_weight is not None else None,
            )
            val_prob_active = apply_calibration(calibrator, val_prob_uncal)
        fold_reports.append(
            evaluate_holdout_metrics(
                y_val_fold,
                val_prob_active,
                fractions=DEFAULT_FRACTIONS,
            )
        )

    if not fold_reports:
        return {"inner_folds": 0}

    aggregated: dict[str, float | int] = {"inner_folds": len(fold_reports)}
    metric_keys = sorted(
        {
            key
            for report in fold_reports
            for key, value in report.items()
            if isinstance(value, (int, float))
        }
    )
    for key in metric_keys:
        values = np.asarray(
            [float(report.get(key, float("nan"))) for report in fold_reports], dtype=float
        )
        aggregated[f"inner_mean_{key}"] = float(np.nanmean(values))
        aggregated[f"inner_std_{key}"] = float(np.nanstd(values))

    return aggregated


def _build_feature_variants(base_features: FeaturesConfig) -> list[tuple[str, FeaturesConfig]]:
    fp_bits = max(0, int(base_features.ligand_morgan_fp_bits))
    vina_flag = bool(base_features.vina_score)
    pocket_feature_flag = bool(base_features.pocket_features)
    baseline = FeaturesConfig(
        ligand_descriptors=True,
        ligand_extra_descriptors=False,
        ligand_morgan_fp_bits=fp_bits,
        pocket_fpocket=True,
        pocket_features=pocket_feature_flag,
        vina_score=vina_flag,
    )
    return [
        ("baseline", baseline),
        ("ligand_only", replace(baseline, pocket_fpocket=False)),
        (
            "pocket_only",
            FeaturesConfig(
                ligand_descriptors=False,
                ligand_extra_descriptors=False,
                ligand_morgan_fp_bits=0,
                pocket_fpocket=True,
                pocket_features=pocket_feature_flag,
                vina_score=False,
            ),
        ),
        ("no_fp", replace(baseline, ligand_morgan_fp_bits=0)),
        ("add_ligand_extras", replace(baseline, ligand_extra_descriptors=True)),
    ]


def _build_consolidated_feature_columns(
    *,
    feature_variants: list[tuple[str, FeaturesConfig]],
    bigbind_root: Path,
    atlas_cfg: dict[str, Any],
    config,
    logger: logging.Logger,
) -> list[str]:
    columns = list(_FEATURE_AUDIT_META_COLUMNS)
    for _variant_name, features in feature_variants:
        featurizer = BigBindFeaturizer(
            bigbind_root=bigbind_root,
            features=features,
            atlas_cfg=atlas_cfg,
            fpocket_center_columns=config.fpocket_center_columns,
            fpocket_variant_column=config.fpocket_variant_column,
            fpocket_ph_column=config.fpocket_ph_column,
            fpocket_centers_by_pdb=config.fpocket_centers_by_pdb,
            logger=logger,
        )
        metadata = featurizer.metadata()
        continuous_names = metadata.get("continuous_feature_names")
        if isinstance(continuous_names, list):
            for name in continuous_names:
                text = str(name)
                if text not in columns:
                    columns.append(text)
    return columns


def _append_feature_rows(
    *,
    source_csv: Path,
    output_csv: Path,
    variant_name: str,
    split_name: str,
    feature_columns: list[str],
) -> None:
    if not source_csv.exists():
        return

    write_header = not output_csv.exists()
    with source_csv.open("r", encoding="utf-8", newline="") as src_handle:
        reader = csv.DictReader(src_handle)
        if not reader.fieldnames:
            return

        with output_csv.open("a", encoding="utf-8", newline="") as dst_handle:
            writer = csv.writer(dst_handle)
            if write_header:
                writer.writerow(["variant", "split", *feature_columns])
            for row in reader:
                writer.writerow([variant_name, split_name, *[row.get(col, "") for col in feature_columns]])


def _append_onbits_rows(
    *,
    source_csv: Path,
    output_csv: Path,
    variant_name: str,
    split_name: str,
) -> None:
    if not source_csv.exists():
        return

    write_header = not output_csv.exists()
    with source_csv.open("r", encoding="utf-8", newline="") as src_handle:
        reader = csv.DictReader(src_handle)
        if not reader.fieldnames:
            return

        with output_csv.open("a", encoding="utf-8", newline="") as dst_handle:
            writer = csv.writer(dst_handle)
            if write_header:
                writer.writerow(["variant", "split", "row_number", "bit"])
            for row in reader:
                writer.writerow(
                    [
                        variant_name,
                        split_name,
                        row.get("row_number", ""),
                        row.get("bit", ""),
                    ]
                )


def _append_variant_audits_to_consolidated(
    *,
    output_dir: Path,
    variant_name: str,
    safe_variant: str,
    feature_columns: list[str],
) -> None:
    for split_name in ("train", "holdout"):
        features_src = output_dir / f"features_{split_name}_{safe_variant}.csv"
        onbits_src = output_dir / f"morgan_onbits_{split_name}_{safe_variant}.csv"
        _append_feature_rows(
            source_csv=features_src,
            output_csv=output_dir / _FEATURES_ALL_VARIANTS_NAME,
            variant_name=variant_name,
            split_name=split_name,
            feature_columns=feature_columns,
        )
        _append_onbits_rows(
            source_csv=onbits_src,
            output_csv=output_dir / _ONBITS_ALL_VARIANTS_NAME,
            variant_name=variant_name,
            split_name=split_name,
        )


def _collect_fpocket_cache_keys(metrics_df: pd.DataFrame | None) -> list[str]:
    if metrics_df is None or metrics_df.empty:
        return []
    keys: list[str] = []
    for _, row in metrics_df.iterrows():
        pocket = str(row.get("pocket", "")).strip()
        receptor = str(row.get("receptor_pdb", "")).strip()
        cx = row.get("pocket_center_x")
        cy = row.get("pocket_center_y")
        cz = row.get("pocket_center_z")
        keys.append(f"{pocket}|{receptor}|{cx}|{cy}|{cz}")
    return sorted(set(keys))


def _run_variant(
    *,
    variant_name: str,
    features: FeaturesConfig,
    bigbind_root: Path,
    train_df: pd.DataFrame,
    holdout_df: pd.DataFrame,
    output_dir: Path,
    consolidated_feature_columns: list[str],
    atlas_cfg: dict[str, Any],
    config,
    logger: logging.Logger,
) -> dict[str, float | int | str]:
    safe_variant = re.sub(r"[^A-Za-z0-9_-]+", "_", variant_name.strip()) or "variant"
    featurizer = BigBindFeaturizer(
        bigbind_root=bigbind_root,
        features=features,
        atlas_cfg=atlas_cfg,
        fpocket_center_columns=config.fpocket_center_columns,
        fpocket_variant_column=config.fpocket_variant_column,
        fpocket_ph_column=config.fpocket_ph_column,
        fpocket_centers_by_pdb=config.fpocket_centers_by_pdb,
        logger=logger,
    )
    train_features = featurizer.transform(
        train_df,
        audit_dir=output_dir,
        audit_tag=f"train_{safe_variant}",
    )
    holdout_features = featurizer.transform(
        holdout_df,
        audit_dir=output_dir,
        audit_tag=f"holdout_{safe_variant}",
    )
    _append_variant_audits_to_consolidated(
        output_dir=output_dir,
        variant_name=variant_name,
        safe_variant=safe_variant,
        feature_columns=consolidated_feature_columns,
    )

    model_family = str(getattr(config, "model_family", "logreg")).strip().lower()
    if model_family == "lgbm":
        model_family = "lightgbm"
    if model_family == "xgb":
        model_family = "xgboost"
    rank_enabled = bool(getattr(config.rank, "enabled", False))
    if rank_enabled and "rank" not in model_family and model_family in {"lightgbm", "xgboost"}:
        model_family = f"{model_family}_rank"
    model_params = dict(getattr(config, "model_params", {}) or {})
    train_used = train_df.loc[train_features.source_index].copy()

    def _group_values(df: pd.DataFrame, key: str) -> np.ndarray:
        key_norm = str(key).strip().lower()
        if key_norm == "target":
            return df["ex_rec_pdb"].astype(str).to_numpy(dtype=object)
        if key_norm == "pocket":
            return df["pocket"].astype(str).to_numpy(dtype=object)
        return (df["ex_rec_pdb"].astype(str) + "::" + df["pocket"].astype(str)).to_numpy(dtype=object)

    train_rank_groups = (
        _group_values(train_used, str(config.rank.group_key))
        if "rank" in model_family
        else None
    )
    sample_weight = (
        derive_sample_weight(train_used, weight_cap=config.labels.weight_cap)
        if config.labels.use_sample_weights
        else None
    )
    y_train_fit = train_features.y.astype(float if config.labels.smoothing_eps > 0 else int)
    if config.labels.smoothing_eps > 0 and "logreg" not in model_family and "rank" not in model_family:
        y_train_fit = apply_label_smoothing(train_features.y, config.labels.smoothing_eps)

    inner_metrics = _inner_grouped_metrics(
        X=train_features.X,
        y=train_features.y,
        murcko_scaffolds=train_features.murcko_scaffolds,
        n_splits_requested=config.inner_scaffold_folds,
        model_family=model_family,
        model_params=model_params,
        calibration_enabled=bool(config.calibration_enabled),
        calibration_method=str(config.calibration_method),
        calibration_cv_folds=int(config.calibration_cv_folds),
        calibration_seed=int(config.calibration_seed),
        random_seed=config.random_seed,
        ranking_query_groups=train_rank_groups,
        sample_weight=sample_weight,
    )

    model = build_model(
        model_family=model_family,
        model_params=model_params,
        random_seed=config.random_seed,
    )
    model.fit(
        train_features.X,
        np.asarray(y_train_fit),
        query_groups=train_rank_groups,
        sample_weight=sample_weight,
    )

    holdout_prob_uncal = np.asarray(model.predict_proba(holdout_features.X)[:, 1], dtype=float)
    holdout_prob_active = holdout_prob_uncal
    if config.calibration_enabled:

        def _builder(seed_value: int):
            return build_model(
                model_family=model_family,
                model_params=model_params,
                random_seed=seed_value,
            )

        calibrator, _ = fit_oof_calibrator(
            X=train_features.X,
            y=train_features.y,
            groups=np.asarray(train_features.murcko_scaffolds, dtype=object),
            base_model_builder=_builder,
            method=config.calibration_method,
            cv_folds=config.calibration_cv_folds,
            seed=config.calibration_seed,
            fit_query_groups=train_rank_groups,
            fit_sample_weights=sample_weight,
        )
        holdout_prob_active = apply_calibration(calibrator, holdout_prob_uncal)

    holdout_metrics = evaluate_holdout_metrics(
        holdout_features.y,
        holdout_prob_active,
        fractions=DEFAULT_FRACTIONS,
    )

    row: dict[str, float | int | str] = {
        "variant": variant_name,
        "train_rows": int(train_features.X.shape[0]),
        "holdout_rows": int(holdout_features.X.shape[0]),
        "feature_count": int(train_features.X.shape[1]),
        "ligand_descriptors": int(features.ligand_descriptors),
        "ligand_extra_descriptors": int(features.ligand_extra_descriptors),
        "ligand_morgan_fp_bits": int(features.ligand_morgan_fp_bits),
        "pocket_fpocket": int(features.pocket_fpocket),
        "pocket_features": int(features.pocket_features),
        "vina_score": int(features.vina_score),
        "model_family": model_family,
        "model_params_json": json.dumps(model_params, sort_keys=True),
        "calibration_enabled": int(bool(config.calibration_enabled)),
        "calibration_method": str(config.calibration_method),
        "rank_enabled": int(rank_enabled),
        "rank_group_key": str(config.rank.group_key),
        "fpocket_loaded_count": int(featurizer.loaded_fpocket_count),
        "fpocket_missing_center_count": int(featurizer.missing_center_count),
        "fpocket_missing_info_count": int(featurizer.missing_info_count),
    }
    row.update(inner_metrics)
    row.update({f"holdout_{key}": value for key, value in holdout_metrics.items()})
    return row


def run_feature_ablation_experiments(
    config_path: str,
) -> tuple[Path, pd.DataFrame]:
    logger = _setup_logging()
    config = load_config(config_path)
    atlas_cfg = load_atlas_cfg(config)
    np.random.seed(config.random_seed)

    dataset = load_train_holdout_from_bigbind(
        bigbind_dir=config.bigbind_dir,
        train_pdb=config.train_pdb,
        splits=config.splits,
        max_rows=config.max_rows,
        random_seed=config.random_seed,
        dataset_split_mode=config.dataset_split_mode,
        exclude_target_prefixes=config.exclude_target_prefixes,
        logger=logger,
    )
    feature_variants = _build_feature_variants(config.features)
    run_id = _resolve_run_id(config.run_id)
    output_dir = Path(__file__).resolve().parent / "outputs" / run_id
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / _FEATURES_ALL_VARIANTS_NAME).unlink(missing_ok=True)
    (output_dir / _ONBITS_ALL_VARIANTS_NAME).unlink(missing_ok=True)
    consolidated_feature_columns = _build_consolidated_feature_columns(
        feature_variants=feature_variants,
        bigbind_root=dataset.bigbind_root,
        atlas_cfg=atlas_cfg,
        config=config,
        logger=logger,
    )

    train_df = dataset.train_df.copy()
    holdout_df = dataset.holdout_df.copy()
    if config.hard_negatives.enabled:
        group_key = str(config.rank.group_key if config.rank.enabled else "target")
        train_df = merge_hard_negatives(
            train_df,
            policy=config.hard_negatives.policy,
            ratio=config.hard_negatives.ratio,
            per_active_k=config.hard_negatives.per_active_k,
            within_group=config.hard_negatives.within_group,
            group_key=group_key if config.hard_negatives.within_group else "target",
            max_candidates=config.hard_negatives.max_candidates,
            seed=config.hard_negatives.seed,
        )
    fpocket_metrics_df: pd.DataFrame | None = None
    if any(
        (features.pocket_fpocket or features.pocket_features)
        for _, features in feature_variants
    ):
        combined = pd.concat([train_df, holdout_df], ignore_index=False)
        fpocket_metrics_df = precompute_fpocket_for_bigbind_df(
            df=combined,
            bigbind_root=dataset.bigbind_root,
            atlas_cfg=atlas_cfg,
            run_dir=output_dir,
        )
        train_df = merge_fpocket_metrics_on_pocket(train_df, fpocket_metrics_df)
        holdout_df = merge_fpocket_metrics_on_pocket(holdout_df, fpocket_metrics_df)

    rows: list[dict[str, float | int | str]] = []
    for variant_name, features in feature_variants:
        logger.info("[ml.experiment] running variant=%s", variant_name)
        rows.append(
            _run_variant(
                variant_name=variant_name,
                features=features,
                bigbind_root=dataset.bigbind_root,
                train_df=train_df,
                holdout_df=holdout_df,
                output_dir=output_dir,
                consolidated_feature_columns=consolidated_feature_columns,
                atlas_cfg=atlas_cfg,
                config=config,
                logger=logger,
            )
        )

    results_df = pd.DataFrame(rows)

    results_path = output_dir / "experiment_results.csv"
    results_df.to_csv(results_path, index=False)

    (output_dir / "config_snapshot.txt").write_text(
        json.dumps(config_to_dict(config), indent=2),
        encoding="utf-8",
    )
    (output_dir / "config_snapshot.json").write_text(
        json.dumps(config_to_dict(config), indent=2),
        encoding="utf-8",
    )
    if config.registry_enabled:
        repo_root = Path(__file__).resolve().parents[1]
        dataset_hash = compute_dataset_hash(
            train_df=train_df,
            holdout_df=holdout_df,
            feature_columns=consolidated_feature_columns,
            config_snapshot=config_to_dict(config),
        )
        registry_record = {
            "run_id": run_id,
            "run_dir": str(output_dir),
            "dataset_hash": dataset_hash,
            "git_sha": get_git_sha(repo_root),
            "num_variants": int(len(results_df)),
            "fpocket_cache_keys_used": _collect_fpocket_cache_keys(fpocket_metrics_df),
            "environment": collect_env_versions(),
        }
        write_registry_record(output_dir, registry_record)
        append_registry_index(
            output_dir.parent / "registry_index.csv",
            {
                "run_id": run_id,
                "run_dir": str(output_dir),
                "dataset_hash": dataset_hash,
                "git_sha": registry_record["git_sha"],
                "num_variants": int(len(results_df)),
            },
        )

    return output_dir, results_df


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    output_dir, results = run_feature_ablation_experiments(args.config)
    print(
        "[ml.experiment] "
        f"saved={output_dir / 'experiment_results.csv'} "
        f"variants={len(results)}"
    )


if __name__ == "__main__":
    main()
