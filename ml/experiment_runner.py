#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
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

from ml.config import FeaturesConfig, config_to_dict, load_atlas_cfg, load_config
from ml.data.bigbind import load_train_holdout_from_bigbind
from ml.evaluate import DEFAULT_FRACTIONS, evaluate_holdout_metrics
from ml.featurize import BigBindFeaturizer
from ml.modeling import train_logistic_regression


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
    model_config: Any,
    random_seed: int,
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

        model = train_logistic_regression(
            X_train=X[train_idx],
            y_train=y_train_fold,
            model_config=model_config,
            random_seed=random_seed + fold_idx,
        )
        val_prob_active = model.predict_proba(X[val_idx])[:, 1]
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
    baseline = FeaturesConfig(
        ligand_descriptors=True,
        ligand_extra_descriptors=False,
        ligand_morgan_fp_bits=fp_bits,
        pocket_fpocket=True,
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
                vina_score=False,
            ),
        ),
        ("no_fp", replace(baseline, ligand_morgan_fp_bits=0)),
        ("add_ligand_extras", replace(baseline, ligand_extra_descriptors=True)),
    ]


def _run_variant(
    *,
    variant_name: str,
    features: FeaturesConfig,
    dataset,
    atlas_cfg: dict[str, Any],
    config,
    logger: logging.Logger,
) -> dict[str, float | int | str]:
    featurizer = BigBindFeaturizer(
        bigbind_root=dataset.bigbind_root,
        features=features,
        atlas_cfg=atlas_cfg,
        fpocket_center_columns=config.fpocket_center_columns,
        fpocket_variant_column=config.fpocket_variant_column,
        fpocket_ph_column=config.fpocket_ph_column,
        fpocket_centers_by_pdb=config.fpocket_centers_by_pdb,
        logger=logger,
    )
    train_features = featurizer.transform(dataset.train_df)
    holdout_features = featurizer.transform(dataset.holdout_df)

    inner_metrics = _inner_grouped_metrics(
        X=train_features.X,
        y=train_features.y,
        murcko_scaffolds=train_features.murcko_scaffolds,
        n_splits_requested=config.inner_scaffold_folds,
        model_config=config.model,
        random_seed=config.random_seed,
    )

    model = train_logistic_regression(
        X_train=train_features.X,
        y_train=train_features.y,
        model_config=config.model,
        random_seed=config.random_seed,
    )
    holdout_prob_active = model.predict_proba(holdout_features.X)[:, 1]
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
        "vina_score": int(features.vina_score),
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
        logger=logger,
    )
    feature_variants = _build_feature_variants(config.features)

    rows: list[dict[str, float | int | str]] = []
    for variant_name, features in feature_variants:
        logger.info("[ml.experiment] running variant=%s", variant_name)
        rows.append(
            _run_variant(
                variant_name=variant_name,
                features=features,
                dataset=dataset,
                atlas_cfg=atlas_cfg,
                config=config,
                logger=logger,
            )
        )

    results_df = pd.DataFrame(rows)

    run_id = _resolve_run_id(config.run_id)
    output_dir = Path(__file__).resolve().parent / "outputs" / run_id
    output_dir.mkdir(parents=True, exist_ok=True)
    results_path = output_dir / "experiment_results.csv"
    results_df.to_csv(results_path, index=False)

    snapshot = config_to_dict(config)
    snapshot["atlas_cfg_keys"] = sorted(str(key) for key in atlas_cfg.keys())
    (output_dir / "config_snapshot.json").write_text(
        json.dumps(snapshot, indent=2),
        encoding="utf-8",
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
