#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from ml.calibration import apply_calibration, fit_oof_calibrator
from ml.config import config_to_dict, load_atlas_cfg, load_config
from ml.data.bigbind import load_train_holdout_from_bigbind
from ml.evaluate import evaluate_holdout_metrics
from ml.featurize import BigBindFeaturizer
from ml.fpocket_bigbind import merge_fpocket_metrics_on_pocket, precompute_fpocket_for_bigbind_df
from ml.hard_negatives import merge_hard_negatives
from ml.labels import apply_label_smoothing, derive_sample_weight
from ml.metrics_extra import (
    brier_score,
    expected_calibration_error,
    reliability_curve,
    save_reliability_plot,
)
from ml.models import build_model
from ml.pipeline import run_two_stage_pipeline
from ml.pocket_features import POCKET_FEATURE_COLUMNS
from ml.registry import (
    append_registry_index,
    collect_env_versions,
    compute_dataset_hash,
    get_git_sha,
    write_registry_record,
)
from ml.splits import assign_family_split_labels, assign_pocket_similarity_clusters


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "CPU-only BigBind training workflow: load activities, featurize, "
            "train classifier/ranker, and evaluate holdout metrics."
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
    return logging.getLogger("ml.train_and_eval")


def _resolve_run_id(config_run_id: str | None) -> str:
    if config_run_id:
        return config_run_id
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def _normalize_family_name(name: str) -> str:
    family = str(name or "").strip().lower()
    if family == "lgbm":
        return "lightgbm"
    if family == "xgb":
        return "xgboost"
    return family


def _is_rank_family(model_family: str, rank_enabled: bool) -> bool:
    family = _normalize_family_name(model_family)
    if "rank" in family:
        return True
    return bool(rank_enabled and family in {"lightgbm", "xgboost"})


def _as_rank_family(model_family: str) -> str:
    family = _normalize_family_name(model_family)
    if family == "lightgbm":
        return "lightgbm_rank"
    if family == "xgboost":
        return "xgboost_rank"
    return family


def _sigmoid(x: np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    return 1.0 / (1.0 + np.exp(-arr))


def _normalize_model_spec(config) -> tuple[str, dict[str, Any]]:
    family = str(getattr(config, "model_family", "") or "").strip().lower()
    if not family:
        family = str(getattr(config.model, "type", "logreg")).strip().lower()
    family = _normalize_family_name(family)

    params = dict(getattr(config, "model_params", {}) or {})
    if family in {"logreg", "logistic", "logistic_regression"}:
        params.setdefault("C", float(config.model.C))
        params.setdefault("class_weight", config.model.class_weight)
        params.setdefault("max_iter", int(config.model.max_iter))
        params.setdefault("solver", "liblinear")
        params.setdefault("penalty", "l2")
    if _is_rank_family(family, config.rank.enabled):
        family = _as_rank_family(family)
    return family, params


def _group_values(df: pd.DataFrame, group_key: str) -> np.ndarray:
    key = str(group_key).strip().lower()
    if key == "target":
        return df["ex_rec_pdb"].astype(str).to_numpy(dtype=object)
    if key == "pocket":
        return df["pocket"].astype(str).to_numpy(dtype=object)
    return (
        df["ex_rec_pdb"].astype(str) + "::" + df["pocket"].astype(str)
    ).to_numpy(dtype=object)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


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


def _collect_pocket_feature_cache_keys(metrics_df: pd.DataFrame | None) -> list[str]:
    if metrics_df is None or metrics_df.empty:
        return []
    keys: list[str] = []
    for raw in metrics_df.get("pocket_abs_path", pd.Series(dtype=object)).astype(str).tolist():
        token = raw.strip()
        if not token:
            continue
        digest = hashlib.sha1(token.encode("utf-8")).hexdigest()  # noqa: S324
        keys.append(digest)
    return sorted(set(keys))


def _adversarial_validation(
    train_X: sparse.csr_matrix,
    holdout_X: sparse.csr_matrix,
    *,
    seed: int,
) -> dict[str, Any]:
    y_adv = np.concatenate(
        [
            np.zeros(train_X.shape[0], dtype=int),
            np.ones(holdout_X.shape[0], dtype=int),
        ]
    )
    X_adv = sparse.vstack([train_X, holdout_X], format="csr")
    model = LogisticRegression(max_iter=1000, random_state=int(seed))
    model.fit(X_adv, y_adv)
    p = model.predict_proba(X_adv)[:, 1]
    auc = float(roc_auc_score(y_adv, p))
    return {
        "adversarial_auc": auc,
        "train_rows": int(train_X.shape[0]),
        "holdout_rows": int(holdout_X.shape[0]),
        "warning": bool(auc > 0.8),
    }


def _feature_block_masks(featurizer: BigBindFeaturizer) -> dict[str, np.ndarray]:
    names = list(featurizer.feature_names)
    cont_names = list(featurizer.metadata().get("continuous_feature_names", []))
    n_features = len(names)
    mask_lig = np.zeros(n_features, dtype=bool)
    mask_fpocket = np.zeros(n_features, dtype=bool)
    mask_pocket = np.zeros(n_features, dtype=bool)
    mask_morgan = np.zeros(n_features, dtype=bool)

    cont_count = len(cont_names)
    for idx, name in enumerate(names):
        if idx >= cont_count:
            mask_morgan[idx] = True
            continue
        text = str(name)
        if text.startswith("lig_"):
            mask_lig[idx] = True
        if text.startswith("fpocket_"):
            mask_fpocket[idx] = True
        if text in POCKET_FEATURE_COLUMNS:
            mask_pocket[idx] = True
    return {
        "ligand_continuous": mask_lig,
        "morgan_bits": mask_morgan,
        "fpocket_metrics": mask_fpocket,
        "pocket_features": mask_pocket,
    }


def _train_with_block_drop(
    *,
    train_X: sparse.csr_matrix,
    train_y: np.ndarray,
    holdout_X: sparse.csr_matrix,
    holdout_y: np.ndarray,
    keep_mask: np.ndarray,
    model_family: str,
    model_params: dict[str, Any],
    seed: int,
    train_query_groups: np.ndarray | None,
    sample_weight: np.ndarray | None,
) -> dict[str, float | int]:
    model = build_model(model_family=model_family, model_params=model_params, random_seed=seed)
    model.fit(
        train_X[:, keep_mask],
        train_y,
        sample_weight=sample_weight,
        query_groups=train_query_groups,
    )
    p = np.asarray(model.predict_proba(holdout_X[:, keep_mask])[:, 1], dtype=float)
    return evaluate_holdout_metrics(holdout_y, p, fractions=(0.01, 0.02, 0.05, 0.10))


def _drop_feature_tests(
    *,
    train_X: sparse.csr_matrix,
    train_y: np.ndarray,
    holdout_X: sparse.csr_matrix,
    holdout_y: np.ndarray,
    featurizer: BigBindFeaturizer,
    base_metrics: dict[str, float | int],
    model_family: str,
    model_params: dict[str, Any],
    seed: int,
    train_query_groups: np.ndarray | None,
    sample_weight: np.ndarray | None,
) -> pd.DataFrame:
    blocks = _feature_block_masks(featurizer)
    rows: list[dict[str, Any]] = []
    for block_name, drop_mask in blocks.items():
        if not np.any(drop_mask):
            continue
        keep_mask = ~drop_mask
        report = _train_with_block_drop(
            train_X=train_X,
            train_y=train_y,
            holdout_X=holdout_X,
            holdout_y=holdout_y,
            keep_mask=keep_mask,
            model_family=model_family,
            model_params=model_params,
            seed=seed,
            train_query_groups=train_query_groups,
            sample_weight=sample_weight,
        )
        rows.append(
            {
                "dropped_block": block_name,
                "holdout_PR_AUC": float(report.get("PR_AUC", float("nan"))),
                "holdout_EF@1%": float(report.get("EF@1%", float("nan"))),
                "delta_PR_AUC": float(report.get("PR_AUC", float("nan")))
                - float(base_metrics.get("PR_AUC", float("nan"))),
                "delta_EF@1%": float(report.get("EF@1%", float("nan")))
                - float(base_metrics.get("EF@1%", float("nan"))),
            }
        )
    return pd.DataFrame(rows)


def _permutation_importance_continuous(
    *,
    model,
    holdout_X: sparse.csr_matrix,
    holdout_y: np.ndarray,
    featurizer: BigBindFeaturizer,
    max_rows: int,
    seed: int,
) -> pd.DataFrame:
    if holdout_X.shape[0] == 0:
        return pd.DataFrame(columns=["feature", "delta_PR_AUC", "delta_EF@1%"])

    rng = np.random.default_rng(int(seed))
    cont_names = list(featurizer.metadata().get("continuous_feature_names", []))
    cont_count = len(cont_names)
    if cont_count == 0:
        return pd.DataFrame(columns=["feature", "delta_PR_AUC", "delta_EF@1%"])

    rows_idx = np.arange(holdout_X.shape[0], dtype=int)
    if rows_idx.size > int(max_rows):
        rows_idx = rng.choice(rows_idx, size=int(max_rows), replace=False)
    X_sub = holdout_X[rows_idx].toarray()
    y_sub = holdout_y[rows_idx]
    base_p = np.asarray(model.predict_proba(sparse.csr_matrix(X_sub))[:, 1], dtype=float)
    base = evaluate_holdout_metrics(y_sub, base_p, fractions=(0.01, 0.02, 0.05, 0.10))
    base_pr = float(base.get("PR_AUC", float("nan")))
    base_ef1 = float(base.get("EF@1%", float("nan")))

    out_rows: list[dict[str, Any]] = []
    for idx, feat_name in enumerate(cont_names):
        X_perm = X_sub.copy()
        shuffled = X_perm[:, idx].copy()
        rng.shuffle(shuffled)
        X_perm[:, idx] = shuffled
        p = np.asarray(model.predict_proba(sparse.csr_matrix(X_perm))[:, 1], dtype=float)
        report = evaluate_holdout_metrics(y_sub, p, fractions=(0.01, 0.02, 0.05, 0.10))
        out_rows.append(
            {
                "feature": feat_name,
                "delta_PR_AUC": base_pr - float(report.get("PR_AUC", float("nan"))),
                "delta_EF@1%": base_ef1 - float(report.get("EF@1%", float("nan"))),
            }
        )
    df = pd.DataFrame(out_rows)
    if not df.empty:
        df.sort_values(["delta_PR_AUC", "delta_EF@1%"], ascending=False, inplace=True)
    return df


def _stress_split_reports(
    *,
    output_dir: Path,
    holdout_df: pd.DataFrame,
    holdout_y: np.ndarray,
    holdout_prob: np.ndarray,
    config,
) -> None:
    mode = str(config.stress_splits.mode).strip().lower()
    base_report = evaluate_holdout_metrics(holdout_y, holdout_prob, fractions=(0.01, 0.02, 0.05, 0.10))
    pd.DataFrame([base_report]).to_csv(output_dir / "metrics_report_standard.csv", index=False)
    _write_json(output_dir / "metrics_report_standard.json", {k: float(v) if isinstance(v, (int, float)) else v for k, v in base_report.items()})

    if mode == "protein_family":
        family = assign_family_split_labels(
            holdout_df,
            method="heuristic",
            identity_threshold=config.stress_splits.protein_family.identity_threshold,
            seed=config.stress_splits.protein_family.seed,
        )
        group_metrics: list[dict[str, Any]] = []
        for fam in sorted(set(family.astype(str).tolist())):
            mask = family.astype(str) == fam
            if int(mask.sum()) <= 1:
                continue
            rep = evaluate_holdout_metrics(holdout_y[mask.values], holdout_prob[mask.values], fractions=(0.01, 0.02, 0.05, 0.10))
            rep["family_id"] = fam
            group_metrics.append(rep)
        df = pd.DataFrame(group_metrics)
        df.to_csv(output_dir / "metrics_report_protein_family.csv", index=False)
        _write_json(
            output_dir / "metrics_report_protein_family.json",
            {
                "mode": "protein_family",
                "groups": int(len(df)),
                "overall": base_report,
            },
        )
        return

    if mode == "pocket_similarity":
        cols = [col for col in POCKET_FEATURE_COLUMNS if col in holdout_df.columns]
        clusters = assign_pocket_similarity_clusters(
            holdout_df,
            feature_columns=cols,
            n_clusters=config.stress_splits.pocket_similarity.n_clusters,
            seed=config.stress_splits.pocket_similarity.seed,
        )
        group_metrics = []
        for cluster_id in sorted(set(clusters.astype(int).tolist())):
            mask = clusters.astype(int) == int(cluster_id)
            if int(mask.sum()) <= 1:
                continue
            rep = evaluate_holdout_metrics(holdout_y[mask.values], holdout_prob[mask.values], fractions=(0.01, 0.02, 0.05, 0.10))
            rep["cluster_id"] = int(cluster_id)
            group_metrics.append(rep)
        df = pd.DataFrame(group_metrics)
        df.to_csv(output_dir / "metrics_report_pocket_similarity.csv", index=False)
        _write_json(
            output_dir / "metrics_report_pocket_similarity.json",
            {
                "mode": "pocket_similarity",
                "groups": int(len(df)),
                "overall": base_report,
            },
        )


def run_pipeline(config_path: str) -> tuple[Path, dict[str, float | int]]:
    logger = _setup_logging()
    config = load_config(config_path)
    config_payload = config_to_dict(config)
    atlas_cfg = load_atlas_cfg(config)
    np.random.seed(config.random_seed)
    run_id = _resolve_run_id(config.run_id)
    output_dir = Path(__file__).resolve().parent / "outputs" / run_id
    output_dir.mkdir(parents=True, exist_ok=True)

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
    train_df = dataset.train_df.copy()
    holdout_df = dataset.holdout_df.copy()

    model_family, model_params = _normalize_model_spec(config)
    is_ranker = _is_rank_family(model_family, config.rank.enabled)
    group_key = str(config.rank.group_key if config.rank.enabled else "target_pocket")

    if config.hard_negatives.enabled:
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
        logger.info("[ml] Hard-negative mining applied: train rows now=%d", len(train_df))

    fpocket_metrics_df: pd.DataFrame | None = None
    if config.features.pocket_fpocket or config.features.pocket_features:
        combined = pd.concat([train_df, holdout_df], ignore_index=False)
        fpocket_metrics_df = precompute_fpocket_for_bigbind_df(
            df=combined,
            bigbind_root=dataset.bigbind_root,
            atlas_cfg=atlas_cfg,
            run_dir=output_dir,
        )
        train_df = merge_fpocket_metrics_on_pocket(train_df, fpocket_metrics_df)
        holdout_df = merge_fpocket_metrics_on_pocket(holdout_df, fpocket_metrics_df)

    featurizer = BigBindFeaturizer(
        bigbind_root=dataset.bigbind_root,
        features=config.features,
        atlas_cfg=atlas_cfg,
        fpocket_center_columns=config.fpocket_center_columns,
        fpocket_variant_column=config.fpocket_variant_column,
        fpocket_ph_column=config.fpocket_ph_column,
        fpocket_centers_by_pdb=config.fpocket_centers_by_pdb,
        logger=logger,
    )
    train_features = featurizer.transform(train_df, audit_dir=output_dir, audit_tag="train")
    holdout_features = featurizer.transform(holdout_df, audit_dir=output_dir, audit_tag="holdout")

    train_used = train_df.loc[train_features.source_index].copy()
    holdout_used = holdout_df.loc[holdout_features.source_index].copy()
    train_query_groups = (
        _group_values(train_used, group_key=group_key) if is_ranker else None
    )

    sample_weight = (
        derive_sample_weight(train_used, weight_cap=config.labels.weight_cap)
        if config.labels.use_sample_weights
        else None
    )
    y_train_fit = train_features.y.astype(float if config.labels.smoothing_eps > 0 else int)
    if config.labels.smoothing_eps > 0 and model_family not in {"logreg", "logistic", "logistic_regression"}:
        y_train_fit = apply_label_smoothing(train_features.y, config.labels.smoothing_eps)

    def _builder(seed_value: int):
        return build_model(
            model_family=model_family,
            model_params=model_params,
            random_seed=seed_value,
        )

    calibration_oof_df = pd.DataFrame(
        columns=[
            "row_id",
            "source_index",
            "fold_id",
            "group",
            "y_true",
            "p_uncalibrated",
            "used_for_calibration",
            "lig_smiles",
            "active",
            "ex_rec_pdb",
            "pocket",
            "murcko_scaffold",
        ]
    )
    calibrator = None
    if config.calibration_enabled:
        calibrator, oof_df = fit_oof_calibrator(
            X=train_features.X,
            y=train_features.y,
            groups=np.asarray(train_features.murcko_scaffolds, dtype=object),
            base_model_builder=_builder,
            method=config.calibration_method,
            cv_folds=config.calibration_cv_folds,
            seed=config.calibration_seed,
            fit_query_groups=train_query_groups,
            fit_sample_weights=sample_weight,
        )
        calibration_oof_df = oof_df.copy()
        calibration_oof_df["source_index"] = np.asarray(train_features.source_index, dtype=int)
        calibration_oof_df["lig_smiles"] = train_used["lig_smiles"].astype(str).values
        calibration_oof_df["active"] = train_used["active"].astype(int).values
        calibration_oof_df["ex_rec_pdb"] = train_used["ex_rec_pdb"].astype(str).values
        calibration_oof_df["pocket"] = train_used["pocket"].astype(str).values
        calibration_oof_df["murcko_scaffold"] = np.asarray(train_features.murcko_scaffolds, dtype=object)

    model = _builder(config.random_seed)
    model.fit(
        train_features.X,
        np.asarray(y_train_fit),
        sample_weight=sample_weight,
        query_groups=train_query_groups,
    )

    train_prob_uncal = np.asarray(model.predict_proba(train_features.X)[:, 1], dtype=float)
    holdout_prob_uncal = np.asarray(model.predict_proba(holdout_features.X)[:, 1], dtype=float)
    holdout_raw_score = np.asarray(model.raw_score(holdout_features.X), dtype=float)

    train_prob_cal = apply_calibration(calibrator, train_prob_uncal) if calibrator is not None else train_prob_uncal
    holdout_prob_cal = apply_calibration(calibrator, holdout_prob_uncal) if calibrator is not None else holdout_prob_uncal

    stage1_df = pd.DataFrame()
    stage2_df = pd.DataFrame()
    pipeline_summary: dict[str, Any] | None = None
    if config.pipeline.enabled:
        stage1_family = _normalize_family_name(config.pipeline.stage1.model_family)
        stage2_family = _normalize_family_name(config.pipeline.stage2.model_family)
        if _is_rank_family(stage2_family, True):
            stage2_family = _as_rank_family(stage2_family)
        two_stage = run_two_stage_pipeline(
            train_X=train_features.X,
            train_y=train_features.y,
            holdout_X=holdout_features.X,
            train_df=train_used,
            holdout_df=holdout_used,
            stage1_model_family=stage1_family,
            stage1_model_params=config.pipeline.stage1.model_params,
            stage2_model_family=stage2_family,
            stage2_model_params=config.pipeline.stage2.model_params,
            group_key=config.pipeline.stage2.group_key,
            keep_top_pct=config.pipeline.stage1.keep_top_pct,
            keep_prob_ge=config.pipeline.stage1.keep_prob_ge,
            keep_within_group=config.pipeline.stage1.keep_within_group,
            random_seed=config.random_seed,
        )
        holdout_raw_score = np.asarray(two_stage.stage2_rank_score_holdout, dtype=float)
        holdout_prob_cal = _sigmoid(holdout_raw_score)
        stage1_df = pd.DataFrame(
            {
                "row_id": np.arange(len(holdout_used), dtype=int),
                "lig_smiles": holdout_used["lig_smiles"].astype(str).values,
                "active": holdout_features.y.astype(int),
                "ex_rec_pdb": holdout_used["ex_rec_pdb"].astype(str).values,
                "pocket": holdout_used["pocket"].astype(str).values,
                "stage1_prob": two_stage.stage1_prob_holdout.astype(float),
                "kept_for_stage2": two_stage.stage1_keep_holdout.astype(int),
            }
        )
        stage2_df = pd.DataFrame(
            {
                "row_id": np.arange(len(holdout_used), dtype=int),
                "stage2_rank_score": holdout_raw_score.astype(float),
                "prob_calibrated": holdout_prob_cal.astype(float),
            }
        )
        pipeline_summary = {
            "enabled": True,
            "stage1_kept_train_fraction": float(np.mean(two_stage.stage1_keep_train)),
            "stage1_kept_holdout_fraction": float(np.mean(two_stage.stage1_keep_holdout)),
            "stage1_keep_top_pct": config.pipeline.stage1.keep_top_pct,
            "stage1_keep_prob_ge": config.pipeline.stage1.keep_prob_ge,
            "stage2_model_family": stage2_family,
        }

    metrics = evaluate_holdout_metrics(
        holdout_features.y,
        holdout_prob_cal,
        fractions=(0.01, 0.02, 0.05, 0.10),
    )

    calibration_metrics: dict[str, Any] = {
        "calibration_enabled": bool(config.calibration_enabled),
        "calibration_method": str(config.calibration_method),
        "calibration_cv_folds": int(config.calibration_cv_folds),
    }
    train_curve = reliability_curve(train_features.y, train_prob_cal, n_bins=config.metrics_ece_bins)
    holdout_curve = reliability_curve(holdout_features.y, holdout_prob_cal, n_bins=config.metrics_ece_bins)
    save_reliability_plot(output_dir / "reliability_plot.png", holdout_curve)

    train_brier = brier_score(train_features.y, train_prob_cal)
    holdout_brier = brier_score(holdout_features.y, holdout_prob_cal)
    train_ece = expected_calibration_error(train_features.y, train_prob_cal, n_bins=config.metrics_ece_bins)
    holdout_ece = expected_calibration_error(holdout_features.y, holdout_prob_cal, n_bins=config.metrics_ece_bins)
    calibration_metrics["train"] = {
        "brier_score": float(train_brier),
        "ece": float(train_ece),
        "bin_centers": train_curve[0].astype(float).tolist(),
        "fraction_positives": train_curve[1].astype(float).tolist(),
        "mean_predicted": train_curve[2].astype(float).tolist(),
        "bin_counts": train_curve[3].astype(int).tolist(),
    }
    calibration_metrics["holdout"] = {
        "brier_score": float(holdout_brier),
        "ece": float(holdout_ece),
        "bin_centers": holdout_curve[0].astype(float).tolist(),
        "fraction_positives": holdout_curve[1].astype(float).tolist(),
        "mean_predicted": holdout_curve[2].astype(float).tolist(),
        "bin_counts": holdout_curve[3].astype(int).tolist(),
    }
    if config.metrics_report_brier:
        metrics["Brier"] = float(holdout_brier)
    if config.metrics_report_ece:
        metrics["ECE"] = float(holdout_ece)

    holdout_predictions = holdout_used[
        [col for col in ("split", "lig_smiles", "active", "pocket", "ex_rec_pdb") if col in holdout_used.columns]
    ].copy()
    holdout_predictions["murcko_scaffold"] = holdout_features.murcko_scaffolds
    holdout_predictions["p_active"] = holdout_prob_cal
    holdout_predictions["raw_score"] = holdout_raw_score
    holdout_predictions["p_active_uncalibrated"] = holdout_prob_uncal
    holdout_predictions.to_csv(output_dir / "holdout_predictions.csv", index=False)

    calibrated_predictions = pd.DataFrame(
        {
            "row_id": np.arange(holdout_features.X.shape[0], dtype=int),
            "lig_smiles": holdout_used["lig_smiles"].astype(str).values,
            "active": holdout_features.y.astype(int),
            "ex_rec_pdb": holdout_used["ex_rec_pdb"].astype(str).values,
            "pocket": holdout_used["pocket"].astype(str).values,
            "scaffold": np.asarray(holdout_features.murcko_scaffolds, dtype=object),
            "raw_score": holdout_raw_score,
            "prob_uncalibrated": holdout_prob_uncal,
            "prob_calibrated": holdout_prob_cal,
        }
    )

    # Core outputs
    joblib.dump(model, output_dir / "model.joblib")
    snapshot_text = json.dumps(config_payload, indent=2)
    (output_dir / "config_snapshot.txt").write_text(snapshot_text, encoding="utf-8")
    (output_dir / "config_snapshot.json").write_text(snapshot_text, encoding="utf-8")
    _write_json(output_dir / "featurizer_metadata.json", featurizer.metadata())
    _write_json(output_dir / "metrics_report.json", metrics)
    pd.DataFrame([metrics]).to_csv(output_dir / "metrics_report.csv", index=False)
    calibrated_predictions.to_csv(output_dir / "calibrated_holdout_predictions.csv", index=False)
    calibration_oof_df.to_csv(output_dir / "calibration_oof_train.csv", index=False)
    _write_json(output_dir / "metrics_calibration.json", calibration_metrics)

    if config.pipeline.enabled:
        stage1_df.to_csv(output_dir / "stage1_predictions_holdout.csv", index=False)
        stage2_df.to_csv(output_dir / "stage2_rank_scores_holdout.csv", index=False)
        summary_payload = dict(pipeline_summary or {})
        summary_payload["metrics"] = metrics
        _write_json(output_dir / "pipeline_summary.json", summary_payload)

    # Diagnostics
    if config.diagnostics.adversarial_validation:
        adv = _adversarial_validation(train_features.X, holdout_features.X, seed=config.random_seed)
        _write_json(output_dir / "adversarial_validation.json", adv)
    if config.diagnostics.drop_feature_tests:
        drop_df = _drop_feature_tests(
            train_X=train_features.X,
            train_y=train_features.y,
            holdout_X=holdout_features.X,
            holdout_y=holdout_features.y,
            featurizer=featurizer,
            base_metrics=metrics,
            model_family=model_family,
            model_params=model_params,
            seed=config.random_seed,
            train_query_groups=train_query_groups,
            sample_weight=sample_weight,
        )
        drop_df.to_csv(output_dir / "drop_feature_report.csv", index=False)
        _write_json(
            output_dir / "drop_feature_report.json",
            {"rows": drop_df.to_dict(orient="records")},
        )
    if config.diagnostics.permutation_importance and "rank" not in model_family:
        perm_df = _permutation_importance_continuous(
            model=model,
            holdout_X=holdout_features.X,
            holdout_y=holdout_features.y,
            featurizer=featurizer,
            max_rows=config.diagnostics.permutation_max_rows,
            seed=config.random_seed,
        )
        perm_df.to_csv(output_dir / "permutation_importance.csv", index=False)

    _stress_split_reports(
        output_dir=output_dir,
        holdout_df=holdout_used,
        holdout_y=holdout_features.y,
        holdout_prob=holdout_prob_cal,
        config=config,
    )

    if config.registry_enabled:
        repo_root = Path(__file__).resolve().parents[1]
        dataset_hash = compute_dataset_hash(
            train_df=train_used,
            holdout_df=holdout_used,
            feature_columns=featurizer.feature_names,
            config_snapshot=config_payload,
        )
        registry_record: dict[str, Any] = {
            "run_id": run_id,
            "run_dir": str(output_dir),
            "dataset_hash": dataset_hash,
            "git_sha": get_git_sha(repo_root),
            "model_family": model_family,
            "model_params": model_params,
            "rank_enabled": bool(config.rank.enabled),
            "rank_group_key": str(config.rank.group_key),
            "calibration_enabled": bool(config.calibration_enabled),
            "calibration_method": str(config.calibration_method),
            "hard_negatives": config_payload.get("hard_negatives", {}),
            "labels": config_payload.get("labels", {}),
            "metrics": metrics,
            "fpocket_cache_keys_used": _collect_fpocket_cache_keys(fpocket_metrics_df),
            "pocket_feature_cache_keys": _collect_pocket_feature_cache_keys(fpocket_metrics_df),
            "pocket_feature_columns": list(POCKET_FEATURE_COLUMNS),
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
                "model_family": model_family,
                "rank_enabled": int(bool(config.rank.enabled)),
                "calibration_enabled": int(bool(config.calibration_enabled)),
                "calibration_method": str(config.calibration_method),
                "holdout_PR_AUC": float(metrics.get("PR_AUC", float("nan"))),
                "holdout_EF@1%": float(metrics.get("EF@1%", float("nan"))),
            },
        )

    return output_dir, metrics


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    output_dir, metrics = run_pipeline(args.config)
    print(
        "[ml] "
        f"saved={output_dir} "
        f"N={metrics['N']} "
        f"n_actives={metrics['n_actives']} "
        f"PR_AUC={metrics['PR_AUC']:.6f} "
        f"EF@1%={metrics['EF@1%']:.6f} "
        f"EF@5%={metrics['EF@5%']:.6f}"
    )


if __name__ == "__main__":
    main()
