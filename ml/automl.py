#!/usr/bin/env python3
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import itertools
import logging
import threading
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any, Iterable, cast
import json

import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import GroupKFold

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from ml.calibration import apply_calibration, fit_oof_calibrator
from ml.config import (
    FeaturesConfig,
    config_to_dict,
    load_atlas_cfg,
    load_config,
    load_config_payload,
)
from ml.data.bigbind import load_train_holdout_from_bigbind
from ml.evaluate import DEFAULT_FRACTIONS, evaluate_holdout_metrics
from ml.featurize import BigBindFeaturizer, FeaturizedRows
from ml.fpocket_bigbind import merge_fpocket_metrics_on_pocket, precompute_fpocket_for_bigbind_df
from ml.hard_negatives import merge_hard_negatives
from ml.labels import derive_sample_weight
from ml.models import build_model
from ml.output_views import write_truncated_views_for_run_dir
from ml.registry import (
    append_registry_index,
    collect_env_versions,
    compute_dataset_hash,
    get_git_sha,
    write_registry_record,
)


_FEATURE_VARIANT_ORDER = (
    "baseline",
    "ligand_only",
    "pocket_only",
    "no_fp",
    "add_ligand_extras",
    "add_pocket_fpocket",
    "add_both",
)
_DEFAULT_MODEL_FAMILIES = ("logreg", "lightgbm", "xgboost")

FeatureCacheEntry = tuple[
    FeaturesConfig,
    FeaturizedRows,  # train
    FeaturizedRows,  # validation
    FeaturizedRows,  # test
    dict[str, int],  # fpocket stats
    np.ndarray | None,  # train query groups
    np.ndarray | None,  # train sample weights
]


def _setup_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    return logging.getLogger("ml.automl")


def _resolve_run_id(config_run_id: str | None) -> str:
    if config_run_id:
        return config_run_id
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def _load_automl_payload(config_path: str | Path) -> dict[str, Any]:
    payload = load_config_payload(config_path)
    automl_payload = payload.get("automl")
    return automl_payload if isinstance(automl_payload, dict) else {}


def _normalize_name_list(raw: Any, *, default: Iterable[str]) -> list[str]:
    if raw is None:
        return list(default)
    if isinstance(raw, str):
        names = [part.strip() for part in raw.split(",") if part.strip()]
        return names or list(default)
    if isinstance(raw, (list, tuple)):
        names = [str(part).strip() for part in raw if str(part).strip()]
        return names or list(default)
    return list(default)


def _normalize_grid_space(raw: Any, *, fallback: dict[str, list[Any]]) -> dict[str, list[Any]]:
    if not isinstance(raw, dict):
        return {key: list(values) for key, values in fallback.items()}
    merged: dict[str, list[Any]] = {key: list(values) for key, values in fallback.items()}
    for key, values in raw.items():
        if isinstance(values, list):
            merged[str(key)] = list(values)
        elif isinstance(values, tuple):
            merged[str(key)] = list(values)
        else:
            merged[str(key)] = [values]
    return merged


def _default_grid_spaces() -> dict[str, dict[str, list[Any]]]:
    spaces: dict[str, dict[str, list[Any]]] = {
        "logreg": {
            "C": [0.01, 0.1, 1.0, 10.0],
            "penalty": ["l2"],
            "solver": ["liblinear"],
            "class_weight": ["balanced"],
        },
        "lightgbm": {
            "learning_rate": [0.03, 0.1],
            "num_leaves": [31, 63],
            "max_depth": [-1, 6],
            "min_child_samples": [20, 50],
            "n_estimators": [200],
            "subsample": [1.0],
            "colsample_bytree": [1.0],
            "class_weight": ["balanced"],
        },
        "xgboost": {
            "max_depth": [4, 6],
            "eta": [0.03, 0.1],
            "n_estimators": [200],
            "subsample": [0.8, 1.0],
            "colsample_bytree": [0.8, 1.0],
        },
        "lightgbm_rank": {
            "objective": ["lambdarank"],
            "learning_rate": [0.03, 0.1],
            "num_leaves": [31, 63],
            "min_child_samples": [20],
            "n_estimators": [200],
        },
        "xgboost_rank": {
            "objective": ["rank:pairwise"],
            "max_depth": [4, 6],
            "eta": [0.03, 0.1],
            "n_estimators": [200],
            "subsample": [0.8, 1.0],
            "colsample_bytree": [0.8, 1.0],
        },
        # Legacy compatibility family
        "hgb": {
            "max_depth": [3, 6, None],
            "learning_rate": [0.03, 0.1],
            "max_leaf_nodes": [31, 63],
            "min_samples_leaf": [20, 50],
        },
    }
    return spaces


def _coerce_positive_int(value: Any) -> int | None:
    if value in (None, "", "none", "null"):
        return None
    try:
        parsed = int(value)
    except Exception:
        return None
    return parsed if parsed > 0 else None


def _resolve_parallel_workers(
    *,
    automl_payload: dict[str, Any],
    atlas_cfg: dict[str, Any],
) -> int:
    for key in ("max_workers", "workers", "n_jobs"):
        parsed = _coerce_positive_int(automl_payload.get(key))
        if parsed is not None:
            return parsed
    for key in ("ML_AUTOML_WORKERS", "ML_WORKERS", "CPU"):
        parsed = _coerce_positive_int(atlas_cfg.get(key))
        if parsed is not None:
            return parsed
    return 1


def _resolve_fpocket_workers(atlas_cfg: dict[str, Any]) -> int:
    for key in ("ML_FPOCKET_WORKERS", "ML_WORKERS", "CPU"):
        parsed = _coerce_positive_int(atlas_cfg.get(key))
        if parsed is not None:
            return parsed
    return 1


def _build_feature_variants(base_features: FeaturesConfig) -> dict[str, FeaturesConfig]:
    baseline = replace(base_features)
    return {
        "baseline": baseline,
        "ligand_only": replace(
            baseline,
            pocket_fpocket=False,
            ligand_extra_descriptors=False,
        ),
        "pocket_only": FeaturesConfig(
            ligand_descriptors=False,
            ligand_extra_descriptors=False,
            ligand_morgan_fp_bits=0,
            pocket_fpocket=True,
            pocket_features=baseline.pocket_features,
            vina_score=False,
        ),
        "no_fp": replace(baseline, ligand_morgan_fp_bits=0),
        "add_ligand_extras": replace(baseline, ligand_extra_descriptors=True),
        "add_pocket_fpocket": replace(baseline, pocket_fpocket=True),
        "add_both": replace(baseline, ligand_extra_descriptors=True, pocket_fpocket=True),
    }


def _model_families_from_payload(automl_payload: dict[str, Any]) -> list[str]:
    default = list(_DEFAULT_MODEL_FAMILIES) + ["lightgbm_rank", "xgboost_rank"]
    names = _normalize_name_list(automl_payload.get("model_families"), default=default)
    normalized: list[str] = []
    for raw_name in names:
        family = str(raw_name).strip().lower()
        if family == "lgbm":
            family = "lightgbm"
        if family == "xgb":
            family = "xgboost"
        if family == "lgbm_rank":
            family = "lightgbm_rank"
        if family == "xgb_rank":
            family = "xgboost_rank"
        normalized.append(family)
    valid = set(_default_grid_spaces().keys())
    selected = [name for name in normalized if name in valid]
    available: list[str] = []
    for family in selected:
        try:
            build_model(model_family=family, model_params={}, random_seed=0)
            available.append(family)
        except RuntimeError:
            continue
        except Exception:
            available.append(family)
    return available or ["logreg"]


def _feature_variant_names_from_payload(automl_payload: dict[str, Any]) -> list[str]:
    names = _normalize_name_list(
        automl_payload.get("feature_variants"),
        default=_FEATURE_VARIANT_ORDER,
    )
    valid = set(_FEATURE_VARIANT_ORDER)
    selected = [name for name in names if name in valid]
    return selected or list(_FEATURE_VARIANT_ORDER)


def _calibration_options_from_payload(
    automl_payload: dict[str, Any],
    *,
    default_enabled: bool,
    default_method: str,
) -> list[tuple[bool, str]]:
    enabled_raw = automl_payload.get("calibration_enabled")
    method_raw = automl_payload.get("calibration_methods")

    if enabled_raw is None:
        enabled_values = [bool(default_enabled)]
    elif isinstance(enabled_raw, (list, tuple)):
        enabled_values = [bool(value) for value in enabled_raw]
    else:
        enabled_values = [bool(enabled_raw)]

    methods = _normalize_name_list(method_raw, default=[default_method])
    method_values = []
    for method in methods:
        name = str(method).strip().lower()
        if name in {"sigmoid", "isotonic"} and name not in method_values:
            method_values.append(name)
    if not method_values:
        method_values = [str(default_method).strip().lower() or "sigmoid"]

    options: list[tuple[bool, str]] = []
    for enabled in enabled_values:
        if enabled:
            for method in method_values:
                options.append((True, method))
        else:
            options.append((False, "sigmoid"))
    dedup = []
    seen: set[tuple[bool, str]] = set()
    for option in options:
        if option not in seen:
            seen.add(option)
            dedup.append(option)
    return dedup


def _group_splits(
    murcko_scaffolds: list[str],
    y: np.ndarray,
    n_splits_requested: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    groups = np.asarray(
        [
            text if text else f"__empty_scaffold_{idx}"
            for idx, text in enumerate(str(x).strip() for x in murcko_scaffolds)
        ],
        dtype=object,
    )
    unique_groups = np.unique(groups)
    n_splits = min(max(2, int(n_splits_requested)), int(unique_groups.shape[0]), int(y.shape[0]))
    if n_splits < 2:
        return []
    splitter = GroupKFold(n_splits=n_splits)
    return list(splitter.split(np.zeros((y.shape[0], 1)), y, groups))


def _is_finite_number(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except Exception:
        return False


def _fit_model(
    *,
    model_family: str,
    model_params: dict[str, Any],
    X_train: sparse.csr_matrix,
    y_train: np.ndarray,
    query_groups: np.ndarray | None = None,
    sample_weight: np.ndarray | None = None,
    random_seed: int,
) -> Any:
    family = str(model_family).strip().lower()
    if family == "lgbm":
        family = "lightgbm"
    if family == "xgb":
        family = "xgboost"
    model = build_model(
        model_family=family,
        model_params=model_params,
        random_seed=random_seed,
    )
    model.fit(
        X_train,
        y_train,
        query_groups=query_groups,
        sample_weight=sample_weight,
    )
    return model


def _predict_p_active(model: Any, X: sparse.csr_matrix) -> np.ndarray:
    prob = model.predict_proba(X)
    if prob.ndim != 2 or prob.shape[1] < 2:
        raise ValueError("Expected predict_proba output with at least 2 columns.")
    return np.asarray(prob[:, 1], dtype=float)


def _aggregate_fold_reports(fold_reports: list[dict[str, float | int]]) -> dict[str, float | int]:
    if not fold_reports:
        return {"inner_folds": 0}

    aggregated: dict[str, float | int] = {"inner_folds": len(fold_reports)}
    keys = sorted(
        {
            key
            for report in fold_reports
            for key, value in report.items()
            if isinstance(value, (int, float))
        }
    )
    for key in keys:
        values = np.asarray(
            [float(report.get(key, float("nan"))) for report in fold_reports],
            dtype=float,
        )
        aggregated[f"inner_mean_{key}"] = float(np.nanmean(values))
        aggregated[f"inner_std_{key}"] = float(np.nanstd(values))
    return aggregated


def _inner_cv_metrics(
    *,
    X: sparse.csr_matrix,
    y: np.ndarray,
    murcko_scaffolds: list[str],
    model_family: str,
    model_params: dict[str, Any],
    calibration_enabled: bool,
    calibration_method: str,
    calibration_cv_folds: int,
    calibration_seed: int,
    random_seed: int,
    n_splits_requested: int,
    ranking_query_groups: np.ndarray | None = None,
    sample_weight: np.ndarray | None = None,
) -> dict[str, float | int]:
    splits = _group_splits(murcko_scaffolds, y, n_splits_requested)
    fold_reports: list[dict[str, float | int]] = []
    for fold_idx, (train_idx, val_idx) in enumerate(splits):
        y_train_fold = y[train_idx]
        y_val_fold = y[val_idx]
        if np.unique(y_train_fold).size < 2:
            continue
        if np.unique(y_val_fold).size < 2:
            continue
        model = _fit_model(
            model_family=model_family,
            model_params=model_params,
            X_train=X[train_idx],
            y_train=y_train_fold,
            query_groups=(
                ranking_query_groups[train_idx] if ranking_query_groups is not None else None
            ),
            sample_weight=sample_weight[train_idx] if sample_weight is not None else None,
            random_seed=random_seed + fold_idx,
        )
        p_uncal = _predict_p_active(model, X[val_idx])
        p_active = p_uncal
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
            p_active = apply_calibration(calibrator, p_uncal)
        fold_reports.append(
            evaluate_holdout_metrics(
                y_val_fold,
                p_active,
                fractions=DEFAULT_FRACTIONS,
            )
        )
    return _aggregate_fold_reports(fold_reports)


def _metric_for_sort(value: Any) -> float:
    try:
        as_float = float(value)
    except Exception:
        return float("-inf")
    if not np.isfinite(as_float):
        return float("-inf")
    return as_float


def _trial_sort_key(record: dict[str, Any]) -> tuple[float, float, float]:
    inner_metrics = record.get("inner_metrics", {})
    if isinstance(inner_metrics, dict):
        inner_ef1 = inner_metrics.get("inner_mean_EF@1%")
        inner_pr_auc = inner_metrics.get("inner_mean_PR_AUC")
    else:
        inner_ef1 = None
        inner_pr_auc = None
    return (
        _metric_for_sort(record.get("inner_mean_EF@1%", inner_ef1)),
        _metric_for_sort(record.get("inner_mean_PR_AUC", inner_pr_auc)),
        -float(record.get("trial_id", 0)),
    )


def _resolve_ligand_id_column(df: pd.DataFrame) -> pd.Series:
    for col in ("lig_id", "lig_file", "ligand_id", "lig_smiles"):
        if col in df.columns:
            return df[col].astype(str)
    return pd.Series([str(i) for i in df.index], index=df.index, dtype=str)


def _compute_query_groups(df: pd.DataFrame, group_key: str) -> np.ndarray:
    key = str(group_key).strip().lower()
    if key == "target":
        return df["ex_rec_pdb"].astype(str).to_numpy(dtype=object)
    if key == "pocket":
        return df["pocket"].astype(str).to_numpy(dtype=object)
    return (df["ex_rec_pdb"].astype(str) + "::" + df["pocket"].astype(str)).to_numpy(dtype=object)


def _collect_fpocket_cache_keys_from_df(df: pd.DataFrame | None) -> list[str]:
    if df is None or df.empty:
        return []
    keys: list[str] = []
    for _, row in df.iterrows():
        pocket = str(row.get("pocket", "")).strip()
        receptor = str(row.get("receptor_pdb", "")).strip()
        cx = row.get("pocket_center_x")
        cy = row.get("pocket_center_y")
        cz = row.get("pocket_center_z")
        keys.append(f"{pocket}|{receptor}|{cx}|{cy}|{cz}")
    return sorted(set(keys))


def _evaluate_holdout_trial(
    *,
    trial_record: dict[str, Any],
    train_rows: FeaturizedRows,
    holdout_rows: FeaturizedRows,
    holdout_df: pd.DataFrame,
    random_seed: int,
    calibration_cv_folds: int,
    calibration_seed: int,
    train_query_groups: np.ndarray | None = None,
    sample_weight: np.ndarray | None = None,
) -> tuple[dict[str, float | int], pd.DataFrame]:
    model = _fit_model(
        model_family=str(trial_record["model_family"]),
        model_params=dict(trial_record["model_params"]),
        X_train=train_rows.X,
        y_train=train_rows.y,
        query_groups=train_query_groups,
        sample_weight=sample_weight,
        random_seed=random_seed,
    )
    holdout_p_uncal = _predict_p_active(model, holdout_rows.X)
    holdout_p_active = holdout_p_uncal
    if bool(trial_record.get("calibration_enabled", False)):

        def _builder(seed_value: int):
            return build_model(
                model_family=str(trial_record["model_family"]),
                model_params=dict(trial_record["model_params"]),
                random_seed=seed_value,
            )

        calibrator, _ = fit_oof_calibrator(
            X=train_rows.X,
            y=train_rows.y,
            groups=np.asarray(train_rows.murcko_scaffolds, dtype=object),
            base_model_builder=_builder,
            method=str(trial_record.get("calibration_method", "sigmoid")),
            cv_folds=int(calibration_cv_folds),
            seed=int(calibration_seed),
            fit_query_groups=train_query_groups,
            fit_sample_weights=sample_weight,
        )
        holdout_p_active = apply_calibration(calibrator, holdout_p_uncal)

    holdout_metrics = evaluate_holdout_metrics(
        holdout_rows.y,
        holdout_p_active,
        fractions=DEFAULT_FRACTIONS,
    )

    holdout_used = holdout_df.loc[holdout_rows.source_index].copy()
    ligand_id = _resolve_ligand_id_column(holdout_used)
    pred_df = pd.DataFrame(
        {
            "lig_id": ligand_id.values,
            "p_active": holdout_p_active,
            "p_active_uncalibrated": holdout_p_uncal,
            "active": holdout_rows.y,
        }
    )
    pred_df["ranks"] = (
        pred_df["p_active"].rank(method="first", ascending=False).astype(int)
    )
    pred_df.sort_values("ranks", inplace=True)
    return holdout_metrics, pred_df


def _materialize_trial_row(record: dict[str, Any]) -> dict[str, Any]:
    row = {
        "trial_id": int(record["trial_id"]),
        "mode": str(record["mode"]),
        "status": str(record["status"]),
        "feature_variant": str(record["feature_variant"]),
        "model_family": str(record["model_family"]),
        "calibration_enabled": int(bool(record.get("calibration_enabled", False))),
        "calibration_method": str(record.get("calibration_method", "sigmoid")),
        "model_params_json": json.dumps(record["model_params"], sort_keys=True),
        "features_json": json.dumps(record["features"], sort_keys=True),
        "error": str(record.get("error") or ""),
    }
    for key, value in record.get("inner_metrics", {}).items():
        row[key] = value
    return row


def _iter_grid_param_sets(space: dict[str, list[Any]]) -> Iterable[dict[str, Any]]:
    keys = list(space.keys())
    value_lists = [space[key] for key in keys]
    for combo in itertools.product(*value_lists):
        yield {key: value for key, value in zip(keys, combo)}


def _run_grid_search(
    *,
    feature_variant_names: list[str],
    feature_variants: dict[str, FeaturesConfig],
    model_families: list[str],
    calibration_options: list[tuple[bool, str]],
    grid_spaces: dict[str, dict[str, list[Any]]],
    feature_cache: dict[str, FeatureCacheEntry],
    random_seed: int,
    inner_scaffold_folds: int,
    calibration_cv_folds: int,
    calibration_seed: int,
    parallel_workers: int = 1,
) -> list[dict[str, Any]]:
    trial_specs: list[tuple[int, str, str, dict[str, Any], bool, str]] = []
    trial_id = 0
    for feature_variant_name in feature_variant_names:
        for model_family in model_families:
            if model_family not in grid_spaces:
                continue
            for params in _iter_grid_param_sets(grid_spaces[model_family]):
                for calibration_enabled, calibration_method in calibration_options:
                    trial_specs.append(
                        (
                            trial_id,
                            feature_variant_name,
                            model_family,
                            dict(params),
                            bool(calibration_enabled),
                            str(calibration_method),
                        )
                    )
                    trial_id += 1

    def _evaluate_trial(spec: tuple[int, str, str, dict[str, Any], bool, str]) -> dict[str, Any]:
        (
            trial_id_local,
            feature_variant_name,
            model_family,
            params,
            calibration_enabled,
            calibration_method,
        ) = spec
        train_rows = feature_cache[feature_variant_name][1]
        train_query_groups = feature_cache[feature_variant_name][5]
        train_sample_weight = feature_cache[feature_variant_name][6]
        record: dict[str, Any] = {
            "trial_id": trial_id_local,
            "mode": "grid",
            "status": "ok",
            "feature_variant": feature_variant_name,
            "features": asdict(feature_variants[feature_variant_name]),
            "model_family": model_family,
            "model_params": dict(params),
            "calibration_enabled": calibration_enabled,
            "calibration_method": calibration_method,
            "inner_metrics": {},
            "error": "",
        }
        try:
            inner_metrics = _inner_cv_metrics(
                X=train_rows.X,
                y=train_rows.y,
                murcko_scaffolds=train_rows.murcko_scaffolds,
                model_family=model_family,
                model_params=params,
                calibration_enabled=calibration_enabled,
                calibration_method=calibration_method,
                calibration_cv_folds=calibration_cv_folds,
                calibration_seed=calibration_seed,
                random_seed=random_seed,
                n_splits_requested=inner_scaffold_folds,
                ranking_query_groups=train_query_groups,
                sample_weight=train_sample_weight,
            )
            if int(inner_metrics.get("inner_folds", 0)) <= 0:
                record["status"] = "skipped"
                record["error"] = "no_valid_inner_folds"
            record["inner_metrics"] = inner_metrics
        except Exception as exc:
            record["status"] = "failed"
            record["error"] = f"{type(exc).__name__}: {exc}"
        return record

    workers = min(max(1, int(parallel_workers)), max(1, len(trial_specs)))
    if workers <= 1 or len(trial_specs) <= 1:
        return [_evaluate_trial(spec) for spec in trial_specs]

    ordered_records: list[dict[str, Any] | None] = [None] * len(trial_specs)
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="ml-grid") as pool:
        future_to_trial_id = {
            pool.submit(_evaluate_trial, spec): spec[0]
            for spec in trial_specs
        }
        for future in as_completed(future_to_trial_id):
            trial_id_local = future_to_trial_id[future]
            ordered_records[trial_id_local] = future.result()
    return [record for record in ordered_records if record is not None]


def _run_optuna_search(
    *,
    feature_variant_names: list[str],
    feature_variants: dict[str, FeaturesConfig],
    model_families: list[str],
    calibration_options: list[tuple[bool, str]],
    feature_cache: dict[str, FeatureCacheEntry],
    random_seed: int,
    inner_scaffold_folds: int,
    calibration_cv_folds: int,
    calibration_seed: int,
    max_trials: int,
    parallel_workers: int = 1,
) -> list[dict[str, Any]]:
    try:
        import optuna  # type: ignore[import-untyped]
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "Optuna mode requested but optuna is not installed. "
            "Install with: pip install optuna"
        ) from exc

    trial_records: list[dict[str, Any]] = []
    trial_records_lock = threading.Lock()
    best_inner_ef1 = {"value": float("-inf")}
    best_inner_lock = threading.Lock()
    calibration_choice_map: dict[str, tuple[bool, str]] = {}
    calibration_choice_labels: list[str] = []
    for enabled, method in calibration_options:
        label = f"{'enabled' if enabled else 'disabled'}::{str(method).strip().lower()}"
        if label in calibration_choice_map:
            continue
        calibration_choice_map[label] = (bool(enabled), str(method).strip().lower())
        calibration_choice_labels.append(label)
    if not calibration_choice_labels:
        calibration_choice_labels = ["disabled::sigmoid"]
        calibration_choice_map = {"disabled::sigmoid": (False, "sigmoid")}

    def objective(trial) -> float:  # type: ignore[no-untyped-def]
        feature_variant_name = trial.suggest_categorical(
            "feature_variant",
            feature_variant_names,
        )
        model_family = trial.suggest_categorical("model_family", model_families)
        calibration_choice = trial.suggest_categorical(
            "calibration_option",
            calibration_choice_labels,
        )
        calibration_enabled, calibration_method = calibration_choice_map[calibration_choice]

        if model_family == "logreg":
            params = {
                "C": float(trial.suggest_float("logreg_C", 0.01, 10.0, log=True)),
                "penalty": "l2",
                "solver": "liblinear",
                "class_weight": "balanced",
            }
        elif model_family == "lightgbm":
            params = {
                "max_depth": int(trial.suggest_categorical("lgbm_max_depth", [-1, 6, 10])),
                "learning_rate": float(
                    trial.suggest_float("lgbm_learning_rate", 0.01, 0.2, log=True)
                ),
                "num_leaves": int(
                    trial.suggest_int("lgbm_num_leaves", 31, 127, step=16)
                ),
                "min_child_samples": int(
                    trial.suggest_int("lgbm_min_child_samples", 10, 80, step=10)
                ),
                "n_estimators": int(trial.suggest_int("lgbm_n_estimators", 100, 400, step=100)),
                "class_weight": "balanced",
            }
        elif model_family == "xgboost":
            params = {
                "max_depth": int(trial.suggest_int("xgb_max_depth", 3, 10)),
                "eta": float(trial.suggest_float("xgb_eta", 0.01, 0.2, log=True)),
                "n_estimators": int(trial.suggest_int("xgb_n_estimators", 100, 400, step=100)),
                "subsample": float(trial.suggest_float("xgb_subsample", 0.6, 1.0)),
                "colsample_bytree": float(trial.suggest_float("xgb_colsample", 0.6, 1.0)),
            }
        elif model_family == "hgb":  # legacy option
            params = {
                "max_depth": trial.suggest_categorical("hgb_max_depth", [3, 6, None]),
                "learning_rate": float(
                    trial.suggest_float("hgb_learning_rate", 0.01, 0.2, log=True)
                ),
                "max_leaf_nodes": int(
                    trial.suggest_int("hgb_max_leaf_nodes", 31, 127, step=16)
                ),
                "min_samples_leaf": int(
                    trial.suggest_int("hgb_min_samples_leaf", 10, 80, step=10)
                ),
            }
        elif model_family == "lightgbm_rank":
            params = {
                "objective": "lambdarank",
                "learning_rate": float(
                    trial.suggest_float("lgbmr_learning_rate", 0.01, 0.2, log=True)
                ),
                "num_leaves": int(trial.suggest_int("lgbmr_num_leaves", 31, 127, step=16)),
                "min_child_samples": int(
                    trial.suggest_int("lgbmr_min_child_samples", 10, 80, step=10)
                ),
                "n_estimators": int(trial.suggest_int("lgbmr_n_estimators", 100, 400, step=100)),
            }
        elif model_family == "xgboost_rank":
            params = {
                "objective": "rank:pairwise",
                "max_depth": int(trial.suggest_int("xgbr_max_depth", 3, 10)),
                "eta": float(trial.suggest_float("xgbr_eta", 0.01, 0.2, log=True)),
                "n_estimators": int(trial.suggest_int("xgbr_n_estimators", 100, 400, step=100)),
                "subsample": float(trial.suggest_float("xgbr_subsample", 0.6, 1.0)),
                "colsample_bytree": float(trial.suggest_float("xgbr_colsample", 0.6, 1.0)),
            }
        else:
            raise ValueError(f"Unsupported model_family={model_family!r}")

        record: dict[str, Any] = {
            "trial_id": int(trial.number),
            "mode": "optuna",
            "status": "ok",
            "feature_variant": feature_variant_name,
            "features": asdict(feature_variants[feature_variant_name]),
            "model_family": model_family,
            "model_params": dict(params),
            "calibration_enabled": bool(calibration_enabled),
            "calibration_method": str(calibration_method),
            "inner_metrics": {},
            "error": "",
        }
        with trial_records_lock:
            trial_records.append(record)

        train_rows = feature_cache[feature_variant_name][1]
        train_query_groups = feature_cache[feature_variant_name][5]
        train_sample_weight = feature_cache[feature_variant_name][6]
        splits = _group_splits(
            train_rows.murcko_scaffolds,
            train_rows.y,
            inner_scaffold_folds,
        )
        fold_reports: list[dict[str, float | int]] = []

        for fold_idx, (train_idx, val_idx) in enumerate(splits):
            y_train_fold = train_rows.y[train_idx]
            y_val_fold = train_rows.y[val_idx]
            if np.unique(y_train_fold).size < 2 or np.unique(y_val_fold).size < 2:
                continue
            model = _fit_model(
                model_family=model_family,
                model_params=params,
                X_train=train_rows.X[train_idx],
                y_train=y_train_fold,
                query_groups=(
                    train_query_groups[train_idx] if train_query_groups is not None else None
                ),
                sample_weight=(
                    train_sample_weight[train_idx] if train_sample_weight is not None else None
                ),
                random_seed=random_seed + fold_idx,
            )
            p_active = _predict_p_active(model, train_rows.X[val_idx])
            report = evaluate_holdout_metrics(
                y_val_fold,
                apply_calibration(
                    fit_oof_calibrator(
                        X=train_rows.X[train_idx],
                        y=y_train_fold,
                        groups=np.asarray(
                            [train_rows.murcko_scaffolds[idx] for idx in train_idx.tolist()],
                            dtype=object,
                        ),
                        base_model_builder=lambda fold_seed: build_model(
                            model_family=model_family,
                            model_params=params,
                            random_seed=fold_seed,
                        ),
                        method=str(calibration_method),
                        cv_folds=calibration_cv_folds,
                        seed=calibration_seed + fold_idx,
                        fit_query_groups=(
                            train_query_groups[train_idx] if train_query_groups is not None else None
                        ),
                        fit_sample_weights=(
                            train_sample_weight[train_idx] if train_sample_weight is not None else None
                        ),
                    )[0],
                    p_active,
                )
                if calibration_enabled
                else p_active,
                fractions=DEFAULT_FRACTIONS,
            )
            fold_reports.append(report)

            running_ef1 = float(
                np.nanmean(
                    np.asarray(
                        [
                            float(r.get("EF@1%", float("nan")))
                            for r in fold_reports
                            if _is_finite_number(r.get("EF@1%"))
                        ],
                        dtype=float,
                    )
                )
            )
            if _is_finite_number(running_ef1):
                trial.report(running_ef1, step=fold_idx)
                if (
                    fold_idx >= 1
                    and _is_finite_number(best_inner_ef1["value"])
                    and running_ef1 < (0.5 * float(best_inner_ef1["value"]))
                ):
                    record["status"] = "pruned"
                    record["error"] = "manual_prune_low_ef1"
                    raise optuna.TrialPruned("manual_prune_low_ef1")
                if trial.should_prune():
                    record["status"] = "pruned"
                    record["error"] = "optuna_pruner"
                    raise optuna.TrialPruned("optuna_pruner")

        inner_metrics = _aggregate_fold_reports(fold_reports)
        record["inner_metrics"] = inner_metrics
        if int(inner_metrics.get("inner_folds", 0)) <= 0:
            record["status"] = "skipped"
            record["error"] = "no_valid_inner_folds"
            return float("-inf")

        ef1 = _metric_for_sort(inner_metrics.get("inner_mean_EF@1%"))
        with best_inner_lock:
            if ef1 > best_inner_ef1["value"]:
                best_inner_ef1["value"] = ef1

        pr_auc_mean = _metric_for_sort(inner_metrics.get("inner_mean_PR_AUC"))
        return ef1 + (1e-6 * pr_auc_mean)

    sampler = optuna.samplers.TPESampler(seed=int(random_seed))
    pruner = optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=1)
    study = optuna.create_study(direction="maximize", sampler=sampler, pruner=pruner)
    study.optimize(
        objective,
        n_trials=int(max_trials),
        n_jobs=max(1, int(parallel_workers)),
    )
    return trial_records


def run_automl(
    config_path: str,
    mode: str = "grid",
    *,
    max_trials: int = 50,
    top_k_report: int = 5,
) -> Path:
    logger = _setup_logging()
    config = load_config(config_path)
    raw_automl_payload = _load_automl_payload(config_path)
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
    feature_variant_names = [
        name
        for name in _feature_variant_names_from_payload(raw_automl_payload)
        if name in feature_variants
    ]
    model_families = _model_families_from_payload(raw_automl_payload)
    requires_rank_groups = any("rank" in str(name) for name in model_families)
    calibration_options = _calibration_options_from_payload(
        raw_automl_payload,
        default_enabled=bool(config.calibration_enabled),
        default_method=str(config.calibration_method),
    )
    trial_workers = _resolve_parallel_workers(
        automl_payload=raw_automl_payload,
        atlas_cfg=atlas_cfg,
    )
    fpocket_workers = _resolve_fpocket_workers(atlas_cfg)
    run_id = _resolve_run_id(config.run_id)
    run_dir = Path(__file__).resolve().parent / "outputs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    logger.info(
        "[ml.automl.parallel] trial_workers=%d fpocket_workers=%d",
        trial_workers,
        fpocket_workers,
    )

    train_df = dataset.train_df.copy()
    val_df = dataset.val_df.copy()
    test_df = dataset.test_df.copy()
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
    needs_fpocket_precompute = any(
        (feature_variants[name].pocket_fpocket or feature_variants[name].pocket_features)
        for name in feature_variant_names
    )
    if needs_fpocket_precompute:
        combined = pd.concat([train_df, val_df, test_df], ignore_index=False)
        fpocket_metrics_df = precompute_fpocket_for_bigbind_df(
            df=combined,
            bigbind_root=dataset.bigbind_root,
            atlas_cfg=atlas_cfg,
            run_dir=run_dir,
            max_workers=fpocket_workers,
        )
        train_df = merge_fpocket_metrics_on_pocket(train_df, fpocket_metrics_df)
        val_df = merge_fpocket_metrics_on_pocket(val_df, fpocket_metrics_df)
        test_df = merge_fpocket_metrics_on_pocket(test_df, fpocket_metrics_df)

    grid_spaces = _default_grid_spaces()
    grid_override = raw_automl_payload.get("grid")
    if isinstance(grid_override, dict):
        for model_family, space in grid_override.items():
            family = str(model_family).strip().lower()
            if family == "lgbm":
                family = "lightgbm"
            if family == "xgb":
                family = "xgboost"
            if family == "lgbm_rank":
                family = "lightgbm_rank"
            if family == "xgb_rank":
                family = "xgboost_rank"
            if family in grid_spaces:
                grid_spaces[family] = _normalize_grid_space(
                    space,
                    fallback=grid_spaces[family],
                )

    feature_cache: dict[str, FeatureCacheEntry] = {}
    for feature_variant_name in feature_variant_names:
        features = feature_variants[feature_variant_name]
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
        train_rows = featurizer.transform(train_df)
        val_rows = featurizer.transform(val_df)
        test_rows = featurizer.transform(test_df)
        train_used = train_df.loc[train_rows.source_index].copy()
        train_query_groups = (
            _compute_query_groups(train_used, str(config.rank.group_key))
            if (config.rank.enabled or requires_rank_groups)
            else None
        )
        train_sample_weight = (
            derive_sample_weight(train_used, weight_cap=config.labels.weight_cap)
            if config.labels.use_sample_weights
            else None
        )
        feature_cache[feature_variant_name] = (
            features,
            train_rows,
            val_rows,
            test_rows,
            {
                "fpocket_loaded_count": int(featurizer.loaded_fpocket_count),
                "fpocket_missing_center_count": int(featurizer.missing_center_count),
                "fpocket_missing_info_count": int(featurizer.missing_info_count),
            },
            train_query_groups,
            train_sample_weight,
        )

    mode_norm = str(mode).strip().lower()
    if mode_norm not in {"grid", "optuna"}:
        raise ValueError("mode must be 'grid' or 'optuna'.")

    if mode_norm == "grid":
        trial_records = _run_grid_search(
            feature_variant_names=feature_variant_names,
            feature_variants=feature_variants,
            model_families=model_families,
            calibration_options=calibration_options,
            grid_spaces=grid_spaces,
            feature_cache=feature_cache,
            random_seed=config.random_seed,
            inner_scaffold_folds=config.inner_scaffold_folds,
            calibration_cv_folds=config.calibration_cv_folds,
            calibration_seed=config.calibration_seed,
            parallel_workers=trial_workers,
        )
    else:
        trial_records = _run_optuna_search(
            feature_variant_names=feature_variant_names,
            feature_variants=feature_variants,
            model_families=model_families,
            calibration_options=calibration_options,
            feature_cache=feature_cache,
            random_seed=config.random_seed,
            inner_scaffold_folds=config.inner_scaffold_folds,
            calibration_cv_folds=config.calibration_cv_folds,
            calibration_seed=config.calibration_seed,
            max_trials=max_trials,
            parallel_workers=trial_workers,
        )

    trial_rows = [_materialize_trial_row(record) for record in trial_records]
    trials_df = pd.DataFrame(trial_rows)
    trials_df.to_csv(run_dir / "automl_trials.csv", index=False)

    successful_trials = [r for r in trial_records if r.get("status") == "ok"]
    successful_trials = [
        r
        for r in successful_trials
        if int(r.get("inner_metrics", {}).get("inner_folds", 0)) > 0
    ]
    if not successful_trials:
        raise RuntimeError("No successful AutoML trials with valid inner folds.")

    successful_trials.sort(key=_trial_sort_key, reverse=True)
    best_trial = successful_trials[0]

    best_variant_name = str(best_trial["feature_variant"])
    best_train_rows = feature_cache[best_variant_name][1]
    best_val_rows = feature_cache[best_variant_name][2]
    best_test_rows = feature_cache[best_variant_name][3]
    best_train_query_groups = feature_cache[best_variant_name][5]
    best_train_sample_weight = feature_cache[best_variant_name][6]

    best_model = _fit_model(
        model_family=str(best_trial["model_family"]),
        model_params=dict(best_trial["model_params"]),
        X_train=best_train_rows.X,
        y_train=best_train_rows.y,
        query_groups=best_train_query_groups,
        sample_weight=best_train_sample_weight,
        random_seed=config.random_seed,
    )
    joblib.dump(best_model, run_dir / "best_model.joblib")

    best_val_metrics, best_val_pred_df = _evaluate_holdout_trial(
        trial_record=best_trial,
        train_rows=best_train_rows,
        holdout_rows=best_val_rows,
        holdout_df=val_df,
        random_seed=config.random_seed,
        calibration_cv_folds=config.calibration_cv_folds,
        calibration_seed=config.calibration_seed,
        train_query_groups=best_train_query_groups,
        sample_weight=best_train_sample_weight,
    )
    best_test_metrics, best_pred_df = _evaluate_holdout_trial(
        trial_record=best_trial,
        train_rows=best_train_rows,
        holdout_rows=best_test_rows,
        holdout_df=test_df,
        random_seed=config.random_seed,
        calibration_cv_folds=config.calibration_cv_folds,
        calibration_seed=config.calibration_seed,
        train_query_groups=best_train_query_groups,
        sample_weight=best_train_sample_weight,
    )
    best_val_pred_df.to_csv(run_dir / "best_validation_predictions.csv", index=False)
    pd.DataFrame([best_val_metrics]).to_csv(
        run_dir / "best_validation_report.csv",
        index=False,
    )
    (run_dir / "best_validation_report.json").write_text(
        json.dumps(best_val_metrics, indent=2),
        encoding="utf-8",
    )
    best_pred_df.to_csv(run_dir / "best_holdout_predictions.csv", index=False)
    best_pred_df.to_csv(run_dir / "best_test_predictions.csv", index=False)
    pd.DataFrame([best_test_metrics]).to_csv(
        run_dir / "best_holdout_report.csv",
        index=False,
    )
    pd.DataFrame([best_test_metrics]).to_csv(
        run_dir / "best_test_report.csv",
        index=False,
    )
    (run_dir / "best_holdout_report.json").write_text(
        json.dumps(best_test_metrics, indent=2),
        encoding="utf-8",
    )
    (run_dir / "best_test_report.json").write_text(
        json.dumps(best_test_metrics, indent=2),
        encoding="utf-8",
    )

    best_trial_payload = {
        "trial_id": int(best_trial["trial_id"]),
        "mode": str(best_trial["mode"]),
        "feature_variant": str(best_trial["feature_variant"]),
        "features": dict(best_trial["features"]),
        "model_family": str(best_trial["model_family"]),
        "calibration_enabled": bool(best_trial.get("calibration_enabled", False)),
        "calibration_method": str(best_trial.get("calibration_method", "sigmoid")),
        "model_params": dict(best_trial["model_params"]),
        "inner_metrics": dict(best_trial["inner_metrics"]),
        "selection_metrics": dict(best_val_metrics),
        "holdout_metrics": dict(best_test_metrics),
        "fpocket_stats": dict(feature_cache[best_variant_name][4]),
        "search_space": {
            "feature_variants": feature_variant_names,
            "model_families": model_families,
        },
        "selection_split": "val",
        "final_eval_split": "test",
    }
    (run_dir / "best_trial.json").write_text(
        json.dumps(best_trial_payload, indent=2),
        encoding="utf-8",
    )

    top_k = max(1, min(int(top_k_report), 5))
    top_trials = successful_trials[:top_k]
    top_rows: list[dict[str, Any]] = []
    for rank, trial_record in enumerate(top_trials, start=1):
        variant_name = str(trial_record["feature_variant"])
        train_rows = feature_cache[variant_name][1]
        holdout_rows = feature_cache[variant_name][2]
        train_query_groups = feature_cache[variant_name][5]
        train_sample_weight = feature_cache[variant_name][6]
        holdout_metrics, _ = _evaluate_holdout_trial(
            trial_record=trial_record,
            train_rows=train_rows,
            holdout_rows=holdout_rows,
            holdout_df=val_df,
            random_seed=config.random_seed,
            calibration_cv_folds=config.calibration_cv_folds,
            calibration_seed=config.calibration_seed,
            train_query_groups=train_query_groups,
            sample_weight=train_sample_weight,
        )
        top_row: dict[str, Any] = {
            "rank_by_inner": rank,
            "trial_id": int(trial_record["trial_id"]),
            "feature_variant": str(trial_record["feature_variant"]),
            "model_family": str(trial_record["model_family"]),
            "calibration_enabled": int(bool(trial_record.get("calibration_enabled", False))),
            "calibration_method": str(trial_record.get("calibration_method", "sigmoid")),
            "inner_mean_EF@1%": trial_record["inner_metrics"].get("inner_mean_EF@1%"),
            "inner_mean_PR_AUC": trial_record["inner_metrics"].get("inner_mean_PR_AUC"),
            "model_params_json": json.dumps(trial_record["model_params"], sort_keys=True),
        }
        top_row.update({f"holdout_{k}": v for k, v in holdout_metrics.items()})
        top_rows.append(top_row)
    pd.DataFrame(top_rows).to_csv(run_dir / "top_k_holdout_report.csv", index=False)
    pd.DataFrame(top_rows).to_csv(run_dir / "top_k_validation_report.csv", index=False)

    (run_dir / "config_snapshot.txt").write_text(
        json.dumps(config_to_dict(config), indent=2),
        encoding="utf-8",
    )
    (run_dir / "automl_snapshot.json").write_text(
        json.dumps(
            {
                "mode": mode_norm,
                "max_trials": int(max_trials),
                "top_k_report": int(top_k),
                "top_k_report_split": "val",
                "feature_variants": feature_variant_names,
                "model_families": model_families,
                "calibration_options": calibration_options,
                "raw_payload": raw_automl_payload,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (run_dir / "config_snapshot.json").write_text(
        json.dumps(config_to_dict(config), indent=2),
        encoding="utf-8",
    )
    write_truncated_views_for_run_dir(run_dir)

    if config.registry_enabled:
        repo_root = Path(__file__).resolve().parents[1]
        feature_columns = [f"f{i}" for i in range(int(best_train_rows.X.shape[1]))]
        dataset_hash = compute_dataset_hash(
            train_df=train_df.loc[best_train_rows.source_index].copy(),
            holdout_df=test_df.loc[best_test_rows.source_index].copy(),
            feature_columns=feature_columns,
            config_snapshot=config_to_dict(config),
        )
        registry_record = {
            "run_id": run_id,
            "run_dir": str(run_dir),
            "dataset_hash": dataset_hash,
            "git_sha": get_git_sha(repo_root),
            "best_trial_id": int(best_trial["trial_id"]),
            "best_feature_variant": str(best_trial["feature_variant"]),
            "best_model_family": str(best_trial["model_family"]),
            "best_model_params": dict(best_trial["model_params"]),
            "calibration_enabled": bool(best_trial.get("calibration_enabled", False)),
            "calibration_method": str(best_trial.get("calibration_method", "sigmoid")),
            "rank_enabled": bool(config.rank.enabled),
            "rank_group_key": str(config.rank.group_key),
            "hard_negatives": config_to_dict(config).get("hard_negatives", {}),
            "selection_metrics": dict(best_val_metrics),
            "holdout_metrics": dict(best_test_metrics),
            "fpocket_cache_keys_used": _collect_fpocket_cache_keys_from_df(fpocket_metrics_df),
            "environment": collect_env_versions(),
        }
        write_registry_record(run_dir, registry_record)
        append_registry_index(
            run_dir.parent / "registry_index.csv",
            {
                "run_id": run_id,
                "run_dir": str(run_dir),
                "dataset_hash": dataset_hash,
                "git_sha": registry_record["git_sha"],
                "model_family": registry_record["best_model_family"],
                "calibration_enabled": int(registry_record["calibration_enabled"]),
                "calibration_method": registry_record["calibration_method"],
                "holdout_PR_AUC": float(best_test_metrics.get("PR_AUC", float("nan"))),
                "holdout_EF@1%": float(best_test_metrics.get("EF@1%", float("nan"))),
            },
        )

    inner_metrics_payload = cast(dict[str, Any], best_trial_payload["inner_metrics"])
    logger.info(
        "[ml.automl] saved=%s best_trial=%s best_inner_EF@1%%=%s best_inner_PR_AUC=%s",
        run_dir,
        best_trial_payload["trial_id"],
        inner_metrics_payload.get("inner_mean_EF@1%"),
        inner_metrics_payload.get("inner_mean_PR_AUC"),
    )
    return run_dir
