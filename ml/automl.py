#!/usr/bin/env python3
from __future__ import annotations

import itertools
import json
import logging
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any, Iterable, cast

import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from ml.config import FeaturesConfig, config_to_dict, load_atlas_cfg, load_config
from ml.data.bigbind import load_train_holdout_from_bigbind
from ml.evaluate import DEFAULT_FRACTIONS, evaluate_holdout_metrics
from ml.featurize import BigBindFeaturizer, FeaturizedRows

try:  # pragma: no cover - optional dependency
    from lightgbm import LGBMClassifier  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    LGBMClassifier = None  # type: ignore[assignment]


_FEATURE_VARIANT_ORDER = (
    "baseline",
    "ligand_only",
    "pocket_only",
    "no_fp",
    "add_ligand_extras",
    "add_pocket_fpocket",
    "add_both",
)
_DEFAULT_MODEL_FAMILIES = ("logreg", "hgb")


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


def _load_json_or_yaml(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix in {".json", ".txt"}:
        with path.open("r", encoding="utf-8") as handle:
            loaded = json.load(handle)
        return loaded if isinstance(loaded, dict) else {}
    if suffix in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore[import-untyped]
        except Exception as exc:  # pragma: no cover - optional dependency
            raise RuntimeError(
                "YAML config requested but PyYAML is not installed. "
                "Use JSON config or install PyYAML."
            ) from exc
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
        return loaded if isinstance(loaded, dict) else {}
    raise ValueError(f"Unsupported config extension: {path.suffix}")


def _load_automl_payload(config_path: str | Path) -> dict[str, Any]:
    payload = _load_json_or_yaml(Path(config_path).expanduser().resolve())
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
            "solver": ["lbfgs"],
            "class_weight": ["balanced"],
        },
        "hgb": {
            "max_depth": [3, 6, None],
            "learning_rate": [0.03, 0.1],
            "max_leaf_nodes": [31, 63],
            "min_samples_leaf": [20, 50],
        },
    }
    if LGBMClassifier is not None:
        spaces["lgbm"] = {
            "learning_rate": [0.03, 0.1],
            "num_leaves": [31, 63],
            "min_child_samples": [20, 50],
            "n_estimators": [200],
            "class_weight": ["balanced"],
        }
    return spaces


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
            vina_score=False,
        ),
        "no_fp": replace(baseline, ligand_morgan_fp_bits=0),
        "add_ligand_extras": replace(baseline, ligand_extra_descriptors=True),
        "add_pocket_fpocket": replace(baseline, pocket_fpocket=True),
        "add_both": replace(baseline, ligand_extra_descriptors=True, pocket_fpocket=True),
    }


def _model_families_from_payload(automl_payload: dict[str, Any]) -> list[str]:
    default = list(_DEFAULT_MODEL_FAMILIES)
    if LGBMClassifier is not None:
        default.append("lgbm")
    names = _normalize_name_list(automl_payload.get("model_families"), default=default)
    valid = set(default)
    selected = [name for name in names if name in valid]
    return selected or default


def _feature_variant_names_from_payload(automl_payload: dict[str, Any]) -> list[str]:
    names = _normalize_name_list(
        automl_payload.get("feature_variants"),
        default=_FEATURE_VARIANT_ORDER,
    )
    valid = set(_FEATURE_VARIANT_ORDER)
    selected = [name for name in names if name in valid]
    return selected or list(_FEATURE_VARIANT_ORDER)


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
    random_seed: int,
) -> Any:
    family = model_family.strip().lower()
    if family == "logreg":
        model = LogisticRegression(
            C=float(model_params.get("C", 1.0)),
            penalty=str(model_params.get("penalty", "l2")),
            solver=str(model_params.get("solver", "lbfgs")),
            class_weight=model_params.get("class_weight", "balanced"),
            max_iter=int(model_params.get("max_iter", 1000)),
            random_state=int(random_seed),
        )
        model.fit(X_train, y_train)
        return model

    if family == "hgb":
        model = HistGradientBoostingClassifier(
            max_depth=model_params.get("max_depth"),
            learning_rate=float(model_params.get("learning_rate", 0.1)),
            max_leaf_nodes=int(model_params.get("max_leaf_nodes", 31)),
            min_samples_leaf=int(model_params.get("min_samples_leaf", 20)),
            random_state=int(random_seed),
        )
        x_dense = X_train.toarray()
        model.fit(x_dense, y_train)
        setattr(model, "_atlas_requires_dense", True)
        return model

    if family == "lgbm":
        if LGBMClassifier is None:
            raise RuntimeError("model_family='lgbm' requested but lightgbm is not installed.")
        model = LGBMClassifier(
            learning_rate=float(model_params.get("learning_rate", 0.1)),
            num_leaves=int(model_params.get("num_leaves", 31)),
            min_child_samples=int(model_params.get("min_child_samples", 20)),
            n_estimators=int(model_params.get("n_estimators", 200)),
            class_weight=model_params.get("class_weight", "balanced"),
            random_state=int(random_seed),
            n_jobs=1,
            verbosity=-1,
        )
        model.fit(X_train, y_train)
        return model

    raise ValueError(f"Unsupported model_family={model_family!r}")


def _predict_p_active(model: Any, X: sparse.csr_matrix) -> np.ndarray:
    if getattr(model, "_atlas_requires_dense", False):
        x_used = X.toarray()
    else:
        x_used = X
    prob = model.predict_proba(x_used)
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
    random_seed: int,
    n_splits_requested: int,
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
            random_seed=random_seed + fold_idx,
        )
        p_active = _predict_p_active(model, X[val_idx])
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


def _evaluate_holdout_trial(
    *,
    trial_record: dict[str, Any],
    train_rows: FeaturizedRows,
    holdout_rows: FeaturizedRows,
    holdout_df: pd.DataFrame,
    random_seed: int,
) -> tuple[dict[str, float | int], pd.DataFrame]:
    model = _fit_model(
        model_family=str(trial_record["model_family"]),
        model_params=dict(trial_record["model_params"]),
        X_train=train_rows.X,
        y_train=train_rows.y,
        random_seed=random_seed,
    )
    holdout_p_active = _predict_p_active(model, holdout_rows.X)
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
    grid_spaces: dict[str, dict[str, list[Any]]],
    feature_cache: dict[str, tuple[FeaturesConfig, FeaturizedRows, FeaturizedRows, dict[str, int]]],
    random_seed: int,
    inner_scaffold_folds: int,
) -> list[dict[str, Any]]:
    trial_records: list[dict[str, Any]] = []
    trial_id = 0
    for feature_variant_name in feature_variant_names:
        train_rows = feature_cache[feature_variant_name][1]
        for model_family in model_families:
            if model_family not in grid_spaces:
                continue
            for params in _iter_grid_param_sets(grid_spaces[model_family]):
                record: dict[str, Any] = {
                    "trial_id": trial_id,
                    "mode": "grid",
                    "status": "ok",
                    "feature_variant": feature_variant_name,
                    "features": asdict(feature_variants[feature_variant_name]),
                    "model_family": model_family,
                    "model_params": dict(params),
                    "inner_metrics": {},
                    "error": "",
                }
                trial_id += 1
                try:
                    inner_metrics = _inner_cv_metrics(
                        X=train_rows.X,
                        y=train_rows.y,
                        murcko_scaffolds=train_rows.murcko_scaffolds,
                        model_family=model_family,
                        model_params=params,
                        random_seed=random_seed,
                        n_splits_requested=inner_scaffold_folds,
                    )
                    if int(inner_metrics.get("inner_folds", 0)) <= 0:
                        record["status"] = "skipped"
                        record["error"] = "no_valid_inner_folds"
                    record["inner_metrics"] = inner_metrics
                except Exception as exc:
                    record["status"] = "failed"
                    record["error"] = f"{type(exc).__name__}: {exc}"
                trial_records.append(record)
    return trial_records


def _run_optuna_search(
    *,
    feature_variant_names: list[str],
    feature_variants: dict[str, FeaturesConfig],
    model_families: list[str],
    feature_cache: dict[str, tuple[FeaturesConfig, FeaturizedRows, FeaturizedRows, dict[str, int]]],
    random_seed: int,
    inner_scaffold_folds: int,
    max_trials: int,
) -> list[dict[str, Any]]:
    try:
        import optuna  # type: ignore[import-untyped]
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "Optuna mode requested but optuna is not installed. "
            "Install with: pip install optuna"
        ) from exc

    trial_records: list[dict[str, Any]] = []
    best_inner_ef1 = {"value": float("-inf")}

    def objective(trial) -> float:  # type: ignore[no-untyped-def]
        feature_variant_name = trial.suggest_categorical(
            "feature_variant",
            feature_variant_names,
        )
        model_family = trial.suggest_categorical("model_family", model_families)

        if model_family == "logreg":
            params = {
                "C": float(trial.suggest_float("logreg_C", 0.01, 10.0, log=True)),
                "penalty": "l2",
                "solver": "lbfgs",
                "class_weight": "balanced",
            }
        elif model_family == "hgb":
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
        elif model_family == "lgbm":
            params = {
                "learning_rate": float(
                    trial.suggest_float("lgbm_learning_rate", 0.01, 0.2, log=True)
                ),
                "num_leaves": int(
                    trial.suggest_int("lgbm_num_leaves", 31, 127, step=16)
                ),
                "min_child_samples": int(
                    trial.suggest_int("lgbm_min_child_samples", 10, 80, step=10)
                ),
                "n_estimators": int(
                    trial.suggest_int("lgbm_n_estimators", 100, 400, step=100)
                ),
                "class_weight": "balanced",
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
            "inner_metrics": {},
            "error": "",
        }
        trial_records.append(record)

        train_rows = feature_cache[feature_variant_name][1]
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
                random_seed=random_seed + fold_idx,
            )
            p_active = _predict_p_active(model, train_rows.X[val_idx])
            report = evaluate_holdout_metrics(
                y_val_fold,
                p_active,
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
        if ef1 > best_inner_ef1["value"]:
            best_inner_ef1["value"] = ef1

        pr_auc_mean = _metric_for_sort(inner_metrics.get("inner_mean_PR_AUC"))
        return ef1 + (1e-6 * pr_auc_mean)

    sampler = optuna.samplers.TPESampler(seed=int(random_seed))
    pruner = optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=1)
    study = optuna.create_study(direction="maximize", sampler=sampler, pruner=pruner)
    study.optimize(objective, n_trials=int(max_trials))
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
        logger=logger,
    )

    feature_variants = _build_feature_variants(config.features)
    feature_variant_names = [
        name
        for name in _feature_variant_names_from_payload(raw_automl_payload)
        if name in feature_variants
    ]
    model_families = _model_families_from_payload(raw_automl_payload)

    grid_spaces = _default_grid_spaces()
    grid_override = raw_automl_payload.get("grid")
    if isinstance(grid_override, dict):
        for model_family, space in grid_override.items():
            family = str(model_family).strip().lower()
            if family in grid_spaces:
                grid_spaces[family] = _normalize_grid_space(
                    space,
                    fallback=grid_spaces[family],
                )

    feature_cache: dict[
        str, tuple[FeaturesConfig, FeaturizedRows, FeaturizedRows, dict[str, int]]
    ] = {}
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
        train_rows = featurizer.transform(dataset.train_df)
        holdout_rows = featurizer.transform(dataset.holdout_df)
        feature_cache[feature_variant_name] = (
            features,
            train_rows,
            holdout_rows,
            {
                "fpocket_loaded_count": int(featurizer.loaded_fpocket_count),
                "fpocket_missing_center_count": int(featurizer.missing_center_count),
                "fpocket_missing_info_count": int(featurizer.missing_info_count),
            },
        )

    mode_norm = str(mode).strip().lower()
    if mode_norm not in {"grid", "optuna"}:
        raise ValueError("mode must be 'grid' or 'optuna'.")

    if mode_norm == "grid":
        trial_records = _run_grid_search(
            feature_variant_names=feature_variant_names,
            feature_variants=feature_variants,
            model_families=model_families,
            grid_spaces=grid_spaces,
            feature_cache=feature_cache,
            random_seed=config.random_seed,
            inner_scaffold_folds=config.inner_scaffold_folds,
        )
    else:
        trial_records = _run_optuna_search(
            feature_variant_names=feature_variant_names,
            feature_variants=feature_variants,
            model_families=model_families,
            feature_cache=feature_cache,
            random_seed=config.random_seed,
            inner_scaffold_folds=config.inner_scaffold_folds,
            max_trials=max_trials,
        )

    run_id = _resolve_run_id(config.run_id)
    run_dir = Path(__file__).resolve().parent / "outputs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

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
    best_holdout_rows = feature_cache[best_variant_name][2]

    best_model = _fit_model(
        model_family=str(best_trial["model_family"]),
        model_params=dict(best_trial["model_params"]),
        X_train=best_train_rows.X,
        y_train=best_train_rows.y,
        random_seed=config.random_seed,
    )
    joblib.dump(best_model, run_dir / "best_model.joblib")

    best_holdout_metrics, best_pred_df = _evaluate_holdout_trial(
        trial_record=best_trial,
        train_rows=best_train_rows,
        holdout_rows=best_holdout_rows,
        holdout_df=dataset.holdout_df,
        random_seed=config.random_seed,
    )
    best_pred_df.to_csv(run_dir / "best_holdout_predictions.csv", index=False)
    pd.DataFrame([best_holdout_metrics]).to_csv(
        run_dir / "best_holdout_report.csv",
        index=False,
    )
    (run_dir / "best_holdout_report.json").write_text(
        json.dumps(best_holdout_metrics, indent=2),
        encoding="utf-8",
    )

    best_trial_payload = {
        "trial_id": int(best_trial["trial_id"]),
        "mode": str(best_trial["mode"]),
        "feature_variant": str(best_trial["feature_variant"]),
        "features": dict(best_trial["features"]),
        "model_family": str(best_trial["model_family"]),
        "model_params": dict(best_trial["model_params"]),
        "inner_metrics": dict(best_trial["inner_metrics"]),
        "holdout_metrics": dict(best_holdout_metrics),
        "fpocket_stats": dict(feature_cache[best_variant_name][3]),
        "search_space": {
            "feature_variants": feature_variant_names,
            "model_families": model_families,
        },
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
        holdout_metrics, _ = _evaluate_holdout_trial(
            trial_record=trial_record,
            train_rows=train_rows,
            holdout_rows=holdout_rows,
            holdout_df=dataset.holdout_df,
            random_seed=config.random_seed,
        )
        top_row: dict[str, Any] = {
            "rank_by_inner": rank,
            "trial_id": int(trial_record["trial_id"]),
            "feature_variant": str(trial_record["feature_variant"]),
            "model_family": str(trial_record["model_family"]),
            "inner_mean_EF@1%": trial_record["inner_metrics"].get("inner_mean_EF@1%"),
            "inner_mean_PR_AUC": trial_record["inner_metrics"].get("inner_mean_PR_AUC"),
            "model_params_json": json.dumps(trial_record["model_params"], sort_keys=True),
        }
        top_row.update({f"holdout_{k}": v for k, v in holdout_metrics.items()})
        top_rows.append(top_row)
    pd.DataFrame(top_rows).to_csv(run_dir / "top_k_holdout_report.csv", index=False)

    snapshot = config_to_dict(config)
    snapshot["automl"] = {
        "mode": mode_norm,
        "max_trials": int(max_trials),
        "top_k_report": int(top_k),
        "feature_variants": feature_variant_names,
        "model_families": model_families,
        "raw_payload": raw_automl_payload,
    }
    (run_dir / "config_snapshot.json").write_text(
        json.dumps(snapshot, indent=2),
        encoding="utf-8",
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
