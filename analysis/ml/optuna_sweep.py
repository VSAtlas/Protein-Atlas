from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.split_manifest import load_locked_split_manifest
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier_core import _design_matrix, _model, _model_probabilities, _sample_weights
from analysis.ml.train_classifier_splits import _custom_split, _fixed_validation_split, _temporal_split

SUPPORTED_OPTUNA_PU_MODES = {"standard_binary", "positive_unlabeled_weighted", "case_control_matched"}
MINIMIZE_METRICS = {"Brier", "ECE", "adaptive_ECE", "MCE", "log_loss"}
OPTUNA_METRICS = ("AUPRC", "AUROC", "EF@1%", "EF@5%", "EF@10%", "Brier", "ECE", "adaptive_ECE", "MCE", "log_loss")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def _direction(metric: str) -> str:
    return "minimize" if metric in MINIMIZE_METRICS else "maximize"


def _sample_params(trial: Any, model_type: str) -> dict[str, object]:
    if model_type == "logistic_regression":
        return {"C": trial.suggest_float("C", 1e-3, 100.0, log=True)}
    if model_type in {"elastic_net", "elastic_net_logistic", "elastic_net_logistic_regression", "logistic_elastic_net"}:
        return {
            "C": trial.suggest_float("C", 1e-3, 100.0, log=True),
            "l1_ratio": trial.suggest_float("l1_ratio", 0.05, 0.95),
        }
    if model_type == "random_forest":
        return {
            "n_estimators": trial.suggest_int("n_estimators", 100, 600, step=50),
            "max_depth": trial.suggest_categorical("max_depth", [None, 4, 8, 16, 32]),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 1, 10),
            "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
            "max_features": trial.suggest_categorical("max_features", ["sqrt", "log2", None]),
        }
    if model_type in {"xgboost", "gradient_boosted_trees"}:
        return {
            "n_estimators": trial.suggest_int("n_estimators", 100, 600, step=50),
            "max_depth": trial.suggest_int("max_depth", 2, 8),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
            "subsample": trial.suggest_float("subsample", 0.6, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.6, 1.0),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
            "min_child_weight": trial.suggest_float("min_child_weight", 0.1, 10.0, log=True),
        }
    if model_type in {"lightgbm", "shallow_lightgbm", "lgbm"}:
        return {
            "n_estimators": trial.suggest_int("n_estimators", 100, 500, step=50),
            "max_depth": trial.suggest_int("max_depth", 2, 5),
            "num_leaves": trial.suggest_categorical("num_leaves", [3, 7, 15]),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.2, log=True),
            "subsample": trial.suggest_float("subsample", 0.6, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.6, 1.0),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
        }
    if model_type in {"catboost", "catboost_shallow"}:
        return {
            "iterations": trial.suggest_int("iterations", 100, 500, step=50),
            "depth": trial.suggest_int("depth", 2, 5),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.2, log=True),
            "l2_leaf_reg": trial.suggest_float("l2_leaf_reg", 1.0, 10.0),
        }
    if model_type in {"explainable_boosting_machine", "ebm", "interpret_ebm"}:
        return {
            "max_rounds": trial.suggest_int("max_rounds", 1000, 8000, step=1000),
            "learning_rate": trial.suggest_float("learning_rate", 0.001, 0.05, log=True),
            "interactions": trial.suggest_int("interactions", 0, 20),
            "max_bins": trial.suggest_categorical("max_bins", [64, 128, 256]),
        }
    raise ValueError(f"unsupported Optuna model_type: {model_type}")


def _prepare_validation_data(
    *,
    dataset_path: str | Path,
    label_col: str,
    feature_set: str,
    split_mode: str,
    seed: int,
    split_column: str | None,
    train_values: list[str] | None,
    test_values: list[str] | None,
    temporal_year_col: str | None,
    temporal_cutoff_year: int | None,
    temporal_max_missing_fraction: float,
    validation_fold_col: str | None,
    validation_fold_value: str | None,
    validation_fraction: float,
    exclude_features: list[str] | None,
    split_manifest: str | Path | None,
    allow_label_definition_features: bool,
    allow_partial_rescoring_features: bool,
) -> dict[str, Any]:
    df = pd.read_csv(dataset_path, low_memory=False)
    df["_atlas_observed_label"] = binary_label_series(df[label_col])
    excluded = effective_exclude_features(
        label_col,
        exclude_features,
        allow_label_definition_features=allow_label_definition_features,
    )
    features = [feature for feature in get_feature_set(feature_set) if feature in df.columns and feature not in excluded]
    assert_no_leakage(
        features,
        allow_label_definition_features=allow_label_definition_features,
        allow_partial_rescoring_features=allow_partial_rescoring_features,
    )
    data = df.loc[df["_atlas_observed_label"].notna()].copy()
    data[label_col] = data["_atlas_observed_label"]
    data = data.dropna(subset=[label_col]).copy()
    data[label_col] = data[label_col].astype(int)

    locked_validation_idx = pd.Index([])
    if split_manifest is not None:
        locked = load_locked_split_manifest(data, split_manifest, label_col=label_col, split_mode=split_mode)
        train_idx = locked["train_idx"]
        test_idx = locked["test_idx"]
        locked_validation_idx = locked["validation_idx"]
        split_summary = dict(locked.get("split_summary") or {})
        split_summary.setdefault("locked_split_manifest", str(split_manifest))
    else:
        temporal = _temporal_split(
            data,
            temporal_year_col,
            temporal_cutoff_year,
            max_missing_fraction=temporal_max_missing_fraction,
        )
        custom = _custom_split(data, split_column, train_values, test_values)
        split_summary: dict[str, object] | None
        if temporal is not None:
            train_idx, test_idx, split_summary = temporal
        elif custom is not None:
            train_idx, test_idx, split_summary = custom
        else:
            train_idx, test_idx = make_split(data, split_mode=split_mode, seed=seed)
            split_summary = None

    train = data.loc[train_idx]
    test = data.loc[test_idx]
    if split_summary is None:
        split_summary = split_overlap_summary(train, test, split_mode)
    if split_summary.get("overlaps") and not bool(split_summary.get("passes_holdout", False)):
        raise ValueError(f"holdout split leakage detected: {split_summary}")

    if len(locked_validation_idx):
        model_train = train.copy()
        validation = data.loc[locked_validation_idx].copy()
        validation_idx = locked_validation_idx
    else:
        model_train, validation, validation_idx = _fixed_validation_split(
            train,
            validation_fold_col=validation_fold_col,
            validation_fold_value=validation_fold_value,
            seed=seed,
            validation_fraction=validation_fraction,
        )
    if validation.empty:
        raise ValueError("Optuna sweep requires a non-empty fixed validation split")
    if validation[label_col].nunique() < 2:
        raise ValueError("Optuna validation split must contain both positive and negative labels")

    x_train = _design_matrix(model_train, features)
    x_validation = _design_matrix(validation, features).reindex(columns=x_train.columns, fill_value=0)
    return {
        "data": data,
        "train": train,
        "test": test,
        "model_train": model_train,
        "validation": validation,
        "validation_idx": validation_idx,
        "x_train": x_train,
        "x_validation": x_validation,
        "features": features,
        "excluded_features": sorted(excluded),
        "split_summary": split_summary,
        "split_manifest": str(split_manifest) if split_manifest is not None else None,
    }


def _trial_rows(study: Any) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for trial in study.trials:
        row: dict[str, object] = {
            "number": trial.number,
            "state": trial.state.name,
            "value": trial.value,
        }
        row.update({f"param_{key}": value for key, value in trial.params.items()})
        row.update({f"metric_{key}": value for key, value in trial.user_attrs.items() if key in OPTUNA_METRICS})
        rows.append(row)
    return rows


def _trial_state_counts(study: Any) -> dict[str, int]:
    counts: dict[str, int] = {}
    for trial in study.trials:
        state = str(trial.state.name).lower()
        counts[state] = counts.get(state, 0) + 1
    return counts


def _make_pruner(optuna: Any, name: str, *, seed: int) -> Any:
    if name == "none":
        return optuna.pruners.NopPruner()
    if name == "median":
        return optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=0)
    if name == "successive_halving":
        return optuna.pruners.SuccessiveHalvingPruner()
    if name == "hyperband":
        return optuna.pruners.HyperbandPruner()
    raise ValueError(f"unsupported Optuna pruner: {name}")


def run_optuna_sweep(
    *,
    dataset_path: str | Path,
    label_col: str,
    feature_set: str,
    model_type: str,
    split_mode: str,
    out_dir: str | Path,
    seed: int = 42,
    split_column: str | None = None,
    train_values: list[str] | None = None,
    test_values: list[str] | None = None,
    temporal_year_col: str | None = None,
    temporal_cutoff_year: int | None = None,
    temporal_max_missing_fraction: float = 0.20,
    validation_fold_col: str | None = None,
    validation_fold_value: str | None = None,
    validation_fraction: float = 0.15,
    n_trials: int = 25,
    timeout: int | None = None,
    metric: str = "AUPRC",
    class_weight: str | None = "balanced",
    pu_mode: str = "standard_binary",
    exclude_features: list[str] | None = None,
    split_manifest: str | Path | None = None,
    allow_label_definition_features: bool = False,
    allow_partial_rescoring_features: bool = False,
    study_name: str | None = None,
    storage: str | None = None,
    load_if_exists: bool = False,
    pruner: str = "median",
) -> dict[str, Any]:
    if metric not in OPTUNA_METRICS:
        raise ValueError(f"unsupported Optuna metric {metric!r}; choose one of {', '.join(OPTUNA_METRICS)}")
    if pu_mode not in SUPPORTED_OPTUNA_PU_MODES:
        raise ValueError(f"Optuna sweeps currently support only {sorted(SUPPORTED_OPTUNA_PU_MODES)}")
    if n_trials <= 0:
        raise ValueError("n_trials must be positive")
    try:
        import optuna
    except ImportError as exc:
        raise RuntimeError(
            "Optuna is required only for analysis.cli.run_ml_optuna_sweep; "
            "install the optional ML ops tool bundle before running sweeps."
        ) from exc

    hpo_dir = Path(out_dir) / "hpo" / "optuna"
    hpo_dir.mkdir(parents=True, exist_ok=True)
    prepared = _prepare_validation_data(
        dataset_path=dataset_path,
        label_col=label_col,
        feature_set=feature_set,
        split_mode=split_mode,
        seed=seed,
        split_column=split_column,
        train_values=train_values,
        test_values=test_values,
        temporal_year_col=temporal_year_col,
        temporal_cutoff_year=temporal_cutoff_year,
        temporal_max_missing_fraction=temporal_max_missing_fraction,
        validation_fold_col=validation_fold_col,
        validation_fold_value=validation_fold_value,
        validation_fraction=validation_fraction,
        exclude_features=exclude_features,
        split_manifest=split_manifest,
        allow_label_definition_features=allow_label_definition_features,
        allow_partial_rescoring_features=allow_partial_rescoring_features,
    )
    direction = _direction(metric)
    sampler = optuna.samplers.TPESampler(seed=seed)
    pruner_obj = _make_pruner(optuna, pruner, seed=seed)
    study = optuna.create_study(
        direction=direction,
        study_name=study_name,
        storage=storage,
        load_if_exists=load_if_exists,
        sampler=sampler,
        pruner=pruner_obj,
    )

    def objective(trial: Any) -> float:
        params = _sample_params(trial, model_type)
        model = _model(model_type, seed, class_weight=class_weight, model_params=params)
        sample_weight = _sample_weights(prepared["model_train"], label_col, pu_mode=pu_mode)
        if sample_weight is not None:
            model.fit(prepared["x_train"], prepared["model_train"][label_col], sample_weight=sample_weight)
        else:
            model.fit(prepared["x_train"], prepared["model_train"][label_col])
        probs = _model_probabilities(model, prepared["x_validation"])
        metrics = calibration_metrics(probs, prepared["validation"][label_col].astype(int).tolist())
        for key, value in metrics.items():
            trial.set_user_attr(key, float(value))
        if metric not in metrics:
            raise ValueError(f"metric {metric!r} not produced by calibration_metrics")
        objective_value = float(metrics[metric])
        trial.report(objective_value, step=0)
        if trial.should_prune():
            raise optuna.TrialPruned()
        return objective_value

    study.optimize(objective, n_trials=n_trials, timeout=timeout)
    trials_path = hpo_dir / "optuna_trials.csv"
    pd.DataFrame(_trial_rows(study)).to_csv(trials_path, index=False)
    completed = [trial for trial in study.trials if trial.state.name == "COMPLETE" and trial.value is not None]
    if not completed:
        trial_state_counts = _trial_state_counts(study)
        manifest = {
            "status": "failed",
            "reason": "no_completed_trials",
            "backend": "optuna",
            "study_name": study.study_name,
            "trials_path": str(trials_path),
            "sampler": "TPESampler",
            "pruner": pruner_obj.__class__.__name__,
            "trial_state_counts": trial_state_counts,
            "pruning_status": "pruned_trials_present" if trial_state_counts.get("pruned", 0) else "no_pruned_trials",
        }
        _write_json(hpo_dir / "optuna_sweep_manifest.json", manifest)
        raise RuntimeError("Optuna sweep completed without successful trials")

    best = study.best_trial
    trial_state_counts = _trial_state_counts(study)
    manifest = {
        "status": "run",
        "backend": "optuna",
        "study_name": study.study_name,
        "direction": direction,
        "sampler": "TPESampler",
        "pruner": pruner_obj.__class__.__name__,
        "selection_metric": metric,
        "n_trials_requested": int(n_trials),
        "n_trials": int(len(completed)),
        "trial_state_counts": trial_state_counts,
        "pruning_status": "pruned_trials_present" if trial_state_counts.get("pruned", 0) else "no_pruned_trials",
        "timeout_seconds": timeout,
        "best_trial_number": int(best.number),
        "best_value": float(best.value),
        "best_params": dict(best.params),
        "dataset": str(dataset_path),
        "label_col": label_col,
        "feature_set": feature_set,
        "model_type": model_type,
        "split_mode": split_mode,
        "class_weight": class_weight or "none",
        "pu_mode": pu_mode,
        "seed": int(seed),
        "n_outer_train_rows": int(len(prepared["train"])),
        "n_model_train_rows": int(len(prepared["model_train"])),
        "n_validation_rows": int(len(prepared["validation"])),
        "n_outer_test_rows": int(len(prepared["test"])),
        "validation_fold_col": validation_fold_col,
        "validation_fold_value": validation_fold_value,
        "validation_fraction": float(validation_fraction),
        "features": list(prepared["features"]),
        "excluded_features": list(prepared["excluded_features"]),
        "split_summary": prepared["split_summary"],
        "trials_path": str(trials_path),
        "manifest_path": str(hpo_dir / "optuna_sweep_manifest.json"),
        "policy": [
            "Trials optimize only the fixed inner validation partition.",
            "The outer holdout is not used by Optuna and remains reserved for final train_ml_model scoring.",
            "Optuna is imported only inside run_optuna_sweep so baseline ML imports do not require it.",
            "A single validation checkpoint is reported to the configured Optuna pruner; iterative estimators can extend this to true early stopping.",
        ],
    }
    _write_json(hpo_dir / "best_params.json", dict(best.params))
    _write_json(hpo_dir / "optuna_sweep_manifest.json", manifest)
    return manifest
