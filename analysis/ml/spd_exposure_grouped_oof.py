from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import (
    average_precision_score,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold
from threadpoolctl import threadpool_limits

from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import find_leaky_features
from analysis.ml.model_run_ledger import LEDGER_COLUMNS
from analysis.ml.train_classifier_core import (
    _fit_design_matrix,
    _transform_design_matrix,
)


_RESIDUAL_SIGMA_FLOOR = 0.25
_DEFAULT_TOP_K = (10, 20, 100)
_SUPPORTED_MODEL_TYPES = frozenset({"ridge", "random_forest", "lightgbm"})
_CONTROLLED_MODEL_PARAMS = frozenset(
    {
        "random_state",
        "seed",
        "bagging_seed",
        "feature_fraction_seed",
        "data_random_seed",
        "drop_seed",
        "n_jobs",
        "nthread",
        "num_threads",
        "deterministic",
        "force_col_wise",
        "force_row_wise",
    }
)
_LIGHTGBM_ALLOWED_PARAMS = frozenset(
    {
        "boosting_type",
        "colsample_bytree",
        "learning_rate",
        "max_bin",
        "max_depth",
        "min_child_samples",
        "min_child_weight",
        "min_split_gain",
        "n_estimators",
        "num_leaves",
        "objective",
        "reg_alpha",
        "reg_lambda",
        "subsample",
        "subsample_freq",
    }
)
_AUDIT_ONLY_PROTEIN_BINDING_COLUMNS = ("pk_context_protein_binding_percent",)
_ID_COLUMNS = (
    "canonical_pair_key",
    "drug_id",
    "drug_name",
    "target_id",
    "pdb_id",
    "target_family",
    "protein_class",
    "scaffold_key",
    "chemical_cluster",
    "ligand_chemotype",
    "label_source",
    "source_family",
)
_PERMUTATION_COLUMNS = (
    "stage",
    "feature",
    "fold",
    "repeat",
    "permutation_unit",
    "n_stage_evaluation",
    "n_downstream_evaluation",
    "baseline_stage_mae",
    "permuted_stage_mae",
    "stage_mae_increase",
    "baseline_stage_rmse",
    "permuted_stage_rmse",
    "stage_rmse_increase",
    "baseline_stage_r2",
    "permuted_stage_r2",
    "stage_r2_drop",
    "baseline_AUPRC",
    "permuted_AUPRC",
    "downstream_AUPRC_drop",
    "baseline_AUROC",
    "permuted_AUROC",
    "downstream_AUROC_drop",
    "baseline_Brier",
    "permuted_Brier",
    "downstream_Brier_increase",
    "baseline_ECE",
    "permuted_ECE",
    "downstream_ECE_increase",
)


@dataclass
class _StageFit:
    model: Any
    preprocessing: Mapping[str, Any]
    train_prediction: pd.Series
    test_prediction: pd.Series
    residual_sigma: float
    effective_model_params: Mapping[str, Any]


def _validate_max_threads(max_threads: int) -> int:
    value = int(max_threads)
    if not 1 <= value <= 4:
        raise ValueError("max_threads must be between 1 and 4")
    return value


def _build_regressor(
    model_type: str,
    *,
    seed: int,
    model_params: Mapping[str, Any] | None,
    max_threads: int,
) -> tuple[Any, dict[str, Any]]:
    """Build one explicitly configured stage regressor and return its full settings."""
    if model_type not in _SUPPORTED_MODEL_TYPES:
        allowed = ", ".join(sorted(_SUPPORTED_MODEL_TYPES))
        raise ValueError(f"model_type must be one of: {allowed}")
    threads = _validate_max_threads(max_threads)
    requested = dict(model_params or {})
    controlled = sorted(set(requested) & _CONTROLLED_MODEL_PARAMS)
    if controlled:
        raise ValueError(
            "model_params cannot override reproducibility/thread controls: "
            + ", ".join(controlled)
        )

    if model_type == "ridge":
        params: dict[str, Any] = {"alpha": 1.0}
        params.update(requested)
        model = Ridge(random_state=int(seed), **params)
    elif model_type == "random_forest":
        from sklearn.ensemble import RandomForestRegressor

        params = {"n_estimators": 300, "min_samples_leaf": 3}
        params.update(requested)
        model = RandomForestRegressor(
            random_state=int(seed),
            n_jobs=threads,
            **params,
        )
    else:
        unknown = sorted(set(requested) - _LIGHTGBM_ALLOWED_PARAMS)
        if unknown:
            raise ValueError("unsupported LightGBM model_params: " + ", ".join(unknown))
        try:
            from lightgbm import LGBMRegressor
        except ImportError as exc:
            raise RuntimeError(
                "lightgbm regressor requires the optional lightgbm package"
            ) from exc
        params = {
            "boosting_type": "gbdt",
            "objective": "regression",
            "n_estimators": 200,
            "learning_rate": 0.05,
            "max_depth": 3,
            "num_leaves": 7,
            "min_child_samples": 10,
            "subsample": 0.8,
            "subsample_freq": 1,
            "colsample_bytree": 0.8,
            "reg_alpha": 0.1,
            "reg_lambda": 1.0,
            "verbosity": -1,
        }
        params.update(requested)
        model = LGBMRegressor(
            random_state=int(seed),
            bagging_seed=int(seed),
            feature_fraction_seed=int(seed),
            data_random_seed=int(seed),
            deterministic=True,
            force_col_wise=True,
            n_jobs=threads,
            **params,
        )
    effective = dict(model.get_params(deep=False))
    return model, effective


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_state() -> tuple[str, bool | None]:
    root = Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception:
        return "unknown", None
    git_commit = commit.stdout.strip() if commit.returncode == 0 else "unknown"
    dirty = bool(status.stdout.strip()) if status.returncode == 0 else None
    return git_commit, dirty


def _first_numeric(
    frame: pd.DataFrame, names: Sequence[str]
) -> tuple[pd.Series, str | None]:
    for name in names:
        if name in frame.columns:
            return pd.to_numeric(frame[name], errors="coerce"), name
    return pd.Series(np.nan, index=frame.index, dtype=float), None


def _activity_relation(frame: pd.DataFrame) -> tuple[pd.Series, str | None]:
    for name in ("spd_activity_relation", "activity_relation", "relation"):
        if name in frame.columns:
            values = frame[name].fillna("").astype(str).str.strip()
            return values, name
    return pd.Series("=", index=frame.index, dtype=object), None


def _deduplicate(values: Sequence[str], *, name: str) -> list[str]:
    cleaned = [str(value).strip() for value in values if str(value).strip()]
    if not cleaned:
        raise ValueError(f"at least one {name} is required")
    duplicates = sorted({value for value in cleaned if cleaned.count(value) > 1})
    if duplicates:
        raise ValueError(
            f"duplicate {name} values are not allowed: {', '.join(duplicates)}"
        )
    return cleaned


def _is_fraction_unbound_predictor(feature: str) -> bool:
    return "fraction_unbound" in feature.lower()


def _is_protein_binding_audit_only(feature: str) -> bool:
    return "protein_binding" in feature.lower()


def _is_direct_response_feature(feature: str, label_col: str) -> bool:
    normalized = feature.lower().replace("-", "_")
    if normalized == label_col.lower().replace("-", "_"):
        return True
    if "ac50" in normalized or normalized in {"pic50", "log_pic50"}:
        return True
    if "cmax" in normalized:
        return True
    if "exposure" in normalized and any(
        token in normalized
        for token in ("label", "margin", "relevant", "probability", "prediction")
    ):
        return True
    return "activity_relation" in normalized or "censor" in normalized


def _validate_features(
    frame: pd.DataFrame,
    *,
    label_col: str,
    potency_features: Sequence[str],
    pk_features: Sequence[str],
) -> tuple[list[str], list[str], list[str]]:
    potency = _deduplicate(potency_features, name="--potency-feature")
    pk = _deduplicate(pk_features, name="--pk-feature")
    requested = [*potency, *pk]
    missing = sorted(set(requested) - set(frame.columns))
    if missing:
        raise ValueError(f"requested feature columns are missing: {', '.join(missing)}")
    all_missing = sorted(
        {feature for feature in requested if frame[feature].notna().sum() == 0}
    )
    if all_missing:
        raise ValueError(
            f"requested feature columns are entirely missing: {', '.join(all_missing)}"
        )

    fraction_unbound = sorted(
        {feature for feature in requested if _is_fraction_unbound_predictor(feature)}
    )
    protein_binding = sorted(
        {feature for feature in requested if _is_protein_binding_audit_only(feature)}
    )
    if protein_binding:
        complement_note = (
            " Do not include it with fraction-unbound because the two fields are "
            "deterministic complements."
            if fraction_unbound
            else ""
        )
        raise ValueError(
            "protein-binding percent is audit-only and cannot be a predictive "
            f"feature: {', '.join(protein_binding)}.{complement_note}"
        )
    explicitly_allowed = fraction_unbound
    leaky = set(find_leaky_features(requested)) - set(explicitly_allowed)
    leaky.update(
        feature
        for feature in requested
        if _is_direct_response_feature(feature, label_col)
    )
    if leaky:
        raise ValueError(
            "direct response, label, identity, or provenance features are not allowed: "
            + ", ".join(sorted(leaky))
        )
    return potency, pk, explicitly_allowed


def _prepare_data(
    frame: pd.DataFrame,
    *,
    label_col: str,
    group_col: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = [label_col, "drug_id", group_col]
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"dataset is missing required columns: {', '.join(missing)}")

    ac50, ac50_col = _first_numeric(frame, ("spd_ac50_uM", "ac50_um", "ac50_uM"))
    if ac50.isna().all() and "ac50_nM" in frame.columns:
        ac50 = pd.to_numeric(frame["ac50_nM"], errors="coerce") / 1000.0
        ac50_col = "ac50_nM"
    free_cmax, free_cmax_col = _first_numeric(
        frame,
        (
            "free_cmax_um",
            "free_cmax_uM",
            "free_cmax",
            "pk_context_free_cmax_um",
        ),
    )
    if ac50_col is None:
        raise ValueError("dataset has no supported AC50 target column")
    if free_cmax_col is None:
        raise ValueError("dataset has no supported free-Cmax target column")

    relation, relation_col = _activity_relation(frame)
    labels = binary_label_series(frame[label_col])
    data = frame.copy()
    data["_source_index"] = np.arange(len(data), dtype=int)
    data["_observed_ac50_um"] = ac50.astype(float)
    data["_observed_free_cmax_um"] = free_cmax.astype(float)
    data["_observed_exposure_label"] = labels
    data["_activity_relation"] = relation

    group_present = data[group_col].notna() & data[group_col].astype(
        str
    ).str.strip().ne("")
    drug_present = data["drug_id"].notna() & data["drug_id"].astype(str).str.strip().ne(
        ""
    )
    usable_mask = (
        data["_observed_ac50_um"].gt(0)
        & data["_observed_free_cmax_um"].gt(0)
        & data["_observed_exposure_label"].notna()
        & group_present
        & drug_present
    )
    usable = data.loc[usable_mask].copy()
    if len(usable) < 20:
        raise ValueError(
            "grouped OOF exposure evaluation requires at least 20 usable rows; "
            f"found {len(usable)}"
        )
    if usable["_observed_exposure_label"].nunique() < 2:
        raise ValueError("grouped OOF exposure evaluation requires both label classes")

    groups_per_drug = usable.groupby("drug_id", dropna=False)[group_col].nunique(
        dropna=False
    )
    crossing = groups_per_drug[groups_per_drug.gt(1)]
    if not crossing.empty:
        examples = ", ".join(str(value) for value in crossing.index[:10])
        raise ValueError(
            f"{group_col!r} is not constant within drug_id; this would leak the drug-level "
            f"free-Cmax target across folds. Example drugs: {examples}"
        )

    usable["_log10_ac50_um"] = np.log10(usable["_observed_ac50_um"].astype(float))
    usable["_log10_free_cmax_um"] = np.log10(
        usable["_observed_free_cmax_um"].astype(float)
    )
    contract = {
        "ac50_source_column": ac50_col,
        "free_cmax_source_column": free_cmax_col,
        "activity_relation_source_column": relation_col,
        "n_rows_full": int(len(frame)),
        "n_rows_usable": int(len(usable)),
        "n_rows_excluded": int(len(frame) - len(usable)),
        "n_positive_usable": int(usable["_observed_exposure_label"].eq(1).sum()),
        "n_negative_usable": int(usable["_observed_exposure_label"].eq(0).sum()),
        "n_exact_ac50_usable": int(usable["_activity_relation"].eq("=").sum()),
        "n_right_censored_ac50_usable": int(
            usable["_activity_relation"].str.startswith(">").sum()
        ),
        "n_drugs_usable": int(usable["drug_id"].nunique()),
        "n_groups_usable": int(usable[group_col].nunique()),
    }
    return usable, contract


def _make_folds(
    data: pd.DataFrame,
    *,
    group_col: str,
    n_splits: int,
    seed: int,
) -> tuple[list[tuple[pd.Index, pd.Index]], str, int]:
    if n_splits < 2:
        raise ValueError("n_splits must be at least 2")
    groups = data[group_col].astype(str)
    effective = min(int(n_splits), int(groups.nunique()))
    if effective < 2:
        raise ValueError(f"{group_col!r} must contain at least two groups")
    labels = data["_observed_exposure_label"].astype(int)
    strategy = "StratifiedGroupKFold"
    try:
        splitter = StratifiedGroupKFold(
            n_splits=effective,
            shuffle=True,
            random_state=seed,
        )
        positional = list(splitter.split(data, labels, groups))
    except ValueError:
        strategy = "GroupKFold_fallback"
        splitter = GroupKFold(n_splits=effective)
        positional = list(splitter.split(data, labels, groups))

    folds = [
        (data.index[train_pos], data.index[test_pos])
        for train_pos, test_pos in positional
    ]
    assigned = pd.Series(-1, index=data.index, dtype=int)
    for fold, (train_idx, test_idx) in enumerate(folds):
        if assigned.loc[test_idx].ne(-1).any():
            raise AssertionError("an OOF row was assigned to more than one test fold")
        train_groups = set(groups.loc[train_idx])
        test_groups = set(groups.loc[test_idx])
        if train_groups & test_groups:
            raise AssertionError(f"group overlap detected in fold {fold}")
        assigned.loc[test_idx] = fold
    if assigned.lt(0).any():
        raise AssertionError("not every usable row received an OOF fold")
    drug_fold_counts = assigned.groupby(data["drug_id"].astype(str)).nunique()
    if drug_fold_counts.gt(1).any():
        raise AssertionError("a drug was assigned to more than one OOF test fold")
    return folds, strategy, effective


def _fold_assignment_frame(
    data: pd.DataFrame,
    folds: Sequence[tuple[pd.Index, pd.Index]],
    *,
    group_col: str,
    label_col: str,
) -> pd.DataFrame:
    fold_by_index = pd.Series(-1, index=data.index, dtype=int)
    for fold, (_, test_idx) in enumerate(folds):
        fold_by_index.loc[test_idx] = int(fold)
    if fold_by_index.lt(0).any():
        raise AssertionError("not every usable row received a fold assignment")
    columns = [column for column in _ID_COLUMNS if column in data.columns]
    assignments = data[columns].copy()
    assignments.insert(0, "_source_index", data["_source_index"].astype(int))
    assignments["fold"] = fold_by_index.astype(int)
    assignments["group_col"] = group_col
    assignments["group_value"] = data[group_col].astype(str)
    assignments[label_col] = data["_observed_exposure_label"].astype(int)
    return assignments.sort_values("_source_index").reset_index(drop=True)


def _fold_signature(assignments: pd.DataFrame) -> str:
    columns = ["_source_index", "fold", "group_col", "group_value"]
    canonical = assignments[columns].copy().sort_values("_source_index")
    canonical["_source_index"] = canonical["_source_index"].astype(int)
    canonical["fold"] = canonical["fold"].astype(int)
    canonical["group_col"] = canonical["group_col"].astype(str)
    canonical["group_value"] = canonical["group_value"].astype(str)
    payload = canonical.to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _load_saved_folds(
    data: pd.DataFrame,
    assignments_path: Path,
    *,
    group_col: str,
) -> tuple[list[tuple[pd.Index, pd.Index]], str, int]:
    if not assignments_path.is_file():
        raise FileNotFoundError(f"fold assignments not found: {assignments_path}")
    saved = pd.read_csv(assignments_path, low_memory=False)
    required = {"_source_index", "fold", "group_col", "group_value"}
    missing = sorted(required - set(saved.columns))
    if missing:
        raise ValueError(
            "saved fold assignments are missing columns: " + ", ".join(missing)
        )
    if saved["_source_index"].duplicated().any():
        raise ValueError("saved fold assignments contain duplicate _source_index rows")
    source_numeric = pd.to_numeric(saved["_source_index"], errors="coerce")
    fold_numeric = pd.to_numeric(saved["fold"], errors="coerce")
    if source_numeric.isna().any() or (source_numeric % 1).ne(0).any():
        raise ValueError("saved _source_index values must be integers")
    if fold_numeric.isna().any() or (fold_numeric % 1).ne(0).any():
        raise ValueError("saved fold values must be integers")
    saved = saved.copy()
    saved["_source_index"] = source_numeric.astype(int)
    saved["fold"] = fold_numeric.astype(int)
    if saved["fold"].lt(0).any():
        raise ValueError("saved fold values must be non-negative")

    expected_sources = set(data["_source_index"].astype(int))
    saved_sources = set(saved["_source_index"])
    if saved_sources != expected_sources:
        missing_sources = len(expected_sources - saved_sources)
        extra_sources = len(saved_sources - expected_sources)
        raise ValueError(
            "saved folds do not exactly cover current usable rows: "
            f"missing={missing_sources}, extra={extra_sources}"
        )
    declared_group_cols = set(saved["group_col"].astype(str))
    if declared_group_cols != {group_col}:
        raise ValueError(
            "saved folds use a different group column: "
            + ", ".join(sorted(declared_group_cols))
        )

    current = data.set_index(data["_source_index"].astype(int), drop=False)
    saved = saved.set_index("_source_index", drop=False).loc[sorted(expected_sources)]
    expected_groups = current.loc[saved.index, group_col].astype(str)
    observed_groups = saved["group_value"].astype(str)
    mismatch = observed_groups.ne(expected_groups)
    if mismatch.any():
        examples = ", ".join(str(value) for value in saved.index[mismatch][:10])
        raise ValueError(
            "saved fold group values do not match the current dataset at source rows: "
            + examples
        )
    group_fold_counts = saved.groupby("group_value")["fold"].nunique()
    if group_fold_counts.gt(1).any():
        raise ValueError("a saved group is assigned to more than one test fold")

    fold_ids = sorted(int(value) for value in saved["fold"].unique())
    if fold_ids != list(range(len(fold_ids))) or len(fold_ids) < 2:
        raise ValueError(
            "saved fold IDs must be contiguous from 0 and contain at least two folds"
        )
    source_to_index = dict(
        zip(data["_source_index"].astype(int), data.index, strict=True)
    )
    folds: list[tuple[pd.Index, pd.Index]] = []
    for fold in fold_ids:
        test_sources = saved.loc[saved["fold"].eq(fold), "_source_index"]
        test_idx = pd.Index([source_to_index[int(value)] for value in test_sources])
        train_idx = data.index[~data.index.isin(test_idx)]
        train_groups = set(data.loc[train_idx, group_col].astype(str))
        test_groups = set(data.loc[test_idx, group_col].astype(str))
        if train_groups & test_groups:
            raise ValueError(f"saved fold {fold} has train/test group overlap")
        folds.append((train_idx, test_idx))
    return folds, "SavedFoldAssignments", len(folds)


def save_spd_exposure_group_folds(
    dataset_path: str | Path,
    assignments_path: str | Path,
    *,
    label_col: str = "spd_exposure_label",
    group_col: str = "drug_id",
    n_splits: int = 5,
    seed: int = 42,
) -> dict[str, Any]:
    """Materialize one reusable grouped-OOF assignment for settings comparisons."""
    dataset = Path(dataset_path)
    if not dataset.is_file():
        raise FileNotFoundError(f"dataset not found: {dataset}")
    full = pd.read_csv(dataset, low_memory=False)
    data, target_contract = _prepare_data(
        full,
        label_col=label_col,
        group_col=group_col,
    )
    folds, strategy, effective_splits = _make_folds(
        data,
        group_col=group_col,
        n_splits=n_splits,
        seed=seed,
    )
    assignments = _fold_assignment_frame(
        data,
        folds,
        group_col=group_col,
        label_col=label_col,
    )
    destination = Path(assignments_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(destination, index=False)
    manifest_path = destination.with_name(f"{destination.stem}_manifest.json")
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset_path": str(dataset),
        "dataset_file_sha256": _sha256_file(dataset),
        "assignments_path": str(destination),
        "assignments_file_sha256": _sha256_file(destination),
        "fold_signature_sha256": _fold_signature(assignments),
        "label_col": label_col,
        "group_col": group_col,
        "strategy": strategy,
        "requested_n_splits": int(n_splits),
        "effective_n_splits": int(effective_splits),
        "seed": int(seed),
        "n_rows": int(len(assignments)),
        "n_groups": int(assignments["group_value"].nunique()),
        "target_contract": target_contract,
    }
    _write_json(manifest_path, manifest)
    manifest["manifest_path"] = str(manifest_path)
    return manifest


def _fit_stage(
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    features: list[str],
    target_col: str,
    model_type: str,
    seed: int,
    model_params: Mapping[str, Any] | None,
    max_threads: int,
) -> _StageFit:
    x_train, preprocessing = _fit_design_matrix(train, features)
    x_test = _transform_design_matrix(test, preprocessing)
    model, effective_model_params = _build_regressor(
        model_type,
        seed=seed,
        model_params=model_params,
        max_threads=max_threads,
    )
    target = pd.to_numeric(train[target_col], errors="coerce").astype(float)
    with threadpool_limits(limits=max_threads):
        model.fit(x_train, target)
        train_values = model.predict(x_train)
        test_values = model.predict(x_test)
    train_prediction = pd.Series(train_values, index=train.index, dtype=float)
    test_prediction = pd.Series(test_values, index=test.index, dtype=float)
    sigma = float((target - train_prediction).std(ddof=1))
    if not math.isfinite(sigma) or sigma < _RESIDUAL_SIGMA_FLOOR:
        sigma = _RESIDUAL_SIGMA_FLOOR
    return _StageFit(
        model=model,
        preprocessing=preprocessing,
        train_prediction=train_prediction,
        test_prediction=test_prediction,
        residual_sigma=sigma,
        effective_model_params=effective_model_params,
    )


def _unique_drug_table(frame: pd.DataFrame, features: Sequence[str]) -> pd.DataFrame:
    columns = ["drug_id", "_log10_free_cmax_um", *features]
    work = frame[columns].dropna(subset=["drug_id", "_log10_free_cmax_um"]).copy()
    if work.empty:
        return work
    aggregation: dict[str, Any] = {"_log10_free_cmax_um": "median"}
    for feature in features:
        numeric = pd.to_numeric(work[feature], errors="coerce")
        if (
            numeric.notna().sum() == work[feature].notna().sum()
            and numeric.notna().any()
        ):
            work[feature] = numeric
            aggregation[feature] = "median"
        else:
            aggregation[feature] = "first"
    return work.groupby("drug_id", as_index=False, dropna=False).agg(aggregation)


def _combine_probability(
    log_ac50: pd.Series,
    log_free_cmax: pd.Series,
    *,
    ac50_sigma: float,
    free_cmax_sigma: float,
) -> pd.Series:
    sigma = math.sqrt(ac50_sigma**2 + free_cmax_sigma**2)
    z_score = (1.0 + log_free_cmax.astype(float) - log_ac50.astype(float)) / sigma
    probability = z_score.map(
        lambda value: 0.5 * (1.0 + math.erf(float(value) / math.sqrt(2.0)))
    )
    return probability.clip(1e-9, 1.0 - 1e-9)


def _expected_calibration_error(
    labels: pd.Series,
    probabilities: pd.Series,
    *,
    n_bins: int = 10,
) -> float:
    work = pd.DataFrame(
        {
            "label": pd.to_numeric(labels, errors="coerce"),
            "probability": pd.to_numeric(probabilities, errors="coerce"),
        }
    ).dropna()
    if work.empty:
        return math.nan
    clipped = work["probability"].clip(0.0, 1.0)
    bins = np.minimum((clipped.to_numpy() * n_bins).astype(int), n_bins - 1)
    work["bin"] = bins
    error = 0.0
    for _, group in work.groupby("bin"):
        error += (len(group) / len(work)) * abs(
            float(group["probability"].mean() - group["label"].mean())
        )
    return float(error)


def _binary_metrics(labels: pd.Series, probabilities: pd.Series) -> dict[str, float]:
    work = pd.DataFrame(
        {
            "label": pd.to_numeric(labels, errors="coerce"),
            "probability": pd.to_numeric(probabilities, errors="coerce"),
        }
    ).dropna()
    if work.empty:
        return {
            "AUROC": math.nan,
            "AUPRC": math.nan,
            "Brier": math.nan,
            "ECE": math.nan,
            "prevalence": math.nan,
            "null_Brier": math.nan,
            "Brier_skill": math.nan,
        }
    labels_int = work["label"].astype(int)
    scores = work["probability"].astype(float).clip(0.0, 1.0)
    prevalence = float(labels_int.mean())
    brier = float(((scores - labels_int) ** 2).mean())
    null_brier = float(((prevalence - labels_int) ** 2).mean())
    two_classes = labels_int.nunique() == 2
    return {
        "AUROC": float(roc_auc_score(labels_int, scores)) if two_classes else math.nan,
        "AUPRC": (
            float(average_precision_score(labels_int, scores))
            if two_classes
            else math.nan
        ),
        "Brier": brier,
        "ECE": _expected_calibration_error(labels_int, scores),
        "prevalence": prevalence,
        "null_Brier": null_brier,
        "Brier_skill": 1.0 - (brier / null_brier) if null_brier > 0 else math.nan,
    }


def _regression_metrics(
    observed: pd.Series, predicted: pd.Series
) -> dict[str, float | int]:
    work = pd.DataFrame(
        {
            "observed": pd.to_numeric(observed, errors="coerce"),
            "predicted": pd.to_numeric(predicted, errors="coerce"),
        }
    ).dropna()
    if work.empty:
        return {
            "n": 0,
            "mae": math.nan,
            "rmse": math.nan,
            "r2": math.nan,
            "median_fold_error": math.nan,
        }
    residual = work["observed"] - work["predicted"]
    can_score_r2 = len(work) > 1 and work["observed"].nunique() > 1
    return {
        "n": int(len(work)),
        "mae": float(mean_absolute_error(work["observed"], work["predicted"])),
        "rmse": float(
            math.sqrt(mean_squared_error(work["observed"], work["predicted"]))
        ),
        "r2": (
            float(r2_score(work["observed"], work["predicted"]))
            if can_score_r2
            else math.nan
        ),
        "median_fold_error": float((10.0 ** residual.abs()).median()),
    }


def _top_k_rows(
    labels: pd.Series,
    probabilities: pd.Series,
    *,
    top_k: Sequence[int],
    scope: str,
    fold: int | None,
) -> list[dict[str, Any]]:
    work = pd.DataFrame(
        {
            "label": pd.to_numeric(labels, errors="coerce"),
            "probability": pd.to_numeric(probabilities, errors="coerce"),
        }
    ).dropna()
    work = work.sort_values("probability", ascending=False, kind="mergesort")
    prevalence = float(work["label"].mean()) if not work.empty else math.nan
    total_positive = int(work["label"].eq(1).sum())
    rows: list[dict[str, Any]] = []
    for requested in top_k:
        effective = min(int(requested), len(work))
        selected = work.head(effective)
        recovered = int(selected["label"].eq(1).sum())
        precision = float(selected["label"].mean()) if effective else math.nan
        rows.append(
            {
                "scope": scope,
                "fold": fold,
                "requested_k": int(requested),
                "effective_k": int(effective),
                "n_rows": int(len(work)),
                "n_positive": total_positive,
                "positives_recovered": recovered,
                "precision_at_k": precision,
                "recall_at_k": recovered / total_positive
                if total_positive
                else math.nan,
                "enrichment_at_k": (
                    precision / prevalence
                    if effective and prevalence > 0 and math.isfinite(precision)
                    else math.nan
                ),
                "prevalence": prevalence,
            }
        )
    return rows


def _add_top_k_columns(
    metrics: dict[str, Any], rows: Sequence[Mapping[str, Any]]
) -> None:
    for row in rows:
        suffix = str(row["requested_k"])
        metrics[f"precision_at_{suffix}"] = row["precision_at_k"]
        metrics[f"recall_at_{suffix}"] = row["recall_at_k"]
        metrics[f"enrichment_at_{suffix}"] = row["enrichment_at_k"]
        metrics[f"positives_recovered_at_{suffix}"] = row["positives_recovered"]
        metrics[f"effective_k_{suffix}"] = row["effective_k"]


def _permutation_rng(
    *,
    seed: int,
    fold: int,
    stage: str,
    feature: str,
    repeat: int,
) -> np.random.Generator:
    token = f"{seed}|{fold}|{stage}|{feature}|{repeat}".encode()
    derived = int.from_bytes(hashlib.sha256(token).digest()[:8], "little")
    return np.random.default_rng(derived)


def _permute_potency_feature(
    frame: pd.DataFrame,
    *,
    feature: str,
    unit: str,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, str]:
    permuted = frame.copy()
    if unit == "drug_block":
        representatives = frame.groupby("drug_id", sort=False, dropna=False)[
            feature
        ].agg(lambda values: values.iloc[0])
        donor_values = representatives.iloc[
            rng.permutation(len(representatives))
        ].to_numpy()
        mapping = dict(
            zip(representatives.index.astype(str), donor_values, strict=True)
        )
        permuted[feature] = frame["drug_id"].astype(str).map(mapping)
        return permuted, "drug_block"
    if unit != "row":
        raise ValueError(f"unsupported potency permutation unit: {unit}")
    permuted[feature] = rng.permutation(frame[feature].to_numpy(copy=True))
    return permuted, "row"


def _permutation_record(
    *,
    stage: str,
    feature: str,
    fold: int,
    repeat: int,
    permutation_unit: str,
    baseline_stage: Mapping[str, Any],
    permuted_stage: Mapping[str, Any],
    baseline_downstream: Mapping[str, float],
    permuted_downstream: Mapping[str, float],
    n_downstream: int,
) -> dict[str, Any]:
    return {
        "stage": stage,
        "feature": feature,
        "fold": int(fold),
        "repeat": int(repeat),
        "permutation_unit": permutation_unit,
        "n_stage_evaluation": int(baseline_stage.get("n", 0)),
        "n_downstream_evaluation": int(n_downstream),
        "baseline_stage_mae": baseline_stage.get("mae"),
        "permuted_stage_mae": permuted_stage.get("mae"),
        "stage_mae_increase": (
            float(permuted_stage["mae"]) - float(baseline_stage["mae"])
        ),
        "baseline_stage_rmse": baseline_stage.get("rmse"),
        "permuted_stage_rmse": permuted_stage.get("rmse"),
        "stage_rmse_increase": (
            float(permuted_stage["rmse"]) - float(baseline_stage["rmse"])
        ),
        "baseline_stage_r2": baseline_stage.get("r2"),
        "permuted_stage_r2": permuted_stage.get("r2"),
        "stage_r2_drop": float(baseline_stage["r2"]) - float(permuted_stage["r2"]),
        "baseline_AUPRC": baseline_downstream.get("AUPRC"),
        "permuted_AUPRC": permuted_downstream.get("AUPRC"),
        "downstream_AUPRC_drop": (
            float(baseline_downstream["AUPRC"]) - float(permuted_downstream["AUPRC"])
        ),
        "baseline_AUROC": baseline_downstream.get("AUROC"),
        "permuted_AUROC": permuted_downstream.get("AUROC"),
        "downstream_AUROC_drop": (
            float(baseline_downstream["AUROC"]) - float(permuted_downstream["AUROC"])
        ),
        "baseline_Brier": baseline_downstream.get("Brier"),
        "permuted_Brier": permuted_downstream.get("Brier"),
        "downstream_Brier_increase": (
            float(permuted_downstream["Brier"]) - float(baseline_downstream["Brier"])
        ),
        "baseline_ECE": baseline_downstream.get("ECE"),
        "permuted_ECE": permuted_downstream.get("ECE"),
        "downstream_ECE_increase": (
            float(permuted_downstream["ECE"]) - float(baseline_downstream["ECE"])
        ),
    }


def _fold_permutation_rows(
    *,
    fold: int,
    seed: int,
    repeats: int,
    test: pd.DataFrame,
    pk_test: pd.DataFrame,
    potency_features: Sequence[str],
    pk_features: Sequence[str],
    potency_permutation_units: Mapping[str, str],
    potency_fit: _StageFit,
    pk_fit: _StageFit,
    baseline_log_ac50: pd.Series,
    baseline_log_cmax: pd.Series,
    baseline_probability: pd.Series,
) -> list[dict[str, Any]]:
    if repeats <= 0:
        return []
    labels = test["_observed_exposure_label"].astype(int)
    baseline_downstream = _binary_metrics(labels, baseline_probability)
    exact = test["_activity_relation"].eq("=")
    baseline_potency_stage = _regression_metrics(
        test.loc[exact, "_log10_ac50_um"],
        baseline_log_ac50.loc[exact],
    )
    baseline_pk_stage = _regression_metrics(
        pk_test["_log10_free_cmax_um"],
        pk_fit.test_prediction,
    )
    rows: list[dict[str, Any]] = []

    for feature in potency_features:
        for repeat in range(repeats):
            rng = _permutation_rng(
                seed=seed,
                fold=fold,
                stage="potency",
                feature=feature,
                repeat=repeat,
            )
            permuted, unit = _permute_potency_feature(
                test,
                feature=feature,
                unit=potency_permutation_units[feature],
                rng=rng,
            )
            design = _transform_design_matrix(permuted, potency_fit.preprocessing)
            log_ac50 = pd.Series(
                potency_fit.model.predict(design),
                index=test.index,
                dtype=float,
            )
            probability = _combine_probability(
                log_ac50,
                baseline_log_cmax,
                ac50_sigma=potency_fit.residual_sigma,
                free_cmax_sigma=pk_fit.residual_sigma,
            )
            rows.append(
                _permutation_record(
                    stage="potency",
                    feature=feature,
                    fold=fold,
                    repeat=repeat,
                    permutation_unit=unit,
                    baseline_stage=baseline_potency_stage,
                    permuted_stage=_regression_metrics(
                        test.loc[exact, "_log10_ac50_um"],
                        log_ac50.loc[exact],
                    ),
                    baseline_downstream=baseline_downstream,
                    permuted_downstream=_binary_metrics(labels, probability),
                    n_downstream=len(test),
                )
            )

    for feature in pk_features:
        for repeat in range(repeats):
            rng = _permutation_rng(
                seed=seed,
                fold=fold,
                stage="free_cmax",
                feature=feature,
                repeat=repeat,
            )
            permuted = pk_test.copy()
            permuted[feature] = rng.permutation(pk_test[feature].to_numpy(copy=True))
            design = _transform_design_matrix(permuted, pk_fit.preprocessing)
            pk_prediction = pd.Series(
                pk_fit.model.predict(design),
                index=pk_test.index,
                dtype=float,
            )
            prediction_by_drug = dict(
                zip(
                    pk_test["drug_id"].astype(str),
                    pk_prediction.astype(float),
                    strict=True,
                )
            )
            log_cmax = test["drug_id"].astype(str).map(prediction_by_drug).astype(float)
            probability = _combine_probability(
                baseline_log_ac50,
                log_cmax,
                ac50_sigma=potency_fit.residual_sigma,
                free_cmax_sigma=pk_fit.residual_sigma,
            )
            rows.append(
                _permutation_record(
                    stage="free_cmax",
                    feature=feature,
                    fold=fold,
                    repeat=repeat,
                    permutation_unit="drug",
                    baseline_stage=baseline_pk_stage,
                    permuted_stage=_regression_metrics(
                        pk_test["_log10_free_cmax_um"],
                        pk_prediction,
                    ),
                    baseline_downstream=baseline_downstream,
                    permuted_downstream=_binary_metrics(labels, probability),
                    n_downstream=len(test),
                )
            )
    return rows


def _summarize_permutations(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "stage",
                "feature",
                "permutation_unit",
                "n_fold_repeats",
                "n_folds",
                "mean_stage_rmse_increase",
                "median_stage_rmse_increase",
                "mean_stage_r2_drop",
                "mean_downstream_AUPRC_drop",
                "median_downstream_AUPRC_drop",
                "q025_downstream_AUPRC_drop",
                "q975_downstream_AUPRC_drop",
                "mean_downstream_AUROC_drop",
                "mean_downstream_Brier_increase",
                "mean_downstream_ECE_increase",
            ]
        )
    grouped = frame.groupby(
        ["stage", "feature", "permutation_unit"],
        dropna=False,
    )
    summary = grouped.agg(
        n_fold_repeats=("downstream_AUPRC_drop", "size"),
        n_folds=("fold", "nunique"),
        mean_stage_rmse_increase=("stage_rmse_increase", "mean"),
        median_stage_rmse_increase=("stage_rmse_increase", "median"),
        mean_stage_r2_drop=("stage_r2_drop", "mean"),
        mean_downstream_AUPRC_drop=("downstream_AUPRC_drop", "mean"),
        median_downstream_AUPRC_drop=("downstream_AUPRC_drop", "median"),
        q025_downstream_AUPRC_drop=(
            "downstream_AUPRC_drop",
            lambda values: values.quantile(0.025),
        ),
        q975_downstream_AUPRC_drop=(
            "downstream_AUPRC_drop",
            lambda values: values.quantile(0.975),
        ),
        mean_downstream_AUROC_drop=("downstream_AUROC_drop", "mean"),
        mean_downstream_Brier_increase=("downstream_Brier_increase", "mean"),
        mean_downstream_ECE_increase=("downstream_ECE_increase", "mean"),
    ).reset_index()
    return summary.sort_values(
        ["stage", "mean_downstream_AUPRC_drop"],
        ascending=[True, False],
    )


def _input_variable_table(
    full: pd.DataFrame,
    usable: pd.DataFrame,
    *,
    label_col: str,
    group_col: str,
    potency_features: Sequence[str],
    pk_features: Sequence[str],
    explicitly_allowed: Sequence[str],
    target_contract: Mapping[str, Any],
) -> pd.DataFrame:
    allowed = set(explicitly_allowed)
    rows: list[dict[str, Any]] = []
    for stage, features in (("potency", potency_features), ("free_cmax", pk_features)):
        for order, feature in enumerate(features, start=1):
            multi = usable.groupby("drug_id", dropna=False)[feature].nunique(
                dropna=True
            )
            rows.append(
                {
                    "stage": stage,
                    "variable": feature,
                    "role": "predictor",
                    "supplied_order": order,
                    "source_column": feature,
                    "dtype": str(full[feature].dtype),
                    "n_missing_full": int(full[feature].isna().sum()),
                    "missing_fraction_full": float(full[feature].isna().mean()),
                    "n_missing_usable": int(usable[feature].isna().sum()),
                    "missing_fraction_usable": float(usable[feature].isna().mean()),
                    "n_unique_usable": int(usable[feature].nunique(dropna=True)),
                    "n_drugs_with_multiple_nonmissing_values": int(multi.gt(1).sum()),
                    "permutation_unit": (
                        "drug"
                        if stage == "free_cmax"
                        else ("row" if multi.gt(1).any() else "drug_block")
                    ),
                    "policy_status": (
                        "explicitly_allowed_label_related_mechanistic_sensitivity_predictor"
                        if feature in allowed
                        else "allowed_predictor"
                    ),
                    "policy_note": (
                        "Label-related mechanistic sensitivity only: fraction-unbound is an observed upstream predictor used in free-Cmax construction. It is not part of the clean structure-only PK model, and inference-time availability must be disclosed."
                        if feature in allowed
                        else "No direct AC50, Cmax/free-Cmax, exposure label, or exposure-margin value is supplied to this stage."
                    ),
                }
            )

    for audit_feature in _AUDIT_ONLY_PROTEIN_BINDING_COLUMNS:
        if audit_feature not in full.columns:
            continue
        multi = usable.groupby("drug_id", dropna=False)[audit_feature].nunique(
            dropna=True
        )
        rows.append(
            {
                "stage": "free_cmax_audit",
                "variable": audit_feature,
                "role": "audit_only_not_predictor",
                "supplied_order": pd.NA,
                "source_column": audit_feature,
                "dtype": str(full[audit_feature].dtype),
                "n_missing_full": int(full[audit_feature].isna().sum()),
                "missing_fraction_full": float(full[audit_feature].isna().mean()),
                "n_missing_usable": int(usable[audit_feature].isna().sum()),
                "missing_fraction_usable": float(usable[audit_feature].isna().mean()),
                "n_unique_usable": int(usable[audit_feature].nunique(dropna=True)),
                "n_drugs_with_multiple_nonmissing_values": int(multi.gt(1).sum()),
                "permutation_unit": "not_permuted",
                "policy_status": "audit_only_deterministic_complement",
                "policy_note": (
                    "Protein-binding percent is the deterministic complement of "
                    "fraction-unbound and is excluded from every predictive design "
                    "matrix to prevent double weighting."
                ),
            }
        )

    target_rows = (
        (
            "potency",
            "_observed_ac50_um",
            target_contract.get("ac50_source_column"),
            "training_target_not_predictor",
        ),
        (
            "free_cmax",
            "_observed_free_cmax_um",
            target_contract.get("free_cmax_source_column"),
            "training_target_not_predictor",
        ),
        ("evaluation", label_col, label_col, "evaluation_label_not_predictor"),
        ("split", group_col, group_col, "split_group_not_predictor"),
        (
            "combiner",
            "predicted_AC50 / predicted_free_Cmax <= 10",
            None,
            "derived_probability_rule",
        ),
    )
    for stage, variable, source, policy in target_rows:
        source_values = (
            full[source] if source in full.columns else pd.Series(dtype=float)
        )
        rows.append(
            {
                "stage": stage,
                "variable": variable,
                "role": policy,
                "supplied_order": pd.NA,
                "source_column": source,
                "dtype": str(source_values.dtype) if source is not None else "derived",
                "n_missing_full": int(source_values.isna().sum())
                if source is not None
                else pd.NA,
                "missing_fraction_full": (
                    float(source_values.isna().mean()) if source is not None else pd.NA
                ),
                "n_missing_usable": pd.NA,
                "missing_fraction_usable": pd.NA,
                "n_unique_usable": pd.NA,
                "n_drugs_with_multiple_nonmissing_values": pd.NA,
                "permutation_unit": "not_applicable",
                "policy_status": policy,
                "policy_note": "This variable is never included in either predictive design matrix.",
            }
        )
    return pd.DataFrame(rows)


def _model_run_record(
    *,
    out: Path,
    dataset_path: Path,
    dataset_version: str,
    label_col: str,
    potency_features: Sequence[str],
    pk_features: Sequence[str],
    group_col: str,
    split_strategy: str,
    n_splits: int,
    model_type: str,
    model_params: Mapping[str, Any],
    seed: int,
    max_threads: int,
    fold_signature_sha256: str,
    censored_policy: str,
    permutation_repeats: int,
    usable: pd.DataFrame,
    per_fold: pd.DataFrame,
    pooled_metrics: Mapping[str, Any],
    top_k_rows: Sequence[Mapping[str, Any]],
    git_commit: str,
) -> dict[str, Any]:
    primary_top_k = top_k_rows[0] if top_k_rows else {}
    record = {
        "run_id": out.name,
        "date": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit,
        "dataset_version": dataset_version,
        "dataset_path": str(dataset_path),
        "label_used": label_col,
        "feature_set": json.dumps(
            {"potency": list(potency_features), "free_cmax": list(pk_features)},
            sort_keys=True,
        ),
        "excluded_columns": "AC50;Cmax/free-Cmax;exposure label;exposure margin;identity/provenance fields",
        "split_method": f"pooled_{n_splits}_fold_{split_strategy}:{group_col}",
        "PU_strategy": "standard_binary",
        "model_type": f"two_stage_{model_type}",
        "selected_model": False,
        "hyperparameters": json.dumps(
            {
                "censored_policy": censored_policy,
                "fold_signature_sha256": fold_signature_sha256,
                "max_threads": int(max_threads),
                "model_params": dict(model_params),
                "permutation_repeats": int(permutation_repeats),
                "residual_sigma_floor_log10": _RESIDUAL_SIGMA_FLOOR,
                "resolved_model_settings_path": str(out / "model_settings.csv"),
                "seed": int(seed),
                "stage_seed_schedule": {
                    "potency": "seed + fold",
                    "free_cmax": "seed + 10000 + fold",
                },
            },
            sort_keys=True,
        ),
        "calibration_method": "analytic_independent_normal_residual_propagation",
        "number_of_rows": int(len(usable)),
        "number_of_positives": int(usable["_observed_exposure_label"].eq(1).sum()),
        "n_train": int(round(float(per_fold["n_train_rows"].median()))),
        "n_test": int(len(usable)),
        "AUROC": pooled_metrics.get("AUROC"),
        "PR_AUC": pooled_metrics.get("AUPRC"),
        "precision_at_K": primary_top_k.get("precision_at_k"),
        "precision_K": primary_top_k.get("effective_k"),
        "enrichment_at_K": primary_top_k.get("enrichment_at_k"),
        "enrichment_K": primary_top_k.get("effective_k"),
        "Brier_score": pooled_metrics.get("Brier"),
        "notes": (
            "Pooled grouped OOF two-stage exposure model. Potency and free-Cmax regressors are fit inside each fold; free-Cmax uses one training row per drug. Fraction-unbound is allowed only as a label-related mechanistic sensitivity predictor; protein-binding percent is audit-only."
        ),
        "model_dir": str(out),
    }
    return {column: record.get(column) for column in LEDGER_COLUMNS}


def run_spd_exposure_grouped_oof(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    label_col: str = "spd_exposure_label",
    potency_features: Sequence[str],
    pk_features: Sequence[str],
    group_col: str = "drug_id",
    n_splits: int = 5,
    model_type: str = "ridge",
    model_params: Mapping[str, Any] | None = None,
    seed: int = 42,
    max_threads: int = 4,
    fold_assignments_path: str | Path | None = None,
    fold_seed: int | None = None,
    censored_policy: str = "exclude",
    permutation_repeats: int = 10,
    top_k: Sequence[int] = _DEFAULT_TOP_K,
) -> dict[str, Any]:
    """Run pooled grouped OOF evaluation for the two-stage SPD exposure model."""
    if model_type not in _SUPPORTED_MODEL_TYPES:
        allowed = ", ".join(sorted(_SUPPORTED_MODEL_TYPES))
        raise ValueError(f"model_type must be one of: {allowed}")
    threads = _validate_max_threads(max_threads)
    requested_model_params = dict(model_params or {})
    _build_regressor(
        model_type,
        seed=seed,
        model_params=requested_model_params,
        max_threads=threads,
    )
    if censored_policy not in {"exclude", "bound"}:
        raise ValueError("censored_policy must be exclude or bound")
    if permutation_repeats < 0:
        raise ValueError("permutation_repeats must be non-negative")
    top_k_values = [int(value) for value in top_k]
    if not top_k_values or any(value <= 0 for value in top_k_values):
        raise ValueError("top_k values must be positive integers")
    top_k_values = list(dict.fromkeys(top_k_values))

    dataset = Path(dataset_path)
    if not dataset.is_file():
        raise FileNotFoundError(f"dataset not found: {dataset}")
    dataset_file_sha256 = _sha256_file(dataset)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    full = pd.read_csv(dataset, low_memory=False)
    potency, pk, explicitly_allowed = _validate_features(
        full,
        label_col=label_col,
        potency_features=potency_features,
        pk_features=pk_features,
    )
    data, target_contract = _prepare_data(
        full,
        label_col=label_col,
        group_col=group_col,
    )
    split_seed = int(seed if fold_seed is None else fold_seed)
    saved_assignments = (
        Path(fold_assignments_path) if fold_assignments_path is not None else None
    )
    if saved_assignments is None:
        folds, split_strategy, effective_splits = _make_folds(
            data,
            group_col=group_col,
            n_splits=n_splits,
            seed=split_seed,
        )
        fold_source = "generated"
    else:
        folds, split_strategy, effective_splits = _load_saved_folds(
            data,
            saved_assignments,
            group_col=group_col,
        )
        if int(n_splits) != effective_splits:
            raise ValueError(
                "saved fold count does not match n_splits: "
                f"saved={effective_splits}, requested={n_splits}"
            )
        fold_source = "saved"
    assignments = _fold_assignment_frame(
        data,
        folds,
        group_col=group_col,
        label_col=label_col,
    )
    fold_signature_sha256 = _fold_signature(assignments)
    potency_permutation_units = {
        feature: (
            "drug_block"
            if feature.startswith("rdkit_")
            or data.groupby("drug_id", dropna=False)[feature]
            .nunique(dropna=False)
            .le(1)
            .all()
            else "row"
        )
        for feature in potency
    }

    predictions: list[pd.DataFrame] = []
    fold_metric_rows: list[dict[str, Any]] = []
    top_k_rows: list[dict[str, Any]] = []
    permutation_rows: list[dict[str, Any]] = []
    model_setting_rows: list[dict[str, Any]] = []

    for fold, (train_idx, test_idx) in enumerate(folds):
        train = data.loc[train_idx].copy()
        test = data.loc[test_idx].copy()
        potency_train = (
            train.loc[train["_activity_relation"].eq("=")].copy()
            if censored_policy == "exclude"
            else train.copy()
        )
        if len(potency_train) < 5:
            raise ValueError(
                f"fold {fold} has only {len(potency_train)} potency-training rows under "
                f"the {censored_policy!r} censoring policy"
            )
        pk_train = _unique_drug_table(train, pk)
        pk_test = _unique_drug_table(test, pk)
        if len(pk_train) < 5 or pk_test.empty:
            raise ValueError(
                f"fold {fold} has insufficient free-Cmax drugs: "
                f"train={len(pk_train)}, test={len(pk_test)}"
            )

        potency_seed = int(seed + fold)
        free_cmax_seed = int(seed + 10_000 + fold)
        potency_fit = _fit_stage(
            potency_train,
            test,
            features=potency,
            target_col="_log10_ac50_um",
            model_type=model_type,
            seed=potency_seed,
            model_params=requested_model_params,
            max_threads=threads,
        )
        pk_fit = _fit_stage(
            pk_train,
            pk_test,
            features=pk,
            target_col="_log10_free_cmax_um",
            model_type=model_type,
            seed=free_cmax_seed,
            model_params=requested_model_params,
            max_threads=threads,
        )
        for stage, stage_seed, stage_fit in (
            ("potency", potency_seed, potency_fit),
            ("free_cmax", free_cmax_seed, pk_fit),
        ):
            model_setting_rows.append(
                {
                    "fold": int(fold),
                    "stage": stage,
                    "model_type": model_type,
                    "model_seed": stage_seed,
                    "max_threads": threads,
                    "fold_signature_sha256": fold_signature_sha256,
                    "requested_model_params_json": json.dumps(
                        requested_model_params,
                        sort_keys=True,
                        default=str,
                    ),
                    "effective_model_params_json": json.dumps(
                        dict(stage_fit.effective_model_params),
                        sort_keys=True,
                        default=str,
                    ),
                }
            )
        pk_prediction_by_drug = dict(
            zip(
                pk_test["drug_id"].astype(str),
                pk_fit.test_prediction.astype(float),
                strict=True,
            )
        )
        log_cmax = test["drug_id"].astype(str).map(pk_prediction_by_drug)
        if log_cmax.isna().any():
            raise AssertionError(
                f"fold {fold} is missing held-out free-Cmax predictions"
            )
        log_cmax = log_cmax.astype(float)
        log_ac50 = potency_fit.test_prediction.astype(float)
        probability = _combine_probability(
            log_ac50,
            log_cmax,
            ac50_sigma=potency_fit.residual_sigma,
            free_cmax_sigma=pk_fit.residual_sigma,
        )

        labels = test["_observed_exposure_label"].astype(int)
        fold_metrics: dict[str, Any] = {
            "fold": int(fold),
            "model_type": model_type,
            "model_seed": int(seed),
            "fold_signature_sha256": fold_signature_sha256,
            "n_train_rows": int(len(train)),
            "n_test_rows": int(len(test)),
            "n_train_groups": int(train[group_col].nunique()),
            "n_test_groups": int(test[group_col].nunique()),
            "n_train_drugs": int(train["drug_id"].nunique()),
            "n_test_drugs": int(test["drug_id"].nunique()),
            "n_potency_train_rows": int(len(potency_train)),
            "n_pk_train_drugs": int(len(pk_train)),
            "n_exact_test_rows": int(test["_activity_relation"].eq("=").sum()),
            "n_test_positive": int(labels.eq(1).sum()),
            "n_test_negative": int(labels.eq(0).sum()),
            "ac50_residual_sigma_log10": potency_fit.residual_sigma,
            "free_cmax_residual_sigma_log10": pk_fit.residual_sigma,
            **_binary_metrics(labels, probability),
        }
        exact = test["_activity_relation"].eq("=")
        potency_stage = _regression_metrics(
            test.loc[exact, "_log10_ac50_um"],
            log_ac50.loc[exact],
        )
        pk_stage = _regression_metrics(
            pk_test["_log10_free_cmax_um"],
            pk_fit.test_prediction,
        )
        fold_metrics.update(
            {f"potency_exact_{key}": value for key, value in potency_stage.items()}
        )
        fold_metrics.update(
            {f"free_cmax_{key}": value for key, value in pk_stage.items()}
        )
        fold_top_k = _top_k_rows(
            labels,
            probability,
            top_k=top_k_values,
            scope="fold",
            fold=fold,
        )
        _add_top_k_columns(fold_metrics, fold_top_k)
        fold_metric_rows.append(fold_metrics)
        top_k_rows.extend(fold_top_k)

        prediction_cols = [column for column in _ID_COLUMNS if column in test.columns]
        prediction = test[prediction_cols].copy()
        prediction.insert(0, "_source_index", test["_source_index"].astype(int))
        prediction["fold"] = int(fold)
        prediction["group_col"] = group_col
        prediction["group_value"] = test[group_col].astype(str)
        prediction["observed_ac50_um"] = test["_observed_ac50_um"].astype(float)
        prediction["observed_ac50_relation"] = test["_activity_relation"].astype(str)
        prediction["predicted_log10_ac50_um"] = log_ac50
        prediction["predicted_ac50_um"] = 10.0**log_ac50
        prediction["observed_free_cmax_um"] = test["_observed_free_cmax_um"].astype(
            float
        )
        prediction["predicted_log10_free_cmax_um"] = log_cmax
        prediction["predicted_free_cmax_um"] = 10.0**log_cmax
        prediction["estimated_exposure_margin"] = (
            prediction["predicted_ac50_um"] / prediction["predicted_free_cmax_um"]
        )
        prediction["estimated_exposure_probability"] = probability
        prediction[label_col] = labels
        predictions.append(prediction)

        permutation_rows.extend(
            _fold_permutation_rows(
                fold=fold,
                seed=seed,
                repeats=permutation_repeats,
                test=test,
                pk_test=pk_test,
                potency_features=potency,
                pk_features=pk,
                potency_permutation_units=potency_permutation_units,
                potency_fit=potency_fit,
                pk_fit=pk_fit,
                baseline_log_ac50=log_ac50,
                baseline_log_cmax=log_cmax,
                baseline_probability=probability,
            )
        )

    pooled = pd.concat(predictions, ignore_index=True).sort_values("_source_index")
    if len(pooled) != len(data) or pooled["_source_index"].nunique() != len(data):
        raise AssertionError(
            "pooled OOF predictions do not cover each usable row exactly once"
        )
    if len(assignments) != len(data) or assignments["_source_index"].nunique() != len(
        data
    ):
        raise AssertionError(
            "fold assignments do not cover each usable row exactly once"
        )

    exact_pooled = pooled["observed_ac50_relation"].eq("=")
    potency_pooled = _regression_metrics(
        np.log10(pooled.loc[exact_pooled, "observed_ac50_um"].astype(float)),
        pooled.loc[exact_pooled, "predicted_log10_ac50_um"],
    )
    unique_pk_pooled = pooled.drop_duplicates("drug_id", keep="first")
    free_cmax_pooled = _regression_metrics(
        np.log10(unique_pk_pooled["observed_free_cmax_um"].astype(float)),
        unique_pk_pooled["predicted_log10_free_cmax_um"],
    )
    pooled_metrics: dict[str, Any] = {
        "scope": "pooled_oof",
        "n_rows": int(len(pooled)),
        "n_positive": int(pooled[label_col].eq(1).sum()),
        "n_negative": int(pooled[label_col].eq(0).sum()),
        "n_drugs": int(pooled["drug_id"].nunique()),
        "n_folds": int(effective_splits),
        "fold_signature_sha256": fold_signature_sha256,
        **{f"potency_exact_{key}": value for key, value in potency_pooled.items()},
        **{f"free_cmax_{key}": value for key, value in free_cmax_pooled.items()},
        **_binary_metrics(
            pooled[label_col],
            pooled["estimated_exposure_probability"],
        ),
    }
    pooled_top_k = _top_k_rows(
        pooled[label_col],
        pooled["estimated_exposure_probability"],
        top_k=top_k_values,
        scope="pooled_oof",
        fold=None,
    )
    _add_top_k_columns(pooled_metrics, pooled_top_k)
    top_k_rows.extend(pooled_top_k)

    per_fold = pd.DataFrame(fold_metric_rows)
    permutation_frame = pd.DataFrame(permutation_rows, columns=_PERMUTATION_COLUMNS)
    permutation_summary = _summarize_permutations(permutation_frame)
    input_variables = _input_variable_table(
        full,
        data,
        label_col=label_col,
        group_col=group_col,
        potency_features=potency,
        pk_features=pk,
        explicitly_allowed=explicitly_allowed,
        target_contract=target_contract,
    )

    output_paths = {
        "pooled_predictions": out / "pooled_oof_predictions.csv",
        "fold_assignments": out / "fold_assignments.csv",
        "pooled_metrics": out / "pooled_metrics.csv",
        "per_fold_metrics": out / "per_fold_metrics.csv",
        "top_k_metrics": out / "top_k_metrics.csv",
        "model_settings": out / "model_settings.csv",
        "stage_permutation_importance": out / "stage_permutation_importance.csv",
        "stage_permutation_importance_summary": out
        / "stage_permutation_importance_summary.csv",
        "input_variables": out / "input_variables.csv",
        "dataset_version_manifest": out / "dataset_version_manifest.json",
        "model_run_record_json": out / "model_run_record.json",
        "model_run_record_csv": out / "model_run_record.csv",
        "manifest": out / "spd_exposure_grouped_oof_manifest.json",
    }
    pooled.to_csv(output_paths["pooled_predictions"], index=False)
    assignments.to_csv(output_paths["fold_assignments"], index=False)
    pd.DataFrame([pooled_metrics]).to_csv(output_paths["pooled_metrics"], index=False)
    per_fold.to_csv(output_paths["per_fold_metrics"], index=False)
    pd.DataFrame(top_k_rows).to_csv(output_paths["top_k_metrics"], index=False)
    pd.DataFrame(model_setting_rows).to_csv(
        output_paths["model_settings"],
        index=False,
    )
    permutation_frame.to_csv(output_paths["stage_permutation_importance"], index=False)
    permutation_summary.to_csv(
        output_paths["stage_permutation_importance_summary"],
        index=False,
    )
    input_variables.to_csv(output_paths["input_variables"], index=False)

    dataset_version = dataset_file_sha256
    git_commit, git_dirty = _git_state()
    dataset_manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset_path": str(dataset),
        "dataset_resolved_path": str(dataset.resolve()),
        "dataset_file_sha256": dataset_file_sha256,
        "dataset_version_id": dataset_version,
        "dataset_file_size_bytes": dataset.stat().st_size,
        "git_commit": git_commit,
        "git_worktree_dirty": git_dirty,
        "label_col": label_col,
        "group_col": group_col,
        "potency_features": potency,
        "pk_features": pk,
        "explicitly_allowed_upstream_predictors": explicitly_allowed,
        "pk_input_policy_class": (
            "label_related_mechanistic_sensitivity"
            if explicitly_allowed
            else "clean_structure_only_pk"
        ),
        "audit_only_upstream_columns": [
            column
            for column in _AUDIT_ONLY_PROTEIN_BINDING_COLUMNS
            if column in full.columns
        ],
        "deterministic_complement_policy": (
            "Protein-binding percent is audit-only; fraction-unbound may be used "
            "alone in the label-related mechanistic sensitivity branch."
        ),
        "target_contract": target_contract,
        "label_counts_full": {
            "positive": int(binary_label_series(full[label_col]).eq(1).sum()),
            "negative": int(binary_label_series(full[label_col]).eq(0).sum()),
            "unknown": int(binary_label_series(full[label_col]).isna().sum()),
        },
        "split": {
            "strategy": split_strategy,
            "source": fold_source,
            "requested_n_splits": int(n_splits),
            "effective_n_splits": int(effective_splits),
            "seed": split_seed,
            "fold_signature_sha256": fold_signature_sha256,
            "saved_assignments_path": (
                str(saved_assignments) if saved_assignments is not None else None
            ),
            "saved_assignments_file_sha256": (
                _sha256_file(saved_assignments)
                if saved_assignments is not None
                else None
            ),
            "every_usable_row_oof_once": True,
            "every_drug_oof_once": True,
            "train_test_group_overlap": 0,
        },
    }
    _write_json(output_paths["dataset_version_manifest"], dataset_manifest)

    record = _model_run_record(
        out=out,
        dataset_path=dataset,
        dataset_version=dataset_version,
        label_col=label_col,
        potency_features=potency,
        pk_features=pk,
        group_col=group_col,
        split_strategy=split_strategy,
        n_splits=effective_splits,
        model_type=model_type,
        model_params=requested_model_params,
        seed=seed,
        max_threads=threads,
        fold_signature_sha256=fold_signature_sha256,
        censored_policy=censored_policy,
        permutation_repeats=permutation_repeats,
        usable=data,
        per_fold=per_fold,
        pooled_metrics=pooled_metrics,
        top_k_rows=pooled_top_k,
        git_commit=git_commit,
    )
    _write_json(output_paths["model_run_record_json"], record)
    pd.DataFrame([record]).reindex(columns=LEDGER_COLUMNS).to_csv(
        output_paths["model_run_record_csv"],
        index=False,
    )

    manifest: dict[str, Any] = {
        "schema_version": 2,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "architecture": "two_stage_spd_exposure_pooled_grouped_oof",
        "dataset": str(dataset),
        "out_dir": str(out),
        "label_col": label_col,
        "group_col": group_col,
        "model_type": model_type,
        "seed": int(seed),
        "model_configuration": {
            "model_type": model_type,
            "requested_model_params": requested_model_params,
            "max_threads": threads,
            "base_seed": int(seed),
            "stage_seed_schedule": {
                "potency": "base_seed + fold",
                "free_cmax": "base_seed + 10000 + fold",
            },
            "resolved_settings_artifact": str(output_paths["model_settings"]),
        },
        "censored_policy": censored_policy,
        "potency_features": potency,
        "pk_features": pk,
        "permutation_repeats": int(permutation_repeats),
        "top_k": top_k_values,
        "split": dataset_manifest["split"],
        "target_contract": target_contract,
        "permutation_policy": {
            "potency_feature_units": potency_permutation_units,
            "free_cmax_feature_unit": "drug",
            "importance_definition": (
                "Held-out raw-variable permutation without refitting; stage error "
                "degradation and downstream exposure-metric changes are both reported."
            ),
        },
        "probability_contract": {
            "event": "predicted_AC50_um / predicted_free_Cmax_um <= 10",
            "formula": "Phi((1 + predicted_log10_free_Cmax - predicted_log10_AC50) / sqrt(sigma_AC50^2 + sigma_free_Cmax^2))",
            "residual_independence_assumption": True,
            "residual_sigma_floor_log10": _RESIDUAL_SIGMA_FLOOR,
            "posthoc_calibration": "none",
        },
        "feature_policy": {
            "direct_response_features_rejected": [
                "AC50 or pIC50",
                "Cmax or free-Cmax",
                "exposure labels/relevance",
                "exposure margins/probabilities",
                "activity-relation/censor fields",
            ],
            "protein_binding_fraction_unbound_policy": (
                "Fraction-unbound may be explicitly supplied only in the label-related mechanistic sensitivity branch. Protein-binding percent is its deterministic complement and remains audit-only; supplying it as a predictor is rejected."
            ),
            "explicitly_allowed_upstream_predictors": explicitly_allowed,
            "audit_only_upstream_columns": [
                column
                for column in _AUDIT_ONLY_PROTEIN_BINDING_COLUMNS
                if column in full.columns
            ],
            "pk_input_policy_class": (
                "label_related_mechanistic_sensitivity"
                if explicitly_allowed
                else "clean_structure_only_pk"
            ),
        },
        "pooled_metrics": pooled_metrics,
        "outputs": {key: str(path) for key, path in output_paths.items()},
        "warnings": [
            "The bound censoring policy treats censor bounds as exact values and is sensitivity-only."
        ]
        if censored_policy == "bound"
        else [],
    }
    _write_json(output_paths["manifest"], manifest)
    return manifest
