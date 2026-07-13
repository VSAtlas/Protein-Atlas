from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
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

from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.model_run_ledger import LEDGER_COLUMNS
from analysis.ml.train_classifier_core import (
    _fit_design_matrix,
    _transform_design_matrix,
)


MISSING_GROUP = "<MISSING>"
SUPPORTED_VARIANTS = (
    "pair_only",
    "ligand_only",
    "combined",
    "shortcut_reduced",
)
CONFIG_COLUMNS = ["outer_group", "architecture_variant", "residual_model"]


@dataclass(frozen=True)
class _MatrixTransform:
    preprocessing: Mapping[str, Any]
    mean: pd.Series
    scale: pd.Series

    def transform(self, frame: pd.DataFrame) -> np.ndarray:
        matrix = _transform_design_matrix(frame, self.preprocessing)
        scaled = (matrix - self.mean) / self.scale
        return scaled.fillna(0.0).to_numpy(dtype=float)


@dataclass(frozen=True)
class _PropensityModel:
    matrix: _MatrixTransform
    coefficients: np.ndarray
    constant_logit: float | None
    n_fit_drugs: int

    def predict_logit(self, frame: pd.DataFrame) -> np.ndarray:
        if self.constant_logit is not None:
            return np.full(len(frame), self.constant_logit, dtype=float)
        design = self.matrix.transform(frame)
        return _linear_predictor(design, self.coefficients)


def _unique(values: Iterable[str]) -> list[str]:
    return list(
        dict.fromkeys(str(value).strip() for value in values if str(value).strip())
    )


def _normalize_variants(variants: Iterable[str]) -> list[str]:
    normalized = _unique(value.lower().replace("-", "_") for value in variants)
    unknown = sorted(set(normalized) - set(SUPPORTED_VARIANTS))
    if unknown:
        raise ValueError(f"unsupported architecture variants: {', '.join(unknown)}")
    if not normalized:
        raise ValueError("at least one architecture variant is required")
    return normalized


def _normalize_group_specs(
    group_specs: Sequence[str | Sequence[str]],
) -> list[tuple[str, tuple[str, ...]]]:
    normalized: list[tuple[str, tuple[str, ...]]] = []
    seen: set[tuple[str, ...]] = set()
    for raw in group_specs:
        if isinstance(raw, str):
            columns = tuple(
                _unique(part for token in raw.split("+") for part in token.split(","))
            )
        else:
            columns = tuple(_unique(raw))
        if not columns:
            raise ValueError(f"empty outer group specification: {raw!r}")
        if columns in seen:
            continue
        seen.add(columns)
        normalized.append(("+".join(columns), columns))
    if not normalized:
        raise ValueError("at least one outer group specification is required")
    return normalized


def _group_key(frame: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    values = frame[list(columns)].astype("string").fillna(MISSING_GROUP)
    if len(columns) == 1:
        return values.iloc[:, 0].astype(str)
    return values.apply(
        lambda row: json.dumps(row.tolist(), ensure_ascii=True, separators=(",", ":")),
        axis=1,
    )


def _fit_matrix(
    frame: pd.DataFrame,
    features: Sequence[str],
) -> tuple[np.ndarray, _MatrixTransform]:
    matrix, preprocessing = _fit_design_matrix(frame, list(features))
    mean = matrix.mean(axis=0).fillna(0.0)
    scale = matrix.std(axis=0, ddof=0).replace(0.0, 1.0).fillna(1.0)
    transform = _MatrixTransform(preprocessing=preprocessing, mean=mean, scale=scale)
    scaled = ((matrix - mean) / scale).fillna(0.0)
    return scaled.to_numpy(dtype=float), transform


def _linear_predictor(design: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    return coefficients[0] + design @ coefficients[1:]


def _expit(values: np.ndarray) -> np.ndarray:
    from scipy.special import expit

    return np.asarray(expit(np.clip(values, -40.0, 40.0)), dtype=float)


def _fit_binomial_coefficients(
    design: np.ndarray,
    successes: np.ndarray,
    trials: np.ndarray,
    *,
    l2: float,
) -> np.ndarray:
    from scipy.optimize import minimize
    from scipy.special import expit

    if len(design) == 0 or float(trials.sum()) <= 0:
        raise ValueError("cannot fit propensity model without labeled drug trials")
    rates = successes / trials
    denominator = float(trials.sum())

    def objective(coefficients: np.ndarray) -> tuple[float, np.ndarray]:
        logits = _linear_predictor(design, coefficients)
        residual = trials * (expit(logits) - rates)
        loss = float(np.sum(trials * (np.logaddexp(0.0, logits) - rates * logits)))
        loss = loss / denominator + 0.5 * l2 * float(
            coefficients[1:] @ coefficients[1:]
        )
        gradient = np.empty_like(coefficients)
        gradient[0] = float(residual.sum()) / denominator
        gradient[1:] = design.T @ residual / denominator + l2 * coefficients[1:]
        return loss, gradient

    total_rate = float(successes.sum() / trials.sum())
    initial = np.zeros(design.shape[1] + 1, dtype=float)
    initial[0] = math.log(
        np.clip(total_rate, 1e-6, 1.0 - 1e-6) / np.clip(1.0 - total_rate, 1e-6, 1.0)
    )
    result = minimize(objective, initial, method="L-BFGS-B", jac=True)
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(f"propensity optimization failed: {result.message}")
    return np.asarray(result.x, dtype=float)


def _fit_offset_coefficients(
    design: np.ndarray,
    labels: np.ndarray,
    offsets: np.ndarray,
    *,
    l2: float,
) -> np.ndarray:
    from scipy.optimize import minimize
    from scipy.special import expit

    if len(labels) == 0 or len(np.unique(labels)) < 2:
        raise ValueError("fixed-offset residual training requires both label classes")
    denominator = float(len(labels))

    def objective(coefficients: np.ndarray) -> tuple[float, np.ndarray]:
        logits = offsets + _linear_predictor(design, coefficients)
        residual = expit(logits) - labels
        loss = float(np.sum(np.logaddexp(0.0, logits) - labels * logits))
        loss = loss / denominator + 0.5 * l2 * float(
            coefficients[1:] @ coefficients[1:]
        )
        gradient = np.empty_like(coefficients)
        gradient[0] = float(residual.sum()) / denominator
        gradient[1:] = design.T @ residual / denominator + l2 * coefficients[1:]
        return loss, gradient

    initial = np.zeros(design.shape[1] + 1, dtype=float)
    result = minimize(objective, initial, method="L-BFGS-B", jac=True)
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(
            f"fixed-offset residual optimization failed: {result.message}"
        )
    return np.asarray(result.x, dtype=float)


def _drug_table(
    rows: pd.DataFrame,
    *,
    drug_col: str,
    label_col: str,
    ligand_features: Sequence[str],
) -> pd.DataFrame:
    descriptors = rows.groupby(drug_col, sort=False, dropna=False)[
        list(ligand_features)
    ].first()
    for feature in ligand_features:
        if pd.api.types.is_numeric_dtype(rows[feature]):
            descriptors[feature] = rows.groupby(drug_col, sort=False, dropna=False)[
                feature
            ].median()
    labels = rows.groupby(drug_col, sort=False, dropna=False)[label_col].agg(
        successes="sum",
        trials="size",
    )
    return descriptors.join(labels).reset_index()


def _fit_propensity_model(
    drug_table: pd.DataFrame,
    *,
    ligand_features: Sequence[str],
    l2: float,
) -> _PropensityModel:
    if drug_table.empty:
        raise ValueError("cannot fit propensity model without outer-training drugs")
    design, matrix = _fit_matrix(drug_table, ligand_features)
    successes = drug_table["successes"].to_numpy(dtype=float)
    trials = drug_table["trials"].to_numpy(dtype=float)
    total_successes = float(successes.sum())
    total_trials = float(trials.sum())
    if total_successes <= 0.0 or total_successes >= total_trials:
        smoothed = (total_successes + 0.5) / (total_trials + 1.0)
        constant_logit = math.log(smoothed / (1.0 - smoothed))
        coefficients = np.zeros(design.shape[1] + 1, dtype=float)
        return _PropensityModel(matrix, coefficients, constant_logit, len(drug_table))
    coefficients = _fit_binomial_coefficients(
        design,
        successes,
        trials,
        l2=l2,
    )
    return _PropensityModel(matrix, coefficients, None, len(drug_table))


def _inner_drug_folds(
    drug_table: pd.DataFrame, *, n_splits: int, seed: int
) -> list[np.ndarray]:
    from sklearn.model_selection import KFold, StratifiedKFold

    n_drugs = len(drug_table)
    if n_drugs < 2:
        raise ValueError(
            "inner cross-fitting requires at least two outer-training drugs"
        )
    split_count = min(max(2, int(n_splits)), n_drugs)
    strata = drug_table["successes"].gt(0).astype(int).to_numpy()
    counts = pd.Series(strata).value_counts()
    if len(counts) == 2 and int(counts.min()) >= split_count:
        splitter = StratifiedKFold(
            n_splits=split_count, shuffle=True, random_state=seed
        )
        return [test for _, test in splitter.split(np.zeros(n_drugs), strata)]
    splitter = KFold(n_splits=split_count, shuffle=True, random_state=seed)
    return [test for _, test in splitter.split(np.zeros(n_drugs))]


def _cross_fitted_training_offsets(
    drug_table: pd.DataFrame,
    *,
    drug_col: str,
    ligand_features: Sequence[str],
    n_splits: int,
    seed: int,
    l2: float,
    audit_context: Mapping[str, Any],
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    offsets: dict[str, float] = {}
    audits: list[dict[str, Any]] = []
    for inner_fold, heldout_pos in enumerate(
        _inner_drug_folds(drug_table, n_splits=n_splits, seed=seed)
    ):
        heldout = drug_table.iloc[heldout_pos]
        fit = drug_table.drop(drug_table.index[heldout_pos])
        fit_ids = set(fit[drug_col].astype(str))
        heldout_ids = set(heldout[drug_col].astype(str))
        overlap = fit_ids & heldout_ids
        if overlap:
            raise AssertionError(
                f"inner propensity drug leakage: {sorted(overlap)[:5]}"
            )
        model = _fit_propensity_model(fit, ligand_features=ligand_features, l2=l2)
        logits = model.predict_logit(heldout)
        for drug, logit in zip(heldout[drug_col].astype(str), logits, strict=True):
            offsets[drug] = float(np.clip(logit, -12.0, 12.0))
        audits.append(
            {
                **audit_context,
                "audit_stage": "inner_training_cross_fit",
                "inner_fold": int(inner_fold),
                "test_drug": None,
                "n_fit_drugs": int(len(fit_ids)),
                "n_heldout_drugs": int(len(heldout_ids)),
                "fit_heldout_drug_overlap": int(len(overlap)),
                "exclusion_verified": not overlap,
            }
        )
    expected = set(drug_table[drug_col].astype(str))
    missing = expected - set(offsets)
    if missing:
        raise AssertionError(
            f"inner propensity offsets missing drugs: {sorted(missing)[:5]}"
        )
    return offsets, audits


def _outer_test_offsets(
    train_drugs: pd.DataFrame,
    test_drugs: pd.DataFrame,
    *,
    drug_col: str,
    ligand_features: Sequence[str],
    l2: float,
    audit_context: Mapping[str, Any],
) -> tuple[dict[str, float], dict[str, int], list[dict[str, Any]]]:
    offsets: dict[str, float] = {}
    fit_sizes: dict[str, int] = {}
    audits: list[dict[str, Any]] = []
    train_ids = set(train_drugs[drug_col].astype(str))
    full_model: _PropensityModel | None = None
    for raw_drug in test_drugs[drug_col]:
        drug = str(raw_drug)
        fit = train_drugs.loc[train_drugs[drug_col].astype(str).ne(drug)].copy()
        fit_ids = set(fit[drug_col].astype(str))
        if drug in fit_ids:
            raise AssertionError(
                f"outer-test drug remained in propensity training: {drug}"
            )
        if drug not in train_ids:
            if full_model is None:
                full_model = _fit_propensity_model(
                    train_drugs,
                    ligand_features=ligand_features,
                    l2=l2,
                )
            model = full_model
        else:
            model = _fit_propensity_model(fit, ligand_features=ligand_features, l2=l2)
        descriptor_row = test_drugs.loc[test_drugs[drug_col].astype(str).eq(drug)]
        logit = float(model.predict_logit(descriptor_row)[0])
        offsets[drug] = float(np.clip(logit, -12.0, 12.0))
        fit_sizes[drug] = int(model.n_fit_drugs)
        audits.append(
            {
                **audit_context,
                "audit_stage": "outer_test_drug_exclusion",
                "inner_fold": None,
                "test_drug": drug,
                "n_fit_drugs": int(model.n_fit_drugs),
                "n_heldout_drugs": 1,
                "fit_heldout_drug_overlap": 0,
                "exclusion_verified": drug not in fit_ids,
            }
        )
    return offsets, fit_sizes, audits


def _outer_folds(
    frame: pd.DataFrame,
    *,
    label_col: str,
    group_columns: Sequence[str],
    n_splits: int,
    seed: int,
) -> tuple[list[dict[str, Any]], pd.Series]:
    from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

    if len(group_columns) not in {1, 2}:
        raise ValueError("outer group specifications support one or two columns")

    def axis_assignments(axis: str, axis_seed: int) -> tuple[np.ndarray, int]:
        groups = _group_key(frame, [axis])
        n_groups = int(groups.nunique())
        if n_groups < 2:
            raise ValueError(f"outer axis {axis!r} has fewer than two distinct groups")
        split_count = min(max(2, int(n_splits)), n_groups)
        labels = frame[label_col].to_numpy(dtype=int)
        try:
            splitter = StratifiedGroupKFold(
                n_splits=split_count,
                shuffle=True,
                random_state=axis_seed,
            )
            raw_folds = list(splitter.split(frame, labels, groups))
        except ValueError:
            splitter = GroupKFold(n_splits=split_count)
            raw_folds = list(splitter.split(frame, labels, groups))
        assignments = np.full(len(frame), -1, dtype=int)
        for fold, (_, test_pos) in enumerate(raw_folds):
            assignments[np.asarray(test_pos)] = fold
        if (assignments < 0).any():
            raise AssertionError(f"unassigned outer-axis groups remain for {axis}")
        return assignments, split_count

    if len(group_columns) == 2:
        axis_a, axis_b = group_columns
        assignment_a, count_a = axis_assignments(axis_a, seed)
        assignment_b, count_b = axis_assignments(axis_b, seed + 1009)
        double_cold_folds: list[dict[str, Any]] = []
        for fold_a in range(count_a):
            for fold_b in range(count_b):
                heldout_a = assignment_a == fold_a
                heldout_b = assignment_b == fold_b
                test_mask = heldout_a & heldout_b
                train_mask = ~heldout_a & ~heldout_b
                embargo_mask = ~(test_mask | train_mask)
                double_cold_folds.append(
                    {
                        "train_pos": np.flatnonzero(train_mask),
                        "test_pos": np.flatnonzero(test_mask),
                        "embargo_pos": np.flatnonzero(embargo_mask),
                        "axis_a_fold": int(fold_a),
                        "axis_b_fold": int(fold_b),
                        "split_semantics": "double_cold_cartesian",
                        "test_cell": f"{fold_a}:{fold_b}",
                    }
                )
        row_cells = pd.Series(
            [
                f"{fold_a}:{fold_b}"
                for fold_a, fold_b in zip(assignment_a, assignment_b, strict=True)
            ],
            index=frame.index,
            dtype="string",
        )
        return double_cold_folds, row_cells

    groups = _group_key(frame, group_columns)
    n_groups = int(groups.nunique())
    if n_groups < 2:
        raise ValueError(
            f"outer group {group_columns[0]!r} has fewer than two distinct groups"
        )
    split_count = min(max(2, int(n_splits)), n_groups)
    labels = frame[label_col].to_numpy(dtype=int)
    try:
        splitter = StratifiedGroupKFold(
            n_splits=split_count,
            shuffle=True,
            random_state=seed,
        )
        raw_folds = list(splitter.split(frame, labels, groups))
    except ValueError:
        splitter = GroupKFold(n_splits=split_count)
        raw_folds = list(splitter.split(frame, labels, groups))
    single_cold_folds = [
        {
            "train_pos": np.asarray(train),
            "test_pos": np.asarray(test),
            "embargo_pos": np.asarray([], dtype=int),
            "axis_a_fold": int(fold),
            "axis_b_fold": None,
            "split_semantics": "single_cold",
            "test_cell": str(fold),
        }
        for fold, (train, test) in enumerate(raw_folds)
    ]
    return single_cold_folds, groups


def _fit_residual(
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    features: Sequence[str],
    labels: np.ndarray,
    train_offsets: np.ndarray,
    test_offsets: np.ndarray,
    model_type: str,
    l2: float,
    seed: int,
    model_n_jobs: int,
    lightgbm_params: Mapping[str, Any] | None,
) -> tuple[np.ndarray, int]:
    design_train, transform = _fit_matrix(train, features)
    design_test = transform.transform(test)
    if model_type == "logistic":
        coefficients = _fit_offset_coefficients(
            design_train,
            labels,
            train_offsets,
            l2=l2,
        )
        return _linear_predictor(design_test, coefficients), int(len(coefficients))
    if model_type != "lightgbm":
        raise ValueError(f"unsupported residual model: {model_type}")
    try:
        from lightgbm import LGBMClassifier
    except ImportError as exc:
        raise RuntimeError(
            "lightgbm residual requested but the optional lightgbm package is not installed"
        ) from exc
    params: dict[str, Any] = {
        "objective": "binary",
        "random_state": int(seed),
        "n_estimators": 200,
        "max_depth": 3,
        "num_leaves": 7,
        "learning_rate": 0.05,
        "reg_lambda": float(l2),
        "verbosity": -1,
        "n_jobs": int(model_n_jobs),
        "boost_from_average": False,
    }
    params.update(dict(lightgbm_params or {}))
    model = LGBMClassifier(**params)
    model.fit(design_train, labels, init_score=train_offsets)
    residual = np.asarray(
        model.booster_.predict(design_test, raw_score=True), dtype=float
    )
    return residual, int(model.booster_.num_trees())


def _ece(labels: np.ndarray, probabilities: np.ndarray, *, n_bins: int = 10) -> float:
    bins = np.minimum(
        (np.clip(probabilities, 0.0, 1.0) * n_bins).astype(int), n_bins - 1
    )
    value = 0.0
    for bin_index in range(n_bins):
        mask = bins == bin_index
        if not mask.any():
            continue
        value += float(mask.mean()) * abs(
            float(probabilities[mask].mean()) - float(labels[mask].mean())
        )
    return value


def _metrics(labels: pd.Series, probabilities: pd.Series) -> dict[str, Any]:
    from sklearn.metrics import average_precision_score, roc_auc_score

    y = labels.to_numpy(dtype=int)
    p = probabilities.to_numpy(dtype=float)
    has_both = len(np.unique(y)) == 2
    has_positive = bool((y == 1).any())
    return {
        "n": int(len(y)),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
        "prevalence": float(y.mean()) if len(y) else math.nan,
        "PR_AUC": float(average_precision_score(y, p)) if has_positive else math.nan,
        "AUROC": float(roc_auc_score(y, p)) if has_both else math.nan,
        "Brier": float(np.mean((p - y) ** 2)) if len(y) else math.nan,
        "ECE": _ece(y, p) if len(y) else math.nan,
    }


def _global_topk(
    labels: pd.Series, probabilities: pd.Series, *, k: int
) -> dict[str, Any]:
    ranked = pd.DataFrame({"label": labels, "score": probabilities}).sort_values(
        "score", ascending=False, kind="mergesort"
    )
    top_n = min(max(1, int(k)), len(ranked))
    precision = float(ranked.head(top_n)["label"].mean()) if top_n else math.nan
    prevalence = float(ranked["label"].mean()) if len(ranked) else math.nan
    return {
        "precision_at_K": precision,
        "precision_K": int(top_n),
        "enrichment_at_K": precision / prevalence if prevalence > 0 else math.nan,
        "enrichment_K": int(top_n),
    }


def _pooled_metrics(
    predictions: pd.DataFrame, *, label_col: str, top_k: int
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for key, group in predictions.groupby(CONFIG_COLUMNS, dropna=False):
        row = _metrics(group[label_col], group["diagnostic_probability"])
        row.update(
            _global_topk(group[label_col], group["diagnostic_probability"], k=top_k)
        )
        row.update(dict(zip(CONFIG_COLUMNS, key, strict=True)))
        row["n_unique_source_rows"] = int(group["_source_index"].nunique())
        rows.append(row)
    return pd.DataFrame(rows)


def _within_drug_topk(
    predictions: pd.DataFrame,
    *,
    drug_col: str,
    target_col: str,
    label_col: str,
    top_ks: Sequence[int],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    grouping = [*CONFIG_COLUMNS, drug_col]
    for key, group in predictions.groupby(grouping, dropna=False):
        pairs = (
            group.groupby(target_col, dropna=False)
            .agg(
                pair_label=(label_col, "max"),
                diagnostic_probability=("diagnostic_probability", "max"),
                n_rows=(label_col, "size"),
            )
            .reset_index()
            .sort_values("diagnostic_probability", ascending=False, kind="mergesort")
        )
        n_positive = int(pairs["pair_label"].eq(1).sum())
        base = dict(zip(grouping, key, strict=True))
        for requested_k in top_ks:
            effective_k = min(max(1, int(requested_k)), len(pairs))
            found = int(pairs.head(effective_k)["pair_label"].eq(1).sum())
            rows.append(
                {
                    **base,
                    "requested_K": int(requested_k),
                    "K": int(effective_k),
                    "n_targets": int(len(pairs)),
                    "n_positive_targets": n_positive,
                    "positives_recovered": found,
                    "recovery_at_K": found / n_positive if n_positive else math.nan,
                    "precision_at_K": found / effective_k if effective_k else math.nan,
                }
            )
    return pd.DataFrame(rows)


def _error_summary(
    predictions: pd.DataFrame,
    *,
    group_col: str,
    label_col: str,
) -> pd.DataFrame:
    output_columns = [
        *CONFIG_COLUMNS,
        "error_group_col",
        "error_group_value",
        "n",
        "n_positive",
        "n_negative",
        "prevalence",
        "mean_diagnostic_probability",
        "mean_signed_error",
        "mean_absolute_error",
        "PR_AUC",
        "AUROC",
        "Brier",
        "ECE",
    ]
    if group_col not in predictions.columns:
        return pd.DataFrame(columns=output_columns)
    rows: list[dict[str, Any]] = []
    for key, group in predictions.groupby([*CONFIG_COLUMNS, group_col], dropna=False):
        values = dict(zip([*CONFIG_COLUMNS, group_col], key, strict=True))
        labels = group[label_col].astype(int)
        probability = group["diagnostic_probability"].astype(float)
        row = _metrics(labels, probability)
        row.update(
            {
                **{column: values[column] for column in CONFIG_COLUMNS},
                "error_group_col": group_col,
                "error_group_value": values[group_col],
                "mean_diagnostic_probability": float(probability.mean()),
                "mean_signed_error": float((probability - labels).mean()),
                "mean_absolute_error": float((probability - labels).abs().mean()),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows).reindex(columns=output_columns)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit(repo_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return "unknown"
    value = result.stdout.strip()
    return value if result.returncode == 0 and value else "unknown"


def _ledger_rows(
    pooled: pd.DataFrame,
    *,
    dataset: Path,
    dataset_version: str,
    out_dir: Path,
    run_id: str,
    label_col: str,
    pair_features: Sequence[str],
    shortcut_features: Sequence[str],
    n_rows: int,
    n_positive: int,
    seed: int,
    propensity_l2: float,
    residual_l2: float,
) -> pd.DataFrame:
    repo_root = Path(__file__).resolve().parents[2]
    commit = _git_commit(repo_root)
    rows: list[dict[str, Any]] = []
    for metric in pooled.itertuples(index=False):
        record: dict[str, Any] = {column: None for column in LEDGER_COLUMNS}
        record.update(
            {
                "run_id": run_id,
                "date": datetime.now(timezone.utc).isoformat(),
                "git_commit": commit,
                "dataset_version": dataset_version,
                "dataset_path": str(dataset),
                "label_used": label_col,
                "feature_set": metric.architecture_variant,
                "excluded_columns": ";".join(shortcut_features)
                if metric.architecture_variant == "shortcut_reduced"
                else "",
                "split_method": f"grouped_oof:{metric.outer_group}",
                "PU_strategy": "standard_binary",
                "model_type": (
                    "drug_binomial_logistic"
                    if metric.architecture_variant == "ligand_only"
                    else f"drug_binomial_logistic+fixed_offset_{metric.residual_model}"
                ),
                "selected_model": metric.residual_model,
                "hyperparameters": json.dumps(
                    {
                        "pair_features": list(pair_features),
                        "propensity_l2": propensity_l2,
                        "residual_l2": residual_l2,
                        "seed": seed,
                    },
                    sort_keys=True,
                ),
                "calibration_method": "none_diagnostic_sigmoid_only",
                "number_of_rows": n_rows,
                "number_of_positives": n_positive,
                "n_test": metric.n,
                "AUROC": metric.AUROC,
                "PR_AUC": metric.PR_AUC,
                "precision_at_K": metric.precision_at_K,
                "precision_K": metric.precision_K,
                "enrichment_at_K": metric.enrichment_at_K,
                "enrichment_K": metric.enrichment_K,
                "Brier_score": metric.Brier,
                "notes": (
                    "Exploratory grouped OOF binding-architecture diagnostic only; "
                    "sigmoid outputs are not calibrated probability or publication claims."
                ),
                "model_dir": str(out_dir),
                "ECE": metric.ECE,
                "outer_group": metric.outer_group,
                "architecture_variant": metric.architecture_variant,
            }
        )
        rows.append(record)
    columns = [*LEDGER_COLUMNS, "ECE", "outer_group", "architecture_variant"]
    return pd.DataFrame(rows).reindex(columns=columns)


def run_pair_residual_binding(
    dataset: str | Path,
    out_dir: str | Path,
    *,
    ligand_features: Sequence[str],
    pair_features: Sequence[str],
    outer_group_specs: Sequence[str | Sequence[str]] = ("drug_id",),
    variants: Sequence[str] = ("pair_only", "ligand_only", "combined"),
    shortcut_features: Sequence[str] = (),
    residual_models: Sequence[str] = ("logistic",),
    label_col: str = "spd_binding_label",
    drug_col: str = "drug_id",
    target_col: str = "target_id",
    family_col: str = "target_family",
    source_col: str = "label_source",
    outer_splits: int = 5,
    inner_splits: int = 5,
    top_ks: Sequence[int] = (5, 10, 20),
    propensity_l2: float = 1.0,
    residual_l2: float = 1.0,
    seed: int = 42,
    model_n_jobs: int = 4,
    lightgbm_params: Mapping[str, Any] | None = None,
    run_id: str | None = None,
) -> dict[str, Any]:
    dataset_path = Path(dataset).resolve()
    output = Path(out_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    ligand_features = _unique(ligand_features)
    pair_features = _unique(pair_features)
    shortcut_features = _unique(shortcut_features)
    normalized_variants = _normalize_variants(variants)
    group_specs = _normalize_group_specs(outer_group_specs)
    residual_models = _unique(model.lower() for model in residual_models)
    unknown_models = sorted(set(residual_models) - {"logistic", "lightgbm"})
    if unknown_models:
        raise ValueError(f"unsupported residual models: {', '.join(unknown_models)}")
    if not ligand_features:
        raise ValueError("at least one ligand descriptor feature is required")
    if not pair_features:
        raise ValueError("at least one pair feature is required")
    if not residual_models and any(
        variant != "ligand_only" for variant in normalized_variants
    ):
        raise ValueError(
            "at least one residual model is required for pair-containing variants"
        )
    invalid_shortcuts = sorted(set(shortcut_features) - set(ligand_features))
    if invalid_shortcuts:
        raise ValueError(
            "shortcut-reduced features must be selected ligand features: "
            + ", ".join(invalid_shortcuts)
        )
    if "shortcut_reduced" in normalized_variants and not shortcut_features:
        raise ValueError(
            "shortcut-reduced variant requires at least one --shortcut-feature"
        )
    reduced_ligand_features = [
        feature for feature in ligand_features if feature not in set(shortcut_features)
    ]
    if "shortcut_reduced" in normalized_variants and not reduced_ligand_features:
        raise ValueError("shortcut-reduced variant removed every ligand feature")
    if outer_splits < 2 or inner_splits < 2:
        raise ValueError("outer_splits and inner_splits must both be at least two")
    if propensity_l2 < 0 or residual_l2 < 0:
        raise ValueError("regularization strengths must be non-negative")
    top_ks = sorted({int(value) for value in top_ks if int(value) > 0})
    if not top_ks:
        raise ValueError("at least one positive Top-K value is required")

    raw = pd.read_csv(dataset_path, low_memory=False)
    group_columns = _unique(column for _, columns in group_specs for column in columns)
    required = {
        label_col,
        drug_col,
        target_col,
        *ligand_features,
        *pair_features,
        *group_columns,
    }
    missing = sorted(required - set(raw.columns))
    missing_pair = sorted(set(pair_features) - set(raw.columns))
    if missing_pair:
        raise ValueError(f"missing required pair features: {', '.join(missing_pair)}")
    if missing:
        raise ValueError(f"missing required columns: {', '.join(missing)}")
    all_missing_pair = [
        feature for feature in pair_features if raw[feature].notna().sum() == 0
    ]
    if all_missing_pair:
        raise ValueError(
            f"all-missing required pair features: {', '.join(all_missing_pair)}"
        )
    all_missing_ligand = [
        feature for feature in ligand_features if raw[feature].notna().sum() == 0
    ]
    if all_missing_ligand:
        raise ValueError(
            f"all-missing required ligand features: {', '.join(all_missing_ligand)}"
        )
    assert_no_leakage([*ligand_features, *pair_features])
    raw[label_col] = binary_label_series(raw[label_col])
    excluded_unlabeled_missing_drug_rows = int(
        (raw[label_col].isna() & raw[drug_col].isna()).sum()
    )
    frame = raw.loc[raw[label_col].notna()].copy().reset_index(names="_source_index")
    if frame[drug_col].isna().any():
        raise ValueError(
            f"missing {drug_col} values in labeled rows prevent leak-safe drug cross-fitting"
        )
    frame[label_col] = frame[label_col].astype(int)
    if frame.empty or frame[label_col].nunique() < 2:
        raise ValueError(
            "binding architecture evaluation requires both binary label classes"
        )
    frame[drug_col] = frame[drug_col].astype(str)

    prediction_frames: list[pd.DataFrame] = []
    fold_rows: list[dict[str, Any]] = []
    fold_integrity_rows: list[dict[str, Any]] = []
    skipped_fold_rows: list[dict[str, Any]] = []
    leakage_audits: list[dict[str, Any]] = []
    metadata_columns = _unique(
        [
            "_source_index",
            drug_col,
            target_col,
            family_col,
            source_col,
            *group_columns,
        ]
    )
    metadata_columns = [
        column for column in metadata_columns if column in frame.columns
    ]

    for group_offset, (group_name, columns) in enumerate(group_specs):
        folds, group_keys = _outer_folds(
            frame,
            label_col=label_col,
            group_columns=columns,
            n_splits=outer_splits,
            seed=seed + group_offset * 1009,
        )
        seen_test_positions: set[int] = set()
        for fold, fold_info in enumerate(folds):
            train_pos = np.asarray(fold_info["train_pos"], dtype=int)
            test_pos = np.asarray(fold_info["test_pos"], dtype=int)
            embargo_pos = np.asarray(fold_info["embargo_pos"], dtype=int)
            duplicate_test_positions = seen_test_positions & set(test_pos.tolist())
            if duplicate_test_positions:
                raise AssertionError(
                    f"outer test rows assigned to multiple {group_name} cells: "
                    f"{sorted(duplicate_test_positions)[:5]}"
                )
            seen_test_positions.update(test_pos.tolist())
            train = frame.iloc[train_pos].copy()
            test = frame.iloc[test_pos].copy()
            component_overlap = {
                f"train_test_overlap_{column}": int(
                    len(
                        set(train[column].astype("string").fillna(MISSING_GROUP))
                        & set(test[column].astype("string").fillna(MISSING_GROUP))
                    )
                )
                for column in columns
            }
            if any(component_overlap.values()):
                raise AssertionError(
                    f"outer axis leakage in {group_name}/{fold}: {component_overlap}"
                )
            invalid_reasons: list[str] = []
            if train.empty:
                invalid_reasons.append("empty_train")
            elif train[label_col].nunique() < 2:
                invalid_reasons.append("one_class_train")
            if test.empty:
                invalid_reasons.append("empty_test")
            elif test[label_col].nunique() < 2:
                invalid_reasons.append("one_class_test")
            if train[drug_col].nunique() < 2:
                invalid_reasons.append("fewer_than_two_train_drugs")
            integrity_row = {
                "outer_group": group_name,
                "split_semantics": fold_info["split_semantics"],
                "outer_fold": int(fold),
                "test_cell": fold_info["test_cell"],
                "axis_a": columns[0],
                "axis_a_fold": fold_info["axis_a_fold"],
                "axis_b": columns[1] if len(columns) == 2 else None,
                "axis_b_fold": fold_info["axis_b_fold"],
                "n_train": int(len(train_pos)),
                "n_test": int(len(test_pos)),
                "n_embargo": int(len(embargo_pos)),
                "n_train_positive": int(train[label_col].eq(1).sum()),
                "n_train_negative": int(train[label_col].eq(0).sum()),
                "n_test_positive": int(test[label_col].eq(1).sum()),
                "n_test_negative": int(test[label_col].eq(0).sum()),
                "embargo_source_index_sample": ";".join(
                    frame.iloc[embargo_pos]["_source_index"].astype(str).head(20)
                ),
                "valid": not invalid_reasons,
                "skip_reason": ";".join(invalid_reasons),
                **component_overlap,
            }
            fold_integrity_rows.append(integrity_row)
            if invalid_reasons:
                skipped_fold_rows.append(integrity_row.copy())
                continue
            audit_context = {
                "outer_group": group_name,
                "outer_fold": int(fold),
                "split_semantics": fold_info["split_semantics"],
                "test_cell": fold_info["test_cell"],
                "n_embargo": int(len(embargo_pos)),
            }
            full_propensity_audit_context = {
                **audit_context,
                "propensity_variant": "full_ligand",
            }
            train_drugs = _drug_table(
                train,
                drug_col=drug_col,
                label_col=label_col,
                ligand_features=ligand_features,
            )
            test_drugs = _drug_table(
                test,
                drug_col=drug_col,
                label_col=label_col,
                ligand_features=ligand_features,
            )
            train_offset_map, inner_audit = _cross_fitted_training_offsets(
                train_drugs,
                drug_col=drug_col,
                ligand_features=ligand_features,
                n_splits=inner_splits,
                seed=seed + group_offset * 1009 + fold * 101,
                l2=propensity_l2,
                audit_context=full_propensity_audit_context,
            )
            test_offset_map, test_fit_sizes, test_audit = _outer_test_offsets(
                train_drugs,
                test_drugs,
                drug_col=drug_col,
                ligand_features=ligand_features,
                l2=propensity_l2,
                audit_context=full_propensity_audit_context,
            )
            leakage_audits.extend(inner_audit)
            leakage_audits.extend(test_audit)
            train_propensity = (
                train[drug_col].map(train_offset_map).to_numpy(dtype=float)
            )
            test_propensity = test[drug_col].map(test_offset_map).to_numpy(dtype=float)
            test_fit_n = test[drug_col].map(test_fit_sizes).to_numpy(dtype=int)
            if (
                not np.isfinite(train_propensity).all()
                or not np.isfinite(test_propensity).all()
            ):
                raise AssertionError("non-finite cross-fitted propensity offsets")

            shortcut_train_propensity = train_propensity
            shortcut_test_propensity = test_propensity
            shortcut_fit_n = test_fit_n
            if "shortcut_reduced" in normalized_variants:
                reduced_audit_context = {
                    **audit_context,
                    "propensity_variant": "reduced_ligand",
                }
                reduced_train_map, reduced_inner_audit = _cross_fitted_training_offsets(
                    train_drugs,
                    drug_col=drug_col,
                    ligand_features=reduced_ligand_features,
                    n_splits=inner_splits,
                    seed=seed + group_offset * 1009 + fold * 101 + 17,
                    l2=propensity_l2,
                    audit_context=reduced_audit_context,
                )
                reduced_test_map, reduced_fit_sizes, reduced_test_audit = (
                    _outer_test_offsets(
                        train_drugs,
                        test_drugs,
                        drug_col=drug_col,
                        ligand_features=reduced_ligand_features,
                        l2=propensity_l2,
                        audit_context=reduced_audit_context,
                    )
                )
                leakage_audits.extend(reduced_inner_audit)
                leakage_audits.extend(reduced_test_audit)
                shortcut_train_propensity = (
                    train[drug_col].map(reduced_train_map).to_numpy(dtype=float)
                )
                shortcut_test_propensity = (
                    test[drug_col].map(reduced_test_map).to_numpy(dtype=float)
                )
                shortcut_fit_n = (
                    test[drug_col].map(reduced_fit_sizes).to_numpy(dtype=int)
                )
                if (
                    not np.isfinite(shortcut_train_propensity).all()
                    or not np.isfinite(shortcut_test_propensity).all()
                ):
                    raise AssertionError(
                        "non-finite shortcut-reduced propensity offsets"
                    )

            for variant in normalized_variants:
                if variant == "ligand_only":
                    model_names = ["none"]
                    selected_pair_features: list[str] = []
                else:
                    model_names = residual_models
                    selected_pair_features = pair_features
                if variant == "shortcut_reduced":
                    variant_train_propensity = shortcut_train_propensity
                    variant_test_propensity = shortcut_test_propensity
                    variant_fit_n = shortcut_fit_n
                    propensity_feature_variant = "reduced_ligand"
                else:
                    variant_train_propensity = train_propensity
                    variant_test_propensity = test_propensity
                    variant_fit_n = test_fit_n
                    propensity_feature_variant = "full_ligand"
                for residual_model in model_names:
                    if variant == "pair_only":
                        train_offsets = np.zeros(len(train), dtype=float)
                        test_offsets = np.zeros(len(test), dtype=float)
                    else:
                        train_offsets = variant_train_propensity
                        test_offsets = variant_test_propensity
                    if variant == "ligand_only":
                        residual = np.zeros(len(test), dtype=float)
                        n_parameters = 0
                    else:
                        residual, n_parameters = _fit_residual(
                            train,
                            test,
                            features=selected_pair_features,
                            labels=train[label_col].to_numpy(dtype=int),
                            train_offsets=train_offsets,
                            test_offsets=test_offsets,
                            model_type=residual_model,
                            l2=residual_l2,
                            seed=seed + group_offset * 1009 + fold * 101,
                            model_n_jobs=model_n_jobs,
                            lightgbm_params=lightgbm_params,
                        )
                    final_logit = test_offsets + residual
                    probability = _expit(final_logit)
                    pred = test[metadata_columns + [label_col]].copy()
                    pred["outer_group"] = group_name
                    pred["outer_group_key"] = group_keys.iloc[test_pos].to_numpy()
                    pred["outer_split_semantics"] = fold_info["split_semantics"]
                    pred["outer_test_cell"] = fold_info["test_cell"]
                    pred["outer_embargo_count"] = int(len(embargo_pos))
                    pred["outer_fold"] = int(fold)
                    pred["architecture_variant"] = variant
                    pred["residual_model"] = residual_model
                    pred["propensity_feature_variant"] = propensity_feature_variant
                    pred["propensity_logit"] = variant_test_propensity
                    pred["applied_offset_logit"] = test_offsets
                    pred["residual_logit"] = residual
                    pred["final_logit"] = final_logit
                    pred["diagnostic_probability"] = probability
                    pred["propensity_fit_n_drugs"] = variant_fit_n
                    pred["propensity_excluded_test_drug"] = True
                    prediction_frames.append(pred)

                    fold_metric = _metrics(
                        pred[label_col], pred["diagnostic_probability"]
                    )
                    fold_rows.append(
                        {
                            **fold_metric,
                            "outer_group": group_name,
                            "split_semantics": fold_info["split_semantics"],
                            "test_cell": fold_info["test_cell"],
                            "outer_fold": int(fold),
                            "architecture_variant": variant,
                            "residual_model": residual_model,
                            "n_train": int(len(train)),
                            "n_test": int(len(test)),
                            "n_embargo": int(len(embargo_pos)),
                            "n_train_drugs": int(train[drug_col].nunique()),
                            "n_test_drugs": int(test[drug_col].nunique()),
                            "outer_group_overlap": 0,
                            "n_residual_parameters_or_trees": int(n_parameters),
                            **component_overlap,
                        }
                    )

    if prediction_frames:
        predictions = pd.concat(prediction_frames, ignore_index=True)
    else:
        prediction_columns = _unique(
            [
                *metadata_columns,
                label_col,
                *CONFIG_COLUMNS,
                "diagnostic_probability",
            ]
        )
        predictions = pd.DataFrame(columns=prediction_columns)
    duplicate_keys = [*CONFIG_COLUMNS, "_source_index"]
    if not predictions.empty and predictions.duplicated(duplicate_keys).any():
        raise AssertionError("pooled OOF rows are duplicated within an architecture")
    fold_metrics = pd.DataFrame(fold_rows)
    fold_integrity = pd.DataFrame(fold_integrity_rows)
    skipped_folds = pd.DataFrame(
        skipped_fold_rows,
        columns=fold_integrity.columns,
    )
    pooled = _pooled_metrics(predictions, label_col=label_col, top_k=top_ks[0])
    within_drug = _within_drug_topk(
        predictions,
        drug_col=drug_col,
        target_col=target_col,
        label_col=label_col,
        top_ks=top_ks,
    )
    error_outputs = {
        "per_drug_errors": _error_summary(
            predictions, group_col=drug_col, label_col=label_col
        ),
        "per_target_errors": _error_summary(
            predictions, group_col=target_col, label_col=label_col
        ),
        "per_family_errors": _error_summary(
            predictions, group_col=family_col, label_col=label_col
        ),
        "per_source_errors": _error_summary(
            predictions, group_col=source_col, label_col=label_col
        ),
    }
    dataset_version = _file_sha256(dataset_path)
    resolved_run_id = run_id or output.name
    ledger = _ledger_rows(
        pooled,
        dataset=dataset_path,
        dataset_version=dataset_version,
        out_dir=output,
        run_id=resolved_run_id,
        label_col=label_col,
        pair_features=pair_features,
        shortcut_features=shortcut_features,
        n_rows=len(frame),
        n_positive=int(frame[label_col].eq(1).sum()),
        seed=seed,
        propensity_l2=propensity_l2,
        residual_l2=residual_l2,
    )

    output_paths = {
        "predictions": output / "pair_residual_oof_predictions.csv",
        "fold_metrics": output / "pair_residual_fold_metrics.csv",
        "pooled_metrics": output / "pair_residual_pooled_metrics.csv",
        "within_drug_topk": output / "within_drug_topk_recovery.csv",
        "propensity_leakage_audit": output / "propensity_leakage_audit.csv",
        "fold_integrity": output / "outer_fold_integrity.csv",
        "skipped_folds": output / "outer_skipped_folds.csv",
        "ledger": output / "ml_model_run_ledger.csv",
        **{name: output / f"{name}.csv" for name in error_outputs},
    }
    predictions.to_csv(output_paths["predictions"], index=False)
    fold_metrics.to_csv(output_paths["fold_metrics"], index=False)
    pooled.to_csv(output_paths["pooled_metrics"], index=False)
    within_drug.to_csv(output_paths["within_drug_topk"], index=False)
    pd.DataFrame(leakage_audits).to_csv(
        output_paths["propensity_leakage_audit"], index=False
    )
    fold_integrity.to_csv(output_paths["fold_integrity"], index=False)
    skipped_folds.to_csv(output_paths["skipped_folds"], index=False)
    ledger.to_csv(output_paths["ledger"], index=False)
    for name, errors in error_outputs.items():
        errors.to_csv(output_paths[name], index=False)

    readiness = {
        "status": "exploratory_blocked",
        "probability_claim": False,
        "publication_claim": False,
        "reason": (
            "Cross-fitted logits are converted with a sigmoid only for Brier/ECE diagnostics; "
            "they are not independently calibrated probabilities and no final holdout was consumed."
        ),
    }
    readiness_path = output / "model_claim_readiness.json"
    readiness_path.write_text(
        json.dumps(readiness, indent=2, sort_keys=True), encoding="utf-8"
    )
    descriptor_conflicts = {
        feature: int(
            frame.groupby(drug_col, dropna=False)[feature]
            .nunique(dropna=False)
            .gt(1)
            .sum()
        )
        for feature in ligand_features
    }
    manifest: dict[str, Any] = {
        "status": "exploratory",
        "claim_status": "no_probability_or_publication_claim",
        "run_id": resolved_run_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": str(dataset_path),
        "dataset_sha256": dataset_version,
        "label_col": label_col,
        "n_labeled_rows": int(len(frame)),
        "n_positive": int(frame[label_col].eq(1).sum()),
        "excluded_unlabeled_missing_drug_rows": excluded_unlabeled_missing_drug_rows,
        "drug_col": drug_col,
        "target_col": target_col,
        "family_col": family_col,
        "source_col": source_col,
        "ligand_features": ligand_features,
        "pair_features": pair_features,
        "shortcut_features_removed": shortcut_features,
        "shortcut_ligand_features_removed": shortcut_features,
        "shortcut_reduced_ligand_features": reduced_ligand_features,
        "ligand_descriptor_conflicted_drug_counts": descriptor_conflicts,
        "outer_split_specs": [
            {
                "name": name,
                "axes": list(columns),
                "semantics": (
                    "single_cold" if len(columns) == 1 else "double_cold_cartesian"
                ),
                "test_rule": (
                    "held-out axis groups"
                    if len(columns) == 1
                    else "rows held out on both independently assigned axes"
                ),
                "train_rule": (
                    "all non-held-out axis groups"
                    if len(columns) == 1
                    else "rows held out on neither axis"
                ),
                "embargo_rule": (
                    "none" if len(columns) == 1 else "rows held out on exactly one axis"
                ),
            }
            for name, columns in group_specs
        ],
        "outer_splits_requested": int(outer_splits),
        "inner_splits_requested": int(inner_splits),
        "variants": normalized_variants,
        "residual_models": residual_models,
        "top_ks": top_ks,
        "propensity_l2": float(propensity_l2),
        "residual_l2": float(residual_l2),
        "seed": int(seed),
        "leakage_controls": {
            "outer_train_propensities": "inner cross-fitted with whole drugs held out",
            "outer_test_propensities": (
                "one fit per test drug using outer-train drugs after removing that drug entirely"
            ),
            "stage_two_offsets": "fixed; stage-two training does not refit propensity logits",
            "all_exclusions_verified": all(
                bool(row["exclusion_verified"]) for row in leakage_audits
            ),
        },
        "diagnostic_probability_warning": readiness["reason"],
        "outputs": {name: str(path) for name, path in output_paths.items()},
        "claim_readiness": str(readiness_path),
    }
    manifest_path = output / "architecture_manifest.json"
    manifest["outputs"]["architecture_manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest
