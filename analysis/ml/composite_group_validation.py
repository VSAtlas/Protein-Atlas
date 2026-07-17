from __future__ import annotations

import json
from datetime import datetime, timezone
from itertools import combinations, product
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from analysis.ml.dataset_manifest import write_dataset_version_manifest
from analysis.ml.descriptor_constellation_audit import DESCRIPTOR_RESOLUTIONS
from analysis.ml.feature_sets import PHYSICHEM_DESCRIPTOR_FEATURES
from analysis.ml.grouped_cv_ablation_suite import _git_commit
from analysis.ml.grouped_cv_stability import (
    _binary_metrics,
    _split_valid,
    _stratified_group_folds,
    _topk_metrics,
)
from analysis.ml.labels import binary_label_series
from analysis.ml.model_run_ledger import LEDGER_COLUMNS
from analysis.ml.train_classifier_core import (
    _fit_design_matrix,
    _fit_model_with_optional_weights,
    _model,
    _model_probabilities,
    _sample_weights,
    _transform_design_matrix,
)

MISSING_GROUP = "<MISSING>"


def _group_values(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame[column].fillna(MISSING_GROUP).astype(str)


def _practical_descriptor_signature(frame: pd.DataFrame) -> pd.Series:
    parts = []
    for feature in PHYSICHEM_DESCRIPTOR_FEATURES:
        values = pd.to_numeric(frame[feature], errors="coerce")
        rounded = (values / DESCRIPTOR_RESOLUTIONS[feature]).round().astype("Int64")
        parts.append(rounded.astype("string").fillna("<NA>").rename(feature))
    return pd.concat(parts, axis=1).agg("|".join, axis=1)


def _identity_conflicts(frame: pd.DataFrame) -> pd.DataFrame:
    work = pd.DataFrame(
        {
            "drug_id": _group_values(frame, "drug_id"),
            "practical_descriptor_signature": _practical_descriptor_signature(frame),
        },
        index=frame.index,
    )
    counts = (
        work.groupby("drug_id", dropna=False)["practical_descriptor_signature"]
        .nunique(dropna=False)
        .rename("n_practical_signatures")
        .reset_index()
    )
    return counts.loc[counts["n_practical_signatures"].gt(1)].sort_values(
        ["n_practical_signatures", "drug_id"], ascending=[False, True]
    )


def _axis_assignments(
    frame: pd.DataFrame,
    *,
    axis: str,
    label_col: str,
    n_splits: int,
    seed: int,
) -> tuple[pd.Series, int]:
    folds = _stratified_group_folds(
        frame,
        label_col=label_col,
        group_col=axis,
        n_splits=n_splits,
        seed=seed,
    )
    if not folds:
        raise ValueError(f"cannot construct grouped folds for {axis}")
    values = _group_values(frame, axis)
    group_to_fold: dict[str, int] = {}
    for fold, (_, test_idx) in enumerate(folds):
        for value in values.loc[test_idx].unique():
            if value in group_to_fold:
                raise AssertionError(f"{axis} group assigned to multiple folds: {value}")
            group_to_fold[value] = fold
    assigned = values.map(group_to_fold)
    if assigned.isna().any():
        raise AssertionError(f"unassigned {axis} groups remain")
    return assigned.astype(int), len(folds)


def _single_cold_folds(
    frame: pd.DataFrame,
    *,
    axis: str,
    label_col: str,
    n_splits: int,
    seed: int,
) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    assignment, split_count = _axis_assignments(
        frame, axis=axis, label_col=label_col, n_splits=n_splits, seed=seed
    )
    folds = []
    for fold in range(split_count):
        test_mask = assignment.eq(fold)
        folds.append(
            {
                "fold": fold,
                "axis_a_fold": fold,
                "axis_b_fold": pd.NA,
                "train_idx": frame.index[~test_mask],
                "test_idx": frame.index[test_mask],
                "embargo_idx": frame.index[:0],
            }
        )
    assignments = pd.DataFrame(
        {
            "_source_index": frame.index.astype(str),
            "split_axes": axis,
            "axis_a": axis,
            "axis_a_fold": assignment,
            "axis_b": pd.NA,
            "axis_b_fold": pd.NA,
        }
    )
    return folds, assignments


def _double_cold_folds(
    frame: pd.DataFrame,
    *,
    axis_a: str,
    axis_b: str,
    label_col: str,
    n_splits: int,
    seed: int,
) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    assignment_a, count_a = _axis_assignments(
        frame, axis=axis_a, label_col=label_col, n_splits=n_splits, seed=seed
    )
    assignment_b, count_b = _axis_assignments(
        frame, axis=axis_b, label_col=label_col, n_splits=n_splits, seed=seed + 1009
    )
    folds = []
    for composite_fold, (fold_a, fold_b) in enumerate(product(range(count_a), range(count_b))):
        heldout_a = assignment_a.eq(fold_a)
        heldout_b = assignment_b.eq(fold_b)
        test_mask = heldout_a & heldout_b
        train_mask = ~heldout_a & ~heldout_b
        folds.append(
            {
                "fold": composite_fold,
                "axis_a_fold": fold_a,
                "axis_b_fold": fold_b,
                "train_idx": frame.index[train_mask],
                "test_idx": frame.index[test_mask],
                "embargo_idx": frame.index[~train_mask & ~test_mask],
            }
        )
    assignments = pd.DataFrame(
        {
            "_source_index": frame.index.astype(str),
            "split_axes": f"{axis_a}+{axis_b}",
            "axis_a": axis_a,
            "axis_a_fold": assignment_a,
            "axis_b": axis_b,
            "axis_b_fold": assignment_b,
        }
    )
    return folds, assignments


def _fit_fold(
    frame: pd.DataFrame,
    *,
    train_idx: pd.Index,
    test_idx: pd.Index,
    label_col: str,
    features: list[str],
    model_type: str,
    seed: int,
    model_params: Mapping[str, object] | None,
) -> tuple[Any, Any, pd.DataFrame, pd.DataFrame, np.ndarray]:
    train = frame.loc[train_idx].copy()
    test = frame.loc[test_idx].copy()
    x_train, preprocessing = _fit_design_matrix(train, features)
    x_test = _transform_design_matrix(test, preprocessing)
    model = _model(model_type, seed, class_weight="balanced", model_params=model_params)
    weights = _sample_weights(train, label_col, pu_mode="standard_binary")
    _fit_model_with_optional_weights(model, x_train, train[label_col], weights)
    scores = np.asarray(_model_probabilities(model, x_test), dtype=float)
    return model, preprocessing, train, test, scores


def _permutation_map(groups: pd.Series, rng: np.random.Generator) -> dict[str, str]:
    unique = np.asarray(sorted(groups.unique()), dtype=object)
    if len(unique) < 2:
        return {str(value): str(value) for value in unique}
    donors = rng.permutation(unique)
    return {str(recipient): str(donor) for recipient, donor in zip(unique, donors, strict=True)}


def _permute_drug_block(
    test: pd.DataFrame,
    *,
    block: tuple[str, ...],
    donor_map: Mapping[str, str],
) -> pd.DataFrame:
    out = test.copy()
    groups = _group_values(test, "drug_id")
    representatives = test.assign(_drug_group=groups).groupby("_drug_group", dropna=False)[
        list(block)
    ].first()
    donor_groups = groups.map(donor_map)
    for feature in block:
        out[feature] = donor_groups.map(representatives[feature])
    return out


def _permutation_rows(
    *,
    model: Any,
    preprocessing: Any,
    test: pd.DataFrame,
    baseline_scores: np.ndarray,
    label_col: str,
    split_axes: str,
    fold: int,
    repeats: int,
    seed: int,
) -> list[dict[str, Any]]:
    labels = test[label_col].astype(int)
    baseline = _binary_metrics(labels, pd.Series(baseline_scores, index=test.index))
    singles = [(feature,) for feature in PHYSICHEM_DESCRIPTOR_FEATURES]
    pairs = list(combinations(PHYSICHEM_DESCRIPTOR_FEATURES, 2))
    blocks = [*singles, *pairs, tuple(PHYSICHEM_DESCRIPTOR_FEATURES)]
    groups = _group_values(test, "drug_id")
    rows = []
    for repeat in range(repeats):
        rng = np.random.default_rng(seed + fold * 1009 + repeat)
        donor_map = _permutation_map(groups, rng)
        for block in blocks:
            permuted = _permute_drug_block(test, block=block, donor_map=donor_map)
            x_permuted = _transform_design_matrix(permuted, preprocessing)
            scores = pd.Series(_model_probabilities(model, x_permuted), index=test.index)
            metrics = _binary_metrics(labels, scores)
            if len(block) == 1:
                block_type = "single"
            elif len(block) == len(PHYSICHEM_DESCRIPTOR_FEATURES):
                block_type = "all_descriptors"
            else:
                block_type = "pair"
            rows.append(
                {
                    "split_axes": split_axes,
                    "fold": fold,
                    "repeat": repeat,
                    "block_type": block_type,
                    "feature_1": block[0] if block else pd.NA,
                    "feature_2": block[1] if len(block) == 2 else pd.NA,
                    "features": ";".join(block),
                    "n_test": int(len(test)),
                    "n_test_drugs": int(groups.nunique()),
                    "baseline_AUPRC": baseline["AUPRC"],
                    "permuted_AUPRC": metrics["AUPRC"],
                    "AUPRC_drop": baseline["AUPRC"] - metrics["AUPRC"],
                    "baseline_AUROC": baseline["AUROC"],
                    "permuted_AUROC": metrics["AUROC"],
                    "AUROC_drop": baseline["AUROC"] - metrics["AUROC"],
                }
            )
    return rows


def _summarize_permutations(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return rows
    keys = ["split_axes", "block_type", "feature_1", "feature_2", "features"]
    summary = (
        rows.groupby(keys, dropna=False)
        .agg(
            n_fold_repeats=("AUPRC_drop", "size"),
            mean_AUPRC_drop=("AUPRC_drop", "mean"),
            median_AUPRC_drop=("AUPRC_drop", "median"),
            mean_AUROC_drop=("AUROC_drop", "mean"),
            positive_AUPRC_drop_fraction=("AUPRC_drop", lambda values: float((values > 0).mean())),
        )
        .reset_index()
    )
    singles = summary.loc[summary["block_type"].eq("single"), ["split_axes", "feature_1", "mean_AUPRC_drop"]]
    single_a = singles.rename(columns={"feature_1": "feature_1", "mean_AUPRC_drop": "single_1_AUPRC_drop"})
    single_b = singles.rename(columns={"feature_1": "feature_2", "mean_AUPRC_drop": "single_2_AUPRC_drop"})
    summary = summary.merge(single_a, on=["split_axes", "feature_1"], how="left")
    summary = summary.merge(single_b, on=["split_axes", "feature_2"], how="left")
    pair = summary["block_type"].eq("pair")
    summary.loc[pair, "interaction_above_max_single"] = summary.loc[pair, "mean_AUPRC_drop"] - summary.loc[
        pair, ["single_1_AUPRC_drop", "single_2_AUPRC_drop"]
    ].max(axis=1)
    summary.loc[pair, "interaction_above_sum_singles"] = summary.loc[pair, "mean_AUPRC_drop"] - summary.loc[
        pair, ["single_1_AUPRC_drop", "single_2_AUPRC_drop"]
    ].sum(axis=1)
    return summary.sort_values(["split_axes", "mean_AUPRC_drop"], ascending=[True, False])


def _prediction_frame(
    test: pd.DataFrame,
    *,
    scores: np.ndarray,
    label_col: str,
    split_axes: str,
    fold_info: Mapping[str, Any],
    model_type: str,
    feature_variant: str,
    identity_mode: str,
) -> pd.DataFrame:
    meta = [
        "drug_id",
        "target_id",
        "pdb_id",
        "target_family",
        "scaffold_key",
        "chemical_cluster",
        label_col,
        *PHYSICHEM_DESCRIPTOR_FEATURES,
    ]
    out = test[[column for column in meta if column in test.columns]].copy()
    out["ml_prediction_score"] = scores
    out["_source_index"] = test.index.astype(str)
    out["split_axes"] = split_axes
    out["fold"] = int(fold_info["fold"])
    out["axis_a_fold"] = fold_info["axis_a_fold"]
    out["axis_b_fold"] = fold_info["axis_b_fold"]
    out["model_type"] = model_type
    out["feature_variant"] = feature_variant
    out["identity_mode"] = identity_mode
    return out


def _failure_attribution(predictions: pd.DataFrame, *, label_col: str) -> pd.DataFrame:
    rows = []
    keys = ["identity_mode", "split_axes", "model_type", "feature_variant"]
    for key, group in predictions.groupby(keys, dropna=False):
        ranked = group.sort_values("ml_prediction_score", ascending=False).copy()
        k = int(ranked[label_col].eq(1).sum())
        ranked["called_positive"] = False
        if k:
            ranked.iloc[:k, ranked.columns.get_loc("called_positive")] = True
        ranked["error_type"] = "TN"
        ranked.loc[ranked["called_positive"] & ranked[label_col].eq(1), "error_type"] = "TP"
        ranked.loc[ranked["called_positive"] & ranked[label_col].eq(0), "error_type"] = "FP"
        ranked.loc[~ranked["called_positive"] & ranked[label_col].eq(1), "error_type"] = "FN"
        for drug_id, drug in ranked.groupby("drug_id", dropna=False):
            row = {column: value for column, value in zip(keys, key, strict=True)}
            row.update(
                {
                    "drug_id": drug_id,
                    "n": int(len(drug)),
                    "n_positive": int(drug[label_col].eq(1).sum()),
                    "positive_rate": float(drug[label_col].mean()),
                    "TP": int(drug["error_type"].eq("TP").sum()),
                    "FP": int(drug["error_type"].eq("FP").sum()),
                    "FN": int(drug["error_type"].eq("FN").sum()),
                    "TN": int(drug["error_type"].eq("TN").sum()),
                    "mean_score": float(drug["ml_prediction_score"].mean()),
                }
            )
            for feature in PHYSICHEM_DESCRIPTOR_FEATURES:
                row[feature] = pd.to_numeric(drug[feature], errors="coerce").median()
            rows.append(row)
    return pd.DataFrame(rows)


def _pooled_metrics(predictions: pd.DataFrame, *, label_col: str, top_k: int) -> pd.DataFrame:
    rows = []
    keys = ["identity_mode", "split_axes", "model_type", "feature_variant"]
    for key, group in predictions.groupby(keys, dropna=False):
        metrics = _binary_metrics(group[label_col], group["ml_prediction_score"])
        metrics.update(_topk_metrics(group[label_col], group["ml_prediction_score"], k=top_k))
        metrics.update({column: value for column, value in zip(keys, key, strict=True)})
        metrics["n"] = int(len(group))
        metrics["n_positive"] = int(group[label_col].eq(1).sum())
        rows.append(metrics)
    return pd.DataFrame(rows)


def _within_drug_score_audit(predictions: pd.DataFrame) -> pd.DataFrame:
    keys = ["identity_mode", "split_axes", "model_type", "feature_variant", "fold", "drug_id"]
    ranges = (
        predictions.groupby(keys, dropna=False)["ml_prediction_score"]
        .agg(lambda values: float(values.max() - values.min()))
        .rename("within_fold_drug_score_range")
        .reset_index()
    )
    summary_keys = ["identity_mode", "split_axes", "model_type", "feature_variant"]
    return (
        ranges.groupby(summary_keys, dropna=False)
        .agg(
            n_fold_drugs=("drug_id", "size"),
            n_fold_drugs_with_variable_scores=(
                "within_fold_drug_score_range",
                lambda values: int((values > 1e-12).sum()),
            ),
            max_within_fold_drug_score_range=("within_fold_drug_score_range", "max"),
        )
        .reset_index()
    )


def _ledger_rows(
    *,
    metrics: pd.DataFrame,
    integrity: pd.DataFrame,
    dataset: Path,
    dataset_manifest: Mapping[str, Any],
    out_dir: Path,
    run_id: str,
    repo_root: Path,
    seed: int,
    model_params: Mapping[str, object] | None,
) -> pd.DataFrame:
    rows = []
    for row in metrics.itertuples(index=False):
        split_integrity = integrity.loc[
            integrity["identity_mode"].eq(row.identity_mode)
            & integrity["split_axes"].eq(row.split_axes)
        ]
        record: dict[str, Any] = {column: None for column in LEDGER_COLUMNS}
        record.update(
            {
                "run_id": run_id,
                "date": datetime.now(timezone.utc).isoformat(),
                "git_commit": _git_commit(repo_root),
                "dataset_version": dataset_manifest.get("dataset_version_id"),
                "dataset_path": str(dataset),
                "label_used": "spd_binding_label",
                "feature_set": row.feature_variant,
                "excluded_columns": "",
                "split_method": f"pooled_grouped_oof:{row.split_axes}",
                "PU_strategy": "standard_binary",
                "model_type": row.model_type,
                "selected_model": row.model_type,
                "hyperparameters": json.dumps(
                    {"seed": seed, "model_params": dict(model_params or {})}, sort_keys=True
                ),
                "calibration_method": "none",
                "number_of_rows": dataset_manifest.get("n_rows"),
                "number_of_positives": row.n_positive,
                "n_train": int(split_integrity["n_train"].median()),
                "n_test": int(split_integrity["n_test"].median()),
                "AUROC": row.AUROC,
                "PR_AUC": row.AUPRC,
                "precision_at_K": row.precision_at_K,
                "precision_K": row.K,
                "enrichment_at_K": row.enrichment_at_K,
                "enrichment_K": row.K,
                "Brier_score": row.Brier,
                "notes": (
                    "Exploratory RDKit descriptor transfer audit; mapping repair is pending; "
                    f"identity_mode={row.identity_mode}."
                ),
                "model_dir": str(out_dir),
            }
        )
        rows.append(record)
    return pd.DataFrame(rows, columns=LEDGER_COLUMNS)


def run_composite_descriptor_validation(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    label_col: str = "spd_binding_label",
    n_splits: int = 5,
    model_types: Iterable[str] = ("logistic_regression", "lightgbm"),
    feature_variants: Mapping[str, Iterable[str]] | None = None,
    composite_axes: Iterable[tuple[str, str]] = (
        ("drug_id", "target_id"),
        ("drug_id", "scaffold_key"),
        ("chemical_cluster", "target_id"),
    ),
    min_test_positives: int = 10,
    min_test_negatives: int = 10,
    permutation_repeats: int = 5,
    top_k: int = 20,
    seed: int = 42,
    model_params: Mapping[str, object] | None = None,
    run_id: str | None = None,
) -> dict[str, Any]:
    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(dataset, low_memory=False)
    feature_variants = dict(
        feature_variants
        or {
            "rdkit_full": tuple(PHYSICHEM_DESCRIPTOR_FEATURES),
            "rdkit_without_qed": tuple(
                feature for feature in PHYSICHEM_DESCRIPTOR_FEATURES if feature != "rdkit_qed"
            ),
        }
    )
    required = {"drug_id", label_col, *PHYSICHEM_DESCRIPTOR_FEATURES}
    required.update(axis for pair in composite_axes for axis in pair)
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"missing required columns: {', '.join(missing)}")
    raw[label_col] = binary_label_series(raw[label_col])
    labeled = raw.loc[raw[label_col].notna()].copy()
    labeled[label_col] = labeled[label_col].astype(int)
    conflicts = _identity_conflicts(labeled)
    conflicts.to_csv(out / "identity_conflicted_drugs.csv", index=False)
    conflict_ids = set(conflicts["drug_id"].astype(str))

    modes = {
        "identity_consistent_primary": labeled.loc[
            ~_group_values(labeled, "drug_id").isin(conflict_ids)
        ].copy(),
        "all_rows_sensitivity": labeled.copy(),
    }
    prediction_frames = []
    metric_rows = []
    skipped_rows = []
    integrity_rows = []
    assignment_frames = []
    permutation_rows = []

    for identity_mode, frame in modes.items():
        split_specs: list[tuple[str, list[dict[str, Any]], pd.DataFrame]] = []
        cold_folds, cold_assignments = _single_cold_folds(
            frame,
            axis="drug_id",
            label_col=label_col,
            n_splits=n_splits,
            seed=seed,
        )
        split_specs.append(("drug_id", cold_folds, cold_assignments))
        for offset, (axis_a, axis_b) in enumerate(composite_axes, start=1):
            folds, assignments = _double_cold_folds(
                frame,
                axis_a=axis_a,
                axis_b=axis_b,
                label_col=label_col,
                n_splits=n_splits,
                seed=seed + offset * 10007,
            )
            split_specs.append((f"{axis_a}+{axis_b}", folds, assignments))

        for split_axes, folds, assignments in split_specs:
            assignments["identity_mode"] = identity_mode
            assignment_frames.append(assignments)
            axes = split_axes.split("+")
            for fold_info in folds:
                train_idx = fold_info["train_idx"]
                test_idx = fold_info["test_idx"]
                embargo_idx = fold_info["embargo_idx"]
                overlap = {
                    axis: len(
                        set(_group_values(frame.loc[train_idx], axis))
                        & set(_group_values(frame.loc[test_idx], axis))
                    )
                    for axis in axes
                }
                if any(overlap.values()):
                    raise AssertionError(f"group leakage in {split_axes} fold {fold_info['fold']}: {overlap}")
                valid, counts = _split_valid(
                    frame,
                    train_idx,
                    test_idx,
                    label_col=label_col,
                    min_test_positives=min_test_positives,
                    min_test_negatives=min_test_negatives,
                )
                integrity = {
                    "identity_mode": identity_mode,
                    "split_axes": split_axes,
                    "fold": int(fold_info["fold"]),
                    "axis_a_fold": fold_info["axis_a_fold"],
                    "axis_b_fold": fold_info["axis_b_fold"],
                    "n_embargo": int(len(embargo_idx)),
                    "valid": bool(valid),
                    **counts,
                }
                for axis, count in overlap.items():
                    integrity[f"train_test_overlap_{axis}"] = int(count)
                integrity_rows.append(integrity)
                if not valid:
                    skipped_rows.append({**integrity, "reason": "insufficient train/test class counts"})
                    continue
                for feature_variant, variant_features_raw in feature_variants.items():
                    variant_features = list(variant_features_raw)
                    for model_type in model_types:
                        model_seed = seed + int(fold_info["fold"]) + len(split_axes) * 1009
                        model, preprocessing, _, test, scores = _fit_fold(
                            frame,
                            train_idx=train_idx,
                            test_idx=test_idx,
                            label_col=label_col,
                            features=variant_features,
                            model_type=model_type,
                            seed=model_seed,
                            model_params=model_params,
                        )
                        pred = _prediction_frame(
                            test,
                            scores=scores,
                            label_col=label_col,
                            split_axes=split_axes,
                            fold_info=fold_info,
                            model_type=model_type,
                            feature_variant=feature_variant,
                            identity_mode=identity_mode,
                        )
                        prediction_frames.append(pred)
                        row: dict[str, Any] = _binary_metrics(
                            pred[label_col], pred["ml_prediction_score"]
                        )
                        row.update(
                            _topk_metrics(pred[label_col], pred["ml_prediction_score"], k=top_k)
                        )
                        row.update(
                            {
                                "identity_mode": identity_mode,
                                "split_axes": split_axes,
                                "fold": int(fold_info["fold"]),
                                "model_type": model_type,
                                "feature_variant": feature_variant,
                                **counts,
                            }
                        )
                        metric_rows.append(row)
                        if (
                            identity_mode == "identity_consistent_primary"
                            and feature_variant == "rdkit_full"
                            and model_type == "lightgbm"
                            and permutation_repeats > 0
                        ):
                            permutation_rows.extend(
                                _permutation_rows(
                                    model=model,
                                    preprocessing=preprocessing,
                                    test=test,
                                    baseline_scores=scores,
                                    label_col=label_col,
                                    split_axes=split_axes,
                                    fold=int(fold_info["fold"]),
                                    repeats=permutation_repeats,
                                    seed=model_seed,
                                )
                            )

    predictions = pd.concat(prediction_frames, ignore_index=True) if prediction_frames else pd.DataFrame()
    fold_metrics = pd.DataFrame(metric_rows)
    integrity = pd.DataFrame(integrity_rows)
    assignments = pd.concat(assignment_frames, ignore_index=True) if assignment_frames else pd.DataFrame()
    permutations = pd.DataFrame(permutation_rows)
    permutation_summary = _summarize_permutations(permutations)
    pooled = _pooled_metrics(predictions, label_col=label_col, top_k=top_k)
    attribution = _failure_attribution(predictions, label_col=label_col)
    within_drug_scores = _within_drug_score_audit(predictions)

    predictions.to_csv(out / "grouped_cv_oof_predictions.csv", index=False)
    fold_metrics.to_csv(out / "grouped_cv_fold_metrics.csv", index=False)
    pd.DataFrame(skipped_rows).to_csv(out / "grouped_cv_skipped_folds.csv", index=False)
    pooled.to_csv(out / "pooled_oof_metrics.csv", index=False)
    integrity.to_csv(out / "composite_fold_integrity.csv", index=False)
    assignments.to_csv(out / "composite_fold_assignments.csv", index=False)
    permutations.to_csv(out / "grouped_permutation_fold_results.csv", index=False)
    permutation_summary.to_csv(out / "grouped_permutation_summary.csv", index=False)
    attribution.to_csv(out / "descriptor_constellation_failure_attribution.csv", index=False)
    within_drug_scores.to_csv(out / "within_drug_score_variation.csv", index=False)

    dataset_manifest = write_dataset_version_manifest(
        dataset_path=dataset,
        frame=labeled,
        out_path=out / "dataset_version_manifest.json",
        label_col=label_col,
        feature_set="ligand_physchem_descriptors",
        features=list(PHYSICHEM_DESCRIPTOR_FEATURES),
        requested_features=list(PHYSICHEM_DESCRIPTOR_FEATURES),
        provenance={"identity_conflicts_excluded_primary": sorted(conflict_ids)},
    )
    repo_root = Path(__file__).resolve().parents[2]
    resolved_run_id = run_id or out.name
    ledger = _ledger_rows(
        metrics=pooled,
        integrity=integrity,
        dataset=dataset,
        dataset_manifest=dataset_manifest,
        out_dir=out,
        run_id=resolved_run_id,
        repo_root=repo_root,
        seed=seed,
        model_params=model_params,
    )
    ledger.to_csv(out / "ml_model_run_ledger.csv", index=False)
    readiness = {
        "status": "exploratory_blocked",
        "claim": "Descriptor-transfer diagnostic only",
        "blockers": [
            "FDA identity mapping repair is pending",
            "pair-score contract is unresolved; only RDKit descriptors were trained",
            "no independent final holdout was consumed",
        ],
    }
    (out / "model_claim_readiness.json").write_text(
        json.dumps(readiness, indent=2, sort_keys=True), encoding="utf-8"
    )
    manifest = {
        "dataset": str(dataset),
        "dataset_version_id": dataset_manifest.get("dataset_version_id"),
        "label_col": label_col,
        "features": list(PHYSICHEM_DESCRIPTOR_FEATURES),
        "feature_variants": {
            name: list(features) for name, features in feature_variants.items()
        },
        "n_labeled_rows": int(len(labeled)),
        "n_positive": int(labeled[label_col].eq(1).sum()),
        "n_identity_conflicted_drugs": int(len(conflicts)),
        "identity_modes": list(modes),
        "n_splits": int(n_splits),
        "composite_axes": [list(pair) for pair in composite_axes],
        "model_types": list(model_types),
        "permutation_repeats": int(permutation_repeats),
        "seed": int(seed),
        "status": "exploratory",
        "outputs": {
            "predictions": str(out / "grouped_cv_oof_predictions.csv"),
            "fold_metrics": str(out / "grouped_cv_fold_metrics.csv"),
            "pooled_metrics": str(out / "pooled_oof_metrics.csv"),
            "integrity": str(out / "composite_fold_integrity.csv"),
            "assignments": str(out / "composite_fold_assignments.csv"),
            "permutation_summary": str(out / "grouped_permutation_summary.csv"),
            "failure_attribution": str(out / "descriptor_constellation_failure_attribution.csv"),
            "within_drug_score_variation": str(out / "within_drug_score_variation.csv"),
            "ledger": str(out / "ml_model_run_ledger.csv"),
        },
    }
    (out / "grouped_cv_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest
