from __future__ import annotations

import importlib
import json
import math
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.api import types as ptypes

from analysis.ml.audit_utils import label_balance, load_table
from analysis.ml.data_gap_report import run_data_gap_report
from analysis.ml.dataset_eda import OptionalDependencyMissing, run_dataset_eda
from analysis.ml.feature_sets import FORBIDDEN_FEATURES, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.score_scale_audit import write_score_scale_audit
from analysis.ml.splits import make_split
from analysis.ml.train_classifier_core import (
    _fit_design_matrix,
    _model,
    _model_probabilities,
    _transform_design_matrix,
)


DEFAULT_LABEL_COLS = [
    "spd_binding_label",
    "spd_exposure_relevant",
    "combined_activity_ml_label",
    "external_four_state_ml_label",
    "tissue_site_label",
]

DEFAULT_GROUP_COLS = [
    "target_family",
    "protein_class",
    "scaffold_key",
    "ligand_chemotype",
    "pdb_id",
    "drug_id",
]

DEFAULT_EXTENDED_GROUP_COLS = [
    *DEFAULT_GROUP_COLS,
    "chemical_cluster",
    "ligand_chemotype_broad",
    "label_source",
    "source_family",
    "upstream_source",
    "assay_type",
    "endpoint_type",
    "activity_type",
    "source_label_policy",
]

PROVENANCE_COLUMNS: Mapping[str, list[str]] = {
    "identity": [
        "drug_id",
        "target_id",
        "pdb_id",
        "target_gene",
        "target_uniprot",
        "canonical_pair_key",
    ],
    "assay_source": [
        "label_source",
        "source_family",
        "upstream_source",
        "source_label_policy",
        "assay_type",
        "endpoint_type",
        "activity_type",
        "database_release_year",
        "availability_year",
    ],
    "target_context": [
        "target_family",
        "target_family_source",
        "protein_class",
        "target_adr_evidence_count",
        "ot_safety_liability_count",
    ],
    "ligand_context": [
        "scaffold_key",
        "scaffold_source",
        "chemical_cluster",
        "chemical_cluster_source",
        "ligand_chemotype",
        "ligand_chemotype_source",
        "ligand_chemotype_broad",
        "rdkit_descriptor_source",
        "rdkit_descriptor_status",
    ],
    "structure_score_context": [
        "structure_quality",
        "structure_quality_source",
        "pdb_resolution",
        "pdb_rank_score",
        "pdb_selection_method",
        "banana_score_feature_source",
        "banana_binding_probability_feature_source",
        "SCORCH_score_used_feature_source",
        "atlas_score",
        "atlas_score_source_for_ml",
        "final_score",
        "final_score_source",
        "z_selected",
        "z_selected_source",
        "z_selected_feature_source",
        "z_stage1",
        "z_stage2",
        "consensus_z_score",
        "consensus_z_score_source",
        "z_vs_decoys_consensus",
        "z_vs_compare_run_consensus",
        "z_vs_decoys_blend",
        "consensus_z_decoy_n",
        "consensus_z_decoy_unique",
        "consensus_z_decoy_zero_fraction",
        "consensus_z_decoy_mu",
        "consensus_z_decoy_sigma",
    ],
    "external_evidence": [
        "external_four_state_label_status",
        "external_evidence_sources",
        "external_parent_sources",
        "external_n_raw_evidence_rows",
        "external_n_production_evidence_rows",
        "external_n_benchmark_only_rows",
        "external_max_confidence",
        "external_source_family",
        "external_upstream_source",
        "combined_activity_label_status",
        "mechanism_label_status",
        "mechanism_label_source",
    ],
    "tissue_site_context": [
        "adr_site_group",
        "adr_site_mapping_source",
        "adr_site_mapping_confidence",
        "site_name",
        "tissue_site_label_source",
        "tissue_site_label_evidence_type",
        "tissue_site_label_confidence",
    ],
}

DEFAULT_LABEL_FEATURE_SETS = {
    "spd_binding_label": "spd_binding_nonleaky",
    "spd_exposure_relevant": "spd_exposure_nonleaky",
    "combined_activity_ml_label": "spd_binding_nonleaky",
    "external_four_state_ml_label": "flat_mechanism_baseline",
    "tissue_site_label": "tissue_site_score_only",
}

HIGH_RISK_AUDIT_COLUMN_TOKENS = [
    "_positive_rate",
    "_assay_density",
    "_labelable_rate",
    "label_source_positive_rate",
    "source_family_positive_rate",
]


def run_provenance_profile_audit(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    labels: Sequence[str] | None = None,
    group_cols: Sequence[str] | None = None,
    profile_label: str | None = None,
    run_profile: bool = True,
    run_associations: bool = True,
    run_permutation_importance: bool = True,
    run_drift: bool = True,
    run_validation: bool = True,
    write_augmented_dataset: bool = True,
    sample_size: int = 10_000,
    profile_sample_size: int = 5_000,
    max_association_columns: int = 80,
    max_category_levels: int = 200,
    max_mi_features: int = 300,
    permutation_sample_size: int = 5_000,
    permutation_repeats: int = 5,
    permutation_model: str = "fast_logistic",
    drift_splits: Sequence[str] | None = None,
    seed: int = 42,
    test_fraction: float = 0.2,
    fail_on_missing: bool = False,
) -> dict[str, Any]:
    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    df = load_table(dataset)
    active_labels = _available(labels or DEFAULT_LABEL_COLS, df)
    active_groups = _available(group_cols or DEFAULT_EXTENDED_GROUP_COLS, df)
    primary_label = profile_label if profile_label in df.columns else (active_labels[0] if active_labels else None)

    outputs: dict[str, str] = {}
    skipped: dict[str, str] = {}
    errors: dict[str, str] = {}
    optional_status = _optional_dependency_status(
        ["ydata_profiling", "phik", "dython", "evidently", "pandera", "great_expectations", "sklearn"]
    )
    optional_status.to_csv(out / "optional_dependency_status.csv", index=False)
    outputs["optional_dependency_status"] = str(out / "optional_dependency_status.csv")

    outputs["provenance_columns_present"] = _write_provenance_column_summary(df, out / "provenance_columns_present.csv")
    outputs["provenance_feature_summary"] = _write_provenance_feature_summary(df, out / "provenance_feature_summary.csv")
    selected_audit_features: list[str] = []
    for label in active_labels:
        feature_set = DEFAULT_LABEL_FEATURE_SETS.get(label)
        if feature_set:
            selected_audit_features.extend(get_feature_set(feature_set))
    score_scale = write_score_scale_audit(
        df,
        out / "score_scale",
        feature_names=selected_audit_features,
    )
    outputs.update(
        {f"score_scale_{name}": path for name, path in score_scale["outputs"].items()}
    )
    outputs["all_label_summary"] = _write_all_label_summaries(df, active_labels, out / "label_summary_by_label.csv")
    outputs["label_group_positive_rates"] = _write_group_positive_rates(
        df,
        active_labels,
        active_groups,
        out / "label_group_positive_rates.csv",
    )
    _run_step(
        "data_gap_report",
        outputs,
        skipped,
        errors,
        fail_on_missing,
        lambda: {
            f"data_gap_{key}": value
            for key, value in run_data_gap_report(
                dataset,
                out / "data_gap_report",
                labels=active_labels,
                group_cols=active_groups,
            ).get("outputs", {}).items()
        },
    )

    augmented_path: Path | None = None
    if write_augmented_dataset:
        augmented = add_audit_provenance_columns(df, active_labels, active_groups)
        augmented_path = out / f"{dataset.stem}.audit_provenance.csv"
        augmented.to_csv(augmented_path, index=False)
        outputs["audit_augmented_dataset"] = str(augmented_path)
        outputs["audit_augmented_column_dictionary"] = _write_audit_column_dictionary(
            augmented,
            df.columns,
            out / "audit_augmented_column_dictionary.csv",
        )

    if run_associations:
        outputs.update(
            _run_multi_label_mutual_information(
                df,
                active_labels,
                out,
                max_features=max_mi_features,
                max_category_levels=max_category_levels,
                sample_size=sample_size,
                seed=seed,
            )
        )

    if run_permutation_importance:
        _run_step(
            "permutation_importance",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: _run_permutation_importance_suite(
                df,
                active_labels,
                out,
                sample_size=permutation_sample_size,
                repeats=permutation_repeats,
                model_type=permutation_model,
                seed=seed,
                test_fraction=test_fraction,
            ),
        )

    if run_validation:
        _run_step(
            "pandera_validation",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: _run_pandera_validation(df, active_labels, out),
        )

    if run_drift:
        _run_step(
            "evidently_drift",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: _run_evidently_drift_suite(
                df,
                out,
                labels=active_labels,
                splits=drift_splits or ["target_holdout", "chemical_cluster_holdout", "target_family_holdout"],
                sample_size=sample_size,
                seed=seed,
                test_fraction=test_fraction,
            ),
        )

    if run_profile:
        profile_tools = ["profile"]
        if run_associations:
            profile_tools.extend(["phik", "dython", "mutual-info", "network"])
        profile_out = out / "dataset_eda"
        _run_step(
            "dataset_eda",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: {
                f"eda_{name}": value
                for name, value in run_dataset_eda(
                    dataset,
                    profile_out,
                    tools=profile_tools,
                    label_col=primary_label,
                    sample_rows=profile_sample_size,
                    random_state=seed,
                    max_association_columns=max_association_columns,
                    max_category_levels=max_category_levels,
                    max_mi_features=max_mi_features,
                    ydata_minimal=True,
                    include_id_like=False,
                    include_high_cardinality=False,
                    fail_on_missing=fail_on_missing,
                ).get("outputs", {}).items()
            },
        )

    manifest = {
        "dataset": str(dataset),
        "out_dir": str(out),
        "n_rows": int(len(df)),
        "n_columns": int(df.shape[1]),
        "labels": active_labels,
        "group_cols": active_groups,
        "profile_label": primary_label,
        "policy": {
            "audit_columns_are_not_clean_predictive_features": True,
            "positive_rate_and_assay_density_columns_are_label_derived": True,
            "clean_models_should_continue_to_use_feature_sets_from_analysis_ml_feature_sets": True,
        },
        "high_risk_audit_column_tokens": HIGH_RISK_AUDIT_COLUMN_TOKENS,
        "score_scale": score_scale,
        "outputs": outputs,
        "skipped": skipped,
        "errors": errors,
    }
    (out / "provenance_profile_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    _write_markdown_summary(out / "provenance_profile_audit.md", manifest)
    return manifest


def add_audit_provenance_columns(
    df: pd.DataFrame,
    labels: Sequence[str],
    group_cols: Sequence[str],
) -> pd.DataFrame:
    derived_columns: dict[str, pd.Series] = {}
    for label in labels:
        label_series = binary_label_series(df[label])
        label_token = _safe_column_token(label)
        for group_col in group_cols:
            if group_col not in df.columns:
                continue
            group = df[group_col].fillna("<missing>").astype(str)
            stats = _group_label_stats(group, label_series)
            density_col = f"_audit_{label_token}_{_safe_column_token(group_col)}_assay_density"
            rate_col = f"_audit_{label_token}_{_safe_column_token(group_col)}_positive_rate"
            labelable_col = f"_audit_{label_token}_{_safe_column_token(group_col)}_labelable_rate"
            derived_columns[density_col] = group.map(stats["labelable"]).fillna(0).astype(int)
            derived_columns[rate_col] = group.map(stats["positive_rate"])
            derived_columns[labelable_col] = group.map(stats["labelable_rate"])
    if not derived_columns:
        return df.copy()
    return pd.concat([df.copy(), pd.DataFrame(derived_columns, index=df.index)], axis=1)


def _group_label_stats(group: pd.Series, label_series: pd.Series) -> pd.DataFrame:
    data = pd.DataFrame({"group": group, "label": label_series})
    rows: list[dict[str, Any]] = []
    total_counts = data.groupby("group", dropna=False).size()
    for key, sub in data.groupby("group", dropna=False):
        labelable = sub["label"].dropna()
        positives = int((labelable == 1).sum())
        negatives = int((labelable == 0).sum())
        n_labelable = positives + negatives
        rows.append(
            {
                "group": key,
                "labelable": n_labelable,
                "positive_rate": positives / n_labelable if n_labelable else math.nan,
                "labelable_rate": n_labelable / int(total_counts.loc[key]) if int(total_counts.loc[key]) else math.nan,
            }
        )
    return pd.DataFrame(rows).set_index("group") if rows else pd.DataFrame(columns=["labelable", "positive_rate", "labelable_rate"])


def _write_provenance_column_summary(df: pd.DataFrame, path: Path) -> str:
    rows: list[dict[str, Any]] = []
    for category, columns in PROVENANCE_COLUMNS.items():
        for col in columns:
            present = col in df.columns
            row: dict[str, Any] = {"category": category, "column": col, "present": present}
            if present:
                series = df[col]
                row.update(
                    {
                        "non_null": int(series.notna().sum()),
                        "missing": int(series.isna().sum()),
                        "missing_fraction": float(series.isna().mean()),
                        "n_unique": int(series.nunique(dropna=True)),
                        "dtype": str(series.dtype),
                    }
                )
            rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _write_provenance_feature_summary(df: pd.DataFrame, path: Path) -> str:
    rows: list[dict[str, Any]] = []
    provenance_cols = sorted({col for cols in PROVENANCE_COLUMNS.values() for col in cols if col in df.columns})
    for col in provenance_cols:
        series = df[col]
        row: dict[str, Any] = {
            "column": col,
            "dtype": str(series.dtype),
            "n_rows": int(len(series)),
            "non_null": int(series.notna().sum()),
            "missing_fraction": float(series.isna().mean()),
            "n_unique": int(series.nunique(dropna=True)),
        }
        if ptypes.is_numeric_dtype(series):
            numeric = pd.to_numeric(series, errors="coerce")
            row.update(
                {
                    "min": _safe_float(numeric.min()),
                    "median": _safe_float(numeric.median()),
                    "max": _safe_float(numeric.max()),
                    "mean": _safe_float(numeric.mean()),
                }
            )
        else:
            counts = series.dropna().astype(str).value_counts().head(5)
            row["top_values_json"] = json.dumps(counts.to_dict(), sort_keys=True)
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _write_all_label_summaries(df: pd.DataFrame, labels: Sequence[str], path: Path) -> str:
    rows: list[dict[str, Any]] = []
    for label in labels:
        balance = label_balance(df, label)
        rows.append({"label": label, **balance})
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _write_group_positive_rates(
    df: pd.DataFrame,
    labels: Sequence[str],
    group_cols: Sequence[str],
    path: Path,
) -> str:
    frames: list[pd.DataFrame] = []
    for label in labels:
        for group_col in group_cols:
            if group_col not in df.columns:
                continue
            table = _positive_rate_by_group(df, label, group_col)
            if not table.empty:
                frames.append(table)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(path, index=False)
    return str(path)


def _positive_rate_by_group(df: pd.DataFrame, label: str, group_col: str) -> pd.DataFrame:
    labels = binary_label_series(df[label])
    raw_numeric = pd.to_numeric(df[label], errors="coerce")
    excluded = raw_numeric.eq(-1)
    work = pd.DataFrame(
        {
            "group_value": df[group_col].fillna("<missing>").astype(str),
            "label": labels,
            "excluded": excluded,
        }
    )
    rows: list[dict[str, Any]] = []
    for group_value, group in work.groupby("group_value", dropna=False):
        positives = int(group["label"].eq(1).sum())
        negatives = int(group["label"].eq(0).sum())
        excluded_count = int(group["excluded"].sum())
        labelable = positives + negatives
        unknown = int(len(group) - labelable - excluded_count)
        rows.append(
            {
                "label": label,
                "group_col": group_col,
                "group_value": group_value,
                "n_rows": int(len(group)),
                "n_labelable": labelable,
                "n_positive": positives,
                "n_negative": negatives,
                "n_unknown": unknown,
                "n_excluded_or_conflicting": excluded_count,
                "positive_rate": positives / labelable if labelable else math.nan,
                "labelable_rate": labelable / len(group) if len(group) else math.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["label", "group_col", "n_rows"], ascending=[True, True, False])


def _run_multi_label_mutual_information(
    df: pd.DataFrame,
    labels: Sequence[str],
    out: Path,
    *,
    max_features: int,
    max_category_levels: int,
    sample_size: int,
    seed: int,
) -> dict[str, str]:
    if not labels:
        return {}
    work = _sample(df, sample_size, seed)
    rows: list[dict[str, Any]] = []
    excluded_rows: list[dict[str, Any]] = []
    for label in labels:
        label_series = binary_label_series(work[label])
        keep = label_series.notna()
        if int(keep.sum()) < 10 or label_series[keep].nunique() < 2:
            excluded_rows.append({"label": label, "feature": "", "reason": "label_has_too_few_two-class_rows"})
            continue
        features, excluded = _select_mi_features(
            work.loc[keep],
            label,
            max_features=max_features,
            max_category_levels=max_category_levels,
        )
        excluded_rows.extend({"label": label, **row} for row in excluded)
        if not features:
            continue
        encoded, discrete_mask = _encode_features_for_mi(work.loc[keep, features])
        try:
            feature_selection = importlib.import_module("sklearn.feature_selection")
            scores = feature_selection.mutual_info_classif(
                encoded,
                label_series.loc[keep].astype(int).to_numpy(),
                discrete_features=discrete_mask,
                random_state=seed,
            )
        except Exception as exc:
            excluded_rows.append({"label": label, "feature": "", "reason": f"mutual_info_failed: {type(exc).__name__}: {exc}"})
            continue
        for feature, score, discrete in zip(features, scores, discrete_mask):
            rows.append(
                {
                    "label": label,
                    "feature": feature,
                    "mutual_information": float(score),
                    "discrete": bool(discrete),
                    "dtype": str(work[feature].dtype),
                    "n_unique": int(work[feature].nunique(dropna=True)),
                    "missing_fraction": float(work[feature].isna().mean()),
                }
            )
    ranking = pd.DataFrame(rows).sort_values(["label", "mutual_information"], ascending=[True, False]) if rows else pd.DataFrame(rows)
    excluded = pd.DataFrame(excluded_rows)
    ranking_path = out / "mutual_information_by_label.csv"
    excluded_path = out / "mutual_information_by_label_excluded_features.csv"
    ranking.to_csv(ranking_path, index=False)
    excluded.to_csv(excluded_path, index=False)
    return {
        "mutual_information_by_label": str(ranking_path),
        "mutual_information_by_label_excluded_features": str(excluded_path),
    }


def _select_mi_features(
    frame: pd.DataFrame,
    label: str,
    *,
    max_features: int,
    max_category_levels: int,
) -> tuple[list[str], list[dict[str, Any]]]:
    candidate_cols = sorted(
        {
            *[col for cols in PROVENANCE_COLUMNS.values() for col in cols],
            *DEFAULT_EXTENDED_GROUP_COLS,
            "atlas_score",
            "consensus_score",
            "consensus_z_score",
            "banana_score_normalized",
            "binding_expert_score",
            "structure_quality",
            "rdkit_mol_wt",
            "rdkit_mol_logp",
            "rdkit_tpsa",
            "rdkit_hbd",
            "rdkit_hba",
            "rdkit_rotatable_bonds",
            "rdkit_formal_charge",
            "rdkit_aromatic_rings",
            "rdkit_fraction_csp3",
            "rdkit_qed",
        }
    )
    selected: list[tuple[tuple[int, float, int, str], str]] = []
    excluded: list[dict[str, Any]] = []
    for col in candidate_cols:
        if col == label or col not in frame.columns:
            continue
        series = frame[col]
        n_unique = int(series.nunique(dropna=True))
        if n_unique <= 1:
            excluded.append({"feature": col, "reason": "constant_or_empty", "n_unique": n_unique})
            continue
        if col in FORBIDDEN_FEATURES:
            excluded.append({"feature": col, "reason": "forbidden_predictive_feature_id", "n_unique": n_unique})
            continue
        if not ptypes.is_numeric_dtype(series) and n_unique > max_category_levels:
            excluded.append({"feature": col, "reason": "high_cardinality", "n_unique": n_unique})
            continue
        selected.append(((0 if ptypes.is_numeric_dtype(series) else 1, float(series.isna().mean()), n_unique, col), col))
    ordered = [col for _, col in sorted(selected, key=lambda item: item[0])]
    for col in ordered[max_features:]:
        excluded.append({"feature": col, "reason": "over_max_features", "n_unique": int(frame[col].nunique(dropna=True))})
    return ordered[:max_features], excluded


def _encode_features_for_mi(frame: pd.DataFrame) -> tuple[np.ndarray, list[bool]]:
    arrays: list[np.ndarray] = []
    discrete: list[bool] = []
    for col in frame.columns:
        series = frame[col]
        if ptypes.is_numeric_dtype(series):
            numeric = pd.to_numeric(series, errors="coerce")
            median = numeric.median()
            arrays.append(numeric.fillna(0.0 if pd.isna(median) else float(median)).to_numpy(dtype=float))
            discrete.append(False)
        else:
            encoded, _ = pd.factorize(series.fillna("<missing>").astype(str), sort=True)
            arrays.append(encoded.astype(float))
            discrete.append(True)
    return np.column_stack(arrays), discrete


def _run_permutation_importance_suite(
    df: pd.DataFrame,
    labels: Sequence[str],
    out: Path,
    *,
    sample_size: int,
    repeats: int,
    model_type: str,
    seed: int,
    test_fraction: float,
) -> dict[str, str]:
    inspection = importlib.import_module("sklearn.inspection")
    metrics = importlib.import_module("sklearn.metrics")
    outputs: dict[str, str] = {}
    frames: list[pd.DataFrame] = []
    raw_frames: list[pd.DataFrame] = []
    skipped_rows: list[dict[str, Any]] = []
    for label in labels:
        feature_set = DEFAULT_LABEL_FEATURE_SETS.get(label)
        if not feature_set:
            skipped_rows.append({"label": label, "reason": "no_default_feature_set"})
            continue
        try:
            features = [feature for feature in get_feature_set(feature_set) if feature in df.columns and feature not in FORBIDDEN_FEATURES]
        except ValueError as exc:
            skipped_rows.append({"label": label, "reason": str(exc)})
            continue
        if not features:
            skipped_rows.append({"label": label, "reason": "no_available_features", "feature_set": feature_set})
            continue
        labels_binary = binary_label_series(df[label])
        data = df.loc[labels_binary.notna()].copy()
        data[label] = labels_binary.loc[labels_binary.notna()].astype(int)
        if len(data) > sample_size > 0:
            data = data.sample(n=sample_size, random_state=seed)
        if data[label].nunique() < 2 or len(data) < 25:
            skipped_rows.append({"label": label, "reason": "too_few_two_class_rows", "feature_set": feature_set, "n": int(len(data))})
            continue
        try:
            train_idx, test_idx = make_split(data, "random", seed=seed, test_fraction=test_fraction)
        except Exception as exc:
            skipped_rows.append({"label": label, "reason": f"split_failed: {exc}", "feature_set": feature_set})
            continue
        train = data.loc[train_idx]
        test = data.loc[test_idx]
        if train[label].nunique() < 2 or test[label].nunique() < 2:
            skipped_rows.append({"label": label, "reason": "train_or_test_one_class", "feature_set": feature_set})
            continue
        try:
            x_train, preprocessing = _fit_design_matrix(train, features)
            x_test = _transform_design_matrix(test, preprocessing, strict_raw_features=False, strict_categories=False)
            model = _model(model_type, seed, class_weight="balanced", model_params={"n_jobs": 4})
            model.fit(x_train, train[label].astype(int))
            scorer = "average_precision" if int(test[label].sum()) > 0 else "roc_auc"
            result = inspection.permutation_importance(
                model,
                x_test,
                test[label].astype(int),
                scoring=scorer,
                n_repeats=repeats,
                random_state=seed,
                n_jobs=1,
            )
            probabilities = _model_probabilities(model, x_test)
            baseline_ap = metrics.average_precision_score(test[label].astype(int), probabilities)
            baseline_auc = (
                metrics.roc_auc_score(test[label].astype(int), probabilities)
                if test[label].nunique() == 2
                else math.nan
            )
        except Exception as exc:
            skipped_rows.append({"label": label, "reason": f"permutation_failed: {type(exc).__name__}: {exc}", "feature_set": feature_set})
            continue
        design_rows = []
        for col, mean, std in zip(x_test.columns, result.importances_mean, result.importances_std):
            design_rows.append(
                {
                    "label": label,
                    "feature_set": feature_set,
                    "model_type": model_type,
                    "scoring": scorer,
                    "design_feature": str(col),
                    "importance_mean": float(mean),
                    "importance_std": float(std),
                    "baseline_average_precision": float(baseline_ap),
                    "baseline_auroc": _safe_float(baseline_auc),
                    "n_train": int(len(train)),
                    "n_test": int(len(test)),
                    "test_positive_rate": float(test[label].mean()),
                }
            )
        design = pd.DataFrame(design_rows).sort_values("importance_mean", ascending=False)
        frames.append(design)
        raw_frames.append(_aggregate_permutation_importance_to_raw_features(design, features))
    design_all = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    raw_all = pd.concat(raw_frames, ignore_index=True) if raw_frames else pd.DataFrame()
    skipped = pd.DataFrame(skipped_rows)
    design_path = out / "permutation_importance_design_features.csv"
    raw_path = out / "permutation_importance_raw_features.csv"
    skipped_path = out / "permutation_importance_skipped.csv"
    design_all.to_csv(design_path, index=False)
    raw_all.to_csv(raw_path, index=False)
    skipped.to_csv(skipped_path, index=False)
    outputs.update(
        {
            "permutation_importance_design_features": str(design_path),
            "permutation_importance_raw_features": str(raw_path),
            "permutation_importance_skipped": str(skipped_path),
        }
    )
    return outputs


def _aggregate_permutation_importance_to_raw_features(design: pd.DataFrame, raw_features: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feature in raw_features:
        mask = design["design_feature"].eq(feature)
        mask |= design["design_feature"].str.startswith(f"{feature}_", na=False)
        mask |= design["design_feature"].eq(f"{feature}__missing")
        subset = design.loc[mask]
        if subset.empty:
            continue
        first = subset.iloc[0]
        rows.append(
            {
                "label": first["label"],
                "feature_set": first["feature_set"],
                "model_type": first["model_type"],
                "raw_feature": feature,
                "importance_sum": float(subset["importance_mean"].sum()),
                "importance_max": float(subset["importance_mean"].max()),
                "n_design_features": int(len(subset)),
                "baseline_average_precision": first["baseline_average_precision"],
                "baseline_auroc": first["baseline_auroc"],
                "n_train": first["n_train"],
                "n_test": first["n_test"],
                "test_positive_rate": first["test_positive_rate"],
            }
        )
    return pd.DataFrame(rows).sort_values(["label", "importance_sum"], ascending=[True, False]) if rows else pd.DataFrame()


def _run_pandera_validation(df: pd.DataFrame, labels: Sequence[str], out: Path) -> dict[str, str]:
    try:
        import pandera.pandas as pa
        from pandera.pandas import Check, Column, DataFrameSchema
        from pandera.errors import SchemaErrors
    except Exception as exc:  # pragma: no cover - optional dependency
        raise OptionalDependencyMissing(f"install optional package pandera: {exc}") from exc

    columns: dict[str, Any] = {}
    for col in ["drug_id", "target_id", "pdb_id"]:
        if col in df.columns:
            columns[col] = Column(pa.String, nullable=False, coerce=True)
    for label in labels:
        if label in df.columns:
            columns[label] = Column(
                float,
                checks=Check(lambda series: _label_values_are_valid(series), element_wise=False),
                nullable=True,
                required=True,
                coerce=True,
            )
    for col in [
        "consensus_score",
        "consensus_z_score",
        "z_selected",
        "banana_score_normalized",
        "structure_quality",
        "rdkit_mol_wt",
        "rdkit_mol_logp",
    ]:
        if col in df.columns:
            columns[col] = Column(float, nullable=True, coerce=True)
    for col in ["target_family", "protein_class", "scaffold_key", "ligand_chemotype", "source_family", "endpoint_type"]:
        if col in df.columns:
            columns[col] = Column(pa.String, nullable=True, coerce=True)
    schema = DataFrameSchema(columns, strict=False, coerce=False)
    status: dict[str, Any] = {"status": "passed", "n_failure_cases": 0}
    failures = pd.DataFrame()
    try:
        schema.validate(df, lazy=True)
    except SchemaErrors as exc:
        status["status"] = "failed"
        failures = exc.failure_cases
        status["n_failure_cases"] = int(len(failures))
    status_path = out / "pandera_validation_summary.json"
    failure_path = out / "pandera_validation_failures.csv"
    status_path.write_text(json.dumps(status, indent=2, sort_keys=True, default=str), encoding="utf-8")
    failures.to_csv(failure_path, index=False)
    return {
        "pandera_validation_summary": str(status_path),
        "pandera_validation_failures": str(failure_path),
    }


def _label_values_are_valid(series: pd.Series) -> pd.Series:
    missing = series.isna()
    text = series.astype("object").where(~missing, "").astype(str).str.strip().str.lower()
    valid_text = text.isin({"", "0", "0.0", "1", "1.0", "-1", "-1.0", "true", "false", "yes", "no", "positive", "negative", "active", "inactive"})
    numeric = pd.to_numeric(series, errors="coerce")
    valid_numeric = numeric.isin([-1, 0, 1])
    return missing | valid_text | valid_numeric


def _run_evidently_drift_suite(
    df: pd.DataFrame,
    out: Path,
    *,
    labels: Sequence[str],
    splits: Sequence[str],
    sample_size: int,
    seed: int,
    test_fraction: float,
) -> dict[str, str]:
    try:
        from evidently import Report
        from evidently.presets import DataDriftPreset
    except Exception as exc:  # pragma: no cover - optional dependency
        raise OptionalDependencyMissing(f"install optional package evidently: {exc}") from exc

    numeric_cols = _evidently_numeric_columns(df, labels)
    if not numeric_cols:
        raise ValueError("no numeric columns available for Evidently drift")
    drift_dir = out / "evidently_drift"
    drift_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    outputs: dict[str, str] = {}
    for split in splits:
        try:
            data = df[numeric_cols].copy()
            train_idx, test_idx = make_split(df, split, seed=seed, test_fraction=test_fraction)
            reference = data.loc[train_idx]
            current = data.loc[test_idx]
            reference = _sample(reference, sample_size, seed)
            current = _sample(current, sample_size, seed + 1)
            if len(reference) < 10 or len(current) < 10:
                rows.append({"split": split, "status": "skipped", "reason": "too_few_rows"})
                continue
            report = Report([DataDriftPreset()])
            snapshot = report.run(current, reference)
            html_path = drift_dir / f"{split}_data_drift.html"
            json_path = drift_dir / f"{split}_data_drift.json"
            snapshot.save_html(str(html_path))
            snapshot.save_json(str(json_path))
            rows.append(
                {
                    "split": split,
                    "status": "ok",
                    "n_reference": int(len(reference)),
                    "n_current": int(len(current)),
                    "n_columns": int(len(numeric_cols)),
                    "html": str(html_path),
                    "json": str(json_path),
                }
            )
            outputs[f"evidently_drift_{split}_html"] = str(html_path)
            outputs[f"evidently_drift_{split}_json"] = str(json_path)
        except Exception as exc:
            rows.append({"split": split, "status": "failed", "reason": f"{type(exc).__name__}: {exc}"})
    summary_path = drift_dir / "evidently_drift_summary.csv"
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    outputs["evidently_drift_summary"] = str(summary_path)
    return outputs


def _evidently_numeric_columns(df: pd.DataFrame, labels: Sequence[str]) -> list[str]:
    preferred = [
        "atlas_score",
        "consensus_score",
        "consensus_z_score",
        "banana_score_normalized",
        "binding_expert_score",
        "structure_quality",
        "rdkit_mol_wt",
        "rdkit_mol_logp",
        "rdkit_tpsa",
        "rdkit_hbd",
        "rdkit_hba",
        "rdkit_rotatable_bonds",
        "rdkit_formal_charge",
        "rdkit_aromatic_rings",
        "rdkit_fraction_csp3",
        "rdkit_qed",
        "site_relevance_score",
        "expression_presence_score",
        "site_specificity_score",
        "expression_concordance_score",
        "hpa_consensus_ntpm",
        "gtex_median_tpm",
        "bgee_expression_score",
        "ot_expression_value",
        "ot_expression_zscore",
        *labels,
    ]
    cols = []
    for col in preferred:
        if col in df.columns and ptypes.is_numeric_dtype(df[col]) and df[col].nunique(dropna=True) > 1:
            cols.append(col)
    return cols[:80]


def _write_audit_column_dictionary(augmented: pd.DataFrame, original_columns: Iterable[str], path: Path) -> str:
    original = set(original_columns)
    rows: list[dict[str, Any]] = []
    for col in augmented.columns:
        if col in original:
            continue
        rows.append(
            {
                "column": col,
                "kind": "audit_provenance_label_derived",
                "allowed_as_clean_model_feature": False,
                "reason": "group positive rates and assay densities are diagnostics for shortcut/source bias, not nonleaky predictive features",
            }
        )
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _optional_dependency_status(modules: Sequence[str]) -> pd.DataFrame:
    rows = []
    for module in modules:
        try:
            imported = importlib.import_module(module)
            rows.append(
                {
                    "module": module,
                    "available": True,
                    "version": str(getattr(imported, "__version__", "")),
                    "error": "",
                }
            )
        except Exception as exc:
            rows.append({"module": module, "available": False, "version": "", "error": f"{type(exc).__name__}: {exc}"})
    return pd.DataFrame(rows)


def _run_step(
    name: str,
    outputs: dict[str, str],
    skipped: dict[str, str],
    errors: dict[str, str],
    fail_on_missing: bool,
    callback: Any,
) -> None:
    try:
        produced = callback()
        if produced:
            outputs.update(produced)
    except OptionalDependencyMissing as exc:
        skipped[name] = str(exc)
        if fail_on_missing:
            raise
    except Exception as exc:
        errors[name] = f"{type(exc).__name__}: {exc}"
        if fail_on_missing:
            raise


def _available(values: Sequence[str], df: pd.DataFrame) -> list[str]:
    return [value for value in values if value in df.columns]


def _safe_column_token(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in str(value).strip()).strip("_").lower()


def _sample(frame: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    if n <= 0 or len(frame) <= n:
        return frame.copy()
    return frame.sample(n=n, random_state=seed).copy()


def _safe_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def _write_markdown_summary(path: Path, manifest: Mapping[str, Any]) -> None:
    lines = [
        "# ML Provenance Profile Audit",
        "",
        f"- Dataset: `{manifest['dataset']}`",
        f"- Rows: `{manifest['n_rows']}`",
        f"- Columns: `{manifest['n_columns']}`",
        f"- Labels: `{', '.join(manifest['labels'])}`",
        "",
        "## Policy",
        "- Positive-rate and assay-density columns are label-derived audit diagnostics.",
        "- Do not add `_audit_*_positive_rate`, `_audit_*_assay_density`, or `_audit_*_labelable_rate` columns to clean predictive feature sets.",
        "- Use these outputs to decide where grouped holdouts, source balancing, or additional data acquisition are required.",
        "",
        "## Outputs",
    ]
    for name, output in sorted(dict(manifest.get("outputs", {})).items()):
        lines.append(f"- `{name}`: `{output}`")
    if manifest.get("skipped"):
        lines.extend(["", "## Skipped"])
        for name, reason in sorted(dict(manifest["skipped"]).items()):
            lines.append(f"- `{name}`: {reason}")
    if manifest.get("errors"):
        lines.extend(["", "## Errors"])
        for name, reason in sorted(dict(manifest["errors"]).items()):
            lines.append(f"- `{name}`: {reason}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
