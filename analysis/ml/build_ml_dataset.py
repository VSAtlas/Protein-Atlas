from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.ml.feature_sets import (
    RETAINED_CONTEXT_AUDIT_COLUMNS,
    effective_exclude_features,
    get_feature_set,
    label_definition_exclude_features,
)
from analysis.ml.labels import binary_label_series, truthy_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.materialize_final_score import (
    canonical_atlas_score_source,
    canonical_final_score_source,
)


def _score_sort_col(df: pd.DataFrame) -> str | None:
    for col in (
        "consensus_z_score",
        "atlas_score",
        "z_selected",
        "priority_score",
        "final_score",
        "consensus_score",
    ):
        if col in df.columns:
            return col
    return None


def _score_source(
    frame: pd.DataFrame, column: str, *, default: str = ""
) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype="object")
    return frame[column].fillna("").astype(str).str.strip()


def _derive_atlas_score(
    frame: pd.DataFrame,
    *,
    value_column: str,
    source_column: str,
    default_source: str = "",
    final_score: bool = False,
) -> None:
    values = pd.to_numeric(frame[value_column], errors="coerce")
    sources = _score_source(frame, source_column, default=default_source)
    canonicalizer = (
        canonical_final_score_source
        if final_score
        else canonical_atlas_score_source
    )
    accepted = sources.map(canonicalizer).ne("")
    frame["atlas_score"] = values.where(accepted)
    frame["atlas_score_source_for_ml"] = sources


def add_standard_ml_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "drug_id" not in out.columns:
        for source in ("dedup_drug_key", "drugcentral_id", "mapped_drug_name", "ligand_base"):
            if source in out.columns:
                out["drug_id"] = out[source]
                break
    if "target_id" not in out.columns:
        for source in ("dedup_target_key", "target_uniprot", "target_gene", "pdb_id"):
            if source in out.columns:
                out["target_id"] = out[source]
                break
    if "atlas_score" not in out.columns:
        if "consensus_z_score" in out.columns:
            _derive_atlas_score(
                out,
                value_column="consensus_z_score",
                source_column="consensus_z_score_source",
                default_source="consensus_z_score",
            )
        elif "z_selected" in out.columns:
            _derive_atlas_score(
                out,
                value_column="z_selected",
                source_column="z_selected_source",
            )
        elif "final_score" in out.columns:
            _derive_atlas_score(
                out,
                value_column="final_score",
                source_column="final_score_source",
                final_score=True,
            )
    if "free_cmax_um" not in out.columns and "free_cmax_uM" in out.columns:
        out["free_cmax_um"] = pd.to_numeric(out["free_cmax_uM"], errors="coerce")
    if "free_cmax" not in out.columns and "free_cmax_um" in out.columns:
        out["free_cmax"] = pd.to_numeric(out["free_cmax_um"], errors="coerce")
    if "chemical_cluster" not in out.columns:
        for source in ("chemical_cluster_id", "ecfp_cluster", "butina_cluster", "umap_cluster", "ligand_cluster"):
            if source in out.columns:
                out["chemical_cluster"] = out[source]
                break
    if "target_family" not in out.columns:
        for source in ("protein_family", "target_class", "protein_class"):
            if source in out.columns:
                out["target_family"] = out[source]
                break
    return out


def deduplicate_ml_rows(
    df: pd.DataFrame,
    label_col: str,
    group_cols: list[str] | None = None,
) -> pd.DataFrame:
    groups = group_cols or [col for col in ("drug_id", "target_id") if col in df.columns]
    if len(groups) < 2:
        return df.copy()
    work = df.copy()
    labels = binary_label_series(work[label_col])
    work["_atlas_ml_label_sort"] = labels.fillna(0).astype(int)
    score_col = _score_sort_col(work)
    if score_col is not None:
        work["_atlas_ml_score_sort"] = pd.to_numeric(work[score_col], errors="coerce")
    else:
        work["_atlas_ml_score_sort"] = 0.0
    work = work.sort_values(
        [*groups, "_atlas_ml_label_sort", "_atlas_ml_score_sort"],
        ascending=[*(True for _ in groups), False, False],
    )
    return work.drop_duplicates(groups, keep="first").drop(
        columns=["_atlas_ml_label_sort", "_atlas_ml_score_sort"]
    )


def _unique_existing_columns(columns: list[str], df: pd.DataFrame) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for col in columns:
        if col in df.columns and col not in seen:
            out.append(col)
            seen.add(col)
    return out


def build_ml_dataset(
    pair_table_path: str | Path,
    label_col: str,
    feature_set: str,
    out_path: str | Path,
    *,
    deduplicate_drug_target: bool = False,
    group_cols: list[str] | None = None,
    eligible_col: str | None = None,
    eligible_value: str = "1",
    include_unlabeled_as_background: bool = False,
    unlabeled_weight: float = 0.25,
    label_source_include: list[str] | None = None,
    require_year_col: str | None = None,
    require_nonmissing_cols: list[str] | None = None,
    exclude_flag_cols: list[str] | None = None,
    exclude_features: list[str] | None = None,
    allow_label_definition_features: bool = False,
    keep_all_columns: bool = False,
) -> pd.DataFrame:
    df = add_standard_ml_columns(pd.read_csv(pair_table_path, low_memory=False))
    if eligible_col and eligible_col in df.columns:
        df = df[df[eligible_col].astype(str).str.strip().eq(str(eligible_value))].copy()
    if label_source_include and "label_source" in df.columns:
        allowed = {str(value).strip() for value in label_source_include if str(value).strip()}
        source_text = df["label_source"].fillna("").astype(str)
        mask = pd.Series(False, index=df.index)
        for source in allowed:
            mask |= source_text.eq(source) | source_text.str.split(";").map(lambda parts: source in parts)
        df = df[mask].copy()
    if require_year_col:
        if require_year_col not in df.columns:
            raise ValueError(f"required year column {require_year_col!r} not found")
        df = df[pd.to_numeric(df[require_year_col], errors="coerce").notna()].copy()
    if require_nonmissing_cols:
        for col in require_nonmissing_cols:
            if col not in df.columns:
                raise ValueError(f"required nonmissing column {col!r} not found")
            df = df[pd.to_numeric(df[col], errors="coerce").notna()].copy()
    if exclude_flag_cols:
        for col in exclude_flag_cols:
            if col not in df.columns:
                continue
            df = df[~truthy_series(df[col])].copy()
    parsed_label = binary_label_series(df[label_col])
    was_unlabeled = parsed_label.isna()
    df[label_col] = parsed_label
    if include_unlabeled_as_background:
        df.loc[was_unlabeled, label_col] = 0
        df["_sample_weight"] = 1.0
        df.loc[was_unlabeled, "_sample_weight"] = float(unlabeled_weight)
        df["_unlabeled_background"] = was_unlabeled
        df["negative_evidence_type"] = df.get("negative_evidence_type", pd.Series("", index=df.index))
        df.loc[was_unlabeled, "negative_evidence_type"] = "pu_background_unknown_not_confirmed_negative"
        df["negative_confidence"] = df.get("negative_confidence", pd.Series(pd.NA, index=df.index))
        df.loc[was_unlabeled, "negative_confidence"] = float(unlabeled_weight)
    label_definition_fields = label_definition_exclude_features(label_col)
    auto_label_definition_excluded = (
        set()
        if allow_label_definition_features
        else label_definition_fields
    )
    excluded = effective_exclude_features(
        label_col,
        exclude_features,
        allow_label_definition_features=allow_label_definition_features,
    )
    retained_excluded_features = sorted(feature for feature in excluded if feature in df.columns)
    features = [feature for feature in get_feature_set(feature_set) if feature in df.columns and feature not in excluded]
    retained_context_not_trained = sorted(
        col for col in RETAINED_CONTEXT_AUDIT_COLUMNS if col in df.columns and col not in features
    )
    assert_no_leakage(
        features,
        allow_label_definition_features=allow_label_definition_features,
    )
    if not features:
        raise ValueError(f"feature set {feature_set!r} has no columns in {pair_table_path}")
    metadata_cols = [
        "drug_id",
        "target_id",
        "pdb_id",
        "ligand_chemotype",
        "scaffold_key",
        "chemical_cluster",
        "ligand_cluster",
        "protein_class",
        "target_family",
        "protein_family",
        "label_source",
        "source_family",
        "upstream_source",
        "source_release",
        "source_version",
        "assay_type",
        "assay_mode",
        "endpoint_type",
        "activity_type",
        "standard_type",
        "label_source_count",
        "database_release_year",
        "activity_publication_year",
        "activity_publication_year_source",
        "free_cmax_missing_policy",
        "chembl_label_available",
        "papyrus_label_available",
        "toxcast_label_available",
        "spd_label_available",
        "mechanism_label_status",
        "negative_evidence_type",
        "negative_source",
        "negative_confidence",
        "negative_label_status",
        "negative_selection_fold",
        "negative_selection_features",
        "_unlabeled_background",
        "excluded_reason",
        "_sample_weight",
        "label_definition_leakage_policy",
        "retained_excluded_features",
        "retained_excluded_policy",
        "retained_context_not_trained",
        "retained_context_not_trained_policy",
    ]
    source_metadata_cols = [
        col
        for col in df.columns
        if col.endswith("_activity_publication_year_min")
        or col.endswith("_activity_publication_year_max")
        or col.endswith("_database_release_year")
        or col.endswith("_activity_document_ids")
    ]
    df["label_definition_leakage_policy"] = (
        "allowed_for_explicit_sensitivity_model"
        if allow_label_definition_features and label_definition_fields
        else (
            "auto_excluded:" + ",".join(sorted(auto_label_definition_excluded))
            if auto_label_definition_excluded
            else "no_label_definition_features_configured"
        )
    )
    df["retained_excluded_features"] = ";".join(retained_excluded_features)
    df["retained_excluded_policy"] = (
        "present_in_dataset_explicitly_excluded_from_training"
        if retained_excluded_features
        else "no_retained_explicit_exclusions"
    )
    df["retained_context_not_trained"] = ";".join(retained_context_not_trained)
    df["retained_context_not_trained_policy"] = (
        "present_in_dataset_not_selected_by_feature_set"
        if retained_context_not_trained
        else "no_retained_context_columns_outside_feature_set"
    )
    if keep_all_columns:
        out = df.dropna(subset=[label_col]).copy()
    else:
        keep = _unique_existing_columns([*metadata_cols, *source_metadata_cols, label_col, *features], df)
        out = df[keep].dropna(subset=[label_col]).copy()
    if deduplicate_drug_target:
        out = deduplicate_ml_rows(out, label_col, group_cols)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out
