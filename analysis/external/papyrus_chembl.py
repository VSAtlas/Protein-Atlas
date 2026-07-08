from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import apply_pair_mapping, normalize_columns, read_optional_mapping
from analysis.external.source_tables import read_source_table


BIOACTIVITY_ALIASES = {
    "drug_id": ["drug_id", "compound_id", "chembl_id", "molecule_chembl_id", "connectivity"],
    "target_id": ["target_id", "uniprot", "uniprotid", "target_chembl_id", "accession", "protein_id"],
    "activity_nM": ["activity_nM", "activity_nm", "standard_value", "value"],
    "pchembl_value": ["pchembl_value", "pchembl", "pchembl_value_mean", "pchembl_value_median"],
    "activity_relation": ["activity_relation", "standard_relation", "relation", "standard_relation"],
    "activity_type": ["activity_type", "standard_type", "type", "standard_type"],
    "activity_units": ["activity_units", "standard_units", "units"],
    "assay_id": ["assay_id", "assay_chembl_id", "aid"],
    "assay_type": ["assay_type", "assay_category"],
    "source": ["source", "dataset"],
    "activity_publication_year": ["activity_publication_year", "document_year", "year", "Year"],
    "activity_document_ids": ["activity_document_ids", "document_chembl_id", "doc_id", "all_doc_ids"],
    "database_release_year": ["database_release_year", "release_year"],
}


def _source_name(path: str | Path) -> str:
    text = str(path).lower()
    if "papyrus" in text:
        return "Papyrus"
    if "chembl" in text:
        return "ChEMBL"
    return "bioactivity"


def _activity_to_nm(df: pd.DataFrame) -> pd.Series:
    values = pd.to_numeric(df["activity_nM"], errors="coerce")
    units = df.get("activity_units", pd.Series("", index=df.index)).astype(str).str.lower()
    values = values.where(~units.isin({"um", "µm", "micromolar"}), values * 1000.0)
    values = values.where(~units.isin({"mm", "millimolar"}), values * 1_000_000.0)
    values = values.where(~units.isin({"pm", "picomolar"}), values / 1000.0)
    pchembl = pd.to_numeric(df.get("pchembl_value", pd.Series(pd.NA, index=df.index)), errors="coerce")
    return values.fillna((10.0 ** (-pchembl)) * 1_000_000_000.0)


def _nullable_any_true(series: pd.Series) -> object:
    nonmissing = series.dropna()
    if nonmissing.empty:
        return pd.NA
    return bool(nonmissing.eq(True).any())


def load_bioactivity(
    path: str | Path,
    mapping_path: str | Path | None = None,
    *,
    active_threshold_nM: float = 10000.0,
) -> pd.DataFrame:
    raw = read_source_table(path)
    df = normalize_columns(raw, BIOACTIVITY_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    df["activity_nM"] = _activity_to_nm(df)
    df["source"] = df["source"].fillna(_source_name(path))
    relation = df.get("activity_relation", pd.Series("", index=df.index)).astype(str).str.strip()
    df["bioactivity_relation"] = relation
    exact = relation.eq("") | relation.eq("=")
    upper_bound = relation.isin({"<", "<="})
    lower_bound = relation.str.startswith(">")
    measured = df["activity_nM"].notna()
    active_mask = (exact | upper_bound) & measured & df["activity_nM"].le(active_threshold_nM)
    inactive_mask = (exact | lower_bound) & measured & df["activity_nM"].gt(active_threshold_nM)
    labelable_mask = active_mask | inactive_mask
    df["bioactivity_active"] = pd.Series(pd.NA, index=df.index, dtype="boolean")
    df.loc[labelable_mask, "bioactivity_active"] = active_mask.loc[labelable_mask]
    df["bioactivity_inactive"] = pd.Series(pd.NA, index=df.index, dtype="boolean")
    df.loc[labelable_mask, "bioactivity_inactive"] = inactive_mask.loc[labelable_mask]
    df["bioactivity_label_status"] = "labeled"
    df.loc[df["activity_nM"].isna(), "bioactivity_label_status"] = "unknown_missing_activity"
    df.loc[upper_bound & measured & df["activity_nM"].gt(active_threshold_nM), "bioactivity_label_status"] = (
        "unknown_censored_activity_lt"
    )
    df.loc[lower_bound & measured & df["activity_nM"].le(active_threshold_nM), "bioactivity_label_status"] = (
        "unknown_censored_activity_gt"
    )
    if "activity_publication_year" in df.columns:
        df["activity_publication_year"] = pd.to_numeric(df["activity_publication_year"], errors="coerce")
    if "database_release_year" in df.columns:
        df["database_release_year"] = pd.to_numeric(df["database_release_year"], errors="coerce")
    return df


def aggregate_bioactivity_assays(df: pd.DataFrame, active_threshold_nM: float = 10000.0) -> pd.DataFrame:
    sortable = df.copy()
    sortable["_active_sort"] = sortable["bioactivity_active"].eq(True)
    sortable = sortable.sort_values(["drug_id", "target_id", "_active_sort", "activity_nM"], ascending=[True, True, False, True])
    grouped = (
        sortable.groupby(["drug_id", "target_id"], dropna=False)
        .agg(
            activity_nM=("activity_nM", "min"),
            activity_type=("activity_type", "first"),
            activity_relation=("bioactivity_relation", "first"),
            assay_id=("assay_id", lambda s: ";".join(sorted({str(v) for v in s.dropna()}))),
            assay_count=("assay_id", "count"),
            source=("source", lambda s: ";".join(sorted({str(v) for v in s.dropna()}))),
            activity_document_ids=(
                "activity_document_ids",
                lambda s: ";".join(sorted({str(v) for v in s.dropna() if str(v).strip()})),
            ),
            activity_publication_year_min=("activity_publication_year", "min"),
            activity_publication_year_max=("activity_publication_year", "max"),
            database_release_year=("database_release_year", "max"),
            bioactivity_active=("bioactivity_active", _nullable_any_true),
            bioactivity_inactive=("bioactivity_inactive", _nullable_any_true),
            bioactivity_label_status=("bioactivity_label_status", "first"),
        )
        .reset_index()
    )
    grouped["bioactivity_conflict"] = (
        grouped["bioactivity_active"].eq(True).fillna(False) & grouped["bioactivity_inactive"].eq(True).fillna(False)
    )
    grouped.loc[grouped["bioactivity_conflict"], "bioactivity_label_status"] = "conflicting_active_inactive"
    grouped["bioactivity_threshold_nM"] = active_threshold_nM
    return grouped


def build_bioactivity_benchmark(
    pair_table_path: str | Path,
    bioactivity_path: str | Path,
    mapping_path: str | Path | None,
    out_path: str | Path,
    *,
    active_threshold_nM: float = 10000.0,
    label_prefix: str | None = None,
) -> pd.DataFrame:
    pair = read_source_table(pair_table_path)
    for key in ("drug_id", "target_id"):
        if key in pair.columns:
            pair[key] = pair[key].fillna("").astype(str)
    label = label_prefix or _source_name(bioactivity_path).lower()
    bioactivity = aggregate_bioactivity_assays(
        load_bioactivity(bioactivity_path, mapping_path, active_threshold_nM=active_threshold_nM),
        active_threshold_nM=active_threshold_nM,
    )
    for key in ("drug_id", "target_id"):
        if key in bioactivity.columns:
            bioactivity[key] = bioactivity[key].fillna("").astype(str)
    rename = {
        "activity_nM": f"{label}_activity_nM",
        "activity_type": f"{label}_activity_type",
        "activity_relation": f"{label}_activity_relation",
        "assay_id": f"{label}_assay_ids",
        "assay_count": f"{label}_assay_count",
        "source": f"{label}_source",
        "activity_document_ids": f"{label}_activity_document_ids",
        "activity_publication_year_min": f"{label}_activity_publication_year_min",
        "activity_publication_year_max": f"{label}_activity_publication_year_max",
        "database_release_year": f"{label}_database_release_year",
        "bioactivity_active": f"{label}_active",
        "bioactivity_inactive": f"{label}_inactive",
        "bioactivity_conflict": f"{label}_conflict",
        "bioactivity_label_status": f"{label}_label_status",
        "bioactivity_threshold_nM": f"{label}_threshold_nM",
    }
    merged = pair.merge(bioactivity.rename(columns=rename), on=["drug_id", "target_id"], how="left")
    status_col = f"{label}_label_status"
    missing_col = f"{label}_missing_reason"
    no_match = merged[status_col].isna()
    merged.loc[no_match, status_col] = f"unknown_no_{label}_pair_match"
    merged.loc[no_match, missing_col] = f"unknown_no_{label}_pair_match"
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    summary = pd.DataFrame(
        [
            {
                "source": label,
                "n_pair_rows": len(pair),
                "n_source_pairs": len(bioactivity),
                "n_joined_pairs": int((~no_match).sum()),
                "n_labelable_pairs": int(merged[f"{label}_active"].notna().sum()) if f"{label}_active" in merged else 0,
                "n_active": int(merged[f"{label}_active"].eq(True).sum()) if f"{label}_active" in merged else 0,
            }
        ]
    )
    summary.to_csv(out.with_name(f"{label}_mapping_summary.csv"), index=False)
    return merged
