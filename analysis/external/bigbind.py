from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import apply_pair_mapping, normalize_columns, read_optional_mapping
from analysis.external.source_tables import read_source_table


BIGBIND_ALIASES = {
    "drug_id": ["drug_id", "compound_id", "ligand_id", "molecule_id", "chembl_id", "molecule_chembl_id"],
    "target_id": ["target_id", "uniprot", "protein_id", "target_chembl_id", "pocket_id"],
    "activity_nM": ["activity_nM", "activity_nm", "standard_value", "value", "affinity_nM"],
    "activity_relation": ["activity_relation", "standard_relation", "relation"],
    "activity_type": ["activity_type", "standard_type", "type"],
    "bigbind_split": ["bigbind_split", "split", "pocket_split"],
    "bigbind_pocket_id": ["bigbind_pocket_id", "pocket_id", "crossdocked_pocket_id"],
    "bigbind_is_putative_inactive": ["bigbind_is_putative_inactive", "putative_inactive", "is_putative_inactive"],
    "source": ["source", "dataset"],
}


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower().isin({"1", "1.0", "true", "yes", "y"})


def load_bigbind(
    path: str | Path,
    mapping_path: str | Path | None = None,
    *,
    active_threshold_nM: float = 10000.0,
) -> pd.DataFrame:
    """Load a local BigBind-style activity table.

    BigBind is treated as a docking/bioactivity benchmark source. It is not an
    exposure or ADR mechanism label source.
    """

    df = normalize_columns(read_source_table(path), BIGBIND_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    df["activity_nM"] = pd.to_numeric(df["activity_nM"], errors="coerce")
    relation = df.get("activity_relation", pd.Series("", index=df.index)).fillna("").astype(str).str.strip()
    exact = relation.eq("") | relation.eq("=")
    upper_bound = relation.isin({"<", "<="})
    lower_bound = relation.str.startswith(">")
    putative_inactive = _truthy(df.get("bigbind_is_putative_inactive", pd.Series(False, index=df.index)))
    active = (exact | upper_bound) & df["activity_nM"].notna() & df["activity_nM"].le(active_threshold_nM)
    inactive = (
        ((exact | lower_bound) & df["activity_nM"].notna() & df["activity_nM"].gt(active_threshold_nM))
        | putative_inactive
    )
    label = pd.Series(pd.NA, index=df.index, dtype="Int64")
    label = label.mask(inactive, 0)
    label = label.mask(active, 1)
    df["bigbind_active"] = label
    df["bigbind_label_status"] = "labeled"
    df.loc[label.isna(), "bigbind_label_status"] = "unknown_unlabelable_relation_or_activity"
    df.loc[putative_inactive, "bigbind_label_status"] = "labeled_putative_inactive"
    df["source"] = df["source"].fillna("BigBind")
    return df


def aggregate_bigbind(df: pd.DataFrame, active_threshold_nM: float = 10000.0) -> pd.DataFrame:
    work = df.copy()
    work["_active_sort"] = work["bigbind_active"].eq(1)
    work["_known_sort"] = work["bigbind_active"].notna()
    work = work.sort_values(
        ["drug_id", "target_id", "_known_sort", "_active_sort", "activity_nM"],
        ascending=[True, True, False, False, True],
    )
    grouped = (
        work.groupby(["drug_id", "target_id"], dropna=False)
        .agg(
            bigbind_activity_nM=("activity_nM", "min"),
            bigbind_active=("bigbind_active", "first"),
            bigbind_activity_type=("activity_type", "first"),
            bigbind_relation=("activity_relation", "first"),
            bigbind_pocket_ids=("bigbind_pocket_id", lambda s: ";".join(sorted({str(v) for v in s.dropna()}))),
            bigbind_splits=("bigbind_split", lambda s: ";".join(sorted({str(v) for v in s.dropna()}))),
            bigbind_assay_count=("activity_nM", "size"),
            bigbind_label_status=("bigbind_label_status", "first"),
            bigbind_source=("source", lambda s: ";".join(sorted({str(v) for v in s.dropna()}))),
        )
        .reset_index()
    )
    grouped["bigbind_threshold_nM"] = active_threshold_nM
    return grouped


def build_bigbind_benchmark(
    pair_table_path: str | Path,
    bigbind_path: str | Path,
    mapping_path: str | Path | None,
    out_path: str | Path,
    *,
    active_threshold_nM: float = 10000.0,
) -> pd.DataFrame:
    pair = read_source_table(pair_table_path)
    for key in ("drug_id", "target_id"):
        if key in pair.columns:
            pair[key] = pair[key].fillna("").astype(str)
    bigbind = aggregate_bigbind(
        load_bigbind(bigbind_path, mapping_path, active_threshold_nM=active_threshold_nM),
        active_threshold_nM=active_threshold_nM,
    )
    for key in ("drug_id", "target_id"):
        if key in bigbind.columns:
            bigbind[key] = bigbind[key].fillna("").astype(str)
    merged = pair.merge(bigbind, on=["drug_id", "target_id"], how="left")
    no_match = merged["bigbind_label_status"].isna()
    merged.loc[no_match, "bigbind_label_status"] = "unknown_no_bigbind_pair_match"
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    summary = pd.DataFrame(
        [
            {
                "n_pair_rows": len(pair),
                "n_bigbind_pairs": len(bigbind),
                "n_joined_rows": int((~no_match).sum()),
                "n_labelable_rows": int(merged["bigbind_active"].notna().sum()),
                "n_active": int(merged["bigbind_active"].eq(1).sum()),
                "n_inactive": int(merged["bigbind_active"].eq(0).sum()),
            }
        ]
    )
    summary.to_csv(out.with_name("bigbind_mapping_summary.csv"), index=False)
    return merged
