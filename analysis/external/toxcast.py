from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import apply_pair_mapping, normalize_columns, read_optional_mapping
from analysis.external.source_tables import read_source_table


TOXCAST_ALIASES = {
    "drug_id": ["drug_id", "chemical_id", "chid", "dsstox_substance_id", "compound"],
    "target_id": ["target_id", "gene", "gene_symbol", "uniprot", "assay_target"],
    "assay_id": ["assay_id", "aeid", "assay"],
    "assay_component": ["assay_component", "component", "assay_component_name"],
    "ac50_nM": ["ac50_nM", "ac50_nm", "ac50", "modl_ac50"],
    "hit_call": ["hit_call", "hitc", "active", "hit"],
    "modl_ga": ["modl_ga", "ga"],
    "source": ["source", "dataset"],
}


def _read_table(path: str | Path) -> pd.DataFrame:
    return read_source_table(path)


def normalize_toxcast_table(toxcast_path: str | Path, mapping_path: str | Path | None = None) -> pd.DataFrame:
    df = normalize_columns(_read_table(toxcast_path), TOXCAST_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    df["ac50_nM"] = pd.to_numeric(df["ac50_nM"], errors="coerce")
    return df


def build_toxcast_benchmark(
    pair_table_path: str | Path,
    toxcast_path: str | Path,
    mapping_path: str | Path | None,
    out_path: str | Path,
    potent_threshold_nM: float = 10000.0,
) -> pd.DataFrame:
    pair = _read_table(pair_table_path)
    tox = normalize_toxcast_table(toxcast_path, mapping_path)
    active_text = tox["hit_call"].astype(str).str.lower()
    tox["toxcast_active"] = active_text.isin({"1", "true", "active", "hit", "positive"})
    tox["toxcast_potent"] = pd.to_numeric(tox["ac50_nM"], errors="coerce") <= potent_threshold_nM
    tox["toxcast_label_status"] = "labeled"
    tox.loc[tox["hit_call"].isna(), "toxcast_label_status"] = "unknown_missing_hit_call"
    grouped = (
        tox.sort_values(["drug_id", "target_id", "toxcast_active", "ac50_nM"], ascending=[True, True, False, True])
        .groupby(["drug_id", "target_id"], dropna=False)
        .agg(
            assay_id=("assay_id", lambda s: ";".join(sorted({str(v) for v in s.dropna()}))),
            assay_component=("assay_component", "first"),
            ac50_nM=("ac50_nM", "min"),
            hit_call=("hit_call", "first"),
            modl_ga=("modl_ga", "first"),
            source=("source", "first"),
            toxcast_active=("toxcast_active", "max"),
            toxcast_potent=("toxcast_potent", "max"),
            toxcast_assayed=("assay_id", "count"),
            toxcast_label_status=("toxcast_label_status", "first"),
        )
        .reset_index()
    )
    merged = pair.merge(grouped, on=["drug_id", "target_id"], how="left", suffixes=("", "_toxcast"))
    no_match = merged["toxcast_assayed"].isna() if "toxcast_assayed" in merged else pd.Series(False, index=merged.index)
    merged.loc[no_match, "toxcast_label_status"] = "unknown_no_toxcast_pair_match"
    merged.loc[no_match, "toxcast_missing_reason"] = "unknown_no_toxcast_pair_match"
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    return merged
