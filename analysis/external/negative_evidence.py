from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.mapping import apply_pair_mapping, normalize_columns, read_optional_mapping
from analysis.external.source_tables import read_source_table


BIOACTIVITY_NEGATIVE_ALIASES = {
    "drug_id": [
        "drug_id",
        "compound_id",
        "chembl_id",
        "molecule_chembl_id",
        "cid",
        "pubchem_cid",
        "substance_id",
        "sid",
        "connectivity",
    ],
    "target_id": [
        "target_id",
        "uniprot",
        "uniprotid",
        "target_chembl_id",
        "accession",
        "protein_id",
        "gene",
        "gene_symbol",
        "target",
    ],
    "assay_id": ["assay_id", "aid", "pubchem_aid", "assay_chembl_id", "assay"],
    "activity_outcome": [
        "activity_outcome",
        "PUBCHEM_ACTIVITY_OUTCOME",
        "outcome",
        "activity",
        "label",
        "curated_bioactivity_label",
        "toxcast_hts_label",
        "source_specific_activity_label",
        "bioactivity_ml_label",
        "ml_binary_label",
        "class",
        "active",
        "inactive",
    ],
    "activity_value_nM": [
        "activity_value_nM",
        "activity_nM",
        "activity_nm",
        "standard_value",
        "value",
        "affinity",
        "kd",
        "ki",
        "ic50",
        "ec50",
    ],
    "activity_relation": ["activity_relation", "standard_relation", "relation", "operator"],
    "activity_type": ["activity_type", "standard_type", "type", "endpoint"],
    "activity_units": ["activity_units", "standard_units", "units"],
    "source": ["source", "dataset", "label_source", "source_objective"],
    "publication_year": ["publication_year", "document_year", "year", "Year"],
    "document_ids": ["document_ids", "document_chembl_id", "doc_id", "pmid", "pubmed_id"],
}

DRUG_ADR_NEGATIVE_ALIASES = {
    "drug_id": [
        "drug_id",
        "exposureName",
        "exposure_name",
        "targetName",
        "target_name",
        "drug",
        "drug_name",
    ],
    "adr_id": ["adr_id", "outcomeName", "outcome_name", "event", "adr", "meddra_term", "pt"],
    "ground_truth": ["groundTruth", "ground_truth", "label", "negative_control", "is_negative"],
    "source": ["source", "dataset"],
}

FAERS_NONSIGNAL_ALIASES = {
    "drug_id": ["drug_id", "drug", "drug_name", "exposureName", "exposure_name"],
    "adr_id": ["adr_id", "event", "adr", "meddra_term", "pt", "outcomeName", "outcome_name"],
    "pair_reports": ["pair_reports", "n_pair", "count", "a"],
    "drug_reports": ["drug_reports", "n_drug", "drug_count", "a_b"],
    "event_reports": ["event_reports", "n_event", "event_count", "a_c"],
    "ror": ["ror", "ROR"],
    "ror_lower": ["ror_lower", "ror025", "ror_lower_ci", "ROR025"],
    "prr": ["prr", "PRR"],
    "ebgm": ["ebgm", "EBGM"],
    "ebgm05": ["ebgm05", "EB05", "eb05"],
    "source": ["source", "dataset"],
}

POSITIVE_EVIDENCE_ALIASES = {
    "drug_id": ["drug_id", "drug", "drug_name", "_source_drug_name"],
    "adr_id": ["adr_id", "event", "adr", "meddra_term", "pt", "outcomeName", "outcome_name"],
}

TRUE_TEXT = {"1", "1.0", "true", "yes", "y", "positive", "active", "hit"}
FALSE_TEXT = {"0", "0.0", "false", "no", "n", "negative", "inactive", "non-hit", "nonhit"}
PUBCHEM_INACTIVE_CODES = {"1", "inactive"}
PUBCHEM_ACTIVE_CODES = {"2", "active", "5", "probe"}


def _source_name(path: str | Path, fallback: str) -> str:
    text = str(path).lower()
    for name in (
        "pubchem",
        "lit-pcba",
        "lit_pcba",
        "excape",
        "bindingdb",
        "davis",
        "kiba",
        "omop",
        "ohdsi",
        "faers",
    ):
        if name in text:
            return name.replace("_", "-")
    return fallback


def _clean_key(series: pd.Series) -> pd.Series:
    return series.astype("object").where(series.notna(), pd.NA).astype(str).str.strip()


def _activity_to_nm(df: pd.DataFrame) -> pd.Series:
    values = pd.to_numeric(df["activity_value_nM"], errors="coerce")
    units = df.get("activity_units", pd.Series("", index=df.index)).astype(str).str.lower()
    values = values.where(~units.isin({"um", "µm", "micromolar"}), values * 1000.0)
    values = values.where(~units.isin({"mm", "millimolar"}), values * 1_000_000.0)
    values = values.where(~units.isin({"pm", "picomolar"}), values / 1000.0)
    return values


def _outcome_masks(outcome: pd.Series, source_format: str) -> tuple[pd.Series, pd.Series]:
    text = outcome.astype("object").where(outcome.notna(), "").astype(str).str.strip().str.lower()
    if source_format == "pubchem_bioassay":
        return text.isin(PUBCHEM_ACTIVE_CODES), text.isin(PUBCHEM_INACTIVE_CODES)
    return text.isin(TRUE_TEXT), text.isin(FALSE_TEXT)


def normalize_measured_negative_source(
    path: str | Path,
    mapping_path: str | Path | None = None,
    *,
    source_format: str = "generic_bioactivity",
    inactive_threshold_nM: float = 10000.0,
    active_threshold_nM: float | None = None,
    confidence: float = 0.9,
) -> pd.DataFrame:
    """Normalize measured active/inactive sources into the negative evidence schema.

    This handles PubChem BioAssay, LIT-PCBA, ExCAPE-DB, BindingDB, Davis, KIBA,
    and generic active/inactive bioactivity tables through aliases rather than
    one parser per source.
    """

    threshold = float(active_threshold_nM if active_threshold_nM is not None else inactive_threshold_nM)
    raw = read_source_table(path)
    df = normalize_columns(raw, BIOACTIVITY_NEGATIVE_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    df["drug_id"] = _clean_key(df["drug_id"])
    df["target_id"] = _clean_key(df["target_id"])
    df["activity_value_nM"] = _activity_to_nm(df)
    relation = df.get("activity_relation", pd.Series("", index=df.index)).astype(str).str.strip()
    df["activity_relation"] = relation
    outcome_active, outcome_inactive = _outcome_masks(df["activity_outcome"], source_format)
    measured = df["activity_value_nM"].notna()
    exact = relation.eq("") | relation.eq("=") | relation.str.lower().eq("nan")
    upper_bound = relation.isin({"<", "<="})
    lower_bound = relation.str.startswith(">")
    numeric_active = (exact | upper_bound) & measured & df["activity_value_nM"].le(threshold)
    numeric_inactive = (exact | lower_bound) & measured & df["activity_value_nM"].gt(threshold)
    active = outcome_active | numeric_active
    inactive = outcome_inactive | numeric_inactive
    ambiguous = active & inactive
    label = pd.Series(pd.NA, index=df.index, dtype="Int64")
    label = label.mask(inactive & ~ambiguous, 0)
    label = label.mask(ambiguous, -1)
    status = pd.Series("unknown_unlabelable_activity", index=df.index)
    status = status.mask(active & ~inactive, "positive_activity_not_negative")
    status = status.mask(inactive & ~active, "measured_inactive")
    status = status.mask(ambiguous, "ambiguous_active_inactive")
    source = df["source"].fillna(_source_name(path, source_format)).astype(str)
    out = pd.DataFrame(
        {
            "drug_id": df["drug_id"],
            "target_id": df["target_id"],
            "adr_id": pd.NA,
            "evidence_scope": "drug_target",
            "mechanism_label": label,
            "negative_evidence_type": "measured_inactive_assay",
            "negative_source": source,
            "negative_confidence": confidence,
            "negative_label_status": status,
            "assay_id": df["assay_id"],
            "activity_type": df["activity_type"],
            "activity_relation": relation,
            "activity_value_nM": df["activity_value_nM"],
            "activity_threshold_nM": threshold,
            "publication_year": pd.to_numeric(df["publication_year"], errors="coerce"),
            "document_ids": df["document_ids"],
            "_source_file": str(path),
        }
    )
    return out[out["mechanism_label"].notna()].copy()


def normalize_omop_ohdsi_negative_source(
    path: str | Path,
    mapping_path: str | Path | None = None,
    *,
    confidence: float = 0.65,
) -> pd.DataFrame:
    raw = read_source_table(path)
    df = normalize_columns(raw, DRUG_ADR_NEGATIVE_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    truth = pd.to_numeric(df["ground_truth"], errors="coerce")
    text = df["ground_truth"].astype("object").where(df["ground_truth"].notna(), "").astype(str).str.lower()
    negative = truth.eq(0) | text.isin({"0", "false", "negative", "negative control"})
    if truth.isna().all() and "negative" in str(path).lower():
        negative = pd.Series(True, index=df.index)
    ambiguous = text.isin({"ambiguous", "conflict", "conflicting", "unknown"})
    label = pd.Series(pd.NA, index=df.index, dtype="Int64").mask(negative, 0).mask(ambiguous, -1)
    source = df["source"].fillna(_source_name(path, "omop_ohdsi")).astype(str)
    out = pd.DataFrame(
        {
            "drug_id": _clean_key(df["drug_id"]),
            "target_id": pd.NA,
            "adr_id": _clean_key(df["adr_id"]),
            "evidence_scope": "drug_adr",
            "mechanism_label": label,
            "negative_evidence_type": "omop_ohdsi_negative_control",
            "negative_source": source,
            "negative_confidence": confidence,
            "negative_label_status": pd.Series("drug_outcome_negative_control", index=df.index).where(
                label.eq(0), "unknown_or_positive_control"
            ),
            "_source_file": str(path),
        }
    )
    return out[out["mechanism_label"].notna()].copy()


def _positive_evidence_keys(path: str | Path | None) -> set[tuple[str, str]]:
    if path is None or not Path(path).exists():
        return set()
    df = normalize_columns(read_source_table(path), POSITIVE_EVIDENCE_ALIASES)
    drugs = _clean_key(df["drug_id"]).str.lower()
    adrs = _clean_key(df["adr_id"]).str.lower()
    return set(zip(drugs, adrs))


def normalize_faers_nonsignal_source(
    path: str | Path,
    mapping_path: str | Path | None = None,
    *,
    positive_evidence_path: str | Path | None = None,
    min_drug_reports: int = 100,
    min_event_reports: int = 100,
    min_pair_reports: int = 0,
    ror_threshold: float = 2.0,
    prr_threshold: float = 2.0,
    ebgm_threshold: float = 2.0,
    confidence: float = 0.4,
) -> pd.DataFrame:
    raw = read_source_table(path)
    df = normalize_columns(raw, FAERS_NONSIGNAL_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    for col in ("pair_reports", "drug_reports", "event_reports", "ror", "ror_lower", "prr", "ebgm", "ebgm05"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    observed = (
        df["drug_reports"].ge(min_drug_reports)
        & df["event_reports"].ge(min_event_reports)
        & df["pair_reports"].fillna(0).ge(min_pair_reports)
    )
    no_signal = (
        df["ror"].fillna(0).lt(ror_threshold)
        & df["ror_lower"].fillna(0).lt(ror_threshold)
        & df["prr"].fillna(0).lt(prr_threshold)
        & df["ebgm"].fillna(0).lt(ebgm_threshold)
        & df["ebgm05"].fillna(0).lt(ebgm_threshold)
    )
    positive_keys = _positive_evidence_keys(positive_evidence_path)
    keys = list(zip(_clean_key(df["drug_id"]).str.lower(), _clean_key(df["adr_id"]).str.lower()))
    has_positive_evidence = pd.Series([key in positive_keys for key in keys], index=df.index)
    negative = observed & no_signal & ~has_positive_evidence
    label = pd.Series(pd.NA, index=df.index, dtype="Int64").mask(negative, 0)
    source = df["source"].fillna(_source_name(path, "faers")).astype(str)
    out = pd.DataFrame(
        {
            "drug_id": _clean_key(df["drug_id"]),
            "target_id": pd.NA,
            "adr_id": _clean_key(df["adr_id"]),
            "evidence_scope": "drug_adr",
            "mechanism_label": label,
            "negative_evidence_type": "faers_observed_nonsignal",
            "negative_source": source,
            "negative_confidence": confidence,
            "negative_label_status": pd.Series("unknown_not_adequately_observed_or_signal", index=df.index).mask(
                negative, "observed_nonsignal_no_positive_evidence"
            ),
            "pair_reports": df["pair_reports"],
            "drug_reports": df["drug_reports"],
            "event_reports": df["event_reports"],
            "ror": df["ror"],
            "ror_lower": df["ror_lower"],
            "prr": df["prr"],
            "ebgm": df["ebgm"],
            "ebgm05": df["ebgm05"],
            "_source_file": str(path),
        }
    )
    return out[out["mechanism_label"].notna()].copy()


def summarize_negative_evidence(df: pd.DataFrame) -> dict[str, Any]:
    if df.empty:
        return {"rows": 0}
    summary: dict[str, Any] = {
        "rows": int(len(df)),
        "label_counts": {
            str(k): int(v) for k, v in df["mechanism_label"].value_counts(dropna=False).items()
        },
        "evidence_type_counts": {
            str(k): int(v) for k, v in df["negative_evidence_type"].value_counts(dropna=False).items()
        },
        "scope_counts": {str(k): int(v) for k, v in df["evidence_scope"].value_counts(dropna=False).items()},
    }
    for key in ("drug_id", "target_id", "adr_id"):
        if key in df.columns:
            summary[f"n_{key}s"] = int(df[key].dropna().astype(str).str.strip().replace("", pd.NA).dropna().nunique())
    return summary


def write_negative_evidence(df: pd.DataFrame, out_path: str | Path, manifest: dict[str, Any]) -> pd.DataFrame:
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    payload = {**manifest, "summary": summarize_negative_evidence(df)}
    out.with_suffix(".summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return df
