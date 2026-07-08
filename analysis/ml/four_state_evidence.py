from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from analysis.external.mapping import apply_pair_mapping, normalize_columns, read_optional_mapping
from analysis.external.source_tables import read_source_table


ACTIVITY_ALIASES = {
    "drug_id": [
        "drug_id",
        "compound_id",
        "chembl_id",
        "molecule_chembl_id",
        "bindingdb_monomerid",
        "monomerid",
        "cid",
        "pubchem_cid",
        "sid",
        "connectivity",
        "ligand_id",
        "gtopdb_ligand_id",
        "iuphar_ligand_id",
        "ligand",
        "ligand_name",
        "Ligand",
        "Ligand ID",
    ],
    "target_id": [
        "target_id",
        "uniprot",
        "uniprot_id",
        "target_uniprot",
        "Target UniProt ID",
        "uniprot_accession",
        "primary_uniprot",
        "target_chembl_id",
        "gtopdb_target_id",
        "iuphar_target_id",
        "accession",
        "protein_id",
        "gene",
        "gene_symbol",
        "Target Gene Symbol",
        "target",
        "target_name",
    ],
    "assay_id": ["assay_id", "aid", "pubchem_aid", "assay_chembl_id", "assay"],
    "activity_outcome": [
        "activity_outcome",
        "pubchem_activity_outcome",
        "outcome",
        "activity",
        "activity_class",
        "label",
        "hit_call",
        "toxcast_hit_call",
        "source_specific_activity_label",
        "curated_bioactivity_label",
        "toxcast_hts_label",
        "bioactivity_ml_label",
        "ml_binary_label",
        "class",
    ],
    "activity_value_nM": [
        "activity_value_nM",
        "activity_nM",
        "activity_nm",
        "standard_value",
        "value",
        "affinity",
        "ki",
        "kd",
        "ic50",
        "ec50",
        "ac50_nm",
        "ac50_nM",
        "affinity",
        "affinity_nm",
        "affinity_nM",
        "affinity_value",
        "value_nm",
        "value_nM",
        "Original Affinity Median nm",
        "Original Affinity Low nm",
        "Original Affinity High nm",
    ],
    "pchembl_value": [
        "pchembl_value",
        "pchembl_value_mean",
        "pchembl_value_Mean",
        "pxc50",
        "pXC50",
        "pki",
        "pKi",
        "pkd",
        "pKd",
        "pic50",
        "pIC50",
        "pec50",
        "pEC50",
        "pa2",
        "pA2",
        "Affinity Median",
        "Affinity Low",
        "Affinity High",
    ],
    "activity_relation": [
        "activity_relation",
        "standard_relation",
        "relation",
        "operator",
        "affinity_relation",
        "Original Affinity Relation",
    ],
    "activity_type": [
        "activity_type",
        "standard_type",
        "type",
        "endpoint",
        "affinity_type",
        "action",
        "ligand_action",
        "Original Affinity Units",
        "Type",
        "Action",
    ],
    "activity_units": ["activity_units", "standard_units", "units", "unit", "affinity_units"],
    "source": ["source", "source_name", "dataset", "label_source", "source_objective"],
    "publication_year": [
        "publication_year",
        "activity_publication_year",
        "document_year",
        "year",
        "Year",
        "reference_year",
    ],
    "document_ids": [
        "document_ids",
        "activity_document_ids",
        "document_chembl_id",
        "doc_id",
        "pmid",
        "pubmed_id",
        "pubmed_ids",
        "reference",
        "references",
        "PubMed ID",
    ],
    "data_validity_comment": [
        "data_validity_comment",
        "chembl_data_validity_comment",
        "standard_text_value",
        "comment",
    ],
    "potential_duplicate": ["potential_duplicate", "chembl_potential_duplicate"],
}

DRUG_ADR_ALIASES = {
    "drug_id": ["drug_id", "exposureName", "exposure_name", "drug", "drug_name", "targetName", "target_name"],
    "adr_id": ["adr_id", "outcomeName", "outcome_name", "event", "adr", "meddra_term", "pt"],
    "label": ["groundTruth", "ground_truth", "label", "negative_control", "is_negative"],
    "source": ["source", "dataset"],
    "confidence": ["confidence", "negative_confidence"],
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

TARGET_ADR_ALIASES = {
    "target_id": ["target_id", "uniprot", "uniprot_ac", "Uniprot AC"],
    "adr_id": ["adr_id", "adr", "adr_term", "ADR Term", "meddra_term"],
    "label": ["label_state", "target_adr_label", "mechanism_label", "label"],
    "source": ["source", "source_name", "dataset"],
    "confidence": ["confidence", "label_confidence"],
    "document_ids": ["document_ids", "pubmed_ids", "pmid", "pubmed_id"],
}

NEGATIVE_EVIDENCE_ALIASES = {
    "drug_id": ["drug_id"],
    "target_id": ["target_id"],
    "adr_id": ["adr_id"],
    "evidence_scope": ["evidence_scope"],
    "mechanism_label": ["mechanism_label"],
    "negative_evidence_type": ["negative_evidence_type"],
    "negative_source": ["negative_source"],
    "negative_confidence": ["negative_confidence"],
    "negative_label_status": ["negative_label_status"],
    "assay_id": ["assay_id"],
    "activity_value_nM": ["activity_value_nM"],
    "activity_relation": ["activity_relation"],
    "activity_type": ["activity_type"],
    "publication_year": ["publication_year"],
    "document_ids": ["document_ids"],
}

COLLAPSED_LABEL_COLUMNS = [
    "canonical_pair_key",
    "pair_type",
    "drug_id",
    "target_id",
    "adr_id",
    "four_state_label",
    "four_state_ml_label",
    "four_state_label_status",
    "evidence_sources",
    "parent_sources",
    "n_raw_evidence_rows",
    "n_production_evidence_rows",
    "n_benchmark_only_rows",
    "evidence_state_counts_json",
    "max_confidence",
]

ACTIVE_TEXT = {"1", "1.0", "true", "yes", "active", "hit", "positive", "2"}
INACTIVE_TEXT = {"0", "0.0", "false", "no", "inactive", "non-hit", "nonhit", "negative"}
INCONCLUSIVE_TEXT = {"inconclusive", "unspecified", "undetermined", "ambiguous", "probe"}
INVALID_TEXT = {"outside typical range", "potential transcription error", "non standard unit"}


def _clean(series: pd.Series) -> pd.Series:
    return series.astype("object").where(series.notna(), pd.NA).astype(str).str.strip()


def _norm_text(series: pd.Series) -> pd.Series:
    return _clean(series).str.lower()


def _to_nm(values: pd.Series, units: pd.Series) -> pd.Series:
    out = pd.to_numeric(values, errors="coerce")
    unit = units.astype("object").where(units.notna(), "").astype(str).str.lower()
    out = out.where(~unit.isin({"um", "µm", "micromolar"}), out * 1000.0)
    out = out.where(~unit.isin({"mm", "millimolar"}), out * 1_000_000.0)
    out = out.where(~unit.isin({"pm", "picomolar"}), out / 1000.0)
    return out


def _canonical_key(pair_type: str, drug: pd.Series, target: pd.Series, adr: pd.Series) -> pd.Series:
    drug_key = _norm_text(drug)
    target_key = _norm_text(target)
    adr_key = _norm_text(adr)
    if pair_type == "drug_target":
        return "drug_target|" + drug_key + "|" + target_key
    if pair_type == "drug_adr":
        return "drug_adr|" + drug_key + "|" + adr_key
    if pair_type == "target_adr":
        return "target_adr|" + target_key + "|" + adr_key
    raise ValueError(f"unsupported pair_type: {pair_type}")


def _provenance(df: pd.DataFrame, skip: Iterable[str]) -> pd.Series:
    skip_set = set(skip)

    def encode(row: pd.Series) -> str:
        payload = {
            str(key): value
            for key, value in row.items()
            if str(key) not in skip_set and pd.notna(value)
        }
        return json.dumps(payload, sort_keys=True)

    return df.apply(encode, axis=1)


def activity_source_to_evidence(
    path: str | Path,
    *,
    source_name: str,
    source_priority: int,
    source_role: str,
    mapping_path: str | Path | None = None,
    benchmark_only: bool = False,
    positive_threshold_nM: float = 1000.0,
    negative_threshold_nM: float = 10000.0,
    use_pchembl: bool = True,
    use_numeric_thresholds: bool = True,
    allow_numeric_negatives: bool = True,
    confidence: float = 0.9,
) -> pd.DataFrame:
    """Normalize direct assay/training/benchmark sources into raw four-state evidence rows."""

    raw = read_source_table(path)
    df = normalize_columns(raw, ACTIVITY_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    drug = _clean(df["drug_id"])
    target = _clean(df["target_id"])
    value_nm = _to_nm(df["activity_value_nM"], df["activity_units"])
    pchembl = pd.to_numeric(df["pchembl_value"], errors="coerce")
    relation = _clean(df["activity_relation"])
    relation_text = relation.str.lower()
    outcome = _norm_text(df["activity_outcome"])
    validity = _norm_text(df["data_validity_comment"])
    duplicate = _norm_text(df["potential_duplicate"]).isin({"1", "1.0", "true", "yes"})

    active_call = outcome.isin(ACTIVE_TEXT)
    inactive_call = outcome.isin(INACTIVE_TEXT)
    inconclusive = outcome.isin(INCONCLUSIVE_TEXT)
    invalid = duplicate | validity.apply(lambda text: any(flag in text for flag in INVALID_TEXT))

    exact = relation_text.isin({"", "=", "nan"})
    upper_bound = relation_text.isin({"<", "<="})
    lower_bound = relation_text.isin({">", ">="})

    numeric_positive = pd.Series(False, index=df.index)
    numeric_negative = pd.Series(False, index=df.index)
    gray_zone = pd.Series(False, index=df.index)
    if use_numeric_thresholds:
        numeric_positive = (exact | upper_bound) & value_nm.le(positive_threshold_nM)
        if allow_numeric_negatives:
            numeric_negative = (exact | lower_bound) & value_nm.ge(negative_threshold_nM)
        gray_zone = value_nm.gt(positive_threshold_nM) & value_nm.lt(negative_threshold_nM)
    if use_pchembl:
        numeric_positive |= pchembl.ge(6.0)
        if allow_numeric_negatives:
            numeric_negative |= pchembl.le(5.0)
        gray_zone |= pchembl.gt(5.0) & pchembl.lt(6.0)

    positive = active_call | numeric_positive
    negative = inactive_call | numeric_negative
    conflicting = positive & negative

    state = pd.Series(pd.NA, index=df.index, dtype="Int64")
    state = state.mask(positive & ~negative & ~invalid, 1)
    state = state.mask(negative & ~positive & ~invalid, 0)
    state = state.mask(gray_zone & ~positive & ~negative & ~invalid, pd.NA)
    state = state.mask(conflicting | invalid | inconclusive, -1)

    status = pd.Series("unknown_or_gray_zone", index=df.index)
    status = status.mask(state.eq(1).fillna(False), "strict_positive")
    status = status.mask(state.eq(0).fillna(False), "measured_or_reliable_negative")
    status = status.mask(inconclusive, "inconclusive")
    status = status.mask(invalid, "invalid_or_low_confidence")
    status = status.mask(conflicting, "conflicting_source_row")

    production_state = pd.Series(pd.NA, index=df.index, dtype="Int64") if benchmark_only else state
    out = pd.DataFrame(
        {
            "canonical_pair_key": _canonical_key("drug_target", drug, target, pd.Series(pd.NA, index=df.index)),
            "pair_type": "drug_target",
            "drug_id": drug,
            "target_id": target,
            "adr_id": pd.NA,
            "source_name": df["source"].fillna(source_name).astype(str),
            "parent_source": df["source"].fillna(source_name).astype(str),
            "source_priority": source_priority,
            "source_role": source_role,
            "evidence_namespace": "benchmark" if benchmark_only else "production",
            "raw_label_state": state,
            "production_label_state": production_state,
            "evidence_value_nM": value_nm,
            "evidence_relation": relation,
            "evidence_type": df["activity_type"],
            "assay_id": df["assay_id"],
            "publication_year": pd.to_numeric(df["publication_year"], errors="coerce"),
            "document_ids": df["document_ids"],
            "confidence": confidence,
            "label_status": status,
            "is_benchmark_only": benchmark_only,
            "_source_file": str(path),
        }
    )
    out["provenance_json"] = _provenance(df, skip=set())
    return out[out["drug_id"].ne("") & out["target_id"].ne("")].copy()


def spd_source_to_evidence(
    path: str | Path,
    *,
    mapping_path: str | Path | None = None,
    confidence: float = 0.95,
) -> pd.DataFrame:
    return activity_source_to_evidence(
        path,
        source_name="SPD",
        source_priority=1,
        source_role="direct_assay_exposure_panel",
        mapping_path=mapping_path,
        positive_threshold_nM=1000.0,
        negative_threshold_nM=10000.0,
        use_pchembl=False,
        confidence=confidence,
    )


def drug_adr_control_source_to_evidence(
    path: str | Path,
    *,
    source_name: str,
    source_priority: int,
    source_role: str,
    mapping_path: str | Path | None = None,
    confidence: float = 0.65,
) -> pd.DataFrame:
    raw = read_source_table(path)
    df = normalize_columns(raw, DRUG_ADR_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    text = _norm_text(df["label"])
    truth = pd.to_numeric(df["label"], errors="coerce")
    positive = truth.eq(1) | text.isin(ACTIVE_TEXT)
    negative = truth.eq(0) | text.isin(INACTIVE_TEXT)
    conflict = positive & negative
    state = pd.Series(pd.NA, index=df.index, dtype="Int64")
    state = state.mask(positive & ~conflict, 1)
    state = state.mask(negative & ~conflict, 0)
    state = state.mask(conflict, -1)
    drug = _clean(df["drug_id"])
    adr = _clean(df["adr_id"])
    out = pd.DataFrame(
        {
            "canonical_pair_key": _canonical_key("drug_adr", drug, pd.Series(pd.NA, index=df.index), adr),
            "pair_type": "drug_adr",
            "drug_id": drug,
            "target_id": pd.NA,
            "adr_id": adr,
            "source_name": df["source"].fillna(source_name).astype(str),
            "parent_source": df["source"].fillna(source_name).astype(str),
            "source_priority": source_priority,
            "source_role": source_role,
            "evidence_namespace": "production",
            "raw_label_state": state,
            "production_label_state": state,
            "evidence_value_nM": pd.NA,
            "evidence_relation": pd.NA,
            "evidence_type": "drug_adr_control",
            "assay_id": pd.NA,
            "publication_year": pd.NA,
            "document_ids": pd.NA,
            "confidence": pd.to_numeric(df["confidence"], errors="coerce").fillna(confidence),
            "label_status": pd.Series("drug_adr_control", index=df.index).mask(conflict, "conflicting_control_row"),
            "is_benchmark_only": False,
            "_source_file": str(path),
        }
    )
    out["provenance_json"] = _provenance(df, skip=set())
    return out[out["drug_id"].ne("") & out["adr_id"].ne("")].copy()


def drug_adr_positive_source_to_evidence(
    path: str | Path,
    *,
    source_name: str,
    source_priority: int,
    source_role: str,
    mapping_path: str | Path | None = None,
    confidence: float = 0.7,
) -> pd.DataFrame:
    """Normalize positive-only drug-ADR resources such as SIDER label rows.

    Positive-only clinical label sources can support known drug-ADR evidence,
    but they do not create negatives. Absence from these tables remains unknown.
    """

    raw = read_source_table(path)
    df = normalize_columns(raw, DRUG_ADR_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    drug = _clean(df["drug_id"])
    adr = _clean(df["adr_id"])
    state = pd.Series(1, index=df.index, dtype="Int64")
    out = pd.DataFrame(
        {
            "canonical_pair_key": _canonical_key("drug_adr", drug, pd.Series(pd.NA, index=df.index), adr),
            "pair_type": "drug_adr",
            "drug_id": drug,
            "target_id": pd.NA,
            "adr_id": adr,
            "source_name": source_name,
            "parent_source": df["source"].fillna(source_name).astype(str),
            "source_priority": source_priority,
            "source_role": source_role,
            "evidence_namespace": "production",
            "raw_label_state": state,
            "production_label_state": state,
            "evidence_value_nM": pd.NA,
            "evidence_relation": pd.NA,
            "evidence_type": "drug_adr_positive_label",
            "assay_id": pd.NA,
            "publication_year": pd.NA,
            "document_ids": pd.NA,
            "confidence": confidence,
            "label_status": "positive_only_drug_adr_label",
            "is_benchmark_only": False,
            "_source_file": str(path),
        }
    )
    out["provenance_json"] = _provenance(df, skip=set())
    return out[out["drug_id"].ne("") & out["adr_id"].ne("")].copy()


def faers_nonsignal_source_to_evidence(
    path: str | Path,
    *,
    mapping_path: str | Path | None = None,
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
    usable_negative = observed & no_signal
    unstable_or_signal = observed & ~no_signal
    state = pd.Series(pd.NA, index=df.index, dtype="Int64")
    state = state.mask(usable_negative, 0)
    state = state.mask(unstable_or_signal, -1)
    drug = _clean(df["drug_id"])
    adr = _clean(df["adr_id"])
    out = pd.DataFrame(
        {
            "canonical_pair_key": _canonical_key("drug_adr", drug, pd.Series(pd.NA, index=df.index), adr),
            "pair_type": "drug_adr",
            "drug_id": drug,
            "target_id": pd.NA,
            "adr_id": adr,
            "source_name": "FAERS non-signal",
            "parent_source": df["source"].fillna("FAERS").astype(str),
            "source_priority": 3,
            "source_role": "observed_faers_nonsignal_control",
            "evidence_namespace": "production",
            "raw_label_state": state,
            "production_label_state": state,
            "evidence_value_nM": pd.NA,
            "evidence_relation": pd.NA,
            "evidence_type": "drug_adr_faers_disproportionality",
            "assay_id": pd.NA,
            "publication_year": pd.NA,
            "document_ids": pd.NA,
            "confidence": confidence,
            "label_status": pd.Series("unknown_not_adequately_observed", index=df.index)
            .mask(usable_negative, "observed_nonsignal")
            .mask(unstable_or_signal, "signal_or_unstable_counts"),
            "is_benchmark_only": False,
            "_source_file": str(path),
            "pair_reports": df["pair_reports"],
            "drug_reports": df["drug_reports"],
            "event_reports": df["event_reports"],
            "ror": df["ror"],
            "ror_lower": df["ror_lower"],
            "prr": df["prr"],
            "ebgm": df["ebgm"],
            "ebgm05": df["ebgm05"],
        }
    )
    out["provenance_json"] = _provenance(df, skip=set())
    return out[out["drug_id"].ne("") & out["adr_id"].ne("")].copy()


def target_adr_source_to_evidence(
    path: str | Path,
    *,
    source_name: str,
    source_priority: int,
    source_role: str,
    confidence: float = 0.9,
) -> pd.DataFrame:
    """Normalize curated target-ADR positives into raw four-state evidence."""

    raw = read_source_table(path)
    df = normalize_columns(raw, TARGET_ADR_ALIASES)
    target = _clean(df["target_id"])
    adr = _clean(df["adr_id"])
    text = _norm_text(df["label"])
    numeric = pd.to_numeric(df["label"], errors="coerce")
    positive = numeric.eq(1) | text.isin(ACTIVE_TEXT) | text.eq("")
    ambiguous = numeric.eq(-1) | text.isin(INCONCLUSIVE_TEXT)
    state = pd.Series(pd.NA, index=df.index, dtype="Int64")
    state = state.mask(positive & ~ambiguous, 1)
    state = state.mask(ambiguous, -1)
    out = pd.DataFrame(
        {
            "canonical_pair_key": _canonical_key("target_adr", pd.Series(pd.NA, index=df.index), target, adr),
            "pair_type": "target_adr",
            "drug_id": pd.NA,
            "target_id": target,
            "adr_id": adr,
            "source_name": source_name,
            "parent_source": df["source"].fillna(source_name).astype(str),
            "source_priority": source_priority,
            "source_role": source_role,
            "evidence_namespace": "production",
            "raw_label_state": state,
            "production_label_state": state,
            "evidence_value_nM": pd.NA,
            "evidence_relation": pd.NA,
            "evidence_type": "target_adr_curated_positive",
            "assay_id": pd.NA,
            "publication_year": pd.NA,
            "document_ids": df["document_ids"],
            "confidence": pd.to_numeric(df["confidence"], errors="coerce").fillna(confidence),
            "label_status": pd.Series("strict_positive_target_adr", index=df.index).mask(
                ambiguous, "excluded_ambiguous_or_conflicting"
            ),
            "is_benchmark_only": False,
            "_source_file": str(path),
        }
    )
    out["provenance_json"] = _provenance(df, skip=set())
    return out[out["target_id"].ne("") & out["adr_id"].ne("")].copy()


def negative_evidence_to_raw(path: str | Path) -> pd.DataFrame:
    raw = read_source_table(path)
    df = normalize_columns(raw, NEGATIVE_EVIDENCE_ALIASES)
    state = pd.to_numeric(df["mechanism_label"], errors="coerce").astype("Int64")
    scope = df["evidence_scope"].fillna("drug_target").astype(str)
    drug = _clean(df["drug_id"])
    target = _clean(df["target_id"])
    adr = _clean(df["adr_id"])
    keys = pd.Series(pd.NA, index=df.index, dtype="object")
    keys = keys.mask(scope.eq("drug_target"), _canonical_key("drug_target", drug, target, adr))
    keys = keys.mask(scope.eq("drug_adr"), _canonical_key("drug_adr", drug, target, adr))
    out = pd.DataFrame(
        {
            "canonical_pair_key": keys,
            "pair_type": scope,
            "drug_id": drug,
            "target_id": target.mask(scope.ne("drug_target"), pd.NA),
            "adr_id": adr.mask(scope.ne("drug_adr"), pd.NA),
            "source_name": df["negative_source"].fillna("negative_evidence").astype(str),
            "parent_source": df["negative_source"].fillna("negative_evidence").astype(str),
            "source_priority": scope.map({"drug_target": 1, "drug_adr": 3}).fillna(4).astype(int),
            "source_role": df["negative_evidence_type"].fillna("negative_evidence").astype(str),
            "evidence_namespace": "production",
            "raw_label_state": state,
            "production_label_state": state,
            "evidence_value_nM": pd.to_numeric(df["activity_value_nM"], errors="coerce"),
            "evidence_relation": df["activity_relation"],
            "evidence_type": df["activity_type"].fillna(df["negative_evidence_type"]),
            "assay_id": df["assay_id"],
            "publication_year": pd.to_numeric(df["publication_year"], errors="coerce"),
            "document_ids": df["document_ids"],
            "confidence": pd.to_numeric(df["negative_confidence"], errors="coerce").fillna(0.5),
            "label_status": df["negative_label_status"],
            "is_benchmark_only": False,
            "_source_file": str(path),
        }
    )
    out["provenance_json"] = _provenance(df, skip=set())
    return out[out["canonical_pair_key"].notna()].copy()


def collapse_four_state_evidence(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame(columns=COLLAPSED_LABEL_COLUMNS)
    work = raw.copy()
    work["production_label_state"] = pd.to_numeric(work["production_label_state"], errors="coerce")
    work["raw_label_state"] = pd.to_numeric(work["raw_label_state"], errors="coerce")
    key_col = "canonical_pair_key"
    first = work.drop_duplicates(key_col, keep="first").set_index(key_col, drop=False)
    group = work.groupby(key_col, sort=False, dropna=False)
    out = pd.DataFrame(index=first.index)
    out["canonical_pair_key"] = first[key_col]
    for column in ("pair_type", "drug_id", "target_id", "adr_id"):
        out[column] = first[column]
    out["n_raw_evidence_rows"] = group.size()
    out["n_production_evidence_rows"] = work["production_label_state"].notna().groupby(work[key_col], sort=False).sum()
    out["n_benchmark_only_rows"] = work["is_benchmark_only"].astype(bool).groupby(work[key_col], sort=False).sum()
    out["max_confidence"] = pd.to_numeric(work["confidence"], errors="coerce").groupby(work[key_col], sort=False).max()

    out["evidence_sources"] = first["source_name"].fillna("").astype(str)
    out["parent_sources"] = first["parent_source"].fillna("").astype(str)
    duplicate_keys = out.index[out["n_raw_evidence_rows"].gt(1)]
    if len(duplicate_keys) > 0:
        duplicate_work = work[work[key_col].isin(duplicate_keys)]
        out.loc[duplicate_keys, "evidence_sources"] = duplicate_work.groupby(key_col, sort=False)["source_name"].agg(
            lambda values: ";".join(sorted(set(values.dropna().astype(str))))
        )
        out.loc[duplicate_keys, "parent_sources"] = duplicate_work.groupby(key_col, sort=False)["parent_source"].agg(
            lambda values: ";".join(sorted(set(values.dropna().astype(str))))
        )

    out["evidence_state_counts_json"] = "{}"
    raw_counts = work.groupby([key_col, "raw_label_state"], sort=False, dropna=False).size()
    if not raw_counts.empty:
        state_frames = raw_counts.rename("n").reset_index()
        state_json = state_frames.groupby(key_col, sort=False).apply(
            lambda frame: json.dumps(
                {str(row["raw_label_state"]): int(row["n"]) for _, row in frame.iterrows()},
                sort_keys=True,
            ),
            include_groups=False,
        )
        out.loc[state_json.index, "evidence_state_counts_json"] = state_json

    production = work[work["production_label_state"].notna()].copy()
    out["four_state_label"] = pd.NA
    out["four_state_label_status"] = "unknown_or_no_production_truth"
    out.loc[out["n_benchmark_only_rows"].gt(0) & out["n_production_evidence_rows"].eq(0), "four_state_label_status"] = (
        "benchmark_only_not_production_truth"
    )
    if not production.empty:
        prod_key = production[key_col]
        invalid_any = production["production_label_state"].eq(-1).groupby(prod_key, sort=False).any()
        pos = production[production["production_label_state"].eq(1)]
        neg = production[production["production_label_state"].eq(0)]
        pos_any = pos.groupby(key_col, sort=False).size().gt(0) if not pos.empty else pd.Series(dtype=bool)
        neg_any = neg.groupby(key_col, sort=False).size().gt(0) if not neg.empty else pd.Series(dtype=bool)
        min_pos = pos.groupby(key_col, sort=False)["source_priority"].min() if not pos.empty else pd.Series(dtype=float)
        min_neg = neg.groupby(key_col, sort=False)["source_priority"].min() if not neg.empty else pd.Series(dtype=float)

        all_prod_keys = out.index[out["n_production_evidence_rows"].gt(0)]
        invalid = invalid_any.reindex(all_prod_keys, fill_value=False)
        has_pos = pos_any.reindex(all_prod_keys, fill_value=False)
        has_neg = neg_any.reindex(all_prod_keys, fill_value=False)
        pos_priority = min_pos.reindex(all_prod_keys).fillna(999)
        neg_priority = min_neg.reindex(all_prod_keys).fillna(999)

        invalid_only = invalid & ~has_pos & ~has_neg
        out.loc[invalid_only[invalid_only].index, ["four_state_label", "four_state_label_status"]] = [
            -1,
            "excluded_ambiguous_or_conflicting",
        ]
        same_priority_conflict = has_pos & has_neg & pos_priority.eq(neg_priority)
        out.loc[same_priority_conflict[same_priority_conflict].index, ["four_state_label", "four_state_label_status"]] = [
            -1,
            "conflicting_comparable_measured_evidence",
        ]
        pos_overrides = has_pos & has_neg & pos_priority.lt(neg_priority)
        out.loc[pos_overrides[pos_overrides].index, ["four_state_label", "four_state_label_status"]] = [
            1,
            "positive_overrides_lower_priority_negative",
        ]
        neg_overrides = has_pos & has_neg & neg_priority.lt(pos_priority)
        out.loc[neg_overrides[neg_overrides].index, ["four_state_label", "four_state_label_status"]] = [
            0,
            "negative_overrides_lower_priority_positive",
        ]
        positive_only = has_pos & ~has_neg
        out.loc[positive_only[positive_only].index, ["four_state_label", "four_state_label_status"]] = [
            1,
            "strict_positive",
        ]
        negative_only = ~has_pos & has_neg
        out.loc[negative_only[negative_only].index, ["four_state_label", "four_state_label_status"]] = [
            0,
            "measured_or_reliable_negative",
        ]

    out["four_state_ml_label"] = out["four_state_label"].where(~out["four_state_label"].isin([-1]), pd.NA)
    out = out.reset_index(drop=True)[COLLAPSED_LABEL_COLUMNS]
    if not out.empty:
        out["four_state_label"] = pd.to_numeric(out["four_state_label"], errors="coerce").astype("Int64")
        out["four_state_ml_label"] = pd.to_numeric(out["four_state_ml_label"], errors="coerce").astype("Int64")
    return out


def join_collapsed_labels(
    feature_table: pd.DataFrame,
    collapsed: pd.DataFrame,
    *,
    pair_type: str = "drug_target",
) -> pd.DataFrame:
    out = feature_table.copy()
    if pair_type == "drug_target":
        key = _canonical_key("drug_target", _clean(out["drug_id"]), _clean(out["target_id"]), pd.Series(pd.NA, index=out.index))
    elif pair_type == "drug_adr":
        adr_col = "adr_id" if "adr_id" in out.columns else "adr_term"
        key = _canonical_key("drug_adr", _clean(out["drug_id"]), pd.Series(pd.NA, index=out.index), _clean(out[adr_col]))
    else:
        raise ValueError(f"unsupported pair_type: {pair_type}")
    out["canonical_pair_key"] = key
    replace_columns = [
        column
        for column in COLLAPSED_LABEL_COLUMNS
        if column not in {"canonical_pair_key", "pair_type", "drug_id", "target_id", "adr_id"}
    ]
    out = out.drop(columns=[column for column in replace_columns if column in out.columns], errors="ignore")
    if collapsed.empty:
        for column in replace_columns:
            if column not in out.columns:
                out[column] = pd.NA
        return out
    keep = [
        "canonical_pair_key",
        "four_state_label",
        "four_state_ml_label",
        "four_state_label_status",
        "evidence_sources",
        "parent_sources",
        "n_raw_evidence_rows",
        "n_production_evidence_rows",
        "n_benchmark_only_rows",
        "evidence_state_counts_json",
        "max_confidence",
    ]
    return out.merge(collapsed[keep], on="canonical_pair_key", how="left")


def write_four_state_outputs(
    raw: pd.DataFrame,
    collapsed: pd.DataFrame,
    out_dir: str | Path,
    *,
    joined: pd.DataFrame | None = None,
) -> dict[str, Any]:
    path = Path(out_dir)
    path.mkdir(parents=True, exist_ok=True)
    raw_path = path / "four_state_raw_evidence.csv"
    collapsed_path = path / "four_state_collapsed_labels.csv"
    raw.to_csv(raw_path, index=False)
    collapsed.to_csv(collapsed_path, index=False)
    manifest: dict[str, Any] = {
        "raw_evidence_rows": int(len(raw)),
        "collapsed_label_rows": int(len(collapsed)),
        "raw_sources": {str(k): int(v) for k, v in raw["source_name"].value_counts(dropna=False).items()} if not raw.empty else {},
        "collapsed_status_counts": {
            str(k): int(v) for k, v in collapsed["four_state_label_status"].value_counts(dropna=False).items()
        }
        if not collapsed.empty
        else {},
        "label_counts": {
            str(k): int(v) for k, v in collapsed["four_state_label"].value_counts(dropna=False).items()
        }
        if not collapsed.empty
        else {},
    }
    if joined is not None:
        joined_path = path / "four_state_joined_feature_table.csv"
        joined.to_csv(joined_path, index=False)
        manifest["joined_feature_table_rows"] = int(len(joined))
        manifest["joined_labelable_rows"] = (
            int(joined["four_state_ml_label"].notna().sum()) if "four_state_ml_label" in joined.columns else 0
        )
    (path / "four_state_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
