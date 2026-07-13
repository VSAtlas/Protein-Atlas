from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable

import pandas as pd

from analysis.external.mapping import normalize_columns
from analysis.external.source_tables import read_source_table


PK_CONTEXT_COLUMNS = [
    "drug_id", "drug_name", "inchikey", "source_name", "source_version",
    "source_record_id", "source_url", "study_id", "reference", "population",
    "species", "dose_value", "dose_unit", "dose_text", "route", "regimen",
    "formulation", "steady_state", "parent_or_metabolite", "cmax_value_raw",
    "cmax_unit_raw", "cmax_um", "protein_binding_percent",
    "fraction_unbound_plasma", "free_cmax_um", "free_cmax_method",
    "clearance_value", "clearance_unit", "bioavailability_value",
    "bioavailability_unit", "extraction_method", "source_confidence",
    "context_status", "license_note", "missing_reason",
]

ALIASES = {
    "drug_id": ["drug_id", "atlas_drug_id", "compound_id", "drugbank_id", "substance_id", "id"],
    "drug_name": ["drug_name", "display_name", "generic_name", "compound_name", "substance", "name"],
    "inchikey": ["inchikey", "inchi_key", "drugcentral_inchikey"],
    "source_record_id": ["source_record_id", "record_id", "row_id", "pk_id", "id"],
    "study_id": ["study_id", "study", "study_sid", "trial_id"],
    "reference": ["reference", "citation", "publication", "source_reference", "pmid", "doi"],
    "population": ["population", "subject_population", "group", "cohort", "age_group"],
    "species": ["species", "organism"],
    "dose_value": ["dose_value", "normalized_dose_value", "dose", "dose_amount"],
    "dose_unit": ["dose_unit", "normalized_dose_unit"],
    "dose_text": ["dose_text", "dosage", "dose_description"],
    "route": ["route", "route_of_administration", "dosage_route", "application"],
    "regimen": ["regimen", "dosing_regimen", "schedule", "frequency"],
    "formulation": ["formulation", "dosage_form", "form"],
    "steady_state": ["steady_state", "is_steady_state"],
    "parent_or_metabolite": ["parent_or_metabolite", "analyte_type", "cmax_from_parent"],
    "cmax_value_raw": ["cmax_value_raw", "cmax", "value", "normalized_value"],
    "cmax_unit_raw": ["cmax_unit_raw", "cmax_unit", "unit", "normalized_unit"],
    "cmax_um": ["cmax_um", "cmax_uM", "total_cmax_um", "combined_cmax_um"],
    "protein_binding_percent": ["protein_binding_percent", "protein_binding", "ppb_percent", "ppb", "bound_percent"],
    "fraction_unbound_plasma": ["fraction_unbound_plasma", "fraction_unbound", "fu", "fu_plasma"],
    "free_cmax_um": ["free_cmax_um", "free_cmax_uM", "cmax_free_um", "combined_free_cmax_um"],
    "clearance_value": ["clearance_value", "clearance", "cl"],
    "clearance_unit": ["clearance_unit", "cl_unit"],
    "bioavailability_value": ["bioavailability_value", "bioavailability", "absolute_bioavailability"],
    "bioavailability_unit": ["bioavailability_unit"],
}

SOURCE_PRIORITY = {
    "SPD": 0,
    "existing_phase1": 1,
    "NCATS_Inxight_FRDB": 2,
    "PK-DB": 3,
    "DrugBank_Cmax": 4,
    "DrugBank_protein_binding": 4,
    "DailyMed_openFDA_SPL": 5,
    "existing_pk": 6,
}

NUMERIC_COLUMNS = [
    "dose_value", "cmax_value_raw", "cmax_um", "protein_binding_percent",
    "fraction_unbound_plasma", "free_cmax_um", "clearance_value",
    "bioavailability_value",
]


def normalize_key(value: Any) -> str:
    if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
        return ""
    text = str(value).strip().casefold()
    if text in {"", "nan", "none", "null", "unknown"}:
        return ""
    return " ".join(re.sub(r"[^a-z0-9]+", " ", text).split())


def _text(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame.get(column, pd.Series("", index=frame.index)).fillna("").astype(str).str.strip()


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame.get(column, pd.Series(pd.NA, index=frame.index)), errors="coerce")


def _numeric_coverage(values: pd.Series) -> int:
    return int(pd.to_numeric(values, errors="coerce").notna().sum())


def _text_coverage(values: pd.Series) -> int:
    normalized = values.fillna("").astype(str).str.strip()
    return int(normalized.ne("").sum())


def _joined_text(*values: pd.Series, separator: str = "; ") -> pd.Series:
    if not values:
        return pd.Series(dtype="object")
    frame = pd.concat(
        [value.fillna("").astype(str).str.strip() for value in values], axis=1
    )
    return frame.apply(
        lambda row: separator.join(value for value in row if value), axis=1
    )


def _first_numeric(values: pd.Series) -> pd.Series:
    extracted = (
        values.fillna("")
        .astype(str)
        .str.extract(
            r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
            expand=False,
        )
    )
    return pd.to_numeric(extracted, errors="coerce")


def _empty(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame({column: pd.Series(pd.NA, index=index, dtype="object") for column in PK_CONTEXT_COLUMNS})


def _context_hash(row: pd.Series) -> str:
    fields = [
        row.get("source_name"), row.get("source_record_id"), row.get("study_id"),
        row.get("drug_id"), row.get("drug_name"), row.get("dose_value"),
        row.get("dose_unit"), row.get("route"), row.get("formulation"),
        row.get("cmax_um"),
    ]
    payload = "|".join("" if pd.isna(value) else str(value) for value in fields)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def finalize_context(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.reindex(columns=PK_CONTEXT_COLUMNS).copy()
    for column in NUMERIC_COLUMNS:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    derived_fu = 1.0 - out["protein_binding_percent"] / 100.0
    out["fraction_unbound_plasma"] = out["fraction_unbound_plasma"].fillna(
        derived_fu.where(derived_fu.between(0.0, 1.0, inclusive="both"))
    )
    reported = out["free_cmax_um"].notna()
    out["free_cmax_um"] = out["free_cmax_um"].fillna(
        out["cmax_um"] * out["fraction_unbound_plasma"]
    )
    method = pd.Series(pd.NA, index=out.index, dtype="object")
    method.loc[reported] = "reported"
    method.loc[~reported & out["free_cmax_um"].notna()] = "derived_cmax_times_fu"
    out["free_cmax_method"] = out["free_cmax_method"].fillna(method)
    valid_id = _text(out, "drug_id").str.len().gt(0) | _text(out, "drug_name").str.len().gt(0)
    out = out.loc[valid_id].copy()
    missing = out["cmax_um"].isna() & out["free_cmax_um"].isna()
    out.loc[missing & out["missing_reason"].isna(), "missing_reason"] = "no_cmax_or_free_cmax"
    out.insert(0, "pk_context_id", out.apply(_context_hash, axis=1))
    return out.drop_duplicates("pk_context_id", keep="first").reset_index(drop=True)


def load_spd_pk_context(path: str | Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="S Data 2", header=2)
    out = _empty(raw.index)
    out["drug_id"] = _text(raw, "drugcentral struct id")
    out["drug_name"] = _text(raw, "drugcentral name")
    out["inchikey"] = _text(raw, "drugcentral inchikey")
    out["source_name"] = "SPD"
    out["source_version"] = "Sutherland_2023_Supplementary_Data_2"
    out["source_record_id"] = out["drug_id"]
    out["source_url"] = "https://doi.org/10.1038/s41467-023-40064-9"
    out["cmax_value_raw"] = _numeric(raw, "Cmax tot (uM)")
    out["cmax_unit_raw"] = "uM"
    out["cmax_um"] = _numeric(raw, "Cmax tot (uM)")
    out["protein_binding_percent"] = _numeric(raw, "PPB %")
    out["free_cmax_um"] = _numeric(raw, "Cmax, free (uM)")
    out["parent_or_metabolite"] = _text(raw, "Cmax from parent?")
    out["reference"] = (
        "Cmax: " + _text(raw, "Cmax source") + "; PPB: " + _text(raw, "PPB source")
    ).str.strip("; ")
    out["species"] = "Homo sapiens"
    out["extraction_method"] = "structured_supplement"
    out["source_confidence"] = "high"
    out["context_status"] = "clinical_context_not_preserved_in_spd_supplement"
    out["license_note"] = "CC-BY-4.0 article supplement; upstream source restrictions may apply"
    return finalize_context(out)


def load_existing_phase1_pk_context(model_table: pd.DataFrame) -> pd.DataFrame:
    raw = model_table.drop_duplicates("drug_id").copy()
    out = _empty(raw.index)
    out["drug_id"] = _text(raw, "drug_id")
    out["drug_name"] = _text(raw, "generic_name").where(
        _text(raw, "generic_name").str.len().gt(0),
        _text(raw, "display_name"),
    )
    out["inchikey"] = _text(raw, "inchikey")
    out["source_name"] = "existing_phase1"
    out["source_version"] = "AtlasSPD_phase1_pre_context_join"
    out["source_record_id"] = out["drug_id"]
    out["cmax_value_raw"] = _numeric(raw, "cmax_um")
    out["cmax_unit_raw"] = "uM"
    out["cmax_um"] = _numeric(raw, "cmax_um")
    out["fraction_unbound_plasma"] = _numeric(raw, "fraction_unbound_plasma")
    out["free_cmax_um"] = _numeric(raw, "free_cmax_um")
    out["reference"] = _text(raw, "spd_source")
    out["species"] = "Homo sapiens"
    out["extraction_method"] = "existing_phase1_join"
    out["source_confidence"] = "high"
    out["context_status"] = "legacy_value_preserved; clinical_context_may_be_missing"
    out["license_note"] = "inherits source restrictions from existing Phase 1 table"
    return finalize_context(out)


def _frdb_cmax_um(frame: pd.DataFrame) -> pd.Series:
    value = _numeric(frame, "pk_cmax_value")
    analyte_mw = _numeric(frame, "pk_analyte_mw")
    molecular_weight = analyte_mw.where(
        analyte_mw.gt(0), _numeric(frame, "pk_application_mw")
    )
    unit = (
        _text(frame, "pk_cmax_units")
        .str.casefold()
        .str.replace("μ", "u", regex=False)
        .str.replace("µ", "u", regex=False)
        .str.replace(" ", "", regex=False)
    )
    out = pd.Series(float("nan"), index=frame.index, dtype="float64")

    for normalized_unit, factor in {
        "pm": 1e-6,
        "nm": 1e-3,
        "um": 1.0,
        "mm": 1e3,
        "nmol/ml": 1.0,
        "pmol/ml": 1e-3,
    }.items():
        mask = unit.eq(normalized_unit)
        out.loc[mask] = value.loc[mask] * factor

    valid_mw = molecular_weight.gt(0)
    for normalized_unit, factor in {
        "ng/l": 1e-3,
        "ng/dl": 1e-2,
        "pg/ml": 1e-3,
        "ng/ml": 1.0,
        "ug/l": 1.0,
        "ug/dl": 10.0,
        "ug/ml": 1e3,
        "mg/l": 1e3,
        "mg/dl": 1e4,
        "mg/ml": 1e6,
    }.items():
        mask = unit.eq(normalized_unit) & valid_mw
        out.loc[mask] = value.loc[mask] * factor / molecular_weight.loc[mask]
    return out


def load_ncats_frdb_pk_context(
    path: str | Path,
    *,
    source_version: str = "2024-12-30",
) -> pd.DataFrame:
    """Normalize NCATS FRDB without promoting free-text PK comments."""

    raw = read_source_table(path)
    out = _empty(raw.index)
    application_name = _text(raw, "pk_application_pt")
    analyte_name = _text(raw, "pk_analyte_pt")
    compound_id = _text(raw, "compound_id")
    out["drug_id"] = compound_id.map(
        lambda value: f"frdb:{value}" if value else ""
    )
    out["drug_name"] = application_name.where(
        application_name.str.len().gt(0), analyte_name
    )
    out["source_name"] = "NCATS_Inxight_FRDB"
    out["source_version"] = source_version
    out["source_record_id"] = _text(raw, "id")
    out["source_url"] = _text(raw, "pk_source_uri")
    out["reference"] = _text(raw, "pk_source_uri")
    out["population"] = _joined_text(
        _text(raw, "pk_age_group"), _text(raw, "pk_health_status")
    )
    out["species"] = _text(raw, "pk_species")
    out["dose_value"] = _numeric(raw, "pk_dose_value")
    out["dose_unit"] = _text(raw, "pk_dose_units").where(
        _text(raw, "pk_dose_units").str.len().gt(0),
        _text(raw, "pk_dose_other_units"),
    )
    out["dose_text"] = _text(raw, "pk_dose_type").where(
        _text(raw, "pk_dose_type").str.len().gt(0),
        _text(raw, "pk_dose_other_types"),
    )
    out["route"] = _text(raw, "pk_routes").where(
        _text(raw, "pk_routes").str.len().gt(0),
        _text(raw, "pk_adm_other_routes"),
    )
    frequency = _joined_text(
        _text(raw, "pk_frequency_times"),
        _text(raw, "pk_frequency_period"),
        separator=" per ",
    )
    out["regimen"] = _joined_text(
        _text(raw, "pk_experiment_type"), frequency
    )
    out["formulation"] = _text(raw, "pk_annotated_form")
    out["steady_state"] = _text(raw, "pk_experiment_type").str.contains(
        "steady", case=False, na=False
    )
    application_unii = _text(raw, "pk_application_unii")
    analyte_unii = _text(raw, "pk_analyte_unii")
    out["parent_or_metabolite"] = "unknown"
    same_analyte = application_unii.ne("") & application_unii.eq(analyte_unii)
    distinct_analyte = application_unii.ne("") & analyte_unii.ne("") & ~same_analyte
    out.loc[same_analyte, "parent_or_metabolite"] = "parent"
    out.loc[distinct_analyte, "parent_or_metabolite"] = "analyte_or_metabolite"
    out["cmax_value_raw"] = _numeric(raw, "pk_cmax_value")
    out["cmax_unit_raw"] = _text(raw, "pk_cmax_units")
    out["cmax_um"] = _frdb_cmax_um(raw)
    unsupported_cmax = out["cmax_value_raw"].notna() & out["cmax_um"].isna()
    fraction_unbound_percent = _numeric(raw, "pk_funbound_value")
    valid_fraction = fraction_unbound_percent.between(
        0.0, 100.0, inclusive="both"
    )
    out["fraction_unbound_plasma"] = (fraction_unbound_percent / 100.0).where(
        valid_fraction
    )
    out["protein_binding_percent"] = (100.0 - fraction_unbound_percent).where(
        valid_fraction
    )
    out["extraction_method"] = "structured_ncats_frdb_table"
    out["source_confidence"] = "medium"
    out["context_status"] = (
        "source_context_preserved; clearance_and_bioavailability_comments_not_promoted"
    )
    out.loc[unsupported_cmax, "context_status"] += (
        "; cmax_unit_not_safely_convertible"
    )
    out["license_note"] = "public NCATS FRDB release; retain source citation"
    return finalize_context(out)


def load_reviewed_openfda_pk_context(path: str | Path) -> pd.DataFrame:
    """Load only source-text-adjudicated DailyMed/openFDA PK scenarios."""

    raw = read_source_table(path)
    accepted = raw.get(
        "acceptable_for_model_training", pd.Series(False, index=raw.index)
    )
    accepted = accepted.fillna(False).astype(str).str.casefold().isin(
        {"1", "true", "yes", "y"}
    )
    accepted &= _text(raw, "adjudication_status").eq("accept_model_context")
    raw = raw.loc[accepted].copy()
    out = _empty(raw.index)
    out["drug_id"] = _text(raw, "drug_id")
    out["drug_name"] = _text(raw, "drug_id")
    out["source_name"] = "DailyMed_openFDA_SPL"
    out["source_version"] = _text(raw, "spl_version")
    out["source_record_id"] = _text(raw, "source_record_id")
    out["source_url"] = _text(raw, "source_url")
    out["reference"] = _text(raw, "source_url")
    out["dose_value"] = _numeric(raw, "adjudicated_dose_value").fillna(
        _numeric(raw, "dose_value")
    )
    out["dose_unit"] = _text(raw, "adjudicated_dose_unit").where(
        _text(raw, "adjudicated_dose_unit").ne(""), _text(raw, "dose_unit")
    )
    out["route"] = _text(raw, "adjudicated_route").where(
        _text(raw, "adjudicated_route").ne(""), _text(raw, "route")
    )
    out["regimen"] = _text(raw, "adjudicated_regimen").where(
        _text(raw, "adjudicated_regimen").ne(""), _text(raw, "regimen")
    )
    out["steady_state"] = _text(raw, "steady_state").str.casefold().isin(
        {"1", "true", "yes", "y"}
    )
    out["parent_or_metabolite"] = _text(raw, "adjudicated_context")
    out["cmax_value_raw"] = _first_numeric(_text(raw, "cmax_values_raw"))
    out["cmax_unit_raw"] = _text(raw, "cmax_units_raw")
    out["cmax_um"] = _first_numeric(
        _text(raw, "cmax_converted_um_candidates")
    )
    out["protein_binding_percent"] = _first_numeric(
        _text(raw, "protein_binding_values_pct")
    )
    out["extraction_method"] = "manual_source_text_adjudication"
    out["source_confidence"] = "high"
    out["context_status"] = (
        "source_text_adjudicated; contextual_scenario_not_universal_drug_pk"
    )
    out["license_note"] = (
        "public FDA labeling; retain SPL set/version provenance"
    )
    return finalize_context(out)


def load_flat_pk_context(
    path: str | Path,
    *,
    source_name: str,
    source_version: str = "",
) -> pd.DataFrame:
    frame = normalize_columns(read_source_table(path), ALIASES)
    out = _empty(frame.index)
    for column in ALIASES:
        out[column] = frame[column]
    out["source_name"] = source_name
    out["source_version"] = source_version
    out["extraction_method"] = "structured_table"
    out["source_confidence"] = "medium"
    out["context_status"] = "source_context_preserved_where_available"
    return finalize_context(out)


def combine_pk_context(parts: Iterable[pd.DataFrame]) -> pd.DataFrame:
    usable = [part for part in parts if part is not None and not part.empty]
    if not usable:
        return pd.DataFrame(columns=["pk_context_id", *PK_CONTEXT_COLUMNS])
    raw = pd.concat(usable, ignore_index=True).drop(columns=["pk_context_id"], errors="ignore")
    return finalize_context(raw)


def select_representative_pk_context(context: pd.DataFrame) -> pd.DataFrame:
    if context.empty:
        return context.copy()
    ranked = context.copy()
    unreviewed_spl = ranked["source_name"].eq("DailyMed_openFDA_SPL") & ranked[
        "source_confidence"
    ].eq("low")
    numeric_columns = [
        "cmax_value_raw",
        "cmax_unit_raw",
        "cmax_um",
        "protein_binding_percent",
        "fraction_unbound_plasma",
        "free_cmax_um",
        "free_cmax_method",
    ]
    ranked.loc[unreviewed_spl, numeric_columns] = pd.NA
    ranked.loc[unreviewed_spl, "context_status"] = (
        ranked.loc[unreviewed_spl, "context_status"].fillna("").astype(str)
        + "; numeric_quarantined_pending_source_text_review"
    ).str.lstrip("; ")
    ranked["_source_rank"] = ranked["source_name"].map(SOURCE_PRIORITY).fillna(99)
    ranked["_measurement_rank"] = 3
    ranked.loc[
        pd.to_numeric(ranked["dose_value"], errors="coerce").notna(),
        "_measurement_rank",
    ] = 2
    ranked.loc[ranked["cmax_um"].notna(), "_measurement_rank"] = 1
    ranked.loc[ranked["free_cmax_um"].notna(), "_measurement_rank"] = 0
    ranked["_confidence_rank"] = ranked["source_confidence"].map(
        {"high": 0, "medium": 1, "low": 2}
    ).fillna(3)
    context_present = pd.DataFrame(
        {
            "dose": pd.to_numeric(
                ranked["dose_value"], errors="coerce"
            ).notna(),
            "route": _text(ranked, "route").ne(""),
            "formulation": _text(ranked, "formulation").ne(""),
            "regimen": _text(ranked, "regimen").ne(""),
        },
        index=ranked.index,
    )
    ranked["_context_rank"] = context_present.sum(axis=1) * -1
    ranked["_drug_key"] = _text(ranked, "drug_name").map(normalize_key)
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(ranked.loc[missing], "inchikey").str.upper()
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(ranked.loc[missing], "drug_id").map(normalize_key)
    ranked = ranked.sort_values(
        [
            "_measurement_rank",
            "_source_rank",
            "_confidence_rank",
            "_context_rank",
            "pk_context_id",
        ],
        kind="stable",
    )
    out = ranked.drop_duplicates("_drug_key", keep="first").copy()
    return out.drop(columns=[column for column in out if column.startswith("_")]).reset_index(drop=True)


def _mapping_keys(row: pd.Series) -> set[str]:
    keys: set[str] = set()
    raw_inchikey = row.get("inchikey")
    inchikey = (
        ""
        if raw_inchikey is None or pd.isna(raw_inchikey)
        else str(raw_inchikey).strip().upper()
    )
    if inchikey and inchikey not in {"NAN", "NONE"}:
        keys.add(f"inchikey:{inchikey}")
        keys.add(f"inchikey14:{inchikey[:14]}")
    for column in ("drug_id", "drug_name", "display_name", "generic_name", "ligand_base"):
        value = normalize_key(row.get(column))
        if value:
            keys.add(f"name:{value}")
    return keys


def join_representative_pk_context(
    model_table: pd.DataFrame,
    representative: pd.DataFrame,
) -> pd.DataFrame:
    out = model_table.reset_index(drop=True).copy()
    stale_context = [
        column for column in out.columns if column.startswith("pk_context_")
    ]
    if stale_context:
        out = out.drop(columns=stale_context)
    if representative.empty:
        out["pk_context_join_status"] = "no_context_sources_available"
        return out
    lookup: dict[str, set[int]] = {}
    for idx, row in representative.iterrows():
        for key in _mapping_keys(row):
            lookup.setdefault(key, set()).add(idx)

    def best_candidate(candidates: set[int]) -> int:
        return min(
            candidates,
            key=lambda idx: (
                pd.isna(representative.loc[idx, "free_cmax_um"]),
                pd.isna(representative.loc[idx, "cmax_um"]),
                pd.isna(representative.loc[idx, "dose_value"]),
                SOURCE_PRIORITY.get(str(representative.loc[idx, "source_name"]), 99),
                idx,
            ),
        )

    selected: list[int | None] = []
    methods: list[str] = []
    for _, row in out.iterrows():
        keys = _mapping_keys(row)
        choice: int | None = None
        method = "unmatched"
        exact_identity = set().union(
            *(lookup.get(key, set()) for key in keys if key.startswith("inchikey:"))
        )
        if exact_identity:
            choice = best_candidate(exact_identity)
            method = "exact_inchikey"
        if choice is None:
            name_candidates = set().union(
                *(lookup.get(key, set()) for key in keys if key.startswith("name:"))
            )
            if name_candidates:
                choice = best_candidate(name_candidates)
                method = "exact_name"
        if choice is None:
            candidates = set().union(
                *(lookup.get(key, set()) for key in keys if key.startswith("inchikey14:"))
            )
            if len(candidates) == 1:
                choice = next(iter(candidates))
                method = "inchikey14"
            elif len(candidates) > 1:
                method = "ambiguous"
        if choice is not None:
            selected.append(choice)
            methods.append(method)
        else:
            selected.append(None)
            methods.append(method)
    context_columns = [column for column in representative.columns if column != "drug_id"]
    joined = representative.reindex(selected)[context_columns].reset_index(drop=True)
    out = pd.concat([out, joined.add_prefix("pk_context_")], axis=1)
    out["pk_context_join_status"] = methods
    return out


def join_pk_context_summary(
    model_table: pd.DataFrame,
    context: pd.DataFrame,
) -> pd.DataFrame:
    out = model_table.copy()
    lookup: dict[str, set[int]] = {}
    for idx, row in context.iterrows():
        for key in _mapping_keys(row):
            lookup.setdefault(key, set()).add(idx)

    def joined_values(indices: set[int], column: str) -> str:
        values = {
            str(value).strip()
            for value in context.loc[list(indices), column].dropna()
            if str(value).strip()
        }
        return ";".join(sorted(values))

    summaries: list[dict[str, Any]] = []
    for _, row in out.iterrows():
        keys = _mapping_keys(row)
        indices = set().union(
            *(lookup.get(key, set()) for key in keys if key.startswith("inchikey:"))
        )
        if not indices:
            indices = set().union(
                *(lookup.get(key, set()) for key in keys if key.startswith("name:"))
            )
        if not indices:
            prefix = set().union(
                *(lookup.get(key, set()) for key in keys if key.startswith("inchikey14:"))
            )
            if len({str(context.loc[idx, "inchikey"]) for idx in prefix}) == 1:
                indices = prefix
        summaries.append(
            {
                "pk_context_scenario_count": len(indices),
                "pk_context_all_sources": joined_values(indices, "source_name") if indices else "",
                "pk_context_all_routes": joined_values(indices, "route") if indices else "",
                "pk_context_all_formulations": joined_values(indices, "formulation") if indices else "",
                "pk_context_all_regimens": joined_values(indices, "regimen") if indices else "",
            }
        )
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(summaries)], axis=1)


def write_pk_context_outputs(
    *,
    context: pd.DataFrame,
    model_table: pd.DataFrame,
    out_dir: str | Path,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    representative = select_representative_pk_context(context)
    enriched = join_representative_pk_context(model_table, representative)
    enriched = join_pk_context_summary(enriched, context)
    context.to_csv(out / "pk_context_long.csv", index=False)
    representative.to_csv(out / "pk_context_representative.csv", index=False)
    enriched.to_csv(out / "AtlasSPD_phase1_pk_enriched.csv", index=False)
    summary = context.groupby("source_name", dropna=False).agg(
        rows=("pk_context_id", "size"),
        drugs=("drug_name", "nunique"),
        cmax=("cmax_um", _numeric_coverage),
        free_cmax=("free_cmax_um", _numeric_coverage),
        dose=("dose_value", _numeric_coverage),
        route=("route", _text_coverage),
        formulation=("formulation", _text_coverage),
        protein_binding=("protein_binding_percent", _numeric_coverage),
        fraction_unbound=("fraction_unbound_plasma", _numeric_coverage),
        clearance=("clearance_value", _numeric_coverage),
        bioavailability=("bioavailability_value", _numeric_coverage),
    ).reset_index()
    summary.to_csv(out / "pk_context_source_coverage.csv", index=False)
    numeric_context = enriched[
        [
            "pk_context_cmax_um",
            "pk_context_free_cmax_um",
            "pk_context_dose_value",
        ]
    ].apply(pd.to_numeric, errors="coerce")
    text_context = pd.DataFrame(
        {
            "route": _text(enriched, "pk_context_route").ne(""),
            "formulation": _text(enriched, "pk_context_formulation").ne(""),
            "regimen": _text(enriched, "pk_context_regimen").ne(""),
        },
        index=enriched.index,
    )
    usable_context = numeric_context.notna().any(axis=1) | text_context.any(axis=1)
    external_context = ~_text(enriched, "pk_context_source_name").isin(
        {"", "SPD", "existing_phase1"}
    )
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy": "SPD labels unchanged; external PK is stored in separate pk_context_* columns",
        "context_rows": int(len(context)),
        "representative_rows": int(len(representative)),
        "model_rows": int(len(enriched)),
        "identity_matched_rows": int(enriched["pk_context_join_status"].isin(
            ["exact_inchikey", "exact_name", "inchikey14"]
        ).sum()),
        "exact_inchikey_rows": int(
            enriched["pk_context_join_status"].eq("exact_inchikey").sum()
        ),
        "exact_name_rows": int(
            enriched["pk_context_join_status"].eq("exact_name").sum()
        ),
        "inchikey14_rows": int(
            enriched["pk_context_join_status"].eq("inchikey14").sum()
        ),
        "usable_context_rows": int(usable_context.sum()),
        "external_representative_rows": int((usable_context & external_context).sum()),
        "ambiguous_rows": int(enriched["pk_context_join_status"].eq("ambiguous").sum()),
        "sources": summary.to_dict("records"),
    }
    (out / "pk_context_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest
