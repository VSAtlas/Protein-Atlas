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


PK_MECHANISTIC_CONTEXT_COLUMNS = [
    "pka_value",
    "pka_type",
    "pka_context",
    "logd_7_4",
    "logd_7_4_context",
    "aqueous_solubility_value",
    "aqueous_solubility_unit",
    "aqueous_solubility_context",
    "biorelevant_solubility_value",
    "biorelevant_solubility_unit",
    "biorelevant_solubility_medium",
    "caco2_papp_value",
    "caco2_papp_unit",
    "caco2_papp_direction",
    "peff_value",
    "peff_unit",
    "peff_context",
    "absorption_rate_ka_value",
    "absorption_rate_ka_unit",
    "absorption_rate_ka_context",
    "tmax_value",
    "tmax_unit",
    "tmax_context",
    "blood_to_plasma_ratio",
    "blood_to_plasma_context",
    "fraction_unbound_plasma_context",
    "intrinsic_clearance_value",
    "intrinsic_clearance_unit",
    "intrinsic_clearance_context",
    "hepatic_clearance_value",
    "hepatic_clearance_unit",
    "hepatic_clearance_context",
    "renal_clearance_value",
    "renal_clearance_unit",
    "renal_clearance_context",
    "systemic_clearance_value",
    "systemic_clearance_unit",
    "systemic_clearance_context",
    "vss_value",
    "vss_unit",
    "vss_context",
    "vd_value",
    "vd_unit",
    "vd_context",
    "half_life_value",
    "half_life_unit",
    "half_life_context",
    "body_weight_value",
    "body_weight_unit",
    "body_weight_context",
]

PK_CONTEXT_COLUMNS = [
    "drug_id",
    "drug_name",
    "inchikey",
    "source_name",
    "source_version",
    "source_record_id",
    "source_url",
    "study_id",
    "reference",
    "population",
    "species",
    "dose_value",
    "dose_unit",
    "dose_text",
    "route",
    "regimen",
    "formulation",
    "steady_state",
    "dose_context_type",
    "parent_or_metabolite",
    "measurement_context",
    "reference_route",
    "cmax_value_raw",
    "cmax_unit_raw",
    "cmax_um",
    "protein_binding_percent",
    "fraction_unbound_plasma",
    "free_cmax_um",
    "free_cmax_method",
    "clearance_value",
    "clearance_unit",
    "bioavailability_value",
    "bioavailability_unit",
    "extraction_method",
    "source_confidence",
    "context_status",
    "training_allowed",
    "license_note",
    "missing_reason",
    *PK_MECHANISTIC_CONTEXT_COLUMNS,
]

ALIASES = {
    "drug_id": [
        "drug_id",
        "atlas_drug_id",
        "compound_id",
        "drugbank_id",
        "substance_id",
        "id",
    ],
    "drug_name": [
        "drug_name",
        "display_name",
        "generic_name",
        "compound_name",
        "substance",
        "name",
    ],
    "inchikey": ["inchikey", "inchi_key", "drugcentral_inchikey"],
    "source_record_id": ["source_record_id", "record_id", "row_id", "pk_id", "id"],
    "study_id": ["study_id", "study", "study_sid", "trial_id"],
    "reference": [
        "reference",
        "citation",
        "publication",
        "source_reference",
        "pmid",
        "doi",
    ],
    "population": ["population", "subject_population", "group", "cohort", "age_group"],
    "species": ["species", "organism"],
    "dose_value": ["dose_value", "normalized_dose_value", "dose_amount_value"],
    "dose_unit": ["dose_unit", "normalized_dose_unit"],
    "dose_text": ["dose_text", "dosage", "dose_description"],
    "route": ["route", "route_of_administration", "dosage_route", "application"],
    "regimen": ["regimen", "dosing_regimen", "schedule", "frequency"],
    "formulation": ["formulation", "dosage_form", "form"],
    "steady_state": ["steady_state", "is_steady_state"],
    "dose_context_type": ["dose_context_type", "dose_provenance_type"],
    "measurement_context": ["measurement_context", "endpoint_context", "endpoint_type"],
    "reference_route": ["reference_route", "comparator_route"],
    "parent_or_metabolite": [
        "parent_or_metabolite",
        "analyte_type",
        "cmax_from_parent",
    ],
    "cmax_value_raw": ["cmax_value_raw", "cmax", "cmax_value"],
    "cmax_unit_raw": ["cmax_unit_raw", "cmax_unit"],
    "cmax_um": ["cmax_um", "cmax_uM", "total_cmax_um", "combined_cmax_um"],
    "protein_binding_percent": [
        "protein_binding_percent",
        "protein_binding",
        "ppb_percent",
        "ppb",
        "bound_percent",
    ],
    "fraction_unbound_plasma": [
        "fraction_unbound_plasma",
        "fu_plasma",
        "fup",
    ],
    "free_cmax_um": [
        "free_cmax_um",
        "free_cmax_uM",
        "cmax_free_um",
        "combined_free_cmax_um",
    ],
    "clearance_value": ["clearance_value", "pk_clearance_value"],
    "clearance_unit": ["clearance_unit", "cl_unit"],
    "bioavailability_value": [
        "bioavailability_value",
        "bioavailability",
        "absolute_bioavailability",
    ],
    "bioavailability_unit": ["bioavailability_unit"],
    "training_allowed": ["training_allowed", "license_allows_ml_training"],
    "pka_value": ["pka_value", "pka"],
    "pka_type": ["pka_type", "pka_class", "ionization_type"],
    "pka_context": ["pka_context", "ionization_context"],
    "logd_7_4": ["logd_7_4", "logd7_4", "logd_ph_7_4", "logd_ph7_4"],
    "logd_7_4_context": ["logd_7_4_context", "logd_context"],
    "aqueous_solubility_value": [
        "aqueous_solubility_value",
        "water_solubility_value",
    ],
    "aqueous_solubility_unit": [
        "aqueous_solubility_unit",
        "water_solubility_unit",
    ],
    "aqueous_solubility_context": [
        "aqueous_solubility_context",
        "water_solubility_context",
    ],
    "biorelevant_solubility_value": ["biorelevant_solubility_value"],
    "biorelevant_solubility_unit": ["biorelevant_solubility_unit"],
    "biorelevant_solubility_medium": [
        "biorelevant_solubility_medium",
        "biorelevant_medium",
    ],
    "caco2_papp_value": [
        "caco2_papp_value",
        "caco_2_papp_value",
        "caco2_apparent_permeability_value",
    ],
    "caco2_papp_unit": [
        "caco2_papp_unit",
        "caco_2_papp_unit",
        "caco2_apparent_permeability_unit",
    ],
    "caco2_papp_direction": [
        "caco2_papp_direction",
        "caco_2_papp_direction",
        "caco2_direction",
    ],
    "peff_value": ["peff_value", "effective_permeability_value"],
    "peff_unit": ["peff_unit", "effective_permeability_unit"],
    "peff_context": ["peff_context", "effective_permeability_context"],
    "absorption_rate_ka_value": [
        "absorption_rate_ka_value",
        "absorption_rate_constant_value",
        "ka_absorption_value",
    ],
    "absorption_rate_ka_unit": [
        "absorption_rate_ka_unit",
        "absorption_rate_constant_unit",
        "ka_absorption_unit",
    ],
    "absorption_rate_ka_context": ["absorption_rate_ka_context"],
    "tmax_value": ["tmax_value", "t_max_value", "tmax"],
    "tmax_unit": ["tmax_unit", "t_max_unit"],
    "tmax_context": ["tmax_context", "t_max_context"],
    "blood_to_plasma_ratio": [
        "blood_to_plasma_ratio",
        "blood_plasma_ratio",
        "b_p_ratio",
    ],
    "blood_to_plasma_context": [
        "blood_to_plasma_context",
        "blood_plasma_ratio_context",
    ],
    "fraction_unbound_plasma_context": [
        "fraction_unbound_plasma_context",
        "fu_plasma_context",
    ],
    "intrinsic_clearance_value": [
        "intrinsic_clearance_value",
        "clint_value",
        "cl_int_value",
    ],
    "intrinsic_clearance_unit": [
        "intrinsic_clearance_unit",
        "clint_unit",
        "cl_int_unit",
    ],
    "intrinsic_clearance_context": ["intrinsic_clearance_context"],
    "hepatic_clearance_value": ["hepatic_clearance_value", "clh_value"],
    "hepatic_clearance_unit": ["hepatic_clearance_unit", "clh_unit"],
    "hepatic_clearance_context": ["hepatic_clearance_context"],
    "renal_clearance_value": ["renal_clearance_value", "clr_value"],
    "renal_clearance_unit": ["renal_clearance_unit", "clr_unit"],
    "renal_clearance_context": ["renal_clearance_context"],
    "systemic_clearance_value": ["systemic_clearance_value"],
    "systemic_clearance_unit": ["systemic_clearance_unit"],
    "systemic_clearance_context": ["systemic_clearance_context"],
    "vss_value": ["vss_value", "steady_state_volume_of_distribution_value"],
    "vss_unit": ["vss_unit", "steady_state_volume_of_distribution_unit"],
    "vss_context": ["vss_context", "steady_state_volume_context"],
    "vd_value": ["vd_value", "volume_of_distribution_value"],
    "vd_unit": ["vd_unit", "volume_of_distribution_unit"],
    "vd_context": ["vd_context", "volume_of_distribution_context"],
    "half_life_value": [
        "half_life_value",
        "terminal_half_life_value",
        "elimination_half_life_value",
    ],
    "half_life_unit": [
        "half_life_unit",
        "terminal_half_life_unit",
        "elimination_half_life_unit",
    ],
    "half_life_context": ["half_life_context", "half_life_type"],
    "body_weight_value": ["body_weight_value", "subject_body_weight_value"],
    "body_weight_unit": ["body_weight_unit", "subject_body_weight_unit"],
    "body_weight_context": ["body_weight_context", "subject_body_weight_context"],
}

SOURCE_PRIORITY = {
    "SPD": 0,
    "existing_phase1": 1,
    "NCATS_Inxight_FRDB": 2,
    "PK-DB": 3,
    "PK-DB_recovered_source_TSV": 3,
    "DrugBank_Cmax": 4,
    "DrugBank_protein_binding": 4,
    "DailyMed_openFDA_SPL": 5,
    "DailyMed_openFDA_SPL_semantic_review": 5,
    "existing_pk": 6,
}

CLEARANCE_MEASUREMENT_CONTEXTS = frozenset(
    {
        "apparent_oral_clearance",
        "apparent_plasma_clearance",
        "blood_clearance",
        "hepatic_clearance",
        "intrinsic_clearance",
        "plasma_clearance",
        "renal_clearance",
        "systemic_clearance",
        "total_body_clearance",
    }
)

NUMERIC_COLUMNS = [
    "dose_value",
    "cmax_value_raw",
    "cmax_um",
    "protein_binding_percent",
    "fraction_unbound_plasma",
    "free_cmax_um",
    "clearance_value",
    "bioavailability_value",
    "pka_value",
    "logd_7_4",
    "aqueous_solubility_value",
    "biorelevant_solubility_value",
    "caco2_papp_value",
    "peff_value",
    "absorption_rate_ka_value",
    "tmax_value",
    "blood_to_plasma_ratio",
    "intrinsic_clearance_value",
    "hepatic_clearance_value",
    "renal_clearance_value",
    "systemic_clearance_value",
    "vss_value",
    "vd_value",
    "half_life_value",
    "body_weight_value",
]

PK_MEASUREMENT_VALUE_COLUMNS = [
    "cmax_value_raw",
    "cmax_um",
    "protein_binding_percent",
    "fraction_unbound_plasma",
    "free_cmax_um",
    "clearance_value",
    "bioavailability_value",
    "pka_value",
    "logd_7_4",
    "aqueous_solubility_value",
    "biorelevant_solubility_value",
    "caco2_papp_value",
    "peff_value",
    "absorption_rate_ka_value",
    "tmax_value",
    "blood_to_plasma_ratio",
    "intrinsic_clearance_value",
    "hepatic_clearance_value",
    "renal_clearance_value",
    "systemic_clearance_value",
    "vss_value",
    "vd_value",
    "half_life_value",
]

_CONTEXT_HASH_COLUMNS = tuple(
    column
    for column in PK_CONTEXT_COLUMNS
    if column not in {"context_status", "license_note", "missing_reason"}
)


def normalize_key(value: Any) -> str:
    if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
        return ""
    text = str(value).strip().casefold()
    if text in {"", "nan", "none", "null", "unknown"}:
        return ""
    return " ".join(re.sub(r"[^a-z0-9]+", " ", text).split())


def _text(frame: pd.DataFrame, column: str) -> pd.Series:
    return (
        frame.get(column, pd.Series("", index=frame.index))
        .fillna("")
        .astype(str)
        .str.strip()
    )


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(
        frame.get(column, pd.Series(pd.NA, index=frame.index)), errors="coerce"
    )


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
    return pd.DataFrame(
        {
            column: pd.Series(pd.NA, index=index, dtype="object")
            for column in PK_CONTEXT_COLUMNS
        }
    )


def _context_hash(row: pd.Series) -> str:
    fields = [row.get(column) for column in _CONTEXT_HASH_COLUMNS]
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
    valid_id = _text(out, "drug_id").str.len().gt(0) | _text(
        out, "drug_name"
    ).str.len().gt(0)
    out = out.loc[valid_id].copy()
    missing_measurement = out[PK_MEASUREMENT_VALUE_COLUMNS].isna().all(axis=1)
    out["missing_reason"] = out["missing_reason"].astype("object")
    out.loc[missing_measurement & out["missing_reason"].isna(), "missing_reason"] = (
        "no_recognized_pk_measurement"
    )
    out.insert(0, "pk_context_id", out.apply(_context_hash, axis=1))
    return out.drop_duplicates("pk_context_id", keep="first").reset_index(drop=True)


def load_spd_pk_context(
    path: str | Path,
    *,
    administration_context: str | Path | pd.DataFrame | None = None,
) -> pd.DataFrame:
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
    out["license_note"] = (
        "CC-BY-4.0 article supplement; upstream source restrictions may apply"
    )
    context = finalize_context(out)
    if administration_context is None:
        return context
    from analysis.external.spd_cmax_context import (
        merge_verified_spd_administration_context,
    )

    merged = merge_verified_spd_administration_context(
        context, administration_context
    ).drop(columns="pk_context_id", errors="ignore")
    return finalize_context(merged)


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
    out["drug_id"] = compound_id.map(lambda value: f"frdb:{value}" if value else "")
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
    out["regimen"] = _joined_text(_text(raw, "pk_experiment_type"), frequency)
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
    cmax_dose = out["cmax_um"].notna() & out["dose_value"].notna()
    out.loc[cmax_dose, "dose_context_type"] = "cmax_study_matched"
    out.loc[cmax_dose, "measurement_context"] = "observed_cmax"
    unsupported_cmax = out["cmax_value_raw"].notna() & out["cmax_um"].isna()
    fraction_unbound_percent = _numeric(raw, "pk_funbound_value")
    valid_fraction = fraction_unbound_percent.between(0.0, 100.0, inclusive="both")
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
    out.loc[unsupported_cmax, "context_status"] += "; cmax_unit_not_safely_convertible"
    out["license_note"] = "public NCATS FRDB release; retain source citation"
    return finalize_context(out)


def load_reviewed_openfda_pk_context(path: str | Path) -> pd.DataFrame:
    """Load only source-text-adjudicated DailyMed/openFDA PK scenarios."""

    raw = read_source_table(path)
    accepted = raw.get(
        "acceptable_for_model_training", pd.Series(False, index=raw.index)
    )
    accepted = (
        accepted.fillna(False)
        .astype(str)
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
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
    out["dose_value"] = _numeric(raw, "adjudicated_dose_value")
    out["dose_unit"] = _text(raw, "adjudicated_dose_unit")
    out["route"] = _text(raw, "adjudicated_route").where(
        _text(raw, "adjudicated_route").ne(""), _text(raw, "route")
    )
    out["regimen"] = _text(raw, "adjudicated_regimen").where(
        _text(raw, "adjudicated_regimen").ne(""), _text(raw, "regimen")
    )
    out["steady_state"] = (
        _text(raw, "steady_state").str.casefold().isin({"1", "true", "yes", "y"})
    )
    out["measurement_context"] = "observed_cmax"
    out["parent_or_metabolite"] = _text(raw, "adjudicated_context")
    out["cmax_value_raw"] = _first_numeric(_text(raw, "cmax_values_raw"))
    out["cmax_unit_raw"] = _text(raw, "cmax_units_raw")
    out["cmax_um"] = _first_numeric(_text(raw, "cmax_converted_um_candidates"))
    matched_dose = out["dose_value"].notna() & out["cmax_um"].notna()
    out.loc[matched_dose, "dose_context_type"] = "cmax_study_matched"
    out["protein_binding_percent"] = _first_numeric(
        _text(raw, "protein_binding_values_pct")
    )
    out["extraction_method"] = "manual_source_text_adjudication"
    out["source_confidence"] = "high"
    out["training_allowed"] = True
    out["context_status"] = (
        "source_text_adjudicated; contextual_scenario_not_universal_drug_pk"
    )
    out["license_note"] = "public FDA labeling; retain SPL set/version provenance"
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
    explicit_match = out["dose_context_type"].eq("cmax_study_matched")
    out.loc[explicit_match & out["cmax_um"].notna(), "measurement_context"] = (
        "observed_cmax"
    )
    measurement_context = _text(out, "measurement_context").str.casefold()
    clearance_candidate = pd.to_numeric(out["clearance_value"], errors="coerce").notna()
    valid_clearance = measurement_context.isin(CLEARANCE_MEASUREMENT_CONTEXTS) & _text(
        out, "clearance_unit"
    ).ne("")
    invalid_clearance = clearance_candidate & ~valid_clearance
    out.loc[
        invalid_clearance,
        ["clearance_value", "clearance_unit"],
    ] = pd.NA
    bioavailability_candidate = pd.to_numeric(
        out["bioavailability_value"], errors="coerce"
    ).notna()
    absolute_bioavailability = (
        measurement_context.eq("absolute_bioavailability")
        & _text(out, "reference_route").str.casefold().eq("intravenous")
        & _text(out, "route").ne("")
        & ~_text(out, "route").str.casefold().eq("intravenous")
    )
    invalid_bioavailability = bioavailability_candidate & ~absolute_bioavailability
    out.loc[
        invalid_bioavailability,
        ["bioavailability_value", "bioavailability_unit"],
    ] = pd.NA
    out["source_confidence"] = "medium"
    out["context_status"] = "source_context_preserved_where_available"
    out.loc[invalid_clearance, "context_status"] += (
        "; generic_clearance_quarantined_without_endpoint_context_and_unit"
    )
    out.loc[invalid_bioavailability, "context_status"] += (
        "; generic_bioavailability_quarantined_without_absolute_iv_reference"
    )
    return finalize_context(out)


def combine_pk_context(parts: Iterable[pd.DataFrame]) -> pd.DataFrame:
    usable = [part for part in parts if part is not None and not part.empty]
    if not usable:
        return pd.DataFrame(columns=["pk_context_id", *PK_CONTEXT_COLUMNS])
    raw = pd.concat(usable, ignore_index=True).drop(
        columns=["pk_context_id"], errors="ignore"
    )
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
        "dose_value",
        "dose_unit",
    ]
    ranked.loc[unreviewed_spl, numeric_columns] = pd.NA
    ranked["context_status"] = ranked["context_status"].astype("object")
    training_denied = (
        ranked["training_allowed"]
        .astype("string")
        .fillna("")
        .str.casefold()
        .isin({"0", "false", "no", "n"})
    )
    ranked.loc[training_denied, numeric_columns] = pd.NA
    unmatched_dose = ranked["dose_value"].notna() & ~ranked["dose_context_type"].eq(
        "cmax_study_matched"
    )
    ranked.loc[unmatched_dose, ["dose_value", "dose_unit"]] = pd.NA
    ranked.loc[unreviewed_spl, "context_status"] = (
        ranked.loc[unreviewed_spl, "context_status"].fillna("").astype(str)
        + "; numeric_quarantined_pending_source_text_review"
    ).str.lstrip("; ")
    ranked.loc[training_denied, "context_status"] = (
        ranked.loc[training_denied, "context_status"].fillna("").astype(str)
        + "; numeric_quarantined_training_not_allowed"
    ).str.lstrip("; ")
    ranked["_source_rank"] = ranked["source_name"].map(SOURCE_PRIORITY).fillna(99)
    ranked["_measurement_rank"] = 3
    ranked.loc[
        pd.to_numeric(ranked["dose_value"], errors="coerce").notna(),
        "_measurement_rank",
    ] = 2
    ranked.loc[ranked["cmax_um"].notna(), "_measurement_rank"] = 1
    ranked.loc[ranked["free_cmax_um"].notna(), "_measurement_rank"] = 0
    ranked["_confidence_rank"] = (
        ranked["source_confidence"].map({"high": 0, "medium": 1, "low": 2}).fillna(3)
    )
    context_present = pd.DataFrame(
        {
            "dose": pd.to_numeric(ranked["dose_value"], errors="coerce").notna(),
            "route": _text(ranked, "route").ne(""),
            "formulation": _text(ranked, "formulation").ne(""),
            "regimen": _text(ranked, "regimen").ne(""),
        },
        index=ranked.index,
    )
    ranked["_context_rank"] = context_present.sum(axis=1) * -1
    ranked["_drug_key"] = _text(ranked, "inchikey").str.upper()
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(ranked.loc[missing], "drug_id").map(
        normalize_key
    )
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(
        ranked.loc[missing], "drug_name"
    ).map(normalize_key)
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
    return out.drop(
        columns=[column for column in out if column.startswith("_")]
    ).reset_index(drop=True)


def select_endpoint_pk_context(
    context: pd.DataFrame,
    *,
    value_column: str,
    unit_column: str,
    required_contexts: set[str] | None = None,
) -> pd.DataFrame:
    """Select context-specific endpoint rows without pooling endpoint families."""

    if context.empty or value_column not in context:
        return context.iloc[0:0].copy()
    ranked = context.loc[
        pd.to_numeric(context[value_column], errors="coerce").notna()
    ].copy()
    if required_contexts is not None:
        ranked = ranked.loc[
            _text(ranked, "measurement_context").isin(required_contexts)
        ].copy()
    if ranked.empty:
        return ranked
    allowed = (
        ranked["training_allowed"]
        .astype("string")
        .fillna("")
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
    )
    ranked = ranked.loc[allowed].copy()
    if ranked.empty:
        return ranked
    ranked["_source_rank"] = ranked["source_name"].map(SOURCE_PRIORITY).fillna(99)
    ranked["_confidence_rank"] = (
        ranked["source_confidence"].map({"high": 0, "medium": 1, "low": 2}).fillna(3)
    )
    ranked["_drug_key"] = _text(ranked, "inchikey").str.upper()
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(ranked.loc[missing], "drug_id").map(
        normalize_key
    )
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(
        ranked.loc[missing], "drug_name"
    ).map(normalize_key)
    ranked["_endpoint_group"] = (
        ranked["_drug_key"] + "|" + _text(ranked, "measurement_context").str.casefold()
    )
    ranked["_endpoint_signature"] = (
        pd.to_numeric(ranked[value_column], errors="coerce").map(
            lambda value: f"{value:.12g}"
        )
        + "|"
        + _text(ranked, unit_column).str.casefold()
        + "|"
        + _text(ranked, "route").str.casefold()
        + "|"
        + _text(ranked, "reference_route").str.casefold()
        + "|"
        + _text(ranked, "population").str.casefold()
    )
    signature_counts = ranked.groupby("_endpoint_group")[
        "_endpoint_signature"
    ].transform("nunique")
    ranked["endpoint_context_count"] = signature_counts
    ranked["endpoint_selection_status"] = "unique_context"
    ranked.loc[signature_counts.gt(1), "endpoint_selection_status"] = (
        "ambiguous_multiple_contexts_excluded"
    )
    ranked = ranked.loc[signature_counts.eq(1)].sort_values(
        ["_source_rank", "_confidence_rank", "pk_context_id"],
        kind="stable",
    )
    selected = ranked.drop_duplicates("_endpoint_group", keep="first").copy()
    return selected.drop(
        columns=[column for column in selected if column.startswith("_")]
    ).reset_index(drop=True)


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
    for column in (
        "drug_id",
        "drug_name",
        "display_name",
        "generic_name",
        "ligand_base",
    ):
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
                *(
                    lookup.get(key, set())
                    for key in keys
                    if key.startswith("inchikey14:")
                )
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
    context_columns = [
        column for column in representative.columns if column != "drug_id"
    ]
    joined = representative.reindex(selected)[context_columns].reset_index(drop=True)
    out = pd.concat([out, joined.add_prefix("pk_context_")], axis=1)
    out["pk_context_join_status"] = methods
    return out


def join_endpoint_pk_context(
    model_table: pd.DataFrame,
    representative: pd.DataFrame,
    *,
    prefix: str,
) -> pd.DataFrame:
    """Join a separately selected PK endpoint while retaining its provenance."""

    out = model_table.reset_index(drop=True).copy()
    stale = [column for column in out if column.startswith(prefix)]
    if stale:
        out = out.drop(columns=stale)
    if representative.empty:
        out[f"{prefix}join_status"] = "no_context_sources_available"
        return out
    lookup: dict[str, set[int]] = {}
    for idx, row in representative.iterrows():
        for key in _mapping_keys(row):
            lookup.setdefault(key, set()).add(idx)

    selected: list[int | None] = []
    methods: list[str] = []
    for _, row in out.iterrows():
        keys = _mapping_keys(row)
        exact_identity = set().union(
            *(lookup.get(key, set()) for key in keys if key.startswith("inchikey:"))
        )
        name_matches = set().union(
            *(lookup.get(key, set()) for key in keys if key.startswith("name:"))
        )
        prefix_matches = set().union(
            *(lookup.get(key, set()) for key in keys if key.startswith("inchikey14:"))
        )
        if len(exact_identity) == 1:
            selected.append(next(iter(exact_identity)))
            methods.append("exact_inchikey")
        elif len(name_matches) == 1:
            selected.append(next(iter(name_matches)))
            methods.append("exact_name")
        elif len(prefix_matches) == 1:
            selected.append(next(iter(prefix_matches)))
            methods.append("inchikey14")
        else:
            selected.append(None)
            methods.append(
                "ambiguous"
                if exact_identity or name_matches or prefix_matches
                else "unmatched"
            )
    context_columns = [
        column for column in representative.columns if column != "drug_id"
    ]
    joined = representative.reindex(selected)[context_columns].reset_index(drop=True)
    out = pd.concat([out, joined.add_prefix(prefix)], axis=1)
    out[f"{prefix}join_status"] = methods
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
                *(
                    lookup.get(key, set())
                    for key in keys
                    if key.startswith("inchikey14:")
                )
            )
            if len({str(context.loc[idx, "inchikey"]) for idx in prefix}) == 1:
                indices = prefix
        summaries.append(
            {
                "pk_context_scenario_count": len(indices),
                "pk_context_all_sources": joined_values(indices, "source_name")
                if indices
                else "",
                "pk_context_all_routes": joined_values(indices, "route")
                if indices
                else "",
                "pk_context_all_formulations": joined_values(indices, "formulation")
                if indices
                else "",
                "pk_context_all_regimens": joined_values(indices, "regimen")
                if indices
                else "",
            }
        )
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(summaries)], axis=1)


def validate_pk_context_selection(
    *,
    context: pd.DataFrame,
    representative: pd.DataFrame,
    clearance: pd.DataFrame,
    bioavailability: pd.DataFrame,
    enriched: pd.DataFrame,
) -> dict[str, Any]:
    """Validate that contextual PK evidence cannot leak into universal fields."""

    dose = pd.to_numeric(representative["dose_value"], errors="coerce").notna()
    cmax = pd.to_numeric(representative["cmax_um"], errors="coerce").notna()
    low_spl = representative["source_name"].eq("DailyMed_openFDA_SPL") & representative[
        "source_confidence"
    ].eq("low")
    low_spl_numeric = (
        representative[["cmax_um", "free_cmax_um", "dose_value"]]
        .apply(pd.to_numeric, errors="coerce")
        .notna()
        .any(axis=1)
    )

    clearance_allowed = (
        clearance["training_allowed"]
        .astype("string")
        .fillna("")
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
    )
    bioavailability_allowed = (
        bioavailability["training_allowed"]
        .astype("string")
        .fillna("")
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
    )
    bioavailability_value = pd.to_numeric(
        bioavailability["bioavailability_value"], errors="coerce"
    )
    bioavailability_route = _text(bioavailability, "route").str.casefold()
    bioavailability_reference = _text(bioavailability, "reference_route").str.casefold()

    violations = {
        "duplicate_context_ids": int(context["pk_context_id"].duplicated().sum()),
        "duplicate_output_columns": int(enriched.columns.duplicated().sum()),
        "primary_dose_not_cmax_matched": int(
            (dose & ~representative["dose_context_type"].eq("cmax_study_matched")).sum()
        ),
        "primary_dose_without_cmax": int((dose & ~cmax).sum()),
        "maximum_dose_in_primary_context": int(
            _text(context, "measurement_context")
            .str.casefold()
            .str.startswith("maximum_")
            .sum()
        ),
        "unreviewed_spl_numeric_selected": int((low_spl & low_spl_numeric).sum()),
        "training_denied_numeric_selected": int(
            (
                representative["training_allowed"]
                .astype("string")
                .fillna("")
                .str.casefold()
                .isin({"0", "false", "no", "n"})
                & low_spl_numeric
            ).sum()
        ),
        "clearance_not_training_allowed": int((~clearance_allowed).sum()),
        "clearance_missing_unit_or_context": int(
            (
                _text(clearance, "clearance_unit").eq("")
                | _text(clearance, "measurement_context").eq("")
            ).sum()
        ),
        "clearance_context_not_partitioned": int(
            (
                ~_text(clearance, "measurement_context")
                .str.casefold()
                .isin(CLEARANCE_MEASUREMENT_CONTEXTS)
            ).sum()
        ),
        "bioavailability_not_training_allowed": int((~bioavailability_allowed).sum()),
        "bioavailability_invalid_value": int(
            (
                ~bioavailability_value.between(
                    0.0,
                    100.0,
                    inclusive="right",
                )
            ).sum()
        ),
        "bioavailability_not_absolute": int(
            (
                ~_text(bioavailability, "measurement_context").eq(
                    "absolute_bioavailability"
                )
            ).sum()
        ),
        "bioavailability_route_not_extravascular": int(
            (
                bioavailability_route.eq("") | bioavailability_route.eq("intravenous")
            ).sum()
        ),
        "bioavailability_reference_not_intravenous": int(
            (~bioavailability_reference.eq("intravenous")).sum()
        ),
    }
    return {
        "status": "pass" if not any(violations.values()) else "fail",
        "violations": violations,
        "policy": {
            "primary_dose": "dose matched to the selected Cmax context only",
            "maximum_recommended_dose": "separate sensitivity artifact only",
            "clearance": "context-specific endpoint with unit and provenance",
            "bioavailability": (
                "absolute extravascular value versus intravenous reference"
            ),
        },
    }


def write_pk_context_outputs(
    *,
    context: pd.DataFrame,
    model_table: pd.DataFrame,
    out_dir: str | Path,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    representative = select_representative_pk_context(context)
    clearance = select_endpoint_pk_context(
        context,
        value_column="clearance_value",
        unit_column="clearance_unit",
        required_contexts=set(CLEARANCE_MEASUREMENT_CONTEXTS),
    )
    bioavailability = select_endpoint_pk_context(
        context,
        value_column="bioavailability_value",
        unit_column="bioavailability_unit",
        required_contexts={"absolute_bioavailability"},
    )
    stale_endpoint_columns = [
        column
        for column in model_table
        if column.startswith(("pk_clearance_", "pk_bioavailability_"))
    ]
    clean_model_table = model_table.drop(
        columns=stale_endpoint_columns,
        errors="ignore",
    )
    enriched = join_representative_pk_context(clean_model_table, representative)
    for measurement_context in sorted(_text(clearance, "measurement_context").unique()):
        endpoint_rows = clearance.loc[
            _text(clearance, "measurement_context").eq(measurement_context)
        ].copy()
        suffix = measurement_context.removesuffix("_clearance")
        enriched = join_endpoint_pk_context(
            enriched,
            endpoint_rows,
            prefix=f"pk_clearance_{suffix}_",
        )
    enriched = join_endpoint_pk_context(
        enriched,
        bioavailability,
        prefix="pk_bioavailability_absolute_",
    )
    enriched = join_pk_context_summary(enriched, context)
    context.to_csv(out / "pk_context_long.csv", index=False)
    representative.to_csv(out / "pk_context_representative.csv", index=False)
    clearance.to_csv(out / "pk_clearance_context_representative.csv", index=False)
    bioavailability.to_csv(
        out / "pk_bioavailability_context_representative.csv", index=False
    )
    enriched.to_csv(out / "AtlasSPD_phase1_pk_enriched.csv", index=False)
    validation = validate_pk_context_selection(
        context=context,
        representative=representative,
        clearance=clearance,
        bioavailability=bioavailability,
        enriched=enriched,
    )
    (out / "pk_context_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n",
        encoding="utf-8",
    )
    if validation["status"] != "pass":
        raise ValueError(
            "PK context validation failed: "
            + json.dumps(validation["violations"], sort_keys=True)
        )
    summary = (
        context.groupby("source_name", dropna=False)
        .agg(
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
        )
        .reset_index()
    )
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
        "identity_matched_rows": int(
            enriched["pk_context_join_status"]
            .isin(["exact_inchikey", "exact_name", "inchikey14"])
            .sum()
        ),
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
        "clearance_model_rows": int(
            enriched[
                [
                    column
                    for column in enriched
                    if column.startswith("pk_clearance_")
                    and column.endswith("clearance_value")
                ]
            ]
            .apply(pd.to_numeric, errors="coerce")
            .notna()
            .any(axis=1)
            .sum()
        )
        if any(
            column.startswith("pk_clearance_") and column.endswith("clearance_value")
            for column in enriched
        )
        else 0,
        "bioavailability_model_rows": int(
            enriched[
                [
                    column
                    for column in enriched
                    if column.startswith("pk_bioavailability_")
                    and column.endswith("bioavailability_value")
                ]
            ]
            .apply(pd.to_numeric, errors="coerce")
            .notna()
            .any(axis=1)
            .sum()
        )
        if any(
            column.startswith("pk_bioavailability_")
            and column.endswith("bioavailability_value")
            for column in enriched
        )
        else 0,
        "clearance_contexts": int(len(clearance)),
        "bioavailability_contexts": int(len(bioavailability)),
        "validation": validation,
        "ambiguous_rows": int(enriched["pk_context_join_status"].eq("ambiguous").sum()),
        "sources": summary.to_dict("records"),
    }
    (out / "pk_context_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest
