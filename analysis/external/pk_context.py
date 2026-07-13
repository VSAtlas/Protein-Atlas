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
    text = str(value or "").strip().casefold()
    if text in {"", "nan", "none", "null", "unknown"}:
        return ""
    return " ".join(re.sub(r"[^a-z0-9]+", " ", text).split())


def _text(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame.get(column, pd.Series("", index=frame.index)).fillna("").astype(str).str.strip()


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame.get(column, pd.Series(pd.NA, index=frame.index)), errors="coerce")


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
    ranked["_free_rank"] = ranked["free_cmax_um"].isna().astype(int)
    ranked["_confidence_rank"] = ranked["source_confidence"].map(
        {"high": 0, "medium": 1, "low": 2}
    ).fillna(3)
    ranked["_context_rank"] = (
        ranked[["dose_value", "route", "formulation", "regimen"]].notna().sum(axis=1) * -1
    )
    ranked["_drug_key"] = _text(ranked, "drug_name").map(normalize_key)
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(ranked.loc[missing], "inchikey").str.upper()
    missing = ranked["_drug_key"].str.len().eq(0)
    ranked.loc[missing, "_drug_key"] = _text(ranked.loc[missing], "drug_id").map(normalize_key)
    ranked = ranked.sort_values(
        ["_free_rank", "_source_rank", "_confidence_rank", "_context_rank", "pk_context_id"],
        kind="stable",
    )
    out = ranked.drop_duplicates("_drug_key", keep="first").copy()
    return out.drop(columns=[column for column in out if column.startswith("_")]).reset_index(drop=True)


def _mapping_keys(row: pd.Series) -> set[str]:
    keys: set[str] = set()
    inchikey = str(row.get("inchikey") or "").strip().upper()
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
        candidates = set().union(
            *(
                lookup.get(key, set())
                for key in keys
                if key.startswith(("inchikey:", "name:"))
            )
        )
        if candidates:
            choice = best_candidate(candidates)
            method = "exact_identity_or_name"
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
            *(lookup.get(key, set()) for key in keys if key.startswith(("inchikey:", "name:")))
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
        cmax=("cmax_um", lambda values: int(values.notna().sum())),
        free_cmax=("free_cmax_um", lambda values: int(values.notna().sum())),
        dose=("dose_value", lambda values: int(values.notna().sum())),
        route=("route", lambda values: int(values.notna().sum())),
        formulation=("formulation", lambda values: int(values.notna().sum())),
    ).reset_index()
    summary.to_csv(out / "pk_context_source_coverage.csv", index=False)
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy": "SPD labels unchanged; external PK is stored in separate pk_context_* columns",
        "context_rows": int(len(context)),
        "representative_rows": int(len(representative)),
        "model_rows": int(len(enriched)),
        "matched_rows": int(enriched["pk_context_join_status"].isin(
            ["exact_identity_or_name", "inchikey14"]
        ).sum()),
        "ambiguous_rows": int(enriched["pk_context_join_status"].eq("ambiguous").sum()),
        "sources": summary.to_dict("records"),
    }
    (out / "pk_context_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest
