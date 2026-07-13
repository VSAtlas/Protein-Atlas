from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd

from analysis.external.openfda_pk import (
    _dose_context,
    _first_text,
    _openfda_values,
    _values,
)
from analysis.reporting.ligand_side_effect_cache import (
    _CMAX_VALUE_RE,
    _PROTEIN_BINDING_RE,
    _cmax_to_um,
    _flatten_label_text,
)


_SYSTEMIC_ROUTES = {
    "buccal",
    "intramuscular",
    "intravenous",
    "oral",
    "subcutaneous",
    "sublingual",
    "transdermal",
}
_NON_SYSTEMIC_ROUTES = {
    "dental",
    "intraocular",
    "ophthalmic",
    "otic",
    "topical",
}
_METABOLITE_RE = re.compile(
    r"\b(?:active\s+metabolite|metabolite|nor[a-z0-9-]+)\b",
    re.IGNORECASE,
)
_RANGE_OR_TABLE_RE = re.compile(
    r"\bfrom\s+[0-9.]+(?:\s*[a-zA-Zµ/]+)?\s+(?:to|through)\s+[0-9.]+|"
    r"\b[0-9.]+\s*(?:ng/mL|mcg/mL|ug/mL|mg/L|nM|uM)\s*"
    r"(?:-|to|through)\s*[0-9.]+|"
    r"\b[0-9.]+\s*(?:-|to)\s*[0-9.]+\s*(?:ng/mL|mcg/mL|ug/mL|mg/L|nM|uM)|"
    r"\btable\s+[0-9]+",
    re.IGNORECASE,
)


def _normalize_name(value: Any) -> str:
    return " ".join(re.findall(r"[a-z0-9]+", str(value or "").casefold()))


def _name_match_status(drug_id: str, names: list[str]) -> str:
    query = _normalize_name(drug_id)
    normalized = [_normalize_name(name) for name in names if _normalize_name(name)]
    if not query or not normalized:
        return "unknown"
    if query in normalized:
        return "exact"
    if any(query in name.split() or name.startswith(f"{query} ") for name in normalized):
        return "parent_or_salt_variant"
    return "mismatch"


def _route_status(routes: list[str]) -> str:
    tokens = {
        token
        for route in routes
        for token in re.findall(r"[a-z]+", str(route).casefold())
    }
    if tokens & _NON_SYSTEMIC_ROUTES:
        return "non_systemic"
    if tokens & _SYSTEMIC_ROUTES:
        return "systemic"
    return "unknown"


def _join(values: list[Any]) -> str:
    return ";".join(dict.fromkeys(str(value) for value in values if str(value)))


def _review_record(
    *,
    cache_file: Path,
    payload: dict[str, Any],
    record: dict[str, Any],
    molecular_weight: float | None,
) -> dict[str, Any]:
    drug_id = str(payload.get("drug_id") or "")
    generic_names = _openfda_values(record, "generic_name")
    substance_names = _openfda_values(record, "substance_name")
    routes = [*_values(record.get("route")), *_openfda_values(record, "route")]
    forms = [
        *_values(record.get("dosage_form")),
        *_openfda_values(record, "dosage_form"),
    ]
    pk_text = _first_text(
        record,
        ("pharmacokinetics", "clinical_pharmacology"),
        limit=50000,
    )
    dose_text = _first_text(
        record,
        ("dosage_and_administration", "dosage_forms_and_strengths"),
        limit=20000,
    )
    full_text = _flatten_label_text(record)
    cmax_matches = list(_CMAX_VALUE_RE.finditer(pk_text))
    binding_matches = list(_PROTEIN_BINDING_RE.finditer(full_text))
    cmax_values = [float(match.group(1)) for match in cmax_matches]
    cmax_units = [match.group(2) for match in cmax_matches]
    converted = [
        _cmax_to_um(value, unit, molecular_weight)
        for value, unit in zip(cmax_values, cmax_units)
    ]
    converted = [value for value in converted if value is not None]
    bound_values = [match.group(1).replace(" ", "") for match in binding_matches]
    dose_value, dose_unit, regimen, steady_state = _dose_context(
        f"{pk_text} {dose_text}"
    )
    cmax_context = ""
    if cmax_matches:
        first = cmax_matches[0]
        cmax_context = pk_text[max(0, first.start() - 300) : first.end() + 500]
    range_or_table = bool(_RANGE_OR_TABLE_RE.search(cmax_context))
    metabolite_ambiguity = bool(_METABOLITE_RE.search(cmax_context))
    name_status = _name_match_status(drug_id, [*generic_names, *substance_names])
    route_status = _route_status(routes)

    reasons: list[str] = []
    if not cmax_matches and not binding_matches:
        review_status = "no_numeric_candidate"
    else:
        if name_status in {"mismatch", "unknown"}:
            reasons.append(f"name_{name_status}")
        if route_status != "systemic":
            reasons.append(f"route_{route_status}")
        if len(cmax_matches) > 1:
            reasons.append("multiple_cmax_candidates")
        if range_or_table:
            reasons.append("range_or_table_context")
        if metabolite_ambiguity:
            reasons.append("parent_metabolite_ambiguity")
        if cmax_matches and not converted:
            reasons.append("cmax_unit_or_molecular_weight_not_convertible")
        if cmax_matches and dose_value is None:
            reasons.append("dose_context_missing")
        if name_status == "mismatch" or route_status == "non_systemic":
            review_status = "reject_mapping_or_route"
        elif reasons:
            review_status = "manual_review_required"
        else:
            review_status = "machine_review_candidate"

    set_id = str(record.get("set_id") or record.get("id") or "")
    return {
        "cache_file": str(cache_file),
        "drug_id": drug_id,
        "search": str(payload.get("search") or ""),
        "source_record_id": set_id,
        "source_url": (
            f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={set_id}"
            if set_id
            else ""
        ),
        "spl_version": record.get("version"),
        "effective_time": record.get("effective_time"),
        "label_generic_names": _join(generic_names),
        "label_substance_names": _join(substance_names),
        "route": _join(routes),
        "dosage_form": _join(forms),
        "name_match_status": name_status,
        "route_status": route_status,
        "molecular_weight": molecular_weight,
        "pk_text_present": bool(pk_text),
        "cmax_candidate_count": len(cmax_matches),
        "cmax_values_raw": _join(cmax_values),
        "cmax_units_raw": _join(cmax_units),
        "cmax_converted_um_candidates": _join(
            [f"{value:.9g}" for value in converted]
        ),
        "protein_binding_candidate_count": len(binding_matches),
        "protein_binding_values_pct": _join(bound_values),
        "dose_value": dose_value,
        "dose_unit": dose_unit,
        "regimen": regimen,
        "steady_state": steady_state,
        "range_or_table_context": range_or_table,
        "metabolite_ambiguity": metabolite_ambiguity,
        "review_status": review_status,
        "review_reasons": _join(reasons),
        "source_excerpt": " ".join(cmax_context.split())[:1600],
        "manual_decision": "",
        "reviewer": "",
        "review_date": "",
        "review_notes": "",
    }


def audit_openfda_pk_cache(
    *,
    cache_dir: str | Path,
    model_table: pd.DataFrame,
    out_dir: str | Path,
) -> dict[str, Any]:
    cache = Path(cache_dir)
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    mw_lookup = (
        model_table.assign(
            _mw=pd.to_numeric(model_table.get("rdkit_mol_wt"), errors="coerce")
        )
        .dropna(subset=["drug_id"])
        .drop_duplicates("drug_id")
        .set_index("drug_id")["_mw"]
        .to_dict()
    )
    rows: list[dict[str, Any]] = []
    parse_failures = 0
    cache_files = sorted(cache.glob("*.json")) if cache.exists() else []
    for cache_file in cache_files:
        try:
            payload = json.loads(cache_file.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            parse_failures += 1
            continue
        drug_id = str(payload.get("drug_id") or "")
        raw_mw = mw_lookup.get(drug_id)
        molecular_weight = (
            float(raw_mw) if pd.notna(raw_mw) and float(raw_mw) > 0 else None
        )
        for record in payload.get("records") or []:
            if isinstance(record, dict):
                rows.append(
                    _review_record(
                        cache_file=cache_file,
                        payload=payload,
                        record=record,
                        molecular_weight=molecular_weight,
                    )
                )
    review = pd.DataFrame(rows)
    review_path = output / "openfda_spl_source_text_review.csv"
    review.to_csv(review_path, index=False)
    if review.empty:
        queue = review.copy()
        summary = pd.DataFrame(columns=["review_status", "records", "drugs"])
    else:
        queue = review.loc[
            review["review_status"].isin(
                ["machine_review_candidate", "manual_review_required"]
            )
        ].copy()
        summary = (
            review.groupby("review_status", dropna=False)
            .agg(records=("source_record_id", "size"), drugs=("drug_id", "nunique"))
            .reset_index()
        )
    queue.to_csv(output / "openfda_spl_manual_review_queue.csv", index=False)
    summary.to_csv(output / "openfda_spl_source_text_summary.csv", index=False)
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cache_dir": str(cache),
        "cache_files": len(cache_files),
        "cache_parse_failures": parse_failures,
        "records_reviewed": int(len(review)),
        "drugs_reviewed": int(review["drug_id"].nunique()) if not review.empty else 0,
        "records_with_cmax_candidates": int(
            review["cmax_candidate_count"].gt(0).sum()
        )
        if not review.empty
        else 0,
        "records_with_binding_candidates": int(
            review["protein_binding_candidate_count"].gt(0).sum()
        )
        if not review.empty
        else 0,
        "machine_review_candidates": int(
            review["review_status"].eq("machine_review_candidate").sum()
        )
        if not review.empty
        else 0,
        "manual_review_required": int(
            review["review_status"].eq("manual_review_required").sum()
        )
        if not review.empty
        else 0,
        "rejected_mapping_or_route": int(
            review["review_status"].eq("reject_mapping_or_route").sum()
        )
        if not review.empty
        else 0,
        "policy": (
            "machine_review_candidate is a review queue state, not trusted PK; "
            "no SPL numeric value is promoted automatically"
        ),
        "review_output": str(review_path),
    }
    (output / "openfda_spl_source_text_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest
