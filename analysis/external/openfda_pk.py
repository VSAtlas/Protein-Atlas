from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import re
import time
from typing import Any

import pandas as pd
import requests  # type: ignore[import-untyped]

from analysis.external.pk_context import PK_CONTEXT_COLUMNS, finalize_context
from analysis.reporting.ligand_side_effect_cache import (
    LigandQuery,
    OPENFDA_LABEL_SOURCE,
    OPENFDA_LABEL_URL,
    _extract_label_exposure,
    _label_query_strings,
    _openfda_label_records,
)


LOG = logging.getLogger(__name__)

_DOSE_RE = re.compile(
    r"(?:following|after|received|administered|dose(?:d)?(?:\s+at|\s+of)?)[^.;\n]{0,100}?"
    r"([0-9]+(?:\.[0-9]+)?)\s*(mg/kg|mcg/kg|ug/kg|µg/kg|g|mg|mcg|ug|µg)\b",
    re.IGNORECASE,
)
_REGIMEN_RE = re.compile(
    r"\b(single dose|once daily|twice daily|three times daily|every\s+[0-9]+\s+hours|"
    r"[0-9]+\s+times daily|steady[ -]state)\b",
    re.IGNORECASE,
)


def _values(value: Any) -> list[str]:
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value or "").strip()
    return [text] if text and text.lower() not in {"nan", "none"} else []


def _first_text(record: dict[str, Any], fields: tuple[str, ...], limit: int = 2000) -> str:
    parts: list[str] = []
    for field in fields:
        parts.extend(_values(record.get(field)))
    return " ".join(" ".join(parts).split())[:limit]


def _openfda_values(record: dict[str, Any], field: str) -> list[str]:
    nested = record.get("openfda")
    if not isinstance(nested, dict):
        return []
    return _values(nested.get(field))


def _dose_context(text: str) -> tuple[float | None, str, str, str]:
    dose_match = _DOSE_RE.search(text)
    regimen_match = _REGIMEN_RE.search(text)
    dose_value = float(dose_match.group(1)) if dose_match else None
    dose_unit = dose_match.group(2).replace("µ", "u") if dose_match else ""
    regimen = regimen_match.group(1) if regimen_match else ""
    steady_state = (
        "yes"
        if "steady" in regimen.casefold()
        else "no" if regimen.casefold() == "single dose" else "unknown"
    )
    return dose_value, dose_unit, regimen, steady_state


def _query_from_row(row: pd.Series) -> LigandQuery:
    names: list[str] = []
    for column in ("generic_name", "display_name", "drug_id"):
        value = str(row.get(column) or "").strip()
        if value and value.lower() not in {"nan", "none"} and value.casefold() not in {
            item.casefold() for item in names
        }:
            names.append(value)
    mw = pd.to_numeric(pd.Series([row.get("rdkit_mol_wt")]), errors="coerce").iloc[0]
    return LigandQuery(
        ligand_label=str(row.get("drug_id") or ""),
        ligand_base=str(row.get("drug_id") or ""),
        ligand_display=str(row.get("display_name") or row.get("generic_name") or ""),
        query_names=tuple(names[:4]),
        brand_names=(),
        rxcuis=(),
        uniis=(),
        molecular_weight=float(mw) if pd.notna(mw) and float(mw) > 0 else None,
        mapping_used=True,
    )


def _record_context(
    *,
    drug_row: pd.Series,
    record: dict[str, Any],
    last_updated: str,
    search: str,
) -> dict[str, Any]:
    query = _query_from_row(drug_row)
    exposure = _extract_label_exposure([record], molecular_weight=query.molecular_weight)
    routes = [*_values(record.get("route")), *_openfda_values(record, "route")]
    forms = [
        *_values(record.get("dosage_form")),
        *_openfda_values(record, "dosage_form"),
        *_values(record.get("dosage_forms_and_strengths")),
    ]
    set_id = str(record.get("set_id") or record.get("id") or "")
    dose_text = _first_text(
        record,
        ("dosage_and_administration", "dosage_forms_and_strengths"),
    )
    pk_text = _first_text(record, ("pharmacokinetics", "clinical_pharmacology"))
    source_excerpt = pk_text[:1000]
    dose_value, dose_unit, regimen, steady_state = _dose_context(pk_text + " " + dose_text)
    row = {column: pd.NA for column in PK_CONTEXT_COLUMNS}
    row.update(
        {
            "drug_id": drug_row.get("drug_id"),
            "drug_name": drug_row.get("generic_name")
            or drug_row.get("display_name")
            or drug_row.get("drug_id"),
            "inchikey": drug_row.get("inchikey"),
            "source_name": "DailyMed_openFDA_SPL",
            "source_version": last_updated or str(record.get("effective_time") or ""),
            "source_record_id": set_id,
            "source_url": f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={set_id}" if set_id else OPENFDA_LABEL_URL,
            "reference": f"openFDA search={search}; SPL version={record.get('version', '')}; excerpt={source_excerpt}",
            "species": "Homo sapiens",
            "dose_value": dose_value,
            "dose_unit": dose_unit,
            "dose_text": dose_text,
            "route": ";".join(dict.fromkeys(routes)),
            "regimen": regimen,
            "formulation": ";".join(dict.fromkeys(forms)),
            "steady_state": steady_state,
            "cmax_um": exposure.get("cmax_um"),
            "fraction_unbound_plasma": exposure.get("fraction_unbound_plasma"),
            "free_cmax_um": exposure.get("free_cmax_um"),
            "extraction_method": "regex_from_structured_product_label",
            "source_confidence": "low",
            "context_status": "label_context_preserved; numeric extraction_requires_review",
            "license_note": "US government public labeling; preserve SPL set_id/version",
            "missing_reason": "" if exposure else "no_machine_extractable_cmax_or_fraction_unbound",
        }
    )
    return row


def _openfda_eligible_mask(frame: pd.DataFrame) -> tuple[pd.Series, str]:
    status = frame.get(
        "canonical_identity_regulatory_status",
        pd.Series("", index=frame.index),
    ).fillna("").astype(str).str.casefold()
    approved = status.isin(
        {"fda_approved_current_or_historical", "drugcentral_fda_approved"}
    )
    claim_flag = frame.get(
        "primary_fda_claim_allowed", pd.Series(False, index=frame.index)
    )
    claim_flag = claim_flag.astype("string").fillna("").str.casefold().isin(
        {"1", "true", "yes", "y"}
    )
    probe_flag = frame.get(
        "probe_sensitivity_only", pd.Series(False, index=frame.index)
    )
    probe_flag = probe_flag.astype("string").fillna("").str.casefold().isin(
        {"1", "true", "yes", "y"}
    )
    eligible = (approved | claim_flag) & ~probe_flag
    if eligible.any():
        return eligible, "canonical_fda_status_or_primary_claim"
    return pd.Series(True, index=frame.index), "no_populated_eligibility_axis"


def fetch_openfda_pk_context(
    model_table: pd.DataFrame,
    *,
    cache_dir: str | Path,
    max_drugs: int = 0,
    label_limit: int = 5,
    timeout: int = 30,
    sleep_sec: float = 0.25,
    reuse_cache: bool = True,
    eligible_only: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    all_drugs = model_table.drop_duplicates("drug_id").copy()
    eligibility_policy = "disabled"
    if eligible_only:
        eligible, eligibility_policy = _openfda_eligible_mask(all_drugs)
        drugs = all_drugs.loc[eligible].copy()
    else:
        drugs = all_drugs
    if max_drugs > 0:
        drugs = drugs.head(max_drugs)
    session = requests.Session()
    rows: list[dict[str, Any]] = []
    statuses: list[dict[str, Any]] = []
    for position, (_, drug_row) in enumerate(drugs.iterrows(), start=1):
        query = _query_from_row(drug_row)
        cache_key = hashlib.sha256(
            "|".join(query.query_names).casefold().encode("utf-8")
        ).hexdigest()[:20]
        cache_path = cache / f"{cache_key}.json"
        records: list[dict[str, Any]] = []
        last_updated = ""
        selected_search = ""
        status = "no_results"
        if reuse_cache and cache_path.exists():
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            records = payload.get("records") or []
            last_updated = str(payload.get("last_updated") or "")
            selected_search = str(payload.get("search") or "")
            status = str(payload.get("status") or "cached")
        else:
            for search in _label_query_strings(query):
                if sleep_sec:
                    time.sleep(sleep_sec)
                try:
                    records, last_updated = _openfda_label_records(
                        session,
                        search,
                        limit=label_limit,
                        timeout=timeout,
                    )
                except Exception as exc:
                    status = f"query_failed:{type(exc).__name__}"
                    LOG.warning("openFDA PK query failed drug=%s error=%s", query.ligand_label, exc)
                    continue
                if records:
                    selected_search = search
                    status = "ok"
                    break
            cache_path.write_text(
                json.dumps(
                    {
                        "drug_id": query.ligand_label,
                        "search": selected_search,
                        "last_updated": last_updated,
                        "status": status,
                        "records": records,
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
        for record in records:
            rows.append(
                _record_context(
                    drug_row=drug_row,
                    record=record,
                    last_updated=last_updated,
                    search=selected_search,
                )
            )
        statuses.append(
            {
                "drug_id": query.ligand_label,
                "status": status,
                "records": len(records),
                "search": selected_search,
            }
        )
        if position % 50 == 0:
            LOG.info("openFDA PK progress %d/%d", position, len(drugs))
    status_frame = pd.DataFrame(statuses)
    status_frame.to_csv(cache.parent / "openfda_pk_query_status.csv", index=False)
    context = finalize_context(pd.DataFrame(rows, columns=PK_CONTEXT_COLUMNS))
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_display": OPENFDA_LABEL_SOURCE,
        "candidate_drugs": int(len(all_drugs)),
        "eligibility_policy": eligibility_policy,
        "drugs_excluded_by_eligibility": int(len(all_drugs) - len(drugs)),
        "drugs_queried": int(len(drugs)),
        "drugs_with_records": int(status_frame["records"].gt(0).sum()) if not status_frame.empty else 0,
        "records": int(len(rows)),
        "rows_with_cmax": int(context["cmax_um"].notna().sum()) if not context.empty else 0,
        "rows_with_free_cmax": int(context["free_cmax_um"].notna().sum()) if not context.empty else 0,
        "confidence_policy": "low until numeric extraction is reviewed against SPL context",
    }
    (cache.parent / "openfda_pk_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return context, manifest
