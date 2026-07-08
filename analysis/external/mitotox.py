from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd


MITOTOX_API = "https://www.mitotox.org/api"


def _get_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "AtlasAnalysis/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def _collect_paginated(endpoint: str, *, sleep_sec: float = 0.1) -> pd.DataFrame:
    url: str | None = f"{MITOTOX_API}/{endpoint.strip('/')}"
    rows: list[dict[str, Any]] = []
    while url:
        payload = _get_json(url)
        rows.extend(payload.get("results", []))
        url = payload.get("next")
        if url:
            time.sleep(max(0.0, sleep_sec))
    return pd.DataFrame(rows)


def stage_mitotox(out_dir: str | Path, *, sleep_sec: float = 0.1) -> dict[str, Any]:
    """Download MitoTox API tables and create Atlas-normalized mechanism evidence."""

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    compounds = _collect_paginated("compounds/list", sleep_sec=sleep_sec)
    targets = _collect_paginated("targets/list", sleep_sec=sleep_sec)
    records = _collect_paginated("records/list", sleep_sec=sleep_sec)

    compounds.to_csv(out / "mitotox_compounds.csv", index=False)
    targets.to_csv(out / "mitotox_targets.csv", index=False)
    records.to_csv(out / "mitotox_records.csv", index=False)

    norm = normalize_mitotox_records(compounds, targets, records)
    norm.to_csv(out / "mitotox_mechanism_evidence.tsv", sep="\t", index=False)
    manifest = {
        "source": "MitoTox",
        "license": "CC BY-NC 4.0 per https://www.mitotox.org/api/",
        "api": MITOTOX_API,
        "compound_rows": int(len(compounds)),
        "target_rows": int(len(targets)),
        "record_rows": int(len(records)),
        "normalized_rows": int(len(norm)),
        "outputs": {
            "compounds": str(out / "mitotox_compounds.csv"),
            "targets": str(out / "mitotox_targets.csv"),
            "records": str(out / "mitotox_records.csv"),
            "normalized": str(out / "mitotox_mechanism_evidence.tsv"),
        },
        "label_policy": (
            "MitoTox result=true creates positive mitochondrial mechanism evidence. "
            "result=false is retained as measured negative evidence for the specific assay context, not as a global ADR negative."
        ),
    }
    (out / "mitotox_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def _split_id_label(value: Any) -> tuple[str, str]:
    text = str(value or "").strip()
    if ":" not in text:
        return text, ""
    left, right = text.split(":", 1)
    return left.strip(), right.strip()


def normalize_mitotox_records(compounds: pd.DataFrame, targets: pd.DataFrame, records: pd.DataFrame) -> pd.DataFrame:
    compound_map = compounds.copy()
    if not compound_map.empty:
        compound_map["compound_id"] = compound_map.get("compound_ID", pd.Series(dtype=object)).astype(str)
        compound_map = compound_map.set_index("compound_id", drop=False)
    target_map = targets.copy()
    if not target_map.empty:
        target_map["target_source_id"] = target_map.get("targ_id", pd.Series(dtype=object)).astype(str)
        target_map = target_map.set_index("target_source_id", drop=False)

    rows: list[dict[str, Any]] = []
    for row in records.to_dict("records"):
        compound_id, compound_name = _split_id_label(row.get("compound"))
        target_id, target_name = _split_id_label(row.get("targ"))
        compound = compound_map.loc[compound_id].to_dict() if compound_id in compound_map.index else {}
        target = target_map.loc[target_id].to_dict() if target_id in target_map.index else {}
        result = row.get("result")
        if result is True:
            label_state: int | None = 1
            status = "mitochondrial_toxicity_positive"
        elif result is False:
            label_state = 0
            status = "mitochondrial_toxicity_measured_negative"
        else:
            label_state = None
            status = "unknown_or_unreported"
        rows.append(
            {
                "drug_id": compound.get("cid") or compound_id,
                "compound_namespace": "PubChemCID" if compound.get("cid") else "MitoTox",
                "drug_name": compound.get("name") or compound_name,
                "target_id": target.get("Uniport_ID") or target.get("gene_symbol") or target_id,
                "target_namespace": "UniProt" if target.get("Uniport_ID") else "GeneSymbol" if target.get("gene_symbol") else "MitoTox",
                "gene_symbol": target.get("gene_symbol"),
                "target_name": target.get("full_name") or target_name,
                "assay_id": row.get("rec_id"),
                "activity_outcome": label_state,
                "activity_type": row.get("method"),
                "activity_relation": "",
                "activity_units": "",
                "source": "MitoTox",
                "source_row_id": row.get("rec_id"),
                "label_status": status,
                "evidence_namespace": "drug_target_adr_mechanism",
                "mitotox_function": row.get("func"),
                "mitotox_target_action": row.get("targ_act"),
                "mitotox_function_action": row.get("func_act"),
                "mitotox_dose": row.get("dose"),
                "mitotox_time": row.get("time"),
                "mitotox_species": row.get("species"),
                "mitotox_model": row.get("animal_model"),
                "document_ids": row.get("paper"),
            }
        )
    return pd.DataFrame(rows)
