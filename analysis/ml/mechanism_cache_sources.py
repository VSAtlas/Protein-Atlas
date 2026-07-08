from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


def _json_files(cache_dir: Path) -> list[Path]:
    if not cache_dir.exists():
        return []
    return sorted(path for path in cache_dir.glob("*.json") if path.is_file())


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def export_opentargets_safety_cache(cache_path: str | Path, out_path: str | Path) -> pd.DataFrame:
    """Normalize Atlas report Open Targets safety cache into graph-ready rows."""

    path = Path(cache_path)
    payload = _load_json(path) if path.exists() else {}
    by_uniprot = payload.get("by_uniprot", {}) if isinstance(payload, dict) else {}
    by_gene = payload.get("by_gene", {}) if isinstance(payload, dict) else {}
    rows: list[dict[str, object]] = []

    def add_profile(target_id: str, profile: dict[str, Any], target_kind: str) -> None:
        evidence = profile.get("evidence") or []
        if not evidence and profile.get("safety_buckets"):
            evidence = [
                {
                    "event": bucket,
                    "bucket": bucket,
                    "source": source,
                    "weight": 1.0,
                }
                for bucket in profile.get("safety_buckets", [])
                for source in (profile.get("safety_sources") or ["Open Targets safety cache"])
            ]
        for item in evidence:
            event = str(item.get("event") or item.get("bucket") or "").strip()
            if not event:
                continue
            rows.append(
                {
                    "target_id": target_id,
                    "target_kind": target_kind,
                    "adr": event,
                    "confidence": pd.to_numeric(pd.Series([item.get("weight", 1.0)]), errors="coerce").fillna(1.0).iloc[0],
                    "pubmed_ids": "",
                    "source": item.get("source") or "Open Targets safety cache",
                }
            )

    if isinstance(by_uniprot, dict):
        for target_id, profile in by_uniprot.items():
            if isinstance(profile, dict):
                add_profile(str(target_id), profile, "uniprot")
    if isinstance(by_gene, dict):
        for target_id, profile in by_gene.items():
            if isinstance(profile, dict):
                add_profile(str(target_id), profile, "gene")

    out = pd.DataFrame(rows).drop_duplicates()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def export_reactome_cache(cache_dir: str | Path, out_path: str | Path) -> pd.DataFrame:
    """Normalize cached Reactome UniProt pathway lookups into graph-ready rows."""

    rows: list[dict[str, object]] = []
    for path in _json_files(Path(cache_dir)):
        if not path.name.startswith("pathway_"):
            continue
        try:
            payload = _load_json(path)
        except (OSError, ValueError):
            continue
        request = payload.get("request", {}) if isinstance(payload, dict) else {}
        data = payload.get("data", {}) if isinstance(payload, dict) else {}
        if request.get("source") != "reactome_uniprot_to_pathways":
            continue
        target_id = data.get("uniprot") or request.get("query")
        if not target_id:
            continue
        for pathway in data.get("pathways") or []:
            if not isinstance(pathway, dict):
                continue
            name = str(pathway.get("name") or "").strip()
            if not name:
                continue
            rows.append(
                {
                    "target_id": target_id,
                    "pathway": name,
                    "pathway_id": pathway.get("stId", ""),
                    "adr": "",
                    "source": "Reactome cache",
                }
            )
    out = pd.DataFrame(rows).drop_duplicates()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def export_openfda_ligand_side_effect_cache(cache_path: str | Path, out_path: str | Path) -> pd.DataFrame:
    """Normalize Atlas ligand side-effect cache into SIDER-like drug-ADR rows."""

    path = Path(cache_path)
    payload = _load_json(path) if path.exists() else {}
    entries = payload.get("entries", {}) if isinstance(payload, dict) else {}
    rows: list[dict[str, object]] = []
    iterable: Iterable[Any]
    if isinstance(entries, dict):
        iterable = entries.values()
    elif isinstance(entries, list):
        iterable = entries
    else:
        iterable = []
    for entry in iterable:
        if not isinstance(entry, dict):
            continue
        drug_id = str(entry.get("ligand_base") or entry.get("drug_id") or "").strip().replace(" ", "_")
        if not drug_id:
            continue
        for adr in entry.get("side_effects") or []:
            term = str(adr).strip()
            if term:
                rows.append({"drug_id": drug_id, "adr": term, "source": "openFDA ligand cache"})
        for bucket in entry.get("safety_buckets") or []:
            term = str(bucket).strip()
            if term:
                rows.append({"drug_id": drug_id, "adr": term, "source": "openFDA ligand cache safety bucket"})
    out = pd.DataFrame(rows).drop_duplicates()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def export_cached_mechanism_sources(
    cache_dir: str | Path,
    out_dir: str | Path,
) -> dict[str, str]:
    cache = Path(cache_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    safety_cache = cache / "target_safety_aggregated.json"
    if safety_cache.exists():
        outputs["opentargets_safety"] = str(out / "opentargets_safety_from_cache.tsv")
        export_opentargets_safety_cache(safety_cache, outputs["opentargets_safety"])
    ligand_cache = cache / "ligand_side_effects_openfda.json"
    if ligand_cache.exists():
        outputs["sider"] = str(out / "sider_like_openfda_from_cache.tsv")
        export_openfda_ligand_side_effect_cache(ligand_cache, outputs["sider"])
    outputs["reactome"] = str(out / "reactome_from_cache.tsv")
    export_reactome_cache(cache, outputs["reactome"])
    return outputs
