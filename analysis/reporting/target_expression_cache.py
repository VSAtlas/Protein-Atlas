from __future__ import annotations

import csv
import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, cast

import requests  # type: ignore[import-untyped]

HPA_SEARCH_URL = "https://www.proteinatlas.org/api/search_download.php"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
TARGET_EXPRESSION_CACHE_VERSION = 1
DEFAULT_CACHE_REL = Path("pathways") / "cache" / "target_expression_hpa.json"

LOG = logging.getLogger(__name__)


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _normalize_text(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def _normalize_upper(value: Any) -> str:
    return _normalize_text(value).upper()


def _normalize_key(value: Any) -> str:
    return _normalize_text(value).casefold()


def cache_path(repo_root: Path) -> Path:
    return repo_root / DEFAULT_CACHE_REL


def load_target_expression_cache(repo_root: Path) -> Dict[str, Any]:
    path = cache_path(repo_root)
    if not path.exists():
        return {
            "version": TARGET_EXPRESSION_CACHE_VERSION,
            "source": "Human Protein Atlas",
            "entries": {},
            "by_uniprot": {},
            "by_gene": {},
        }
    payload = json.loads(path.read_text(encoding="utf-8"))
    return cast(Dict[str, Any], payload if isinstance(payload, dict) else {})


def lookup_target_expression_entry(
    repo_root: Path,
    *,
    uniprots: Iterable[str],
    gene_symbols: Iterable[str] | None = None,
) -> Dict[str, Any]:
    payload = load_target_expression_cache(repo_root)
    entries = cast(Dict[str, Any], payload.get("entries") or {})
    by_uniprot = cast(Dict[str, Any], payload.get("by_uniprot") or {})
    by_gene = cast(Dict[str, Any], payload.get("by_gene") or {})
    for uniprot in sorted({_normalize_upper(value) for value in uniprots if _normalize_text(value)}):
        entry_key = _normalize_key(by_uniprot.get(uniprot))
        if entry_key and isinstance(entries.get(entry_key), dict):
            return cast(Dict[str, Any], entries[entry_key])
    for gene in sorted({_normalize_upper(value) for value in (gene_symbols or []) if _normalize_text(value)}):
        entry_key = _normalize_key(by_gene.get(gene))
        if entry_key and isinstance(entries.get(entry_key), dict):
            return cast(Dict[str, Any], entries[entry_key])
    return {}


def _dedupe(values: Iterable[Any], *, limit: int = 0) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        clean = _normalize_text(value)
        key = clean.casefold()
        if not clean or key in seen:
            continue
        seen.add(key)
        out.append(clean)
        if limit and len(out) >= limit:
            break
    return out


def _extract_uniprot_gene_symbols(payload: Mapping[str, Any]) -> list[str]:
    genes: list[str] = []
    for gene in cast(Sequence[Any], payload.get("genes") or []):
        if not isinstance(gene, Mapping):
            continue
        gene_name = cast(Mapping[str, Any], gene.get("geneName") or {})
        genes.append(_normalize_text(gene_name.get("value")))
        for synonym in cast(Sequence[Any], gene.get("synonyms") or []):
            if isinstance(synonym, Mapping):
                genes.append(_normalize_text(synonym.get("value")))
    return _dedupe(genes, limit=12)


def resolve_uniprot_gene_symbols(
    accession: str,
    *,
    session: Optional[requests.Session] = None,
    timeout: int = 20,
) -> list[str]:
    clean = _normalize_upper(accession)
    if not clean:
        return []
    client = session or requests.Session()
    response = client.get(UNIPROT_URL.format(accession=clean), timeout=timeout)
    if response.status_code == 404:
        return []
    response.raise_for_status()
    payload = response.json()
    organism = _normalize_text(cast(Mapping[str, Any], payload.get("organism") or {}).get("scientificName"))
    if organism and organism.casefold() != "homo sapiens":
        return []
    return _extract_uniprot_gene_symbols(cast(Mapping[str, Any], payload))


def _hpa_search(
    gene_symbol: str,
    *,
    session: Optional[requests.Session] = None,
    timeout: int = 20,
) -> list[Dict[str, Any]]:
    client = session or requests.Session()
    response = client.get(
        HPA_SEARCH_URL,
        params={
            "search": gene_symbol,
            "format": "json",
            "columns": "g,gs,rnats,rnatd,rnatss,rnatsm",
            "compress": "no",
        },
        timeout=timeout,
    )
    if response.status_code == 404:
        return []
    response.raise_for_status()
    payload = response.json()
    return cast(list[Dict[str, Any]], payload if isinstance(payload, list) else [])


def _select_hpa_record(records: Sequence[Mapping[str, Any]], gene_symbols: Sequence[str]) -> Mapping[str, Any]:
    wanted = {_normalize_upper(symbol) for symbol in gene_symbols if _normalize_text(symbol)}
    for record in records:
        gene = _normalize_upper(record.get("Gene"))
        if gene in wanted:
            return record
    for record in records:
        synonyms = {
            _normalize_upper(value)
            for value in cast(Sequence[Any], record.get("Gene synonym") or [])
            if _normalize_text(value)
        }
        if synonyms & wanted:
            return record
    return records[0] if records else {}


def _expression_level_and_score(record: Mapping[str, Any]) -> tuple[str, float]:
    specificity = _normalize_text(record.get("RNA tissue specificity")).casefold()
    distribution = _normalize_text(record.get("RNA tissue distribution")).casefold()
    ntpm = cast(Mapping[str, Any], record.get("RNA tissue specific nTPM") or {})
    max_ntpm = 0.0
    for raw_value in ntpm.values():
        try:
            max_ntpm = max(max_ntpm, float(str(raw_value)))
        except (TypeError, ValueError):
            continue
    if "not detected" in distribution:
        return "not detected", 0.05
    if max_ntpm >= 50 or "detected in all" in distribution or "detected in many" in distribution:
        return "high", 0.85
    if "detected in some" in distribution or max_ntpm >= 10:
        return "medium", 0.67
    if "detected in single" in distribution or "tissue enriched" in specificity or max_ntpm > 0:
        return "low", 0.33
    return "not detected", 0.05


def _expression_summary(record: Mapping[str, Any]) -> str:
    specificity = _normalize_text(record.get("RNA tissue specificity"))
    distribution = _normalize_text(record.get("RNA tissue distribution"))
    ntpm = cast(Mapping[str, Any], record.get("RNA tissue specific nTPM") or {})
    enriched = ", ".join(
        f"{tissue}: {value} nTPM"
        for tissue, value in sorted(ntpm.items(), key=lambda item: str(item[0]).lower())[:5]
    )
    parts = [part for part in (specificity, distribution, enriched) if part]
    return "; ".join(parts)


def fetch_target_expression_entry(
    *,
    uniprots: Iterable[str],
    gene_symbols: Iterable[str] | None = None,
    session: Optional[requests.Session] = None,
    timeout: int = 20,
    sleep_sec: float = 0.0,
) -> Dict[str, Any]:
    client = session or requests.Session()
    clean_uniprots = _dedupe((_normalize_upper(value) for value in uniprots), limit=12)
    genes = _dedupe(gene_symbols or [], limit=12)
    for uniprot in clean_uniprots:
        if sleep_sec > 0:
            time.sleep(sleep_sec)
        try:
            genes.extend(resolve_uniprot_gene_symbols(uniprot, session=client, timeout=timeout))
        except Exception as exc:
            LOG.warning("UniProt gene lookup failed accession=%s error=%s", uniprot, exc)
    genes = _dedupe(genes, limit=12)

    selected_record: Mapping[str, Any] = {}
    searched_gene = ""
    for gene in genes:
        if sleep_sec > 0:
            time.sleep(sleep_sec)
        try:
            records = _hpa_search(gene, session=client, timeout=timeout)
        except Exception as exc:
            LOG.warning("HPA expression query failed gene=%s error=%s", gene, exc)
            continue
        record = _select_hpa_record(records, [gene, *genes])
        if record:
            selected_record = record
            searched_gene = gene
            break

    if not selected_record:
        return {
            "uniprot_accessions": clean_uniprots,
            "gene_symbols": genes,
            "hpa_query_status": "no_results" if genes else "no_query",
        }

    level, score = _expression_level_and_score(selected_record)
    gene = _normalize_upper(selected_record.get("Gene")) or _normalize_upper(searched_gene)
    return {
        "gene": gene,
        "gene_symbols": genes,
        "uniprot_accessions": clean_uniprots,
        "target_tissue_expression": level,
        "target_tissue_expression_label": level,
        "target_tissue_expression_score": f"{score:.3f}",
        "target_tissue_expression_source": "Human Protein Atlas",
        "target_tissue_expression_summary": _expression_summary(selected_record),
        "hpa_query_gene": searched_gene,
        "hpa_query_status": "ok",
    }


def build_cache_payload(entries: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    payload_entries: Dict[str, Any] = {}
    by_uniprot: Dict[str, str] = {}
    by_gene: Dict[str, str] = {}
    for entry in entries:
        gene = _normalize_upper(entry.get("gene") or next(iter(entry.get("gene_symbols") or []), ""))
        primary = _normalize_key(gene or next(iter(entry.get("uniprot_accessions") or []), ""))
        if not primary:
            continue
        payload_entries[primary] = dict(entry)
        for uniprot in cast(Sequence[Any], entry.get("uniprot_accessions") or []):
            clean = _normalize_upper(uniprot)
            if clean:
                by_uniprot[clean] = primary
        for symbol in cast(Sequence[Any], entry.get("gene_symbols") or []):
            clean = _normalize_upper(symbol)
            if clean:
                by_gene[clean] = primary
        if gene:
            by_gene[gene] = primary
    return {
        "version": TARGET_EXPRESSION_CACHE_VERSION,
        "source": "Human Protein Atlas",
        "generated_at": _now_iso(),
        "entries": payload_entries,
        "by_uniprot": by_uniprot,
        "by_gene": by_gene,
    }


def write_cache(repo_root: Path, payload: Mapping[str, Any]) -> Path:
    path = cache_path(repo_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def pdb_ids_from_heatmap_csv(path: Path) -> list[str]:
    pdb_ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            pdb_id = _normalize_upper(row.get("pdb_id"))
            if not pdb_id:
                target_id = _normalize_text(row.get("target_id"))
                pdb_id = _normalize_upper(target_id.split("|", 1)[0])
            if pdb_id:
                pdb_ids.add(pdb_id)
    return sorted(pdb_ids)
