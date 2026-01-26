from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import logging
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, cast

import requests  # type: ignore[import-untyped]

SCHEMA_VERSION = 1
DEFAULT_ORGANISM = "Homo sapiens"

REACTOME_SEARCH_URL = "https://reactome.org/ContentService/search/query"
REACTOME_PARTICIPANTS_URL = "https://reactome.org/ContentService/data/participants/{pathway_id}"

KEGG_BASE_URL = "https://rest.kegg.jp"

WIKIPATHWAYS_BASE_URL = "https://webservice.wikipathways.org"

PDBE_MAPPING_BASES = [
    "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures",
    "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures",
]
PDBE_LIGAND_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers"

_UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-Z0-9]{10})$"
)
_UNIPROT_TAG_RE = re.compile(r"UniProt:([A-Za-z0-9]+)")


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _normalize_text(text: str) -> str:
    return " ".join(text.strip().lower().split())


def _clean_name(text: str) -> str:
    return re.sub(r"<[^>]+>", "", text)


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _normalize_uniprot(uniprot: str) -> str:
    return uniprot.strip().upper()


def _is_uniprot_accession(value: str) -> bool:
    return bool(_UNIPROT_RE.match(value.strip().upper()))


def _extract_uniprot_from_text(text: str) -> str | None:
    match = _UNIPROT_TAG_RE.search(text)
    if not match:
        return None
    acc = match.group(1).strip().upper()
    if _is_uniprot_accession(acc):
        return acc
    return None

def _load_config(path: str) -> dict[str, Any]:
    module = importlib.import_module("input_and_export_functions")
    loader = cast(Callable[[str], dict[str, Any]], getattr(module, "load_config"))
    return loader(path)


@dataclass
class CachePayload:
    created_at: str
    schema_version: int
    request: dict[str, Any]
    data: dict[str, Any]


class Cache:
    def __init__(
        self, cache_dir: Path, refresh: bool = False, logger: logging.Logger | None = None
    ) -> None:
        self.cache_dir = cache_dir
        self.refresh = refresh
        self.logger = logger or logging.getLogger(__name__)

    def pathway_cache_path(self, source: str, organism: str, query: str) -> Path:
        key = f"{source}|{_normalize_text(organism)}|{_normalize_text(query)}"
        digest = _sha256(key)
        return self.cache_dir / f"pathway_{digest}.json"

    def uniprot_cache_path(self, uniprot: str) -> Path:
        return self.cache_dir / f"uniprot_{_normalize_uniprot(uniprot)}.json"

    def _read_json(self, path: Path, cache_label: str) -> dict[str, Any] | None:
        if self.refresh:
            self.logger.info("Cache bypassed (--refresh) for %s", cache_label)
            return None
        if not path.exists():
            self.logger.info("Cache miss for %s", cache_label)
            return None
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:  # pragma: no cover - corrupted cache
            self.logger.warning("Cache read failed for %s: %s", cache_label, exc)
            return None
        if data.get("schema_version") != SCHEMA_VERSION:
            self.logger.info("Cache schema mismatch for %s", cache_label)
            return None
        self.logger.info("Cache hit for %s", cache_label)
        return data

    def _write_json(self, path: Path, payload: CachePayload) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        payload_dict = {
            "created_at": payload.created_at,
            "schema_version": payload.schema_version,
            "request": payload.request,
            "data": payload.data,
        }
        tmp_path = path.with_suffix(".json.tmp")
        tmp_path.write_text(
            json.dumps(payload_dict, indent=2, sort_keys=True), encoding="utf-8"
        )
        tmp_path.replace(path)

    def read_pathway(self, source: str, organism: str, query: str) -> dict[str, Any] | None:
        path = self.pathway_cache_path(source, organism, query)
        return self._read_json(path, "pathway")

    def write_pathway(
        self,
        source: str,
        organism: str,
        query: str,
        data: dict[str, Any],
    ) -> None:
        path = self.pathway_cache_path(source, organism, query)
        payload = CachePayload(
            created_at=_now_iso(),
            schema_version=SCHEMA_VERSION,
            request={"source": source, "organism": organism, "query": query},
            data=data,
        )
        self._write_json(path, payload)

    def read_uniprot(self, uniprot: str) -> dict[str, Any] | None:
        path = self.uniprot_cache_path(uniprot)
        return self._read_json(path, f"uniprot:{uniprot}")

    def write_uniprot(self, uniprot: str, data: dict[str, Any]) -> None:
        path = self.uniprot_cache_path(uniprot)
        payload = CachePayload(
            created_at=_now_iso(),
            schema_version=SCHEMA_VERSION,
            request={"uniprot": _normalize_uniprot(uniprot)},
            data=data,
        )
        self._write_json(path, payload)


class HttpClient:
    def __init__(self, timeout: int = 30) -> None:
        self.timeout = timeout

    def get_json(self, url: str, params: dict[str, str] | None = None) -> Any:
        resp = requests.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    def get_text(self, url: str, params: dict[str, str] | None = None) -> str:
        resp = requests.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        return resp.text


def _matches_organism(result: dict[str, Any], organism: str) -> bool:
    org_norm = _normalize_text(organism)
    species = result.get("species") or result.get("speciesName")
    if isinstance(species, str):
        return _normalize_text(species) == org_norm
    if isinstance(species, list):
        for entry in species:
            if isinstance(entry, dict):
                name = entry.get("name")
                if name and _normalize_text(str(name)) == org_norm:
                    return True
            elif isinstance(entry, str) and _normalize_text(entry) == org_norm:
                return True
    return False


def _pick_best_pathway(
    results: Iterable[dict[str, Any]], query: str
) -> dict[str, Any] | None:
    query_norm = _normalize_text(query)
    best = None
    best_key: tuple[int, float, str] | None = None
    for res in results:
        name_raw = str(res.get("name") or res.get("displayName") or res.get("id") or "")
        name = _clean_name(name_raw)
        name_norm = _normalize_text(name)
        exact = 0 if name_norm == query_norm else 1
        score = float(res.get("score") or 0.0)
        key = (exact, -score, name_norm)
        if best_key is None or key < best_key:
            best_key = key
            best = res
    return best


def _extract_uniprots_from_reactome(
    participants: Iterable[dict[str, Any]],
) -> list[str]:
    uniprots: list[str] = []
    for item in participants:
        if not isinstance(item, dict):
            continue
        db = str(item.get("databaseName") or item.get("database") or "")
        identifier = item.get("identifier") or item.get("identifierId") or item.get("id")
        if identifier is None:
            identifier = ""
        acc = str(identifier).strip()
        if acc and ("uniprot" in db.lower() or _is_uniprot_accession(acc)):
            if _is_uniprot_accession(acc):
                uniprots.append(_normalize_uniprot(acc))
                continue

        display = str(item.get("displayName") or "")
        from_display = _extract_uniprot_from_text(display)
        if from_display:
            uniprots.append(from_display)

        ref_entities = item.get("refEntities")
        if isinstance(ref_entities, list):
            for ref in ref_entities:
                if not isinstance(ref, dict):
                    continue
                ref_id = ref.get("identifier")
                if isinstance(ref_id, str) and _is_uniprot_accession(ref_id):
                    uniprots.append(_normalize_uniprot(ref_id))
                    continue
                ref_display = str(ref.get("displayName") or "")
                from_ref = _extract_uniprot_from_text(ref_display)
                if from_ref:
                    uniprots.append(from_ref)
    return uniprots


def _resolve_reactome(
    pathway_query: str, organism: str, http: HttpClient, logger: logging.Logger
) -> tuple[str, str, list[str]]:
    params = {"query": pathway_query, "species": organism, "types": "Pathway"}
    resp = http.get_json(REACTOME_SEARCH_URL, params=params)
    raw_results: list[Any] = []
    if isinstance(resp, dict):
        raw_results = resp.get("results") or resp.get("entries") or []
    elif isinstance(resp, list):
        raw_results = resp

    results: list[dict[str, Any]] = []
    if isinstance(raw_results, list):
        for item in raw_results:
            if not isinstance(item, dict):
                continue
            entries = item.get("entries")
            if isinstance(entries, list):
                results.extend([e for e in entries if isinstance(e, dict)])
            else:
                results.append(item)

    filtered = [
        r for r in results if isinstance(r, dict) and _matches_organism(r, organism)
    ]
    best = _pick_best_pathway(filtered or results, pathway_query)
    if not best:
        return ("", "", [])
    pathway_id = str(best.get("stId") or best.get("id") or "")
    if not pathway_id:
        logger.warning("Reactome pathway ID missing for query: %s", pathway_query)
        return ("", "", [])
    pathway_name = _clean_name(
        str(best.get("name") or best.get("displayName") or pathway_id)
    )
    logger.info("Reactome pathway match: %s (%s)", pathway_name, pathway_id)

    participants_url = REACTOME_PARTICIPANTS_URL.format(pathway_id=pathway_id)
    participants_resp = http.get_json(participants_url)
    participants: list[dict[str, Any]] = []
    if isinstance(participants_resp, list):
        participants = participants_resp
    elif isinstance(participants_resp, dict):
        participants = participants_resp.get("participants") or []

    uniprots = _extract_uniprots_from_reactome(participants)
    logger.info(
        "Reactome participants: %d (UniProt: %d)",
        len(participants),
        len(uniprots),
    )
    return (pathway_id, pathway_name, uniprots)


def _kegg_organism_code(
    organism: str, http: HttpClient, logger: logging.Logger
) -> str | None:
    if organism.islower() and len(organism) <= 4:
        return organism
    resp = http.get_text(f"{KEGG_BASE_URL}/list/organism")
    org_norm = _normalize_text(organism)
    best_code = None
    for line in resp.splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        code = parts[1].strip() if len(parts) > 1 else ""
        name = parts[2].strip() if len(parts) > 2 else ""
        if _normalize_text(name).startswith(org_norm):
            best_code = code
            break
    if not best_code:
        logger.warning("KEGG organism code not found for %s", organism)
    return best_code


def _kegg_pick_pathway(
    lines: list[str], query: str, org_code: str | None
) -> tuple[str, str]:
    query_norm = _normalize_text(query)
    best_id = ""
    best_name = ""
    best_key: tuple[int, int, str] | None = None
    for line in lines:
        parts = line.split("\t", 1)
        if len(parts) < 2:
            continue
        pathway_id, name = parts[0].strip(), parts[1].strip()
        if org_code and not pathway_id.startswith(f"path:{org_code}"):
            continue
        clean_name = name.split(" - ")[0].strip()
        clean_norm = _normalize_text(clean_name)
        exact = 0 if clean_norm == query_norm else 1
        contains = 0 if query_norm in clean_norm else 1
        key = (exact, contains, clean_norm)
        if best_key is None or key < best_key:
            best_key = key
            best_id = pathway_id
            best_name = clean_name
    return (best_id, best_name)


def _extract_uniprots_from_kegg_link(text: str) -> list[str]:
    uniprots = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        target = parts[1].strip()
        if ":" in target:
            target = target.split(":", 1)[1]
        if _is_uniprot_accession(target):
            uniprots.append(_normalize_uniprot(target))
    return uniprots


def _resolve_kegg(
    pathway_query: str, organism: str, http: HttpClient, logger: logging.Logger
) -> tuple[str, str, list[str]]:
    org_code = _kegg_organism_code(organism, http, logger)
    search_resp = http.get_text(f"{KEGG_BASE_URL}/find/pathway/{pathway_query}")
    lines = [line for line in search_resp.splitlines() if line.strip()]
    pathway_id, pathway_name = _kegg_pick_pathway(lines, pathway_query, org_code)
    if not pathway_id:
        return ("", "", [])
    logger.info("KEGG pathway match: %s (%s)", pathway_name, pathway_id)

    link_resp = http.get_text(f"{KEGG_BASE_URL}/link/uniprot/{pathway_id}")
    uniprots = _extract_uniprots_from_kegg_link(link_resp)
    if uniprots:
        logger.info("KEGG participants: UniProt %d", len(uniprots))
        return (pathway_id, pathway_name, uniprots)

    gene_resp = http.get_text(f"{KEGG_BASE_URL}/link/genes/{pathway_id}")
    gene_lines = [line for line in gene_resp.splitlines() if line.strip()]
    for line in gene_lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        gene_id = parts[1].strip()
        conv_resp = http.get_text(f"{KEGG_BASE_URL}/conv/uniprot/{gene_id}")
        if not conv_resp.strip():
            logger.warning("KEGG gene %s missing UniProt mapping", gene_id)
            continue
        uniprots.extend(_extract_uniprots_from_kegg_link(conv_resp))

    logger.info("KEGG participants: UniProt %d", len(uniprots))
    return (pathway_id, pathway_name, uniprots)


def _extract_uniprots_from_wikipathways(pathway: dict[str, Any]) -> list[str]:
    uniprots: list[str] = []
    elements = pathway.get("elements") or pathway.get("element") or []
    if isinstance(elements, dict):
        elements = [elements]
    for element in elements:
        if not isinstance(element, dict):
            continue
        source = str(element.get("xrefDataSource") or "")
        xref_id = element.get("xrefId")
        if source.lower().startswith("uniprot") and xref_id:
            acc = str(xref_id).strip()
            if _is_uniprot_accession(acc):
                uniprots.append(_normalize_uniprot(acc))
        xrefs = element.get("xref") or element.get("xrefs")
        if isinstance(xrefs, list):
            for xref in xrefs:
                if not isinstance(xref, dict):
                    continue
                source = str(xref.get("datasource") or xref.get("dataSource") or "")
                xref_id = xref.get("id")
                if source.lower().startswith("uniprot") and xref_id:
                    acc = str(xref_id).strip()
                    if _is_uniprot_accession(acc):
                        uniprots.append(_normalize_uniprot(acc))
    return uniprots


def _resolve_wikipathways(
    pathway_query: str, organism: str, http: HttpClient, logger: logging.Logger
) -> tuple[str, str, list[str]]:
    params = {"query": pathway_query, "species": organism, "format": "json"}
    resp = http.get_json(f"{WIKIPATHWAYS_BASE_URL}/findPathwaysByText", params=params)
    results: list[dict[str, Any]] = []
    if isinstance(resp, dict):
        raw_results = resp.get("pathways") or resp.get("result", {}).get("pathways") or []
        if isinstance(raw_results, list):
            results = [r for r in raw_results if isinstance(r, dict)]
    elif isinstance(resp, list):
        results = [r for r in resp if isinstance(r, dict)]
    filtered = [r for r in results if isinstance(r, dict) and _matches_organism(r, organism)]
    best = _pick_best_pathway(filtered or results, pathway_query)
    if not best:
        return ("", "", [])
    pathway_id = str(best.get("id") or best.get("pwId") or "")
    pathway_name = str(best.get("name") or pathway_id)
    logger.info("WikiPathways match: %s (%s)", pathway_name, pathway_id)

    pathway_resp = http.get_json(
        f"{WIKIPATHWAYS_BASE_URL}/getPathway", params={"pwId": pathway_id, "format": "json"}
    )
    pathway: dict[str, Any] = {}
    if isinstance(pathway_resp, dict):
        pathway = pathway_resp.get("pathway") or pathway_resp.get("result") or {}

    uniprots = _extract_uniprots_from_wikipathways(pathway)
    logger.info("WikiPathways participants: UniProt %d", len(uniprots))
    return (pathway_id, pathway_name, uniprots)


def resolve_uniprots(
    pathway_query: str,
    source: str,
    organism: str,
    cache: Cache,
    http: HttpClient,
) -> list[str]:
    cached = cache.read_pathway(source, organism, pathway_query)
    if cached:
        data = cached.get("data", {})
        return sorted({_normalize_uniprot(p) for p in data.get("uniprots", [])})

    logger = cache.logger
    source_norm = source.lower()
    if source_norm == "reactome":
        pathway_id, pathway_name, uniprots = _resolve_reactome(
            pathway_query, organism, http, logger
        )
    elif source_norm == "kegg":
        pathway_id, pathway_name, uniprots = _resolve_kegg(
            pathway_query, organism, http, logger
        )
    elif source_norm == "wikipathways":
        pathway_id, pathway_name, uniprots = _resolve_wikipathways(
            pathway_query, organism, http, logger
        )
    else:
        raise ValueError(f"Unknown source: {source}")

    unique = sorted({_normalize_uniprot(p) for p in uniprots})
    cache.write_pathway(
        source,
        organism,
        pathway_query,
        {"pathway_id": pathway_id, "pathway_name": pathway_name, "uniprots": unique},
    )
    logger.info("Resolved UniProts: %d", len(unique))
    return unique


def _pdbe_best_structures(
    uniprot: str, http: HttpClient, logger: logging.Logger
) -> list[dict[str, Any]]:
    last_error: Exception | None = None
    for base in PDBE_MAPPING_BASES:
        url = f"{base}/{_normalize_uniprot(uniprot)}"
        try:
            resp = http.get_json(url)
        except Exception as exc:  # pragma: no cover - network errors
            last_error = exc
            continue
        if isinstance(resp, dict):
            for key in (
                _normalize_uniprot(uniprot),
                _normalize_uniprot(uniprot).lower(),
            ):
                if key in resp:
                    return resp.get(key) or []
            if "best_structures" in resp:
                return resp.get("best_structures") or []
        elif isinstance(resp, list):
            return resp
    if last_error:
        logger.warning("PDBe mapping failed for %s: %s", uniprot, last_error)
    return []


def _normalize_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    pdb_id = str(candidate.get("pdb_id") or candidate.get("pdb") or "").strip()
    pdb_id = pdb_id.upper()
    resolution = candidate.get("resolution")
    if resolution is not None:
        try:
            resolution = float(resolution)
        except Exception:
            resolution = None
    method = candidate.get("experimental_method") or candidate.get("method")
    coverage = candidate.get("coverage")
    if coverage is not None:
        try:
            coverage = float(coverage)
        except Exception:
            coverage = None
    unp_start = candidate.get("unp_start") or candidate.get("start")
    unp_end = candidate.get("unp_end") or candidate.get("end")
    if unp_start is not None:
        try:
            unp_start = int(unp_start)
        except Exception:
            unp_start = None
    if unp_end is not None:
        try:
            unp_end = int(unp_end)
        except Exception:
            unp_end = None
    mutation_count = candidate.get("mutation_count") or candidate.get("mutations")
    if mutation_count is not None:
        try:
            mutation_count = int(mutation_count)
        except Exception:
            mutation_count = None
    has_ligand = candidate.get("has_ligand")
    if isinstance(has_ligand, str):
        has_ligand = has_ligand.lower() in {"true", "yes", "1"}
    elif not isinstance(has_ligand, bool):
        has_ligand = None
    return {
        "pdb_id": pdb_id,
        "resolution": resolution,
        "method": method,
        "coverage": coverage,
        "unp_start": unp_start,
        "unp_end": unp_end,
        "mutation_count": mutation_count,
        "has_ligand": has_ligand,
    }


def _pdbe_has_ligand(pdb_id: str, http: HttpClient) -> bool | None:
    url = f"{PDBE_LIGAND_URL}/{pdb_id.lower()}"
    try:
        resp = http.get_json(url)
    except Exception:  # pragma: no cover - network errors
        return None
    if not isinstance(resp, dict):
        return None
    data = resp.get(pdb_id.lower()) or resp.get(pdb_id.upper())
    if data is None:
        return None
    return bool(data)


def _candidate_coverage_score(candidate: dict[str, Any]) -> float:
    coverage = candidate.get("coverage")
    if isinstance(coverage, (int, float)):
        return float(coverage)
    unp_start = candidate.get("unp_start")
    unp_end = candidate.get("unp_end")
    if isinstance(unp_start, int) and isinstance(unp_end, int):
        return float(max(0, unp_end - unp_start + 1))
    return -1.0


def _method_rank(method: str | None) -> int:
    if not method:
        return 3
    method_norm = method.lower()
    if "x-ray" in method_norm or "xray" in method_norm:
        return 0
    if "cryo" in method_norm:
        return 1
    if "nmr" in method_norm:
        return 2
    return 3


def _candidate_sort_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    has_ligand = candidate.get("has_ligand") is True
    ligand_rank = 0 if has_ligand else 1
    coverage_score = _candidate_coverage_score(candidate)
    resolution = candidate.get("resolution")
    if isinstance(resolution, (int, float)):
        resolution_present = 0
        resolution_val = float(resolution)
        method_rank = 0
    else:
        resolution_present = 1
        resolution_val = 99.0
        method_rank = _method_rank(candidate.get("method"))
    mutation = candidate.get("mutation_count")
    mutation_rank = int(mutation) if isinstance(mutation, int) else 999
    pdb_id = str(candidate.get("pdb_id") or "")
    return (
        ligand_rank,
        -coverage_score,
        resolution_present,
        resolution_val,
        method_rank,
        mutation_rank,
        pdb_id,
    )


def map_uniprot_to_pdb(
    uniprot: str, cache: Cache, http: HttpClient
) -> list[dict[str, Any]]:
    cached = cache.read_uniprot(uniprot)
    if cached:
        data = cached.get("data", {})
        return data.get("candidates", [])

    logger = cache.logger
    raw_candidates = _pdbe_best_structures(uniprot, http, logger)
    candidates = [
        _normalize_candidate(c)
        for c in raw_candidates
        if isinstance(c, dict) and (c.get("pdb_id") or c.get("pdb"))
    ]
    candidates = [c for c in candidates if c.get("pdb_id")]

    to_enrich = [c for c in candidates if c.get("has_ligand") is None]
    if to_enrich:
        for candidate in sorted(to_enrich, key=_candidate_sort_key)[:3]:
            pdb_id = candidate.get("pdb_id")
            if not pdb_id:
                continue
            has_ligand = _pdbe_has_ligand(str(pdb_id), http)
            if has_ligand is not None:
                candidate["has_ligand"] = has_ligand

    cache.write_uniprot(
        uniprot,
        {
            "snapshot_at": _now_iso(),
            "candidates": candidates,
        },
    )
    logger.info("PDBe candidates for %s: %d", _normalize_uniprot(uniprot), len(candidates))
    return candidates


def select_representatives(
    uniprot: str, candidates: list[dict[str, Any]], max_n: int
) -> list[str]:
    max_n = max(1, min(3, max_n))
    sorted_candidates = sorted(candidates, key=_candidate_sort_key)
    selected: list[str] = []
    seen: set[str] = set()
    for candidate in sorted_candidates:
        pdb_id = str(candidate.get("pdb_id") or "").strip().upper()
        if len(pdb_id) != 4:
            continue
        if pdb_id in seen:
            continue
        selected.append(pdb_id)
        seen.add(pdb_id)
        if len(selected) >= max_n:
            break
    return selected


def write_resolved_pdbs(pdb_ids: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = "\n".join([p.strip().upper() for p in pdb_ids if p.strip()]) + "\n"
    output_path.write_text(payload, encoding="utf-8")


def _configure_logging(level: str) -> logging.Logger:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    return logging.getLogger("pathway_resolver")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Resolve pathway query to PDB IDs.")
    parser.add_argument("pathway", nargs="?", help="Pathway query text.")
    parser.add_argument("--pathway", dest="pathway_flag", help="Pathway query text.")
    parser.add_argument(
        "--source",
        default="reactome",
        choices=["reactome", "kegg", "wikipathways"],
        help="Pathway data source.",
    )
    parser.add_argument("--organism", help="Organism name override.")
    parser.add_argument(
        "--output",
        default="pathways/resolved_pdbs.txt",
        help="Output file path.",
    )
    parser.add_argument(
        "--cache-dir",
        default="pathways/cache",
        help="Cache directory.",
    )
    parser.add_argument(
        "--max-pdbs-per-uniprot",
        type=int,
        default=1,
        choices=[1, 2, 3],
        help="Max PDBs per UniProt (1-3).",
    )
    parser.add_argument("--refresh", action="store_true", help="Bypass cache.")
    parser.add_argument(
        "--config",
        default="config.txt",
        help="Config file path.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level.",
    )
    return parser


def run(argv: list[str], http_client: HttpClient | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    pathway_query = args.pathway_flag or args.pathway
    if not pathway_query:
        parser.error("Pathway query is required (positional or --pathway).")

    logger = _configure_logging(args.log_level)
    cfg = _load_config(args.config)
    organism = args.organism or cfg.get("PATHWAY_ORGANISM") or DEFAULT_ORGANISM

    output_path = Path(args.output)
    cache_dir = Path(args.cache_dir)
    root = Path(cfg.get("OVERALL_DIR") or Path.cwd())
    if not output_path.is_absolute():
        output_path = root / output_path
    if not cache_dir.is_absolute():
        cache_dir = root / cache_dir

    cache = Cache(cache_dir=cache_dir, refresh=args.refresh, logger=logger)
    http = http_client or HttpClient()

    uniprots = resolve_uniprots(pathway_query, args.source, organism, cache, http)
    if not uniprots:
        logger.warning("No UniProt accessions resolved for %s", pathway_query)

    selected_pdbs: list[str] = []
    seen_pdbs: set[str] = set()
    for uniprot in sorted(uniprots):
        candidates = map_uniprot_to_pdb(uniprot, cache, http)
        selected = select_representatives(
            uniprot, candidates, args.max_pdbs_per_uniprot
        )
        for pdb_id in selected:
            if pdb_id in seen_pdbs:
                logger.warning(
                    "Duplicate PDB %s for UniProt %s (already selected)", pdb_id, uniprot
                )
                continue
            selected_pdbs.append(pdb_id)
            seen_pdbs.add(pdb_id)

    write_resolved_pdbs(selected_pdbs, output_path)
    logger.info("Resolved PDBs written: %s", output_path)
    return 0


def main() -> None:
    raise SystemExit(run(sys.argv[1:]))


if __name__ == "__main__":
    main()
