"""Query and rank RCSB PDB entries for docking panel selection.

Defaults:
- cache dir: ``chemdb/.cache/pdb_query`` (relative to this module)
- cache TTL: 7 days
- request rate limit: 5 req/s

Optional ligand lists:
- ``--trivial-ligands-file``: line/CSV-like list of comp_ids to classify as trivial
- ``--keep-ligands-file``: line/CSV-like list of comp_ids to always keep as non-trivial
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import random
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import requests  # type: ignore[import-untyped]

LOGGER = logging.getLogger("pdb_query")

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_NONPOLYMER_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{}/{}"

DEFAULT_CACHE_TTL_SECONDS = 7 * 24 * 60 * 60
DEFAULT_RATE_LIMIT_RPS = 5.0
DEFAULT_RESOLUTION_MAX = 3.2
DEFAULT_TIMEOUT_SECONDS = 30.0
DEFAULT_MAX_RETRIES = 4
DEFAULT_SPECIES = "Homo sapiens"
DEFAULT_TAXONOMY_ID = 9606
DEFAULT_ENTITY_TYPE = "Protein"
DEFAULT_EXPERIMENTAL_METHODS: tuple[str, ...] = ()
DEFAULT_METHOD_FILTER = "any"
DEFAULT_DEDUPE_SEQUENCE_IDENTITY: int | None = 90
SEQUENCE_IDENTITY_LEVELS = (100, 95, 90, 70, 50, 30)
QUALITY_PRESETS = ("any", "publication")
LIGAND_FILTER_MODES = ("strict", "relaxed")
DEFAULT_CACHE_DIR = Path(__file__).resolve().parents[1] / "chemdb" / ".cache" / "pdb_query"


@dataclass(frozen=True)
class SpeciesAlias:
    """Common organism shortcut resolved to RCSB source organism filters."""

    canonical: str
    taxonomy_id: int
    aliases: tuple[str, ...]


SPECIES_ALIASES: tuple[SpeciesAlias, ...] = (
    SpeciesAlias("Homo sapiens", 9606, ("human", "homo sapiens", "h sapiens")),
    SpeciesAlias("Mus musculus", 10090, ("mouse", "mice", "mus musculus", "m musculus")),
    SpeciesAlias("Rattus norvegicus", 10116, ("rat", "rattus norvegicus", "r norvegicus")),
    SpeciesAlias("Danio rerio", 7955, ("zebrafish", "danio rerio", "d rerio")),
    SpeciesAlias(
        "Drosophila melanogaster",
        7227,
        ("fly", "fruit fly", "drosophila", "drosophila melanogaster", "d melanogaster"),
    ),
    SpeciesAlias(
        "Caenorhabditis elegans",
        6239,
        ("worm", "c elegans", "celegans", "caenorhabditis elegans"),
    ),
    SpeciesAlias(
        "Saccharomyces cerevisiae",
        4932,
        ("yeast", "budding yeast", "s cerevisiae", "saccharomyces cerevisiae"),
    ),
    SpeciesAlias(
        "Escherichia coli",
        562,
        ("ecoli", "e coli", "escherichia coli", "e. coli"),
    ),
    SpeciesAlias(
        "Arabidopsis thaliana",
        3702,
        ("arabidopsis", "arabidopsis thaliana", "a thaliana"),
    ),
    SpeciesAlias("Gallus gallus", 9031, ("chicken", "gallus gallus", "g gallus")),
    SpeciesAlias("Bos taurus", 9913, ("cow", "cattle", "bovine", "bos taurus")),
    SpeciesAlias("Sus scrofa", 9823, ("pig", "porcine", "sus scrofa")),
    SpeciesAlias("Canis lupus familiaris", 9615, ("dog", "canine", "canis familiaris")),
    SpeciesAlias("Macaca mulatta", 9544, ("rhesus", "rhesus macaque", "macaca mulatta")),
    SpeciesAlias("Cricetulus griseus", 10029, ("hamster", "cho", "chinese hamster")),
    SpeciesAlias("Oryctolagus cuniculus", 9986, ("rabbit", "oryctolagus cuniculus")),
    SpeciesAlias(
        "Severe acute respiratory syndrome coronavirus 2",
        2697049,
        ("sars-cov-2", "sars cov 2", "covid", "covid-19"),
    ),
    SpeciesAlias(
        "Mycobacterium tuberculosis",
        1773,
        ("tb", "tuberculosis", "m tuberculosis", "mycobacterium tuberculosis"),
    ),
    SpeciesAlias(
        "Plasmodium falciparum",
        5833,
        ("malaria", "p falciparum", "plasmodium falciparum"),
    ),
    SpeciesAlias(
        "Candida albicans",
        5476,
        ("candida", "c albicans", "candida albicans"),
    ),
    SpeciesAlias(
        "Staphylococcus aureus",
        1280,
        ("staph", "staph aureus", "s aureus", "staphylococcus aureus"),
    ),
)


def _species_alias_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


SPECIES_ALIAS_BY_KEY: dict[str, SpeciesAlias] = {
    _species_alias_key(alias): item
    for item in SPECIES_ALIASES
    for alias in (item.canonical, *item.aliases)
}

TRIVIAL_NONPOLYMERS = {
    "ACT",
    "BR",
    "CA",
    "CL",
    "CO",
    "CU",
    "DOD",
    "EDO",
    "FE",
    "FMT",
    "GOL",
    "HOH",
    "IOD",
    "K",
    "MG",
    "MN",
    "MPD",
    "NA",
    "NI",
    "PEG",
    "PO4",
    "SO4",
    "ZN",
}

COMMON_CRYSTALLIZATION_ADDITIVES = {
    "BME",
    "CAC",
    "CIT",
    "DMS",
    "DTT",
    "HEP",
    "MES",
    "MRD",
    "PGE",
    "PG4",
    "TRS",
}

GLYCAN_NONPOLYMERS = {
    "A2G",
    "A2M",
    "BMA",
    "FUC",
    "GAL",
    "GLC",
    "MAN",
    "NAG",
    "NDG",
    "NNN",
    "SIA",
    "SIR",
}

# Hard-include LSD-bound HTR2B structures for safety panel / validation.
FORCED_PDBS_BY_GENE: dict[str, list[str]] = {"HTR2B": ["5TVN", "7SRS"]}


@dataclass
class RankedEntry:
    """Ranked entry payload used for programmatic usage and CSV output."""

    gene: str
    pdb_id: str
    resolution: Optional[float]
    method: str
    nonpolymer_comp_ids: list[str] = field(default_factory=list)
    nontrivial_comp_ids: list[str] = field(default_factory=list)
    rank_score: float = 0.0
    rank_reason: dict[str, Any] = field(default_factory=dict)
    forced_included: bool = False
    valid_entry: bool = True


@dataclass
class SelectionResult:
    """Per-gene selection result with summary statistics."""

    entries: list[RankedEntry]
    stats: dict[str, int]


@dataclass(frozen=True)
class QueryFilters:
    """Optional RCSB/Atlas filters shared by gene, UniProt, and text queries."""

    ligand_comp_ids: tuple[str, ...] = ()
    dedupe_sequence_identity: int | None = DEFAULT_DEDUPE_SEQUENCE_IDENTITY
    quality: str = "any"
    ligand_filter_mode: str = "strict"


class PDBQueryClient:
    """HTTP client with caching, polite rate limiting, and retry/backoff."""

    def __init__(
        self,
        *,
        session: Optional[requests.Session] = None,
        cache_dir: Optional[Path] = None,
        use_cache: bool = True,
        cache_ttl_seconds: int = DEFAULT_CACHE_TTL_SECONDS,
        rate_limit_rps: float = DEFAULT_RATE_LIMIT_RPS,
        max_retries: int = DEFAULT_MAX_RETRIES,
        timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    ) -> None:
        self.session = session or requests.Session()
        self.cache_dir = (cache_dir or DEFAULT_CACHE_DIR).expanduser()
        self.use_cache = use_cache
        self.cache_ttl_seconds = cache_ttl_seconds
        self.rate_limit_rps = rate_limit_rps
        self.max_retries = max_retries
        self.timeout_seconds = timeout_seconds
        self._last_request_ts = 0.0
        self._rng = random.Random()
        self._logger = logging.getLogger("pdb_query.http")

    def search(self, payload: dict[str, Any]) -> Optional[dict[str, Any]]:
        return self._request_json(
            method="POST",
            url=RCSB_SEARCH_URL,
            payload=payload,
            timeout=max(self.timeout_seconds, 60.0),
        )

    def fetch_entry(self, entry_id: str) -> Optional[dict[str, Any]]:
        return self._request_json(method="GET", url=RCSB_ENTRY_URL.format(entry_id))

    def fetch_nonpolymer_entity(
        self, entry_id: str, entity_id: str
    ) -> Optional[dict[str, Any]]:
        return self._request_json(
            method="GET",
            url=RCSB_NONPOLYMER_ENTITY_URL.format(entry_id, entity_id),
        )

    def _request_json(
        self,
        *,
        method: str,
        url: str,
        payload: Optional[dict[str, Any]] = None,
        timeout: Optional[float] = None,
    ) -> Optional[dict[str, Any]]:
        method = method.upper()
        use_cache = method == "GET" and self.use_cache
        if use_cache:
            cached = self._load_cache(url)
            if cached is not None:
                return cached

        transient_statuses = {429, 500, 502, 503, 504}
        request_timeout = timeout if timeout is not None else self.timeout_seconds

        for attempt in range(self.max_retries + 1):
            self._wait_for_rate_limit()
            try:
                if method == "POST":
                    response = self.session.post(
                        url, json=payload, timeout=request_timeout
                    )
                elif method == "GET":
                    response = self.session.get(url, timeout=request_timeout)
                else:
                    raise ValueError(f"Unsupported HTTP method: {method}")
            except requests.RequestException as exc:
                if attempt >= self.max_retries:
                    self._logger.warning(
                        "HTTP %s failed for %s after retries: %s", method, url, exc
                    )
                    return None
                self._sleep_backoff(attempt, url, status_code=None)
                continue

            if response.status_code == 200:
                try:
                    parsed = response.json()
                except ValueError:
                    self._logger.warning("Invalid JSON from %s", url)
                    return None
                if use_cache:
                    self._write_cache(url, parsed)
                return parsed

            if response.status_code not in transient_statuses:
                self._logger.warning("HTTP %s for %s", response.status_code, url)
                return None

            if attempt >= self.max_retries:
                self._logger.warning(
                    "HTTP %s for %s after retries exhausted",
                    response.status_code,
                    url,
                )
                return None

            self._sleep_backoff(attempt, url, status_code=response.status_code)
        return None

    def _cache_path(self, url: str) -> Path:
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def _load_cache(self, url: str) -> Optional[dict[str, Any]]:
        path = self._cache_path(url)
        if not path.exists():
            return None
        try:
            blob = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        fetched_at = blob.get("fetched_at")
        payload = blob.get("payload")
        if not isinstance(fetched_at, (int, float)) or not isinstance(payload, dict):
            return None
        age = time.time() - float(fetched_at)
        if age > self.cache_ttl_seconds:
            return None
        return payload

    def _write_cache(self, url: str, payload: dict[str, Any]) -> None:
        path = self._cache_path(url)
        path.parent.mkdir(parents=True, exist_ok=True)
        blob = {"fetched_at": time.time(), "payload": payload}
        path.write_text(json.dumps(blob), encoding="utf-8")

    def _wait_for_rate_limit(self) -> None:
        if self.rate_limit_rps <= 0:
            return
        min_interval = 1.0 / self.rate_limit_rps
        elapsed = time.time() - self._last_request_ts
        if elapsed < min_interval:
            time.sleep(min_interval - elapsed)
        self._last_request_ts = time.time()

    def _sleep_backoff(
        self, attempt: int, url: str, status_code: Optional[int]
    ) -> None:
        wait_s = min(8.0, 0.5 * (2**attempt)) + self._rng.uniform(0.0, 0.2)
        if status_code is None:
            self._logger.info("Retrying %s after network error in %.2fs", url, wait_s)
        else:
            self._logger.info(
                "Retrying %s after HTTP %s in %.2fs", url, status_code, wait_s
            )
        time.sleep(wait_s)


_DEFAULT_CLIENT: Optional[PDBQueryClient] = None


def _default_client() -> PDBQueryClient:
    global _DEFAULT_CLIENT
    if _DEFAULT_CLIENT is None:
        _DEFAULT_CLIENT = PDBQueryClient()
    return _DEFAULT_CLIENT


def _upper_list(values: Sequence[Any]) -> list[str]:
    return [str(value).strip().upper() for value in values if str(value).strip()]


# Public API: species alias catalog for qol_cli and downstream callers.
def known_species_aliases() -> list[SpeciesAlias]:
    return sorted(SPECIES_ALIASES, key=lambda item: item.canonical)


def resolve_species_filter(
    species_name: str | None = None,
    taxonomy_id: int | None = None,
) -> tuple[str | None, int | None]:
    """Resolve common organism shortcuts for RCSB source organism filters."""

    raw_species = str(species_name or "").strip()
    if not raw_species and taxonomy_id is None:
        return DEFAULT_SPECIES, DEFAULT_TAXONOMY_ID

    token = raw_species.lower().strip()
    if token in {"any", "all", "none", "unfiltered", "all-organisms", "all organisms"}:
        return None, taxonomy_id

    alias = SPECIES_ALIAS_BY_KEY.get(_species_alias_key(raw_species))
    if alias is not None:
        return alias.canonical, taxonomy_id if taxonomy_id is not None else alias.taxonomy_id

    if raw_species:
        return raw_species, taxonomy_id
    return None, taxonomy_id


def _get_resolution(entry_json: dict[str, Any]) -> Optional[float]:
    """Extract the best (lowest) combined resolution from an entry payload."""
    info = entry_json.get("rcsb_entry_info", {})
    res = info.get("resolution_combined")
    if isinstance(res, list) and res:
        parsed = [float(x) for x in res if isinstance(x, (int, float))]
        return min(parsed) if parsed else None
    if isinstance(res, (int, float)):
        return float(res)
    return None


def _get_experimental_method(entry_json: dict[str, Any]) -> str:
    """Return normalized experimental method (upper-case) from ``exptl``."""
    expt = entry_json.get("exptl", [])
    if isinstance(expt, list) and expt:
        method = expt[0].get("method", "")
        return str(method or "").upper()
    return ""


def _first_number_by_key(payload: Any, key_names: set[str]) -> float | None:
    if isinstance(payload, Mapping):
        for key, value in payload.items():
            if key in key_names or str(key).lower() in key_names:
                try:
                    return float(value)
                except (TypeError, ValueError):
                    pass
            found = _first_number_by_key(value, key_names)
            if found is not None:
                return found
    elif isinstance(payload, list):
        for item in payload:
            found = _first_number_by_key(item, key_names)
            if found is not None:
                return found
    return None


def _quality_metrics(entry_json: dict[str, Any]) -> dict[str, float | None]:
    return {
        "r_free": _first_number_by_key(entry_json, {"ls_R_factor_R_free", "ls_r_factor_r_free"}),
        "r_work": _first_number_by_key(entry_json, {"ls_R_factor_R_work", "ls_r_factor_r_work"}),
        "clashscore": _first_number_by_key(entry_json, {"clashscore"}),
        "q_score": _first_number_by_key(entry_json, {"Q_score", "q_score", "q_score_combined"}),
    }


def _quality_allowed(
    entry_json: dict[str, Any],
    quality: str | None,
    metrics: dict[str, float | None] | None = None,
) -> bool:
    if _normalize_quality_preset(quality) == "any":
        return True
    metrics = metrics if metrics is not None else _quality_metrics(entry_json)
    checks = (
        ("r_free", metrics["r_free"], 0.30, "max"),
        ("r_work", metrics["r_work"], 0.27, "max"),
        ("clashscore", metrics["clashscore"], 30.0, "max"),
        ("q_score", metrics["q_score"], 0.40, "min"),
    )
    for _name, value, threshold, mode in checks:
        if value is None:
            continue
        if mode == "max" and value > threshold:
            return False
        if mode == "min" and value < threshold:
            return False
    return True


def _normalize_method_filters(methods: str | Sequence[str] | None) -> tuple[str, ...]:
    if methods is None:
        return DEFAULT_EXPERIMENTAL_METHODS
    if isinstance(methods, str):
        raw_methods = re.split(r"[,;|]+", methods)
    else:
        raw_methods = list(methods)
    normalized: list[str] = []
    aliases = {
        "": "",
        "ANY": "",
        "ALL": "",
        "NONE": "",
        "X-RAY": "X-RAY",
        "XRAY": "X-RAY",
        "X RAY": "X-RAY",
        "X-RAY DIFFRACTION": "X-RAY",
        "EM": "ELECTRON MICROSCOPY",
        "CRYOEM": "ELECTRON MICROSCOPY",
        "CRYO-EM": "ELECTRON MICROSCOPY",
        "ELECTRON MICROSCOPY": "ELECTRON MICROSCOPY",
        "NMR": "NMR",
        "SOLUTION NMR": "NMR",
    }
    for method in raw_methods:
        token = str(method).strip().upper()
        mapped = aliases.get(token, token)
        if mapped and mapped not in normalized:
            normalized.append(mapped)
    return tuple(normalized)


def _method_allowed(method: str, methods: str | Sequence[str] | None) -> bool:
    filters = _normalize_method_filters(methods)
    return not filters or any(token in method.upper() for token in filters)


def _normalize_comp_ids(values: str | Sequence[str] | None) -> tuple[str, ...]:
    if values is None:
        return ()
    if isinstance(values, str):
        raw_values = re.split(r"[,;|\s]+", values)
    else:
        raw_values = [
            token
            for value in values
            for token in re.split(r"[,;|\s]+", str(value).strip())
        ]
    normalized: list[str] = []
    for raw_value in raw_values:
        token = str(raw_value).strip().upper()
        if token and token not in normalized:
            normalized.append(token)
    return tuple(normalized)


def _normalize_quality_preset(value: str | None) -> str:
    preset = str(value or "any").strip().lower()
    if preset in {"", "none", "off", "all"}:
        return "any"
    if preset not in QUALITY_PRESETS:
        allowed = ", ".join(QUALITY_PRESETS)
        raise ValueError(f"unknown quality preset '{value}'. Known presets: {allowed}")
    return preset


def _normalize_ligand_filter_mode(value: str | None) -> str:
    token = str(value or "strict").strip().lower().replace("_", "-")
    if token in {"", "default", "auto"}:
        return "strict"
    if token not in LIGAND_FILTER_MODES:
        allowed = ", ".join(LIGAND_FILTER_MODES)
        raise ValueError(f"unknown ligand filter mode '{value}'. Known modes: {allowed}")
    return token


def _normalize_sequence_identity(value: int | str | None) -> int | None:
    if value is None:
        return None
    token = str(value).strip().lower().replace("%", "")
    if token in {"", "none", "off", "any", "no"}:
        return None
    cutoff = int(float(token))
    if cutoff not in SEQUENCE_IDENTITY_LEVELS:
        allowed = ", ".join(str(level) for level in SEQUENCE_IDENTITY_LEVELS)
        raise ValueError(f"--dedupe sequence identity must be one of: {allowed}")
    return cutoff


def _query_filters(
    *,
    ligand_comp_ids: str | Sequence[str] | None = None,
    dedupe_sequence_identity: int | str | None = None,
    quality: str | None = None,
    ligand_filter_mode: str | None = None,
) -> QueryFilters:
    return QueryFilters(
        ligand_comp_ids=_normalize_comp_ids(ligand_comp_ids),
        dedupe_sequence_identity=_normalize_sequence_identity(dedupe_sequence_identity),
        quality=_normalize_quality_preset(quality),
        ligand_filter_mode=_normalize_ligand_filter_mode(ligand_filter_mode),
    )


def normalize_comp_ids(values: str | Sequence[str] | None) -> tuple[str, ...]:
    """Normalize one or more PDB chemical component IDs."""

    return _normalize_comp_ids(values)


def normalize_quality_preset(value: str | None) -> str:
    """Normalize an Atlas target-query quality preset."""

    return _normalize_quality_preset(value)


def parse_dedupe_option(
    dedupe: Sequence[str] | None,
    dedupe_sequence_identity: int | str | None = None,
    *,
    no_dedupe: bool = False,
) -> int | None:
    """Parse CLI de-duplication options such as ``--dedupe sequence-identity 90``."""

    if no_dedupe:
        return None
    if dedupe_sequence_identity is not None:
        return _normalize_sequence_identity(dedupe_sequence_identity)
    if not dedupe:
        return DEFAULT_DEDUPE_SEQUENCE_IDENTITY
    if len(dedupe) != 2:
        raise ValueError("--dedupe expects: sequence-identity <cutoff>")
    mode = str(dedupe[0]).strip().lower().replace("_", "-")
    if mode not in {"sequence-identity", "sequence", "seqid", "seq-id"}:
        raise ValueError("--dedupe currently supports only `sequence-identity <cutoff>`")
    return _normalize_sequence_identity(dedupe[1])


def _entry_matches_species(
    entry_json: dict[str, Any],
    *,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
) -> bool:
    species_token = str(species_name or "").upper()
    taxonomy_token = str(taxonomy_id or "")
    serialized = json.dumps(entry_json).upper()
    if (species_token and species_token in serialized) or (
        taxonomy_token and taxonomy_token in serialized
    ):
        return True

    # Entry-level payloads may not include source-organism blocks; in that case,
    # trust the upstream search query constraints that already enforce organism hits.
    has_any_organism_block = (
        "RCSB_ENTITY_SOURCE_ORGANISM" in serialized
        or "ENTITY_SRC_GEN" in serialized
        or "ENTITY_SRC_NAT" in serialized
    )
    return not has_any_organism_block


def _method_score(method: str) -> float:
    if "X-RAY" in method:
        return 40.0
    if "ELECTRON MICROSCOPY" in method:
        return 25.0
    if "NMR" in method:
        return 10.0
    return 5.0


def _score_entry(
    entry_json: dict[str, Any], nontrivial_count: int, resolution: Optional[float]
) -> tuple[float, dict[str, Any]]:
    method = _get_experimental_method(entry_json)
    info = entry_json.get("rcsb_entry_info", {})
    polymer_count = int(info.get("polymer_entity_count", 0) or 0)
    protein_count = int(info.get("polymer_entity_count_protein", 0) or 0)
    nonprotein_polymer_count = max(0, polymer_count - protein_count)
    assembly_count = int(info.get("assembly_count", 1) or 1)

    method_component = _method_score(method)
    resolution_component = 0.0
    if resolution is not None:
        resolution_component = max(0.0, 30.0 - (resolution * 8.0))

    ligand_component = 15.0 + min(nontrivial_count, 4) * 2.0 if nontrivial_count else -20.0
    nonprotein_penalty = float(nonprotein_polymer_count * 4)
    size_penalty = max(0.0, (assembly_count - 1) * 1.5) + max(
        0.0, (polymer_count - 4) * 1.0
    )
    total = (
        method_component
        + resolution_component
        + ligand_component
        - nonprotein_penalty
        - size_penalty
    )

    reason = {
        "method": method,
        "method_component": round(method_component, 3),
        "resolution": resolution,
        "resolution_component": round(resolution_component, 3),
        "nontrivial_ligand_count": nontrivial_count,
        "ligand_component": round(ligand_component, 3),
        "nonprotein_polymer_count": nonprotein_polymer_count,
        "nonprotein_penalty": round(nonprotein_penalty, 3),
        "assembly_count": assembly_count,
        "polymer_entity_count": polymer_count,
        "size_penalty": round(size_penalty, 3),
        "co_crystallized_ligand": nontrivial_count > 0,
        "score": round(total, 3),
    }
    return total, reason


def _read_comp_id_file(path: Optional[str]) -> set[str]:
    if not path:
        return set()
    comp_ids: set[str] = set()
    text = Path(path).read_text(encoding="utf-8")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [token.strip() for token in re.split(r"[,\s;|]+", line) if token.strip()]
        comp_ids.update(value.upper() for value in parts)
    return comp_ids


def _classify_nonpolymer_comp_ids(
    comp_ids: Sequence[str],
    *,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
    ligand_filter_mode: str = "strict",
) -> tuple[list[str], list[str]]:
    normalized = _upper_list(comp_ids)
    blocked = set(TRIVIAL_NONPOLYMERS)
    blocked.update(COMMON_CRYSTALLIZATION_ADDITIVES)
    if ligand_filter_mode == "strict":
        blocked.update(GLYCAN_NONPOLYMERS)
    if trivial_ligands:
        blocked.update(_upper_list(sorted(trivial_ligands)))
    keep = set(_upper_list(sorted(keep_ligands or set())))
    nontrivial = [comp for comp in normalized if comp in keep or comp not in blocked]
    return normalized, nontrivial


def _list_nonpolymer_comp_ids(
    entry_id: str, entry_json: dict[str, Any], client: Optional[PDBQueryClient] = None
) -> list[str]:
    """Return nonpolymer comp_ids from ``/nonpolymer_entity`` endpoints."""
    active_client = client or _default_client()
    container = entry_json.get("rcsb_entry_container_identifiers", {})
    nonpolymer_entity_ids = container.get("non_polymer_entity_ids")
    if not isinstance(nonpolymer_entity_ids, list) or not nonpolymer_entity_ids:
        info = entry_json.get("rcsb_entry_info", {})
        nonpolymer_entity_ids = info.get("nonpolymer_entity_ids")
    if not isinstance(nonpolymer_entity_ids, list) or not nonpolymer_entity_ids:
        return []

    comp_ids: list[str] = []
    for raw_ent_id in nonpolymer_entity_ids:
        entity = active_client.fetch_nonpolymer_entity(entry_id, str(raw_ent_id))
        if not entity:
            continue
        comp_id = (
            entity.get("rcsb_nonpolymer_entity_container_identifiers", {}).get(
                "nonpolymer_comp_id"
            )
            or entity.get("rcsb_nonpolymer_entity_container_identifiers", {}).get(
                "chem_ref_def_id"
            )
            or entity.get("pdbx_entity_nonpoly", {}).get("comp_id")
            or entity.get("chem_comp", {}).get("id")
        )
        if isinstance(comp_id, str) and comp_id:
            comp_ids.append(comp_id.upper())
    return comp_ids


def _fetch_entry(
    entry_id: str, client: Optional[PDBQueryClient] = None
) -> Optional[dict[str, Any]]:
    """Fetch RCSB entry metadata for a PDB ID."""
    active_client = client or _default_client()
    return active_client.fetch_entry(entry_id.upper())


def _search_rcsb(
    payload: dict[str, Any],
    *,
    client: Optional[PDBQueryClient] = None,
    max_candidates: int = 200,
    dedupe_sequence_identity: int | None = None,
) -> list[str]:
    active_client = client or _default_client()
    query_payload = _with_sequence_identity_grouping(payload, dedupe_sequence_identity)
    response = active_client.search(query_payload)
    if not response:
        return []
    result_set = response.get("result_set", [])
    if not isinstance(result_set, list):
        return []
    candidate_ids = _upper_list(
        [
            _entry_id_from_search_identifier(str(item.get("identifier") or ""))
            for item in result_set
            if isinstance(item, dict)
        ]
    )
    deduped = list(dict.fromkeys(candidate_ids))
    return deduped[:max_candidates]


def _entry_id_from_search_identifier(identifier: str) -> str:
    token = identifier.strip().upper()
    if not token:
        return ""
    for separator in ("_", ".", "-"):
        if separator in token:
            return token.split(separator, 1)[0]
    return token


def _with_sequence_identity_grouping(
    payload: dict[str, Any],
    cutoff: int | None,
) -> dict[str, Any]:
    if cutoff is None:
        return payload
    grouped = json.loads(json.dumps(payload))
    grouped["return_type"] = "polymer_entity"
    request_options = grouped.setdefault("request_options", {})
    request_options["group_by"] = {
        "aggregation_method": "sequence_identity",
        "similarity_cutoff": cutoff,
    }
    request_options["group_by_return_type"] = "representatives"
    return grouped


def _organism_and_polymer_nodes(
    *,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
) -> list[dict[str, Any]]:
    nodes: list[dict[str, Any]] = []
    organism_nodes = []
    if species_name:
        organism_nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                    "operator": "exact_match",
                    "value": species_name,
                },
            }
        )
    if taxonomy_id is not None:
        organism_nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                    "operator": "exact_match",
                    "value": str(taxonomy_id),
                },
            }
        )
    if len(organism_nodes) == 1:
        nodes.append(organism_nodes[0])
    elif organism_nodes:
        nodes.append(
            {
                "type": "group",
                "logical_operator": "or",
                "nodes": organism_nodes,
            }
        )

    normalized_entity_type = (entity_type or "").strip()
    if normalized_entity_type.lower() == "protein":
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                    "operator": "greater_or_equal",
                    "value": 1,
                },
            }
        )
    elif normalized_entity_type:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                    "operator": "exact_match",
                    "value": normalized_entity_type,
                },
            }
        )
    return nodes


def _ligand_nodes(ligand_comp_ids: Sequence[str]) -> list[dict[str, Any]]:
    comp_ids = _normalize_comp_ids(ligand_comp_ids)
    if not comp_ids:
        return []
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                "operator": "exact_match",
                "value": comp_id,
            },
        }
        for comp_id in comp_ids
    ]
    if len(nodes) == 1:
        return nodes
    return [{"type": "group", "logical_operator": "or", "nodes": nodes}]


def _rcsb_entry_search_payload(nodes: list[dict[str, Any]]) -> dict[str, Any]:
    """Shared AND-group query + ``return_all_hits`` for RCSB search ``/v2/query`` (entry)."""
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "request_options": {"return_all_hits": True},
        "return_type": "entry",
    }


def _build_gene_query(
    gene_name: str,
    *,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    ligand_comp_ids: Sequence[str] = (),
) -> dict[str, Any]:
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                "operator": "exact_match",
                "value": gene_name.strip().upper(),
            },
        }
    ]
    nodes.extend(
        _organism_and_polymer_nodes(
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=entity_type,
        )
    )
    nodes.extend(_ligand_nodes(ligand_comp_ids))
    return _rcsb_entry_search_payload(nodes)


def _build_uniprot_query(
    uniprot_id: str,
    attribute: str,
    *,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    ligand_comp_ids: Sequence[str] = (),
) -> dict[str, Any]:
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": attribute,
                "operator": "exact_match",
                "value": uniprot_id,
            },
        }
    ]
    nodes.extend(
        _organism_and_polymer_nodes(
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=entity_type,
        )
    )
    nodes.extend(_ligand_nodes(ligand_comp_ids))
    return _rcsb_entry_search_payload(nodes)


def _build_text_query(
    text: str,
    *,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    ligand_comp_ids: Sequence[str] = (),
) -> dict[str, Any]:
    nodes = [
        {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": text.strip(),
            },
        }
    ]
    nodes.extend(
        _organism_and_polymer_nodes(
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=entity_type,
        )
    )
    nodes.extend(_ligand_nodes(ligand_comp_ids))
    return _rcsb_entry_search_payload(nodes)


def _rank_entries(entries: list[RankedEntry], *, max_return: int) -> list[RankedEntry]:
    ordered = sorted(entries, key=lambda row: (-row.rank_score, row.pdb_id))
    return ordered[:max_return]


def _evaluate_entry(
    *,
    gene_name: str,
    entry_id: str,
    resolution_max: float,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    experimental_methods: str | Sequence[str] | None = DEFAULT_EXPERIMENTAL_METHODS,
    require_ligand: bool = True,
    filters: QueryFilters = QueryFilters(),
    client: Optional[PDBQueryClient] = None,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
    ligand_filter_mode: str = "strict",
) -> Optional[RankedEntry]:
    entry_json = _fetch_entry(entry_id, client=client)
    if not entry_json or not _entry_matches_species(
        entry_json,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
    ):
        return None

    resolution = _get_resolution(entry_json)
    if resolution is not None and resolution > resolution_max:
        return None

    method = _get_experimental_method(entry_json)
    if not _method_allowed(method, experimental_methods):
        return None

    all_comp_ids, nontrivial_comp_ids = _classify_nonpolymer_comp_ids(
        _list_nonpolymer_comp_ids(entry_id, entry_json, client=client),
        trivial_ligands=trivial_ligands,
        keep_ligands=keep_ligands,
        ligand_filter_mode=ligand_filter_mode,
    )
    if filters.ligand_comp_ids and not set(filters.ligand_comp_ids).intersection(all_comp_ids):
        return None
    if require_ligand and not nontrivial_comp_ids:
        return None
    quality_metrics = _quality_metrics(entry_json)
    if not _quality_allowed(entry_json, filters.quality, quality_metrics):
        return None

    score, reason = _score_entry(entry_json, len(nontrivial_comp_ids), resolution)
    reason["quality"] = filters.quality
    reason["quality_metrics"] = quality_metrics
    return RankedEntry(
        gene=gene_name,
        pdb_id=entry_id.upper(),
        resolution=resolution,
        method=method,
        nonpolymer_comp_ids=all_comp_ids,
        nontrivial_comp_ids=nontrivial_comp_ids,
        rank_score=score,
        rank_reason=reason,
    )


# Public API: legacy gene search entry point; returns ranked PDB IDs only.
def search_rcsb_by_gene(
    gene_name: str,
    max_candidates: int = 200,
    max_return: int = 10,
    *,
    resolution_max: float = DEFAULT_RESOLUTION_MAX,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    experimental_methods: str | Sequence[str] | None = DEFAULT_EXPERIMENTAL_METHODS,
    require_ligand: bool = True,
    ligand_comp_ids: str | Sequence[str] | None = None,
    dedupe_sequence_identity: int | str | None = DEFAULT_DEDUPE_SEQUENCE_IDENTITY,
    quality: str | None = "any",
    ligand_filter_mode: str | None = None,
    client: Optional[PDBQueryClient] = None,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
) -> list[str]:
    """Search by gene symbol and return ranked PDB IDs."""
    filters = _query_filters(
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        quality=quality,
        ligand_filter_mode=ligand_filter_mode,
    )
    query = _build_gene_query(
        gene_name,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
        entity_type=entity_type,
        ligand_comp_ids=filters.ligand_comp_ids,
    )
    candidate_ids = _search_rcsb(
        query,
        client=client,
        max_candidates=max_candidates,
        dedupe_sequence_identity=filters.dedupe_sequence_identity,
    )
    ranked: list[RankedEntry] = []
    for entry_id in candidate_ids:
        row = _evaluate_entry(
            gene_name=gene_name,
            entry_id=entry_id,
            resolution_max=resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            experimental_methods=experimental_methods,
            require_ligand=require_ligand,
            filters=filters,
            client=client,
            trivial_ligands=trivial_ligands,
            keep_ligands=keep_ligands,
            ligand_filter_mode=filters.ligand_filter_mode,
        )
        if row is not None:
            ranked.append(row)
    return [row.pdb_id for row in _rank_entries(ranked, max_return=max_return)]


# Public API: legacy UniProt search entry point; returns ranked PDB IDs only.
def search_rcsb_by_uniprot(
    uniprot_id: str,
    max_candidates: int = 200,
    max_return: int = 10,
    *,
    resolution_max: float = DEFAULT_RESOLUTION_MAX,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    experimental_methods: str | Sequence[str] | None = DEFAULT_EXPERIMENTAL_METHODS,
    require_ligand: bool = True,
    ligand_comp_ids: str | Sequence[str] | None = None,
    dedupe_sequence_identity: int | str | None = DEFAULT_DEDUPE_SEQUENCE_IDENTITY,
    quality: str | None = "any",
    ligand_filter_mode: str | None = None,
    client: Optional[PDBQueryClient] = None,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
) -> list[str]:
    """Search by UniProt accession and return ranked PDB IDs."""
    filters = _query_filters(
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        quality=quality,
        ligand_filter_mode=ligand_filter_mode,
    )
    attributes = [
        "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
    ]
    candidate_ids: list[str] = []
    for attribute in attributes:
        query = _build_uniprot_query(
            uniprot_id,
            attribute=attribute,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=entity_type,
            ligand_comp_ids=filters.ligand_comp_ids,
        )
        candidate_ids.extend(
            _search_rcsb(
                query,
                client=client,
                max_candidates=max_candidates,
                dedupe_sequence_identity=filters.dedupe_sequence_identity,
            )
        )
    deduped_ids = list(dict.fromkeys(_upper_list(candidate_ids)))[:max_candidates]

    ranked: list[RankedEntry] = []
    for entry_id in deduped_ids:
        row = _evaluate_entry(
            gene_name=uniprot_id,
            entry_id=entry_id,
            resolution_max=resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            experimental_methods=experimental_methods,
            require_ligand=require_ligand,
            filters=filters,
            client=client,
            trivial_ligands=trivial_ligands,
            keep_ligands=keep_ligands,
            ligand_filter_mode=filters.ligand_filter_mode,
        )
        if row is not None:
            ranked.append(row)
    return [row.pdb_id for row in _rank_entries(ranked, max_return=max_return)]


def get_ranked_entries_for_gene(
    gene_name: str,
    *,
    max_candidates: int = 200,
    max_return: int = 10,
    resolution_max: float = DEFAULT_RESOLUTION_MAX,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    experimental_methods: str | Sequence[str] | None = DEFAULT_EXPERIMENTAL_METHODS,
    require_ligand: bool = True,
    ligand_comp_ids: str | Sequence[str] | None = None,
    dedupe_sequence_identity: int | str | None = DEFAULT_DEDUPE_SEQUENCE_IDENTITY,
    quality: str | None = "any",
    ligand_filter_mode: str | None = None,
    gene_to_uniprots: Optional[Mapping[str, Sequence[str]]] = None,
    forced_pdbs_by_gene: Optional[Mapping[str, Sequence[str]]] = None,
    include_forced: bool = True,
    client: Optional[PDBQueryClient] = None,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
) -> SelectionResult:
    """Return ranked entry records (with scoring explanations) for a gene."""
    filters = _query_filters(
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        quality=quality,
        ligand_filter_mode=ligand_filter_mode,
    )
    gene_upper = gene_name.upper()
    stats = {
        "candidates_before_filters": 0,
        "selected_after_filters": 0,
        "selected_final": 0,
    }

    gene_query = _build_gene_query(
        gene_name,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
        entity_type=entity_type,
        ligand_comp_ids=filters.ligand_comp_ids,
    )
    candidate_ids = _search_rcsb(
        gene_query,
        client=client,
        max_candidates=max_candidates,
        dedupe_sequence_identity=filters.dedupe_sequence_identity,
    )
    uniprots = _upper_list((gene_to_uniprots or {}).get(gene_upper, []))
    for uniprot_id in uniprots:
        for attribute in (
            "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        ):
            query = _build_uniprot_query(
                uniprot_id,
                attribute=attribute,
                species_name=species_name,
                taxonomy_id=taxonomy_id,
                entity_type=entity_type,
                ligand_comp_ids=filters.ligand_comp_ids,
            )
            candidate_ids.extend(
                _search_rcsb(
                    query,
                    client=client,
                    max_candidates=max_candidates,
                    dedupe_sequence_identity=filters.dedupe_sequence_identity,
                )
            )
    deduped_candidates = list(dict.fromkeys(_upper_list(candidate_ids)))[:max_candidates]
    stats["candidates_before_filters"] = len(deduped_candidates)

    evaluated: list[RankedEntry] = []
    for entry_id in deduped_candidates:
        row = _evaluate_entry(
            gene_name=gene_name,
            entry_id=entry_id,
            resolution_max=resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            experimental_methods=experimental_methods,
            require_ligand=require_ligand,
            filters=filters,
            client=client,
            trivial_ligands=trivial_ligands,
            keep_ligands=keep_ligands,
            ligand_filter_mode=filters.ligand_filter_mode,
        )
        if row is not None:
            evaluated.append(row)
    selected = _rank_entries(evaluated, max_return=max_return)
    stats["selected_after_filters"] = len(selected)

    forced_map: Mapping[str, Sequence[str]] = forced_pdbs_by_gene or FORCED_PDBS_BY_GENE
    forced_ids = _upper_list(forced_map.get(gene_upper, []))
    by_id: dict[str, RankedEntry] = {row.pdb_id: row for row in selected}

    if include_forced:
        for forced_id in forced_ids:
            if forced_id in by_id:
                by_id[forced_id].forced_included = True
                continue

            forced_entry_json = _fetch_entry(forced_id, client=client)
            if not forced_entry_json:
                LOGGER.warning(
                    "Forced PDB %s for gene %s not found in RCSB entry API",
                    forced_id,
                    gene_name,
                )
                by_id[forced_id] = RankedEntry(
                    gene=gene_name,
                    pdb_id=forced_id,
                    resolution=None,
                    method="",
                    nonpolymer_comp_ids=[],
                    nontrivial_comp_ids=[],
                    rank_score=float("-inf"),
                    rank_reason={
                        "forced": True,
                        "valid_entry": False,
                        "reason": "entry_not_found",
                    },
                    forced_included=True,
                    valid_entry=False,
                )
                continue

            method = _get_experimental_method(forced_entry_json)
            resolution = _get_resolution(forced_entry_json)
            all_comp_ids, nontrivial_comp_ids = _classify_nonpolymer_comp_ids(
                _list_nonpolymer_comp_ids(forced_id, forced_entry_json, client=client),
                trivial_ligands=trivial_ligands,
                keep_ligands=keep_ligands,
                ligand_filter_mode=filters.ligand_filter_mode,
            )
            if filters.ligand_comp_ids and not set(filters.ligand_comp_ids).intersection(all_comp_ids):
                continue
            quality_metrics = _quality_metrics(forced_entry_json)
            if not _quality_allowed(forced_entry_json, filters.quality, quality_metrics):
                continue
            score, reason = _score_entry(
                forced_entry_json,
                nontrivial_count=len(nontrivial_comp_ids),
                resolution=resolution,
            )
            reason["forced"] = True
            reason["quality"] = filters.quality
            reason["quality_metrics"] = quality_metrics
            by_id[forced_id] = RankedEntry(
                gene=gene_name,
                pdb_id=forced_id,
                resolution=resolution,
                method=method,
                nonpolymer_comp_ids=all_comp_ids,
                nontrivial_comp_ids=nontrivial_comp_ids,
                rank_score=score,
                rank_reason=reason,
                forced_included=True,
                valid_entry=True,
            )

    merged = sorted(by_id.values(), key=lambda row: (-row.rank_score, row.pdb_id))
    stats["selected_final"] = len(merged)
    return SelectionResult(entries=merged, stats=stats)


def get_ranked_entries_for_uniprot(
    uniprot_id: str,
    *,
    max_candidates: int = 200,
    max_return: int = 10,
    resolution_max: float = DEFAULT_RESOLUTION_MAX,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    experimental_methods: str | Sequence[str] | None = DEFAULT_EXPERIMENTAL_METHODS,
    require_ligand: bool = True,
    ligand_comp_ids: str | Sequence[str] | None = None,
    dedupe_sequence_identity: int | str | None = DEFAULT_DEDUPE_SEQUENCE_IDENTITY,
    quality: str | None = "any",
    ligand_filter_mode: str | None = None,
    client: Optional[PDBQueryClient] = None,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
) -> SelectionResult:
    filters = _query_filters(
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        quality=quality,
        ligand_filter_mode=ligand_filter_mode,
    )
    attributes = [
        "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
    ]
    candidate_ids: list[str] = []
    for attribute in attributes:
        query = _build_uniprot_query(
            uniprot_id,
            attribute=attribute,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=entity_type,
            ligand_comp_ids=filters.ligand_comp_ids,
        )
        candidate_ids.extend(
            _search_rcsb(
                query,
                client=client,
                max_candidates=max_candidates,
                dedupe_sequence_identity=filters.dedupe_sequence_identity,
            )
        )
    return _rank_candidate_entries(
        label=uniprot_id,
        candidate_ids=candidate_ids,
        max_candidates=max_candidates,
        max_return=max_return,
        resolution_max=resolution_max,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
        experimental_methods=experimental_methods,
        require_ligand=require_ligand,
        filters=filters,
        client=client,
        trivial_ligands=trivial_ligands,
        keep_ligands=keep_ligands,
    )


# Public API: free-text RCSB query with full ranked SelectionResult payload.
def get_ranked_entries_for_text(
    text: str,
    *,
    max_candidates: int = 200,
    max_return: int = 10,
    resolution_max: float = DEFAULT_RESOLUTION_MAX,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    experimental_methods: str | Sequence[str] | None = DEFAULT_EXPERIMENTAL_METHODS,
    require_ligand: bool = True,
    ligand_comp_ids: str | Sequence[str] | None = None,
    dedupe_sequence_identity: int | str | None = DEFAULT_DEDUPE_SEQUENCE_IDENTITY,
    quality: str | None = "any",
    ligand_filter_mode: str | None = None,
    client: Optional[PDBQueryClient] = None,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
) -> SelectionResult:
    filters = _query_filters(
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        quality=quality,
        ligand_filter_mode=ligand_filter_mode,
    )
    query = _build_text_query(
        text,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
        entity_type=entity_type,
        ligand_comp_ids=filters.ligand_comp_ids,
    )
    candidate_ids = _search_rcsb(
        query,
        client=client,
        max_candidates=max_candidates,
        dedupe_sequence_identity=filters.dedupe_sequence_identity,
    )
    return _rank_candidate_entries(
        label=text,
        candidate_ids=candidate_ids,
        max_candidates=max_candidates,
        max_return=max_return,
        resolution_max=resolution_max,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
        experimental_methods=experimental_methods,
        require_ligand=require_ligand,
        filters=filters,
        client=client,
        trivial_ligands=trivial_ligands,
        keep_ligands=keep_ligands,
    )


def _rank_candidate_entries(
    *,
    label: str,
    candidate_ids: Sequence[str],
    max_candidates: int,
    max_return: int,
    resolution_max: float,
    species_name: str | None,
    taxonomy_id: int | None,
    experimental_methods: str | Sequence[str] | None,
    require_ligand: bool,
    filters: QueryFilters,
    client: Optional[PDBQueryClient],
    trivial_ligands: Optional[set[str]],
    keep_ligands: Optional[set[str]],
) -> SelectionResult:
    deduped_candidates = list(dict.fromkeys(_upper_list(candidate_ids)))[:max_candidates]
    stats = {
        "candidates_before_filters": len(deduped_candidates),
        "selected_after_filters": 0,
        "selected_final": 0,
    }
    evaluated: list[RankedEntry] = []
    for entry_id in deduped_candidates:
        row = _evaluate_entry(
            gene_name=label,
            entry_id=entry_id,
            resolution_max=resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            experimental_methods=experimental_methods,
            require_ligand=require_ligand,
            filters=filters,
            client=client,
            trivial_ligands=trivial_ligands,
            keep_ligands=keep_ligands,
            ligand_filter_mode=filters.ligand_filter_mode,
        )
        if row is not None:
            evaluated.append(row)
    selected = _rank_entries(evaluated, max_return=max_return)
    stats["selected_after_filters"] = len(selected)
    stats["selected_final"] = len(selected)
    return SelectionResult(entries=selected, stats=stats)


# Public API: gene-to-PDB helper used by tests and selection workflows.
def get_pdbs_for_gene(
    gene_name: str,
    max_candidates: int = 200,
    max_return: int = 10,
    *,
    gene_to_uniprots: Optional[Mapping[str, Sequence[str]]] = None,
    include_forced: bool = True,
    forced_pdbs_by_gene: Optional[Mapping[str, Sequence[str]]] = None,
    resolution_max: float = DEFAULT_RESOLUTION_MAX,
    species_name: str | None = DEFAULT_SPECIES,
    taxonomy_id: int | None = DEFAULT_TAXONOMY_ID,
    entity_type: str | None = DEFAULT_ENTITY_TYPE,
    experimental_methods: str | Sequence[str] | None = DEFAULT_EXPERIMENTAL_METHODS,
    require_ligand: bool = True,
    ligand_comp_ids: str | Sequence[str] | None = None,
    dedupe_sequence_identity: int | str | None = DEFAULT_DEDUPE_SEQUENCE_IDENTITY,
    quality: str | None = "any",
    client: Optional[PDBQueryClient] = None,
    trivial_ligands: Optional[set[str]] = None,
    keep_ligands: Optional[set[str]] = None,
) -> list[str]:
    """Return top PDB IDs for a gene, optionally merging forced inclusions."""
    result = get_ranked_entries_for_gene(
        gene_name,
        max_candidates=max_candidates,
        max_return=max_return,
        resolution_max=resolution_max,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
        entity_type=entity_type,
        experimental_methods=experimental_methods,
        require_ligand=require_ligand,
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        quality=quality,
        gene_to_uniprots=gene_to_uniprots,
        forced_pdbs_by_gene=forced_pdbs_by_gene,
        include_forced=include_forced,
        client=client,
        trivial_ligands=trivial_ligands,
        keep_ligands=keep_ligands,
    )
    return [entry.pdb_id for entry in result.entries]


def _read_gene_list(path: str) -> list[str]:
    genes: list[str] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames and "gene" in [name.lower() for name in reader.fieldnames]:
            gene_field = next(
                name for name in reader.fieldnames if name and name.lower() == "gene"
            )
            for row in reader:
                value = (row.get(gene_field) or "").strip()
                if value:
                    genes.append(value)
            return genes

    with Path(path).open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "," in line:
                value = line.split(",", 1)[0].strip()
            else:
                value = line
            if value and value.lower() != "gene":
                genes.append(value)
    return genes


def _read_gene_to_uniprot_map(path: Optional[str]) -> dict[str, list[str]]:
    if not path:
        return {}
    mapping: dict[str, list[str]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        gene_col = next((c for c in reader.fieldnames if c and c.lower() == "gene"), None)
        uniprot_col = next(
            (
                c
                for c in reader.fieldnames
                if c and c.lower() in {"uniprot", "uniprot_id", "uniprots"}
            ),
            None,
        )
        if gene_col is None or uniprot_col is None:
            return {}
        for row in reader:
            gene = (row.get(gene_col) or "").strip().upper()
            raw = (row.get(uniprot_col) or "").strip()
            if not gene or not raw:
                continue
            parts = [token for token in re.split(r"[,\s;|]+", raw) if token]
            mapping.setdefault(gene, [])
            for part in _upper_list(parts):
                if part not in mapping[gene]:
                    mapping[gene].append(part)
    return mapping


def _read_forced_map(path: Optional[str]) -> dict[str, list[str]]:
    if not path:
        return {}
    forced: dict[str, list[str]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        gene_col = next((c for c in reader.fieldnames if c and c.lower() == "gene"), None)
        pdb_col = next(
            (c for c in reader.fieldnames if c and c.lower() in {"pdb_id", "pdb", "entry"}),
            None,
        )
        if gene_col is None or pdb_col is None:
            return {}
        for row in reader:
            gene = (row.get(gene_col) or "").strip().upper()
            pdb = (row.get(pdb_col) or "").strip().upper()
            if not gene or not pdb:
                continue
            forced.setdefault(gene, [])
            if pdb not in forced[gene]:
                forced[gene].append(pdb)
    return forced


def _merge_forced_maps(*maps: Mapping[str, Sequence[str]]) -> dict[str, list[str]]:
    merged: dict[str, list[str]] = {}
    for mapping in maps:
        for gene, ids in mapping.items():
            gene_upper = gene.upper()
            merged.setdefault(gene_upper, [])
            for entry_id in _upper_list(ids):
                if entry_id not in merged[gene_upper]:
                    merged[gene_upper].append(entry_id)
    return merged


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Query RCSB PDB entries for docking panels")
    parser.add_argument("genes", nargs="*", help="Gene symbols to query")
    parser.add_argument("--genes-file", help="CSV/TXT file with genes (column `gene` preferred)")
    parser.add_argument(
        "--uniprot",
        action="append",
        default=[],
        help="Direct UniProt accession to query. Repeat for multiple accessions.",
    )
    parser.add_argument("--gene-uniprot-file", help="CSV with columns gene,uniprot")
    parser.add_argument("--forced-csv", help="CSV with columns gene,pdb_id")
    parser.add_argument("--trivial-ligands-file", help="Override/additional trivial ligand comp_ids")
    parser.add_argument("--keep-ligands-file", help="Allowlist comp_ids always treated non-trivial")
    parser.add_argument("--max-candidates", type=int, default=200)
    parser.add_argument("--max-return", type=int, default=10)
    parser.add_argument("--resolution-max", type=float, default=DEFAULT_RESOLUTION_MAX)
    parser.add_argument(
        "--species",
        "--organism",
        default=None,
        help="Organism name or alias, e.g. human, mouse, rat, zebrafish, or 'any'.",
    )
    parser.add_argument(
        "--taxonomy-id",
        "--taxon",
        dest="taxonomy_id",
        type=int,
        default=None,
        help="NCBI taxonomy id used in the RCSB query; inferred for known aliases.",
    )
    parser.add_argument(
        "--entity-type",
        default=DEFAULT_ENTITY_TYPE,
        help="RCSB polymer entity type; defaults to Protein.",
    )
    parser.add_argument(
        "--methods",
        default=DEFAULT_METHOD_FILTER,
        help="Comma-separated experimental methods; use 'any' to disable method filtering.",
    )
    parser.add_argument(
        "--allow-apo",
        action="store_true",
        help="Keep structures without a non-trivial bound ligand.",
    )
    parser.add_argument(
        "--ligand-filter-mode",
        default="strict",
        choices=LIGAND_FILTER_MODES,
        help="How to treat ambiguous ligand-like components in filter checks.",
    )
    parser.add_argument(
        "--ligand",
        "--contains-ligand",
        action="append",
        default=[],
        dest="ligands",
        help="Require a co-crystallized PDB chemical component ID, e.g. ATP or HEM.",
    )
    parser.add_argument(
        "--dedupe",
        nargs=2,
        metavar=("MODE", "CUTOFF"),
        help="De-duplicate search hits, e.g. `--dedupe sequence-identity 90`.",
    )
    parser.add_argument(
        "--dedupe-sequence-identity",
        type=int,
        choices=SEQUENCE_IDENTITY_LEVELS,
        help="Expert alias for `--dedupe sequence-identity <cutoff>`.",
    )
    parser.add_argument(
        "--no-dedupe",
        action="store_true",
        help="Disable default sequence-identity de-duplication.",
    )
    parser.add_argument(
        "--quality",
        choices=QUALITY_PRESETS,
        default="any",
        help="Target quality preset; `publication` filters bad R-free/R-work/clashscore/Q-score when reported.",
    )
    parser.add_argument("--include-forced", action="store_true", default=True)
    parser.add_argument("--exclude-forced", dest="include_forced", action="store_false")
    parser.add_argument("--out", required=True, help="Output CSV path")
    parser.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR))
    parser.add_argument("--cache-ttl-days", type=float, default=7.0)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--rate-limit", type=float, default=DEFAULT_RATE_LIMIT_RPS)
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT_SECONDS)
    parser.add_argument("--max-retries", type=int, default=DEFAULT_MAX_RETRIES)
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser


def _entries_to_rows(entries: Sequence[RankedEntry]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for entry in entries:
        rows.append(
            {
                "gene": entry.gene,
                "pdb_id": entry.pdb_id,
                "resolution": "" if entry.resolution is None else f"{entry.resolution:.3f}",
                "method": entry.method,
                "nonpolymer_comp_ids": ";".join(entry.nonpolymer_comp_ids),
                "nontrivial_comp_ids": ";".join(entry.nontrivial_comp_ids),
                "rank_score": f"{entry.rank_score:.3f}",
                "rank_reason": json.dumps(entry.rank_reason, sort_keys=True),
                "forced_included": str(entry.forced_included).lower(),
                "valid_entry": str(entry.valid_entry).lower(),
            }
        )
    return rows


def _write_output_csv(path: str, rows: Sequence[dict[str, str]]) -> None:
    fieldnames = [
        "gene",
        "pdb_id",
        "resolution",
        "method",
        "nonpolymer_comp_ids",
        "nontrivial_comp_ids",
        "rank_score",
        "rank_reason",
        "forced_included",
        "valid_entry",
    ]
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    genes = list(args.genes)
    if args.genes_file:
        genes.extend(_read_gene_list(args.genes_file))
    genes = [gene.strip() for gene in genes if gene.strip()]
    genes = list(dict.fromkeys(genes))
    uniprots = list(dict.fromkeys(_upper_list(args.uniprot)))
    if not genes and not uniprots:
        parser.error("Provide gene names, --genes-file, or --uniprot")

    trivial_ligands = _read_comp_id_file(args.trivial_ligands_file)
    keep_ligands = _read_comp_id_file(args.keep_ligands_file)
    gene_to_uniprots = _read_gene_to_uniprot_map(args.gene_uniprot_file)
    forced_map = _merge_forced_maps(FORCED_PDBS_BY_GENE, _read_forced_map(args.forced_csv))
    species_name, taxonomy_id = resolve_species_filter(args.species, args.taxonomy_id)
    try:
        ligand_comp_ids = normalize_comp_ids(args.ligands)
        dedupe_sequence_identity = parse_dedupe_option(
            args.dedupe,
            args.dedupe_sequence_identity,
            no_dedupe=args.no_dedupe,
        )
        quality = normalize_quality_preset(args.quality)
        ligand_filter_mode = _normalize_ligand_filter_mode(args.ligand_filter_mode)
    except ValueError as exc:
        parser.error(str(exc))

    client = PDBQueryClient(
        cache_dir=Path(args.cache_dir),
        use_cache=not args.no_cache,
        cache_ttl_seconds=max(1, int(args.cache_ttl_days * 24 * 60 * 60)),
        rate_limit_rps=args.rate_limit,
        timeout_seconds=args.timeout,
        max_retries=args.max_retries,
    )

    output_rows: list[dict[str, str]] = []
    for gene in genes:
        result = get_ranked_entries_for_gene(
            gene,
            max_candidates=args.max_candidates,
            max_return=args.max_return,
            resolution_max=args.resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=args.entity_type,
            experimental_methods=args.methods,
            require_ligand=not args.allow_apo,
            ligand_comp_ids=ligand_comp_ids,
            dedupe_sequence_identity=dedupe_sequence_identity,
            quality=quality,
            ligand_filter_mode=ligand_filter_mode,
            gene_to_uniprots=gene_to_uniprots,
            forced_pdbs_by_gene=forced_map,
            include_forced=args.include_forced,
            client=client,
            trivial_ligands=trivial_ligands,
            keep_ligands=keep_ligands,
        )
        LOGGER.info(
            "Gene %s: candidates=%d, post_filter=%d, final=%d",
            gene,
            result.stats["candidates_before_filters"],
            result.stats["selected_after_filters"],
            result.stats["selected_final"],
        )
        output_rows.extend(_entries_to_rows(result.entries))

    for uniprot_id in uniprots:
        result = get_ranked_entries_for_uniprot(
            uniprot_id,
            max_candidates=args.max_candidates,
            max_return=args.max_return,
            resolution_max=args.resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=args.entity_type,
            experimental_methods=args.methods,
            require_ligand=not args.allow_apo,
            ligand_comp_ids=ligand_comp_ids,
            dedupe_sequence_identity=dedupe_sequence_identity,
            quality=quality,
            ligand_filter_mode=ligand_filter_mode,
            client=client,
            trivial_ligands=trivial_ligands,
            keep_ligands=keep_ligands,
        )
        LOGGER.info(
            "UniProt %s: candidates=%d, post_filter=%d, final=%d",
            uniprot_id,
            result.stats["candidates_before_filters"],
            result.stats["selected_after_filters"],
            result.stats["selected_final"],
        )
        output_rows.extend(_entries_to_rows(result.entries))

    _write_output_csv(args.out, output_rows)
    LOGGER.info("Wrote %d rows to %s", len(output_rows), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
