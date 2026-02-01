#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import logging
import os
import random
import re
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence
from urllib.parse import quote, urlencode

import requests  # type: ignore[import-untyped]

LOGGER = logging.getLogger(__name__)

FDA_PLACEHOLDER_RE = re.compile(r"^fda_\d+$", re.IGNORECASE)
UNK_PLACEHOLDER_RE = re.compile(r"^unk_", re.IGNORECASE)
IUPAC_PAREN_DIGIT_RE = re.compile(r"\(\d")
IUPAC_DIGIT_LOCANT_RE = re.compile(r"\b\d{1,3}[a-z]?\s*-\s*[A-Za-z]")
IUPAC_NOS_LOCANT_RE = re.compile(r"\bN-\b|\bO-\b|\bS-\b")
IUPAC_COMMA_LOCANT_RE = re.compile(r"\b\d,\d")
IUPAC_TOKEN_RE = re.compile(
    r"(azanium|ylazanium|methoxy|propyl|butyl|penta|hexa|cyclo|carboxylate|benzamide)",
    re.IGNORECASE,
)
IUPAC_LONG_PUNCT_RE = re.compile(r"[,;()]")
PUBCHEM_PUNCT_RE = re.compile(r"[()\[\],;:/]")
PUBCHEM_ALPHA_RE = re.compile(r"^[A-Za-z][A-Za-z\\s-]*$")
PUBCHEM_DIGIT_RE = re.compile(r"\d")
UNICHEM_CHEMBL_SOURCE_ID = 1


def _clean_text(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _collapse_spaces(value: Any) -> str:
    text = _clean_text(value)
    return " ".join(text.split()) if text else ""


def _normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.strip().lower())


def _choose_header(headers: Sequence[str], *candidates: str) -> str:
    normalized = {_normalize_header(name): name for name in headers}
    for candidate in candidates:
        key = _normalize_header(candidate)
        if key in normalized:
            return normalized[key]
    return ""


_RDKit_AVAILABLE: bool | None = None


def _rdkit_available() -> bool:
    global _RDKit_AVAILABLE
    if _RDKit_AVAILABLE is not None:
        return _RDKit_AVAILABLE
    try:
        pass  # type: ignore[import-untyped]
    except Exception:
        _RDKit_AVAILABLE = False
        return False
    _RDKit_AVAILABLE = True
    return True


def rdkit_inchikey_from_inchi(inchi_value: str) -> str:
    text = _clean_text(inchi_value)
    if not text:
        return ""
    try:
        from rdkit import Chem  # type: ignore[import-untyped]
        from rdkit.Chem import inchi as rdkit_inchi  # type: ignore[import-untyped]
    except Exception:
        return ""
    try:
        mol = Chem.MolFromInchi(text)
    except Exception:
        mol = None
    if mol is None:
        return ""
    try:
        return _clean_text(rdkit_inchi.MolToInchiKey(mol))
    except Exception:
        return ""


def rdkit_inchikey_from_smiles(smiles_value: str) -> str:
    text = _clean_text(smiles_value)
    if not text:
        return ""
    try:
        from rdkit import Chem  # type: ignore[import-untyped]
        from rdkit.Chem import inchi as rdkit_inchi  # type: ignore[import-untyped]
    except Exception:
        return ""
    try:
        mol = Chem.MolFromSmiles(text)
    except Exception:
        mol = None
    if mol is None:
        return ""
    try:
        return _clean_text(rdkit_inchi.MolToInchiKey(mol))
    except Exception:
        return ""


def load_drugcentral_index(path: Path) -> dict[str, str]:
    if not path.exists():
        raise ValueError(f"DrugCentral TSV not found: {path}")
    inchikey_index: dict[str, str] = {}

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        if not headers:
            raise ValueError(f"DrugCentral TSV has no headers: {path}")

        inn_col = _choose_header(headers, "inn", "inn_name")
        name_col = _choose_header(headers, "name")
        drug_name_col = _choose_header(headers, "drug_name", "drugname")
        inchikey_col = _choose_header(
            headers,
            "inchikey",
            "inchi_key",
            "inchikey_2d",
            "inchikey2d",
        )
        inchi_col = _choose_header(headers, "inchi", "inchi_string")
        smiles_col = _choose_header(headers, "smiles", "canonical_smiles")
        rdkit_available = _rdkit_available()
        if not inchikey_col and not rdkit_available:
            LOGGER.warning(
                "DrugCentral TSV missing INCHIKEY column; RDKit not available so cannot "
                "compute InChIKeys. Install RDKit or provide a TSV with INCHIKEY. "
                "DrugCentral stage will have low coverage."
            )
        LOGGER.info(
            "[drugcentral] columns detected: inn=%s name=%s drug_name=%s inchikey=%s inchi=%s smiles=%s",
            inn_col or "-",
            name_col or "-",
            drug_name_col or "-",
            inchikey_col or "-",
            inchi_col or "-",
            smiles_col or "-",
        )

        for row in reader:
            name = _pick_first(
                row.get(inn_col),
                row.get(name_col),
                row.get(drug_name_col),
            )
            if not name:
                continue
            name = _collapse_spaces(name)

            inchikey = _clean_text(row.get(inchikey_col)) if inchikey_col else ""
            if not inchikey:
                inchi_value = _clean_text(row.get(inchi_col)) if inchi_col else ""
                smiles_value = _clean_text(row.get(smiles_col)) if smiles_col else ""
                if not rdkit_available:
                    continue
                if inchi_value:
                    inchikey = rdkit_inchikey_from_inchi(inchi_value)
                elif smiles_value:
                    inchikey = rdkit_inchikey_from_smiles(smiles_value)
            if not inchikey:
                continue
            inchikey = inchikey.strip().upper()
            if inchikey and inchikey not in inchikey_index:
                inchikey_index[inchikey] = name

    LOGGER.info("[drugcentral] indexed %s structures", len(inchikey_index))
    return inchikey_index


def is_bad_display_name(name: str, pubchem_iupac_name: str | None = None) -> bool:
    text = _clean_text(name)
    if not text:
        return True
    if FDA_PLACEHOLDER_RE.match(text):
        return True
    if UNK_PLACEHOLDER_RE.match(text):
        return True
    if pubchem_iupac_name:
        if text.lower() == _clean_text(pubchem_iupac_name).lower():
            return True
    if IUPAC_PAREN_DIGIT_RE.search(text):
        return True
    if IUPAC_DIGIT_LOCANT_RE.search(text):
        return True
    if IUPAC_NOS_LOCANT_RE.search(text):
        return True
    if IUPAC_COMMA_LOCANT_RE.search(text):
        return True
    if IUPAC_TOKEN_RE.search(text):
        return True
    if len(text) > 60 and IUPAC_LONG_PUNCT_RE.search(text):
        return True
    return False


def is_good_display_name(name: str, pubchem_iupac_name: str | None = None) -> bool:
    return not is_bad_display_name(name, pubchem_iupac_name)


@dataclass(frozen=True)
class PubChemResult:
    cid: str | None = None
    title: str | None = None
    iupac_name: str | None = None
    synonyms: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class RxNormResult:
    rxcui: str | None = None
    name: str | None = None


@dataclass(frozen=True)
class ChEMBLResolution:
    name: str = ""
    approved_by_max_phase: bool = False
    approved_by_first_approval: bool = False
    rejected_withdrawn: bool = False


class RateLimiter:
    def __init__(self, qps: float) -> None:
        if qps <= 0:
            raise ValueError("--qps must be > 0.")
        self._interval = 1.0 / qps
        self._lock = threading.Lock()
        self._next_allowed = 0.0

    def wait(self) -> None:
        if self._interval <= 0:
            return
        while True:
            with self._lock:
                now = time.monotonic()
                if now >= self._next_allowed:
                    self._next_allowed = now + self._interval
                    return
                wait_time = self._next_allowed - now
            time.sleep(wait_time)


class CachedJsonClient:
    def __init__(
        self,
        cache_dir: Path,
        sleep: float,
        *,
        cache_mode: str,
        rate_limiter: RateLimiter | None,
        max_retries: int,
        backoff_base: float,
        backoff_max: float,
    ) -> None:
        self.cache_dir = cache_dir
        self.cache_mode = cache_mode
        self.sleep = max(float(sleep), 0.0)
        self.rate_limiter = rate_limiter
        self.max_retries = max_retries
        self.backoff_base = max(float(backoff_base), 0.0)
        self.backoff_max = max(float(backoff_max), 0.0)
        self._lock = threading.Lock()
        if self.cache_mode != "off":
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _cache_path(self, key: str) -> Path:
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def _load_cache(self, cache_path: Path) -> Optional[Any]:
        with self._lock:
            if not cache_path.exists():
                return None
            try:
                cached_payload = json.loads(cache_path.read_text(encoding="utf-8"))
            except Exception:
                return None
        if isinstance(cached_payload, dict) and "payload" in cached_payload:
            return cached_payload.get("payload")
        if isinstance(cached_payload, (dict, list)):
            return cached_payload
        return None

    def _write_cache(
        self, cache_path: Path, url: str, payload: Any, *, overwrite: bool
    ) -> None:
        with self._lock:
            if cache_path.exists() and not overwrite:
                return
            try:
                tmp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
                tmp_path.write_text(
                    json.dumps({"url": url, "payload": payload}, ensure_ascii=True),
                    encoding="utf-8",
                )
                tmp_path.replace(cache_path)
            except Exception:
                pass

    def _parse_retry_after(
        self, response: Optional[requests.Response]
    ) -> Optional[float]:
        if response is None:
            return None
        try:
            value = response.headers.get("Retry-After")
        except Exception:
            return None
        if not value:
            return None
        try:
            return float(value)
        except Exception:
            return None

    def _compute_backoff(
        self, retry_count: int, response: Optional[requests.Response]
    ) -> float:
        base = self.backoff_base * (2 ** (retry_count - 1))
        jitter = 0.0
        if base > 0:
            jitter = random.uniform(0.0, base * 0.1)
        wait_time = base + jitter
        retry_after = self._parse_retry_after(response)
        if retry_after is not None:
            wait_time = max(wait_time, retry_after)
        if self.backoff_max > 0:
            wait_time = min(wait_time, self.backoff_max)
        return max(wait_time, 0.0)

    def _wait_before_request(self) -> None:
        if self.rate_limiter is not None:
            self.rate_limiter.wait()
            return
        if self.sleep:
            time.sleep(self.sleep)

    def _request_json(
        self,
        url: str,
        *,
        method: str = "get",
        json_payload: Optional[dict[str, Any]] = None,
    ) -> Optional[Any]:
        retries = 0
        while True:
            self._wait_before_request()
            try:
                if method == "post":
                    response = requests.post(url, json=json_payload, timeout=10)
                else:
                    response = requests.get(url, timeout=10)
            except requests.RequestException as exc:
                if retries >= self.max_retries:
                    LOGGER.warning("[http] exc=%s url=%s", exc.__class__.__name__, url)
                    return None
                retries += 1
                wait_time = self._compute_backoff(retries, None)
                LOGGER.warning(
                    "[http] exc=%s retry=%s wait=%.2f url=%s",
                    exc.__class__.__name__,
                    retries,
                    wait_time,
                    url,
                )
                time.sleep(wait_time)
                continue

            status = response.status_code
            if status == 200:
                try:
                    payload = response.json()
                except Exception:
                    LOGGER.warning("[http] status=200 url=%s reason=json_error", url)
                    return None
                if isinstance(payload, (dict, list)):
                    return payload
                LOGGER.warning("[http] status=200 url=%s reason=non_object_json", url)
                return None

            if status in (429, 503):
                if retries >= self.max_retries:
                    LOGGER.warning(
                        "[http] status=%s url=%s reason=max_retries", status, url
                    )
                    return None
                retries += 1
                wait_time = self._compute_backoff(retries, response)
                LOGGER.warning(
                    "[http] status=%s retry=%s wait=%.2f url=%s",
                    status,
                    retries,
                    wait_time,
                    url,
                )
                time.sleep(wait_time)
                continue

            LOGGER.warning("[http] status=%s url=%s reason=non_200", status, url)
            return None

    def get_json(self, url: str) -> Optional[Any]:
        cache_path = self._cache_path(url) if self.cache_mode != "off" else None
        if self.cache_mode == "use" and cache_path is not None:
            cached_payload = self._load_cache(cache_path)
            if cached_payload is not None:
                return cached_payload

        payload = self._request_json(url)
        if payload is None:
            return None
        if self.cache_mode != "off" and cache_path is not None:
            self._write_cache(
                cache_path,
                url,
                payload,
                overwrite=self.cache_mode == "refresh",
            )
        return payload

    def post_json(self, url: str, payload: dict[str, Any]) -> Optional[Any]:
        cache_key = (
            f"POST {url} {json.dumps(payload, sort_keys=True, separators=(',', ':'))}"
        )
        cache_path = self._cache_path(cache_key) if self.cache_mode != "off" else None
        if self.cache_mode == "use" and cache_path is not None:
            cached_payload = self._load_cache(cache_path)
            if cached_payload is not None:
                return cached_payload

        response_payload = self._request_json(url, method="post", json_payload=payload)
        if response_payload is None:
            return None
        if self.cache_mode != "off" and cache_path is not None:
            self._write_cache(
                cache_path,
                url,
                response_payload,
                overwrite=self.cache_mode == "refresh",
            )
        return response_payload


def _clear_cache_dir(cache_dir: Path) -> None:
    if not cache_dir.exists():
        return
    if not cache_dir.is_dir():
        raise ValueError(f"Cache directory is not a directory: {cache_dir}")
    for child in cache_dir.iterdir():
        try:
            if child.is_symlink() or child.is_file():
                child.unlink()
            elif child.is_dir():
                shutil.rmtree(child)
        except FileNotFoundError:
            continue


def _resolution_key(row: dict[str, Any]) -> tuple[str | None, str | None]:
    for key in (
        "inchikey",
        "remark_inchikey",
        "smiles_neutral",
        "smiles",
        "remark_smiles",
    ):
        value = _clean_text(row.get(key))
        if value:
            return key, value
    return None, None


def _parse_existing_pubchem_cid(row: dict[str, Any]) -> Optional[str]:
    raw = _clean_text(row.get("pubchem_cid_resolved"))
    if not raw:
        return None
    try:
        cid = int(float(raw))
    except Exception:
        return None
    if cid <= 0:
        return None
    return str(cid)


def _parse_pubchem_cid(payload: Optional[dict[str, Any]]) -> Optional[str]:
    if not payload:
        return None
    try:
        cids = payload["IdentifierList"]["CID"]
    except Exception:
        return None
    if not cids:
        return None
    return str(cids[0])


def _parse_pubchem_properties(
    payload: Optional[dict[str, Any]],
) -> tuple[Optional[str], Optional[str]]:
    if not payload:
        return None, None
    try:
        props = payload["PropertyTable"]["Properties"]
        if not props:
            return None, None
        entry = props[0]
        return entry.get("Title"), entry.get("IUPACName")
    except Exception:
        return None, None


def _parse_pubchem_synonyms(payload: Optional[dict[str, Any]]) -> list[str]:
    if not payload:
        return []
    try:
        info = payload["InformationList"]["Information"]
        if not info:
            return []
        entry = info[0]
        syns = entry.get("Synonym", [])
        return [str(s) for s in syns if str(s).strip()]
    except Exception:
        return []


def _row_inchikey(row: dict[str, Any]) -> str:
    return _pick_first(row.get("inchikey"), row.get("remark_inchikey")).upper()


def _dedupe_preserve(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        text = _clean_text(value)
        if not text:
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        ordered.append(text)
    return ordered


def _parse_unichem_items(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        for key in ("compounds", "results", "payload"):
            items = payload.get(key)
            if isinstance(items, list):
                return [item for item in items if isinstance(item, dict)]
    return []


def _parse_unichem_chembl_ids(payload: Any) -> list[str]:
    chembl_ids: list[str] = []
    for item in _parse_unichem_items(payload):
        sources = item.get("sources") if isinstance(item, dict) else None
        if isinstance(sources, list):
            for source in sources:
                if not isinstance(source, dict):
                    continue
                source_id = (
                    source.get("id") or source.get("src_id") or source.get("source_id")
                )
                if source_id is None:
                    continue
                try:
                    source_id_int = int(source_id)
                except Exception:
                    continue
                if source_id_int != UNICHEM_CHEMBL_SOURCE_ID:
                    continue
                chembl_id = (
                    source.get("compoundId")
                    or source.get("compound_id")
                    or source.get("src_compound_id")
                    or source.get("source_compound_id")
                )
                if chembl_id:
                    chembl_ids.append(str(chembl_id))
            continue

        lowered = {str(k).lower(): v for k, v in item.items()}
        source_id = (
            lowered.get("src_id") or lowered.get("source_id") or lowered.get("sourceid")
        )
        if source_id is None:
            continue
        try:
            source_id_int = int(source_id)
        except Exception:
            continue
        if source_id_int != UNICHEM_CHEMBL_SOURCE_ID:
            continue
        chembl_id = (
            lowered.get("src_compound_id")
            or lowered.get("source_compound_id")
            or lowered.get("src_compound")
        )
        if chembl_id:
            chembl_ids.append(str(chembl_id))
    return _dedupe_preserve(chembl_ids)


def _extract_chembl_max_phase(payload: dict[str, Any]) -> Optional[int]:
    for key in ("max_phase", "maxphase"):
        if key in payload:
            value = payload.get(key)
            if value is None:
                continue
            try:
                return int(value)
            except Exception:
                return None
    try:
        props = payload.get("molecule_properties", {})
    except Exception:
        return None
    if isinstance(props, dict):
        value = props.get("max_phase") or props.get("maxphase")
        if value is None:
            return None
        try:
            return int(value)
        except Exception:
            return None
    return None


def _parse_boolish(value: Any) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    text = _clean_text(value).lower()
    if text in {"1", "true", "yes", "y", "t"}:
        return True
    if text in {"0", "false", "no", "n", "f"}:
        return False
    return None


def _extract_chembl_first_approval(payload: dict[str, Any]) -> Optional[int]:
    value = payload.get("first_approval")
    if value is None:
        return None
    try:
        year = int(float(str(value).strip()))
    except Exception:
        return None
    return year if year > 0 else None


def _extract_chembl_withdrawn(payload: dict[str, Any]) -> Optional[bool]:
    for key in ("withdrawn_flag", "withdrawn"):
        if key in payload:
            parsed = _parse_boolish(payload.get(key))
            if parsed is not None:
                return parsed
    return None


def _extract_chembl_synonyms(payload: dict[str, Any]) -> list[str]:
    synonyms: list[str] = []
    raw = payload.get("molecule_synonyms", [])
    if not isinstance(raw, list):
        return synonyms
    for item in raw:
        if not isinstance(item, dict):
            continue
        value = item.get("molecule_synonym") or item.get("synonym") or item.get("name")
        if value:
            synonyms.append(str(value))
    return synonyms


def resolve_chembl_name_from_unichem(
    inchikey: str,
    *,
    unichem_client: CachedJsonClient,
    chembl_client: CachedJsonClient,
) -> ChEMBLResolution:
    key = _clean_text(inchikey).upper()
    if not key:
        return ChEMBLResolution()
    payload = unichem_client.post_json(
        "https://www.ebi.ac.uk/unichem/api/v1/compounds",
        {"compound": key, "type": "inchikey"},
    )
    chembl_ids = _parse_unichem_chembl_ids(payload)
    rejected_withdrawn = False
    current_year = dt.date.today().year
    for chembl_id in chembl_ids:
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
        chembl_payload = chembl_client.get_json(url)
        if not isinstance(chembl_payload, dict):
            continue
        withdrawn = _extract_chembl_withdrawn(chembl_payload)
        if withdrawn is True:
            rejected_withdrawn = True
            continue
        max_phase = _extract_chembl_max_phase(chembl_payload)
        approved_by_max_phase = max_phase == 4
        first_approval = _extract_chembl_first_approval(chembl_payload)
        approved_by_first_approval = False
        if first_approval is not None and 1900 <= first_approval <= current_year + 1:
            approved_by_first_approval = True
        if not (approved_by_max_phase or approved_by_first_approval):
            continue
        pref_name = _clean_text(chembl_payload.get("pref_name"))
        if pref_name:
            return ChEMBLResolution(
                name=pref_name,
                approved_by_max_phase=approved_by_max_phase,
                approved_by_first_approval=approved_by_first_approval,
                rejected_withdrawn=rejected_withdrawn,
            )
        for synonym in _extract_chembl_synonyms(chembl_payload):
            if is_good_display_name(synonym):
                return ChEMBLResolution(
                    name=_clean_text(synonym),
                    approved_by_max_phase=approved_by_max_phase,
                    approved_by_first_approval=approved_by_first_approval,
                    rejected_withdrawn=rejected_withdrawn,
                )
    if rejected_withdrawn:
        return ChEMBLResolution(rejected_withdrawn=True)
    return ChEMBLResolution()


def resolve_pubchem(
    row: dict[str, Any],
    client: CachedJsonClient,
    *,
    require_existing_cid: bool = False,
) -> PubChemResult:
    cid = _parse_existing_pubchem_cid(row)
    if not cid:
        if require_existing_cid:
            return PubChemResult()
        key_type, key_value = _resolution_key(row)
        if not key_type or not key_value:
            return PubChemResult()

        if key_type in ("inchikey", "remark_inchikey"):
            url = (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
                f"{quote(key_value, safe='')}/cids/JSON"
            )
        else:
            url = (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
                f"{quote(key_value, safe='')}/cids/JSON"
            )

        cid = _parse_pubchem_cid(client.get_json(url))
        if not cid:
            return PubChemResult()

    props_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/Title,IUPACName/JSON"
    )
    syn_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" f"{cid}/synonyms/JSON"
    )
    title, iupac_name = _parse_pubchem_properties(client.get_json(props_url))
    synonyms = _parse_pubchem_synonyms(client.get_json(syn_url))
    return PubChemResult(cid=cid, title=title, iupac_name=iupac_name, synonyms=synonyms)


def _parse_rxnorm_candidate(payload: Optional[dict[str, Any]]) -> Optional[str]:
    if not payload:
        return None
    try:
        candidates = payload["approximateGroup"]["candidate"]
    except Exception:
        return None
    if not candidates:
        return None

    def score_of(item: dict[str, Any]) -> int:
        try:
            return int(item.get("score", 0))
        except Exception:
            return 0

    best = max(candidates, key=score_of)
    return best.get("rxcui")


def _parse_rxnorm_name(payload: Optional[dict[str, Any]]) -> Optional[str]:
    if not payload:
        return None
    try:
        props = payload["properties"]
        return props.get("name")
    except Exception:
        return None


def resolve_rxnorm(term: str, client: CachedJsonClient) -> RxNormResult:
    query = _clean_text(term)
    if not query:
        return RxNormResult()
    url = f"https://rxnav.nlm.nih.gov/REST/approximateTerm.json?{urlencode({'term': query, 'maxEntries': 1})}"
    rxcui = _parse_rxnorm_candidate(client.get_json(url))
    if not rxcui:
        return RxNormResult()
    props_url = f"https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/properties.json"
    name = _parse_rxnorm_name(client.get_json(props_url))
    return RxNormResult(rxcui=rxcui, name=name)


def _maybe_set_field(row: dict[str, Any], field: str, value: Optional[str]) -> bool:
    if field not in row:
        return False
    if _clean_text(row.get(field)):
        return False
    if not _clean_text(value):
        return False
    row[field] = value
    return True


def _pick_first(*values: Optional[str]) -> str:
    for value in values:
        text = _clean_text(value)
        if text:
            return text
    return ""


def _split_synonyms(raw: str) -> list[str]:
    if not raw:
        return []
    parts = [p.strip() for p in re.split(r"[;|]", raw) if p.strip()]
    return parts


def _is_standard_casing(text: str) -> bool:
    if not text:
        return False
    if text.isupper():
        return False
    if text.islower():
        return True
    words = [w for w in text.split() if w]
    if words and all(w[:1].isupper() and w[1:].islower() for w in words):
        return True
    return text[:1].isupper()


def _score_pubchem_candidate(
    candidate: str, pubchem_iupac_name: str | None
) -> Optional[int]:
    if is_bad_display_name(candidate, pubchem_iupac_name):
        return None
    score = 0
    text = _clean_text(candidate)
    if len(text) <= 30:
        score += 3
    if PUBCHEM_ALPHA_RE.match(text):
        score += 3
    if _is_standard_casing(text):
        score += 2

    has_punct = bool(PUBCHEM_PUNCT_RE.search(text))
    has_digits = bool(PUBCHEM_DIGIT_RE.search(text))
    if has_punct:
        score -= 1
        if has_digits:
            score -= 3
    elif has_digits:
        score -= 1

    if (
        IUPAC_PAREN_DIGIT_RE.search(text)
        or IUPAC_DIGIT_LOCANT_RE.search(text)
        or IUPAC_NOS_LOCANT_RE.search(text)
        or IUPAC_COMMA_LOCANT_RE.search(text)
    ):
        score -= 2
    if IUPAC_TOKEN_RE.search(text):
        score -= 1
    return score


def _select_pubchem_candidate(
    candidates: Iterable[str], pubchem_iupac_name: str | None
) -> str:
    best: tuple[int, int, str] | None = None
    best_name = ""
    for candidate in candidates:
        score = _score_pubchem_candidate(candidate, pubchem_iupac_name)
        if score is None:
            continue
        text = _clean_text(candidate)
        ranking = (score, -len(text), text.lower())
        if best is None or ranking > best:
            best = ranking
            best_name = text
    return best_name


def _normalize_numeric(text: str) -> str:
    raw = _clean_text(text)
    if not raw:
        return ""
    try:
        return str(int(float(raw)))
    except Exception:
        return ""


def _fallback_display_name(row: dict[str, Any], row_index: int) -> str:
    file_num = _normalize_numeric(_clean_text(row.get("file_num")))
    if file_num:
        return f"UNK_{file_num}"
    sdf_index = _normalize_numeric(_clean_text(row.get("sdf_index")))
    if sdf_index:
        return f"UNK_{sdf_index}"
    inchikey = _pick_first(row.get("inchikey"), row.get("remark_inchikey"))
    if inchikey:
        return f"UNK_{inchikey[:8]}"
    return f"UNK_ROW{row_index}"


def _row_identifier(row: dict[str, Any], row_index: int) -> str:
    return _pick_first(
        row.get("sdf_title"),
        row.get("path"),
        row.get("file_num"),
        row.get("sdf_index"),
        f"row_{row_index}",
    )


def _select_pubchem_display_name(
    row: dict[str, Any],
    pubchem: PubChemResult,
    pubchem_iupac: str,
) -> str:
    candidates: list[str] = []
    title = _pick_first(
        pubchem.title, row.get("pubchem_record_title"), row.get("pubchem_name")
    )
    if title:
        candidates.append(title)
    synonyms = pubchem.synonyms
    if not synonyms:
        synonyms = _split_synonyms(_clean_text(row.get("pubchem_synonyms")))
    candidates.extend(synonyms)
    return _select_pubchem_candidate(_dedupe_preserve(candidates), pubchem_iupac)


@dataclass
class Summary:
    total_rows: int = 0
    rows_evaluated: int = 0
    rows_updated: int = 0
    drugcentral_resolved: int = 0
    chembl_resolved: int = 0
    chembl_approved_by_max_phase: int = 0
    chembl_approved_by_first_approval: int = 0
    chembl_rejected_withdrawn: int = 0
    pubchem_reranked: int = 0
    pubchem_resolved: int = 0
    rxnorm_resolved: int = 0
    rxnorm_promoted: int = 0
    fallbacks_used: int = 0


@dataclass
class RowResult:
    row_index: int
    row: dict[str, Any]
    row_evaluated: bool
    row_updated: bool
    drugcentral_resolved: bool
    chembl_resolved: bool
    chembl_approved_by_max_phase: bool
    chembl_approved_by_first_approval: bool
    chembl_rejected_withdrawn: bool
    pubchem_reranked: bool
    pubchem_resolved: bool
    rxnorm_resolved: bool
    rxnorm_promoted: bool
    used_fallback: bool
    final_bad: bool
    row_identifier: str


def _process_row(
    row_index: int,
    row: dict[str, Any],
    only_fix_bad_display_names: bool,
    pubchem_client: CachedJsonClient | None,
    rxnorm_client: CachedJsonClient | None,
    unichem_client: CachedJsonClient | None,
    chembl_client: CachedJsonClient | None,
    drugcentral_index: dict[str, str],
    enable_unichem: bool,
    enable_pubchem_rerank: bool,
) -> RowResult:
    current_display = _clean_text(row.get("display_name"))
    pubchem_iupac = _clean_text(row.get("pubchem_iupac_name"))
    display_bad = is_bad_display_name(current_display, pubchem_iupac)
    evaluate = (not only_fix_bad_display_names) or display_bad
    row_updated = False
    used_fallback = False
    row_evaluated = False
    pubchem_resolved = False
    rxnorm_resolved = False
    drugcentral_resolved = False
    chembl_resolved = False
    pubchem_reranked = False
    chembl_approved_by_max_phase = False
    chembl_approved_by_first_approval = False
    chembl_rejected_withdrawn = False
    rxnorm_promoted = False

    pubchem = PubChemResult()
    if evaluate:
        row_evaluated = True
        if enable_pubchem_rerank and pubchem_client and rxnorm_client:
            pubchem = resolve_pubchem(row, pubchem_client)
            if pubchem.cid:
                pubchem_resolved = True
            row_updated |= _maybe_set_field(row, "pubchem_cid_resolved", pubchem.cid)
            row_updated |= _maybe_set_field(row, "pubchem_record_title", pubchem.title)
            row_updated |= _maybe_set_field(row, "pubchem_name", pubchem.title)
            row_updated |= _maybe_set_field(
                row, "pubchem_iupac_name", pubchem.iupac_name
            )
            if pubchem.synonyms:
                syn_value = "; ".join(pubchem.synonyms)
                row_updated |= _maybe_set_field(row, "pubchem_synonyms", syn_value)

            rxnorm_name = _clean_text(row.get("rxnorm_generic_name"))
            if not rxnorm_name:
                candidate = _pick_first(
                    pubchem.title,
                    row.get("pubchem_record_title"),
                    row.get("pubchem_name"),
                    row.get("generic_name"),
                    current_display,
                )
                rxnorm = resolve_rxnorm(candidate, rxnorm_client)
                if rxnorm.rxcui or rxnorm.name:
                    rxnorm_resolved = True
                row_updated |= _maybe_set_field(row, "rxnorm_generic_name", rxnorm.name)
                row_updated |= _maybe_set_field(row, "rxnorm_rxcui", rxnorm.rxcui)

    pubchem_iupac = _clean_text(row.get("pubchem_iupac_name"))
    display_bad = is_bad_display_name(current_display, pubchem_iupac)
    if evaluate and display_bad:
        resolved_name = ""
        if drugcentral_index:
            inchikey = _row_inchikey(row)
            candidate = drugcentral_index.get(inchikey, "")
            if candidate and is_good_display_name(candidate, pubchem_iupac):
                if _maybe_set_field(row, "generic_name", candidate):
                    row_updated = True
                if candidate != current_display:
                    row["display_name"] = candidate
                    row_updated = True
                resolved_name = candidate
                drugcentral_resolved = True
                LOGGER.info(
                    "[drugcentral] row=%s inchikey=%s name=%s",
                    row_index,
                    inchikey,
                    candidate,
                )

        if not resolved_name:
            rxnorm_name = _clean_text(row.get("rxnorm_generic_name"))
            if rxnorm_name and is_good_display_name(rxnorm_name, pubchem_iupac):
                if rxnorm_name != current_display:
                    row["display_name"] = rxnorm_name
                    row_updated = True
                resolved_name = rxnorm_name
                rxnorm_promoted = True

        if not resolved_name and enable_unichem and unichem_client and chembl_client:
            inchikey = _row_inchikey(row)
            chembl_resolution = resolve_chembl_name_from_unichem(
                inchikey, unichem_client=unichem_client, chembl_client=chembl_client
            )
            if chembl_resolution.rejected_withdrawn:
                chembl_rejected_withdrawn = True
            candidate = chembl_resolution.name
            if candidate and is_good_display_name(candidate, pubchem_iupac):
                if candidate != current_display:
                    row["display_name"] = candidate
                    row_updated = True
                resolved_name = candidate
                chembl_resolved = True
                chembl_approved_by_max_phase = chembl_resolution.approved_by_max_phase
                chembl_approved_by_first_approval = (
                    chembl_resolution.approved_by_first_approval
                )
                LOGGER.info(
                    "[chembl] row=%s inchikey=%s name=%s",
                    row_index,
                    inchikey,
                    candidate,
                )

        if (
            not resolved_name
            and enable_pubchem_rerank
            and pubchem_client
            and pubchem.cid
        ):
            pubchem_iupac = _pick_first(
                row.get("pubchem_iupac_name"), pubchem.iupac_name
            )
            candidate = _select_pubchem_display_name(row, pubchem, pubchem_iupac)
            if candidate and is_good_display_name(candidate, pubchem_iupac):
                if candidate != current_display:
                    row["display_name"] = candidate
                    row_updated = True
                resolved_name = candidate
                pubchem_reranked = True
                LOGGER.info(
                    "[pubchem] row=%s cid=%s name=%s", row_index, pubchem.cid, candidate
                )

        if not resolved_name:
            fallback = _fallback_display_name(row, row_index)
            if fallback and fallback != current_display:
                row["display_name"] = fallback
                row_updated = True
            used_fallback = True

    final_display = _clean_text(row.get("display_name"))
    final_bad = is_bad_display_name(
        final_display, _clean_text(row.get("pubchem_iupac_name"))
    )
    row_identifier = (
        _row_identifier(row, row_index) if (used_fallback or final_bad) else ""
    )
    return RowResult(
        row_index=row_index,
        row=row,
        row_evaluated=row_evaluated,
        row_updated=row_updated,
        drugcentral_resolved=drugcentral_resolved,
        chembl_resolved=chembl_resolved,
        chembl_approved_by_max_phase=chembl_approved_by_max_phase,
        chembl_approved_by_first_approval=chembl_approved_by_first_approval,
        chembl_rejected_withdrawn=chembl_rejected_withdrawn,
        pubchem_reranked=pubchem_reranked,
        pubchem_resolved=pubchem_resolved,
        rxnorm_resolved=rxnorm_resolved,
        rxnorm_promoted=rxnorm_promoted,
        used_fallback=used_fallback,
        final_bad=final_bad,
        row_identifier=row_identifier,
    )


def _process_row_item(
    item: tuple[int, dict[str, Any]],
    only_fix_bad_display_names: bool,
    pubchem_client: CachedJsonClient | None,
    rxnorm_client: CachedJsonClient | None,
    unichem_client: CachedJsonClient | None,
    chembl_client: CachedJsonClient | None,
    drugcentral_index: dict[str, str],
    enable_unichem: bool,
    enable_pubchem_rerank: bool,
) -> RowResult:
    row_index, row = item
    return _process_row(
        row_index,
        row,
        only_fix_bad_display_names=only_fix_bad_display_names,
        pubchem_client=pubchem_client,
        rxnorm_client=rxnorm_client,
        unichem_client=unichem_client,
        chembl_client=chembl_client,
        drugcentral_index=drugcentral_index,
        enable_unichem=enable_unichem,
        enable_pubchem_rerank=enable_pubchem_rerank,
    )


def fill_mapping_names(
    in_csv: Path,
    out_csv: Path,
    *,
    only_fix_bad_display_names: bool,
    require_all_named: bool,
    report_path: Path,
    drugcentral_structures_tsv: Path | None,
    enable_unichem: bool,
    enable_pubchem_rerank: bool,
    cache_dir: Path,
    cache_mode: str,
    clear_cache: bool,
    sleep: float,
    rate_limiter: RateLimiter | None,
    max_workers: int,
    max_retries: int,
    backoff_base: float,
    backoff_max: float,
) -> int:
    if cache_mode not in {"use", "refresh", "off"}:
        raise ValueError("--cache_mode must be one of: use, refresh, off.")
    if cache_mode != "off":
        cache_dir.mkdir(parents=True, exist_ok=True)
        if clear_cache:
            _clear_cache_dir(cache_dir)
    elif clear_cache:
        LOGGER.warning("[cache] clear_cache ignored because cache_mode=off")

    drugcentral_index: dict[str, str] = {}
    if drugcentral_structures_tsv:
        drugcentral_index = load_drugcentral_index(drugcentral_structures_tsv)
    else:
        LOGGER.warning(
            "[drugcentral] --drugcentral_structures_tsv not provided; skipping DrugCentral stage."
        )

    pubchem_client: CachedJsonClient | None = None
    rxnorm_client: CachedJsonClient | None = None
    if enable_pubchem_rerank:
        pubchem_client = CachedJsonClient(
            cache_dir / "pubchem",
            sleep,
            cache_mode=cache_mode,
            rate_limiter=rate_limiter,
            max_retries=max_retries,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
        )
        rxnorm_client = CachedJsonClient(
            cache_dir / "rxnorm",
            sleep,
            cache_mode=cache_mode,
            rate_limiter=rate_limiter,
            max_retries=max_retries,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
        )

    unichem_client: CachedJsonClient | None = None
    chembl_client: CachedJsonClient | None = None
    if enable_unichem:
        unichem_client = CachedJsonClient(
            cache_dir / "unichem",
            sleep,
            cache_mode=cache_mode,
            rate_limiter=rate_limiter,
            max_retries=max_retries,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
        )
        chembl_client = CachedJsonClient(
            cache_dir / "chembl",
            sleep,
            cache_mode=cache_mode,
            rate_limiter=rate_limiter,
            max_retries=max_retries,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
        )

    unresolved: list[str] = []
    has_bad = False
    summary = Summary()

    in_place = in_csv.resolve() == out_csv.resolve()
    temp_path = out_csv.with_suffix(out_csv.suffix + ".tmp") if in_place else None
    write_path = temp_path if temp_path else out_csv
    write_path.parent.mkdir(parents=True, exist_ok=True)

    if max_workers < 1:
        raise ValueError("--max_workers must be >= 1.")

    with in_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError("Input CSV has no headers.")

        with write_path.open("w", newline="", encoding="utf-8") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
            writer.writeheader()
            process_row = partial(
                _process_row_item,
                only_fix_bad_display_names=only_fix_bad_display_names,
                pubchem_client=pubchem_client,
                rxnorm_client=rxnorm_client,
                unichem_client=unichem_client,
                chembl_client=chembl_client,
                drugcentral_index=drugcentral_index,
                enable_unichem=enable_unichem,
                enable_pubchem_rerank=enable_pubchem_rerank,
            )

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                for result in executor.map(process_row, enumerate(reader, start=1)):
                    summary.total_rows += 1
                    if result.row_evaluated:
                        summary.rows_evaluated += 1
                    if result.drugcentral_resolved:
                        summary.drugcentral_resolved += 1
                    if result.chembl_resolved:
                        summary.chembl_resolved += 1
                    if result.chembl_approved_by_max_phase:
                        summary.chembl_approved_by_max_phase += 1
                    if result.chembl_approved_by_first_approval:
                        summary.chembl_approved_by_first_approval += 1
                    if result.chembl_rejected_withdrawn:
                        summary.chembl_rejected_withdrawn += 1
                    if result.pubchem_reranked:
                        summary.pubchem_reranked += 1
                    if result.pubchem_resolved:
                        summary.pubchem_resolved += 1
                    if result.rxnorm_resolved:
                        summary.rxnorm_resolved += 1
                    if result.rxnorm_promoted:
                        summary.rxnorm_promoted += 1
                    if result.used_fallback:
                        summary.fallbacks_used += 1
                    if result.row_updated:
                        summary.rows_updated += 1
                    if result.row_identifier:
                        unresolved.append(result.row_identifier)
                    if result.final_bad:
                        has_bad = True

                    writer.writerow(
                        {name: result.row.get(name, "") for name in fieldnames}
                    )

    if temp_path:
        temp_path.replace(out_csv)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    ordered_unresolved = list(dict.fromkeys(unresolved))
    report_path.write_text(
        "\n".join(ordered_unresolved) + ("\n" if ordered_unresolved else ""),
        encoding="utf-8",
    )

    print(f"total rows: {summary.total_rows}")
    print(f"rows evaluated: {summary.rows_evaluated}")
    print(f"rows updated: {summary.rows_updated}")
    print(f"drugcentral resolved: {summary.drugcentral_resolved}")
    print(f"chembl resolved: {summary.chembl_resolved}")
    print(f"chembl approved by max_phase: {summary.chembl_approved_by_max_phase}")
    print(
        f"chembl approved by first_approval: {summary.chembl_approved_by_first_approval}"
    )
    if summary.chembl_rejected_withdrawn:
        print(f"chembl rejected withdrawn: {summary.chembl_rejected_withdrawn}")
    print(f"pubchem reranked: {summary.pubchem_reranked}")
    print(f"pubchem resolved: {summary.pubchem_resolved}")
    print(f"rxnorm resolved: {summary.rxnorm_resolved}")
    print(f"rxnorm promoted: {summary.rxnorm_promoted}")
    print(f"fallbacks used: {summary.fallbacks_used}")

    if require_all_named and has_bad:
        return 2
    return 0


def _parse_bool(value: str) -> bool:
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError("Expected a boolean value (true/false).")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fill FDA mapping display names from DrugCentral/ChEMBL/PubChem without overwriting inputs."
        )
    )
    parser.add_argument("--in_csv", required=True, help="Input FDA mapping CSV.")
    parser.add_argument(
        "--out_csv", help="Output CSV path (required unless --inplace)."
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Overwrite --in_csv after creating a timestamped .bak backup.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.25,
        help="Sleep between external requests when --qps is not set (deprecated).",
    )
    parser.add_argument(
        "--qps",
        type=float,
        default=2.0,
        help="Global request rate limit (queries per second) across all threads.",
    )
    parser.add_argument(
        "--only_fix_bad_display_names",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Only update rows with placeholder/IUPAC-like display names (default: true).",
    )
    parser.add_argument(
        "--require_all_named",
        action="store_true",
        help="Exit with code 2 if any row still has a placeholder/IUPAC-like display_name.",
    )
    parser.add_argument(
        "--fail_if_unresolved",
        action="store_true",
        help="Deprecated alias for --require_all_named.",
    )
    parser.add_argument(
        "--drugcentral_structures_tsv",
        help=(
            "DrugCentral structures.smiles.tsv path (download the SMILES and InChI file "
            "from DrugCentral)."
        ),
    )
    parser.add_argument(
        "--enable_unichem",
        type=_parse_bool,
        default=True,
        help="Enable UniChem->ChEMBL max_phase=4 resolution (default: true).",
    )
    parser.add_argument(
        "--enable_pubchem_rerank",
        type=_parse_bool,
        default=True,
        help="Enable PubChem synonym re-ranking for CID-known rows (default: true).",
    )
    parser.add_argument(
        "--report_path",
        help="Path for unresolved report (default: <out_csv>.unresolved.txt).",
    )
    parser.add_argument(
        "--cache_dir",
        default=".cache/fda_name_fill/",
        help="Directory for JSON response caches.",
    )
    parser.add_argument(
        "--cache_mode",
        choices=("use", "refresh", "off"),
        default="use",
        help="Cache behavior: use (default), refresh (overwrite), or off (no cache).",
    )
    parser.add_argument(
        "--clear_cache",
        action="store_true",
        help="Clear cache_dir contents before processing (ignored when cache_mode=off).",
    )
    parser.add_argument(
        "--max_workers",
        type=int,
        default=min(32, os.cpu_count() or 4),
        help="Maximum worker threads for parallel resolution (default: min(32, cpu_count)).",
    )
    parser.add_argument(
        "--max_retries",
        type=int,
        default=5,
        help="Maximum retries for 429/503/network errors (default: 5).",
    )
    parser.add_argument(
        "--backoff_base",
        type=float,
        default=0.5,
        help="Base seconds for exponential backoff (default: 0.5).",
    )
    parser.add_argument(
        "--backoff_max",
        type=float,
        default=8.0,
        help="Maximum seconds for backoff or Retry-After (default: 8.0).",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    argv_list = list(argv) if argv is not None else sys.argv[1:]
    qps_provided = any(arg == "--qps" or arg.startswith("--qps=") for arg in argv_list)
    sleep_provided = any(
        arg == "--sleep" or arg.startswith("--sleep=") for arg in argv_list
    )

    in_path = Path(args.in_csv)
    if not in_path.exists():
        parser.error(f"--in_csv not found: {in_path}")

    if args.inplace:
        if args.out_csv:
            out_path = Path(args.out_csv)
            if out_path.resolve() != in_path.resolve():
                parser.error(
                    "--inplace requires --out_csv to match --in_csv when provided."
                )
        out_path = in_path
        timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = in_path.with_name(f"{in_path.name}.{timestamp}.bak")
        shutil.copy2(in_path, backup_path)
    else:
        if not args.out_csv:
            parser.error("--out_csv is required unless --inplace is set.")
        out_path = Path(args.out_csv)
        if out_path.resolve() == in_path.resolve():
            parser.error("--out_csv must differ from --in_csv unless --inplace is set.")

    report_path = (
        Path(args.report_path)
        if args.report_path
        else Path(f"{out_path}.unresolved.txt")
    )
    drugcentral_path = (
        Path(args.drugcentral_structures_tsv)
        if args.drugcentral_structures_tsv
        else None
    )
    if drugcentral_path and not drugcentral_path.exists():
        parser.error(f"--drugcentral_structures_tsv not found: {drugcentral_path}")
    cache_dir = Path(args.cache_dir)
    cache_mode = str(args.cache_mode)
    require_all_named = bool(args.require_all_named or args.fail_if_unresolved)

    if int(args.max_retries) < 0:
        parser.error("--max_retries must be >= 0.")
    if float(args.backoff_base) < 0:
        parser.error("--backoff_base must be >= 0.")
    if float(args.backoff_max) < 0:
        parser.error("--backoff_max must be >= 0.")

    rate_limiter: RateLimiter | None = None
    if qps_provided or not sleep_provided:
        if float(args.qps) <= 0:
            parser.error("--qps must be > 0.")
        rate_limiter = RateLimiter(float(args.qps))

    try:
        return fill_mapping_names(
            in_path,
            out_path,
            only_fix_bad_display_names=bool(args.only_fix_bad_display_names),
            require_all_named=require_all_named,
            report_path=report_path,
            drugcentral_structures_tsv=drugcentral_path,
            enable_unichem=bool(args.enable_unichem),
            enable_pubchem_rerank=bool(args.enable_pubchem_rerank),
            cache_dir=cache_dir,
            cache_mode=cache_mode,
            clear_cache=bool(args.clear_cache),
            sleep=float(args.sleep),
            rate_limiter=rate_limiter,
            max_workers=int(args.max_workers),
            max_retries=int(args.max_retries),
            backoff_base=float(args.backoff_base),
            backoff_max=float(args.backoff_max),
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
