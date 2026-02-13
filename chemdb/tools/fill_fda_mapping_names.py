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
UNII_TOKEN_RE = re.compile(r"\b[A-Z0-9]{10}\b")
UNII_EXACT_RE = re.compile(r"^[A-Z0-9]{10}$")
IDENTIFIER_LIKE_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^(?:S?CHEMBL)\d+$", re.IGNORECASE),
    re.compile(r"^CAS[-\s]?\d{2,7}-\d{2}-\d$", re.IGNORECASE),
    re.compile(r"^NSC\d+$", re.IGNORECASE),
    re.compile(r"^ORB\d+$", re.IGNORECASE),
    re.compile(r"^CID\d+$", re.IGNORECASE),
    re.compile(r"^CHEBI:\d+$", re.IGNORECASE),
    re.compile(r"^DB\d+$", re.IGNORECASE),
)
INCHIKEY_EXACT_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", re.IGNORECASE)
CAS_ONLY_RE = re.compile(r"^(?:CAS[-\s:]*)?\d{2,7}-\d{2}-\d$", re.IGNORECASE)
CAS_VALUE_RE = re.compile(r"(?:CAS[-\s:]*)?(\d{2,7}-\d{2}-\d)", re.IGNORECASE)
FORMULA_TOKEN_RE = re.compile(r"[A-Z][a-z]?\d*")
SPACED_IUPAC_LOCANT_TOKEN_RE = re.compile(r"^\d+[a-z]?$", re.IGNORECASE)
SPACED_IUPAC_TERM_RE = re.compile(
    r"(hydroxy|phenyl|methyl|ethyl|propyl|butyl|oxo|amino|carboxy|inden|imidazol|benzyl|acetyl|sulfon|amine|amide|acid|one|ol)$",
    re.IGNORECASE,
)
NAME_HINT_RE = re.compile(r"^[A-Za-z][A-Za-z0-9'`\-\s]{1,80}$")
UNICHEM_CHEMBL_SOURCE_ID = 1
DEFAULT_GSRS_BASE_URL = "https://gsrs.ncats.nih.gov/ginas/app/api/v1"
DEFAULT_DRUGCENTRAL_STRUCTURES_TSV = (
    Path(__file__).resolve().parents[2] / "drugcentral" / "structures.smiles.tsv"
)


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
        from rdkit import Chem as _Chem  # type: ignore[import-untyped]

        _ = _Chem
    except Exception:
        _RDKit_AVAILABLE = False
        return False
    _RDKit_AVAILABLE = True
    return True


def inchikey_to_2d(key: str) -> str:
    text = _clean_text(key).upper()
    if not text:
        return ""
    block = text.split("-", 1)[0] if "-" in text else text[:14]
    if len(block) != 14 or not block.isalnum():
        return ""
    return block


def canonicalize_for_key(mol: Any) -> Any:
    if mol is None:
        return None
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore[import-untyped]
    except Exception:
        return mol

    normalized = mol
    try:
        normalized = rdMolStandardize.Cleanup(normalized)
    except Exception:
        pass
    try:
        normalized = rdMolStandardize.FragmentParent(normalized)
    except Exception:
        pass
    try:
        normalized = rdMolStandardize.ChargeParent(normalized)
    except Exception:
        pass
    try:
        if hasattr(normalized, "GetNumAtoms") and normalized.GetNumAtoms() > 120:
            return normalized
        enumerator = rdMolStandardize.TautomerEnumerator()
        if hasattr(enumerator, "SetMaxTautomers"):
            enumerator.SetMaxTautomers(128)
        if hasattr(enumerator, "SetMaxTransforms"):
            enumerator.SetMaxTransforms(256)
        normalized = enumerator.Canonicalize(normalized)
    except Exception:
        pass
    return normalized if normalized is not None else mol


def _rdkit_parent_mol(mol: Any) -> Any:
    return canonicalize_for_key(mol)


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
    mol = _rdkit_parent_mol(mol)
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
    mol = _rdkit_parent_mol(mol)
    try:
        return _clean_text(rdkit_inchi.MolToInchiKey(mol))
    except Exception:
        return ""


def rdkit_inchi_from_smiles(smiles_value: str) -> str:
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
    mol = _rdkit_parent_mol(mol)
    try:
        return _clean_text(rdkit_inchi.MolToInchi(mol))
    except Exception:
        return ""


@dataclass
class DrugCentralIndex:
    full: dict[str, str] = field(default_factory=dict)
    two_d: dict[str, str] = field(default_factory=dict)


def load_drugcentral_index(path: Path) -> DrugCentralIndex:
    if not path.exists():
        raise ValueError(f"DrugCentral TSV not found: {path}")
    full_index: dict[str, str] = {}
    two_d_index: dict[str, str] = {}

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
        )
        inchikey_2d_col = _choose_header(
            headers,
            "inchikey_2d",
            "inchi_key_2d",
            "inchikey2d",
            "inchikeyconnectivity",
            "inchi_key2d",
        )
        inchi_col = _choose_header(headers, "inchi", "inchi_string")
        smiles_col = _choose_header(headers, "smiles", "canonical_smiles")
        rdkit_available = _rdkit_available()
        if not inchikey_col and not inchikey_2d_col and not rdkit_available:
            LOGGER.warning(
                "DrugCentral TSV missing INCHIKEY column; RDKit not available so cannot "
                "compute InChIKeys. Install RDKit or provide a TSV with INCHIKEY. "
                "DrugCentral stage will have low coverage."
            )
        LOGGER.info(
            "[drugcentral] columns detected: inn=%s name=%s drug_name=%s "
            "inchikey=%s inchikey_2d=%s inchi=%s smiles=%s",
            inn_col or "-",
            name_col or "-",
            drug_name_col or "-",
            inchikey_col or "-",
            inchikey_2d_col or "-",
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

            full_key = _clean_text(row.get(inchikey_col)).upper() if inchikey_col else ""
            two_d_key = (
                _clean_text(row.get(inchikey_2d_col)).upper() if inchikey_2d_col else ""
            )
            inchi_value = _clean_text(row.get(inchi_col)) if inchi_col else ""
            smiles_value = _clean_text(row.get(smiles_col)) if smiles_col else ""

            full_candidates: list[str] = []
            if full_key:
                full_candidates.append(full_key)

            if rdkit_available and not full_key:
                if inchi_value:
                    derived_inchi_key = rdkit_inchikey_from_inchi(inchi_value).upper()
                    if derived_inchi_key:
                        full_candidates.append(derived_inchi_key)
                if smiles_value:
                    derived_smiles_key = rdkit_inchikey_from_smiles(smiles_value).upper()
                    if derived_smiles_key:
                        full_candidates.append(derived_smiles_key)

            full_candidates = _dedupe_preserve(full_candidates)
            if not two_d_key and full_candidates:
                two_d_key = inchikey_to_2d(full_candidates[0])

            for candidate in full_candidates:
                if candidate and candidate not in full_index:
                    full_index[candidate] = name
                candidate_2d = inchikey_to_2d(candidate)
                if candidate_2d and candidate_2d not in two_d_index:
                    two_d_index[candidate_2d] = name
            if two_d_key and two_d_key not in two_d_index:
                two_d_index[two_d_key] = name
            if not full_candidates and not two_d_key:
                continue

    sample_full_len = len(next(iter(full_index))) if full_index else 0
    sample_2d_len = len(next(iter(two_d_index))) if two_d_index else 0
    LOGGER.info(
        "[drugcentral] indexed full=%s two_d=%s sample_key_lengths(full=%s,two_d=%s)",
        len(full_index),
        len(two_d_index),
        sample_full_len,
        sample_2d_len,
    )
    return DrugCentralIndex(full=full_index, two_d=two_d_index)


def is_inchikey_string(name: str) -> bool:
    text = _collapse_spaces(name).upper()
    if not text:
        return False
    if INCHIKEY_EXACT_RE.match(text):
        return True
    tokens = text.split()
    if len(tokens) != 3:
        return False
    lengths = [14, 10, 1]
    return all(
        len(token) == expected and token.isalpha()
        for token, expected in zip(tokens, lengths, strict=True)
    )


def is_cas_only_name(name: str) -> bool:
    text = _collapse_spaces(name)
    if not text:
        return False
    return CAS_ONLY_RE.match(text) is not None


def is_molecular_formula_like(name: str) -> bool:
    text = _clean_text(name)
    if not text or " " in text or "-" in text:
        return False
    if not text[:1].isupper():
        return False
    parts = text.split(".")
    if not parts or any(not part for part in parts):
        return False

    total_tokens = 0
    for part in parts:
        tokens = FORMULA_TOKEN_RE.findall(part)
        if not tokens or "".join(tokens) != part:
            return False
        total_tokens += len(tokens)
        if not any(any(ch.isdigit() for ch in token) for token in tokens):
            return False
    return total_tokens >= 2


def is_formula_like_name(name: str) -> bool:
    return is_molecular_formula_like(name)


def is_unii_like_name(name: str) -> bool:
    text = _clean_text(name).upper()
    if not text or " " in text or "-" in text:
        return False
    if UNII_EXACT_RE.match(text) is None:
        return False
    return any(ch.isdigit() for ch in text)


def is_spaced_iupacish_name(name: str) -> bool:
    text = _collapse_spaces(name).lower()
    if not text:
        return False
    tokens = text.split()
    if len(tokens) < 4:
        return False
    locants = [
        token for token in tokens if SPACED_IUPAC_LOCANT_TOKEN_RE.match(token) is not None
    ]
    if len(locants) < 2:
        return False
    if any(SPACED_IUPAC_TERM_RE.search(token) for token in tokens):
        return True
    return False


def is_identifier_like_name(name: str) -> bool:
    text = _clean_text(name)
    if not text:
        return False
    if is_inchikey_string(text):
        return True
    if is_unii_like_name(text):
        return True
    return any(pattern.match(text) is not None for pattern in IDENTIFIER_LIKE_PATTERNS)


def display_name_quality_reasons(
    name: str,
    pubchem_iupac_name: str | None = None,
    *,
    strict: bool = False,
) -> list[str]:
    text = _clean_text(name)
    reasons: list[str] = []
    if not text:
        reasons.append("empty")
        return reasons
    if FDA_PLACEHOLDER_RE.match(text):
        reasons.append("fda_placeholder")
    if UNK_PLACEHOLDER_RE.match(text):
        reasons.append("unk_placeholder")
    if is_inchikey_string(text):
        reasons.append("inchikey_like")
    if is_identifier_like_name(text):
        reasons.append("identifier_like")
    if is_cas_only_name(text):
        reasons.append("cas_only")
    if is_molecular_formula_like(text):
        reasons.append("formula_like")
    if is_unii_like_name(text):
        reasons.append("unii_like")

    if pubchem_iupac_name and text.lower() == _clean_text(pubchem_iupac_name).lower():
        reasons.append("equals_pubchem_iupac")
    if (
        IUPAC_PAREN_DIGIT_RE.search(text)
        or IUPAC_DIGIT_LOCANT_RE.search(text)
        or IUPAC_NOS_LOCANT_RE.search(text)
        or IUPAC_COMMA_LOCANT_RE.search(text)
        or IUPAC_TOKEN_RE.search(text)
        or (len(text) > 60 and IUPAC_LONG_PUNCT_RE.search(text))
    ):
        reasons.append("iupac_like")

    if strict and is_spaced_iupacish_name(text):
        reasons.append("spaced_iupacish")
    return _dedupe_preserve(reasons)


def is_bad_display_name(name: str, pubchem_iupac_name: str | None = None) -> bool:
    return bool(display_name_quality_reasons(name, pubchem_iupac_name))


def is_low_quality_display_name(
    name: str,
    pubchem_iupac_name: str | None = None,
    *,
    strict: bool = False,
) -> bool:
    return bool(
        display_name_quality_reasons(name, pubchem_iupac_name, strict=strict)
    )


def is_good_display_name(name: str, pubchem_iupac_name: str | None = None) -> bool:
    return not is_bad_display_name(name, pubchem_iupac_name)


@dataclass(frozen=True)
class PubChemResult:
    cid: str | None = None
    title: str | None = None
    iupac_name: str | None = None
    synonyms: list[str] = field(default_factory=list)
    smiles_400_inchi_fallback_attempted: bool = False
    inchi_cid_resolved: bool = False


@dataclass(frozen=True)
class PubChemLastMileResult:
    name: str = ""
    record_title: str = ""
    parent_cid: str = ""
    pugview_calls: int = 0
    pugview_success: int = 0
    pugview_fail: int = 0
    pugview_title_resolved: bool = False
    parent_title_resolved: bool = False


@dataclass(frozen=True)
class PubChemCasLastMileResult:
    name: str = ""
    source_cid: str = ""
    record_title: str = ""
    candidates_considered: int = 0
    title_resolved: bool = False
    pugview_calls: int = 0
    pugview_success: int = 0
    pugview_fail: int = 0


@dataclass(frozen=True)
class GSRSResolution:
    name: str = ""
    mode: str = ""


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

    def _request_json_with_status(
        self,
        url: str,
        *,
        method: str = "get",
        json_payload: Optional[dict[str, Any]] = None,
    ) -> tuple[Optional[Any], Optional[int]]:
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
                    return None, None
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
                    return None, status
                if isinstance(payload, (dict, list)):
                    return payload, status
                LOGGER.warning("[http] status=200 url=%s reason=non_object_json", url)
                return None, status

            if status in (429, 503):
                if retries >= self.max_retries:
                    LOGGER.warning(
                        "[http] status=%s url=%s reason=max_retries", status, url
                    )
                    return None, status
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
            return None, status

    def _request_json(
        self,
        url: str,
        *,
        method: str = "get",
        json_payload: Optional[dict[str, Any]] = None,
    ) -> Optional[Any]:
        payload, _ = self._request_json_with_status(
            url,
            method=method,
            json_payload=json_payload,
        )
        return payload

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

    def get_json_with_status(self, url: str) -> tuple[Optional[Any], Optional[int]]:
        cache_path = self._cache_path(url) if self.cache_mode != "off" else None
        if self.cache_mode == "use" and cache_path is not None:
            cached_payload = self._load_cache(cache_path)
            if cached_payload is not None:
                return cached_payload, 200

        payload, status = self._request_json_with_status(url)
        if payload is None:
            return None, status
        if self.cache_mode != "off" and cache_path is not None:
            self._write_cache(
                cache_path,
                url,
                payload,
                overwrite=self.cache_mode == "refresh",
            )
        return payload, status

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


class GSRSClient:
    def __init__(self, base_url: str, client: CachedJsonClient) -> None:
        self.base_url = _clean_text(base_url).rstrip("/")
        self.client = client

    def _get(self, path: str, params: Optional[dict[str, Any]] = None) -> Optional[Any]:
        if not self.base_url:
            return None
        url = f"{self.base_url}/{path.lstrip('/')}"
        if params:
            url = f"{url}?{urlencode(params)}"
        return self.client.get_json(url)

    def _records(self, payload: Optional[Any]) -> list[dict[str, Any]]:
        if payload is None:
            return []
        if isinstance(payload, list):
            return [item for item in payload if isinstance(item, dict)]
        if isinstance(payload, dict):
            if "uuid" in payload:
                return [payload]
            for key in (
                "content",
                "results",
                "data",
                "substances",
                "items",
                "matches",
            ):
                value = payload.get(key)
                if isinstance(value, list):
                    return [item for item in value if isinstance(item, dict)]
        return []

    def _search(self, query: str) -> list[dict[str, Any]]:
        text = _clean_text(query)
        if not text:
            return []
        attempts = [
            ("substances/search", {"q": text}),
            ("substances/search", {"query": text}),
            ("substances", {"q": text}),
            ("substances", {"query": text}),
        ]
        for path, params in attempts:
            payload = self._get(path, params)
            records = self._records(payload)
            if records:
                return records
        return []

    def search_by_unii(self, unii: str) -> list[dict[str, Any]]:
        text = _clean_text(unii).upper()
        if not text:
            return []
        queries = (
            f"approvalID:{text}",
            f"root_approvalID:{text}",
            f"UNII:{text}",
            text,
        )
        for query in queries:
            records = self._search(query)
            if records:
                return records
        return []

    def search_by_inchikey(self, key: str) -> list[dict[str, Any]]:
        text = _clean_text(key).upper()
        if not text:
            return []
        queries = (
            f"inchikey:{text}",
            f"structure.inchikey:{text}",
            text,
        )
        for query in queries:
            records = self._search(query)
            if records:
                return records
        return []

    def search_by_structure(self, structure: str) -> list[dict[str, Any]]:
        text = _clean_text(structure)
        if not text:
            return []
        attempts = [
            ("substances/structureSearch", {"q": text, "type": "exact"}),
            ("substances/structureSearch", {"query": text, "type": "exact"}),
            ("substances/structureSearch", {"smiles": text, "matchType": "exact"}),
            ("substances/structureSearch", {"q": text}),
        ]
        for path, params in attempts:
            payload = self._get(path, params)
            records = self._records(payload)
            if records:
                return records
        return []


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
        "inchi",
        "remark_inchi",
    ):
        value = _clean_text(row.get(key))
        if value:
            return key, value
    return None, None


def _resolution_keys(row: dict[str, Any]) -> list[tuple[str, str]]:
    keys: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for key in (
        "inchikey",
        "remark_inchikey",
        "smiles_neutral",
        "smiles",
        "remark_smiles",
        "inchi",
        "remark_inchi",
    ):
        value = _clean_text(row.get(key))
        if not value:
            continue
        pair = (key, value)
        if pair in seen:
            continue
        seen.add(pair)
        keys.append(pair)
    return keys


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
    try:
        cid = int(cids[0])
    except Exception:
        return None
    if cid <= 0:
        return None
    return str(cid)


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


def _parse_pubchem_cids(payload: Optional[dict[str, Any]]) -> list[str]:
    if not payload:
        return []
    try:
        raw_cids = payload["IdentifierList"]["CID"]
    except Exception:
        return []
    if not isinstance(raw_cids, list):
        return []
    cids: list[str] = []
    for value in raw_cids:
        try:
            cid = int(value)
        except Exception:
            continue
        if cid > 0:
            cids.append(str(cid))
    return cids


def _extract_cas_tokens(value: str) -> list[str]:
    text = _clean_text(value)
    if not text:
        return []
    matches = [match.group(1) for match in CAS_VALUE_RE.finditer(text)]
    if matches:
        return _dedupe_preserve(matches)
    tokens = [
        token.strip()
        for token in re.split(r"[;|,\s]+", text)
        if token.strip()
    ]
    normalized: list[str] = []
    for token in tokens:
        match = CAS_VALUE_RE.fullmatch(token)
        if match:
            normalized.append(match.group(1))
    return _dedupe_preserve(normalized)


def _collect_cas_from_row(row: dict[str, Any]) -> list[str]:
    cas_values: list[str] = []
    for field_name in ("cas", "cas_number", "cas_no", "casrn", "remark_cas"):
        cas_values.extend(_extract_cas_tokens(_clean_text(row.get(field_name))))
    return _dedupe_preserve(cas_values)


def pubchem_cids_from_name(name_str: str, client: CachedJsonClient) -> list[int]:
    query = _clean_text(name_str)
    if not query:
        return []
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{quote(query, safe='')}/cids/JSON"
    )
    payload = client.get_json(url)
    raw_cids = _parse_pubchem_cids(payload if isinstance(payload, dict) else None)
    cids: list[int] = []
    seen: set[int] = set()
    for value in raw_cids:
        try:
            cid = int(value)
        except Exception:
            continue
        if cid <= 0 or cid in seen:
            continue
        seen.add(cid)
        cids.append(cid)
    return cids


def _pubchem_cid_from_inchi(inchi_value: str, client: CachedJsonClient) -> Optional[str]:
    text = _clean_text(inchi_value)
    if not text:
        return None
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/"
        f"{quote(text, safe='')}/cids/JSON"
    )
    payload = client.get_json(url)
    return _parse_pubchem_cid(payload if isinstance(payload, dict) else None)


def _pubchem_fastidentity_cids_from_inchikey(
    inchikey: str, client: CachedJsonClient
) -> list[str]:
    key = _clean_text(inchikey).upper()
    if not key:
        return []
    base = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/inchikey/"
        f"{quote(key, safe='')}/cids/JSON?identity_type="
    )
    for identity_type in ("same_parent_connectivity", "same_connectivity"):
        payload = client.get_json(base + identity_type)
        cids = _parse_pubchem_cids(payload if isinstance(payload, dict) else None)
        if cids:
            return cids
    return []


def _pubchem_cid_from_inchikey(inchikey: str, client: CachedJsonClient) -> Optional[str]:
    key = _clean_text(inchikey).upper()
    if not key:
        return None
    candidate_keys: list[str] = [key]
    key_2d = inchikey_to_2d(key)
    if key_2d and key_2d not in candidate_keys:
        candidate_keys.append(key_2d)

    for candidate in candidate_keys:
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
            f"{quote(candidate, safe='')}/cids/JSON"
        )
        payload = client.get_json(url)
        cid = _parse_pubchem_cid(payload if isinstance(payload, dict) else None)
        if cid:
            return cid

    for candidate in candidate_keys:
        cids = _pubchem_fastidentity_cids_from_inchikey(candidate, client)
        if cids:
            return cids[0]
    return None


def _iter_pugview_sections(section: dict[str, Any]) -> Iterable[dict[str, Any]]:
    yield section
    children = section.get("Section", [])
    if not isinstance(children, list):
        return
    for child in children:
        if isinstance(child, dict):
            yield from _iter_pugview_sections(child)


def _extract_pugview_strings(section: dict[str, Any]) -> list[str]:
    values: list[str] = []
    info_list = section.get("Information", [])
    if not isinstance(info_list, list):
        return values
    for info in info_list:
        if not isinstance(info, dict):
            continue
        value = info.get("Value", {})
        if not isinstance(value, dict):
            continue
        strings = value.get("StringWithMarkup", [])
        if not isinstance(strings, list):
            continue
        for entry in strings:
            if not isinstance(entry, dict):
                continue
            text = _clean_text(entry.get("String"))
            if text:
                values.append(text)
    return values


def _parse_pubchem_pugview_record_title(payload: Optional[dict[str, Any]]) -> str:
    if not payload or not isinstance(payload, dict):
        return ""
    try:
        record = payload.get("Record", {})
    except Exception:
        return ""
    if not isinstance(record, dict):
        return ""
    return _clean_text(record.get("RecordTitle"))


def _pugview_root_heading_allowed(heading: str) -> bool:
    text = _clean_text(heading).lower()
    if not text:
        return False
    return any(
        token in text
        for token in ("name", "identifier", "drug", "medication", "pharmacology")
    )


def _pugview_section_heading_allowed(section_heading: str, root_heading: str) -> bool:
    section_text = _clean_text(section_heading).lower()
    root_text = _clean_text(root_heading).lower()
    combined = f"{root_text} {section_text}".strip()
    if "iupac" in section_text:
        return False
    return any(
        token in combined
        for token in (
            "name",
            "title",
            "synonym",
            "drug",
            "medication",
            "drugbank",
            "fda",
            "inn",
            "usan",
            "generic",
            "international nonproprietary",
            "common name",
        )
    )


def _parse_pubchem_pugview_name_candidates_with_sources(
    payload: Optional[dict[str, Any]],
) -> list[tuple[str, str]]:
    if not payload or not isinstance(payload, dict):
        return []
    try:
        record = payload.get("Record", {})
        root_sections = record.get("Section", []) if isinstance(record, dict) else []
    except Exception:
        return []
    if not isinstance(root_sections, list):
        return []

    candidates: list[tuple[str, str]] = []
    for root in root_sections:
        if not isinstance(root, dict):
            continue
        root_heading = _clean_text(root.get("TOCHeading"))
        if not _pugview_root_heading_allowed(root_heading):
            continue
        for section in _iter_pugview_sections(root):
            section_heading = _clean_text(section.get("TOCHeading"))
            if not _pugview_section_heading_allowed(section_heading, root_heading):
                continue
            source_heading = section_heading or root_heading
            for value in _extract_pugview_strings(section):
                candidates.append((value, source_heading))

    deduped: list[tuple[str, str]] = []
    seen: set[str] = set()
    for value, source_heading in candidates:
        text = _clean_text(value)
        if not text:
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        deduped.append((text, _clean_text(source_heading)))
    return deduped


def _parse_pubchem_pugview_name_candidates(payload: Optional[dict[str, Any]]) -> list[str]:
    return [
        candidate
        for candidate, _source in _parse_pubchem_pugview_name_candidates_with_sources(
            payload
        )
    ]


def _parse_pubchem_pugview_unii(payload: Optional[dict[str, Any]]) -> str:
    if not payload or not isinstance(payload, dict):
        return ""
    try:
        record = payload.get("Record", {})
        root_sections = record.get("Section", []) if isinstance(record, dict) else []
    except Exception:
        return ""
    if not isinstance(root_sections, list):
        return ""
    for root in root_sections:
        if not isinstance(root, dict):
            continue
        for section in _iter_pugview_sections(root):
            heading = _clean_text(section.get("TOCHeading")).lower()
            if "unii" not in heading:
                continue
            for value in _extract_pugview_strings(section):
                cleaned = _clean_text(value).upper()
                if cleaned:
                    return cleaned
    return ""


def pubchem_pugview_compound(cid: int, client: CachedJsonClient) -> Optional[dict[str, Any]]:
    if cid <= 0:
        return None
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
        f"{cid}/JSON/?response_type=display"
    )
    payload = client.get_json(url)
    return payload if isinstance(payload, dict) else None


def pubchem_fastidentity_parent_cids(cid: int, client: CachedJsonClient) -> list[int]:
    if cid <= 0:
        return []
    base = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/cid/"
        f"{cid}/cids/JSON?identity_type="
    )
    candidate_urls = [
        base + "same_parent_connectivity",
        base + "same_connectivity",
    ]
    for url in candidate_urls:
        payload = client.get_json(url)
        cids = _parse_pubchem_cids(payload if isinstance(payload, dict) else None)
        parsed: list[int] = []
        for item in cids:
            try:
                parsed_cid = int(item)
            except Exception:
                continue
            if parsed_cid > 0:
                parsed.append(parsed_cid)
        if parsed:
            return parsed
    return []


def _resolve_pubchem_last_mile(
    *,
    row: dict[str, Any],
    pubchem_cid: str,
    pubchem_iupac_name: str,
    pubchem_client: CachedJsonClient,
) -> PubChemLastMileResult:
    try:
        cid_int = int(pubchem_cid)
    except Exception:
        return PubChemLastMileResult()

    calls = 0
    success = 0
    fail = 0

    payload = pubchem_pugview_compound(cid_int, pubchem_client)
    calls += 1
    if payload is None:
        fail += 1
    else:
        success += 1
    record_title = _parse_pubchem_pugview_record_title(payload)
    if record_title and is_good_display_name(record_title, pubchem_iupac_name):
        _maybe_set_field_overwrite_if_bad(
            row,
            "pubchem_record_title",
            record_title,
            pubchem_iupac_name=pubchem_iupac_name,
        )
    unii = _parse_pubchem_pugview_unii(payload)
    if unii:
        _maybe_set_field(row, "pubchem_unii_list", unii)

    if record_title and is_good_display_name(record_title, pubchem_iupac_name):
        return PubChemLastMileResult(
            name=record_title,
            record_title=record_title,
            pugview_calls=calls,
            pugview_success=success,
            pugview_fail=fail,
            pugview_title_resolved=True,
        )

    candidate_sources = _parse_pubchem_pugview_name_candidates_with_sources(payload)
    source_headings = {name.lower(): source for name, source in candidate_sources}
    selected = _select_pubchem_candidate(
        [name for name, _source in candidate_sources],
        pubchem_iupac_name,
        source_headings=source_headings,
        rxnorm_generic_name=_clean_text(row.get("rxnorm_generic_name")),
    )
    if selected:
        _maybe_set_field_overwrite_if_bad(
            row,
            "pubchem_name",
            selected,
            pubchem_iupac_name=pubchem_iupac_name,
        )
        return PubChemLastMileResult(
            name=selected,
            record_title=record_title,
            pugview_calls=calls,
            pugview_success=success,
            pugview_fail=fail,
            pugview_title_resolved=True,
        )

    parent_cids = pubchem_fastidentity_parent_cids(cid_int, pubchem_client)
    for parent_cid in parent_cids:
        if parent_cid <= 0 or parent_cid == cid_int:
            continue
        parent_payload = pubchem_pugview_compound(parent_cid, pubchem_client)
        calls += 1
        if parent_payload is None:
            fail += 1
            continue
        success += 1
        parent_title = _parse_pubchem_pugview_record_title(parent_payload)
        if parent_title and is_good_display_name(parent_title, pubchem_iupac_name):
            _maybe_set_field_overwrite_if_bad(
                row,
                "pubchem_record_title",
                parent_title,
                pubchem_iupac_name=pubchem_iupac_name,
            )
            return PubChemLastMileResult(
                name=parent_title,
                record_title=record_title,
                parent_cid=str(parent_cid),
                pugview_calls=calls,
                pugview_success=success,
                pugview_fail=fail,
                parent_title_resolved=True,
            )

        parent_candidates = _parse_pubchem_pugview_name_candidates_with_sources(
            parent_payload
        )
        parent_source_headings = {
            name.lower(): source for name, source in parent_candidates
        }
        parent_selected = _select_pubchem_candidate(
            [name for name, _source in parent_candidates],
            pubchem_iupac_name,
            source_headings=parent_source_headings,
            rxnorm_generic_name=_clean_text(row.get("rxnorm_generic_name")),
        )
        if parent_selected:
            _maybe_set_field_overwrite_if_bad(
                row,
                "pubchem_name",
                parent_selected,
                pubchem_iupac_name=pubchem_iupac_name,
            )
            return PubChemLastMileResult(
                name=parent_selected,
                record_title=parent_title,
                parent_cid=str(parent_cid),
                pugview_calls=calls,
                pugview_success=success,
                pugview_fail=fail,
                parent_title_resolved=True,
            )

    return PubChemLastMileResult(
        pugview_calls=calls,
        pugview_success=success,
        pugview_fail=fail,
    )


def _resolve_pubchem_cas_last_mile(
    *,
    row: dict[str, Any],
    pubchem_iupac_name: str,
    pubchem_client: CachedJsonClient,
) -> PubChemCasLastMileResult:
    cas_values = _collect_cas_from_row(row)
    if not cas_values:
        return PubChemCasLastMileResult()

    calls = 0
    success = 0
    fail = 0
    candidates: list[str] = []
    source_headings: dict[str, str] = {}
    source_cid_by_name: dict[str, str] = {}
    record_title_by_name: dict[str, str] = {}
    seen_candidates: set[str] = set()
    cid_queue: list[int] = []
    seen_cids: set[int] = set()

    def enqueue_cid(value: int) -> None:
        if value <= 0 or value in seen_cids:
            return
        seen_cids.add(value)
        cid_queue.append(value)

    def add_candidate(
        value: str,
        source_heading: str,
        source_cid: int,
        record_title: str,
    ) -> None:
        text = _clean_text(value)
        if not text:
            return
        key = text.lower()
        if key in seen_candidates:
            return
        seen_candidates.add(key)
        candidates.append(text)
        source_headings[key] = _clean_text(source_heading)
        source_cid_by_name[key] = str(source_cid)
        record_title_by_name[key] = _clean_text(record_title)

    for cas_value in cas_values:
        for cid_int in pubchem_cids_from_name(cas_value, pubchem_client):
            enqueue_cid(cid_int)

    # PubChem may miss CAS-based lookups for some salts; probe neutral structure as fallback.
    if not cid_queue:
        neutral_smiles = _clean_text(row.get("smiles_neutral"))
        if neutral_smiles:
            neutral_probe = resolve_pubchem(
                {"smiles_neutral": neutral_smiles},
                pubchem_client,
                require_existing_cid=False,
            )
            if neutral_probe.cid:
                try:
                    enqueue_cid(int(neutral_probe.cid))
                except Exception:
                    pass

    for cid_int in cid_queue:
        payload = pubchem_pugview_compound(cid_int, pubchem_client)
        calls += 1
        if payload is None:
            fail += 1
            continue
        success += 1
        record_title = _parse_pubchem_pugview_record_title(payload)
        if record_title:
            add_candidate(record_title, "record_title", cid_int, record_title)

        unii = _parse_pubchem_pugview_unii(payload)
        if unii:
            _maybe_set_field(row, "pubchem_unii_list", unii)

        for value, source in _parse_pubchem_pugview_name_candidates_with_sources(payload):
            add_candidate(value, source or "pugview_name", cid_int, record_title)

        props_url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid_int}/property/Title,IUPACName/JSON"
        )
        title, _iupac_name = _parse_pubchem_properties(pubchem_client.get_json(props_url))
        if title:
            add_candidate(title, "property_title", cid_int, record_title)

        syn_url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid_int}/synonyms/JSON"
        )
        for synonym in _parse_pubchem_synonyms(pubchem_client.get_json(syn_url)):
            add_candidate(synonym, "cid_synonym", cid_int, record_title)

    selected = _select_pubchem_candidate(
        candidates,
        pubchem_iupac_name,
        source_headings=source_headings,
        rxnorm_generic_name=_clean_text(row.get("rxnorm_generic_name")),
    )
    if not selected:
        return PubChemCasLastMileResult(
            candidates_considered=len(candidates),
            pugview_calls=calls,
            pugview_success=success,
            pugview_fail=fail,
        )

    selected_key = selected.lower()
    selected_record_title = _clean_text(record_title_by_name.get(selected_key))
    source_cid = _clean_text(source_cid_by_name.get(selected_key))
    if _row_has_besylate_counterion(row) and not selected.lower().endswith("besylate"):
        selected = f"{selected} besylate"
    return PubChemCasLastMileResult(
        name=selected,
        source_cid=source_cid,
        record_title=selected_record_title,
        candidates_considered=len(candidates),
        title_resolved=selected_record_title.lower() == _clean_text(selected).lower()
        if selected_record_title
        else False,
        pugview_calls=calls,
        pugview_success=success,
        pugview_fail=fail,
    )


def _row_inchikey(row: dict[str, Any]) -> str:
    return _clean_text(_pick_first(row.get("inchikey"), row.get("remark_inchikey"))).upper()


def _drugcentral_lookup_candidates(row: dict[str, Any]) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []

    full = _row_inchikey(row)
    if full:
        candidates.append((full, "row_inchikey_full"))
        full_2d = inchikey_to_2d(full)
        if full_2d:
            candidates.append((full_2d, "row_inchikey_2d"))

    if _rdkit_available():
        for field_name in ("inchi", "remark_inchi"):
            inchi_value = _clean_text(row.get(field_name))
            if not inchi_value:
                continue
            rdkit_key = rdkit_inchikey_from_inchi(inchi_value).upper()
            if not rdkit_key:
                continue
            candidates.append((rdkit_key, f"{field_name}_tautomer_full"))
            rdkit_2d = inchikey_to_2d(rdkit_key)
            if rdkit_2d:
                candidates.append((rdkit_2d, f"{field_name}_tautomer_2d"))

        for field_name in ("smiles_neutral", "smiles", "remark_smiles"):
            smiles_value = _clean_text(row.get(field_name))
            if not smiles_value:
                continue
            rdkit_key = rdkit_inchikey_from_smiles(smiles_value).upper()
            if not rdkit_key:
                continue
            candidates.append((rdkit_key, f"{field_name}_tautomer_full"))
            rdkit_2d = inchikey_to_2d(rdkit_key)
            if rdkit_2d:
                candidates.append((rdkit_2d, f"{field_name}_tautomer_2d"))

    deduped: list[tuple[str, str]] = []
    seen: set[str] = set()
    for key, source in candidates:
        cleaned = _clean_text(key).upper()
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        deduped.append((cleaned, source))
    return deduped


def _match_drugcentral_name(
    index: DrugCentralIndex,
    candidates: Iterable[tuple[str, str]],
) -> tuple[str, str, str, str]:
    for key, source in candidates:
        key_text = _clean_text(key).upper()
        if not key_text:
            continue
        if len(key_text) == 14:
            candidate_name = _clean_text(index.two_d.get(key_text))
            if candidate_name:
                return candidate_name, "2d", source, key_text
            continue

        candidate_name = _clean_text(index.full.get(key_text))
        if candidate_name:
            return candidate_name, "exact", source, key_text
        key_2d = inchikey_to_2d(key_text)
        if not key_2d:
            continue
        candidate_name = _clean_text(index.two_d.get(key_2d))
        if candidate_name:
            return candidate_name, "2d", source, key_2d
    return "", "", "", ""


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
    smiles_400_inchi_fallback_attempted = False
    inchi_cid_resolved = False
    if not cid:
        if require_existing_cid:
            return PubChemResult()
        resolution_keys = _resolution_keys(row)
        if not resolution_keys:
            return PubChemResult()
        for key_type, key_value in resolution_keys:
            if key_type in ("inchikey", "remark_inchikey"):
                cid = _pubchem_cid_from_inchikey(key_value, client)
                if cid:
                    break
                continue

            if key_type in ("inchi", "remark_inchi"):
                cid = _pubchem_cid_from_inchi(key_value, client)
                if not cid and _rdkit_available():
                    derived_inchikey = rdkit_inchikey_from_inchi(key_value)
                    if derived_inchikey:
                        cid = _pubchem_cid_from_inchikey(derived_inchikey, client)
                        if cid:
                            inchi_cid_resolved = True
                if cid:
                    break
                continue

            if key_type in ("smiles_neutral", "smiles", "remark_smiles"):
                url = (
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
                    f"{quote(key_value, safe='')}/cids/JSON"
                )
                cid_payload, status = client.get_json_with_status(url)
                cid = _parse_pubchem_cid(
                    cid_payload if isinstance(cid_payload, dict) else None
                )
                if cid:
                    break

                if status == 400:
                    inchi_value = _clean_text(
                        _pick_first(row.get("inchi"), row.get("remark_inchi"))
                    )
                    if not inchi_value and _rdkit_available():
                        inchi_value = rdkit_inchi_from_smiles(key_value)
                    if inchi_value:
                        smiles_400_inchi_fallback_attempted = True
                        cid = _pubchem_cid_from_inchi(inchi_value, client)
                        if cid:
                            inchi_cid_resolved = True
                            break
                        if _rdkit_available():
                            derived_inchikey = rdkit_inchikey_from_inchi(inchi_value)
                            if not derived_inchikey:
                                derived_inchikey = rdkit_inchikey_from_smiles(key_value)
                            if derived_inchikey:
                                cid = _pubchem_cid_from_inchikey(derived_inchikey, client)
                                if cid:
                                    inchi_cid_resolved = True
                                    break
                continue

            # Unknown resolution key type
            continue

        if not cid:
            return PubChemResult(
                smiles_400_inchi_fallback_attempted=smiles_400_inchi_fallback_attempted,
                inchi_cid_resolved=inchi_cid_resolved,
            )

    props_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/Title,IUPACName/JSON"
    )
    syn_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" f"{cid}/synonyms/JSON"
    )
    title, iupac_name = _parse_pubchem_properties(client.get_json(props_url))
    synonyms = _parse_pubchem_synonyms(client.get_json(syn_url))
    return PubChemResult(
        cid=cid,
        title=title,
        iupac_name=iupac_name,
        synonyms=synonyms,
        smiles_400_inchi_fallback_attempted=smiles_400_inchi_fallback_attempted,
        inchi_cid_resolved=inchi_cid_resolved,
    )


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


def _parse_openfda_substance_name(payload: Optional[dict[str, Any]]) -> str:
    if not payload or not isinstance(payload, dict):
        return ""
    results = payload.get("results")
    if not isinstance(results, list):
        return ""
    for item in results:
        if not isinstance(item, dict):
            continue
        name = _collapse_spaces(
            _pick_first(
                item.get("substance_name"),
                item.get("name"),
                item.get("preferred_name"),
            )
        )
        if name:
            return name
    return ""


def resolve_openfda_substance_name(unii: str, client: CachedJsonClient) -> str:
    token = _clean_text(unii).upper()
    if not token:
        return ""
    query = quote(f"unii:{token}", safe="")
    url = f"https://api.fda.gov/other/substance.json?search={query}&limit=1"
    payload = client.get_json(url)
    return _parse_openfda_substance_name(
        payload if isinstance(payload, dict) else None
    )


def _maybe_set_field(row: dict[str, Any], field: str, value: Optional[str]) -> bool:
    if field not in row:
        return False
    if _clean_text(row.get(field)):
        return False
    if not _clean_text(value):
        return False
    row[field] = value
    return True


def _maybe_set_field_overwrite_if_bad(
    row: dict[str, Any],
    field: str,
    value: Optional[str],
    *,
    pubchem_iupac_name: str | None = None,
    require_good_replacement: bool = True,
) -> bool:
    if field not in row:
        return False
    new_value = _clean_text(value)
    if not new_value:
        return False
    existing = _clean_text(row.get(field))
    if not existing:
        row[field] = value
        return True
    if not is_bad_display_name(existing, pubchem_iupac_name):
        return False
    if require_good_replacement and is_bad_display_name(new_value, pubchem_iupac_name):
        return False
    if existing == new_value:
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


def _looks_like_non_smiles_text(value: str) -> bool:
    text = _clean_text(value)
    if not text:
        return False
    lower = text.lower()
    if ";" in text:
        return True
    if any(token in lower for token in ("reference standard", "pharmacopoeia", "powder")):
        return True
    return False


def _extract_salvage_name_hint(row: dict[str, Any]) -> str:
    for field_name in ("smiles_neutral", "smiles", "remark_smiles"):
        raw = _clean_text(row.get(field_name))
        if not raw:
            continue
        if ";" in raw:
            hint = raw.rsplit(";", 1)[-1]
        elif _looks_like_non_smiles_text(raw):
            hint = raw
        else:
            continue
        hint = _collapse_spaces(hint.strip(" '\"\t\r\n.,"))
        if not hint:
            continue
        if not NAME_HINT_RE.match(hint):
            continue
        if is_identifier_like_name(hint):
            continue
        return hint
    return ""


def _row_has_besylate_counterion(row: dict[str, Any]) -> bool:
    fields = (
        "smiles",
        "remark_smiles",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_name",
    )
    for field_name in fields:
        text = _clean_text(row.get(field_name)).lower()
        if not text:
            continue
        if "benzenesulfonic acid" in text:
            return True
        if "o=s(=o)(o)c1ccccc1" in text:
            return True
    return False


def _is_trueish(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = _clean_text(value).lower()
    return text in {"1", "true", "yes", "y", "t"}


def _extract_unii_tokens(value: str) -> list[str]:
    text = _clean_text(value).upper()
    if not text:
        return []
    matches = [m.group(0) for m in UNII_TOKEN_RE.finditer(text)]
    if matches:
        return _dedupe_preserve(matches)
    tokens = [
        token.strip().upper()
        for token in re.split(r"[;|,\s]+", text)
        if token.strip()
    ]
    return [token for token in _dedupe_preserve(tokens) if len(token) == 10]


def _collect_uniis_from_row(row: dict[str, Any]) -> list[str]:
    field_names = (
        "pubchem_unii_list",
        "pubchem_unii",
        "unii",
        "unii_list",
        "gsrs_unii",
        "approval_id",
    )
    uniis: list[str] = []
    for field_name in field_names:
        value = _clean_text(row.get(field_name))
        if not value:
            continue
        uniis.extend(_extract_unii_tokens(value))
    return _dedupe_preserve(uniis)


def _rdkit_canonical_smiles(smiles_value: str) -> str:
    text = _clean_text(smiles_value)
    if not text:
        return ""
    if not _rdkit_available():
        return ""
    try:
        from rdkit import Chem  # type: ignore[import-untyped]
    except Exception:
        return ""
    try:
        mol = Chem.MolFromSmiles(text)
    except Exception:
        mol = None
    if mol is None:
        return ""
    mol = _rdkit_parent_mol(mol)
    try:
        return _clean_text(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    except Exception:
        return ""


def rdkit_smiles_from_inchi(inchi_value: str) -> str:
    text = _clean_text(inchi_value)
    if not text:
        return ""
    if not _rdkit_available():
        return ""
    try:
        from rdkit import Chem  # type: ignore[import-untyped]
    except Exception:
        return ""
    try:
        mol = Chem.MolFromInchi(text)
    except Exception:
        mol = None
    if mol is None:
        return ""
    mol = _rdkit_parent_mol(mol)
    try:
        return _clean_text(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    except Exception:
        return ""


def _structure_queries_from_row(row: dict[str, Any]) -> list[str]:
    if not _rdkit_available():
        return []
    queries: list[str] = []
    for field_name in ("smiles_neutral", "smiles", "remark_smiles"):
        canonical = _rdkit_canonical_smiles(_clean_text(row.get(field_name)))
        if canonical:
            queries.append(canonical)
    for field_name in ("inchi", "remark_inchi"):
        converted = rdkit_smiles_from_inchi(_clean_text(row.get(field_name)))
        if converted:
            queries.append(converted)
    return _dedupe_preserve(queries)


def _select_gsrs_name(
    records: Iterable[dict[str, Any]], pubchem_iupac_name: str | None
) -> str:
    preferred: list[str] = []
    general: list[str] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        for key in ("displayName", "preferredName", "name"):
            value = _collapse_spaces(record.get(key))
            if value:
                preferred.append(value)
                general.append(value)

        names = record.get("names")
        if not isinstance(names, list):
            continue
        for name_entry in names:
            if not isinstance(name_entry, dict):
                continue
            name_value = _collapse_spaces(
                _pick_first(
                    name_entry.get("displayName"),
                    name_entry.get("name"),
                    name_entry.get("stdName"),
                    name_entry.get("preferredName"),
                )
            )
            if not name_value:
                continue
            general.append(name_value)
            if _is_trueish(name_entry.get("preferred")):
                preferred.append(name_value)
                continue
            if _is_trueish(name_entry.get("displayName")):
                preferred.append(name_value)
                continue
            type_text = _clean_text(name_entry.get("type")).lower()
            if type_text in {"preferred", "display", "primary"}:
                preferred.append(name_value)

    for candidate in _dedupe_preserve(preferred):
        if is_good_display_name(candidate, pubchem_iupac_name):
            return candidate

    good_general = [
        candidate
        for candidate in _dedupe_preserve(general)
        if is_good_display_name(candidate, pubchem_iupac_name)
    ]
    if not good_general:
        return ""
    return min(good_general, key=lambda item: (len(item), item.lower()))


def resolve_gsrs_last_mile(
    row: dict[str, Any], pubchem_iupac_name: str, gsrs_client: GSRSClient
) -> GSRSResolution:
    for unii in _collect_uniis_from_row(row):
        records = gsrs_client.search_by_unii(unii)
        candidate = _select_gsrs_name(records, pubchem_iupac_name)
        if candidate:
            return GSRSResolution(name=candidate, mode="unii")

    full_inchikey = _row_inchikey(row)
    if full_inchikey:
        records = gsrs_client.search_by_inchikey(full_inchikey)
        candidate = _select_gsrs_name(records, pubchem_iupac_name)
        if candidate:
            return GSRSResolution(name=candidate, mode="inchikey_full")

        inchikey_2d = inchikey_to_2d(full_inchikey)
        if inchikey_2d and inchikey_2d != full_inchikey:
            records = gsrs_client.search_by_inchikey(inchikey_2d)
            candidate = _select_gsrs_name(records, pubchem_iupac_name)
            if candidate:
                return GSRSResolution(name=candidate, mode="inchikey_2d")

    for structure_query in _structure_queries_from_row(row):
        records = gsrs_client.search_by_structure(structure_query)
        candidate = _select_gsrs_name(records, pubchem_iupac_name)
        if candidate:
            return GSRSResolution(name=candidate, mode="structure")

    return GSRSResolution()


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


def _pubchem_heading_weight(source_heading: str) -> int:
    heading = _clean_text(source_heading).lower()
    if not heading:
        return 0
    if heading in {"record_title", "record title"}:
        return 6
    if heading == "property_title":
        return 5
    if "removed synonym" in heading:
        return -3
    if (
        "international nonproprietary name" in heading
        or "generic name" in heading
        or heading == "inn"
        or heading == "usan"
    ):
        return 6
    if "drug and medication information" in heading:
        return 5
    if "drugbank" in heading or "fda" in heading:
        return 5
    if "common name" in heading:
        return 4
    if "drug" in heading or "medication" in heading:
        return 4
    if "name" in heading or "synonym" in heading or "title" in heading:
        return 2
    return 0


def _score_pubchem_candidate(
    candidate: str,
    pubchem_iupac_name: str | None,
    *,
    source_heading: str = "",
    rxnorm_generic_name: str = "",
) -> Optional[int]:
    if is_identifier_like_name(candidate):
        return None
    if is_cas_only_name(candidate):
        return None
    if is_molecular_formula_like(candidate):
        return None
    if is_unii_like_name(candidate):
        return None
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
    if is_spaced_iupacish_name(text):
        score -= 8

    score += _pubchem_heading_weight(source_heading)

    rxnorm_name = _clean_text(rxnorm_generic_name).lower()
    if rxnorm_name:
        candidate_lower = text.lower()
        if candidate_lower == rxnorm_name:
            score += 4
        elif candidate_lower in rxnorm_name or rxnorm_name in candidate_lower:
            score += 2

    return score


def _select_pubchem_candidate(
    candidates: Iterable[str],
    pubchem_iupac_name: str | None,
    *,
    source_headings: Optional[dict[str, str]] = None,
    rxnorm_generic_name: str = "",
) -> str:
    best: tuple[int, int, str] | None = None
    best_name = ""
    for candidate in candidates:
        text = _clean_text(candidate)
        if not text:
            continue
        source_heading = ""
        if source_headings:
            source_heading = _clean_text(source_headings.get(text.lower()))
        score = _score_pubchem_candidate(
            text,
            pubchem_iupac_name,
            source_heading=source_heading,
            rxnorm_generic_name=rxnorm_generic_name,
        )
        if score is None:
            continue
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
    return _select_pubchem_candidate(
        _dedupe_preserve(candidates),
        pubchem_iupac,
        rxnorm_generic_name=_clean_text(row.get("rxnorm_generic_name")),
    )


@dataclass
class Summary:
    total_rows: int = 0
    rows_evaluated: int = 0
    rows_updated: int = 0
    drugcentral_resolved: int = 0
    drugcentral_resolved_exact: int = 0
    drugcentral_resolved_2d: int = 0
    chembl_resolved: int = 0
    chembl_approved_by_max_phase: int = 0
    chembl_approved_by_first_approval: int = 0
    chembl_rejected_withdrawn: int = 0
    pubchem_reranked: int = 0
    pubchem_resolved: int = 0
    pubchem_pugview_title_resolved: int = 0
    pubchem_parent_title_resolved: int = 0
    pubchem_pugview_calls: int = 0
    pubchem_pugview_success: int = 0
    pubchem_pugview_fail: int = 0
    pubchem_cas_title_resolved: int = 0
    pubchem_cas_candidates_considered: int = 0
    pubchem_inchi_cid_resolved: int = 0
    pubchem_smiles_400_inchi_fallback_attempted: int = 0
    gsrs_resolved_total: int = 0
    gsrs_resolved_by_unii: int = 0
    gsrs_resolved_by_inchikey_full: int = 0
    gsrs_resolved_by_inchikey_2d: int = 0
    gsrs_resolved_by_structure: int = 0
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
    drugcentral_resolved_exact: bool
    drugcentral_resolved_2d: bool
    chembl_resolved: bool
    chembl_approved_by_max_phase: bool
    chembl_approved_by_first_approval: bool
    chembl_rejected_withdrawn: bool
    pubchem_reranked: bool
    pubchem_resolved: bool
    pubchem_pugview_title_resolved: bool
    pubchem_parent_title_resolved: bool
    pubchem_pugview_calls: int
    pubchem_pugview_success: int
    pubchem_pugview_fail: int
    pubchem_cas_title_resolved: bool
    pubchem_cas_candidates_considered: int
    pubchem_inchi_cid_resolved: bool
    pubchem_smiles_400_inchi_fallback_attempted: bool
    gsrs_resolved_mode: str
    rxnorm_resolved: bool
    rxnorm_promoted: bool
    used_fallback: bool
    final_bad: bool
    final_quality_reasons: list[str]
    row_identifier: str
    unresolved_reason: str


def _process_row(
    row_index: int,
    row: dict[str, Any],
    only_fix_bad_display_names: bool,
    pubchem_client: CachedJsonClient | None,
    rxnorm_client: CachedJsonClient | None,
    unichem_client: CachedJsonClient | None,
    chembl_client: CachedJsonClient | None,
    gsrs_client: GSRSClient | None,
    drugcentral_index: DrugCentralIndex,
    enable_unichem: bool,
    enable_pubchem_rerank: bool,
    enable_gsrs: bool,
    openfda_client: CachedJsonClient | None = None,
    quality_fix_mode: str = "off",
) -> RowResult:
    current_display = _clean_text(row.get("display_name"))
    pubchem_iupac = _clean_text(row.get("pubchem_iupac_name"))
    strict_quality_mode = _clean_text(quality_fix_mode).lower() == "strict"
    display_bad = is_low_quality_display_name(
        current_display,
        pubchem_iupac,
        strict=strict_quality_mode,
    )
    evaluate = (not only_fix_bad_display_names) or display_bad
    row_updated = False
    used_fallback = False
    row_evaluated = False
    pubchem_resolved = False
    rxnorm_resolved = False
    drugcentral_resolved = False
    drugcentral_resolved_exact = False
    drugcentral_resolved_2d = False
    chembl_resolved = False
    pubchem_reranked = False
    pubchem_pugview_title_resolved = False
    pubchem_parent_title_resolved = False
    pubchem_pugview_calls = 0
    pubchem_pugview_success = 0
    pubchem_pugview_fail = 0
    pubchem_cas_title_resolved = False
    pubchem_cas_candidates_considered = 0
    pubchem_inchi_cid_resolved = False
    pubchem_smiles_400_inchi_fallback_attempted = False
    gsrs_resolved_mode = ""
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
            pubchem_inchi_cid_resolved = pubchem.inchi_cid_resolved
            pubchem_smiles_400_inchi_fallback_attempted = (
                pubchem.smiles_400_inchi_fallback_attempted
            )
            row_updated |= _maybe_set_field(row, "pubchem_cid_resolved", pubchem.cid)
            row_updated |= _maybe_set_field(
                row, "pubchem_iupac_name", pubchem.iupac_name
            )
            pubchem_iupac_for_fields = _clean_text(
                _pick_first(row.get("pubchem_iupac_name"), pubchem.iupac_name)
            )
            row_updated |= _maybe_set_field_overwrite_if_bad(
                row,
                "pubchem_record_title",
                pubchem.title,
                pubchem_iupac_name=pubchem_iupac_for_fields,
            )
            row_updated |= _maybe_set_field_overwrite_if_bad(
                row,
                "pubchem_name",
                pubchem.title,
                pubchem_iupac_name=pubchem_iupac_for_fields,
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
    display_bad = is_low_quality_display_name(
        current_display,
        pubchem_iupac,
        strict=strict_quality_mode,
    )
    if evaluate and display_bad:
        resolved_name = ""
        rxnorm_name = _clean_text(row.get("rxnorm_generic_name"))
        if rxnorm_name and is_good_display_name(rxnorm_name, pubchem_iupac):
            if rxnorm_name != current_display:
                row["display_name"] = rxnorm_name
                row_updated = True
            resolved_name = rxnorm_name
            rxnorm_promoted = True

        if not resolved_name and (drugcentral_index.full or drugcentral_index.two_d):
            drugcentral_name, match_mode, match_source, matched_key = _match_drugcentral_name(
                drugcentral_index,
                _drugcentral_lookup_candidates(row),
            )
            if drugcentral_name and is_good_display_name(drugcentral_name, pubchem_iupac):
                row_updated |= _maybe_set_field_overwrite_if_bad(
                    row,
                    "generic_name",
                    drugcentral_name,
                    pubchem_iupac_name=pubchem_iupac,
                )
                if drugcentral_name != current_display:
                    row["display_name"] = drugcentral_name
                    row_updated = True
                resolved_name = drugcentral_name
                drugcentral_resolved = True
                if match_mode == "2d":
                    drugcentral_resolved_2d = True
                else:
                    drugcentral_resolved_exact = True
                LOGGER.info(
                    "[drugcentral] row=%s mode=%s source=%s key=%s name=%s",
                    row_index,
                    match_mode,
                    match_source,
                    matched_key,
                    drugcentral_name,
                )

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

        if not resolved_name and enable_pubchem_rerank and pubchem_client:
            cid_value = _pick_first(pubchem.cid, _parse_existing_pubchem_cid(row))
            if cid_value:
                pubchem_iupac = _pick_first(
                    row.get("pubchem_iupac_name"),
                    pubchem.iupac_name,
                )
                last_mile = _resolve_pubchem_last_mile(
                    row=row,
                    pubchem_cid=cid_value,
                    pubchem_iupac_name=pubchem_iupac,
                    pubchem_client=pubchem_client,
                )
                pubchem_pugview_calls += last_mile.pugview_calls
                pubchem_pugview_success += last_mile.pugview_success
                pubchem_pugview_fail += last_mile.pugview_fail
                if last_mile.name and is_good_display_name(last_mile.name, pubchem_iupac):
                    if last_mile.name != current_display:
                        row["display_name"] = last_mile.name
                        row_updated = True
                    if last_mile.record_title:
                        row_updated |= _maybe_set_field_overwrite_if_bad(
                            row,
                            "pubchem_record_title",
                            last_mile.record_title,
                            pubchem_iupac_name=pubchem_iupac,
                        )
                    row_updated |= _maybe_set_field_overwrite_if_bad(
                        row,
                        "pubchem_name",
                        last_mile.name,
                        pubchem_iupac_name=pubchem_iupac,
                    )
                    if last_mile.parent_cid:
                        row_updated |= _maybe_set_field(
                            row,
                            "pubchem_parent_cid",
                            last_mile.parent_cid,
                        )
                    resolved_name = last_mile.name
                    pubchem_pugview_title_resolved = last_mile.pugview_title_resolved
                    pubchem_parent_title_resolved = last_mile.parent_title_resolved

        if not resolved_name and enable_pubchem_rerank and pubchem_client:
            pubchem_iupac = _pick_first(
                row.get("pubchem_iupac_name"),
                pubchem.iupac_name,
            )
            cas_last_mile = _resolve_pubchem_cas_last_mile(
                row=row,
                pubchem_iupac_name=pubchem_iupac,
                pubchem_client=pubchem_client,
            )
            pubchem_pugview_calls += cas_last_mile.pugview_calls
            pubchem_pugview_success += cas_last_mile.pugview_success
            pubchem_pugview_fail += cas_last_mile.pugview_fail
            pubchem_cas_candidates_considered += cas_last_mile.candidates_considered
            if cas_last_mile.title_resolved:
                pubchem_cas_title_resolved = True
            if cas_last_mile.name and is_good_display_name(
                cas_last_mile.name, pubchem_iupac
            ):
                if cas_last_mile.name != current_display:
                    row["display_name"] = cas_last_mile.name
                    row_updated = True
                if cas_last_mile.record_title:
                    row_updated |= _maybe_set_field_overwrite_if_bad(
                        row,
                        "pubchem_record_title",
                        cas_last_mile.record_title,
                        pubchem_iupac_name=pubchem_iupac,
                    )
                row_updated |= _maybe_set_field_overwrite_if_bad(
                    row,
                    "pubchem_name",
                    cas_last_mile.name,
                    pubchem_iupac_name=pubchem_iupac,
                )
                if cas_last_mile.source_cid:
                    row_updated |= _maybe_set_field(
                        row,
                        "pubchem_cid_resolved",
                        cas_last_mile.source_cid,
                    )
                resolved_name = cas_last_mile.name

        if not resolved_name and enable_gsrs and gsrs_client:
            gsrs_resolution = resolve_gsrs_last_mile(row, pubchem_iupac, gsrs_client)
            if gsrs_resolution.name and is_good_display_name(
                gsrs_resolution.name, pubchem_iupac
            ):
                if gsrs_resolution.name != current_display:
                    row["display_name"] = gsrs_resolution.name
                    row_updated = True
                resolved_name = gsrs_resolution.name
                gsrs_resolved_mode = gsrs_resolution.mode
                LOGGER.info(
                    "[gsrs] row=%s mode=%s name=%s",
                    row_index,
                    gsrs_resolution.mode,
                    gsrs_resolution.name,
                )

        if not resolved_name and openfda_client:
            for unii in _collect_uniis_from_row(row):
                openfda_name = resolve_openfda_substance_name(unii, openfda_client)
                if not openfda_name:
                    continue
                openfda_name = _collapse_spaces(openfda_name)
                if not is_good_display_name(openfda_name, pubchem_iupac):
                    continue
                if openfda_name != current_display:
                    row["display_name"] = openfda_name
                    row_updated = True
                resolved_name = openfda_name
                LOGGER.info(
                    "[openfda] row=%s unii=%s name=%s",
                    row_index,
                    unii,
                    openfda_name,
                )
                break

        if not resolved_name and rxnorm_client:
            name_hint = _extract_salvage_name_hint(row)
            if name_hint:
                rxnorm_hint = resolve_rxnorm(name_hint, rxnorm_client)
                if rxnorm_hint.rxcui or rxnorm_hint.name:
                    rxnorm_resolved = True
                row_updated |= _maybe_set_field(row, "rxnorm_generic_name", rxnorm_hint.name)
                row_updated |= _maybe_set_field(row, "rxnorm_rxcui", rxnorm_hint.rxcui)
                hint_name = _clean_text(rxnorm_hint.name)
                if hint_name and is_good_display_name(hint_name, pubchem_iupac):
                    if hint_name != current_display:
                        row["display_name"] = hint_name
                        row_updated = True
                    resolved_name = hint_name
                    rxnorm_promoted = True
                    LOGGER.info(
                        "[rxnorm-salvage] row=%s hint=%s name=%s",
                        row_index,
                        name_hint,
                        hint_name,
                    )

        if not resolved_name:
            fallback = _fallback_display_name(row, row_index)
            if fallback and fallback != current_display:
                row["display_name"] = fallback
                row_updated = True
            used_fallback = True

    final_display = _clean_text(row.get("display_name"))
    final_quality_reasons = display_name_quality_reasons(
        final_display,
        _clean_text(row.get("pubchem_iupac_name")),
        strict=strict_quality_mode,
    )
    final_bad = bool(final_quality_reasons)
    row_identifier = (
        _row_identifier(row, row_index) if (used_fallback or final_bad) else ""
    )
    unresolved_reason = ""
    if final_bad:
        unresolved_reason = f"low_quality_after_resolution:{','.join(final_quality_reasons)}"
    elif used_fallback:
        unresolved_reason = "used_fallback"
    return RowResult(
        row_index=row_index,
        row=row,
        row_evaluated=row_evaluated,
        row_updated=row_updated,
        drugcentral_resolved=drugcentral_resolved,
        drugcentral_resolved_exact=drugcentral_resolved_exact,
        drugcentral_resolved_2d=drugcentral_resolved_2d,
        chembl_resolved=chembl_resolved,
        chembl_approved_by_max_phase=chembl_approved_by_max_phase,
        chembl_approved_by_first_approval=chembl_approved_by_first_approval,
        chembl_rejected_withdrawn=chembl_rejected_withdrawn,
        pubchem_reranked=pubchem_reranked,
        pubchem_resolved=pubchem_resolved,
        pubchem_pugview_title_resolved=pubchem_pugview_title_resolved,
        pubchem_parent_title_resolved=pubchem_parent_title_resolved,
        pubchem_pugview_calls=pubchem_pugview_calls,
        pubchem_pugview_success=pubchem_pugview_success,
        pubchem_pugview_fail=pubchem_pugview_fail,
        pubchem_cas_title_resolved=pubchem_cas_title_resolved,
        pubchem_cas_candidates_considered=pubchem_cas_candidates_considered,
        pubchem_inchi_cid_resolved=pubchem_inchi_cid_resolved,
        pubchem_smiles_400_inchi_fallback_attempted=pubchem_smiles_400_inchi_fallback_attempted,
        gsrs_resolved_mode=gsrs_resolved_mode,
        rxnorm_resolved=rxnorm_resolved,
        rxnorm_promoted=rxnorm_promoted,
        used_fallback=used_fallback,
        final_bad=final_bad,
        final_quality_reasons=final_quality_reasons,
        row_identifier=row_identifier,
        unresolved_reason=unresolved_reason,
    )


def _process_row_item(
    item: tuple[int, dict[str, Any]],
    only_fix_bad_display_names: bool,
    pubchem_client: CachedJsonClient | None,
    rxnorm_client: CachedJsonClient | None,
    unichem_client: CachedJsonClient | None,
    chembl_client: CachedJsonClient | None,
    gsrs_client: GSRSClient | None,
    drugcentral_index: DrugCentralIndex,
    enable_unichem: bool,
    enable_pubchem_rerank: bool,
    enable_gsrs: bool,
    openfda_client: CachedJsonClient | None = None,
    quality_fix_mode: str = "off",
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
        gsrs_client=gsrs_client,
        openfda_client=openfda_client,
        drugcentral_index=drugcentral_index,
        enable_unichem=enable_unichem,
        enable_pubchem_rerank=enable_pubchem_rerank,
        enable_gsrs=enable_gsrs,
        quality_fix_mode=quality_fix_mode,
    )


def fill_mapping_names(
    in_csv: Path,
    out_csv: Path,
    *,
    only_fix_bad_display_names: bool,
    quality_fix_mode: str,
    require_all_named: bool,
    report_path: Path,
    flagged_report_path: Path,
    drugcentral_structures_tsv: Path | None,
    enable_unichem: bool,
    enable_pubchem_rerank: bool,
    enable_gsrs: bool,
    gsrs_base_url: str,
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
    quality_mode = _clean_text(quality_fix_mode).lower() or "off"
    if quality_mode not in {"off", "strict"}:
        raise ValueError("--quality_fix_mode must be one of: off, strict.")
    strict_quality_mode = quality_mode == "strict"

    if cache_mode != "off":
        cache_dir.mkdir(parents=True, exist_ok=True)
        if clear_cache:
            _clear_cache_dir(cache_dir)
    elif clear_cache:
        LOGGER.warning("[cache] clear_cache ignored because cache_mode=off")

    drugcentral_index = DrugCentralIndex()
    if drugcentral_structures_tsv:
        drugcentral_index = load_drugcentral_index(drugcentral_structures_tsv)
    else:
        LOGGER.warning(
            "[drugcentral] --drugcentral_structures_tsv not provided; skipping DrugCentral stage."
        )

    pubchem_client: CachedJsonClient | None = None
    rxnorm_client: CachedJsonClient | None = None
    if enable_pubchem_rerank or strict_quality_mode:
        rxnorm_client = CachedJsonClient(
            cache_dir / "rxnorm",
            sleep,
            cache_mode=cache_mode,
            rate_limiter=rate_limiter,
            max_retries=max_retries,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
        )
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

    gsrs_client: GSRSClient | None = None
    if enable_gsrs:
        gsrs_http_client = CachedJsonClient(
            cache_dir / "gsrs",
            sleep,
            cache_mode=cache_mode,
            rate_limiter=rate_limiter,
            max_retries=max_retries,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
        )
        gsrs_client = GSRSClient(gsrs_base_url, gsrs_http_client)

    openfda_client: CachedJsonClient | None = None
    if strict_quality_mode:
        openfda_client = CachedJsonClient(
            cache_dir / "openfda",
            sleep,
            cache_mode=cache_mode,
            rate_limiter=rate_limiter,
            max_retries=max_retries,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
        )

    unresolved: list[str] = []
    flagged_quality_rows: list[dict[str, str]] = []
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
        fieldnames = list(reader.fieldnames or [])
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
                gsrs_client=gsrs_client,
                openfda_client=openfda_client,
                drugcentral_index=drugcentral_index,
                enable_unichem=enable_unichem,
                enable_pubchem_rerank=enable_pubchem_rerank,
                enable_gsrs=enable_gsrs,
                quality_fix_mode=quality_mode,
            )

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                for result in executor.map(process_row, enumerate(reader, start=1)):
                    summary.total_rows += 1
                    if result.row_evaluated:
                        summary.rows_evaluated += 1
                    if result.drugcentral_resolved:
                        summary.drugcentral_resolved += 1
                    if result.drugcentral_resolved_exact:
                        summary.drugcentral_resolved_exact += 1
                    if result.drugcentral_resolved_2d:
                        summary.drugcentral_resolved_2d += 1
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
                    if result.pubchem_pugview_title_resolved:
                        summary.pubchem_pugview_title_resolved += 1
                    if result.pubchem_parent_title_resolved:
                        summary.pubchem_parent_title_resolved += 1
                    summary.pubchem_pugview_calls += result.pubchem_pugview_calls
                    summary.pubchem_pugview_success += result.pubchem_pugview_success
                    summary.pubchem_pugview_fail += result.pubchem_pugview_fail
                    if result.pubchem_cas_title_resolved:
                        summary.pubchem_cas_title_resolved += 1
                    summary.pubchem_cas_candidates_considered += (
                        result.pubchem_cas_candidates_considered
                    )
                    if result.pubchem_inchi_cid_resolved:
                        summary.pubchem_inchi_cid_resolved += 1
                    if result.pubchem_smiles_400_inchi_fallback_attempted:
                        summary.pubchem_smiles_400_inchi_fallback_attempted += 1
                    if result.gsrs_resolved_mode:
                        summary.gsrs_resolved_total += 1
                        if result.gsrs_resolved_mode == "unii":
                            summary.gsrs_resolved_by_unii += 1
                        elif result.gsrs_resolved_mode == "inchikey_full":
                            summary.gsrs_resolved_by_inchikey_full += 1
                        elif result.gsrs_resolved_mode == "inchikey_2d":
                            summary.gsrs_resolved_by_inchikey_2d += 1
                        elif result.gsrs_resolved_mode == "structure":
                            summary.gsrs_resolved_by_structure += 1
                    if result.rxnorm_resolved:
                        summary.rxnorm_resolved += 1
                    if result.rxnorm_promoted:
                        summary.rxnorm_promoted += 1
                    if result.used_fallback:
                        summary.fallbacks_used += 1
                    if result.row_updated:
                        summary.rows_updated += 1
                    if result.row_identifier:
                        reason = _clean_text(result.unresolved_reason)
                        entry = (
                            f"{result.row_identifier}\t{reason}"
                            if reason
                            else result.row_identifier
                        )
                        unresolved.append(entry)
                    if result.final_bad:
                        has_bad = True
                    if result.final_quality_reasons:
                        flagged_quality_rows.append(
                            {
                                "sdf_title": _clean_text(result.row.get("sdf_title")),
                                "inchikey": _row_inchikey(result.row),
                                "current_display_name": _clean_text(
                                    result.row.get("display_name")
                                ),
                                "reason_codes": ";".join(result.final_quality_reasons),
                                "pubchem_cid_resolved": _clean_text(
                                    result.row.get("pubchem_cid_resolved")
                                ),
                                "pubchem_record_title": _clean_text(
                                    result.row.get("pubchem_record_title")
                                ),
                                "pubchem_unii_list": _clean_text(
                                    result.row.get("pubchem_unii_list")
                                ),
                            }
                        )

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

    flagged_report_path.parent.mkdir(parents=True, exist_ok=True)
    flagged_headers = [
        "sdf_title",
        "inchikey",
        "current_display_name",
        "reason_codes",
        "pubchem_cid_resolved",
        "pubchem_record_title",
        "pubchem_unii_list",
    ]
    with flagged_report_path.open("w", newline="", encoding="utf-8") as flagged_handle:
        flagged_writer = csv.DictWriter(flagged_handle, fieldnames=flagged_headers)
        flagged_writer.writeheader()
        for row in flagged_quality_rows:
            flagged_writer.writerow({name: row.get(name, "") for name in flagged_headers})

    print(f"total rows: {summary.total_rows}")
    print(f"rows evaluated: {summary.rows_evaluated}")
    print(f"rows updated: {summary.rows_updated}")
    print(f"drugcentral resolved: {summary.drugcentral_resolved}")
    print(f"drugcentral resolved exact: {summary.drugcentral_resolved_exact}")
    print(f"drugcentral resolved 2d: {summary.drugcentral_resolved_2d}")
    print(f"chembl resolved: {summary.chembl_resolved}")
    print(f"chembl approved by max_phase: {summary.chembl_approved_by_max_phase}")
    print(
        f"chembl approved by first_approval: {summary.chembl_approved_by_first_approval}"
    )
    if summary.chembl_rejected_withdrawn:
        print(f"chembl rejected withdrawn: {summary.chembl_rejected_withdrawn}")
    print(f"pubchem reranked: {summary.pubchem_reranked}")
    print(f"pubchem resolved: {summary.pubchem_resolved}")
    print(f"pubchem pugview title resolved: {summary.pubchem_pugview_title_resolved}")
    print(f"pubchem parent title resolved: {summary.pubchem_parent_title_resolved}")
    print(f"pubchem cas title resolved: {summary.pubchem_cas_title_resolved}")
    print(
        "pubchem cas candidates considered: "
        f"{summary.pubchem_cas_candidates_considered}"
    )
    print(f"pubchem pugview calls: {summary.pubchem_pugview_calls}")
    print(f"pubchem pugview success: {summary.pubchem_pugview_success}")
    print(f"pubchem pugview fail: {summary.pubchem_pugview_fail}")
    print(
        f"pubchem smiles 400 inchi fallback attempted: "
        f"{summary.pubchem_smiles_400_inchi_fallback_attempted}"
    )
    print(f"pubchem inchi cid resolved: {summary.pubchem_inchi_cid_resolved}")
    print(f"gsrs resolved total: {summary.gsrs_resolved_total}")
    print(f"gsrs resolved by unii: {summary.gsrs_resolved_by_unii}")
    print(f"gsrs resolved by inchikey full: {summary.gsrs_resolved_by_inchikey_full}")
    print(f"gsrs resolved by inchikey 2d: {summary.gsrs_resolved_by_inchikey_2d}")
    print(f"gsrs resolved by structure: {summary.gsrs_resolved_by_structure}")
    print(f"rxnorm resolved: {summary.rxnorm_resolved}")
    print(f"rxnorm promoted: {summary.rxnorm_promoted}")
    print(f"fallbacks used: {summary.fallbacks_used}")
    print(f"flagged name-quality rows: {len(flagged_quality_rows)}")

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
        "--quality_fix_mode",
        choices=("off", "strict"),
        default="off",
        help=(
            "Additional low-quality name cleanup mode. "
            "'strict' re-processes CAS/formula/UNII/spaced-IUPAC-like display names."
        ),
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
            "from DrugCentral). If omitted, the tool auto-enables "
            f"{DEFAULT_DRUGCENTRAL_STRUCTURES_TSV} when present."
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
        "--enable_gsrs",
        type=_parse_bool,
        default=True,
        help="Enable GSRS last-mile resolution for still-unresolved rows (default: true).",
    )
    parser.add_argument(
        "--gsrs_base_url",
        default=DEFAULT_GSRS_BASE_URL,
        help="GSRS API base URL.",
    )
    parser.add_argument(
        "--report_path",
        help="Path for unresolved report (default: <out_csv>.unresolved.txt).",
    )
    parser.add_argument(
        "--flagged_report_path",
        help=(
            "Path for flagged low-quality name report "
            "(default: <out_csv with .flagged_name_quality.csv suffix>)."
        ),
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
    if args.flagged_report_path:
        flagged_report_path = Path(args.flagged_report_path)
    else:
        if out_path.suffix:
            flagged_report_path = out_path.with_suffix(".flagged_name_quality.csv")
        else:
            flagged_report_path = out_path.with_name(
                f"{out_path.name}.flagged_name_quality.csv"
            )
    if args.drugcentral_structures_tsv:
        drugcentral_path = Path(args.drugcentral_structures_tsv)
    elif DEFAULT_DRUGCENTRAL_STRUCTURES_TSV.exists():
        drugcentral_path = DEFAULT_DRUGCENTRAL_STRUCTURES_TSV
        print(f"[drugcentral] auto-using local TSV: {drugcentral_path}")
    else:
        drugcentral_path = None
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
            quality_fix_mode=_clean_text(args.quality_fix_mode).lower() or "off",
            require_all_named=require_all_named,
            report_path=report_path,
            flagged_report_path=flagged_report_path,
            drugcentral_structures_tsv=drugcentral_path,
            enable_unichem=bool(args.enable_unichem),
            enable_pubchem_rerank=bool(args.enable_pubchem_rerank),
            enable_gsrs=bool(args.enable_gsrs),
            gsrs_base_url=str(args.gsrs_base_url).strip() or DEFAULT_GSRS_BASE_URL,
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
