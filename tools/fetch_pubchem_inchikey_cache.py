#!/usr/bin/env python3
"""Cache exact PubChem InChIKey property lookups with offline provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from rdkit import Chem


SCHEMA = "atlas.pubchem-inchikey-cache.v1"
SOURCE = "PubChem PUG REST exact InChIKey property lookup"
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPERTIES = (
    "CanonicalSMILES",
    "ConnectivitySMILES",
    "IsomericSMILES",
    "InChI",
    "InChIKey",
    "IUPACName",
    "MolecularFormula",
    "Title",
)
INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
RETRYABLE = {408, 425, 429, 500, 502, 503, 504}


class CacheError(RuntimeError):
    """Raised when cache or remote content violates the lookup contract."""


class RateLimiter:
    def __init__(self, requests_per_second: float) -> None:
        if requests_per_second <= 0:
            raise ValueError("requests_per_second must be positive")
        self.interval = 1.0 / requests_per_second
        self.last: float | None = None

    def wait(self) -> None:
        now = time.monotonic()
        if self.last is not None:
            time.sleep(max(0.0, self.interval - (now - self.last)))
        self.last = time.monotonic()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _clean(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def _truth(value: Any) -> bool:
    return _clean(value).casefold() in {"1", "true", "yes", "y"}


def _validate_key(value: str, context: str) -> str:
    key = value.strip().upper()
    if not INCHIKEY_RE.fullmatch(key):
        raise ValueError(f"invalid InChIKey {value!r} at {context}")
    return key


def read_queries(args: argparse.Namespace) -> tuple[list[str], dict[str, list[dict[str, str]]], int]:
    provenance: dict[str, list[dict[str, str]]] = {}
    selected_rows = 0
    for csv_path in args.input_csv:
        selected_rows += _read_query_csv(csv_path, args, provenance)
    for index, raw_key in enumerate(args.inchikey, start=1):
        key = _validate_key(raw_key, f"--inchikey item {index}")
        provenance.setdefault(key, []).append({"input_csv": "", "row_number": "", "key_column": "--inchikey"})
    if not provenance:
        raise ValueError("no InChIKey queries selected")
    return sorted(provenance), provenance, selected_rows


def _read_query_csv(
    csv_path: Path,
    args: argparse.Namespace,
    provenance: dict[str, list[dict[str, str]]],
) -> int:
    selected = 0
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {args.key_column}
        if not args.all_rows:
            required.update({args.validated_column, args.usable_column})
        if args.fallback_key_column:
            required.add(args.fallback_key_column)
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{csv_path} lacks required columns: {', '.join(missing)}")
        for row_number, row in enumerate(reader, start=2):
            if _skip_query_row(row, args):
                continue
            selected += 1
            key, source_column = _row_query_key(row, args, csv_path, row_number)
            provenance.setdefault(key, []).append(
                {
                    "input_csv": str(csv_path.resolve()),
                    "row_number": str(row_number),
                    "mapping_row_number": _clean(row.get("mapping_row_number")),
                    "rdk_id": _clean(row.get("rdk_id")),
                    "key_column": source_column,
                }
            )
    return selected


def _skip_query_row(row: Mapping[str, Any], args: argparse.Namespace) -> bool:
    return not args.all_rows and (
        _truth(row.get(args.validated_column)) or not _truth(row.get(args.usable_column))
    )


def _row_query_key(
    row: Mapping[str, Any], args: argparse.Namespace, csv_path: Path, row_number: int
) -> tuple[str, str]:
    source_column = args.key_column
    raw_key = _clean(row.get(source_column))
    if not raw_key and args.fallback_key_column:
        source_column = args.fallback_key_column
        raw_key = _clean(row.get(source_column))
    return _validate_key(raw_key, f"{csv_path}:{row_number}:{source_column}"), source_column


def _url(key: str) -> str:
    properties = ",".join(PROPERTIES)
    encoded = urllib.parse.quote(key, safe="")
    return f"{BASE_URL}/compound/inchikey/{encoded}/property/{properties}/JSON"


def _query_id(key: str) -> str:
    return _sha(_json_bytes({"inchikey": key, "properties": list(PROPERTIES), "url": _url(key)}))


def _smiles(record: Mapping[str, Any]) -> str:
    for name in ("IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES"):
        value = record.get(name)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def parse_hit(payload: bytes, query_key: str) -> tuple[list[dict[str, Any]], dict[str, list[str]]]:
    try:
        rows = json.loads(payload)["PropertyTable"]["Properties"]
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        raise CacheError(f"invalid PubChem response for {query_key}: {exc}") from exc
    if not isinstance(rows, list):
        raise CacheError(f"PubChem Properties is not a list for {query_key}")
    records: list[dict[str, Any]] = []
    invalid: dict[str, list[str]] = {}
    seen: set[int] = set()
    for row in rows:
        cid, errors = _validate_hit_record(row, query_key, seen)
        if errors:
            invalid[str(cid)] = errors
        records.append(row)
    if not records:
        raise CacheError(f"successful PubChem response contained no records for {query_key}")
    return records, invalid


def _validate_hit_record(row: Any, query_key: str, seen: set[int]) -> tuple[int, list[str]]:
    if not isinstance(row, dict) or not isinstance(row.get("CID"), int):
        raise CacheError(f"PubChem returned a malformed record for {query_key}")
    cid = int(row["CID"])
    if cid in seen:
        raise CacheError(f"PubChem returned duplicate CID {cid} for {query_key}")
    seen.add(cid)
    errors: list[str] = []
    if _clean(row.get("InChIKey")).upper() != query_key:
        errors.append(f"returned InChIKey {_clean(row.get('InChIKey'))!r} is not the exact query key")
    smiles = _smiles(row)
    if not smiles:
        errors.append("no SMILES property")
    elif Chem.MolFromSmiles(smiles) is None:
        errors.append("SMILES is not RDKit-parseable")
    return cid, errors


def fetch(
    key: str,
    *,
    limiter: RateLimiter,
    timeout: float,
    retries: int,
    user_agent: str,
) -> tuple[int, bytes, dict[str, str]]:
    request = urllib.request.Request(_url(key), headers={"Accept": "application/json", "User-Agent": user_agent})
    for attempt in range(retries + 1):
        limiter.wait()
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                headers = {
                    name.lower(): value
                    for name, value in response.headers.items()
                    if name.lower() in {"content-type", "date", "server", "x-throttling-control"}
                }
                return int(response.status), response.read(), headers
        except urllib.error.HTTPError as exc:
            handled = _inchikey_http_error(exc, key, attempt, retries)
            if isinstance(handled, tuple):
                return handled
            delay = handled
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt >= retries:
                raise CacheError(f"PubChem request failed for {key}: {exc}") from exc
            delay = 0.75 * (2**attempt)
        time.sleep(min(30.0, delay) + random.uniform(0.0, 0.25))
    raise AssertionError("retry loop exhausted")


def _inchikey_http_error(
    exc: urllib.error.HTTPError, key: str, attempt: int, retries: int
) -> float | tuple[int, bytes, dict[str, str]]:
    payload = exc.read()
    if exc.code == 404:
        return 404, payload, {"content-type": exc.headers.get("Content-Type", "")}
    if exc.code not in RETRYABLE or attempt >= retries:
        raise CacheError(f"PubChem HTTP {exc.code} for {key}: {payload[:500]!r}") from exc
    retry_after = exc.headers.get("Retry-After")
    try:
        return float(retry_after) if retry_after else 0.75 * (2**attempt)
    except ValueError:
        return 0.75 * (2**attempt)


def _paths(cache_dir: Path, key: str) -> tuple[Path, Path]:
    query_id = _query_id(key)
    query_dir = cache_dir / "queries"
    return query_dir / f"{query_id}.json", query_dir / f"{query_id}.metadata.json"


def write_query(
    cache_dir: Path,
    key: str,
    status: int,
    payload: bytes,
    headers: Mapping[str, str],
) -> dict[str, Any]:
    raw_path, metadata_path = _paths(cache_dir, key)
    if status == 404:
        records: list[dict[str, Any]] = []
        invalid: dict[str, list[str]] = {}
        outcome = "no_hit"
    else:
        records, invalid = parse_hit(payload, key)
        outcome = "one_hit" if len(records) == 1 else "one_to_many"
    metadata = {
        "schema": SCHEMA,
        "query_id": _query_id(key),
        "query_inchikey": key,
        "source": SOURCE,
        "source_url": _url(key),
        "http_status": status,
        "outcome": outcome,
        "returned_cids": [record["CID"] for record in records],
        "invalid_structures": invalid,
        "retrieved_at_utc": _now(),
        "response_headers": dict(headers),
        "content_file": raw_path.name,
        "content_bytes": len(payload),
        "content_sha256": _sha(payload),
    }
    _atomic_write(raw_path, payload)
    _atomic_write(metadata_path, _json_bytes(metadata))
    return {"metadata": metadata, "records": records}


def load_cache(cache_dir: Path) -> dict[str, dict[str, Any]]:
    cached: dict[str, dict[str, Any]] = {}
    query_dir = cache_dir / "queries"
    if not query_dir.exists():
        return cached
    for metadata_path in sorted(query_dir.glob("*.metadata.json")):
        key, metadata, records = _load_cached_query(metadata_path, query_dir)
        cached[key] = {"metadata": metadata, "records": records}
    return cached


def _load_cached_query(metadata_path: Path, query_dir: Path) -> tuple[str, dict[str, Any], list[dict[str, Any]]]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    key = _clean(metadata.get("query_inchikey"))
    _validate_cached_query_identity(metadata, metadata_path, key)
    raw_path = query_dir / _clean(metadata.get("content_file"))
    payload = raw_path.read_bytes()
    if _sha(payload) != metadata.get("content_sha256") or len(payload) != metadata.get("content_bytes"):
        raise CacheError(f"raw response hash/size mismatch for {key}")
    records, invalid, outcome = _cached_query_content(metadata.get("http_status"), payload, key)
    if metadata.get("outcome") != outcome:
        raise CacheError(f"outcome mismatch for cached query {key}")
    if metadata.get("returned_cids") != [record["CID"] for record in records]:
        raise CacheError(f"CID coverage mismatch for cached query {key}")
    if metadata.get("invalid_structures") != invalid:
        raise CacheError(f"structure validation mismatch for cached query {key}")
    return key, metadata, records


def _validate_cached_query_identity(metadata: Mapping[str, Any], metadata_path: Path, key: str) -> None:
    if metadata.get("schema") != SCHEMA or not INCHIKEY_RE.fullmatch(key):
        raise CacheError(f"invalid cache metadata {metadata_path}")
    query_id = _query_id(key)
    if metadata.get("query_id") != query_id or metadata_path.stem.removesuffix(".metadata") != query_id:
        raise CacheError(f"query identity mismatch in {metadata_path}")


def _cached_query_content(
    status: Any, payload: bytes, key: str
) -> tuple[list[dict[str, Any]], dict[str, list[str]], str]:
    if status == 404:
        return [], {}, "no_hit"
    if status != 200:
        raise CacheError(f"unsupported cached HTTP status {status!r} for {key}")
    records, invalid = parse_hit(payload, key)
    return records, invalid, "one_hit" if len(records) == 1 else "one_to_many"


def write_outputs(
    cache_dir: Path,
    keys: Sequence[str],
    provenance: Mapping[str, list[dict[str, str]]],
    results: Mapping[str, dict[str, Any]],
    *,
    selected_rows: int,
    offline: bool,
    fetched: int,
) -> dict[str, Any]:
    lines: list[bytes] = []
    outcomes: dict[str, int] = {"one_hit": 0, "one_to_many": 0, "no_hit": 0}
    invalid_queries = 0
    for key in keys:
        result = results[key]
        metadata = result["metadata"]
        outcome = metadata["outcome"]
        outcomes[outcome] += 1
        if metadata["invalid_structures"]:
            invalid_queries += 1
        normalized = {
            "query_inchikey": key,
            "outcome": outcome,
            "http_status": metadata["http_status"],
            "returned_cids": metadata["returned_cids"],
            "invalid_structures": metadata["invalid_structures"],
            "records": result["records"],
            "input_rows": provenance[key],
            "source_content_sha256": metadata["content_sha256"],
        }
        lines.append(json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode() + b"\n")
    payload = b"".join(lines)
    records_path = cache_dir / "records.jsonl"
    _atomic_write(records_path, payload)
    manifest = {
        "schema": SCHEMA,
        "created_at_utc": _now(),
        "source": SOURCE,
        "offline": offline,
        "selected_input_rows": selected_rows,
        "requested_unique_inchikeys": len(keys),
        "covered_unique_inchikeys": len(results),
        "coverage_percent": round(100.0 * len(results) / len(keys), 6),
        "outcome_counts": outcomes,
        "invalid_structure_query_count": invalid_queries,
        "fetched_query_count": fetched,
        "reused_query_count": len(keys) - fetched,
        "requested_inchikeys": list(keys),
        "records_file": records_path.name,
        "records_bytes": len(payload),
        "records_sha256": _sha(payload),
    }
    _atomic_write(cache_dir / "manifest.json", _json_bytes(manifest))
    return manifest


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input-csv", action="append", type=Path, default=[])
    result.add_argument("--key-column", default="mapping_exact_inchikey")
    result.add_argument("--fallback-key-column", default="identity_exact_inchikey")
    result.add_argument("--validated-column", default="identity_structure_validated")
    result.add_argument("--usable-column", default="legacy_file_usable")
    result.add_argument("--all-rows", action="store_true", help="Do not apply validated=false and usable=true filters")
    result.add_argument("--inchikey", action="append", default=[])
    result.add_argument("--cache-dir", type=Path, required=True)
    result.add_argument("--requests-per-second", type=float, default=3.0)
    result.add_argument("--timeout", type=float, default=60.0)
    result.add_argument("--retries", type=int, default=5)
    result.add_argument("--contact-email")
    result.add_argument("--offline", action="store_true")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        manifest, cache_dir = _run_inchikey_cache(args)
    except (CacheError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        f"PubChem InChIKey cache: rows={manifest['selected_input_rows']} "
        f"keys={manifest['requested_unique_inchikeys']} outcomes={manifest['outcome_counts']} "
        f"invalid={manifest['invalid_structure_query_count']} cache_dir={cache_dir}"
    )
    return 0


def _run_inchikey_cache(args: argparse.Namespace) -> tuple[dict[str, Any], Path]:
    keys, provenance, selected_rows = read_queries(args)
    cache_dir = args.cache_dir.expanduser().resolve()
    results = load_cache(cache_dir)
    missing = [key for key in keys if key not in results]
    if args.offline and missing:
        raise CacheError(f"offline cache lacks {len(missing)} of {len(keys)} exact InChIKeys")
    fetched = _fetch_inchikey_queries(args, cache_dir, missing, results)
    manifest = write_outputs(
        cache_dir,
        keys,
        provenance,
        {key: results[key] for key in keys},
        selected_rows=selected_rows,
        offline=args.offline,
        fetched=fetched,
    )
    return manifest, cache_dir


def _fetch_inchikey_queries(
    args: argparse.Namespace,
    cache_dir: Path,
    missing: Sequence[str],
    results: dict[str, dict[str, Any]],
) -> int:
    if args.offline:
        return 0
    limiter = RateLimiter(args.requests_per_second)
    contact = f"; contact={args.contact_email}" if args.contact_email else ""
    for index, key in enumerate(missing, start=1):
        status, payload, headers = fetch(
            key,
            limiter=limiter,
            timeout=args.timeout,
            retries=args.retries,
            user_agent=f"atlas-fda-pubchem-inchikey-cache/1.0{contact}",
        )
        results[key] = write_query(cache_dir, key, status, payload, headers)
        if index % 25 == 0 or index == len(missing):
            print(f"PubChem exact InChIKey progress: {index}/{len(missing)}", file=sys.stderr)
    return len(missing)


if __name__ == "__main__":
    raise SystemExit(main())
