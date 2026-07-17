#!/usr/bin/env python3
"""Build a validated, provenance-preserving PubChem CID property cache.

This utility is intentionally independent of FDA mapping consumers.  It reads
CIDs, reuses any valid cached batch containing them, fetches only missing CIDs,
and emits a normalized JSONL index plus a run manifest.  Raw PubChem responses
are retained byte-for-byte with SHA-256 sidecars so later runs can be offline.
"""

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
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from rdkit import Chem


CACHE_SCHEMA = "atlas.pubchem-cid-cache.v1"
SOURCE_NAME = "PubChem PUG REST"
SOURCE_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
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
CID_SPLIT_RE = re.compile(r"[\s,;|]+")
RETRYABLE_HTTP_CODES = {408, 425, 429, 500, 502, 503, 504}


class CacheError(RuntimeError):
    """Raised when fetched or cached PubChem content fails validation."""


@dataclass(frozen=True)
class CachedRecord:
    record: dict[str, Any]
    batch_id: str
    content_sha256: str


class RateLimiter:
    def __init__(self, requests_per_second: float) -> None:
        if requests_per_second <= 0:
            raise ValueError("requests_per_second must be positive")
        self._interval = 1.0 / requests_per_second
        self._last_request: float | None = None

    def wait(self) -> None:
        now = time.monotonic()
        if self._last_request is not None:
            delay = self._interval - (now - self._last_request)
            if delay > 0:
                time.sleep(delay)
        self._last_request = time.monotonic()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_json_bytes(value: Any) -> bytes:
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


def _parse_cid_token(token: str, context: str) -> int:
    value = token.strip()
    if not value:
        raise ValueError(f"empty CID token at {context}")
    if re.fullmatch(r"[0-9]+\.0+", value):
        value = value.split(".", 1)[0]
    if not value.isdigit() or int(value) <= 0:
        raise ValueError(f"invalid PubChem CID {token!r} at {context}")
    return int(value)


def _parse_cid_cell(value: str, context: str) -> list[int]:
    stripped = value.strip()
    if not stripped:
        return []
    return [_parse_cid_token(token, context) for token in CID_SPLIT_RE.split(stripped) if token]


def read_cids(csv_paths: Sequence[Path], column: str, direct: Sequence[str]) -> list[int]:
    cids: set[int] = set()
    for path in csv_paths:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or column not in reader.fieldnames:
                raise ValueError(f"{path} has no {column!r} column")
            for row_number, row in enumerate(reader, start=2):
                cids.update(_parse_cid_cell(row.get(column, ""), f"{path}:{row_number}:{column}"))
    for item_number, value in enumerate(direct, start=1):
        cids.update(_parse_cid_cell(value, f"--cid item {item_number}"))
    if not cids:
        raise ValueError("no PubChem CIDs were supplied")
    return sorted(cids)


def _batch_id(cids: Sequence[int]) -> str:
    request_identity = {
        "cids": list(cids),
        "properties": list(PROPERTIES),
        "source_url": SOURCE_BASE_URL,
    }
    return _sha256(_canonical_json_bytes(request_identity))


def _property_url() -> str:
    properties = ",".join(PROPERTIES)
    return f"{SOURCE_BASE_URL}/compound/cid/property/{properties}/JSON"


def _structure_smiles(record: dict[str, Any]) -> str:
    for key in ("IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES"):
        value = record.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def validate_record(record: dict[str, Any], requested: set[int]) -> list[str]:
    errors: list[str] = []
    cid_value = record.get("CID")
    if not isinstance(cid_value, int) or cid_value <= 0:
        return ["CID is absent or is not a positive integer"]
    if cid_value not in requested:
        errors.append("CID was not requested")
    inchikey = record.get("InChIKey")
    if not isinstance(inchikey, str) or not INCHIKEY_RE.fullmatch(inchikey):
        errors.append("InChIKey is absent or malformed")
    smiles = _structure_smiles(record)
    if not smiles:
        errors.append("no SMILES structure property was returned")
    elif Chem.MolFromSmiles(smiles) is None:
        errors.append("returned SMILES cannot be parsed by RDKit")
    return errors


def parse_response(payload: bytes, requested_cids: Sequence[int]) -> tuple[dict[int, dict[str, Any]], dict[int, list[str]]]:
    try:
        document = json.loads(payload)
        rows = document["PropertyTable"]["Properties"]
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        raise CacheError(f"invalid PubChem property response: {exc}") from exc
    if not isinstance(rows, list):
        raise CacheError("PubChem PropertyTable.Properties is not a list")

    requested = set(requested_cids)
    records: dict[int, dict[str, Any]] = {}
    invalid: dict[int, list[str]] = {}
    for row in rows:
        if not isinstance(row, dict):
            raise CacheError("PubChem property response contains a non-object record")
        cid_value = row.get("CID")
        errors = validate_record(row, requested)
        if not isinstance(cid_value, int):
            raise CacheError("PubChem property response contains a record without an integer CID")
        if cid_value in records:
            raise CacheError(f"PubChem returned duplicate CID {cid_value}")
        if cid_value not in requested:
            raise CacheError(f"PubChem returned unrequested CID {cid_value}")
        records[cid_value] = row
        if errors:
            invalid[cid_value] = errors
    return records, invalid


def _retry_after_seconds(headers: Any) -> float | None:
    raw = headers.get("Retry-After") if headers is not None else None
    if raw is None:
        return None
    try:
        return max(0.0, float(raw))
    except (TypeError, ValueError):
        return None


def fetch_batch(
    cids: Sequence[int],
    *,
    timeout: float,
    retries: int,
    limiter: RateLimiter,
    user_agent: str,
) -> tuple[bytes, dict[str, str], int]:
    body = urllib.parse.urlencode({"cid": ",".join(str(cid) for cid in cids)}).encode("ascii")
    request = urllib.request.Request(
        _property_url(),
        data=body,
        method="POST",
        headers={
            "Accept": "application/json",
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": user_agent,
        },
    )
    for attempt in range(retries + 1):
        limiter.wait()
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = response.read()
                headers = {
                    key.lower(): value
                    for key, value in response.headers.items()
                    if key.lower() in {"content-type", "date", "server", "x-throttling-control"}
                }
                return payload, headers, int(response.status)
        except urllib.error.HTTPError as exc:
            server_delay = _cid_http_retry_delay(exc, attempt, retries)
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt >= retries:
                raise CacheError(f"PubChem request failed after {attempt + 1} attempts: {exc}") from exc
            server_delay = None
        delay = server_delay if server_delay is not None else min(30.0, 0.75 * (2**attempt))
        time.sleep(delay + random.uniform(0.0, min(0.25, delay / 4.0)))
    raise AssertionError("retry loop did not return or raise")


def _cid_http_retry_delay(exc: urllib.error.HTTPError, attempt: int, retries: int) -> float | None:
    if exc.code in RETRYABLE_HTTP_CODES and attempt < retries:
        return _retry_after_seconds(exc.headers)
    detail = exc.read(500).decode("utf-8", errors="replace")
    raise CacheError(f"PubChem HTTP {exc.code}: {detail}") from exc


def _batch_paths(cache_dir: Path, batch_id: str) -> tuple[Path, Path]:
    batch_dir = cache_dir / "batches"
    return batch_dir / f"{batch_id}.json", batch_dir / f"{batch_id}.metadata.json"


def write_batch(
    cache_dir: Path,
    cids: Sequence[int],
    payload: bytes,
    records: dict[int, dict[str, Any]],
    invalid: dict[int, list[str]],
    response_headers: dict[str, str],
    http_status: int,
) -> tuple[str, str]:
    batch_id = _batch_id(cids)
    content_hash = _sha256(payload)
    raw_path, metadata_path = _batch_paths(cache_dir, batch_id)
    returned = sorted(records)
    metadata = {
        "schema": CACHE_SCHEMA,
        "batch_id": batch_id,
        "source": SOURCE_NAME,
        "source_url": _property_url(),
        "request_method": "POST",
        "requested_cids": list(cids),
        "returned_cids": returned,
        "missing_cids": sorted(set(cids) - set(returned)),
        "invalid_structures": {str(cid): errors for cid, errors in sorted(invalid.items())},
        "properties": list(PROPERTIES),
        "retrieved_at_utc": _utc_now(),
        "http_status": http_status,
        "response_headers": response_headers,
        "content_file": raw_path.name,
        "content_bytes": len(payload),
        "content_sha256": content_hash,
    }
    _atomic_write(raw_path, payload)
    _atomic_write(metadata_path, _canonical_json_bytes(metadata))
    return batch_id, content_hash


def load_cached_records(cache_dir: Path) -> dict[int, CachedRecord]:
    records: dict[int, CachedRecord] = {}
    batch_dir = cache_dir / "batches"
    if not batch_dir.exists():
        return records
    for metadata_path in sorted(batch_dir.glob("*.metadata.json")):
        parsed, expected_batch_id, actual_hash = _load_cached_batch(metadata_path, batch_dir)
        for cid, record in parsed.items():
            candidate = CachedRecord(record, expected_batch_id, actual_hash)
            previous = records.get(cid)
            if previous is not None and previous.record != candidate.record:
                raise CacheError(f"conflicting cached PubChem records for CID {cid}")
            records[cid] = candidate
    return records


def _load_cached_batch(
    metadata_path: Path, batch_dir: Path
) -> tuple[dict[int, dict[str, Any]], str, str]:
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CacheError(f"cannot read cache metadata {metadata_path}: {exc}") from exc
    requested, expected_batch_id = _validated_batch_identity(metadata, metadata_path)
    raw_path = batch_dir / str(metadata.get("content_file", ""))
    payload = raw_path.read_bytes()
    actual_hash = _sha256(payload)
    if actual_hash != metadata.get("content_sha256"):
        raise CacheError(f"SHA-256 mismatch for cached response {raw_path}")
    parsed, invalid = parse_response(payload, requested)
    current_invalid = {str(cid): errors for cid, errors in sorted(invalid.items())}
    if metadata.get("invalid_structures", {}) != current_invalid:
        raise CacheError(f"structure validation metadata mismatch in {metadata_path}")
    if sorted(parsed) != metadata.get("returned_cids"):
        raise CacheError(f"returned CID coverage mismatch in {metadata_path}")
    return parsed, expected_batch_id, actual_hash


def _validated_batch_identity(metadata: Mapping[str, Any], metadata_path: Path) -> tuple[list[int], str]:
    if metadata.get("schema") != CACHE_SCHEMA:
        raise CacheError(f"unsupported cache schema in {metadata_path}")
    requested = metadata.get("requested_cids")
    if not isinstance(requested, list) or not all(isinstance(cid, int) for cid in requested):
        raise CacheError(f"invalid requested_cids in {metadata_path}")
    expected_batch_id = _batch_id(requested)
    if metadata.get("batch_id") != expected_batch_id:
        raise CacheError(f"request identity mismatch in {metadata_path}")
    return requested, expected_batch_id


def _chunks(values: Sequence[int], size: int) -> Iterable[list[int]]:
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def write_outputs(
    cache_dir: Path,
    requested: Sequence[int],
    records: dict[int, CachedRecord],
    invalid: dict[int, list[str]],
    *,
    offline: bool,
    fetched_batches: int,
    reused_count: int,
) -> dict[str, Any]:
    index_lines: list[bytes] = []
    for cid in sorted(set(requested) & records.keys()):
        cached = records[cid]
        item = dict(cached.record)
        item["_source_batch_id"] = cached.batch_id
        item["_source_content_sha256"] = cached.content_sha256
        index_lines.append(json.dumps(item, sort_keys=True, separators=(",", ":")).encode("utf-8") + b"\n")
    index_payload = b"".join(index_lines)
    index_path = cache_dir / "records.jsonl"
    _atomic_write(index_path, index_payload)

    returned = sorted(set(requested) & records.keys())
    missing = sorted(set(requested) - set(returned))
    manifest = {
        "schema": CACHE_SCHEMA,
        "created_at_utc": _utc_now(),
        "source": SOURCE_NAME,
        "source_url": _property_url(),
        "offline": offline,
        "requested_count": len(requested),
        "returned_count": len(returned),
        "missing_count": len(missing),
        "invalid_structure_count": len(invalid),
        "requested_cids": list(requested),
        "missing_cids": missing,
        "invalid_structures": {str(cid): errors for cid, errors in sorted(invalid.items())},
        "reused_record_count": reused_count,
        "fetched_batch_count": fetched_batches,
        "records_file": index_path.name,
        "records_bytes": len(index_payload),
        "records_sha256": _sha256(index_payload),
    }
    _atomic_write(cache_dir / "manifest.json", _canonical_json_bytes(manifest))
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", action="append", type=Path, default=[], help="CSV containing PubChem CIDs; repeatable")
    parser.add_argument("--cid-column", default="pubchem_cid_resolved", help="CID column in each input CSV")
    parser.add_argument("--cid", action="append", default=[], help="CID or delimited CID list; repeatable")
    parser.add_argument("--cache-dir", type=Path, required=True, help="Output/cache directory")
    parser.add_argument("--batch-size", type=int, default=100, help="CIDs per POST request (default: 100; maximum: 500)")
    parser.add_argument("--requests-per-second", type=float, default=3.0, help="Maximum request start rate (default: 3)")
    parser.add_argument("--timeout", type=float, default=60.0, help="Per-request timeout in seconds")
    parser.add_argument("--retries", type=int, default=5, help="Retries for transient failures")
    parser.add_argument("--contact-email", help="Optional contact included in the HTTP User-Agent")
    parser.add_argument("--offline", action="store_true", help="Forbid network access and use verified cache content only")
    parser.add_argument("--refresh", action="store_true", help="Refetch all requested CIDs (incompatible with --offline)")
    parser.add_argument("--allow-partial", action="store_true", help="Exit successfully while recording missing or invalid records")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not 1 <= args.batch_size <= 500:
        raise SystemExit("--batch-size must be between 1 and 500")
    if args.retries < 0:
        raise SystemExit("--retries cannot be negative")
    if args.offline and args.refresh:
        raise SystemExit("--offline and --refresh cannot be combined")
    try:
        manifest, fetched_batches, cache_dir = _run_cid_cache(args)
    except (CacheError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(
        f"PubChem cache: requested={manifest['requested_count']} "
        f"returned={manifest['returned_count']} missing={manifest['missing_count']} "
        f"invalid={manifest['invalid_structure_count']} fetched_batches={fetched_batches} "
        f"cache_dir={cache_dir}"
    )
    incomplete = bool(manifest["missing_count"] or manifest["invalid_structure_count"])
    return 0 if args.allow_partial or not incomplete else 2


def _run_cid_cache(args: argparse.Namespace) -> tuple[dict[str, Any], int, Path]:
    requested = read_cids(args.input_csv, args.cid_column, args.cid)
    cache_dir = args.cache_dir.expanduser().resolve()
    cached = load_cached_records(cache_dir)
    reused = {} if args.refresh else {cid: cached[cid] for cid in requested if cid in cached}
    missing = list(requested) if args.refresh else sorted(set(requested) - reused.keys())
    if args.offline and missing:
        print(f"offline cache is missing {len(missing)} of {len(requested)} requested CIDs", file=sys.stderr)
    all_records = dict(reused)
    fetched_batches = _fetch_cid_batches(args, cache_dir, missing, all_records)
    requested_set = set(requested)
    invalid = {
        cid: errors
        for cid, cached_record in all_records.items()
        if (errors := validate_record(cached_record.record, requested_set))
    }
    manifest = write_outputs(
        cache_dir,
        requested,
        all_records,
        invalid,
        offline=args.offline,
        fetched_batches=fetched_batches,
        reused_count=len(reused),
    )
    return manifest, fetched_batches, cache_dir


def _fetch_cid_batches(
    args: argparse.Namespace,
    cache_dir: Path,
    missing: Sequence[int],
    records: dict[int, CachedRecord],
) -> int:
    if args.offline:
        return 0
    limiter = RateLimiter(args.requests_per_second)
    contact = f"; contact={args.contact_email}" if args.contact_email else ""
    fetched = 0
    for batch in _chunks(missing, args.batch_size):
        payload, headers, status = fetch_batch(
            batch,
            timeout=args.timeout,
            retries=args.retries,
            limiter=limiter,
            user_agent=f"atlas-fda-pubchem-cache/1.0{contact}",
        )
        parsed, invalid = parse_response(payload, batch)
        batch_id, content_hash = write_batch(cache_dir, batch, payload, parsed, invalid, headers, status)
        records.update({cid: CachedRecord(record, batch_id, content_hash) for cid, record in parsed.items()})
        fetched += 1
    return fetched


if __name__ == "__main__":
    raise SystemExit(main())
