#!/usr/bin/env python3
"""Independently verify a terminal FDA legacy reconciliation bundle.

The verifier accepts artifact paths and schema columns at the CLI boundary and
does not import resolver code.  It always writes a machine-readable report and
returns nonzero if any invariant cannot be established.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


TERMINAL_DISPOSITIONS = (
    "confirmed_fda_active_ingredient",
    "confirmed_fda_salt_parent_docked",
    "retained_salt_counterion",
    "collapsed_combination_product",
    "additive_excipient",
    "identified_non_fda",
    "incorrect_name_structure_mapping_resolved",
    "prepared_file_unusable",
)
FDA_REGULATORY_STATUSES = (
    "fda_approved_current_or_historical",
    "drugcentral_fda_approved",
)
STRUCTURE_EVIDENCE_CODES = (
    "manifest_approved_structure",
    "drugcentral_exact",
    "drugcentral_parent",
    "pubchem_exact",
    "pubchem_parent",
    "pubchem_validated_connectivity_name",
)
USABLE_IDENTITY_EVIDENCE_CODES = STRUCTURE_EVIDENCE_CODES + (
    "mapping_exact_stereoisomer_pdbqt_lineage",
    "mapping_exact_structure_pdbqt_lineage",
)
FDA_SOURCE_CODES = (
    "drugcentral_approved_manifest",
    "drugcentral_approved_id",
    "drugsatfda_active_ingredient",
    "drugsatfda_salt_core",
)
DEFAULT_USABLE_STATUSES = ("selected_unique_bytes", "preserved_usable", "ready", "usable")
UNSAFE_CANONICAL_STATUS = "quarantined_unsafe_parent_collapse"
UNSAFE_REASON_CODE = "unsafe_parent_collapse"
UNSAFE_MAPPING_REASON_CODE = "unsafe_parent_collapse"
UNSAFE_SOURCE_MARKER = "repair_quarantine_source_record_linked"
CODE_SPLIT_RE = re.compile(r"[\s,;|]+")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
MAX_FAILURE_SAMPLES = 20


@dataclass
class Check:
    name: str
    passed: bool = True
    metrics: dict[str, Any] = field(default_factory=dict)
    failures: list[str] = field(default_factory=list)
    failure_count: int = 0

    def fail(self, message: str) -> None:
        self.passed = False
        self.failure_count += 1
        if len(self.failures) < MAX_FAILURE_SAMPLES:
            self.failures.append(message)

    def as_dict(self) -> dict[str, Any]:
        metrics = dict(self.metrics)
        metrics["failure_count"] = self.failure_count
        return {
            "name": self.name,
            "passed": self.passed,
            "metrics": metrics,
            "failure_samples": self.failures,
            "failure_samples_truncated": self.failure_count > len(self.failures),
        }


class VerificationError(RuntimeError):
    """Raised when an artifact cannot be parsed sufficiently to continue."""


class HashCache:
    def __init__(self) -> None:
        self._values: dict[Path, str] = {}

    def sha256(self, path: Path) -> str:
        resolved = path.resolve()
        if resolved not in self._values:
            digest = hashlib.sha256()
            with resolved.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            self._values[resolved] = digest.hexdigest()
        return self._values[resolved]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _atomic_json(path: Path, value: Any) -> None:
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise VerificationError(f"CSV has no header: {path}")
            return list(reader.fieldnames), [dict(row) for row in reader]
    except (OSError, csv.Error) as exc:
        raise VerificationError(f"cannot read CSV {path}: {exc}") from exc


def _require_columns(path: Path, fields: Sequence[str], required: Iterable[str]) -> None:
    missing = sorted(set(required) - set(fields))
    if missing:
        raise VerificationError(f"{path} is missing required columns: {', '.join(missing)}")


def _clean(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def _percent(numerator: int, denominator: int) -> float:
    return round(100.0 * numerator / denominator, 6) if denominator else 100.0


def _codes(value: Any) -> set[str]:
    return {token for token in CODE_SPLIT_RE.split(_clean(value)) if token}


def _row_key(row: Mapping[str, Any], row_column: str, rdk_column: str) -> tuple[str, str]:
    return _clean(row.get(row_column)), _clean(row.get(rdk_column))


def _resolve_artifact_path(raw: str, table_path: Path) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else table_path.parent / path


def _validate_file_hash(
    check: Check,
    hashes: HashCache,
    raw_path: str,
    expected_hash: str,
    table_path: Path,
    context: str,
) -> None:
    path_text = _clean(raw_path)
    expected = _clean(expected_hash).lower()
    if not path_text or not expected:
        check.fail(f"{context}: path/hash pair is incomplete")
        return
    if not SHA256_RE.fullmatch(expected):
        check.fail(f"{context}: malformed SHA-256 {expected!r}")
        return
    path = _resolve_artifact_path(path_text, table_path)
    if not path.is_file():
        check.fail(f"{context}: file does not exist: {path}")
        return
    try:
        actual = hashes.sha256(path)
    except OSError as exc:
        check.fail(f"{context}: cannot hash {path}: {exc}")
        return
    if actual != expected:
        check.fail(f"{context}: SHA-256 mismatch for {path}; expected {expected}, got {actual}")


def check_mapping_uniqueness(
    rows: Sequence[Mapping[str, str]], row_column: str, rdk_column: str
) -> Check:
    check = Check("mapping_row_and_rdk_uniqueness")
    row_ids = [_clean(row.get(row_column)) for row in rows]
    rdk_ids = [_clean(row.get(rdk_column)) for row in rows]
    row_counts = Counter(row_ids)
    rdk_counts = Counter(rdk_ids)
    _check_unique_counts(check, row_column, row_counts)
    _check_unique_counts(check, rdk_column, rdk_counts)
    check.metrics = {
        "mapping_rows": len(rows),
        "unique_mapping_rows": len(set(row_ids) - {""}),
        "unique_rdk_ids": len(set(rdk_ids) - {""}),
        "failure_count": sum(1 for value in row_counts if not value or row_counts[value] != 1)
        + sum(1 for value in rdk_counts if not value or rdk_counts[value] != 1),
    }
    return check


def _check_unique_counts(check: Check, label: str, counts: Counter[str]) -> None:
    for value, count in counts.items():
        if not value:
            check.fail(f"blank {label} occurs {count} time(s)")
        elif count != 1:
            check.fail(f"duplicate {label}={value!r} occurs {count} times")


def check_prepared_files(
    rows: Sequence[Mapping[str, str]],
    mapping_path: Path,
    args: argparse.Namespace,
    hashes: HashCache,
) -> Check:
    check = Check("prepared_path_sha_parity")
    checked = 0
    incomplete = 0
    resolved_paths = [
        _resolve_artifact_path(value, mapping_path).resolve()
        for row in rows
        if (value := _clean(row.get(args.path_column)))
    ]
    path_counts = Counter(resolved_paths)
    for index, row in enumerate(rows, start=2):
        path_value = _clean(row.get(args.path_column))
        hash_value = _clean(row.get(args.sha_column))
        if not path_value and not hash_value:
            continue
        checked += 1
        _validate_file_hash(check, hashes, path_value, hash_value, mapping_path, f"mapping row {index}")
        if not path_value or not hash_value:
            incomplete += 1
        if path_value:
            path = _resolve_artifact_path(path_value, mapping_path).resolve()
            if path_counts[path] != 1:
                check.fail(f"mapping row {index}: prepared path occurs {path_counts[path]} times: {path}")
    check.metrics = {
        "mapping_rows": len(rows),
        "rows_with_prepared_artifact": checked,
        "prepared_artifact_coverage_percent": _percent(checked, len(rows)),
        "incomplete_path_hash_pairs": incomplete,
        "failure_count": max(incomplete, len(check.failures)),
    }
    return check


def check_regression_sentinel(rows: Sequence[Mapping[str, str]], args: argparse.Namespace) -> Check:
    check = Check("deserpidine_regression_sentinel")
    matches = [
        row
        for row in rows
        if _clean(row.get(args.row_column)) == args.regression_mapping_row
        and _clean(row.get(args.rdk_column)) == args.regression_rdk_id
    ]
    if len(matches) != 1:
        check.fail(
            f"expected exactly one sentinel ({args.regression_mapping_row}, {args.regression_rdk_id}), found {len(matches)}"
        )
    else:
        _validate_regression_sentinel(check, matches[0], args)
    check.metrics = {
        "matching_rows": len(matches),
        "expected_mapping_row": args.regression_mapping_row,
        "expected_rdk_id": args.regression_rdk_id,
        "expected_preferred_name": args.regression_preferred_name,
        "expected_selected_sha256": args.regression_selected_sha256,
    }
    return check


def _validate_regression_sentinel(
    check: Check, row: Mapping[str, str], args: argparse.Namespace
) -> None:
    preferred = _clean(row.get(args.preferred_name_column))
    if preferred.casefold() != args.regression_preferred_name.casefold():
        check.fail(f"sentinel preferred name is {preferred!r}, expected {args.regression_preferred_name!r}")
    regulatory = _clean(row.get(args.regulatory_column))
    if regulatory not in set(args.fda_regulatory_status):
        check.fail(f"sentinel is not FDA: regulatory_status={regulatory!r}")
    actual_sha = _clean(row.get(args.sha_column)).lower()
    if actual_sha != args.regression_selected_sha256.lower():
        check.fail(f"sentinel selected SHA-256 is {actual_sha!r}, expected {args.regression_selected_sha256!r}")
    disposition = _clean(row.get(args.disposition_column))
    if disposition not in set(args.terminal_value) or disposition in {
        "identified_non_fda",
        "prepared_file_unusable",
    }:
        check.fail(f"sentinel has unacceptable terminal disposition {disposition!r}")
    structure_hits = _codes(row.get(args.identity_evidence_column)) & set(args.structure_evidence_code)
    approval_hits = _codes(row.get(args.approval_evidence_column)) & set(args.fda_source_code)
    if not structure_hits or not approval_hits:
        check.fail("sentinel lacks required structure and FDA approval evidence")


def check_terminal_partition(
    mapping_rows: Sequence[Mapping[str, str]],
    terminal_rows: Sequence[Mapping[str, str]],
    args: argparse.Namespace,
) -> tuple[Check, dict[tuple[str, str], Mapping[str, str]]]:
    check = Check("terminal_disposition_partition")
    mapping_keys = [_row_key(row, args.row_column, args.rdk_column) for row in mapping_rows]
    terminal_keys = [_row_key(row, args.row_column, args.rdk_column) for row in terminal_rows]
    terminal_counts = Counter(terminal_keys)
    terminal_by_key: dict[tuple[str, str], Mapping[str, str]] = {}
    for key, row in zip(terminal_keys, terminal_rows):
        if key not in terminal_by_key:
            terminal_by_key[key] = row
    _check_terminal_key_counts(check, terminal_counts)
    missing = sorted(set(mapping_keys) - set(terminal_keys))
    extra = sorted(set(terminal_keys) - set(mapping_keys))
    _report_key_difference(check, missing, "mapping key absent from terminal partition")
    _report_key_difference(check, extra, "terminal partition contains unknown mapping key")

    allowed = set(args.terminal_value)
    disposition_counts: Counter[str] = Counter()
    _check_terminal_dispositions(check, terminal_by_key, args.disposition_column, allowed, disposition_counts)
    parity_columns = (
        args.resolution_column,
        args.disposition_column,
        args.regulatory_column,
        args.path_column,
        args.sha_column,
        args.file_status_column,
    )
    mapping_by_key = {_row_key(row, args.row_column, args.rdk_column): row for row in mapping_rows}
    _check_terminal_parity(check, mapping_by_key, terminal_by_key, parity_columns)
    check.metrics = {
        "mapping_rows": len(mapping_rows),
        "terminal_rows": len(terminal_rows),
        "terminal_partition_coverage_percent": _percent(len(set(mapping_keys) & set(terminal_keys)), len(set(mapping_keys))),
        "missing_partition_rows": len(missing),
        "extra_partition_rows": len(extra),
        "disposition_counts": dict(sorted(disposition_counts.items())),
        "failure_count": len(check.failures),
    }
    return check, terminal_by_key


def _check_terminal_key_counts(check: Check, counts: Counter[tuple[str, str]]) -> None:
    for key, count in counts.items():
        if count != 1:
            check.fail(f"terminal key {key!r} occurs {count} times")


def _report_key_difference(check: Check, keys: Sequence[tuple[str, str]], message: str) -> None:
    for key in keys:
        check.fail(f"{message}: {key!r}")


def _check_terminal_dispositions(
    check: Check,
    rows: Mapping[tuple[str, str], Mapping[str, str]],
    column: str,
    allowed: set[str],
    counts: Counter[str],
) -> None:
    for key, row in rows.items():
        disposition = _clean(row.get(column))
        counts[disposition] += 1
        if disposition not in allowed:
            check.fail(f"terminal key {key!r} has unexpected disposition {disposition!r}")


def _check_terminal_parity(
    check: Check,
    mapping_rows: Mapping[tuple[str, str], Mapping[str, str]],
    terminal_rows: Mapping[tuple[str, str], Mapping[str, str]],
    columns: Sequence[str],
) -> None:
    for key in set(mapping_rows) & set(terminal_rows):
        for column in columns:
            if column in mapping_rows[key] and column in terminal_rows[key]:
                if _clean(mapping_rows[key][column]) != _clean(terminal_rows[key][column]):
                    check.fail(f"key {key!r}: {column} differs between mapping and terminal table")


def check_no_nonterminal_usable(
    mapping_rows: Sequence[Mapping[str, str]],
    terminal_by_key: Mapping[tuple[str, str], Mapping[str, str]],
    args: argparse.Namespace,
) -> Check:
    check = Check("zero_nonterminal_usable_rows")
    usable_count = 0
    nonterminal_count = 0
    unusable_count = 0
    for row in mapping_rows:
        classification, error = _classify_usable_row(row, terminal_by_key, args)
        if classification == "unusable":
            unusable_count += 1
            continue
        if classification == "absent":
            continue
        usable_count += 1
        if error:
            nonterminal_count += 1
            check.fail(error)
    if usable_count != args.expected_prepared_usable:
        check.fail(f"prepared usable count is {usable_count}, expected {args.expected_prepared_usable}")
    if unusable_count != args.expected_prepared_unusable:
        check.fail(f"prepared unusable count is {unusable_count}, expected {args.expected_prepared_unusable}")
    check.metrics = {
        "usable_rows": usable_count,
        "expected_usable_rows": args.expected_prepared_usable,
        "prepared_file_unusable_rows": unusable_count,
        "expected_prepared_file_unusable_rows": args.expected_prepared_unusable,
        "usable_terminal_identity_coverage_percent": _percent(usable_count - nonterminal_count, usable_count),
        "nonterminal_usable_rows": nonterminal_count,
        "failure_count": nonterminal_count,
    }
    return check


def _classify_usable_row(
    row: Mapping[str, str],
    terminal_by_key: Mapping[tuple[str, str], Mapping[str, str]],
    args: argparse.Namespace,
) -> tuple[str, str]:
    key = _row_key(row, args.row_column, args.rdk_column)
    terminal = terminal_by_key.get(key)
    disposition = _clean(terminal.get(args.disposition_column)) if terminal else ""
    if disposition == "prepared_file_unusable":
        return "unusable", ""
    has_artifact = bool(_clean(row.get(args.path_column)) or _clean(row.get(args.sha_column)))
    if not has_artifact and _clean(row.get(args.file_status_column)) not in set(args.usable_status):
        return "absent", ""
    resolution = _clean(terminal.get(args.resolution_column)) if terminal else ""
    preferred = _clean(row.get(args.preferred_name_column))
    identity_hits = _codes(row.get(args.identity_evidence_column)) & set(args.usable_identity_evidence_code)
    if _usable_identity_complete(disposition, resolution, preferred, identity_hits, args):
        return "usable", ""
    error = (
        f"usable key {key!r} lacks terminal identity: resolution={resolution!r}, "
        f"disposition={disposition!r}, preferred_name={preferred!r}, "
        f"recognized_identity_evidence={sorted(identity_hits)}"
    )
    return "usable", error


def _usable_identity_complete(
    disposition: str,
    resolution: str,
    preferred: str,
    identity_hits: set[str],
    args: argparse.Namespace,
) -> bool:
    return bool(disposition in set(args.terminal_value) and resolution and preferred and identity_hits)


def check_fda_evidence(rows: Sequence[Mapping[str, str]], args: argparse.Namespace) -> Check:
    check = Check("fda_rows_have_structure_and_approval_evidence")
    fda_statuses = set(args.fda_regulatory_status)
    structure_allowlist = set(args.structure_evidence_code)
    source_allowlist = set(args.fda_source_code)
    fda_rows = 0
    failures = 0
    for index, row in enumerate(rows, start=2):
        status = _clean(row.get(args.regulatory_column))
        if status not in fda_statuses:
            continue
        fda_rows += 1
        identity_codes = _codes(row.get(args.identity_evidence_column))
        approval_codes = _codes(row.get(args.approval_evidence_column))
        structure_hits = sorted(identity_codes & structure_allowlist)
        source_hits = sorted(approval_codes & source_allowlist)
        if not structure_hits or not source_hits:
            failures += 1
            check.fail(
                f"mapping row {index}: FDA status {status!r}, structure hits={structure_hits}, approval hits={source_hits}"
            )
    check.metrics = {
        "fda_rows": fda_rows,
        "fda_required_evidence_coverage_percent": _percent(fda_rows - failures, fda_rows),
        "fda_rows_without_required_evidence": failures,
        "failure_count": failures,
    }
    return check


def check_fda_canonical_parents(rows: Sequence[Mapping[str, str]], args: argparse.Namespace) -> Check:
    check = Check("fda_canonical_parent_readiness")
    fda_statuses = set(args.fda_regulatory_status)
    allowed_statuses = set(args.canonical_status)
    ready_statuses = set(args.canonical_ready_status)
    fda_rows = 0
    valid_rows = 0
    hashes = HashCache()
    for index, row in enumerate(rows, start=2):
        if _clean(row.get(args.regulatory_column)) not in fda_statuses:
            continue
        fda_rows += 1
        before = check.failure_count
        parent_id = _clean(row.get(args.canonical_parent_id_column))
        status = _clean(row.get(args.canonical_status_column))
        if not parent_id:
            check.fail(f"mapping row {index}: FDA row lacks canonical_parent_id")
        if status not in allowed_statuses:
            check.fail(f"mapping row {index}: unexpected canonical parent status {status!r}")
        if status in ready_statuses:
            _validate_file_hash(
                check,
                hashes,
                _clean(row.get(args.canonical_path_column)),
                _clean(row.get(args.canonical_sha_column)),
                args.mapping.expanduser().resolve(),
                f"mapping row {index} canonical parent",
            )
        if check.failure_count == before:
            valid_rows += 1
    check.metrics = {
        "fda_rows": fda_rows,
        "valid_canonical_parent_rows": valid_rows,
        "canonical_parent_coverage_percent": _percent(valid_rows, fda_rows),
    }
    return check


def check_unsafe_parent_quarantine(
    rows: Sequence[Mapping[str, str]],
    repair_path: Path | None,
    repair_fields: Sequence[str],
    repair_rows: Sequence[Mapping[str, str]],
    args: argparse.Namespace,
) -> Check:
    check = Check("unsafe_parent_collapse_repair_quarantine")
    unsafe_rows = [
        (index, row)
        for index, row in enumerate(rows, start=2)
        if _clean(row.get(args.canonical_status_column)) == args.unsafe_canonical_status
    ]
    check.metrics = {
        "unsafe_mapping_rows": len(unsafe_rows),
        "repair_quarantine_path": str(repair_path) if repair_path else "",
        "matched_source_rows": 0,
        "matched_source_row_hashes": [],
    }
    if not unsafe_rows:
        return check
    if repair_path is None:
        check.fail("unsafe canonical rows exist but --repair-quarantine was not supplied")
        return check
    required = {
        args.repair_scope_column,
        args.repair_source_index_column,
        args.repair_parent_column,
        args.repair_reason_column,
    }
    missing = sorted(required - set(repair_fields))
    if missing:
        check.fail(f"repair quarantine is missing columns: {', '.join(missing)}")
        return check

    matched = 0
    row_hashes: list[dict[str, str]] = []
    for index, row in unsafe_rows:
        parent_key, source_index, context = _validate_unsafe_mapping_row(check, row, index, args)
        candidates = _matching_repair_rows(repair_rows, parent_key, source_index, args)
        if len(candidates) != 1:
            check.fail(
                f"{context}: expected one repair-quarantine source row for source={source_index!r}, "
                f"parent={parent_key!r}; found {len(candidates)}"
            )
            continue
        matched += 1
        row_payload = json.dumps(candidates[0], sort_keys=True, separators=(",", ":")).encode("utf-8")
        row_hashes.append(
            {
                "mapping_row_number": _clean(row.get(args.row_column)),
                "source_record_index": source_index,
                "parent_inchikey": parent_key,
                "repair_source_row_sha256": hashlib.sha256(row_payload).hexdigest(),
            }
        )
    check.metrics["matched_source_rows"] = matched
    check.metrics["matched_source_row_hashes"] = row_hashes
    check.metrics["unsafe_source_coverage_percent"] = _percent(matched, len(unsafe_rows))
    return check


def _validate_unsafe_mapping_row(
    check: Check, row: Mapping[str, str], index: int, args: argparse.Namespace
) -> tuple[str, str, str]:
    parent_id = _clean(row.get(args.canonical_parent_id_column))
    parent_key = _clean(row.get(args.identity_parent_column)).upper()
    source_index = _clean(row.get(args.approved_source_index_column))
    context = f"mapping row {index} ({_clean(row.get(args.rdk_column))})"
    if not parent_id or not parent_key or parent_key not in parent_id.upper():
        check.fail(f"{context}: synthetic canonical_parent_id does not encode identity parent {parent_key!r}")
    if _clean(row.get(args.canonical_path_column)) or _clean(row.get(args.canonical_sha_column)):
        check.fail(f"{context}: unsafe canonical quarantine must have blank PDBQT path and SHA-256")
    if not source_index:
        check.fail(f"{context}: approved source record index is blank")
    mapping_reasons = _codes(row.get(args.terminal_reason_column))
    required_reasons = {args.unsafe_mapping_reason_code, args.unsafe_source_marker}
    if not required_reasons <= mapping_reasons:
        check.fail(f"{context}: terminal reasons lack {sorted(required_reasons - mapping_reasons)}")
    return parent_key, source_index, context


def _matching_repair_rows(
    rows: Sequence[Mapping[str, str]], parent_key: str, source_index: str, args: argparse.Namespace
) -> list[Mapping[str, str]]:
    return [
        row
        for row in rows
        if _clean(row.get(args.repair_scope_column)) == args.repair_source_scope
        and source_index in _codes(row.get(args.repair_source_index_column))
        and _clean(row.get(args.repair_parent_column)).upper() == parent_key
        and args.unsafe_reason_code in _codes(row.get(args.repair_reason_column))
    ]


def check_prepared_connectivity_partition(
    rows: Sequence[Mapping[str, str]], args: argparse.Namespace, hashes: HashCache
) -> Check:
    check = Check("prepared_connectivity_and_score_reuse_partition")
    expected_pairs = {
        "coordinate_connectivity_exact": ("reuse_supported_coordinate_exact", args.expected_connectivity_exact),
        "coordinate_connectivity_parent": ("reuse_supported_coordinate_parent", args.expected_connectivity_parent),
        "ambiguous": ("no_reuse_ambiguous_connectivity", args.expected_connectivity_ambiguous),
        "unverifiable": ("no_reuse_missing_file", args.expected_connectivity_unverifiable),
    }
    counts: Counter[str] = Counter()
    reuse_counts: Counter[str] = Counter()
    fda_ambiguous = 0
    delafloxacin_matches = 0
    mapping_path = args.mapping.expanduser().resolve()
    for index, row in enumerate(rows, start=2):
        state = _clean(row.get(args.connectivity_state_column))
        reuse = _clean(row.get(args.score_reuse_column))
        counts[state] += 1
        reuse_counts[reuse] += 1
        expected = expected_pairs.get(state)
        if expected is None:
            check.fail(f"mapping row {index}: unexpected prepared connectivity state {state!r}")
            continue
        fda_ambiguous += _check_connectivity_row(check, hashes, row, index, state, reuse, expected[0], mapping_path, args)
        if _clean(row.get(args.rdk_column)) == args.ambiguous_fda_rdk_id:
            delafloxacin_matches += 1
            if state != "ambiguous" or reuse != "no_reuse_ambiguous_connectivity":
                check.fail(
                    f"required ambiguous FDA sentinel {args.ambiguous_fda_rdk_id} has state={state!r}, reuse={reuse!r}"
                )
    expected_reuse_counts = {
        "reuse_supported_coordinate_exact": args.expected_reuse_coordinate_exact,
        "reuse_supported_coordinate_parent": args.expected_reuse_coordinate_parent,
        "no_reuse_ambiguous_connectivity": args.expected_connectivity_ambiguous,
        "no_reuse_missing_file": args.expected_connectivity_unverifiable,
        "no_reuse_unsafe_parent_collapse": args.expected_no_reuse_unsafe_parent,
    }
    _check_count_expectations(check, counts, {key: value[1] for key, value in expected_pairs.items()}, "prepared connectivity")
    _check_count_expectations(check, reuse_counts, expected_reuse_counts, "legacy score reuse status")
    if delafloxacin_matches != 1:
        check.fail(
            f"expected one ambiguous FDA sentinel {args.ambiguous_fda_rdk_id}, found {delafloxacin_matches}"
        )
    check.metrics = {
        "mapping_rows": len(rows),
        "connectivity_counts": dict(sorted(counts.items())),
        "score_reuse_counts": dict(sorted(reuse_counts.items())),
        "coordinate_safe_reuse_rows": counts["coordinate_connectivity_exact"]
        + counts["coordinate_connectivity_parent"],
        "score_reuse_authorized_rows": reuse_counts["reuse_supported_coordinate_exact"]
        + reuse_counts["reuse_supported_coordinate_parent"],
        "no_reuse_rows": sum(count for status, count in reuse_counts.items() if status.startswith("no_reuse_")),
        "ambiguous_fda_rows_with_required_canonical_check": fda_ambiguous,
        "partition_coverage_percent": _percent(sum(counts.values()), len(rows)),
    }
    return check


def _check_count_expectations(
    check: Check, actual: Counter[str], expected: Mapping[str, int], label: str
) -> None:
    for value, expected_count in expected.items():
        if actual[value] != expected_count:
            check.fail(f"{label} {value!r} count is {actual[value]}, expected {expected_count}")
    unexpected = sorted(set(actual) - set(expected))
    if unexpected:
        check.fail(f"unexpected {label} values: {unexpected}")


def _check_connectivity_row(
    check: Check,
    hashes: HashCache,
    row: Mapping[str, str],
    index: int,
    state: str,
    reuse: str,
    standard_reuse: str,
    mapping_path: Path,
    args: argparse.Namespace,
) -> int:
    unsafe_parent = _clean(row.get(args.canonical_status_column)) == args.unsafe_canonical_status
    expected_reuse = "no_reuse_unsafe_parent_collapse" if unsafe_parent else standard_reuse
    if unsafe_parent and state != "coordinate_connectivity_parent":
        check.fail(f"mapping row {index}: unsafe parent-collapse row must retain coordinate_connectivity_parent, got {state!r}")
    if reuse != expected_reuse:
        check.fail(f"mapping row {index}: connectivity {state!r} has incompatible reuse status {reuse!r}")
    _check_connectivity_artifact(check, row, index, state, args)
    is_fda_ambiguous = state == "ambiguous" and _clean(row.get(args.regulatory_column)) in set(args.fda_regulatory_status)
    if is_fda_ambiguous:
        _check_ambiguous_fda_replacement(check, hashes, row, index, mapping_path, args)
    return int(is_fda_ambiguous)


def _check_connectivity_artifact(
    check: Check, row: Mapping[str, str], index: int, state: str, args: argparse.Namespace
) -> None:
    has_artifact = bool(_clean(row.get(args.path_column)) and _clean(row.get(args.sha_column)))
    if state != "unverifiable" and not has_artifact:
        check.fail(f"mapping row {index}: {state} row lacks prepared path/SHA")
    if state == "unverifiable" and has_artifact:
        check.fail(f"mapping row {index}: unverifiable row unexpectedly has prepared path/SHA")


def _check_ambiguous_fda_replacement(
    check: Check,
    hashes: HashCache,
    row: Mapping[str, str],
    index: int,
    mapping_path: Path,
    args: argparse.Namespace,
) -> None:
    if _clean(row.get(args.canonical_status_column)) != "ready":
        check.fail(f"mapping row {index}: ambiguous FDA row lacks ready canonical replacement")
        return
    _validate_file_hash(
        check,
        hashes,
        _clean(row.get(args.canonical_path_column)),
        _clean(row.get(args.canonical_sha_column)),
        mapping_path,
        f"mapping row {index} ambiguous FDA canonical replacement",
    )


def _integer_tokens(value: Any, context: str) -> set[int]:
    tokens = _codes(value)
    if not all(token.isdigit() and int(token) > 0 for token in tokens):
        raise VerificationError(f"{context} contains invalid source-record indices: {_clean(value)!r}")
    return {int(token) for token in tokens}


def _count_sdf_records(path: Path) -> int:
    count = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.rstrip("\r\n") == "$$$$":
                    count += 1
    except OSError as exc:
        raise VerificationError(f"cannot read FDA source SDF {path}: {exc}") from exc
    return count


def check_fda_source_partition(
    named_rows: Sequence[Mapping[str, str]],
    repair_rows: Sequence[Mapping[str, str]],
    source_sdf: Path | None,
    args: argparse.Namespace,
) -> Check:
    check = Check("fda_source_record_partition")
    if source_sdf is None:
        if repair_rows:
            check.fail("repair quarantine was supplied without --fda-source-sdf")
        check.metrics = {"skipped": True, "reason": "fda_source_sdf_not_supplied"}
        return check
    named_indices = _collect_named_source_indices(check, named_rows, args)
    repair_indices, repair_source_rows = _collect_repair_source_indices(check, repair_rows, args)
    record_count = _count_sdf_records(source_sdf)
    expected_universe = set(range(1, record_count + 1))
    if record_count != args.expected_fda_source_records:
        check.fail(f"FDA source SDF contains {record_count} records, expected {args.expected_fda_source_records}")
    if len(named_indices) != args.expected_named_source_records:
        check.fail(f"named manifest covers {len(named_indices)} source records, expected {args.expected_named_source_records}")
    if len(repair_indices) != args.expected_repair_source_records:
        check.fail(f"repair quarantine covers {len(repair_indices)} source records, expected {args.expected_repair_source_records}")
    if named_indices & repair_indices:
        check.fail(f"named and repair-quarantine source partitions overlap at {sorted(named_indices & repair_indices)}")
    union = named_indices | repair_indices
    if union != expected_universe:
        check.fail(
            f"FDA source partition incomplete: missing={len(expected_universe - union)}, extra={len(union - expected_universe)}"
        )
    check.metrics = {
        "source_sdf_records": record_count,
        "named_source_records": len(named_indices),
        "repair_quarantine_source_records": len(repair_indices),
        "repair_quarantine_source_rows": repair_source_rows,
        "partition_overlap_records": len(named_indices & repair_indices),
        "partition_coverage_percent": _percent(len(union & expected_universe), len(expected_universe)),
    }
    return check


def _collect_named_source_indices(
    check: Check, rows: Sequence[Mapping[str, str]], args: argparse.Namespace
) -> set[int]:
    indices: set[int] = set()
    for index, row in enumerate(rows, start=2):
        values = _integer_tokens(row.get(args.named_source_indices_column), f"named manifest row {index}")
        overlap = indices & values
        if overlap:
            check.fail(f"named manifest row {index}: duplicate source indices {sorted(overlap)}")
        indices.update(values)
    return indices


def _collect_repair_source_indices(
    check: Check, rows: Sequence[Mapping[str, str]], args: argparse.Namespace
) -> tuple[set[int], int]:
    indices: set[int] = set()
    count = 0
    for index, row in enumerate(rows, start=2):
        if _clean(row.get(args.repair_scope_column)) != args.repair_source_scope:
            continue
        count += 1
        values = _integer_tokens(row.get(args.repair_source_index_column), f"repair quarantine row {index}")
        overlap = indices & values
        if overlap:
            check.fail(f"repair quarantine row {index}: duplicate source indices {sorted(overlap)}")
        indices.update(values)
        if args.unsafe_reason_code not in _codes(row.get(args.repair_reason_column)):
            check.fail(f"repair quarantine row {index}: approved source row lacks unsafe_parent_collapse reason")
    return indices, count


def check_named_manifest(
    rows: Sequence[Mapping[str, str]],
    manifest_path: Path,
    args: argparse.Namespace,
    hashes: HashCache,
) -> Check:
    check = Check("named_library_manifest_counts_and_hashes")
    ready_values = set(args.named_ready_status)
    quarantine_values = set(args.named_quarantine_status)
    allowed = ready_values | quarantine_values
    counts: Counter[str] = Counter()
    output_paths: Counter[Path] = Counter()
    sidecar_paths: Counter[Path] = Counter()
    hash_pairs_checked = 0
    for index, row in enumerate(rows, start=2):
        hash_pairs_checked += _check_named_row(
            check, hashes, row, index, manifest_path, args, allowed, ready_values, counts, output_paths, sidecar_paths
        )
    for label, paths in (("output", output_paths), ("sidecar", sidecar_paths)):
        for path, count in paths.items():
            if count != 1:
                check.fail(f"named manifest {label} path occurs {count} times: {path}")
    if counts["ready"] != args.expected_named_ready:
        check.fail(f"named ready count is {counts['ready']}, expected {args.expected_named_ready}")
    if counts["quarantine"] != args.expected_named_quarantine:
        check.fail(f"named quarantine count is {counts['quarantine']}, expected {args.expected_named_quarantine}")
    check.metrics = {
        "manifest_rows": len(rows),
        "ready_rows": counts["ready"],
        "quarantine_rows": counts["quarantine"],
        "unexpected_rows": counts["unexpected"],
        "hash_pairs_checked": hash_pairs_checked,
        "failure_count": len(check.failures),
    }
    return check


def _check_named_row(
    check: Check,
    hashes: HashCache,
    row: Mapping[str, str],
    index: int,
    manifest_path: Path,
    args: argparse.Namespace,
    allowed: set[str],
    ready_values: set[str],
    counts: Counter[str],
    output_paths: Counter[Path],
    sidecar_paths: Counter[Path],
) -> int:
    status = _clean(row.get(args.named_status_column))
    if status not in allowed:
        check.fail(f"named manifest row {index}: unexpected materialization status {status!r}")
        counts["unexpected"] += 1
        return 0
    bucket = "ready" if status in ready_values else "quarantine"
    counts[bucket] += 1
    checked = 0
    pairs = ((args.named_path_column, args.named_sha_column, "PDBQT"), (args.sidecar_path_column, args.sidecar_sha_column, "sidecar"))
    for path_column, sha_column, label in pairs:
        path_value, sha_value = _clean(row.get(path_column)), _clean(row.get(sha_column))
        if bucket == "ready":
            checked += 1
            _validate_file_hash(check, hashes, path_value, sha_value, manifest_path, f"named manifest row {index} {label}")
        elif sha_value:
            check.fail(f"named manifest row {index} quarantined {label} has a hash for a non-materialized artifact")
        if path_value:
            target = output_paths if label == "PDBQT" else sidecar_paths
            target[_resolve_artifact_path(path_value, manifest_path).resolve()] += 1
    return checked


def _parse_pubchem_raw(path: Path) -> set[int]:
    try:
        document = json.loads(path.read_bytes())
        properties = document["PropertyTable"]["Properties"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        raise VerificationError(f"cannot parse PubChem raw batch {path}: {exc}") from exc
    if not isinstance(properties, list):
        raise VerificationError(f"PubChem raw batch has non-list Properties: {path}")
    cids: list[int] = []
    for record in properties:
        cid = record.get("CID") if isinstance(record, dict) else None
        if not isinstance(cid, int) or cid <= 0:
            raise VerificationError(f"PubChem raw batch contains invalid CID: {path}")
        cids.append(cid)
    if len(cids) != len(set(cids)):
        raise VerificationError(f"PubChem raw batch contains duplicate CIDs: {path}")
    return set(cids)


def _mapping_pubchem_cids(rows: Sequence[Mapping[str, str]], column: str) -> set[int]:
    cids: set[int] = set()
    for index, row in enumerate(rows, start=2):
        value = _clean(row.get(column))
        if not value:
            continue
        if re.fullmatch(r"[0-9]+\.0+", value):
            value = value.split(".", 1)[0]
        if not value.isdigit() or int(value) <= 0:
            raise VerificationError(f"mapping row {index} has invalid PubChem CID {value!r}")
        cids.add(int(value))
    return cids


def check_pubchem_cache(
    manifest_path: Path,
    mapping_rows: Sequence[Mapping[str, str]],
    args: argparse.Namespace,
    hashes: HashCache,
) -> Check:
    check = Check("pubchem_manifest_batch_hashes_and_offline_coverage")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot read PubChem manifest {manifest_path}: {exc}") from exc
    cache_dir = manifest_path.parent
    records_path = cache_dir / _clean(manifest.get("records_file"))
    expected_records_hash = _clean(manifest.get("records_sha256")).lower()
    _validate_file_hash(
        check,
        hashes,
        str(records_path),
        expected_records_hash,
        manifest_path,
        "PubChem normalized records",
    )
    record_cids = _read_pubchem_records(check, records_path)
    requested_set = _check_pubchem_manifest_fields(check, manifest, records_path, record_cids)
    metadata_paths, batch_requested, batch_returned = _check_pubchem_batches(check, hashes, cache_dir)
    if not metadata_paths:
        check.fail("PubChem cache has no batch metadata files")
    if not record_cids <= batch_returned:
        check.fail(f"{len(record_cids - batch_returned)} normalized CIDs have no validated raw batch")
    mapping_cids = _mapping_pubchem_cids(mapping_rows, args.pubchem_cid_column)
    uncovered = mapping_cids - record_cids
    if uncovered:
        check.fail(f"offline PubChem cache does not cover {len(uncovered)} mapping CIDs")
    check.metrics = {
        "manifest_requested_cids": len(requested_set),
        "normalized_record_cids": len(record_cids),
        "batch_metadata_files": len(metadata_paths),
        "batch_requested_union": len(batch_requested),
        "batch_returned_union": len(batch_returned),
        "mapping_pubchem_cids": len(mapping_cids),
        "mapping_pubchem_cids_uncovered": len(uncovered),
        "mapping_pubchem_offline_coverage_percent": _percent(len(mapping_cids) - len(uncovered), len(mapping_cids)),
        "failure_count": len(check.failures),
    }
    return check


def _read_pubchem_records(check: Check, path: Path) -> set[int]:
    cids: set[int] = set()
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if line.strip():
                    _add_pubchem_record(check, cids, json.loads(line), line_number)
    except (OSError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot parse PubChem records {path}: {exc}") from exc
    return cids


def _add_pubchem_record(check: Check, cids: set[int], record: Any, line_number: int) -> None:
    cid = record.get("CID") if isinstance(record, dict) else None
    if not isinstance(cid, int) or cid <= 0 or cid in cids:
        check.fail(f"records.jsonl line {line_number}: invalid or duplicate CID {cid!r}")
        return
    has_structure = _clean(record.get("InChIKey")) and any(
        _clean(record.get(key)) for key in ("IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES")
    )
    if not has_structure:
        check.fail(f"records.jsonl line {line_number}: CID {cid} lacks validated structure fields")
        return
    cids.add(cid)


def _check_pubchem_manifest_fields(
    check: Check, manifest: Mapping[str, Any], records_path: Path, record_cids: set[int]
) -> set[int]:
    requested = manifest.get("requested_cids")
    requested_values = requested if isinstance(requested, list) else []
    valid_requested = isinstance(requested, list) and all(isinstance(cid, int) for cid in requested_values)
    requested_set = set(requested_values) if valid_requested else set()
    if not valid_requested:
        check.fail("PubChem manifest requested_cids is invalid")
    elif len(requested_set) != len(requested_values):
        check.fail("PubChem manifest requested_cids contains duplicates")
    _check_pubchem_manifest_counts(check, manifest, records_path, requested_set, record_cids)
    return requested_set


def _check_pubchem_manifest_counts(
    check: Check,
    manifest: Mapping[str, Any],
    records_path: Path,
    requested: set[int],
    records: set[int],
) -> None:
    checks = (
        (manifest.get("missing_cids") in ([], None), "PubChem manifest reports missing CIDs"),
        (manifest.get("invalid_structures") in ({}, None), "PubChem manifest reports invalid structures"),
        (requested == records, f"PubChem normalized coverage differs: missing={len(requested-records)}, extra={len(records-requested)}"),
        (manifest.get("requested_count") == len(requested), "PubChem requested_count does not match requested_cids"),
        (manifest.get("returned_count") == len(records), "PubChem returned_count does not match records.jsonl"),
        (manifest.get("records_bytes") == records_path.stat().st_size, "PubChem records_bytes does not match records.jsonl"),
    )
    for passed, message in checks:
        if not passed:
            check.fail(message)


def _check_pubchem_batches(
    check: Check, hashes: HashCache, cache_dir: Path
) -> tuple[list[Path], set[int], set[int]]:
    batch_dir = cache_dir / "batches"
    paths = sorted(batch_dir.glob("*.metadata.json")) if batch_dir.is_dir() else []
    requested: set[int] = set()
    returned: set[int] = set()
    for path in paths:
        result = _check_pubchem_batch(check, hashes, batch_dir, path)
        if result:
            requested.update(result[0])
            returned.update(result[1])
    return paths, requested, returned


def _check_pubchem_batch(
    check: Check, hashes: HashCache, batch_dir: Path, metadata_path: Path
) -> tuple[list[int], set[int]] | None:
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        raw_path = batch_dir / _clean(metadata.get("content_file"))
        _validate_file_hash(check, hashes, str(raw_path), _clean(metadata.get("content_sha256")), metadata_path, f"PubChem batch {metadata_path.name}")
        raw_cids = _parse_pubchem_raw(raw_path)
    except (OSError, json.JSONDecodeError, VerificationError) as exc:
        check.fail(f"cannot validate PubChem batch {metadata_path}: {exc}")
        return None
    requested = metadata.get("requested_cids")
    returned = metadata.get("returned_cids")
    if not isinstance(requested, list) or not all(isinstance(cid, int) for cid in requested):
        check.fail(f"{metadata_path.name}: invalid requested_cids")
        return None
    if not isinstance(returned, list) or set(returned) != raw_cids:
        check.fail(f"{metadata_path.name}: returned_cids differs from raw content")
        return None
    if metadata.get("missing_cids") != sorted(set(requested) - raw_cids):
        check.fail(f"{metadata_path.name}: missing_cids is inconsistent")
    if metadata.get("content_bytes") != raw_path.stat().st_size:
        check.fail(f"{metadata_path.name}: content_bytes is inconsistent")
    return requested, raw_cids


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True, help="fda_legacy_mapping_v3.csv")
    parser.add_argument("--terminal-dispositions", type=Path, required=True, help="fda_legacy_terminal_dispositions_v3.csv")
    parser.add_argument("--named-manifest", type=Path, required=True)
    parser.add_argument("--pubchem-manifest", type=Path, required=True)
    parser.add_argument(
        "--repair-quarantine",
        type=Path,
        help="fda_repair_quarantine.csv; required when unsafe parent-collapse rows exist",
    )
    parser.add_argument("--fda-source-sdf", type=Path, help="Exact 1,876-record FDA source SDF")
    parser.add_argument("--output", type=Path, required=True, help="JSON verification report")
    parser.add_argument("--row-column", default="mapping_row_number")
    parser.add_argument("--rdk-column", default="rdk_id")
    parser.add_argument("--resolution-column", default="identity_resolution_status")
    parser.add_argument("--disposition-column", default="terminal_disposition")
    parser.add_argument("--regulatory-column", default="regulatory_status")
    parser.add_argument("--preferred-name-column", default="resolved_preferred_name")
    parser.add_argument("--identity-evidence-column", default="identity_evidence_codes")
    parser.add_argument("--approval-evidence-column", default="approval_evidence_codes")
    parser.add_argument("--path-column", default="selected_pdbqt_path")
    parser.add_argument("--sha-column", default="selected_pdbqt_sha256")
    parser.add_argument("--file-status-column", default="legacy_file_status")
    parser.add_argument("--pubchem-cid-column", default="pubchem_cid_resolved")
    parser.add_argument("--terminal-value", action="append", help="Allowed terminal value; repeat to replace defaults")
    parser.add_argument("--fda-regulatory-status", action="append", help="FDA status; repeat to replace defaults")
    parser.add_argument("--structure-evidence-code", action="append", help="Strong structure code; repeat to replace defaults")
    parser.add_argument("--fda-source-code", action="append", help="FDA source code; repeat to replace defaults")
    parser.add_argument("--usable-status", action="append", help="Usable file status; repeat to replace defaults")
    parser.add_argument(
        "--usable-identity-evidence-code",
        action="append",
        help="Recognized identity evidence for any usable row; repeat to replace defaults",
    )
    parser.add_argument("--expected-prepared-usable", type=int, default=4729)
    parser.add_argument("--expected-prepared-unusable", type=int, default=1)
    parser.add_argument("--canonical-parent-id-column", default="canonical_parent_id")
    parser.add_argument("--canonical-status-column", default="canonical_materialization_status")
    parser.add_argument("--canonical-path-column", default="canonical_pdbqt_path")
    parser.add_argument("--canonical-sha-column", default="canonical_pdbqt_sha256")
    parser.add_argument("--canonical-status", action="append", help="Canonical status; repeat to replace ready/quarantine defaults")
    parser.add_argument("--canonical-ready-status", action="append", help="Canonical status requiring a ready PDBQT")
    parser.add_argument("--identity-parent-column", default="identity_parent_inchikey")
    parser.add_argument("--approved-source-index-column", default="approved_source_record_index")
    parser.add_argument("--terminal-reason-column", default="terminal_reason_codes")
    parser.add_argument("--unsafe-canonical-status", default=UNSAFE_CANONICAL_STATUS)
    parser.add_argument("--unsafe-reason-code", default=UNSAFE_REASON_CODE)
    parser.add_argument("--unsafe-mapping-reason-code", default=UNSAFE_MAPPING_REASON_CODE)
    parser.add_argument("--unsafe-source-marker", default=UNSAFE_SOURCE_MARKER)
    parser.add_argument("--repair-scope-column", default="scope")
    parser.add_argument("--repair-source-index-column", default="source_record_indices")
    parser.add_argument("--repair-parent-column", default="parent_inchikey")
    parser.add_argument("--repair-reason-column", default="reason_codes")
    parser.add_argument("--repair-source-scope", default="approved_source_record")
    parser.add_argument("--connectivity-state-column", default="prepared_connectivity_state")
    parser.add_argument("--score-reuse-column", default="legacy_score_reuse_status")
    parser.add_argument("--expected-connectivity-exact", type=int, default=3760)
    parser.add_argument("--expected-connectivity-parent", type=int, default=963)
    parser.add_argument("--expected-connectivity-ambiguous", type=int, default=6)
    parser.add_argument("--expected-connectivity-unverifiable", type=int, default=1)
    parser.add_argument("--expected-reuse-coordinate-exact", type=int, default=3760)
    parser.add_argument("--expected-reuse-coordinate-parent", type=int, default=962)
    parser.add_argument("--expected-no-reuse-unsafe-parent", type=int, default=1)
    parser.add_argument("--ambiguous-fda-rdk-id", default="rdk_0003623")
    parser.add_argument("--named-status-column", default="materialization_status")
    parser.add_argument("--named-path-column", default="output_pdbqt_path")
    parser.add_argument("--named-sha-column", default="output_pdbqt_sha256")
    parser.add_argument("--sidecar-path-column", default="sidecar_path")
    parser.add_argument("--sidecar-sha-column", default="sidecar_sha256")
    parser.add_argument("--named-ready-status", action="append", help="Named ready status; repeat to replace default")
    parser.add_argument("--named-quarantine-status", action="append", help="Named quarantine status; repeat to replace default")
    parser.add_argument("--expected-named-ready", type=int, default=1850)
    parser.add_argument("--expected-named-quarantine", type=int, default=13)
    parser.add_argument("--named-source-indices-column", default="source_record_indices")
    parser.add_argument("--expected-fda-source-records", type=int, default=1876)
    parser.add_argument("--expected-named-source-records", type=int, default=1863)
    parser.add_argument("--expected-repair-source-records", type=int, default=13)
    parser.add_argument("--regression-mapping-row", default="4714")
    parser.add_argument("--regression-rdk-id", default="rdk_0004924")
    parser.add_argument("--regression-preferred-name", default="deserpidine")
    parser.add_argument(
        "--regression-selected-sha256",
        default="0cf969069aa5bcc708b2a9572a8a9cbdb07d6a28d059cbfb63906cce66367005",
    )
    return parser


def _apply_cli_defaults(args: argparse.Namespace) -> None:
    defaults = {
        "terminal_value": list(TERMINAL_DISPOSITIONS),
        "fda_regulatory_status": list(FDA_REGULATORY_STATUSES),
        "structure_evidence_code": list(STRUCTURE_EVIDENCE_CODES),
        "fda_source_code": list(FDA_SOURCE_CODES),
        "usable_status": list(DEFAULT_USABLE_STATUSES),
        "usable_identity_evidence_code": list(USABLE_IDENTITY_EVIDENCE_CODES),
        "canonical_status": ["ready", "quarantined", "quarantined_non_dockable", args.unsafe_canonical_status],
        "canonical_ready_status": ["ready"],
        "named_ready_status": ["ready"],
        "named_quarantine_status": ["quarantine", "quarantined", "quarantined_non_dockable"],
    }
    for name, value in defaults.items():
        setattr(args, name, getattr(args, name) or value)


def verify(args: argparse.Namespace) -> dict[str, Any]:
    mapping_path = args.mapping.expanduser().resolve()
    terminal_path = args.terminal_dispositions.expanduser().resolve()
    named_path = args.named_manifest.expanduser().resolve()
    pubchem_path = args.pubchem_manifest.expanduser().resolve()
    repair_path = args.repair_quarantine.expanduser().resolve() if args.repair_quarantine else None
    source_sdf = args.fda_source_sdf.expanduser().resolve() if args.fda_source_sdf else None
    mapping_fields, mapping_rows = _read_csv(mapping_path)
    terminal_fields, terminal_rows = _read_csv(terminal_path)
    named_fields, named_rows = _read_csv(named_path)
    if repair_path is not None:
        repair_fields, repair_rows = _read_csv(repair_path)
    else:
        repair_fields, repair_rows = [], []
    mapping_required = (
        args.row_column,
        args.rdk_column,
        args.resolution_column,
        args.disposition_column,
        args.regulatory_column,
        args.preferred_name_column,
        args.identity_evidence_column,
        args.approval_evidence_column,
        args.path_column,
        args.sha_column,
        args.file_status_column,
        args.pubchem_cid_column,
        args.canonical_parent_id_column,
        args.canonical_status_column,
        args.canonical_path_column,
        args.canonical_sha_column,
        args.identity_parent_column,
        args.approved_source_index_column,
        args.terminal_reason_column,
        args.connectivity_state_column,
        args.score_reuse_column,
    )
    terminal_required = (
        args.row_column,
        args.rdk_column,
        args.resolution_column,
        args.disposition_column,
    )
    named_required = (
        args.named_status_column,
        args.named_path_column,
        args.named_sha_column,
        args.sidecar_path_column,
        args.sidecar_sha_column,
    )
    _require_columns(mapping_path, mapping_fields, mapping_required)
    _require_columns(terminal_path, terminal_fields, terminal_required)
    _require_columns(named_path, named_fields, named_required)
    if source_sdf is not None and args.named_source_indices_column not in named_fields:
        raise VerificationError(
            f"{named_path} is missing required column: {args.named_source_indices_column}"
        )

    hashes = HashCache()
    checks: list[Check] = []
    checks.append(check_mapping_uniqueness(mapping_rows, args.row_column, args.rdk_column))
    checks.append(check_prepared_files(mapping_rows, mapping_path, args, hashes))
    checks.append(check_prepared_connectivity_partition(mapping_rows, args, hashes))
    checks.append(check_regression_sentinel(mapping_rows, args))
    partition_check, terminal_by_key = check_terminal_partition(mapping_rows, terminal_rows, args)
    checks.append(partition_check)
    checks.append(check_no_nonterminal_usable(mapping_rows, terminal_by_key, args))
    checks.append(check_fda_evidence(mapping_rows, args))
    checks.append(check_fda_canonical_parents(mapping_rows, args))
    checks.append(
        check_unsafe_parent_quarantine(
            mapping_rows,
            repair_path,
            repair_fields,
            repair_rows,
            args,
        )
    )
    checks.append(check_named_manifest(named_rows, named_path, args, hashes))
    checks.append(check_fda_source_partition(named_rows, repair_rows, source_sdf, args))
    checks.append(check_pubchem_cache(pubchem_path, mapping_rows, args, hashes))
    passed = all(check.passed for check in checks)
    return {
        "schema": "atlas.fda-terminal-bundle-verification.v1",
        "created_at_utc": _utc_now(),
        "passed": passed,
        "artifacts": {
            "mapping": {"path": str(mapping_path), "sha256": hashes.sha256(mapping_path)},
            "terminal_dispositions": {"path": str(terminal_path), "sha256": hashes.sha256(terminal_path)},
            "named_manifest": {"path": str(named_path), "sha256": hashes.sha256(named_path)},
            "pubchem_manifest": {"path": str(pubchem_path), "sha256": hashes.sha256(pubchem_path)},
            **(
                {
                    "repair_quarantine": {
                        "path": str(repair_path),
                        "sha256": hashes.sha256(repair_path),
                    }
                }
                if repair_path is not None
                else {}
            ),
            **(
                {
                    "fda_source_sdf": {
                        "path": str(source_sdf),
                        "sha256": hashes.sha256(source_sdf),
                    }
                }
                if source_sdf is not None
                else {}
            ),
        },
        "checks": [check.as_dict() for check in checks],
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _apply_cli_defaults(args)
    try:
        report = verify(args)
    except Exception as exc:
        report = {
            "schema": "atlas.fda-terminal-bundle-verification.v1",
            "created_at_utc": _utc_now(),
            "passed": False,
            "fatal_error": f"{type(exc).__name__}: {exc}",
            "checks": [],
        }
    try:
        _atomic_json(args.output.expanduser().resolve(), report)
    except OSError as exc:
        print(f"failed to write verification report: {exc}", file=sys.stderr)
        return 2
    print(f"FDA terminal bundle verification: passed={report['passed']} report={args.output}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
