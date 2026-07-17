"""Hash-bound linkage between an FDA reconciliation bundle and a launched run."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml  # type: ignore[import-untyped]


BINDING_SCHEMA_VERSION = 1
BINDING_FILENAME = "fda_legacy_run_binding_v2.json"
_RECONCILIATION_VERSION = "fda-legacy-reconciliation-v2"
_REDOCK_LINKAGE_FILENAME = "fda_legacy_redock_linkage_v2.csv"
_RECONCILIATION_SUMMARY_FILENAME = "fda_legacy_reconciliation_v2_summary.json"
_READY_REDOCK_STATUS = "parent_prepared_ready_no_docking_launched"
_READY_DELTA_STATUS = "ready"
_PDB_NAME = re.compile(r"^[0-9A-Za-z]{4}\.pdb$")


@dataclass(frozen=True)
class FDALegacyRunBindingOutput:
    binding_json: Path


@dataclass(frozen=True)
class _Snapshot:
    path: Path
    data: bytes
    sha256: str

    @property
    def size(self) -> int:
        return len(self.data)


def bind_legacy_reconciliation_run(
    *,
    reconciliation_dir: Path,
    run_id: str,
    target_dir: Path,
    target_manifest_csv: Path,
    delta_summary_json: Path,
    run_manifest_yaml: Path,
    cpu_count: int,
    observed_worker_concurrency: int | None = None,
    expected_target_count: int = 93,
    expected_ready_count: int = 586,
    expected_apo_holo_mode: str = "holo",
    expected_run_status: str = "running",
) -> tuple[FDALegacyRunBindingOutput, dict[str, Any]]:
    """Create one post-launch binding artifact without changing pre-launch files."""

    reconciliation_dir = _require_dir(reconciliation_dir, "reconciliation directory")
    target_dir = _require_dir(target_dir, "frozen target directory")
    target_manifest = _snapshot(target_manifest_csv, "frozen target manifest")
    delta_summary = _snapshot(delta_summary_json, "delta summary")
    run_manifest = _snapshot(run_manifest_yaml, "run manifest")
    redock_linkage = _snapshot(
        reconciliation_dir / _REDOCK_LINKAGE_FILENAME,
        "pre-launch redock linkage",
    )
    reconciliation_summary = _snapshot(
        reconciliation_dir / _RECONCILIATION_SUMMARY_FILENAME,
        "reconciliation summary",
    )

    if not run_id.strip():
        raise ValueError("run_id must not be empty")
    if not 1 <= cpu_count <= 32:
        raise ValueError("cpu_count must be between 1 and 32")
    if observed_worker_concurrency is not None and not (
        1 <= observed_worker_concurrency <= 32
    ):
        raise ValueError("observed_worker_concurrency must be between 1 and 32")
    if expected_target_count <= 0 or expected_ready_count <= 0:
        raise ValueError("expected counts must be positive")

    reconciliation_payload = _json_mapping(reconciliation_summary)
    if reconciliation_payload.get("reconciliation_version") != _RECONCILIATION_VERSION:
        raise ValueError(
            "reconciliation summary is not fda-legacy-reconciliation-v2"
        )

    delta_payload = _json_mapping(delta_summary)
    delta_manifest = _delta_manifest_snapshot(delta_payload, delta_summary.path)
    redock_rows = _csv_rows(redock_linkage)
    delta_rows = _csv_rows(delta_manifest)
    ready_binding = _ready_ligand_binding(
        redock_rows=redock_rows,
        delta_rows=delta_rows,
        delta_payload=delta_payload,
        expected_ready_count=expected_ready_count,
    )
    target_binding = _target_binding(
        target_dir, target_manifest, expected_target_count
    )

    run_payload = _yaml_mapping(run_manifest)
    run_binding = _run_binding(
        run_payload=run_payload,
        run_id=run_id,
        target_dir=target_dir,
        target_binding=target_binding,
        delta_payload=delta_payload,
        expected_apo_holo_mode=expected_apo_holo_mode,
        expected_run_status=expected_run_status,
    )

    output = FDALegacyRunBindingOutput(
        binding_json=reconciliation_dir / BINDING_FILENAME
    )
    protected_inputs = {
        redock_linkage.path,
        reconciliation_summary.path,
        target_manifest.path,
        delta_summary.path,
        delta_manifest.path,
        run_manifest.path,
    }
    if output.binding_json.resolve() in protected_inputs:
        raise ValueError("run binding output collides with an input")

    observed_at = datetime.now(timezone.utc).isoformat()
    payload: dict[str, Any] = {
        "schema_version": BINDING_SCHEMA_VERSION,
        "binding_type": "fda-legacy-reconciliation-post-launch-v2",
        "observed_at": observed_at,
        "run": run_binding,
        "launch_resources": {
            "cpu_count": cpu_count,
            "cpu_count_evidence": "operator_launch_declaration",
            "cpu_count_manifest_limitation": (
                "The live run manifest records host cores, not the launch CPU cap."
            ),
            "observed_worker_concurrency": observed_worker_concurrency,
            "observed_worker_concurrency_evidence": (
                "operator_supplied_process_table_observation"
                if observed_worker_concurrency is not None
                else "not_observed"
            ),
            "apo_holo_mode": expected_apo_holo_mode,
        },
        "frozen_target_manifest": target_binding,
        "ready_ligand_binding": ready_binding,
        "source_artifacts": {
            "reconciliation_summary": _metadata(reconciliation_summary),
            "prelaunch_redock_linkage": {
                **_metadata(redock_linkage),
                "semantics": "pre_launch_plan",
                "mutated_by_binding": False,
            },
            "frozen_target_selection_manifest": _metadata(target_manifest),
            "delta_summary": _metadata(delta_summary),
            "delta_manifest": _metadata(delta_manifest),
            "run_manifest_snapshot": {
                **_metadata(run_manifest),
                "hash_scope": "exact_live_manifest_bytes_at_observed_at",
                "mutable_source_expected_to_change_while_run_active": True,
            },
        },
        "validation": {
            "run_exists": True,
            "run_status_matches_expected": True,
            "run_target_set_matches_selection_and_frozen_files": True,
            "run_library_matches_delta_summary": True,
            "run_apo_holo_mode_matches_expected": True,
            "ready_ligand_sets_and_hashes_match": True,
            "prelaunch_linkage_semantics_preserved": True,
        },
        "policy": {
            "binding_artifact_only": True,
            "run_outputs_mutated": False,
            "prelaunch_redock_linkage_rewritten": False,
            "binding_does_not_claim_run_completion_or_score_reuse": True,
        },
    }
    _write_json(output.binding_json, payload)
    return output, payload


def _ready_ligand_binding(
    *,
    redock_rows: Sequence[Mapping[str, str]],
    delta_rows: Sequence[Mapping[str, str]],
    delta_payload: Mapping[str, Any],
    expected_ready_count: int,
) -> dict[str, Any]:
    ready_redock = _rows_with_status(
        redock_rows,
        field="redock_linkage_status",
        status=_READY_REDOCK_STATUS,
        expected_count=expected_ready_count,
        label="ready redock rows",
    )
    ready_delta = _rows_with_status(
        delta_rows,
        field="materialization_status",
        status=_READY_DELTA_STATUS,
        expected_count=expected_ready_count,
        label="ready delta rows",
    )
    expected_inventory = _validated_delta_inventory(
        delta_payload, expected_ready_count
    )
    redock_hashes = _unique_hash_map(
        ready_redock,
        key_field="desired_parent_inchikey",
        hash_field="named_parent_pdbqt_sha256",
        label="pre-launch ready rows",
    )
    delta_hashes = _unique_hash_map(
        ready_delta,
        key_field="parent_inchikey",
        hash_field="output_pdbqt_sha256",
        label="delta ready rows",
    )
    _validate_matching_ready_hashes(redock_hashes, delta_hashes)
    verified_files, verified_bytes = _verify_ready_pdbqts(ready_delta)
    ready_digest = _row_set_digest(redock_hashes)
    return {
        "ready_row_count": expected_ready_count,
        "ready_parent_hash_set_sha256": ready_digest,
        "ready_parent_hash_set_canonicalization": (
            "sorted uppercase_parent_inchikey<TAB>lowercase_pdbqt_sha256<LF>"
        ),
        "delta_library_inventory_sha256": str(
            expected_inventory.get("sha256", "")
        ),
        "verified_current_pdbqt_count": verified_files,
        "verified_current_pdbqt_bytes": verified_bytes,
        "prelaunch_status": _READY_REDOCK_STATUS,
        "delta_status": _READY_DELTA_STATUS,
    }


def _rows_with_status(
    rows: Sequence[Mapping[str, str]],
    *,
    field: str,
    status: str,
    expected_count: int,
    label: str,
) -> list[Mapping[str, str]]:
    selected = [row for row in rows if row.get(field) == status]
    if len(selected) != expected_count:
        raise ValueError(
            f"{label}: expected {expected_count}, observed {len(selected)}"
        )
    return selected


def _validated_delta_inventory(
    delta_payload: Mapping[str, Any], expected_count: int
) -> Mapping[str, Any]:
    if _as_int(delta_payload.get("ready_pdbqt_count")) != expected_count:
        raise ValueError("delta summary ready_pdbqt_count does not match expectation")
    inventory = _mapping(delta_payload.get("library_inventory"))
    if _as_int(inventory.get("count")) != expected_count:
        raise ValueError("delta library inventory count does not match expectation")
    return inventory


def _validate_matching_ready_hashes(
    redock_hashes: Mapping[str, str], delta_hashes: Mapping[str, str]
) -> None:
    if redock_hashes == delta_hashes:
        return
    missing_from_delta = sorted(set(redock_hashes) - set(delta_hashes))
    missing_from_redock = sorted(set(delta_hashes) - set(redock_hashes))
    hash_mismatch = sorted(
        key
        for key in set(redock_hashes) & set(delta_hashes)
        if redock_hashes[key] != delta_hashes[key]
    )
    raise ValueError(
        "ready ligand rows differ between pre-launch linkage and delta manifest: "
        f"missing_from_delta={missing_from_delta[:5]} "
        f"missing_from_redock={missing_from_redock[:5]} "
        f"hash_mismatch={hash_mismatch[:5]}"
    )


def _verify_ready_pdbqts(
    ready_rows: Sequence[Mapping[str, str]],
) -> tuple[int, int]:
    verified_bytes = 0
    for row in ready_rows:
        path = _require_file(
            Path(row.get("output_pdbqt_path", "")), "delta ready PDBQT"
        )
        observed = _sha256_bytes(path.read_bytes())
        expected = row.get("output_pdbqt_sha256", "").lower()
        if observed != expected:
            raise ValueError(f"delta ready PDBQT checksum mismatch: {path}")
        verified_bytes += path.stat().st_size
    return len(ready_rows), verified_bytes


def _target_binding(
    target_dir: Path, target_manifest: _Snapshot, expected_count: int
) -> dict[str, Any]:
    selection_ids = _target_selection_ids(target_manifest, expected_count)
    target_files = _target_paths(target_dir, expected_count)
    rows = _target_file_rows(target_files)
    file_ids = [str(row["pdb_id"]) for row in rows]
    if sorted(selection_ids) != sorted(file_ids):
        raise ValueError(
            "target selection manifest IDs differ from frozen target PDB files"
        )
    digest_source = "".join(
        f"{row['pdb_id']}\t{row['sha256']}\t{row['bytes']}\n" for row in rows
    ).encode("utf-8")
    return {
        "source_directory": str(target_dir),
        "selection_manifest_path": str(target_manifest.path),
        "selection_manifest_sha256": target_manifest.sha256,
        "target_count": len(rows),
        "target_set_sha256": _sha256_bytes(digest_source),
        "target_set_canonicalization": (
            "sorted uppercase_pdb_id<TAB>lowercase_file_sha256<TAB>bytes<LF>"
        ),
        "targets": rows,
    }


def _target_selection_ids(
    target_manifest: _Snapshot, expected_count: int
) -> list[str]:
    selection_rows = _csv_rows(target_manifest)
    if len(selection_rows) != expected_count:
        raise ValueError(
            f"target selection rows: expected {expected_count}, "
            f"observed {len(selection_rows)}"
        )
    selection_ids = [row.get("pdb_id", "").upper() for row in selection_rows]
    if any(not value for value in selection_ids):
        raise ValueError("target selection manifest contains a missing pdb_id")
    if len(set(selection_ids)) != expected_count:
        raise ValueError("target selection manifest contains duplicate pdb_id values")
    return selection_ids


def _target_paths(target_dir: Path, expected_count: int) -> list[Path]:
    paths = sorted(
        (
            path
            for path in target_dir.iterdir()
            if path.is_file() and _PDB_NAME.fullmatch(path.name)
        ),
        key=lambda path: path.name.upper(),
    )
    if len(paths) != expected_count:
        raise ValueError(
            f"frozen target PDB count: expected {expected_count}, "
            f"observed {len(paths)}"
        )
    return paths


def _target_file_rows(target_files: Sequence[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for path in target_files:
        pdb_id = path.stem.upper()
        if pdb_id in seen:
            raise ValueError(f"duplicate frozen target ID: {pdb_id}")
        seen.add(pdb_id)
        snapshot = _snapshot(path, f"frozen target {pdb_id}")
        rows.append(
            {
                "pdb_id": pdb_id,
                "path": str(snapshot.path),
                "sha256": snapshot.sha256,
                "bytes": snapshot.size,
            }
        )
    return rows


def _run_binding(
    *,
    run_payload: Mapping[str, Any],
    run_id: str,
    target_dir: Path,
    target_binding: Mapping[str, Any],
    delta_payload: Mapping[str, Any],
    expected_apo_holo_mode: str,
    expected_run_status: str,
) -> dict[str, Any]:
    observed_status = _validate_run_header(
        run_payload, run_id, expected_run_status
    )
    command = _mapping(run_payload.get("command"))
    paths = _mapping(run_payload.get("paths"))
    observed_mode, observed_library = _validate_run_scope(
        command=command,
        paths=paths,
        target_dir=target_dir,
        delta_payload=delta_payload,
        expected_apo_holo_mode=expected_apo_holo_mode,
    )
    summary = _mapping(run_payload.get("summary"))
    run_target_ids = _run_target_ids(summary)
    frozen_target_ids = _frozen_target_ids(target_binding)
    if run_target_ids != frozen_target_ids:
        raise ValueError("run target universe differs from frozen target manifest")

    resources = _mapping(run_payload.get("resources"))
    timing = _mapping(run_payload.get("timing"))
    git = _mapping(run_payload.get("git"))
    run_dir = Path(str(paths.get("run_dir", ""))).expanduser().resolve()
    config_file = str(paths.get("config_file", "")).strip()
    config_snapshot_path = run_dir / config_file if config_file else run_dir
    return {
        "run_id": run_id,
        "observed_status": observed_status,
        "expected_status": expected_run_status,
        "created_at": timing.get("created_at"),
        "started_at": timing.get("started_at"),
        "finished_at": timing.get("finished_at"),
        "config_hash": command.get("config_hash"),
        "apo_holo_mode": observed_mode,
        "library_token": observed_library,
        "frozen_target_count": len(run_target_ids),
        "scheduler_admitted_target_count_at_observation": _as_int(
            summary.get("total_proteins_scheduled")
        ),
        "host": resources.get("host"),
        "host_n_cores": resources.get("n_cores"),
        "git_commit": git.get("commit"),
        "git_dirty": git.get("dirty"),
        "expected_config_snapshot": {
            "path": str(config_snapshot_path),
            "exists": config_snapshot_path.is_file(),
            "note": (
                "The manifest config_hash is retained even when the advertised "
                "configuration snapshot was not written."
            ),
        },
    }


def _validate_run_header(
    run_payload: Mapping[str, Any], run_id: str, expected_status: str
) -> str:
    observed_run_id = str(run_payload.get("run_id", ""))
    if observed_run_id != run_id:
        raise ValueError(
            f"run manifest ID mismatch: expected {run_id}, observed {observed_run_id}"
        )
    observed_status = str(run_payload.get("status", ""))
    if observed_status != expected_status:
        raise ValueError(
            "run manifest status mismatch: "
            f"expected {expected_status}, observed {observed_status}"
        )
    return observed_status


def _validate_run_scope(
    *,
    command: Mapping[str, Any],
    paths: Mapping[str, Any],
    target_dir: Path,
    delta_payload: Mapping[str, Any],
    expected_apo_holo_mode: str,
) -> tuple[str, str]:
    observed_mode = str(command.get("APO_HOLO_MODE", "")).casefold()
    if observed_mode != expected_apo_holo_mode.casefold():
        raise ValueError(
            "APO_HOLO_MODE mismatch: "
            f"expected {expected_apo_holo_mode}, observed {observed_mode}"
        )
    delta_dir = _require_dir(
        Path(str(delta_payload.get("output_dir", ""))), "delta output directory"
    )
    observed_library = str(command.get("TEST_MODE_ENABLE", ""))
    if observed_library != delta_dir.name:
        raise ValueError(
            "run library mismatch: "
            f"manifest={observed_library!r}, delta={delta_dir.name!r}"
        )
    manifest_target_dir = _require_dir(
        Path(str(paths.get("input_pdb_dir", ""))),
        "run manifest input PDB directory",
    )
    if manifest_target_dir != target_dir:
        raise ValueError(
            "run target directory mismatch: "
            f"manifest={manifest_target_dir}, expected={target_dir}"
        )
    return observed_mode, observed_library


def _run_target_ids(summary: Mapping[str, Any]) -> list[str]:
    raw_target_ids = summary.get("total_protein_list")
    if not isinstance(raw_target_ids, list):
        raise ValueError("run manifest summary.total_protein_list is not a list")
    return sorted({str(value).upper() for value in raw_target_ids})


def _frozen_target_ids(target_binding: Mapping[str, Any]) -> list[str]:
    manifest_rows = target_binding.get("targets")
    if not isinstance(manifest_rows, list):
        raise ValueError("generated frozen target manifest rows are unavailable")
    return sorted(
        str(row.get("pdb_id", "")).upper()
        for row in manifest_rows
        if isinstance(row, Mapping)
    )


def _delta_manifest_snapshot(
    delta_payload: Mapping[str, Any], delta_summary_path: Path
) -> _Snapshot:
    raw_path = str(delta_payload.get("manifest_csv", "")).strip()
    if not raw_path:
        raise ValueError("delta summary does not name manifest_csv")
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = delta_summary_path.parent / path
    manifest = _snapshot(path, "delta manifest")
    stated = str(delta_payload.get("manifest_sha256", "")).lower()
    if manifest.sha256 != stated:
        raise ValueError("delta manifest checksum differs from delta summary")
    return manifest


def _unique_hash_map(
    rows: Sequence[Mapping[str, str]],
    *,
    key_field: str,
    hash_field: str,
    label: str,
) -> dict[str, str]:
    output: dict[str, str] = {}
    for row in rows:
        key = row.get(key_field, "").strip().upper()
        digest = row.get(hash_field, "").strip().lower()
        if not key or not _is_sha256(digest):
            raise ValueError(f"{label} contains a missing key or invalid SHA-256")
        if key in output:
            raise ValueError(f"{label} contains duplicate key {key}")
        output[key] = digest
    return output


def _row_set_digest(rows: Mapping[str, str]) -> str:
    data = "".join(f"{key}\t{rows[key]}\n" for key in sorted(rows)).encode(
        "utf-8"
    )
    return _sha256_bytes(data)


def _snapshot(path: Path, label: str) -> _Snapshot:
    resolved = _require_file(path, label)
    data = resolved.read_bytes()
    if not data:
        raise ValueError(f"{label} is empty: {resolved}")
    return _Snapshot(path=resolved, data=data, sha256=_sha256_bytes(data))


def _require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} is not a file: {resolved}")
    return resolved


def _require_dir(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise FileNotFoundError(f"{label} is not a directory: {resolved}")
    return resolved


def _json_mapping(snapshot: _Snapshot) -> Mapping[str, Any]:
    try:
        payload = json.loads(snapshot.data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid JSON: {snapshot.path}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"JSON root is not an object: {snapshot.path}")
    return payload


def _yaml_mapping(snapshot: _Snapshot) -> Mapping[str, Any]:
    try:
        payload = yaml.safe_load(snapshot.data.decode("utf-8")) or {}
    except (UnicodeDecodeError, yaml.YAMLError) as exc:
        raise ValueError(f"invalid YAML: {snapshot.path}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"YAML root is not an object: {snapshot.path}")
    return payload


def _csv_rows(snapshot: _Snapshot) -> list[dict[str, str]]:
    try:
        text = snapshot.data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"CSV is not UTF-8: {snapshot.path}") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if not reader.fieldnames:
        raise ValueError(f"CSV has no header: {snapshot.path}")
    return [
        {key: str(value or "").strip() for key, value in row.items()}
        for row in reader
    ]


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _as_int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _is_sha256(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9a-f]{64}", value))


def _metadata(snapshot: _Snapshot) -> dict[str, Any]:
    return {
        "path": str(snapshot.path),
        "sha256": snapshot.sha256,
        "bytes": snapshot.size,
    }


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(tmp, path)


__all__ = [
    "BINDING_FILENAME",
    "FDALegacyRunBindingOutput",
    "bind_legacy_reconciliation_run",
]
