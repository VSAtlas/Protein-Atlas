"""Versioned, non-destructive reconciliation of the legacy FDA ligand library."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = 2
RECONCILIATION_VERSION = "fda-legacy-reconciliation-v2"
_HISTORICAL_HASH_FIELDS = (
    "ligand_pdbqt_sha256",
    "pdbqt_sha256",
    "ligand_sha256",
    "input_ligand_sha256",
)
_SCORE_LIGAND_FIELDS = (
    "rdk_id",
    "ligand_base",
    "ligand",
    "ligand_file",
    "Ligand_ID",
)

_MAPPING_FIELDS = [
    "legacy_reconciliation_version",
    "mapping_row_number",
    "rdk_id",
    "legacy_display_name",
    "legacy_generic_name",
    "legacy_drugcentral_id",
    "legacy_drugcentral_generic_name",
    "identity_resolution_status",
    "identity_change_applied",
    "identity_change_reason_codes",
    "resolved_drugcentral_id",
    "resolved_preferred_name",
    "approved_source_record_index",
    "approved_full_form_exact_inchikey",
    "approved_parent_inchikey",
    "selected_pdbqt_path",
    "selected_pdbqt_sha256",
    "legacy_file_status",
    "chemistry_redock_required",
]

_CHANGE_FIELDS = [
    "mapping_row_number",
    "rdk_id",
    "old_display_name",
    "new_display_name",
    "old_generic_name",
    "new_generic_name",
    "old_drugcentral_id",
    "new_drugcentral_id",
    "old_drugcentral_generic_name",
    "new_drugcentral_generic_name",
    "approved_source_record_index",
    "approved_full_form_exact_inchikey",
    "approved_parent_inchikey",
    "selected_pdbqt_path",
    "selected_pdbqt_sha256",
    "change_reason_codes",
]

_INVENTORY_FIELDS = [
    "mapping_row_number",
    "rdk_id",
    "legacy_mapping_path",
    "selected_pdbqt_path",
    "selected_pdbqt_sha256",
    "verified_pdbqt_sha256",
    "pdbqt_bytes",
    "preservation_status",
    "preservation_action",
    "identity_resolution_status",
    "legacy_display_name",
    "resolved_preferred_name",
]

_UNRESOLVED_FIELDS = [
    "mapping_row_number",
    "rdk_id",
    "resolution_status",
    "second_pass_priority",
    "legacy_display_name",
    "legacy_vendor_id",
    "legacy_cas",
    "legacy_exact_inchikey",
    "legacy_parent_inchikey",
    "audit_category",
    "audit_category_basis",
    "selected_pdbqt_path",
    "selected_pdbqt_sha256",
    "legacy_file_status",
    "repair_reason_codes",
    "recommended_second_pass",
]

_SCORE_PREFIX_FIELDS = [
    "score_source_csv",
    "score_source_sha256",
    "score_source_row_number",
]
_SCORE_JOIN_FIELDS = [
    "score_row_stable_id",
    "raw_ligand_identifier",
    "normalized_ligand_identifier",
    "joined_rdk_id",
    "candidate_mapping_rdk_ids",
    "identity_join_status",
    "identity_resolution_status",
    "resolved_drugcentral_id",
    "resolved_preferred_name",
    "approved_parent_inchikey",
    "current_pdbqt_path",
    "current_pdbqt_sha256",
    "observed_current_pdbqt_sha256",
    "current_file_checksum_status",
    "historical_pdbqt_sha256",
    "historical_checksum_status",
    "score_reconciliation_status",
    "score_reconciliation_reason_codes",
    "chemistry_redock_required",
]

_REDOCK_JOIN_FIELDS = [
    "linked_legacy_rdk_ids",
    "redock_trigger_rdk_ids",
    "trigger_mapping_rows",
    "legacy_selected_paths",
    "legacy_selected_sha256s",
    "linked_score_row_ids",
    "named_parent_materialization_status",
    "named_parent_pdbqt_path",
    "named_parent_pdbqt_sha256",
    "named_parent_status_reason_codes",
    "redock_linkage_status",
    "docking_launched",
]


@dataclass(frozen=True)
class FDALegacyReconciliationOutputs:
    versioned_mapping_csv: Path
    identity_changes_csv: Path
    file_inventory_csv: Path
    score_identity_joins_csv: Path
    unresolved_second_pass_csv: Path
    redock_linkage_csv: Path
    summary_json: Path


@dataclass(frozen=True)
class _ScoreResolution:
    repaired: Mapping[str, str] | None
    joined_rdk_id: str
    candidate_rdk_ids: tuple[str, ...]
    join_status: str
    reasons: tuple[str, ...]


def legacy_reconciliation_paths(output_dir: Path) -> FDALegacyReconciliationOutputs:
    return FDALegacyReconciliationOutputs(
        versioned_mapping_csv=output_dir / "fda_legacy_mapping_v2.csv",
        identity_changes_csv=output_dir / "fda_legacy_identity_changes_v2.csv",
        file_inventory_csv=output_dir / "fda_legacy_file_inventory_v2.csv",
        score_identity_joins_csv=output_dir / "fda_legacy_score_identity_joins_v2.csv",
        unresolved_second_pass_csv=output_dir
        / "fda_unresolved_second_pass_queue_v2.csv",
        redock_linkage_csv=output_dir / "fda_legacy_redock_linkage_v2.csv",
        summary_json=output_dir / "fda_legacy_reconciliation_v2_summary.json",
    )


def reconcile_legacy_fda(
    *,
    original_mapping_csv: Path,
    repaired_mapping_csv: Path,
    redock_delta_csv: Path,
    score_csvs: Sequence[Path],
    output_dir: Path,
    named_library_manifest_csv: Path | None = None,
) -> tuple[FDALegacyReconciliationOutputs, dict[str, Any]]:
    """Write a versioned reconciliation bundle without changing legacy inputs."""

    inputs = _resolved_inputs(
        original_mapping_csv,
        repaired_mapping_csv,
        redock_delta_csv,
        score_csvs,
        named_library_manifest_csv,
    )
    outputs = legacy_reconciliation_paths(Path(output_dir).expanduser().resolve())
    _validate_output_collisions(outputs, inputs.values())
    original_fields, original_rows = _read_csv(inputs["original_mapping_csv"])
    _, repaired_rows = _read_csv(inputs["repaired_mapping_csv"])
    _, delta_rows = _read_csv(inputs["redock_delta_csv"])
    _validate_mapping_pair(original_rows, repaired_rows)

    repaired_by_id = {row["rdk_id"]: row for row in repaired_rows}
    inventory = _file_inventory(original_rows, repaired_rows)
    versioned, changes = _mapping_outputs(
        original_rows, repaired_rows, original_fields
    )
    inventory_by_id = {row["rdk_id"]: row for row in inventory}
    unresolved = _unresolved_rows(repaired_rows)
    score_fields, score_joins = _score_identity_joins(
        inputs, repaired_by_id, inventory_by_id
    )
    redock = _redock_linkage(
        delta_rows,
        repaired_rows,
        inputs.get("named_library_manifest_csv"),
        score_joins,
    )
    _validate_current_counts(
        original_rows=original_rows,
        changes=changes,
        inventory=inventory,
        unresolved=unresolved,
        redock=redock,
        score_joins=score_joins,
    )

    outputs.versioned_mapping_csv.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        outputs.versioned_mapping_csv,
        versioned,
        [*original_fields, *_without_duplicates(_MAPPING_FIELDS, original_fields)],
    )
    _write_csv(outputs.identity_changes_csv, changes, _CHANGE_FIELDS)
    _write_csv(outputs.file_inventory_csv, inventory, _INVENTORY_FIELDS)
    _write_csv(outputs.score_identity_joins_csv, score_joins, score_fields)
    _write_csv(outputs.unresolved_second_pass_csv, unresolved, _UNRESOLVED_FIELDS)
    redock_fields = _redock_fields(delta_rows)
    _write_csv(outputs.redock_linkage_csv, redock, redock_fields)

    summary = _summary_payload(
        inputs=inputs,
        outputs=outputs,
        original_rows=original_rows,
        changes=changes,
        inventory=inventory,
        score_joins=score_joins,
        unresolved=unresolved,
        redock=redock,
    )
    _write_json(outputs.summary_json, summary)
    return outputs, summary


def _resolved_inputs(
    original_mapping_csv: Path,
    repaired_mapping_csv: Path,
    redock_delta_csv: Path,
    score_csvs: Sequence[Path],
    named_library_manifest_csv: Path | None,
) -> dict[str, Path]:
    values = {
        "original_mapping_csv": Path(original_mapping_csv).expanduser().resolve(),
        "repaired_mapping_csv": Path(repaired_mapping_csv).expanduser().resolve(),
        "redock_delta_csv": Path(redock_delta_csv).expanduser().resolve(),
    }
    if named_library_manifest_csv is not None:
        values["named_library_manifest_csv"] = (
            Path(named_library_manifest_csv).expanduser().resolve()
        )
    for index, path in enumerate(score_csvs, start=1):
        values[f"score_csv_{index}"] = Path(path).expanduser().resolve()
    missing = [str(path) for path in values.values() if not path.is_file()]
    if missing:
        raise ValueError("missing FDA legacy reconciliation input(s): " + ", ".join(missing))
    return values


def _mapping_outputs(
    original_rows: Sequence[Mapping[str, str]],
    repaired_rows: Sequence[Mapping[str, str]],
    original_fields: Sequence[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    versioned = []
    changes = []
    for original, repaired in zip(original_rows, repaired_rows, strict=True):
        output, change = _versioned_mapping_row(original, repaired, original_fields)
        versioned.append(output)
        if change is not None:
            changes.append(change)
    return versioned, changes


def _versioned_mapping_row(
    original: Mapping[str, str],
    repaired: Mapping[str, str],
    original_fields: Sequence[str],
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    output = {field: original.get(field, "") for field in original_fields}
    changed = _identity_change_applies(repaired)
    legacy = _legacy_identity_fields(original)
    if changed:
        output["display_name"] = repaired["corrected_preferred_name"]
        output["generic_name"] = repaired["corrected_preferred_name"]
        output["drugcentral_id"] = repaired["corrected_drugcentral_id"]
        output["drugcentral_generic_name"] = repaired["corrected_preferred_name"]
    output.update(legacy)
    output.update(_reconciliation_fields(repaired, changed))
    change = _change_row(original, repaired) if changed else None
    return output, change


def _legacy_identity_fields(original: Mapping[str, str]) -> dict[str, str]:
    return {
        "legacy_display_name": original.get("display_name", ""),
        "legacy_generic_name": original.get("generic_name", ""),
        "legacy_drugcentral_id": original.get("drugcentral_id", ""),
        "legacy_drugcentral_generic_name": original.get(
            "drugcentral_generic_name", ""
        ),
    }


def _reconciliation_fields(
    repaired: Mapping[str, str], changed: bool
) -> dict[str, str]:
    return {
        "legacy_reconciliation_version": RECONCILIATION_VERSION,
        "mapping_row_number": repaired["mapping_row_number"],
        "rdk_id": repaired["rdk_id"],
        "identity_resolution_status": _identity_resolution_status(repaired),
        "identity_change_applied": _bool(changed),
        "identity_change_reason_codes": repaired.get("repair_reason_codes", ""),
        "resolved_drugcentral_id": repaired.get("corrected_drugcentral_id", ""),
        "resolved_preferred_name": repaired.get("corrected_preferred_name", ""),
        "approved_source_record_index": repaired.get(
            "approved_source_record_index", ""
        ),
        "approved_full_form_exact_inchikey": repaired.get(
            "approved_full_form_exact_inchikey", ""
        ),
        "approved_parent_inchikey": repaired.get("approved_parent_inchikey", ""),
        "selected_pdbqt_path": repaired.get("selected_pdbqt_path", ""),
        "selected_pdbqt_sha256": repaired.get("selected_pdbqt_sha256", ""),
        "legacy_file_status": repaired.get("selected_pdbqt_status", ""),
        "chemistry_redock_required": repaired.get(
            "chemistry_redock_required", ""
        ),
    }


def _change_row(
    original: Mapping[str, str], repaired: Mapping[str, str]
) -> dict[str, str]:
    new_name = repaired["corrected_preferred_name"]
    return {
        "mapping_row_number": repaired["mapping_row_number"],
        "rdk_id": repaired["rdk_id"],
        "old_display_name": original.get("display_name", ""),
        "new_display_name": new_name,
        "old_generic_name": original.get("generic_name", ""),
        "new_generic_name": new_name,
        "old_drugcentral_id": original.get("drugcentral_id", ""),
        "new_drugcentral_id": repaired["corrected_drugcentral_id"],
        "old_drugcentral_generic_name": original.get(
            "drugcentral_generic_name", ""
        ),
        "new_drugcentral_generic_name": new_name,
        "approved_source_record_index": repaired.get(
            "approved_source_record_index", ""
        ),
        "approved_full_form_exact_inchikey": repaired.get(
            "approved_full_form_exact_inchikey", ""
        ),
        "approved_parent_inchikey": repaired.get("approved_parent_inchikey", ""),
        "selected_pdbqt_path": repaired.get("selected_pdbqt_path", ""),
        "selected_pdbqt_sha256": repaired.get("selected_pdbqt_sha256", ""),
        "change_reason_codes": repaired.get("repair_reason_codes", ""),
    }


def _identity_change_applies(repaired: Mapping[str, str]) -> bool:
    return bool(
        repaired.get("repair_action") == "repair_name_and_identity"
        and repaired.get("corrected_drugcentral_id")
        and repaired.get("corrected_preferred_name")
    )


def _identity_resolution_status(repaired: Mapping[str, str]) -> str:
    if not repaired.get("corrected_drugcentral_id"):
        return "unresolved_requires_second_pass"
    if _identity_change_applies(repaired):
        return "manifest_resolved_corrected"
    return "manifest_resolved_unchanged"


def _file_inventory(
    original_rows: Sequence[Mapping[str, str]],
    repaired_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, Any]]:
    return [
        _file_inventory_row(original, repaired)
        for original, repaired in zip(original_rows, repaired_rows, strict=True)
    ]


def _file_inventory_row(
    original: Mapping[str, str], repaired: Mapping[str, str]
) -> dict[str, Any]:
    path_text = repaired.get("selected_pdbqt_path", "")
    expected_hash = repaired.get("selected_pdbqt_sha256", "")
    output = {
        "mapping_row_number": repaired["mapping_row_number"],
        "rdk_id": repaired["rdk_id"],
        "legacy_mapping_path": original.get("path", ""),
        "selected_pdbqt_path": path_text,
        "selected_pdbqt_sha256": expected_hash,
        "identity_resolution_status": _identity_resolution_status(repaired),
        "legacy_display_name": original.get("display_name", ""),
        "resolved_preferred_name": repaired.get("corrected_preferred_name", ""),
    }
    output.update(_verified_file_fields(path_text, expected_hash))
    return output


def _verified_file_fields(path_text: str, expected_hash: str) -> dict[str, Any]:
    path = Path(path_text).expanduser() if path_text else None
    if path is None or not path.is_file() or path.stat().st_size <= 0:
        return {
            "verified_pdbqt_sha256": "",
            "pdbqt_bytes": "",
            "preservation_status": "unusable_missing_or_empty",
            "preservation_action": "retained_unresolved_no_usable_file",
        }
    verified_hash = _sha256(path)
    if verified_hash != expected_hash:
        raise ValueError(f"legacy FDA PDBQT checksum changed: {path}")
    return {
        "verified_pdbqt_sha256": verified_hash,
        "pdbqt_bytes": path.stat().st_size,
        "preservation_status": "preserved_usable",
        "preservation_action": "retained_in_place_byte_verified",
    }


def _unresolved_rows(
    repaired_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    return [
        _unresolved_row(row)
        for row in repaired_rows
        if not row.get("corrected_drugcentral_id")
    ]


def _unresolved_row(row: Mapping[str, str]) -> dict[str, str]:
    category = row.get("audit_category", "")
    priority = "high" if category == "incorrect name–structure mapping" else "standard"
    return {
        "mapping_row_number": row["mapping_row_number"],
        "rdk_id": row["rdk_id"],
        "resolution_status": "unresolved_not_classified_non_fda",
        "second_pass_priority": priority,
        "legacy_display_name": row.get("original_display_name", ""),
        "legacy_vendor_id": row.get("original_mapping_vendor_id", ""),
        "legacy_cas": row.get("original_mapping_cas", ""),
        "legacy_exact_inchikey": row.get("original_mapping_exact_inchikey", ""),
        "legacy_parent_inchikey": row.get("original_mapping_parent_inchikey", ""),
        "audit_category": category,
        "audit_category_basis": row.get("audit_category_basis", ""),
        "selected_pdbqt_path": row.get("selected_pdbqt_path", ""),
        "selected_pdbqt_sha256": row.get("selected_pdbqt_sha256", ""),
        "legacy_file_status": row.get("selected_pdbqt_status", ""),
        "repair_reason_codes": row.get("repair_reason_codes", ""),
        "recommended_second_pass": (
            "stable_identifier_then_exact_or_full_parent_structure_review"
        ),
    }


def _score_identity_joins(
    inputs: Mapping[str, Path],
    repaired_by_id: Mapping[str, Mapping[str, str]],
    inventory_by_id: Mapping[str, Mapping[str, Any]],
) -> tuple[list[str], list[dict[str, Any]]]:
    score_paths = [
        path for key, path in inputs.items() if key.startswith("score_csv_")
    ]
    original_fields: list[str] = []
    output: list[dict[str, Any]] = []
    repaired_by_hash = _repaired_by_observed_hash(repaired_by_id, inventory_by_id)
    for path in score_paths:
        fields, rows = _read_csv(path)
        original_fields.extend(field for field in fields if field not in original_fields)
        source_hash = _sha256(path)
        output.extend(
            _joined_score_row(
                path,
                source_hash,
                index,
                row,
                repaired_by_id,
                repaired_by_hash,
                inventory_by_id,
            )
            for index, row in enumerate(rows, start=1)
        )
    fields = [
        *_SCORE_PREFIX_FIELDS,
        *_without_duplicates(original_fields, _SCORE_PREFIX_FIELDS),
        *_without_duplicates(
            _SCORE_JOIN_FIELDS, [*_SCORE_PREFIX_FIELDS, *original_fields]
        ),
    ]
    return fields, output


def _joined_score_row(
    path: Path,
    source_hash: str,
    row_number: int,
    score_row: Mapping[str, str],
    repaired_by_id: Mapping[str, Mapping[str, str]],
    repaired_by_hash: Mapping[str, Sequence[Mapping[str, str]]],
    inventory_by_id: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    raw_id = _first(score_row, _SCORE_LIGAND_FIELDS)
    rdk_id = _parse_score_rdk_id(raw_id)
    historical_hash = _first(score_row, _HISTORICAL_HASH_FIELDS)
    resolution = _resolve_score_identity(
        raw_id,
        rdk_id,
        historical_hash,
        repaired_by_id,
        repaired_by_hash,
        inventory_by_id,
    )
    output: dict[str, Any] = {
        "score_source_csv": str(path),
        "score_source_sha256": source_hash,
        "score_source_row_number": row_number,
        **score_row,
    }
    output.update(
        _score_join_fields(
            source_hash,
            row_number,
            raw_id,
            historical_hash,
            resolution,
            inventory_by_id,
        )
    )
    return output


def _score_join_fields(
    source_hash: str,
    row_number: int,
    raw_id: str,
    historical_hash: str,
    resolution: _ScoreResolution,
    inventory_by_id: Mapping[str, Mapping[str, Any]],
) -> dict[str, str]:
    repaired = resolution.repaired
    if repaired is None:
        return _missing_score_join_fields(
            source_hash, row_number, raw_id, historical_hash, resolution
        )
    inventory = inventory_by_id[resolution.joined_rdk_id]
    current_hash = repaired.get("selected_pdbqt_sha256", "")
    observed_hash = str(inventory.get("verified_pdbqt_sha256", ""))
    historical_status = _historical_checksum_status(observed_hash, historical_hash)
    score_status, score_reasons = _score_status(
        repaired, resolution.join_status, historical_status
    )
    return {
        "score_row_stable_id": f"{source_hash}:{row_number}",
        "raw_ligand_identifier": raw_id,
        "normalized_ligand_identifier": resolution.joined_rdk_id,
        "joined_rdk_id": resolution.joined_rdk_id,
        "candidate_mapping_rdk_ids": _joined(resolution.candidate_rdk_ids),
        "identity_join_status": resolution.join_status,
        "identity_resolution_status": _identity_resolution_status(repaired),
        "resolved_drugcentral_id": repaired.get("corrected_drugcentral_id", ""),
        "resolved_preferred_name": repaired.get("corrected_preferred_name", ""),
        "approved_parent_inchikey": repaired.get("approved_parent_inchikey", ""),
        "current_pdbqt_path": repaired.get("selected_pdbqt_path", ""),
        "current_pdbqt_sha256": current_hash,
        "observed_current_pdbqt_sha256": observed_hash,
        "current_file_checksum_status": _current_file_checksum_status(
            current_hash, observed_hash
        ),
        "historical_pdbqt_sha256": historical_hash,
        "historical_checksum_status": historical_status,
        "score_reconciliation_status": score_status,
        "score_reconciliation_reason_codes": ";".join(
            [*resolution.reasons, *score_reasons]
        ),
        "chemistry_redock_required": repaired.get(
            "chemistry_redock_required", ""
        ),
    }


def _missing_score_join_fields(
    source_hash: str,
    row_number: int,
    raw_id: str,
    historical_hash: str,
    resolution: _ScoreResolution,
) -> dict[str, str]:
    return {
        "score_row_stable_id": f"{source_hash}:{row_number}",
        "raw_ligand_identifier": raw_id,
        "normalized_ligand_identifier": "",
        "joined_rdk_id": "",
        "candidate_mapping_rdk_ids": _joined(resolution.candidate_rdk_ids),
        "identity_join_status": resolution.join_status,
        "identity_resolution_status": "unresolved",
        "resolved_drugcentral_id": "",
        "resolved_preferred_name": "",
        "approved_parent_inchikey": "",
        "current_pdbqt_path": "",
        "current_pdbqt_sha256": "",
        "observed_current_pdbqt_sha256": "",
        "current_file_checksum_status": "not_joined",
        "historical_pdbqt_sha256": historical_hash,
        "historical_checksum_status": _historical_checksum_status(
            "", historical_hash
        ),
        "score_reconciliation_status": "unresolved",
        "score_reconciliation_reason_codes": ";".join(resolution.reasons),
        "chemistry_redock_required": "",
    }


def _score_status(
    repaired: Mapping[str, str],
    join_status: str,
    historical_status: str,
) -> tuple[str, list[str]]:
    if join_status in {"current_artifact_missing", "current_artifact_hash_changed"}:
        return "blocked_current_file_unusable", [join_status]
    if not repaired.get("selected_pdbqt_sha256"):
        return "blocked_current_file_unusable", ["current_pdbqt_checksum_unavailable"]
    if join_status.startswith("conflict_") or historical_status == "mismatch":
        return "blocked_historical_checksum_mismatch", [
            "historical_ligand_checksum_mismatch"
        ]
    if historical_status == "missing_historical":
        return "candidate_context_unverified", [
            "rdk_id_joined",
            "current_checksum_attached",
            "historical_ligand_checksum_missing",
        ]
    if historical_status == "invalid_historical":
        return "blocked_historical_checksum_invalid", [
            "historical_ligand_checksum_invalid"
        ]
    return "checksum_match_context_incomplete", [
        "rdk_id_and_ligand_checksum_match",
        "receptor_config_engine_context_not_validated",
    ]


def _resolve_score_identity(
    raw_id: str,
    rdk_id: str,
    historical_hash: str,
    repaired_by_id: Mapping[str, Mapping[str, str]],
    repaired_by_hash: Mapping[str, Sequence[Mapping[str, str]]],
    inventory_by_id: Mapping[str, Mapping[str, Any]],
) -> _ScoreResolution:
    if not raw_id:
        return _hash_only_resolution(
            "unresolved_no_identifier", historical_hash, repaired_by_hash
        )
    if not rdk_id:
        return _hash_only_resolution(
            "unresolved_identifier_unparseable", historical_hash, repaired_by_hash
        )
    repaired = repaired_by_id.get(rdk_id)
    if repaired is None:
        return _hash_only_resolution(
            "unresolved_rdk_not_in_mapping", historical_hash, repaired_by_hash
        )
    return _rdk_score_resolution(
        rdk_id, repaired, historical_hash, repaired_by_hash, inventory_by_id
    )


def _rdk_score_resolution(
    rdk_id: str,
    repaired: Mapping[str, str],
    historical_hash: str,
    repaired_by_hash: Mapping[str, Sequence[Mapping[str, str]]],
    inventory_by_id: Mapping[str, Mapping[str, Any]],
) -> _ScoreResolution:
    inventory = inventory_by_id[rdk_id]
    observed = str(inventory.get("verified_pdbqt_sha256", ""))
    if not observed:
        return _ScoreResolution(
            repaired, rdk_id, (rdk_id,), "current_artifact_missing", ()
        )
    stated = repaired.get("selected_pdbqt_sha256", "")
    if stated != observed:
        return _ScoreResolution(
            repaired, rdk_id, (rdk_id,), "current_artifact_hash_changed", ()
        )
    if not historical_hash:
        return _ScoreResolution(
            repaired,
            rdk_id,
            (rdk_id,),
            "joined_rdk_only_missing_historical_hash",
            ("historical_ligand_checksum_missing",),
        )
    if not _is_sha256(historical_hash):
        return _ScoreResolution(
            repaired,
            rdk_id,
            (rdk_id,),
            "ambiguous_historical_hash",
            ("historical_ligand_checksum_invalid",),
        )
    hash_matches = repaired_by_hash.get(historical_hash.casefold(), ())
    hash_ids = tuple(sorted(row["rdk_id"] for row in hash_matches))
    if historical_hash.casefold() == observed.casefold():
        return _ScoreResolution(
            repaired, rdk_id, (rdk_id,), "joined_rdk_hash_match", ()
        )
    status = "conflict_rdk_vs_hash" if hash_matches else "conflict_rdk_hash_mismatch"
    candidates = tuple(sorted({rdk_id, *hash_ids}))
    return _ScoreResolution(
        repaired,
        rdk_id,
        candidates,
        status,
        ("rdk_id_and_historical_checksum_disagree",),
    )


def _hash_only_resolution(
    fallback_status: str,
    historical_hash: str,
    repaired_by_hash: Mapping[str, Sequence[Mapping[str, str]]],
) -> _ScoreResolution:
    if not historical_hash:
        return _ScoreResolution(
            None, "", (), fallback_status, ("historical_ligand_checksum_missing",)
        )
    if not _is_sha256(historical_hash):
        return _ScoreResolution(
            None,
            "",
            (),
            "ambiguous_historical_hash",
            ("historical_ligand_checksum_invalid",),
        )
    candidates = repaired_by_hash.get(historical_hash.casefold(), ())
    ids = tuple(sorted(row["rdk_id"] for row in candidates))
    if len(candidates) == 1:
        return _ScoreResolution(
            candidates[0],
            ids[0],
            ids,
            "joined_unique_hash_only",
            ("raw_ligand_identifier_unresolved",),
        )
    if len(candidates) > 1:
        return _ScoreResolution(
            None,
            "",
            ids,
            "ambiguous_historical_hash",
            ("historical_checksum_maps_to_multiple_legacy_rows",),
        )
    return _ScoreResolution(
        None,
        "",
        (),
        fallback_status,
        ("historical_checksum_not_in_legacy_inventory",),
    )


def _repaired_by_observed_hash(
    repaired_by_id: Mapping[str, Mapping[str, str]],
    inventory_by_id: Mapping[str, Mapping[str, Any]],
) -> dict[str, list[Mapping[str, str]]]:
    output: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for rdk_id, repaired in repaired_by_id.items():
        digest = str(inventory_by_id[rdk_id].get("verified_pdbqt_sha256", ""))
        if digest:
            output[digest.casefold()].append(repaired)
    return dict(output)


def _current_file_checksum_status(stated: str, observed: str) -> str:
    if not observed:
        return "current_artifact_missing"
    if stated.casefold() != observed.casefold():
        return "current_artifact_hash_changed"
    return "current_checksum_verified"


def _historical_checksum_status(current: str, historical: str) -> str:
    if not historical:
        return "missing_historical"
    if not _is_sha256(historical):
        return "invalid_historical"
    if not current:
        return "missing_current"
    return "match" if current.casefold() == historical.casefold() else "mismatch"


def _redock_linkage(
    delta_rows: Sequence[Mapping[str, str]],
    repaired_rows: Sequence[Mapping[str, str]],
    named_manifest_path: Path | None,
    score_joins: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    named = _named_manifest_by_parent(named_manifest_path)
    repaired_by_parent: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in repaired_rows:
        key = row.get("approved_parent_inchikey", "")
        if key:
            repaired_by_parent[key].append(row)
    scores_by_parent: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in score_joins:
        key = str(row.get("approved_parent_inchikey", ""))
        if key:
            scores_by_parent[key].append(row)
    return [
        _redock_row(
            row,
            repaired_by_parent.get(row["desired_parent_inchikey"], []),
            named,
            scores_by_parent.get(row["desired_parent_inchikey"], []),
        )
        for row in delta_rows
    ]


def _redock_row(
    delta: Mapping[str, str],
    linked: Sequence[Mapping[str, str]],
    named: Mapping[str, Mapping[str, str]],
    score_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    key = delta["desired_parent_inchikey"]
    named_row = named.get(key, {})
    status = named_row.get("materialization_status", "")
    output = dict(delta)
    output.update(
        {
            "linked_legacy_rdk_ids": _joined(row.get("rdk_id", "") for row in linked),
            "redock_trigger_rdk_ids": _joined(
                row.get("rdk_id", "")
                for row in linked
                if row.get("chemistry_redock_required") == "true"
            ),
            "trigger_mapping_rows": _joined(
                row.get("mapping_row_number", "")
                for row in linked
                if row.get("chemistry_redock_required") == "true"
            ),
            "legacy_selected_paths": _joined(
                row.get("selected_pdbqt_path", "") for row in linked
            ),
            "legacy_selected_sha256s": _joined(
                row.get("selected_pdbqt_sha256", "") for row in linked
            ),
            "linked_score_row_ids": _joined(
                row.get("score_row_stable_id", "") for row in score_rows
            ),
            "named_parent_materialization_status": status,
            "named_parent_pdbqt_path": named_row.get("output_pdbqt_path", ""),
            "named_parent_pdbqt_sha256": named_row.get(
                "output_pdbqt_sha256", ""
            ),
            "named_parent_status_reason_codes": named_row.get(
                "status_reason_codes", ""
            ),
            "redock_linkage_status": _redock_status(status),
            "docking_launched": "false",
        }
    )
    return output


def _named_manifest_by_parent(
    path: Path | None,
) -> dict[str, Mapping[str, str]]:
    if path is None:
        return {}
    _, rows = _read_csv(path)
    return {row["parent_inchikey"]: row for row in rows}


def _redock_status(materialization_status: str) -> str:
    if materialization_status == "ready":
        return "parent_prepared_ready_no_docking_launched"
    if materialization_status == "quarantined_non_dockable":
        return "parent_quarantined_non_dockable"
    return "parent_materialization_not_linked"


def _validate_mapping_pair(
    original: Sequence[Mapping[str, str]],
    repaired: Sequence[Mapping[str, str]],
) -> None:
    if len(original) != len(repaired):
        raise ValueError("legacy and repaired FDA mapping cardinalities differ")
    expected = [str(index) for index in range(1, len(repaired) + 1)]
    actual = [row.get("mapping_row_number", "") for row in repaired]
    if actual != expected:
        raise ValueError("repaired FDA mapping row numbers are not exactly 1..N")
    rdk_ids = [row.get("rdk_id", "") for row in repaired]
    if not all(rdk_ids) or len(set(rdk_ids)) != len(rdk_ids):
        raise ValueError("repaired FDA mapping rdk_id values are not unique")


def _validate_current_counts(
    *,
    original_rows: Sequence[Mapping[str, str]],
    changes: Sequence[Mapping[str, Any]],
    inventory: Sequence[Mapping[str, Any]],
    unresolved: Sequence[Mapping[str, Any]],
    redock: Sequence[Mapping[str, Any]],
    score_joins: Sequence[Mapping[str, Any]],
) -> None:
    if len(original_rows) != 4730:
        return
    counts = {
        "changes": len(changes),
        "usable": sum(row["preservation_status"] == "preserved_usable" for row in inventory),
        "unresolved": len(unresolved),
        "redock": len(redock),
        "score_rows": len(score_joins),
    }
    expected = {
        "changes": 945,
        "usable": 4729,
        "unresolved": 3459,
        "redock": 599,
        "score_rows": 342,
    }
    if counts != expected:
        raise RuntimeError(
            f"FDA legacy reconciliation count drift: expected={expected} actual={counts}"
        )


def _summary_payload(
    *,
    inputs: Mapping[str, Path],
    outputs: FDALegacyReconciliationOutputs,
    original_rows: Sequence[Mapping[str, str]],
    changes: Sequence[Mapping[str, Any]],
    inventory: Sequence[Mapping[str, Any]],
    score_joins: Sequence[Mapping[str, Any]],
    unresolved: Sequence[Mapping[str, Any]],
    redock: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    output_paths = _output_path_map(outputs)
    return {
        "schema_version": SCHEMA_VERSION,
        "reconciliation_version": RECONCILIATION_VERSION,
        "input_paths": {key: str(path) for key, path in inputs.items()},
        "input_sha256": {key: _sha256(path) for key, path in inputs.items()},
        "input_metadata": _input_metadata(inputs),
        "counts": _bundle_counts(
            original_rows, changes, inventory, score_joins, unresolved, redock
        ),
        "output_paths": {key: str(path) for key, path in output_paths.items()},
        "output_sha256": _output_hashes(output_paths),
        "policy": {
            "canonical_legacy_mapping_mutated": False,
            "legacy_pdbqt_files_mutated": False,
            "score_tables_mutated": False,
            "unresolved_rows_are_non_fda": False,
            "unresolved_rows_require_second_pass": True,
            "historical_checksum_missing_authorizes_score_reuse": False,
            "score_filename_establishes_fda_scope": False,
            "score_library_scope_columns_preserved": True,
            "cohort_dependent_score_statistics_require_recomputation": True,
            "docking_launched": False,
        },
    }


def _bundle_counts(
    original_rows: Sequence[Mapping[str, str]],
    changes: Sequence[Mapping[str, Any]],
    inventory: Sequence[Mapping[str, Any]],
    score_joins: Sequence[Mapping[str, Any]],
    unresolved: Sequence[Mapping[str, Any]],
    redock: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    preserved = sum(
        row["preservation_status"] == "preserved_usable" for row in inventory
    )
    return {
        "mapping_rows": len(original_rows),
        "identity_changed_rows": len(changes),
        "preserved_usable_files": preserved,
        "unusable_legacy_files": len(inventory) - preserved,
        "score_join_rows": len(score_joins),
        "score_identity_join_status": _field_counts(
            score_joins, "identity_join_status"
        ),
        "score_reconciliation_status": _field_counts(
            score_joins, "score_reconciliation_status"
        ),
        "unresolved_second_pass_rows": len(unresolved),
        "redock_linkage_rows": len(redock),
        "redock_linkage_status": _field_counts(redock, "redock_linkage_status"),
    }


def _field_counts(
    rows: Sequence[Mapping[str, Any]], field: str
) -> dict[str, int]:
    return dict(sorted(Counter(str(row[field]) for row in rows).items()))


def _output_hashes(paths: Mapping[str, Path]) -> dict[str, str]:
    return {
        key: _sha256(path)
        for key, path in paths.items()
        if key != "summary_json"
    }


def _output_path_map(
    outputs: FDALegacyReconciliationOutputs,
) -> dict[str, Path]:
    return {
        "versioned_mapping_csv": outputs.versioned_mapping_csv,
        "identity_changes_csv": outputs.identity_changes_csv,
        "file_inventory_csv": outputs.file_inventory_csv,
        "score_identity_joins_csv": outputs.score_identity_joins_csv,
        "unresolved_second_pass_csv": outputs.unresolved_second_pass_csv,
        "redock_linkage_csv": outputs.redock_linkage_csv,
        "summary_json": outputs.summary_json,
    }


def _input_metadata(inputs: Mapping[str, Path]) -> dict[str, dict[str, Any]]:
    output = {}
    for key, path in inputs.items():
        fields, rows = _read_csv(path)
        output[key] = {
            "path": str(path),
            "sha256": _sha256(path),
            "bytes": path.stat().st_size,
            "row_count": len(rows),
            "header": fields,
        }
    return output


def _redock_fields(delta_rows: Sequence[Mapping[str, str]]) -> list[str]:
    base = list(delta_rows[0]) if delta_rows else []
    return [*base, *_without_duplicates(_REDOCK_JOIN_FIELDS, base)]


def _validate_output_collisions(
    outputs: FDALegacyReconciliationOutputs, inputs: Iterable[Path]
) -> None:
    output_paths = set(_output_path_map(outputs).values())
    collisions = {path.resolve() for path in output_paths} & {
        path.resolve() for path in inputs
    }
    if len(output_paths) != 7 or collisions:
        raise ValueError(
            "FDA legacy reconciliation outputs collide with each other or inputs"
        )


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f"CSV has no header: {path}")
        return fields, [
            {key: _clean(value) for key, value in row.items()} for row in reader
        ]


def _write_csv(
    path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(tmp, path)


def _parse_score_rdk_id(value: object) -> str:
    text = Path(str(value or "").strip()).name
    lowered = text.casefold()
    for suffix in (".pdbqt", ".mol2", ".sdf", ".pdb"):
        if lowered.endswith(suffix):
            text = text[: -len(suffix)]
            break
    match = re.fullmatch(r"rdk[_-]?(\d+)", text, flags=re.IGNORECASE)
    return f"rdk_{int(match.group(1)):07d}" if match else ""


def _is_sha256(value: object) -> bool:
    return bool(re.fullmatch(r"[0-9a-fA-F]{64}", str(value or "").strip()))


def _first(row: Mapping[str, str], fields: Sequence[str]) -> str:
    return next((row.get(field, "") for field in fields if row.get(field)), "")


def _joined(values: Iterable[object]) -> str:
    return ";".join(
        sorted({str(value).strip() for value in values if str(value).strip()})
    )


def _without_duplicates(values: Sequence[str], existing: Sequence[str]) -> list[str]:
    seen = set(existing)
    output = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        output.append(value)
    return output


def _clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text.casefold() == "nan" else text


def _bool(value: bool) -> str:
    return "true" if value else "false"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "FDALegacyReconciliationOutputs",
    "legacy_reconciliation_paths",
    "reconcile_legacy_fda",
]
