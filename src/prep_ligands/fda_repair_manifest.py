"""Non-destructive repair and score-reuse manifests for an FDA ligand audit."""

from __future__ import annotations

import csv
import hashlib
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


_RDK_ID_RE = re.compile(r"rdk[_-]?(\d+)", re.IGNORECASE)
_EXACT_PREPARED_STATES = frozenset(
    {"high_fidelity_exact", "coordinate_connectivity_exact"}
)
_PARENT_PREPARED_STATES = frozenset(
    {"high_fidelity_parent", "coordinate_connectivity_parent"}
)
_CONSISTENT_PREPARED_STATES = _EXACT_PREPARED_STATES | _PARENT_PREPARED_STATES
_UNSAFE_SOURCE_CLASSES = frozenset(
    {
        "additive_or_nonmedication",
        "multi_active_combination",
        "metal_complex",
        "ambiguous_mixture",
    }
)


@dataclass(frozen=True)
class FDARepairOutputs:
    repaired_mapping_csv: Path
    score_reuse_manifest_csv: Path
    redock_delta_csv: Path
    repair_quarantine_csv: Path


_REPAIRED_FIELDS = [
    "mapping_row_number",
    "rdk_id",
    "original_display_name",
    "original_mapping_vendor_id",
    "original_mapping_drugcentral_id",
    "original_mapping_cas",
    "original_mapping_exact_inchikey",
    "original_mapping_parent_inchikey",
    "audit_category",
    "audit_category_basis",
    "repair_action",
    "repair_reason_codes",
    "corrected_drugcentral_id",
    "corrected_preferred_name",
    "name_changed",
    "approved_source_record_index",
    "approved_full_form_preferred_name",
    "approved_salt_name",
    "approved_full_form_canonical_smiles",
    "approved_full_form_exact_inchikey",
    "approved_parent_smiles",
    "approved_parent_inchikey",
    "approved_full_form_fragment_count",
    "approved_source_substance_class",
    "pdbqt_candidate_paths",
    "pdbqt_candidate_sha256s",
    "selected_pdbqt_path",
    "selected_pdbqt_sha256",
    "selected_pdbqt_bytes",
    "selected_pdbqt_status",
    "prepared_identity_method",
    "prepared_connectivity_state",
    "docked_parent_inchikey",
    "approved_docked_relation",
    "chemistry_redock_required",
    "claim_grade_score_reuse_possible",
]

_SCORE_FIELDS = [
    "score_source_csv",
    "score_source_sha256",
    "score_row_number",
    "run_id",
    "pdb_id",
    "variant",
    "ph_label",
    "library",
    "run_mode",
    "engine",
    "stage",
    "raw_ligand_id",
    "rdk_id",
    "score_value",
    "current_pdbqt_path",
    "current_pdbqt_sha256",
    "historical_pdbqt_sha256",
    "historical_checksum_status",
    "identity_join_status",
    "corrected_drugcentral_id",
    "corrected_preferred_name",
    "approved_parent_inchikey",
    "chemistry_redock_required",
    "claim_grade_score_reuse_status",
    "reuse_reason_codes",
    "historical_receptor_sha256",
    "historical_config_sha256",
    "historical_engine_input_sha256",
]

_DELTA_FIELDS = [
    "desired_parent_inchikey",
    "desired_parent_smiles",
    "approved_source_record_indices",
    "approved_drugcentral_ids",
    "approved_preferred_names",
    "approved_full_form_exact_inchikeys",
    "approved_source_substance_classes",
    "redock_action",
    "redock_reason_codes",
]

_QUARANTINE_FIELDS = [
    "scope",
    "mapping_row_number",
    "rdk_id",
    "source_record_indices",
    "display_or_preferred_names",
    "parent_inchikey",
    "reason_codes",
]


def repair_output_paths(output_dir: Path) -> FDARepairOutputs:
    return FDARepairOutputs(
        repaired_mapping_csv=output_dir / "fda_mapping_repaired.csv",
        score_reuse_manifest_csv=output_dir / "fda_score_reuse_manifest.csv",
        redock_delta_csv=output_dir / "fda_redock_delta.csv",
        repair_quarantine_csv=output_dir / "fda_repair_quarantine.csv",
    )


def write_fda_repair_artifacts(
    *,
    audited_rows: Sequence[Any],
    source_index: Any,
    source_class_by_index: Mapping[int, str],
    output_dir: Path,
    score_csvs: Sequence[Path] = (),
    protected_input_paths: Sequence[Path] = (),
) -> tuple[FDARepairOutputs, dict[str, Any]]:
    """Write repaired identities and a chemistry-minimal redocking delta."""

    outputs = repair_output_paths(output_dir)
    score_sources = _score_source_provenance(score_csvs)
    _validate_output_collisions(
        outputs,
        [*protected_input_paths, *(Path(item["path"]) for item in score_sources)],
    )
    repaired = [_repaired_row(row, source_class_by_index) for row in audited_rows]
    _validate_repaired_rows(repaired, len(audited_rows))
    _write_csv(outputs.repaired_mapping_csv, repaired, _REPAIRED_FIELDS)

    score_rows, score_stats = _score_reuse_rows(
        [Path(item["path"]) for item in score_sources], repaired
    )
    _write_csv(outputs.score_reuse_manifest_csv, score_rows, _SCORE_FIELDS)

    delta, source_quarantine, parent_stats = _parent_delta_rows(
        repaired, source_index, source_class_by_index
    )
    _write_csv(outputs.redock_delta_csv, delta, _DELTA_FIELDS)

    mapping_quarantine = _mapping_quarantine_rows(repaired)
    quarantine = [*mapping_quarantine, *source_quarantine]
    _write_csv(outputs.repair_quarantine_csv, quarantine, _QUARANTINE_FIELDS)

    return outputs, {
        "schema_version": 1,
        "mapping_rows": len(repaired),
        "mapping_rows_with_selected_pdbqt": sum(
            bool(row["selected_pdbqt_sha256"]) for row in repaired
        ),
        "mapping_rows_with_corrected_drugcentral_identity": sum(
            bool(row["corrected_drugcentral_id"]) for row in repaired
        ),
        "mapping_name_changes": sum(row["name_changed"] == "true" for row in repaired),
        "mapping_chemistry_redock_required": sum(
            row["chemistry_redock_required"] == "true" for row in repaired
        ),
        "score_sources": score_sources,
        "score_join": score_stats,
        "canonical_parent_inventory": parent_stats,
        "redock_delta_unique_parents": len(delta),
        "quarantine_rows": len(quarantine),
        "outputs": {
            "repaired_mapping_csv": str(outputs.repaired_mapping_csv),
            "score_reuse_manifest_csv": str(outputs.score_reuse_manifest_csv),
            "redock_delta_csv": str(outputs.redock_delta_csv),
            "repair_quarantine_csv": str(outputs.repair_quarantine_csv),
        },
        "policy": {
            "original_mapping_mutated": False,
            "score_tables_mutated": False,
            "corrected_names_require_manifest_backed_structure_link": True,
            "redock_delta_grain": "unique_approved_parent_inchikey",
            "missing_historical_score_hashes_force_chemistry_redock": False,
            "missing_historical_score_hashes_block_claim_grade_reuse": True,
            "library_dependent_statistics_must_be_recomputed": True,
        },
    }


def _repaired_row(row: Any, source_class_by_index: Mapping[int, str]) -> dict[str, Any]:
    source_record = _approved_source_record(row)
    source_structure = _record_structure(source_record)
    source_id = _record_id(source_record)
    source_name = _record_name(source_record)
    source_class = _record_source_class(source_record, source_class_by_index)
    selection = _select_pdbqt(row.rdk_id, row.actual_paths)
    relation = _approved_docked_relation(row, source_structure, selection)
    repair_action, repair_reasons = _repair_action(
        row, source_record, selection, relation
    )
    name_changed = _approved_name_changed(source_name, row.display_name)
    docked_parent = _docked_parent_inchikey(row, selection)
    chemistry_redock = _chemistry_redock_required(
        source_record, source_class, relation
    )
    return {
        "mapping_row_number": row.row_number,
        "rdk_id": row.rdk_id,
        "original_display_name": row.display_name,
        "original_mapping_vendor_id": _clean(row.mapping_row.get("id")),
        "original_mapping_drugcentral_id": _clean(
            row.mapping_row.get("drugcentral_id")
        ),
        "original_mapping_cas": _clean(row.mapping_row.get("cas")),
        "original_mapping_exact_inchikey": row.mapping_structure.exact_inchikey,
        "original_mapping_parent_inchikey": row.mapping_structure.parent_inchikey,
        "audit_category": row.cross_reference_category,
        "audit_category_basis": row.cross_reference_basis,
        "repair_action": repair_action,
        "repair_reason_codes": ";".join(repair_reasons),
        "corrected_drugcentral_id": source_id,
        "corrected_preferred_name": source_name,
        "name_changed": _bool(name_changed),
        "approved_source_record_index": _record_index(source_record),
        "approved_full_form_preferred_name": source_name,
        "approved_salt_name": _approved_salt_name(
            source_name, source_structure, source_class
        ),
        "approved_full_form_canonical_smiles": _structure_value(
            source_structure, "canonical_smiles"
        ),
        "approved_full_form_exact_inchikey": _structure_value(
            source_structure, "exact_inchikey"
        ),
        "approved_parent_smiles": _structure_value(
            source_structure, "parent_smiles"
        ),
        "approved_parent_inchikey": _structure_value(
            source_structure, "parent_inchikey"
        ),
        "approved_full_form_fragment_count": _structure_value(
            source_structure, "fragment_count"
        ),
        "approved_source_substance_class": source_class,
        "pdbqt_candidate_paths": ";".join(selection["paths"]),
        "pdbqt_candidate_sha256s": ";".join(selection["hashes"]),
        "selected_pdbqt_path": selection["path"],
        "selected_pdbqt_sha256": selection["sha256"],
        "selected_pdbqt_bytes": selection["bytes"],
        "selected_pdbqt_status": selection["status"],
        "prepared_identity_method": _prepared_identity_method(row),
        "prepared_connectivity_state": row.prepared_form_state,
        "docked_parent_inchikey": docked_parent,
        "approved_docked_relation": relation,
        "chemistry_redock_required": _bool(chemistry_redock),
        "claim_grade_score_reuse_possible": "unknown_requires_historical_context_hashes",
    }


def _approved_source_record(row: Any) -> Any | None:
    if not row.approval_verified:
        return None
    return row.source_record


def _record_structure(record: Any | None) -> Any | None:
    return record.structure if record is not None else None


def _record_name(record: Any | None) -> str:
    return str(record.name or "") if record is not None else ""


def _record_index(record: Any | None) -> int | str:
    return record.index if record is not None else ""


def _record_source_class(
    record: Any | None, source_class_by_index: Mapping[int, str]
) -> str:
    if record is None:
        return ""
    return source_class_by_index.get(int(record.index), "")


def _structure_value(structure: Any | None, name: str) -> Any:
    return getattr(structure, name) if structure is not None else ""


def _approved_name_changed(source_name: str, display_name: str) -> bool:
    return bool(
        source_name
        and _normalized_name(source_name) != _normalized_name(display_name)
    )


def _docked_parent_inchikey(row: Any, selection: Mapping[str, Any]) -> str:
    if (
        selection["sha256"]
        and row.prepared_form_state in _CONSISTENT_PREPARED_STATES
    ):
        return row.mapping_structure.parent_inchikey
    return ""


def _chemistry_redock_required(
    source_record: Any | None, source_class: str, relation: str
) -> bool:
    return bool(
        source_record is not None
        and source_class not in _UNSAFE_SOURCE_CLASSES
        and relation not in {"approved_full_form_docked", "approved_parent_docked"}
    )


def _approved_salt_name(
    source_name: str, source_structure: Any | None, source_class: str
) -> str:
    if source_structure is None:
        return ""
    if source_structure.fragment_count > 1 or source_class == "salt_or_solvate":
        return source_name
    return ""


def _select_pdbqt(rdk_id: str, raw_paths: Sequence[str]) -> dict[str, Any]:
    paths = sorted({str(Path(path).expanduser().resolve()) for path in raw_paths if path})
    usable: list[tuple[Path, str, int]] = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.is_file():
            continue
        size = path.stat().st_size
        if size <= 0:
            continue
        usable.append((path, _sha256(path), size))
    hashes = sorted({item[1] for item in usable})
    if not usable:
        return {
            "paths": paths,
            "hashes": [],
            "path": "",
            "sha256": "",
            "bytes": "",
            "status": "missing_or_empty",
        }
    if len(hashes) > 1:
        return {
            "paths": paths,
            "hashes": hashes,
            "path": "",
            "sha256": "",
            "bytes": "",
            "status": "ambiguous_different_bytes",
        }
    preferred = sorted(
        usable,
        key=lambda item: (
            item[0].stem.casefold() != rdk_id.casefold(),
            len(str(item[0])),
            str(item[0]),
        ),
    )[0]
    return {
        "paths": paths,
        "hashes": hashes,
        "path": str(preferred[0]),
        "sha256": preferred[1],
        "bytes": preferred[2],
        "status": "selected_unique_bytes"
        if len(usable) == 1
        else "selected_identical_aliases",
    }


def _approved_docked_relation(
    row: Any, source_structure: Any | None, selection: Mapping[str, Any]
) -> str:
    if source_structure is None:
        return "unverifiable"
    if not selection["sha256"]:
        return "ambiguous" if selection["status"].startswith("ambiguous") else "unverifiable"
    if row.prepared_form_state not in _CONSISTENT_PREPARED_STATES:
        return "different"
    mapping = row.mapping_structure
    if (
        row.prepared_form_state in _EXACT_PREPARED_STATES
        and mapping.exact_inchikey
        and mapping.exact_inchikey == source_structure.exact_inchikey
    ):
        return "approved_full_form_docked"
    if (
        mapping.parent_inchikey
        and mapping.parent_inchikey == source_structure.parent_inchikey
    ):
        return "approved_parent_docked"
    return "different"


def _repair_action(
    row: Any,
    source_record: Any | None,
    selection: Mapping[str, Any],
    relation: str,
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    if source_record is None:
        reasons.append("manifest_fda_structure_not_uniquely_linked")
        return "retain_original_for_review", reasons
    if not selection["sha256"]:
        reasons.append(selection["status"])
    if relation not in {"approved_full_form_docked", "approved_parent_docked"}:
        reasons.append(f"approved_docked_relation_{relation}")
    if _normalized_name(source_record.name) != _normalized_name(row.display_name):
        reasons.append("preferred_name_reassigned_from_structure_link")
        action = "repair_name_and_identity"
    else:
        action = "confirm_manifest_identity"
    return action, reasons


def _prepared_identity_method(row: Any) -> str:
    if row.actual_structures:
        return "high_fidelity_pdbqt_or_sidecar_structure"
    if row.actual_connectivity:
        return "coordinate_graph_plus_mapping_structure"
    return "unavailable"


def _record_id(record: Any | None) -> str:
    if record is None:
        return ""
    return _clean(record.identifiers.get("id"))


def _score_reuse_rows(
    score_csvs: Sequence[Path], repaired: Sequence[Mapping[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    by_id = {str(row["rdk_id"]): row for row in repaired}
    output: list[dict[str, Any]] = []
    stats: defaultdict[str, int] = defaultdict(int)
    for raw_path in sorted({Path(path).expanduser().resolve() for path in score_csvs}):
        output.extend(_score_rows_from_csv(raw_path, by_id, stats))
    stats["joined_rows"] = len(output)
    return output, dict(sorted(stats.items()))


def _score_rows_from_csv(
    path: Path,
    by_id: Mapping[str, Mapping[str, Any]],
    stats: defaultdict[str, int],
) -> list[dict[str, Any]]:
    if not path.is_file():
        stats["missing_score_csvs"] += 1
        return []
    output = []
    source_hash = _sha256(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row_number, score_row in enumerate(reader, start=1):
            score_reuse_row = _score_reuse_row(
                path, source_hash, row_number, score_row, by_id
            )
            if score_reuse_row is None:
                stats["rows_without_rdk_id"] += 1
                continue
            output.append(score_reuse_row)
            stats[str(score_reuse_row["claim_grade_score_reuse_status"])] += 1
    return output


def _score_reuse_row(
    path: Path,
    source_hash: str,
    row_number: int,
    score_row: Mapping[str, Any],
    by_id: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any] | None:
    raw_ligand = _score_ligand_value(score_row)
    rdk_id = _rdk_id(raw_ligand)
    if not rdk_id:
        return None
    repair = by_id.get(rdk_id)
    context = _score_repair_context(repair)
    historical_hash = _first_value(
        score_row,
        "ligand_pdbqt_sha256",
        "pdbqt_sha256",
        "ligand_sha256",
        "input_ligand_sha256",
    )
    checksum_status = _checksum_status(context["current_pdbqt_sha256"], historical_hash)
    status, reasons = _score_reuse_status(repair, checksum_status, score_row)
    return {
        "score_source_csv": str(path),
        "score_source_sha256": source_hash,
        "score_row_number": row_number,
        "run_id": _first_value(score_row, "run_id"),
        "pdb_id": _first_value(score_row, "pdb_id", "target_id"),
        "variant": _first_value(score_row, "variant"),
        "ph_label": _first_value(score_row, "ph_label", "ph"),
        "library": _first_value(score_row, "library"),
        "run_mode": _first_value(score_row, "run_mode"),
        "engine": _infer_engine(path, score_row),
        "stage": _first_value(score_row, "stage", "selected_stage"),
        "raw_ligand_id": raw_ligand,
        "rdk_id": rdk_id,
        "score_value": _score_value(score_row),
        **context,
        "historical_pdbqt_sha256": historical_hash,
        "historical_checksum_status": checksum_status,
        "claim_grade_score_reuse_status": status,
        "reuse_reason_codes": ";".join(reasons),
        "historical_receptor_sha256": _first_value(
            score_row, "receptor_pdbqt_sha256", "receptor_sha256"
        ),
        "historical_config_sha256": _first_value(
            score_row, "docking_config_sha256", "config_sha256"
        ),
        "historical_engine_input_sha256": _first_value(
            score_row, "engine_input_sha256", "mol2_sha256"
        ),
    }


def _score_repair_context(
    repair: Mapping[str, Any] | None,
) -> dict[str, str]:
    if repair is None:
        return {
            "current_pdbqt_path": "",
            "current_pdbqt_sha256": "",
            "identity_join_status": "rdk_id_not_in_mapping",
            "corrected_drugcentral_id": "",
            "corrected_preferred_name": "",
            "approved_parent_inchikey": "",
            "chemistry_redock_required": "",
        }
    return {
        "current_pdbqt_path": str(repair["selected_pdbqt_path"]),
        "current_pdbqt_sha256": str(repair["selected_pdbqt_sha256"]),
        "identity_join_status": "joined_unique_rdk_id",
        "corrected_drugcentral_id": str(repair["corrected_drugcentral_id"]),
        "corrected_preferred_name": str(repair["corrected_preferred_name"]),
        "approved_parent_inchikey": str(repair["approved_parent_inchikey"]),
        "chemistry_redock_required": str(repair["chemistry_redock_required"]),
    }


def _score_reuse_status(
    repair: Mapping[str, Any] | None,
    checksum_status: str,
    score_row: Mapping[str, Any],
) -> tuple[str, list[str]]:
    if repair is None:
        return "not_joined", ["rdk_id_not_in_repaired_mapping"]
    if not repair["corrected_drugcentral_id"]:
        return "unresolved_identity", ["manifest_fda_identity_not_resolved"]
    if repair["chemistry_redock_required"] == "true":
        return "redock_required_structure", ["prepared_structure_not_reusable"]
    if checksum_status == "mismatch":
        return "redock_required_ligand_bytes", ["historical_ligand_sha256_mismatch"]
    if checksum_status == "missing_historical":
        return "candidate_context_unverified", ["missing_historical_ligand_sha256"]
    missing = []
    if not _first_value(score_row, "receptor_pdbqt_sha256", "receptor_sha256"):
        missing.append("missing_historical_receptor_sha256")
    if not _first_value(score_row, "docking_config_sha256", "config_sha256"):
        missing.append("missing_historical_config_sha256")
    if missing:
        return "candidate_context_unverified", missing
    return "provenance_candidate_requires_pose_validation", [
        "ligand_receptor_config_hashes_present"
    ]


def _parent_delta_rows(
    repaired: Sequence[Mapping[str, Any]],
    source_index: Any,
    source_class_by_index: Mapping[int, str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    if not bool(source_index.approval_provenance_verified):
        return _unverified_parent_delta(source_index)
    reusable = _reusable_parent_rows(repaired)
    groups, source_quarantine = _approved_parent_groups(
        source_index, source_class_by_index
    )
    delta, covered = _missing_parent_rows(
        groups, reusable, source_class_by_index
    )
    return delta, source_quarantine, {
        "eligible_unique_parents": len(groups),
        "covered_by_existing_pdbqt": covered,
        "missing_parent_redock_delta": len(delta),
        "unsafe_or_unresolved_source_records": len(source_quarantine),
        "approval_manifest_verified": 1,
    }


def _unverified_parent_delta(
    source_index: Any,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    quarantine = [
        {
            "scope": "approved_source_record",
            "mapping_row_number": "",
            "rdk_id": "",
            "source_record_indices": record.index,
            "display_or_preferred_names": record.name,
            "parent_inchikey": record.structure.parent_inchikey,
            "reason_codes": "approval_manifest_unverified",
        }
        for record in source_index.by_index.values()
    ]
    return [], quarantine, {
        "eligible_unique_parents": 0,
        "covered_by_existing_pdbqt": 0,
        "missing_parent_redock_delta": 0,
        "unsafe_or_unresolved_source_records": len(quarantine),
        "approval_manifest_verified": 0,
    }


def _reusable_parent_rows(
    repaired: Sequence[Mapping[str, Any]],
) -> dict[str, list[Mapping[str, Any]]]:
    reusable: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in repaired:
        if (
            row["approved_parent_inchikey"]
            and row["selected_pdbqt_sha256"]
            and row["approved_docked_relation"]
            in {"approved_full_form_docked", "approved_parent_docked"}
        ):
            reusable[str(row["approved_parent_inchikey"])].append(row)
    return reusable


def _approved_parent_groups(
    source_index: Any, source_class_by_index: Mapping[int, str]
) -> tuple[dict[str, list[Any]], list[dict[str, Any]]]:
    groups: dict[str, list[Any]] = defaultdict(list)
    source_quarantine: list[dict[str, Any]] = []
    for record in source_index.by_index.values():
        source_class = source_class_by_index.get(int(record.index), "")
        parent_key = record.structure.parent_inchikey
        if source_class in _UNSAFE_SOURCE_CLASSES or not parent_key:
            source_quarantine.append(
                {
                    "scope": "approved_source_record",
                    "mapping_row_number": "",
                    "rdk_id": "",
                    "source_record_indices": record.index,
                    "display_or_preferred_names": record.name,
                    "parent_inchikey": parent_key,
                    "reason_codes": "unsafe_parent_collapse"
                    if source_class in _UNSAFE_SOURCE_CLASSES
                    else "source_parent_inchikey_unavailable",
                }
            )
            continue
        groups[parent_key].append(record)
    return groups, source_quarantine


def _missing_parent_rows(
    groups: Mapping[str, Sequence[Any]],
    reusable: Mapping[str, Sequence[Mapping[str, Any]]],
    source_class_by_index: Mapping[int, str],
) -> tuple[list[dict[str, Any]], int]:
    delta: list[dict[str, Any]] = []
    covered = 0
    for parent_key, records in sorted(groups.items()):
        if reusable.get(parent_key):
            covered += 1
            continue
        delta.append(_parent_delta_row(parent_key, records, source_class_by_index))
    return delta, covered


def _parent_delta_row(
    parent_key: str,
    records: Sequence[Any],
    source_class_by_index: Mapping[int, str],
) -> dict[str, Any]:
    return {
        "desired_parent_inchikey": parent_key,
        "desired_parent_smiles": _unique_values(
            record.structure.parent_smiles for record in records
        ),
        "approved_source_record_indices": ";".join(
            str(record.index) for record in records
        ),
        "approved_drugcentral_ids": _unique_values(
            _record_id(record) for record in records
        ),
        "approved_preferred_names": _unique_values(
            record.name for record in records
        ),
        "approved_full_form_exact_inchikeys": _unique_values(
            record.structure.exact_inchikey for record in records
        ),
        "approved_source_substance_classes": _unique_values(
            source_class_by_index.get(int(record.index), "") for record in records
        ),
        "redock_action": "prepare_parent_and_dock",
        "redock_reason_codes": "no_compatible_existing_pdbqt_for_approved_parent",
    }


def _mapping_quarantine_rows(
    repaired: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    output = []
    for row in repaired:
        reasons: list[str] = []
        if not row["corrected_drugcentral_id"]:
            reasons.append("manifest_fda_identity_not_resolved")
        if str(row["selected_pdbqt_status"]).startswith(("missing", "ambiguous")):
            reasons.append(str(row["selected_pdbqt_status"]))
        if row["approved_source_substance_class"] in _UNSAFE_SOURCE_CLASSES:
            reasons.append("unsafe_parent_collapse")
        if not reasons:
            continue
        output.append(
            {
                "scope": "mapping_row",
                "mapping_row_number": row["mapping_row_number"],
                "rdk_id": row["rdk_id"],
                "source_record_indices": row["approved_source_record_index"],
                "display_or_preferred_names": row["original_display_name"],
                "parent_inchikey": row["approved_parent_inchikey"],
                "reason_codes": ";".join(dict.fromkeys(reasons)),
            }
        )
    return output


def _validate_repaired_rows(rows: Sequence[Mapping[str, Any]], expected: int) -> None:
    if len(rows) != expected:
        raise RuntimeError("FDA repaired mapping cardinality changed")
    row_numbers = [int(row["mapping_row_number"]) for row in rows]
    if row_numbers != list(range(1, expected + 1)):
        raise RuntimeError("FDA repaired mapping row numbers are not exactly 1..N")
    rdk_ids = [str(row["rdk_id"]) for row in rows]
    if not all(rdk_ids) or len(set(rdk_ids)) != expected:
        raise RuntimeError("FDA repaired mapping requires unique nonempty rdk_id values")


def _score_source_provenance(score_csvs: Sequence[Path]) -> list[dict[str, str]]:
    sources = []
    for path in sorted({Path(item).expanduser().resolve() for item in score_csvs}):
        if not path.is_file():
            raise ValueError(f"FDA score CSV not found: {path}")
        sources.append({"path": str(path), "sha256": _sha256(path)})
    return sources


def _validate_output_collisions(
    outputs: FDARepairOutputs, input_paths: Sequence[Path]
) -> None:
    output_paths = {
        path.resolve()
        for path in (
            outputs.repaired_mapping_csv,
            outputs.score_reuse_manifest_csv,
            outputs.redock_delta_csv,
            outputs.repair_quarantine_csv,
        )
    }
    if len(output_paths) != 4:
        raise ValueError("FDA repair output paths collide with each other")
    collisions = output_paths & {
        Path(path).expanduser().resolve() for path in input_paths
    }
    if collisions:
        joined = ", ".join(str(path) for path in sorted(collisions))
        raise ValueError(f"FDA repair output path collides with an input: {joined}")


def _write_csv(
    path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rdk_id(value: object) -> str:
    match = _RDK_ID_RE.search(str(value or ""))
    return f"rdk_{int(match.group(1)):07d}" if match else ""


def _score_ligand_value(row: Mapping[str, Any]) -> str:
    return _first_value(
        row, "rdk_id", "ligand_base", "ligand", "ligand_file", "Ligand_ID"
    )


def _score_value(row: Mapping[str, Any]) -> str:
    return _first_value(
        row,
        "score",
        "selected_docking_score",
        "final_score",
        "consensus_score",
        "gnina_minimized_affinity_kcal",
        "ledock_best_score_kcal",
        "dock6_best_score_kcal",
        "grid_score",
    )


def _infer_engine(path: Path, row: Mapping[str, Any]) -> str:
    explicit = _first_value(row, "engine", "best_engine")
    if explicit:
        return explicit
    name = path.name.casefold()
    for engine in ("gnina", "ledock", "dock6", "vina"):
        if engine in name:
            return engine
    return "derived_or_unknown"


def _checksum_status(current_hash: str, historical_hash: str) -> str:
    if not historical_hash:
        return "missing_historical"
    if not current_hash:
        return "missing_current"
    return "match" if current_hash.casefold() == historical_hash.casefold() else "mismatch"


def _first_value(row: Mapping[str, Any], *names: str) -> str:
    for name in names:
        value = _clean(row.get(name))
        if value:
            return value
    return ""


def _unique_values(values: Sequence[str] | Any) -> str:
    return ";".join(sorted({_clean(value) for value in values if _clean(value)}))


def _normalized_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", _clean(value).casefold()).strip()


def _clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text.casefold() == "nan" else text


def _bool(value: bool) -> str:
    return "true" if value else "false"


__all__ = [
    "FDARepairOutputs",
    "repair_output_paths",
    "write_fda_repair_artifacts",
]
