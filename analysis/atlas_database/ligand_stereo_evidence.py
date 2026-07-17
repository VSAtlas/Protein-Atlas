"""Fail-closed ingestion of immutable ligand stereochemistry evidence."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import sqlite3
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

AUDIT_STATUSES = frozenset(
    {
        "no_potential_stereo",
        "fully_specified",
        "partially_unspecified",
        "all_unspecified",
        "parse_failed",
    }
)
PARSE_STATUSES = frozenset({"parsed", "parse_failed"})
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_SDF_DELIMITER_RE = re.compile(rb"(?m)^\$\$\$\$[ \t]*(?:\r?\n|\Z)")
_SDF_PROPERTY_RE = re.compile(r"^>\s*<([^>]+)>.*$")


class LigandStereoEvidenceError(ValueError):
    """Raised when stereo evidence is incomplete, ambiguous, or inconsistent."""


def _text(value: Any) -> str:
    return str(value).strip() if value is not None else ""


def _required(record: Mapping[str, Any], field: str, label: str) -> str:
    value = _text(record.get(field))
    if not value:
        raise LigandStereoEvidenceError(f"{label} requires non-empty {field}")
    return value


def _library_name(record: Mapping[str, Any], label: str) -> str:
    library_name = _text(record.get("library_name"))
    library_id = _text(record.get("library_id"))
    if library_name and library_id and library_name != library_id:
        raise LigandStereoEvidenceError(
            f"{label} has conflicting library_name and library_id"
        )
    value = library_name or library_id
    if not value:
        raise LigandStereoEvidenceError(
            f"{label} requires non-empty library_id or library_name"
        )
    return value


def _parent_ligand_id(record: Mapping[str, Any], label: str) -> str | None:
    parent = _text(record.get("parent_ligand_id"))
    canonical_parent = _text(record.get("canonical_parent_id"))
    if parent and canonical_parent and parent != canonical_parent:
        raise LigandStereoEvidenceError(
            f"{label} has conflicting parent_ligand_id and canonical_parent_id"
        )
    return parent or canonical_parent or None


def _sha256_value(record: Mapping[str, Any], field: str, label: str) -> str:
    value = _required(record, field, label).lower()
    if not _SHA256.fullmatch(value):
        raise LigandStereoEvidenceError(
            f"{label}.{field} must be 64 hexadecimal digits"
        )
    return value


def _nonnegative_int(record: Mapping[str, Any], field: str, label: str) -> int:
    value = record.get(field)
    if isinstance(value, bool) or value in (None, ""):
        raise LigandStereoEvidenceError(
            f"{label}.{field} must be a non-negative integer"
        )
    try:
        parsed = int(str(value))
    except (TypeError, ValueError) as exc:
        raise LigandStereoEvidenceError(
            f"{label}.{field} must be a non-negative integer"
        ) from exc
    if parsed < 0 or _text(value).lstrip("+") != str(parsed):
        raise LigandStereoEvidenceError(
            f"{label}.{field} must be a non-negative integer"
        )
    return parsed


def _strict_bool(value: Any, field: str, label: str) -> bool:
    if isinstance(value, bool):
        return value
    token = _text(value).lower()
    if token in {"true", "1"}:
        return True
    if token in {"false", "0"}:
        return False
    raise LigandStereoEvidenceError(
        f"{label}.{field} must be an explicit true or false"
    )


def _stereo_elements(value: Any, label: str) -> str:
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError as exc:
            raise LigandStereoEvidenceError(
                f"{label}.stereo_elements_json is not valid JSON"
            ) from exc
    if not isinstance(value, list):
        raise LigandStereoEvidenceError(
            f"{label}.stereo_elements_json must encode a list"
        )
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _json_mapping(value: Any, field: str, label: str) -> str:
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError as exc:
            raise LigandStereoEvidenceError(
                f"{label}.{field} is not valid JSON"
            ) from exc
    if not isinstance(value, Mapping):
        raise LigandStereoEvidenceError(f"{label}.{field} must encode a mapping")
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def _normalized_csv_row(
    fieldnames: Sequence[str], row: Mapping[str | None, Any]
) -> dict[str, Any]:
    normalized: dict[str, Any] = {
        field: "" if row.get(field) is None else str(row.get(field))
        for field in fieldnames
    }
    extras = row.get(None)
    if extras:
        normalized["__extra_columns__"] = [str(value) for value in extras]
    return normalized


def _csv_record_materials(
    payload: bytes,
    fields: Mapping[str, Any],
    label: str,
) -> dict[int, tuple[bytes, bytes, str]]:
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise LigandStereoEvidenceError(
            f"{label}.source_path is not valid UTF-8 CSV"
        ) from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    fieldnames = reader.fieldnames
    if fieldnames is None or len(fieldnames) != len(set(fieldnames)):
        raise LigandStereoEvidenceError(
            f"{label}.source_path has a missing or duplicate CSV header"
        )
    id_column = str(fields["id_column"])
    structure_column = str(fields["structure_column"])
    missing = [name for name in (id_column, structure_column) if name not in fieldnames]
    if missing:
        raise LigandStereoEvidenceError(
            f"{label}.source_path is missing audited CSV columns: {', '.join(missing)}"
        )
    materials: dict[int, tuple[bytes, bytes, str]] = {}
    for index, row in enumerate(reader, start=1):
        normalized = _normalized_csv_row(fieldnames, row)
        canonical_row = _canonical_json(normalized).encode("utf-8")
        structure = str(normalized[structure_column]).encode("utf-8")
        source_id = str(normalized[id_column]).strip()
        materials[index] = (canonical_row, structure, source_id)
    return materials


def _read_sdf_properties(record_text: str) -> dict[str, str]:
    lines = record_text.splitlines()
    properties: dict[str, str] = {}
    index = 0
    while index < len(lines):
        match = _SDF_PROPERTY_RE.match(lines[index])
        if match is None:
            index += 1
            continue
        key = match.group(1).strip()
        index += 1
        values: list[str] = []
        while index < len(lines) and lines[index].strip():
            if lines[index].strip() == "$$$$":
                break
            values.append(lines[index])
            index += 1
        properties[key] = "\n".join(values)
    return properties


def _sdf_structure_bytes(record_bytes: bytes) -> bytes:
    lines = record_bytes.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if line.strip() == b"M  END":
            return b"".join(lines[: index + 1])
    return record_bytes


def _sdf_record_spans(payload: bytes) -> list[bytes]:
    spans: list[bytes] = []
    start = 0
    for match in _SDF_DELIMITER_RE.finditer(payload):
        spans.append(payload[start : match.end()])
        start = match.end()
    if payload[start:].strip():
        spans.append(payload[start:])
    return spans


def _sdf_record_materials(
    payload: bytes,
    fields: Mapping[str, Any],
) -> dict[int, tuple[bytes, bytes, str]]:
    id_column = str(fields["id_column"])
    materials: dict[int, tuple[bytes, bytes, str]] = {}
    for index, record_bytes in enumerate(_sdf_record_spans(payload), start=1):
        structure = _sdf_structure_bytes(record_bytes)
        record_text = record_bytes.decode("utf-8", errors="replace")
        lines = record_text.splitlines()
        if id_column == "_Name":
            source_id = lines[0].strip() if lines else ""
        else:
            source_id = _read_sdf_properties(record_text).get(id_column, "").strip()
        materials[index] = (record_bytes, structure, source_id)
    return materials


def _source_record_material(
    fields: Mapping[str, Any],
    label: str,
    source_payload_cache: Mapping[Path, tuple[str, bytes]],
    record_material_cache: dict[
        tuple[Path, str, str, str], dict[int, tuple[bytes, bytes, str]]
    ],
) -> tuple[bytes, bytes, str]:
    path = Path(str(fields["ligand_source_path"]))
    expected_file_hash = str(fields["source_file_sha256"])
    observed_file_hash, payload = source_payload_cache[path]
    if observed_file_hash != expected_file_hash:
        raise LigandStereoEvidenceError(
            f"{label}.source_file_sha256 does not match cached source payload: "
            f"expected {expected_file_hash}, observed {observed_file_hash}"
        )
    key = (
        path,
        str(fields["source_kind"]),
        str(fields["id_column"]),
        str(fields.get("structure_column") or ""),
    )
    materials = record_material_cache.get(key)
    if materials is None:
        if fields["source_kind"] == "csv":
            materials = _csv_record_materials(payload, fields, label)
        else:
            materials = _sdf_record_materials(payload, fields)
        record_material_cache[key] = materials
    record_index = int(fields["source_record_index"])
    material = materials.get(record_index)
    if material is None:
        raise LigandStereoEvidenceError(
            f"{label}.source_record_index {record_index} is outside the source "
            f"{str(fields['source_kind']).upper()}"
        )
    return material


def _require_recomputed_hash(
    fields: Mapping[str, Any],
    field: str,
    material: bytes,
    label: str,
) -> None:
    expected = str(fields[field])
    observed = _sha256_bytes(material)
    if observed != expected:
        raise LigandStereoEvidenceError(
            f"{label}.{field} does not match source record {fields['source_record_index']}: "
            f"expected {expected}, observed {observed}"
        )


def _validate_source_record_attachment(
    fields: Mapping[str, Any],
    label: str,
    source_payload_cache: Mapping[Path, tuple[str, bytes]],
    record_material_cache: dict[
        tuple[Path, str, str, str], dict[int, tuple[bytes, bytes, str]]
    ],
) -> None:
    source_record_id = str(fields["source_record_id"])
    canonical_id = str(fields["ligand_canonical_id"])
    if not source_record_id or canonical_id != source_record_id:
        raise LigandStereoEvidenceError(
            f"{label} requires ligand_canonical_id to equal non-empty source_record_id"
        )
    evidence = json.loads(str(fields["source_record_evidence_json"]))
    evidence_id = evidence.get("id_value")
    if not isinstance(evidence_id, str) or evidence_id != source_record_id:
        raise LigandStereoEvidenceError(
            f"{label}.source_record_evidence_json.id_value must equal source_record_id"
        )
    evidence_field_name = "id_column" if fields["source_kind"] == "csv" else "id_field"
    evidence_id_field = evidence.get(evidence_field_name)
    if (
        not isinstance(evidence_id_field, str)
        or evidence_id_field != fields["id_column"]
    ):
        raise LigandStereoEvidenceError(
            f"{label}.source_record_evidence_json.{evidence_field_name} "
            "must equal id_column"
        )
    record_material, structure_material, observed_source_id = _source_record_material(
        fields, label, source_payload_cache, record_material_cache
    )
    if not observed_source_id or observed_source_id != source_record_id:
        raise LigandStereoEvidenceError(
            f"{label}.source_record_id does not match the exact source record"
        )
    _require_recomputed_hash(fields, "source_record_sha256", record_material, label)
    _require_recomputed_hash(
        fields, "audited_structure_sha256", structure_material, label
    )


def _verify_source_file(
    record: Mapping[str, Any],
    annotation_source_path: Path,
    expected_sha256: str,
    label: str,
    cache: dict[Path, tuple[str, bytes]],
) -> str:
    raw = _required(record, "source_path", label)
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = annotation_source_path.parent / path
    path = path.resolve()
    if not path.is_file():
        raise LigandStereoEvidenceError(f"{label}.source_path does not exist: {path}")
    actual = cache.get(path)
    if actual is None:
        payload = path.read_bytes()
        actual = (_sha256_bytes(payload), payload)
        cache[path] = actual
    observed_sha256, _payload = actual
    if observed_sha256 != expected_sha256:
        raise LigandStereoEvidenceError(
            f"{label}.source_file_sha256 does not match {path}: "
            f"expected {expected_sha256}, observed {observed_sha256}"
        )
    return str(path)


_EMPTY = frozenset({"empty"})
_PRESENT = frozenset({"present"})
_ANY_HASH_STATE = frozenset({"empty", "valid", "invalid"})
_PREPARED_ARTIFACT_STATE_MATRIX = {
    "not_declared": {
        "path": _EMPTY,
        "declared": _EMPTY,
        "computed": _EMPTY,
        "suffix": _EMPTY,
        "relation": frozenset({"none"}),
    },
    "missing_path_value": {
        "path": _EMPTY,
        "declared": _ANY_HASH_STATE,
        "computed": _EMPTY,
        "suffix": _EMPTY,
        "relation": frozenset({"none"}),
    },
    "unsupported_artifact_type": {
        "path": _PRESENT,
        "declared": _ANY_HASH_STATE,
        "computed": _EMPTY,
        "suffix": frozenset({"non_pdbqt"}),
        "relation": frozenset({"none"}),
    },
    "missing_file": {
        "path": _PRESENT,
        "declared": _ANY_HASH_STATE,
        "computed": _EMPTY,
        "suffix": frozenset({"pdbqt"}),
        "relation": frozenset({"none"}),
    },
    "computed_only": {
        "path": _PRESENT,
        "declared": _EMPTY,
        "computed": frozenset({"valid"}),
        "suffix": frozenset({"pdbqt"}),
        "relation": frozenset({"none"}),
    },
    "invalid_declared_sha256": {
        "path": _PRESENT,
        "declared": frozenset({"invalid"}),
        "computed": frozenset({"valid"}),
        "suffix": frozenset({"pdbqt"}),
        "relation": frozenset({"none"}),
    },
    "verified_match": {
        "path": _PRESENT,
        "declared": frozenset({"valid"}),
        "computed": frozenset({"valid"}),
        "suffix": frozenset({"pdbqt"}),
        "relation": frozenset({"equal"}),
    },
    "declared_mismatch": {
        "path": _PRESENT,
        "declared": frozenset({"valid"}),
        "computed": frozenset({"valid"}),
        "suffix": frozenset({"pdbqt"}),
        "relation": frozenset({"unequal"}),
    },
}


def _hash_state(value: str) -> str:
    if not value:
        return "empty"
    return "valid" if _SHA256.fullmatch(value) else "invalid"


def _prepared_artifact_state(
    path_text: str,
    declared: str,
    computed: str,
) -> dict[str, str]:
    path_state = "present" if path_text else "empty"
    suffix = "empty"
    if path_text:
        suffix = "pdbqt" if Path(path_text).suffix.lower() == ".pdbqt" else "non_pdbqt"
    declared_state = _hash_state(declared)
    computed_state = _hash_state(computed)
    relation = "none"
    if declared_state == computed_state == "valid":
        relation = "equal" if declared.lower() == computed.lower() else "unequal"
    return {
        "path": path_state,
        "declared": declared_state,
        "computed": computed_state,
        "suffix": suffix,
        "relation": relation,
    }


def _validate_prepared_state_matrix(
    status: str,
    actual: Mapping[str, str],
    label: str,
) -> None:
    expected = _PREPARED_ARTIFACT_STATE_MATRIX.get(status)
    if expected is None:
        raise LigandStereoEvidenceError(
            f"{label} has unknown prepared_artifact_verification_status"
        )
    mismatches = [
        f"{field}={value!r}"
        for field, value in actual.items()
        if value not in expected[field]
    ]
    if mismatches:
        raise LigandStereoEvidenceError(
            f"{label} prepared-artifact state is incompatible with {status}: "
            + ", ".join(mismatches)
        )


def _validate_prepared_artifact(
    fields: Mapping[str, Any],
    label: str,
    annotation_source_path: Path,
    cache: dict[Path, str],
) -> None:
    path_text = _text(fields.get("prepared_artifact_path"))
    declared = _text(fields.get("prepared_artifact_declared_sha256"))
    computed = _text(fields.get("prepared_artifact_computed_sha256"))
    status = _text(fields.get("prepared_artifact_verification_status"))
    actual_state = _prepared_artifact_state(path_text, declared, computed)
    _validate_prepared_state_matrix(status, actual_state, label)
    if status != "verified_match":
        return
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = annotation_source_path.parent / path
    path = path.resolve()
    if not path.is_file():
        raise LigandStereoEvidenceError(
            f"{label} verified prepared artifact is not a regular file: {path}"
        )
    observed = cache.get(path)
    if observed is None:
        observed = _file_sha256(path)
        cache[path] = observed
    if observed != computed:
        raise LigandStereoEvidenceError(
            f"{label} verified prepared-artifact hash does not match {path}: "
            f"recorded {computed}, observed {observed}"
        )


def _status_counts(fields: Mapping[str, Any]) -> tuple[int, int, int, int, int]:
    potential = int(fields["potential_stereo_count"])
    specified = int(fields["specified_stereo_count"])
    unspecified = int(fields["unspecified_stereo_count"])
    unknown = int(fields["unknown_stereo_count"])
    typed = (
        int(fields["tetrahedral_count"])
        + int(fields["double_bond_count"])
        + int(fields["other_stereo_count"])
    )
    return potential, specified, unspecified, unknown, typed


def _validate_status_counts(
    fields: Mapping[str, Any],
    label: str,
    counts: tuple[int, int, int, int, int],
) -> None:
    potential, specified, unspecified, unknown, typed = counts
    if potential != specified + unspecified + unknown:
        raise LigandStereoEvidenceError(
            f"{label} requires potential_stereo_count = specified + unspecified + unknown"
        )
    if potential != typed:
        raise LigandStereoEvidenceError(
            f"{label} requires potential_stereo_count = tetrahedral + double_bond + other"
        )
    if bool(fields["has_unspecified_potential_stereo"]) != (unspecified > 0):
        raise LigandStereoEvidenceError(
            f"{label}.has_unspecified_potential_stereo disagrees with its count"
        )
    if bool(fields["has_unresolved_potential_stereo"]) != (unspecified + unknown > 0):
        raise LigandStereoEvidenceError(
            f"{label}.has_unresolved_potential_stereo disagrees with its counts"
        )


def _parse_failed_status_is_complete(
    parse_status: str,
    audit_status: str,
    potential: int,
    label: str,
) -> bool:
    if parse_status == "parse_failed":
        if audit_status != "parse_failed" or potential != 0:
            raise LigandStereoEvidenceError(
                f"{label} parse_failed evidence must have audit_status=parse_failed "
                "and zero stereo counts"
            )
        return True
    if audit_status == "parse_failed":
        raise LigandStereoEvidenceError(
            f"{label} parsed evidence cannot have audit_status=parse_failed"
        )
    return False


def _validate_no_potential_status(
    potential: int,
    label: str,
) -> None:
    if potential != 0:
        raise LigandStereoEvidenceError(
            f"{label} no_potential_stereo requires zero potential features"
        )


def _validate_fully_specified_status(
    potential: int,
    specified: int,
    unspecified: int,
    unknown: int,
    label: str,
) -> None:
    if potential == 0 or specified != potential or unspecified or unknown:
        raise LigandStereoEvidenceError(
            f"{label} cannot claim fully_specified with unassigned features"
        )


def _validate_partially_unspecified_status(
    potential: int,
    specified: int,
    unresolved: int,
    label: str,
) -> None:
    if not (potential > 0 and specified > 0 and unresolved > 0):
        raise LigandStereoEvidenceError(
            f"{label} partially_unspecified requires assigned and unassigned features"
        )


def _validate_all_unspecified_status(
    potential: int,
    specified: int,
    unresolved: int,
    label: str,
) -> None:
    if not (potential > 0 and specified == 0 and unresolved == potential):
        raise LigandStereoEvidenceError(
            f"{label} all_unspecified requires every potential feature to be unassigned"
        )


def _validate_parsed_audit_status(
    audit_status: str,
    counts: tuple[int, int, int, int, int],
    label: str,
) -> None:
    potential, specified, unspecified, unknown, _typed = counts
    if audit_status == "no_potential_stereo":
        _validate_no_potential_status(potential, label)
    elif audit_status == "fully_specified":
        _validate_fully_specified_status(
            potential, specified, unspecified, unknown, label
        )
    elif audit_status == "partially_unspecified":
        _validate_partially_unspecified_status(
            potential, specified, unspecified + unknown, label
        )
    elif audit_status == "all_unspecified":
        _validate_all_unspecified_status(
            potential, specified, unspecified + unknown, label
        )


def _validate_status(fields: Mapping[str, Any], label: str) -> None:
    counts = _status_counts(fields)
    _validate_status_counts(fields, label, counts)
    parse_status = str(fields["parse_status"])
    audit_status = str(fields["audit_status"])
    if _parse_failed_status_is_complete(parse_status, audit_status, counts[0], label):
        return
    _validate_parsed_audit_status(audit_status, counts, label)


def _validated_record_tokens(
    record: Mapping[str, Any],
    label: str,
) -> tuple[str, str, str, str, str]:
    schema_version = _required(record, "schema_version", label)
    if schema_version != "atlas-ligand-stereo-audit-v1":
        raise LigandStereoEvidenceError(
            f"{label}.schema_version must be atlas-ligand-stereo-audit-v1"
        )
    source_kind = _required(record, "source_kind", label).lower()
    if source_kind not in {"sdf", "csv"}:
        raise LigandStereoEvidenceError(f"{label}.source_kind must be sdf or csv")
    parse_status = _required(record, "parse_status", label).lower()
    audit_status = _required(record, "audit_status", label).lower()
    if parse_status not in PARSE_STATUSES:
        raise LigandStereoEvidenceError(f"{label} has unknown parse_status")
    if audit_status not in AUDIT_STATUSES:
        raise LigandStereoEvidenceError(f"{label} has unknown audit_status")
    source_file_sha256 = _sha256_value(record, "source_file_sha256", label)
    return (
        schema_version,
        source_kind,
        parse_status,
        audit_status,
        source_file_sha256,
    )


def _optional_text(record: Mapping[str, Any], field: str) -> str | None:
    return _text(record.get(field)) or None


def _optional_lower_text(record: Mapping[str, Any], field: str) -> str | None:
    value = _text(record.get(field)).lower()
    return value or None


def _initial_prepared_record(
    record: Mapping[str, Any],
    label: str,
    index: int,
    annotation_source_path: Path,
    annotation_source_sha256: str,
    tokens: tuple[str, str, str, str, str],
) -> dict[str, Any]:
    (
        schema_version,
        source_kind,
        parse_status,
        audit_status,
        source_file_sha256,
    ) = tokens
    return {
        "library_name": _library_name(record, label),
        "parent_ligand_id": _parent_ligand_id(record, label),
        "ligand_canonical_id": _required(record, "ligand_canonical_id", label),
        "evidence_schema_version": schema_version,
        "source_kind": source_kind,
        "source_file_sha256": source_file_sha256,
        "source_record_index": _nonnegative_int(record, "source_record_index", label),
        "source_record_id": _required(record, "source_record_id", label),
        "source_record_sha256": _sha256_value(record, "source_record_sha256", label),
        "audited_structure_sha256": _sha256_value(
            record, "audited_structure_sha256", label
        ),
        "structure_column": _optional_text(record, "structure_column"),
        "id_column": _optional_text(record, "id_column"),
        "structure_format": _required(record, "structure_format", label).lower(),
        "prepared_artifact_path": _optional_text(record, "prepared_artifact_path"),
        "prepared_artifact_declared_sha256": _optional_lower_text(
            record, "prepared_artifact_declared_sha256"
        ),
        "prepared_artifact_computed_sha256": _optional_lower_text(
            record, "prepared_artifact_computed_sha256"
        ),
        "prepared_artifact_verification_status": _required(
            record, "prepared_artifact_verification_status", label
        ),
        "source_record_evidence_json": _json_mapping(
            record.get("source_record_evidence_json"),
            "source_record_evidence_json",
            label,
        ),
        "parse_status": parse_status,
        "audit_status": audit_status,
        "stereo_elements_json": _stereo_elements(
            record.get("stereo_elements_json"), label
        ),
        "rdkit_version": _required(record, "rdkit_version", label),
        "method": _required(record, "method", label),
        "canonical_isomeric_smiles": _optional_text(
            record, "canonical_isomeric_smiles"
        ),
        "canonical_smiles": _optional_text(record, "canonical_smiles"),
        "inchikey": _optional_text(record, "inchikey"),
        "error": _optional_text(record, "error"),
        "annotation_source_path": str(annotation_source_path),
        "annotation_source_sha256": annotation_source_sha256,
        "annotation_record_index": index,
        "source_record_json": json.dumps(
            record, sort_keys=True, separators=(",", ":"), default=str
        ),
    }


def _add_status_metrics(
    prepared: dict[str, Any],
    record: Mapping[str, Any],
    label: str,
) -> None:
    for field in (
        "potential_stereo_count",
        "specified_stereo_count",
        "unspecified_stereo_count",
        "unknown_stereo_count",
        "tetrahedral_count",
        "double_bond_count",
        "other_stereo_count",
    ):
        prepared[field] = _nonnegative_int(record, field, label)
    prepared["has_unresolved_potential_stereo"] = _strict_bool(
        record.get("has_unresolved_potential_stereo"),
        "has_unresolved_potential_stereo",
        label,
    )
    prepared["has_unspecified_potential_stereo"] = _strict_bool(
        record.get("has_unspecified_potential_stereo"),
        "has_unspecified_potential_stereo",
        label,
    )


def _validate_sdf_source_shape(
    prepared: Mapping[str, Any],
    label: str,
) -> None:
    structure_format = str(prepared["structure_format"])
    if structure_format != "sdf":
        raise LigandStereoEvidenceError(
            f"{label} SDF evidence requires structure_format=sdf"
        )
    if not prepared["id_column"] or prepared["structure_column"]:
        raise LigandStereoEvidenceError(
            f"{label} SDF evidence requires id_column and no structure_column"
        )


def _validate_csv_source_shape(
    prepared: Mapping[str, Any],
    label: str,
) -> None:
    structure_format = str(prepared["structure_format"])
    if structure_format not in {"smiles", "inchi", "molblock"}:
        raise LigandStereoEvidenceError(
            f"{label} CSV evidence has unsupported structure_format"
        )
    if not (prepared["structure_column"] and prepared["id_column"]):
        raise LigandStereoEvidenceError(
            f"{label} CSV evidence requires structure_column and id_column"
        )


def _validate_source_shape(prepared: Mapping[str, Any], label: str) -> None:
    if prepared["source_kind"] == "sdf":
        _validate_sdf_source_shape(prepared, label)
    else:
        _validate_csv_source_shape(prepared, label)


def _validate_record_provenance(
    prepared: dict[str, Any],
    record: Mapping[str, Any],
    label: str,
    annotation_source_path: Path,
    source_payload_cache: dict[Path, tuple[str, bytes]],
    record_material_cache: dict[
        tuple[Path, str, str, str], dict[int, tuple[bytes, bytes, str]]
    ],
    artifact_hash_cache: dict[Path, str],
) -> None:
    prepared["ligand_source_path"] = _verify_source_file(
        record,
        annotation_source_path,
        str(prepared["source_file_sha256"]),
        label,
        source_payload_cache,
    )
    _validate_source_record_attachment(
        prepared,
        label,
        source_payload_cache,
        record_material_cache,
    )
    _validate_prepared_artifact(
        prepared,
        label,
        annotation_source_path,
        artifact_hash_cache,
    )


def _add_release_selection(
    prepared: dict[str, Any],
    record: Mapping[str, Any],
    label: str,
) -> None:
    selection_present = "selected_for_release" in record and record.get(
        "selected_for_release"
    ) not in (None, "")
    prepared["selection_present"] = selection_present
    prepared["selected_for_release"] = (
        _strict_bool(record.get("selected_for_release"), "selected_for_release", label)
        if selection_present
        else None
    )


def _prepare_record(
    record: Mapping[str, Any],
    index: int,
    annotation_source_path: Path,
    annotation_source_sha256: str,
    source_payload_cache: dict[Path, tuple[str, bytes]],
    record_material_cache: dict[
        tuple[Path, str, str, str], dict[int, tuple[bytes, bytes, str]]
    ],
    artifact_hash_cache: dict[Path, str],
) -> dict[str, Any]:
    label = f"annotations.ligand_stereo record {index}"
    tokens = _validated_record_tokens(record, label)
    prepared = _initial_prepared_record(
        record,
        label,
        index,
        annotation_source_path,
        annotation_source_sha256,
        tokens,
    )
    if prepared["source_record_index"] < 1:
        raise LigandStereoEvidenceError(
            f"{label}.source_record_index must be at least 1"
        )
    _add_status_metrics(prepared, record, label)
    _validate_source_shape(prepared, label)
    _validate_record_provenance(
        prepared,
        record,
        label,
        annotation_source_path,
        source_payload_cache,
        record_material_cache,
        artifact_hash_cache,
    )
    _validate_status(prepared, label)
    _add_release_selection(prepared, record, label)
    return prepared


def _known_ligand_libraries(
    connection: sqlite3.Connection,
) -> dict[tuple[str, str], int]:
    rows = connection.execute(
        """SELECT DISTINCT rc.library, l.canonical_id, l.ligand_id
        FROM receptor_contexts rc
        JOIN pair_cells p ON p.receptor_context_id=rc.receptor_context_id
        JOIN ligands l ON l.ligand_id=p.ligand_id
        WHERE rc.library IS NOT NULL AND TRIM(rc.library) != ''"""
    ).fetchall()
    return {(str(row[0]), str(row[1])): int(row[2]) for row in rows}


def _prepared_record_key(item: Mapping[str, Any]) -> tuple[str, str]:
    return str(item["library_name"]), str(item["ligand_canonical_id"])


def _prepared_evidence_identity(item: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        item["source_file_sha256"],
        item["source_record_index"],
        item["source_record_sha256"],
        item["audited_structure_sha256"],
    )


def _group_prepared_records(
    prepared: list[dict[str, Any]],
    known: Mapping[tuple[str, str], int],
) -> dict[tuple[str, str], list[dict[str, Any]]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    evidence_attachments: dict[tuple[Any, ...], tuple[str, str]] = {}
    for item in prepared:
        key = _prepared_record_key(item)
        if key not in known:
            raise LigandStereoEvidenceError(
                "annotations.ligand_stereo references an unknown imported "
                f"ligand/library pair: {key!r}"
            )
        identity = _prepared_evidence_identity(item)
        previous_key = evidence_attachments.get(identity)
        if previous_key is not None:
            raise LigandStereoEvidenceError(
                "annotations.ligand_stereo source evidence identity is attached "
                f"more than once: {identity!r}; first={previous_key!r}, next={key!r}"
            )
        evidence_attachments[identity] = key
        groups[key].append(item)
    return groups


def _select_unique_candidate(
    key: tuple[str, str],
    candidate: dict[str, Any],
) -> None:
    if candidate["selection_present"] and not candidate["selected_for_release"]:
        raise LigandStereoEvidenceError(
            f"unique stereo evidence for {key!r} cannot be explicitly unselected"
        )
    candidate["selected_for_release"] = True
    candidate["selection_method"] = (
        "explicit_selected_for_release"
        if candidate["selection_present"]
        else "only_evidence_record"
    )


def _select_duplicate_candidates(
    key: tuple[str, str],
    candidates: list[dict[str, Any]],
) -> None:
    if any(not item["selection_present"] for item in candidates):
        raise LigandStereoEvidenceError(
            f"duplicate stereo evidence for {key!r} requires an explicit "
            "selected_for_release boolean on every candidate"
        )
    selected = [item for item in candidates if item["selected_for_release"]]
    if len(selected) != 1:
        raise LigandStereoEvidenceError(
            f"duplicate stereo evidence for {key!r} must select exactly one "
            f"candidate; observed {len(selected)}"
        )
    for item in candidates:
        item["selection_method"] = "explicit_selected_for_release"


def _select_release_candidates(
    groups: Mapping[tuple[str, str], list[dict[str, Any]]],
) -> None:
    for key, candidates in groups.items():
        if len(candidates) == 1:
            _select_unique_candidate(key, candidates[0])
        else:
            _select_duplicate_candidates(key, candidates)


def _insert_prepared_records(
    connection: sqlite3.Connection,
    prepared: list[dict[str, Any]],
    known: Mapping[tuple[str, str], int],
    sql: str,
) -> None:
    for item in prepared:
        item["ligand_id"] = known[_prepared_record_key(item)]
        item["selected_for_release"] = int(bool(item["selected_for_release"]))
        connection.execute(sql, item)


def ingest_ligand_stereo_evidence(
    connection: sqlite3.Connection,
    records: list[dict[str, Any]],
    annotation_source_path: Path,
    annotation_source_sha256: str,
) -> int:
    """Validate and ingest one explicit ligand stereo evidence source."""
    if not records:
        raise LigandStereoEvidenceError(
            "annotations.ligand_stereo source contains no evidence records"
        )
    known = _known_ligand_libraries(connection)
    source_payload_cache: dict[Path, tuple[str, bytes]] = {}
    record_material_cache: dict[
        tuple[Path, str, str, str], dict[int, tuple[bytes, bytes, str]]
    ] = {}
    artifact_hash_cache: dict[Path, str] = {}
    prepared = [
        _prepare_record(
            record,
            index,
            annotation_source_path,
            annotation_source_sha256,
            source_payload_cache,
            record_material_cache,
            artifact_hash_cache,
        )
        for index, record in enumerate(records, start=1)
    ]
    groups = _group_prepared_records(prepared, known)
    _select_release_candidates(groups)

    sql = """INSERT INTO ligand_stereo_evidence (
        ligand_id, library_name, parent_ligand_id, evidence_schema_version,
        source_kind, ligand_source_path, source_file_sha256, source_record_index,
        source_record_id, source_record_sha256, audited_structure_sha256,
        structure_column, id_column, structure_format,
        prepared_artifact_path, prepared_artifact_declared_sha256,
        prepared_artifact_computed_sha256,
        prepared_artifact_verification_status, source_record_evidence_json,
        parse_status, audit_status,
        potential_stereo_count, specified_stereo_count,
        unspecified_stereo_count, unknown_stereo_count, tetrahedral_count,
        double_bond_count, other_stereo_count,
        has_unspecified_potential_stereo, has_unresolved_potential_stereo,
        canonical_isomeric_smiles,
        canonical_smiles, inchikey, stereo_elements_json, rdkit_version,
        method, error, selected_for_release, selection_method,
        annotation_source_path, annotation_source_sha256,
        annotation_record_index, source_record_json
    ) VALUES (
        :ligand_id, :library_name, :parent_ligand_id, :evidence_schema_version,
        :source_kind, :ligand_source_path, :source_file_sha256,
        :source_record_index, :source_record_id, :source_record_sha256,
        :audited_structure_sha256, :structure_column, :id_column,
        :structure_format, :prepared_artifact_path,
        :prepared_artifact_declared_sha256,
        :prepared_artifact_computed_sha256,
        :prepared_artifact_verification_status, :source_record_evidence_json,
        :parse_status, :audit_status, :potential_stereo_count,
        :specified_stereo_count, :unspecified_stereo_count,
        :unknown_stereo_count, :tetrahedral_count, :double_bond_count,
        :other_stereo_count, :has_unspecified_potential_stereo,
        :has_unresolved_potential_stereo,
        :canonical_isomeric_smiles, :canonical_smiles, :inchikey,
        :stereo_elements_json, :rdkit_version, :method, :error,
        :selected_for_release, :selection_method, :annotation_source_path,
        :annotation_source_sha256, :annotation_record_index,
        :source_record_json
    )"""
    _insert_prepared_records(connection, prepared, known, sql)
    return len(prepared)


__all__ = [
    "AUDIT_STATUSES",
    "LigandStereoEvidenceError",
    "ingest_ligand_stereo_evidence",
]
