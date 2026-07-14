"""Explicit scientific annotation ingestion for Atlas release snapshots."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sqlite3
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

ANNOTATION_SOURCE_NAMES = ("proteins", "receptors", "known_pairs")
EXPLICIT_QUALIFIED_STATUSES = frozenset(
    {"qualified", "explicitly_qualified", "qualification_passed"}
)
CONTROLLED_RECEPTOR_CLASSIFICATIONS = frozenset({"APO", "HOLO"})
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")


class ScientificAnnotationError(ValueError):
    """Raised when explicit scientific annotations are ambiguous or invalid."""


def _text(value: Any) -> str:
    return str(value).strip() if value is not None else ""


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_explicitly_qualified(value: Any) -> bool:
    """Return true only for a controlled, explicitly stored pass token."""
    return _text(value).lower() in EXPLICIT_QUALIFIED_STATUSES


def _source_path(spec: Any, manifest_path: Path, name: str) -> Path:
    if isinstance(spec, str):
        raw = spec.strip()
    elif isinstance(spec, Mapping):
        raw = _text(spec.get("path"))
    else:
        raw = ""
    if not raw:
        raise ScientificAnnotationError(
            f"annotations.{name} must be a path string or mapping with path"
        )
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path.resolve()


def _load_records(path: Path, name: str) -> list[dict[str, Any]]:
    if not path.is_file():
        raise ScientificAnnotationError(
            f"annotations.{name} source does not exist: {path}"
        )
    suffix = path.suffix.lower()
    try:
        if suffix == ".csv":
            with path.open("r", encoding="utf-8", newline="") as handle:
                records: Any = list(csv.DictReader(handle))
        elif suffix == ".json":
            records = json.loads(path.read_text(encoding="utf-8"))
        elif suffix in {".yaml", ".yml"}:
            records = yaml.safe_load(path.read_text(encoding="utf-8"))
        else:
            raise ScientificAnnotationError(
                f"annotations.{name} must use .csv, .json, .yaml, or .yml"
            )
    except (OSError, csv.Error, json.JSONDecodeError, yaml.YAMLError) as exc:
        raise ScientificAnnotationError(
            f"could not read annotations.{name} source {path}: {exc}"
        ) from exc
    if isinstance(records, Mapping):
        records = records.get("records")
    if not isinstance(records, list):
        raise ScientificAnnotationError(
            f"annotations.{name} source must contain a list or records list"
        )
    result: list[dict[str, Any]] = []
    for index, record in enumerate(records, start=1):
        if not isinstance(record, Mapping):
            raise ScientificAnnotationError(
                f"annotations.{name} record {index} must be a mapping"
            )
        result.append(dict(record))
    return result


def _required(record: Mapping[str, Any], field: str, label: str) -> str:
    value = _text(record.get(field))
    if not value:
        raise ScientificAnnotationError(f"{label} requires non-empty {field}")
    return value


def _normalize_ph(value: Any) -> str:
    token = _text(value)
    if token.lower() in {"base", "none", "null"} or token.lower().endswith(".csv"):
        return ""
    return token


def _context_key(record: Mapping[str, Any], label: str) -> tuple[str, str, str, str]:
    run_id = _required(record, "run_id", label)
    pdb_id = _required(record, "pdb_id", label).upper()
    if "variant" not in record:
        raise ScientificAnnotationError(f"{label} requires explicit variant")
    if "ph_label" not in record:
        raise ScientificAnnotationError(f"{label} requires explicit ph_label")
    return (
        run_id,
        pdb_id,
        _text(record.get("variant")).upper(),
        _normalize_ph(record.get("ph_label")),
    )


def _unique_records(
    records: list[dict[str, Any]],
    key_fn: Any,
    name: str,
) -> list[tuple[Any, int, dict[str, Any]]]:
    seen: dict[Any, int] = {}
    result: list[tuple[Any, int, dict[str, Any]]] = []
    for index, record in enumerate(records, start=1):
        key = key_fn(record, f"annotations.{name} record {index}")
        if key in seen:
            raise ScientificAnnotationError(
                f"duplicate annotations.{name} key {key!r} in records "
                f"{seen[key]} and {index}"
            )
        seen[key] = index
        result.append((key, index, record))
    return result


def _provenance(record: Mapping[str, Any]) -> str | None:
    value = record.get("provenance")
    return _json(value) if value not in (None, "", {}, []) else None


def _context_id(
    connection: sqlite3.Connection, key: tuple[str, str, str, str], label: str
) -> int:
    rows = connection.execute(
        """SELECT receptor_context_id FROM receptor_contexts
        WHERE run_id=? AND pdb_id=? AND variant=? AND ph_label=?""",
        key,
    ).fetchall()
    if not rows:
        raise ScientificAnnotationError(
            f"{label} does not match an imported receptor context: {key!r}"
        )
    if len(rows) != 1:
        raise ScientificAnnotationError(
            f"{label} ambiguously matches {len(rows)} receptor contexts: {key!r}"
        )
    return int(rows[0][0])


def _optional_float(value: Any, field: str, label: str) -> float | None:
    if value in (None, ""):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ScientificAnnotationError(
            f"{label} has invalid {field}: {value!r}"
        ) from exc
    if parsed != parsed or parsed in (float("inf"), float("-inf")):
        raise ScientificAnnotationError(f"{label} has non-finite {field}: {value!r}")
    return parsed


def _optional_nonnegative_int(value: Any, field: str, label: str) -> int | None:
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        raise ScientificAnnotationError(f"{label} has invalid {field}: {value!r}")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ScientificAnnotationError(
            f"{label} has invalid {field}: {value!r}"
        ) from exc
    if parsed < 0:
        raise ScientificAnnotationError(f"{label} requires non-negative {field}")
    return parsed


def _optional_bool(value: Any, field: str, label: str) -> int | None:
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        return int(value)
    token = _text(value).lower()
    if token in {"1", "true", "yes"}:
        return 1
    if token in {"0", "false", "no"}:
        return 0
    raise ScientificAnnotationError(f"{label} has invalid {field}: {value!r}")


def _optional_sha256(value: Any, field: str, label: str) -> str | None:
    token = _text(value).lower()
    if not token:
        return None
    if not _SHA256.fullmatch(token):
        raise ScientificAnnotationError(
            f"{label}.{field} must be 64 hexadecimal digits"
        )
    return token


def _evidence_value(
    field: str,
    record: Mapping[str, Any],
    chemistry: Mapping[str, Any],
    provenance: Mapping[str, Any],
) -> Any:
    for source in (record, chemistry, provenance):
        value = source.get(field)
        if value not in (None, ""):
            return value
    return None


def _receptor_evidence_contract(
    record: Mapping[str, Any], label: str
) -> tuple[str | None, dict[str, Any], list[dict[str, Any]]]:
    provenance_value = record.get("provenance")
    if provenance_value in (None, "", {}):
        provenance: Mapping[str, Any] = {}
    elif isinstance(provenance_value, Mapping):
        provenance = provenance_value
    else:
        raise ScientificAnnotationError(f"{label}.provenance must be a mapping")

    chemistry_value = provenance.get("prepared_receptor_chemistry")
    if chemistry_value in (None, "", {}):
        chemistry: Mapping[str, Any] = {}
    elif isinstance(chemistry_value, Mapping):
        chemistry = chemistry_value
    else:
        raise ScientificAnnotationError(
            f"{label}.provenance.prepared_receptor_chemistry must be a mapping"
        )

    classification = _text(record.get("receptor_classification")).upper()
    if classification and classification not in CONTROLLED_RECEPTOR_CLASSIFICATIONS:
        allowed = ", ".join(sorted(CONTROLLED_RECEPTOR_CLASSIFICATIONS))
        raise ScientificAnnotationError(
            f"{label}.receptor_classification must be one of: {allowed}"
        )
    chemistry_classification = _text(
        chemistry.get("receptor_classification")
        or provenance.get("observed_receptor_classification")
    ).upper()
    if (
        classification
        and chemistry_classification
        and classification != chemistry_classification
    ):
        raise ScientificAnnotationError(
            f"{label} receptor classification disagrees with chemistry evidence"
        )

    evidence_status = _text(
        _evidence_value("chemistry_evidence_status", record, chemistry, provenance)
    )
    if evidence_status == "observed" and not classification:
        raise ScientificAnnotationError(
            f"{label} has observed chemistry evidence but no controlled classification"
        )

    fields = {
        "classification_method": _text(
            _evidence_value("classification_method", record, chemistry, provenance)
        )
        or None,
        "chemistry_evidence_status": evidence_status or None,
        "prepared_receptor_sha256": _optional_sha256(
            _evidence_value("prepared_receptor_sha256", record, chemistry, provenance),
            "prepared_receptor_sha256",
            label,
        ),
        "canonical_chemistry_policy_sha256": _optional_sha256(
            _evidence_value(
                "canonical_chemistry_policy_sha256", record, chemistry, provenance
            ),
            "canonical_chemistry_policy_sha256",
            label,
        ),
        "retained_metal_atom_count": _optional_nonnegative_int(
            _evidence_value("retained_metal_atom_count", record, chemistry, provenance),
            "retained_metal_atom_count",
            label,
        ),
        "retained_cofactor_residue_count": _optional_nonnegative_int(
            _evidence_value(
                "retained_cofactor_residue_count", record, chemistry, provenance
            ),
            "retained_cofactor_residue_count",
            label,
        ),
        "requested_observed_conflict": _optional_bool(
            _evidence_value(
                "requested_observed_conflict", record, chemistry, provenance
            ),
            "requested_observed_conflict",
            label,
        ),
    }

    audits: list[dict[str, Any]] = []
    identities: set[tuple[str, str]] = set()
    audit_sources = (
        ("structured_receptor_audits", False),
        ("structured_receptor_audit_rejections", True),
    )
    for source_name, rejected in audit_sources:
        audits_value = provenance.get(source_name)
        if audits_value in (None, ""):
            continue
        if not isinstance(audits_value, list):
            raise ScientificAnnotationError(
                f"{label}.provenance.{source_name} must be a list"
            )
        for audit_index, audit_value in enumerate(audits_value, start=1):
            audit_label = (
                f"{label}.provenance.{source_name}[{audit_index}]"
            )
            if not isinstance(audit_value, Mapping):
                raise ScientificAnnotationError(f"{audit_label} must be a mapping")
            audit = dict(audit_value)
            audit_kind = _required(audit, "audit_kind", audit_label)
            source_path = _required(audit, "source_path", audit_label)
            source_sha256 = _optional_sha256(
                audit.get("source_sha256"), "source_sha256", audit_label
            )
            if source_sha256 is None:
                raise ScientificAnnotationError(
                    f"{audit_label}.source_sha256 is required"
                )
            if rejected:
                rejection_reason = _required(
                    audit, "rejection_reason", audit_label
                )
                parse_status = f"rejected:{rejection_reason}"
            else:
                parse_status = _required(audit, "parse_status", audit_label)
            identity = (audit_kind, source_sha256)
            if identity in identities:
                raise ScientificAnnotationError(
                    f"{audit_label} duplicates audit kind/hash {identity!r}"
                )
            identities.add(identity)
            audits.append(
                {
                    "audit_kind": audit_kind,
                    "source_path": source_path,
                    "source_sha256": source_sha256,
                    "parse_status": parse_status,
                    "audit_json": _json(audit),
                }
            )
    return classification or None, fields, audits


def _ingest_proteins(
    connection: sqlite3.Connection,
    records: list[dict[str, Any]],
    source_path: Path,
    source_sha256: str,
) -> int:
    prepared = _unique_records(
        records,
        lambda record, label: _required(record, "protein_key", label),
        "proteins",
    )
    for protein_key, index, record in prepared:
        identity_values = (
            _text(record.get("uniprot_id")),
            _text(record.get("gene_symbol")),
            _text(record.get("display_name")),
        )
        if not any(identity_values):
            raise ScientificAnnotationError(
                f"annotations.proteins record {index} requires at least one of "
                "uniprot_id, gene_symbol, or display_name"
            )
        connection.execute(
            """INSERT INTO protein_identities
            (protein_key, uniprot_id, gene_symbol, display_name, source_path,
             source_sha256, source_record_index, source_record_json, provenance_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                protein_key,
                identity_values[0] or None,
                identity_values[1] or None,
                identity_values[2] or None,
                str(source_path),
                source_sha256,
                index,
                _json(record),
                _provenance(record),
            ),
        )
    return len(prepared)


def _ingest_receptors(
    connection: sqlite3.Connection,
    records: list[dict[str, Any]],
    source_path: Path,
    source_sha256: str,
) -> int:
    prepared = _unique_records(records, _context_key, "receptors")
    fields = (
        "protein_key",
        "receptor_classification",
        "qualification_status",
        "qualification_reason",
        "native_redock_status",
        "native_redock_reason",
        "native_redock_rmsd",
    )
    for key, index, record in prepared:
        label = f"annotations.receptors record {index}"
        receptor_classification, evidence_fields, receptor_audits = (
            _receptor_evidence_contract(record, label)
        )
        has_legacy_field = any(
            record.get(field) not in (None, "") for field in fields
        )
        if not has_legacy_field and not any(evidence_fields.values()) and not receptor_audits:
            raise ScientificAnnotationError(
                f"{label} requires at least one explicit annotation field"
            )
        context_id = _context_id(connection, key, label)
        protein_key = _text(record.get("protein_key"))
        protein_identity_id: int | None = None
        if protein_key:
            protein_row = connection.execute(
                "SELECT protein_identity_id FROM protein_identities WHERE protein_key=?",
                (protein_key,),
            ).fetchone()
            if protein_row is None:
                raise ScientificAnnotationError(
                    f"{label} references unknown protein_key {protein_key!r}"
                )
            protein_identity_id = int(protein_row[0])
        native_redock_rmsd = _optional_float(
            record.get("native_redock_rmsd"), "native_redock_rmsd", label
        )
        if native_redock_rmsd is not None and native_redock_rmsd < 0:
            raise ScientificAnnotationError(
                f"{label} requires non-negative native_redock_rmsd"
            )
        connection.execute(
            """INSERT INTO receptor_annotations
            (receptor_context_id, protein_identity_id, receptor_classification,
             classification_method, chemistry_evidence_status,
             prepared_receptor_sha256, canonical_chemistry_policy_sha256,
             retained_metal_atom_count, retained_cofactor_residue_count,
             requested_observed_conflict,
             qualification_status, qualification_reason, native_redock_status,
             native_redock_reason, native_redock_rmsd, source_path,
             source_sha256, source_record_index, source_record_json,
             provenance_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                context_id,
                protein_identity_id,
                receptor_classification,
                evidence_fields["classification_method"],
                evidence_fields["chemistry_evidence_status"],
                evidence_fields["prepared_receptor_sha256"],
                evidence_fields["canonical_chemistry_policy_sha256"],
                evidence_fields["retained_metal_atom_count"],
                evidence_fields["retained_cofactor_residue_count"],
                evidence_fields["requested_observed_conflict"],
                _text(record.get("qualification_status")) or None,
                _text(record.get("qualification_reason")) or None,
                _text(record.get("native_redock_status")) or None,
                _text(record.get("native_redock_reason")) or None,
                native_redock_rmsd,
                str(source_path),
                source_sha256,
                index,
                _json(record),
                _provenance(record),
            ),
        )
        for audit in receptor_audits:
            connection.execute(
                """INSERT INTO receptor_audits
                (receptor_context_id, audit_kind, source_path, source_sha256,
                 parse_status, audit_json)
                VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    context_id,
                    audit["audit_kind"],
                    audit["source_path"],
                    audit["source_sha256"],
                    audit["parse_status"],
                    audit["audit_json"],
                ),
            )
    return len(prepared)


def _ingest_known_pairs(
    connection: sqlite3.Connection,
    records: list[dict[str, Any]],
    source_path: Path,
    source_sha256: str,
) -> int:
    prepared = _unique_records(records, _context_key, "known_pairs")
    for key, index, record in prepared:
        label = f"annotations.known_pairs record {index}"
        context_id = _context_id(connection, key, label)
        ligand = _required(record, "ligand_canonical_id", label)
        pair_row = connection.execute(
            """SELECT p.pair_cell_id FROM pair_cells p
            JOIN ligands l ON l.ligand_id=p.ligand_id
            WHERE p.receptor_context_id=? AND l.canonical_id=?""",
            (context_id, ligand),
        ).fetchone()
        connection.execute(
            """INSERT INTO known_pair_selections
            (receptor_context_id, ligand_canonical_id, pair_cell_id,
             selection_label, evidence_reference, source_path, source_sha256,
             source_record_index, source_record_json, provenance_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                context_id,
                ligand,
                int(pair_row[0]) if pair_row is not None else None,
                _text(record.get("selection_label")) or None,
                _text(record.get("evidence_reference")) or None,
                str(source_path),
                source_sha256,
                index,
                _json(record),
                _provenance(record),
            ),
        )
    return len(prepared)


def ingest_scientific_annotations(
    connection: sqlite3.Connection,
    manifest: Mapping[str, Any],
    manifest_path: Path,
) -> dict[str, int]:
    """Ingest only manifest-declared annotation records; never infer values."""
    specs = manifest.get("annotations")
    counts = {name: 0 for name in ANNOTATION_SOURCE_NAMES}
    if specs in (None, {}):
        return counts
    if not isinstance(specs, Mapping):
        raise ScientificAnnotationError("annotations must be a mapping")
    loaded: dict[str, tuple[list[dict[str, Any]], Path, str]] = {}
    for name in ANNOTATION_SOURCE_NAMES:
        if name not in specs:
            continue
        path = _source_path(specs[name], manifest_path, name)
        loaded[name] = (_load_records(path, name), path, _sha256(path))

    if "proteins" in loaded:
        records, path, digest = loaded["proteins"]
        counts["proteins"] = _ingest_proteins(connection, records, path, digest)
    if "receptors" in loaded:
        records, path, digest = loaded["receptors"]
        counts["receptors"] = _ingest_receptors(connection, records, path, digest)
    if "known_pairs" in loaded:
        records, path, digest = loaded["known_pairs"]
        counts["known_pairs"] = _ingest_known_pairs(connection, records, path, digest)
    return counts


__all__ = [
    "ANNOTATION_SOURCE_NAMES",
    "EXPLICIT_QUALIFIED_STATUSES",
    "ScientificAnnotationError",
    "ingest_scientific_annotations",
    "is_explicitly_qualified",
]
