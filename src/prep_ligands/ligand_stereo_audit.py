"""Deterministic, non-mutating stereochemistry evidence audits.

The audit reports stereochemistry encoded by source SDF or delimited
structure records.  It never assigns missing stereochemistry and deliberately
does not accept PDBQT as a stereochemistry source.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from rdkit import Chem, rdBase
from rdkit.Chem import inchi


SCHEMA_VERSION = "atlas-ligand-stereo-audit-v1"
METHOD = "rdkit.FindPotentialStereo(cleanIt=True,flagPossible=True)"

PARSED = "parsed"
PARSE_FAILED = "parse_failed"

NO_POTENTIAL_STEREO = "no_potential_stereo"
FULLY_SPECIFIED = "fully_specified"
PARTIALLY_UNSPECIFIED = "partially_unspecified"
ALL_UNSPECIFIED = "all_unspecified"

_SDF_DELIMITER_RE = re.compile(rb"(?m)^\$\$\$\$[ \t]*(?:\r?\n|\Z)")
# SDF data headers may include a record annotation after the field name.
_SDF_PROPERTY_RE = re.compile(r"^>\s*<([^>]+)>.*$")
_NO_ATOM = (1 << 32) - 1

_CSV_FIELDS = (
    "schema_version",
    "library_id",
    "ligand_canonical_id",
    "parent_ligand_id",
    "source_kind",
    "source_path",
    "source_file_sha256",
    "source_record_index",
    "source_record_id",
    "source_record_sha256",
    "audited_structure_sha256",
    "structure_column",
    "id_column",
    "structure_format",
    "prepared_artifact_path",
    "prepared_artifact_declared_sha256",
    "prepared_artifact_computed_sha256",
    "prepared_artifact_verification_status",
    "source_record_evidence_json",
    "parse_status",
    "audit_status",
    "potential_stereo_count",
    "specified_stereo_count",
    "unspecified_stereo_count",
    "unknown_stereo_count",
    "tetrahedral_count",
    "double_bond_count",
    "other_stereo_count",
    "has_unspecified_potential_stereo",
    "has_unresolved_potential_stereo",
    "canonical_isomeric_smiles",
    "canonical_smiles",
    "inchikey",
    "stereo_elements_json",
    "rdkit_version",
    "method",
    "error",
)


@dataclass(frozen=True)
class LigandStereoSource:
    """One source input and its structure/identifier fields."""

    path: Path
    kind: Literal["sdf", "csv"]
    id_field: str
    parent_id_field: str | None = None
    structure_field: str | None = None
    structure_format: Literal["sdf", "smiles", "inchi", "molblock"] = "sdf"
    prepared_path_field: str | None = None
    prepared_sha256_field: str | None = None
    prepared_path_root: Path | None = None


@dataclass(frozen=True)
class LigandStereoAuditConfig:
    library_id: str
    output_dir: Path
    sources: tuple[LigandStereoSource, ...]


@dataclass(frozen=True)
class LigandStereoAuditResult:
    records_csv: Path
    evidence_jsonl: Path
    summary_json: Path
    total_records: int
    parsed_records: int
    parse_failed_records: int
    records_with_unspecified_stereo: int
    records_with_unresolved_stereo: int


@dataclass(frozen=True)
class _SourceRecord:
    source: LigandStereoSource
    source_path: Path
    source_file_sha256: str
    source_record_index: int
    source_record_id: str
    source_record_sha256: str
    parent_ligand_id: str
    audited_structure_sha256: str
    structure_text: str
    prepared_artifact_path: str
    prepared_artifact_declared_sha256: str
    prepared_artifact_computed_sha256: str
    prepared_artifact_verification_status: str
    evidence: Mapping[str, Any]


@dataclass(frozen=True)
class _AuditRecord:
    schema_version: str
    library_id: str
    ligand_canonical_id: str
    source_kind: str
    parent_ligand_id: str
    source_path: str
    source_file_sha256: str
    source_record_index: int
    source_record_id: str
    source_record_sha256: str
    audited_structure_sha256: str
    structure_column: str
    id_column: str
    structure_format: str
    prepared_artifact_path: str
    prepared_artifact_declared_sha256: str
    prepared_artifact_computed_sha256: str
    prepared_artifact_verification_status: str
    source_record_evidence_json: str
    parse_status: str
    audit_status: str
    potential_stereo_count: int
    specified_stereo_count: int
    unspecified_stereo_count: int
    unknown_stereo_count: int
    tetrahedral_count: int
    double_bond_count: int
    other_stereo_count: int
    has_unspecified_potential_stereo: bool
    has_unresolved_potential_stereo: bool
    canonical_isomeric_smiles: str
    canonical_smiles: str
    inchikey: str
    stereo_elements_json: str
    rdkit_version: str
    method: str
    error: str


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def _enum_name(value: Any) -> str:
    name = getattr(value, "name", None)
    return str(name if name is not None else value)


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


def _iter_sdf_records(source: LigandStereoSource) -> list[_SourceRecord]:
    path = source.path.resolve()
    payload = path.read_bytes()
    file_sha256 = _sha256_bytes(payload)
    spans: list[bytes] = []
    start = 0
    for match in _SDF_DELIMITER_RE.finditer(payload):
        spans.append(payload[start : match.end()])
        start = match.end()
    if payload[start:].strip():
        spans.append(payload[start:])

    records: list[_SourceRecord] = []
    for index, record_bytes in enumerate(spans, start=1):
        record_text = record_bytes.decode("utf-8", errors="replace")
        structure_bytes = _sdf_structure_bytes(record_bytes)
        structure_text = structure_bytes.decode("utf-8", errors="replace")
        properties = _read_sdf_properties(record_text)
        header_name = (
            record_text.splitlines()[0].strip() if record_text.splitlines() else ""
        )
        if source.id_field == "_Name":
            source_id = header_name
        else:
            source_id = properties.get(source.id_field, "").strip()
        parent_ligand_id = (
            properties.get(source.parent_id_field, "").strip()
            if source.parent_id_field
            else ""
        )
        evidence = {
            "header_name": header_name,
            "id_field": source.id_field,
            "id_value": source_id,
            "properties": properties,
        }
        records.append(
            _SourceRecord(
                source=source,
                source_path=path,
                source_file_sha256=file_sha256,
                source_record_index=index,
                source_record_id=source_id,
                source_record_sha256=_sha256_bytes(record_bytes),
                parent_ligand_id=parent_ligand_id,
                audited_structure_sha256=_sha256_bytes(structure_bytes),
                structure_text=structure_text,
                prepared_artifact_path="",
                prepared_artifact_declared_sha256="",
                prepared_artifact_computed_sha256="",
                prepared_artifact_verification_status="not_declared",
                evidence=evidence,
            )
        )
    return records


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


def _validate_csv_fieldnames(
    source: LigandStereoSource,
    path: Path,
    fieldnames: Sequence[str] | None,
) -> Sequence[str]:
    if fieldnames is None:
        raise ValueError(f"CSV source has no header: {path}")
    if len(fieldnames) != len(set(fieldnames)):
        raise ValueError(f"CSV source has duplicate column names: {path}")
    structure_field = source.structure_field
    assert structure_field is not None
    declared_fields = [source.id_field, structure_field]
    optional_fields = (
        source.parent_id_field,
        source.prepared_path_field,
        source.prepared_sha256_field,
    )
    declared_fields.extend(field for field in optional_fields if field)
    missing = [field for field in declared_fields if field not in fieldnames]
    if missing:
        raise ValueError(f"CSV source {path} is missing columns: {', '.join(missing)}")
    return fieldnames


def _optional_csv_value(normalized: Mapping[str, Any], field: str | None) -> str:
    return str(normalized[field]).strip() if field else ""


def _csv_source_record(
    *,
    source: LigandStereoSource,
    path: Path,
    file_sha256: str,
    fieldnames: Sequence[str],
    index: int,
    row: Mapping[str | None, Any],
) -> _SourceRecord:
    normalized = _normalized_csv_row(fieldnames, row)
    canonical_row = _canonical_json(normalized)
    source_id = str(normalized[source.id_field]).strip()
    parent_ligand_id = _optional_csv_value(normalized, source.parent_id_field)
    assert source.structure_field is not None
    structure_text = str(normalized[source.structure_field])
    prepared_path_value = _optional_csv_value(normalized, source.prepared_path_field)
    prepared_declared_sha256 = _optional_csv_value(
        normalized, source.prepared_sha256_field
    ).lower()
    (
        prepared_artifact_path,
        prepared_computed_sha256,
        prepared_verification_status,
    ) = _prepared_artifact_evidence(
        source=source,
        source_path=path,
        declared_path=prepared_path_value,
        declared_sha256=prepared_declared_sha256,
    )
    evidence = {
        "id_column": source.id_field,
        "id_value": source_id,
        "parent_id_column": source.parent_id_field or "",
        "parent_id_value": parent_ligand_id,
        "record": normalized,
        "structure_column": source.structure_field,
        "prepared_artifact": {
            "computed_sha256": prepared_computed_sha256,
            "declared_path": prepared_path_value,
            "declared_sha256": prepared_declared_sha256,
            "path_column": source.prepared_path_field or "",
            "resolved_path": prepared_artifact_path,
            "sha256_column": source.prepared_sha256_field or "",
            "verification_status": prepared_verification_status,
        },
    }
    return _SourceRecord(
        source=source,
        source_path=path,
        source_file_sha256=file_sha256,
        parent_ligand_id=parent_ligand_id,
        source_record_index=index,
        source_record_id=source_id,
        source_record_sha256=_sha256_bytes(canonical_row.encode("utf-8")),
        audited_structure_sha256=_sha256_bytes(structure_text.encode("utf-8")),
        structure_text=structure_text,
        prepared_artifact_path=prepared_artifact_path,
        prepared_artifact_declared_sha256=prepared_declared_sha256,
        prepared_artifact_computed_sha256=prepared_computed_sha256,
        prepared_artifact_verification_status=prepared_verification_status,
        evidence=evidence,
    )


def _iter_csv_records(source: LigandStereoSource) -> list[_SourceRecord]:
    if source.structure_field is None:
        raise ValueError(f"CSV source {source.path} requires structure_field")
    path = source.path.resolve()
    payload = path.read_bytes()
    file_sha256 = _sha256_bytes(payload)
    text = payload.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text, newline=""))
    fieldnames = _validate_csv_fieldnames(source, path, reader.fieldnames)

    return [
        _csv_source_record(
            source=source,
            path=path,
            file_sha256=file_sha256,
            fieldnames=fieldnames,
            index=index,
            row=row,
        )
        for index, row in enumerate(reader, start=1)
    ]


def _prepared_artifact_evidence(
    *,
    source: LigandStereoSource,
    source_path: Path,
    declared_path: str,
    declared_sha256: str,
) -> tuple[str, str, str]:
    if source.prepared_path_field is None:
        return "", "", "not_declared"
    if not declared_path:
        return "", "", "missing_path_value"
    candidate = Path(declared_path).expanduser()
    if not candidate.is_absolute():
        root = source.prepared_path_root or source_path.parent
        candidate = root / candidate
    candidate = candidate.resolve()
    candidate_text = candidate.as_posix()
    if candidate.suffix.lower() != ".pdbqt":
        return candidate_text, "", "unsupported_artifact_type"
    if not candidate.is_file():
        return candidate_text, "", "missing_file"
    computed = _sha256_bytes(candidate.read_bytes())
    if not declared_sha256:
        return candidate_text, computed, "computed_only"
    if not re.fullmatch(r"[0-9a-f]{64}", declared_sha256):
        return candidate_text, computed, "invalid_declared_sha256"
    if computed == declared_sha256:
        return candidate_text, computed, "verified_match"
    return candidate_text, computed, "declared_mismatch"


def _parse_molecule(record: _SourceRecord) -> Chem.Mol | None:
    structure_format = record.source.structure_format
    if structure_format in {"sdf", "molblock"}:
        return Chem.MolFromMolBlock(
            record.structure_text,
            sanitize=True,
            removeHs=False,
            strictParsing=True,
        )
    if structure_format == "smiles":
        return Chem.MolFromSmiles(record.structure_text, sanitize=True)
    if structure_format == "inchi":
        return Chem.MolFromInchi(record.structure_text, sanitize=True, removeHs=False)
    raise ValueError(f"Unsupported structure format: {structure_format}")


def _stereo_elements(molecule: Chem.Mol) -> list[dict[str, Any]]:
    elements: list[dict[str, Any]] = []
    for info in Chem.FindPotentialStereo(
        molecule,
        cleanIt=True,
        flagPossible=True,
    ):
        controlling_atoms = [
            None if int(atom_index) == _NO_ATOM else int(atom_index)
            for atom_index in info.controllingAtoms
        ]
        elements.append(
            {
                "centered_on": int(info.centeredOn),
                "controlling_atoms": controlling_atoms,
                "descriptor": _enum_name(info.descriptor),
                "specified": _enum_name(info.specified),
                "type": _enum_name(info.type),
            }
        )
    return sorted(
        elements,
        key=lambda element: (
            element["type"],
            element["centered_on"],
            element["specified"],
            element["descriptor"],
        ),
    )


def _audit_status(potential: int, specified: int, unresolved: int) -> str:
    if potential == 0:
        return NO_POTENTIAL_STEREO
    if specified == 0 and unresolved == potential:
        return ALL_UNSPECIFIED
    if specified > 0 and unresolved > 0:
        return PARTIALLY_UNSPECIFIED
    if unresolved == 0:
        return FULLY_SPECIFIED
    return PARTIALLY_UNSPECIFIED


def _failed_record(
    config: LigandStereoAuditConfig, record: _SourceRecord, error: str
) -> _AuditRecord:
    ligand_id = record.source_record_id or f"__record_{record.source_record_index:06d}"
    return _AuditRecord(
        schema_version=SCHEMA_VERSION,
        library_id=config.library_id,
        ligand_canonical_id=ligand_id,
        parent_ligand_id=record.parent_ligand_id,
        source_kind=record.source.kind,
        source_path=record.source_path.as_posix(),
        source_file_sha256=record.source_file_sha256,
        source_record_index=record.source_record_index,
        source_record_id=record.source_record_id,
        source_record_sha256=record.source_record_sha256,
        audited_structure_sha256=record.audited_structure_sha256,
        structure_column=record.source.structure_field or "",
        id_column=record.source.id_field,
        structure_format=record.source.structure_format,
        prepared_artifact_path=record.prepared_artifact_path,
        prepared_artifact_declared_sha256=record.prepared_artifact_declared_sha256,
        prepared_artifact_computed_sha256=record.prepared_artifact_computed_sha256,
        prepared_artifact_verification_status=record.prepared_artifact_verification_status,
        source_record_evidence_json=_canonical_json(record.evidence),
        parse_status=PARSE_FAILED,
        audit_status=PARSE_FAILED,
        potential_stereo_count=0,
        specified_stereo_count=0,
        unspecified_stereo_count=0,
        unknown_stereo_count=0,
        tetrahedral_count=0,
        double_bond_count=0,
        other_stereo_count=0,
        has_unspecified_potential_stereo=False,
        has_unresolved_potential_stereo=False,
        canonical_isomeric_smiles="",
        canonical_smiles="",
        inchikey="",
        stereo_elements_json="[]",
        rdkit_version=rdBase.rdkitVersion,
        method=METHOD,
        error=error,
    )


def _audit_record(
    config: LigandStereoAuditConfig, record: _SourceRecord
) -> _AuditRecord:
    ligand_id = record.source_record_id or f"__record_{record.source_record_index:06d}"
    if not record.structure_text.strip():
        return _failed_record(config, record, "empty_structure")
    try:
        molecule = _parse_molecule(record)
    except Exception as exc:  # RDKit parsers can raise format-specific exceptions.
        return _failed_record(
            config, record, f"parse_exception:{type(exc).__name__}:{exc}"
        )
    if molecule is None:
        return _failed_record(config, record, "rdkit_parse_returned_none")

    elements = _stereo_elements(molecule)
    specified = sum(element["specified"] == "Specified" for element in elements)
    unspecified = sum(element["specified"] == "Unspecified" for element in elements)
    unknown = len(elements) - specified - unspecified
    tetrahedral = sum(element["type"] == "Atom_Tetrahedral" for element in elements)
    double_bond = sum(element["type"] == "Bond_Double" for element in elements)
    other = len(elements) - tetrahedral - double_bond

    return _AuditRecord(
        schema_version=SCHEMA_VERSION,
        library_id=config.library_id,
        ligand_canonical_id=ligand_id,
        parent_ligand_id=record.parent_ligand_id,
        source_kind=record.source.kind,
        source_path=record.source_path.as_posix(),
        source_file_sha256=record.source_file_sha256,
        source_record_index=record.source_record_index,
        source_record_id=record.source_record_id,
        source_record_sha256=record.source_record_sha256,
        audited_structure_sha256=record.audited_structure_sha256,
        structure_column=record.source.structure_field or "",
        id_column=record.source.id_field,
        structure_format=record.source.structure_format,
        prepared_artifact_path=record.prepared_artifact_path,
        prepared_artifact_declared_sha256=record.prepared_artifact_declared_sha256,
        prepared_artifact_computed_sha256=record.prepared_artifact_computed_sha256,
        prepared_artifact_verification_status=record.prepared_artifact_verification_status,
        source_record_evidence_json=_canonical_json(record.evidence),
        parse_status=PARSED,
        audit_status=_audit_status(len(elements), specified, unspecified + unknown),
        potential_stereo_count=len(elements),
        specified_stereo_count=specified,
        unspecified_stereo_count=unspecified,
        unknown_stereo_count=unknown,
        tetrahedral_count=tetrahedral,
        double_bond_count=double_bond,
        other_stereo_count=other,
        has_unspecified_potential_stereo=unspecified > 0,
        has_unresolved_potential_stereo=(unspecified + unknown) > 0,
        canonical_isomeric_smiles=Chem.MolToSmiles(
            molecule, canonical=True, isomericSmiles=True
        ),
        canonical_smiles=Chem.MolToSmiles(
            molecule, canonical=True, isomericSmiles=False
        ),
        inchikey=inchi.MolToInchiKey(molecule),
        stereo_elements_json=_canonical_json(elements),
        rdkit_version=rdBase.rdkitVersion,
        method=METHOD,
        error="",
    )


def _validate_sdf_source(source: LigandStereoSource) -> None:
    if source.structure_format != "sdf" or source.structure_field is not None:
        raise ValueError(
            "SDF sources require structure_format='sdf' and no structure_field"
        )
    if source.prepared_path_field or source.prepared_sha256_field:
        raise ValueError("Prepared-artifact columns are supported only for CSV sources")


def _validate_csv_source(source: LigandStereoSource) -> None:
    if source.structure_format not in {"smiles", "inchi", "molblock"}:
        raise ValueError("CSV structure_format must be smiles, inchi, or molblock")
    if not source.structure_field:
        raise ValueError("CSV sources require a structure_field")
    if source.prepared_sha256_field and not source.prepared_path_field:
        raise ValueError("CSV prepared_sha256_field requires prepared_path_field")


def _validate_source(
    source: LigandStereoSource,
    output_paths: set[Path],
) -> None:
    path = source.path.resolve()
    if path.suffix.lower() == ".pdbqt":
        raise ValueError("PDBQT is not accepted as source stereochemistry evidence")
    if not path.is_file():
        raise FileNotFoundError(path)
    if path in output_paths:
        raise ValueError(f"Output would overwrite source input: {path}")
    if source.kind == "sdf":
        _validate_sdf_source(source)
    elif source.kind == "csv":
        _validate_csv_source(source)
    else:
        raise ValueError(f"Unsupported source kind: {source.kind}")


def _validate_config(config: LigandStereoAuditConfig) -> None:
    if not config.library_id.strip():
        raise ValueError("library_id must not be empty")
    if not config.sources:
        raise ValueError("At least one SDF or CSV source is required")
    output_paths = {
        (config.output_dir / "ligand_stereo_audit.csv").resolve(),
        (config.output_dir / "ligand_stereo_evidence.jsonl").resolve(),
        (config.output_dir / "ligand_stereo_summary.json").resolve(),
    }
    for source in config.sources:
        _validate_source(source, output_paths)


def _write_csv(path: Path, records: Sequence[_AuditRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        for record in records:
            row = asdict(record)
            row["has_unspecified_potential_stereo"] = str(
                record.has_unspecified_potential_stereo
            ).lower()
            row["has_unresolved_potential_stereo"] = str(
                record.has_unresolved_potential_stereo
            ).lower()
            writer.writerow(row)


def _write_jsonl(path: Path, records: Sequence[_AuditRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(_canonical_json(asdict(record)))
            handle.write("\n")


def _collect_source_records(
    config: LigandStereoAuditConfig,
) -> tuple[list[_SourceRecord], list[dict[str, Any]]]:
    source_records: list[_SourceRecord] = []
    source_summaries: list[dict[str, Any]] = []
    sources = sorted(config.sources, key=lambda item: item.path.resolve().as_posix())
    for source in sources:
        records = (
            _iter_sdf_records(source)
            if source.kind == "sdf"
            else _iter_csv_records(source)
        )
        path = source.path.resolve()
        source_records.extend(records)
        source_summaries.append(
            {
                "kind": source.kind,
                "path": path.as_posix(),
                "record_count": len(records),
                "sha256": _sha256_bytes(path.read_bytes()),
            }
        )
    return source_records, source_summaries


def _validate_audited_record(record: _AuditRecord) -> None:
    if record.potential_stereo_count != (
        record.specified_stereo_count
        + record.unspecified_stereo_count
        + record.unknown_stereo_count
    ):
        raise AssertionError(
            "Stereo specification counts do not sum to potential count"
        )
    if record.potential_stereo_count != (
        record.tetrahedral_count + record.double_bond_count + record.other_stereo_count
    ):
        raise AssertionError("Stereo type counts do not sum to potential count")


def _audit_source_records(
    config: LigandStereoAuditConfig,
    source_records: Sequence[_SourceRecord],
) -> list[_AuditRecord]:
    audited = [_audit_record(config, record) for record in source_records]
    audited.sort(
        key=lambda record: (
            record.library_id,
            record.ligand_canonical_id,
            record.source_path,
            record.source_record_index,
            record.source_record_sha256,
        )
    )
    for record in audited:
        _validate_audited_record(record)
    return audited


def _write_audit_outputs(
    output_dir: Path,
    audited: Sequence[_AuditRecord],
) -> tuple[Path, Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records_csv = output_dir / "ligand_stereo_audit.csv"
    evidence_jsonl = output_dir / "ligand_stereo_evidence.jsonl"
    summary_json = output_dir / "ligand_stereo_summary.json"
    _write_csv(records_csv, audited)
    _write_jsonl(evidence_jsonl, audited)
    return records_csv, evidence_jsonl, summary_json


def _record_totals(audited: Sequence[_AuditRecord]) -> Counter[str]:
    totals: Counter[str] = Counter()
    for record in audited:
        totals["records"] += 1
        totals["potential"] += record.potential_stereo_count
        totals["specified"] += record.specified_stereo_count
        totals["unspecified"] += record.unspecified_stereo_count
        totals["unknown"] += record.unknown_stereo_count
        totals["records_unspecified"] += record.has_unspecified_potential_stereo
        totals["records_unresolved"] += record.has_unresolved_potential_stereo
    return totals


def _summary_payload(
    *,
    config: LigandStereoAuditConfig,
    audited: Sequence[_AuditRecord],
    source_summaries: Sequence[Mapping[str, Any]],
    records_csv: Path,
    evidence_jsonl: Path,
) -> tuple[dict[str, Any], Counter[str], Counter[str]]:
    status_counts = Counter(record.audit_status for record in audited)
    prepared_status_counts = Counter(
        record.prepared_artifact_verification_status for record in audited
    )
    duplicate_counts = Counter(record.ligand_canonical_id for record in audited)
    totals = _record_totals(audited)
    summary = {
        "audit_status_counts": dict(sorted(status_counts.items())),
        "duplicate_ligand_id_count": sum(
            count > 1 for count in duplicate_counts.values()
        ),
        "library_id": config.library_id,
        "method": METHOD,
        "outputs": {
            "evidence_jsonl": {
                "path": evidence_jsonl.name,
                "sha256": _sha256_bytes(evidence_jsonl.read_bytes()),
            },
            "records_csv": {
                "path": records_csv.name,
                "sha256": _sha256_bytes(records_csv.read_bytes()),
            },
        },
        "parse_failed_records": status_counts[PARSE_FAILED],
        "parsed_records": totals["records"] - status_counts[PARSE_FAILED],
        "prepared_artifact_verification_status_counts": dict(
            sorted(prepared_status_counts.items())
        ),
        "rdkit_version": rdBase.rdkitVersion,
        "records_with_unspecified_stereo": totals["records_unspecified"],
        "records_with_unresolved_stereo": totals["records_unresolved"],
        "schema_version": SCHEMA_VERSION,
        "sources": list(source_summaries),
        "total_potential_stereo_elements": totals["potential"],
        "total_specified_stereo_elements": totals["specified"],
        "total_unspecified_stereo_elements": totals["unspecified"],
        "total_unknown_stereo_elements": totals["unknown"],
        "total_records": totals["records"],
        "unique_ligand_id_count": len(duplicate_counts),
    }
    return summary, status_counts, totals


def _audit_result(
    *,
    records_csv: Path,
    evidence_jsonl: Path,
    summary_json: Path,
    status_counts: Mapping[str, int],
    totals: Mapping[str, int],
) -> LigandStereoAuditResult:
    parse_failed = status_counts[PARSE_FAILED]
    return LigandStereoAuditResult(
        records_csv=records_csv,
        evidence_jsonl=evidence_jsonl,
        summary_json=summary_json,
        total_records=totals["records"],
        parsed_records=totals["records"] - parse_failed,
        parse_failed_records=parse_failed,
        records_with_unspecified_stereo=totals["records_unspecified"],
        records_with_unresolved_stereo=totals["records_unresolved"],
    )


def audit_ligand_stereochemistry(
    config: LigandStereoAuditConfig,
) -> LigandStereoAuditResult:
    """Audit potential stereo encoded in source records without assigning it."""

    _validate_config(config)
    source_records, source_summaries = _collect_source_records(config)
    audited = _audit_source_records(config, source_records)
    records_csv, evidence_jsonl, summary_json = _write_audit_outputs(
        config.output_dir, audited
    )
    summary, status_counts, totals = _summary_payload(
        config=config,
        audited=audited,
        source_summaries=source_summaries,
        records_csv=records_csv,
        evidence_jsonl=evidence_jsonl,
    )
    summary_json.write_text(_canonical_json(summary) + "\n", encoding="utf-8")
    return _audit_result(
        records_csv=records_csv,
        evidence_jsonl=evidence_jsonl,
        summary_json=summary_json,
        status_counts=status_counts,
        totals=totals,
    )
