from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ARTIFACT_NAME = "consensus_reranked_scorch.csv"
_KNOWN_SUFFIX = re.compile(r"\.(?:pdbqt|mol2|sdf)(?:\.gz)?$", re.IGNORECASE)
_LEGACY_STAGE_SUFFIX = re.compile(
    r"(?:_(?:dud_)?(?:gnina|dock6|ledock)?_?stage\d+)+$", re.IGNORECASE
)
_DECOY_PREFIX = re.compile(r"^(?:dud_|decoys?_)", re.IGNORECASE)
_MISSING_TEXT = {"", "na", "n/a", "nan", "none", "null"}
_TRUTHY_TEXT = {"1", "true", "t", "yes", "y", "on"}

PROVENANCE_FIELDS: tuple[str, ...] = (
    "final_score_source",
    "z_vs_decoys_blend",
    "z_vs_decoys_consensus",
    "ml_blend_mode",
    "ml_blend_scorch_weight_effective",
    "ml_blend_cnn_weight_effective",
    "scorch_composite",
    "SCORCH_score_used",
    "SCORCH_certainty_used",
    "scorch_pct",
    "blend_mu_decoy",
    "blend_sigma_decoy",
    "blend_n_decoys",
    "consensus_mu_decoy",
    "consensus_sigma_decoy",
    "consensus_n_decoys",
    "selected_stage",
    "rescored_stage",
    "rescored_flag",
    "stage_match_flag",
    "stage_fallback_reason",
    "scorch_source_used",
    "selected_docking_score",
    "best_engine",
    "best_signal",
)


def _score_source_fields(output_score_column: str) -> tuple[str, str, str]:
    return (
        f"{output_score_column}_reference_source_file",
        f"{output_score_column}_reference_source_sha256",
        f"{output_score_column}_source_provenance",
    )


def _reference_fields(output_score_column: str) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            (
                output_score_column,
                *PROVENANCE_FIELDS,
                *_score_source_fields(output_score_column),
            )
        )
    )


REFERENCE_FIELDS: tuple[str, ...] = _reference_fields("final_score")

_NUMERIC_FIELDS = frozenset(
    {
        "final_score",
        "z_vs_decoys_blend",
        "z_vs_decoys_consensus",
        "ml_blend_scorch_weight_effective",
        "ml_blend_cnn_weight_effective",
        "scorch_composite",
        "SCORCH_score_used",
        "SCORCH_certainty_used",
        "scorch_pct",
        "blend_mu_decoy",
        "blend_sigma_decoy",
        "blend_n_decoys",
        "consensus_mu_decoy",
        "consensus_sigma_decoy",
        "consensus_n_decoys",
    }
)


@dataclass(frozen=True)
class ReferenceRecord:
    key: tuple[str, str]
    values: Mapping[str, str]
    source_file: str
    source_sha256: str


@dataclass
class SourceInventory:
    source_file: str
    source_sha256: str
    pdb_id: str
    rows_total: int = 0
    rows_scored: int = 0
    rows_without_reference_score: int = 0
    rows_decoy: int = 0
    rows_control: int = 0
    rows_missing_ligand: int = 0
    rows_invalid_numeric_provenance: int = 0
    inferred_score_sources: int = 0


@dataclass
class ReferenceLoadResult:
    records: dict[tuple[str, str], ReferenceRecord] = field(default_factory=dict)
    inventory: list[SourceInventory] = field(default_factory=list)
    conflicts: list[dict[str, Any]] = field(default_factory=list)
    duplicate_identical: int = 0
    discovered_files: int = 0


def _is_missing(value: Any) -> bool:
    return str(value or "").strip().lower() in _MISSING_TEXT


def _truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in _TRUTHY_TEXT


def _finite_text(value: Any) -> str:
    text = str(value or "").strip()
    if _is_missing(text):
        return ""
    try:
        number = float(text)
    except (TypeError, ValueError):
        return ""
    return text if math.isfinite(number) else ""


def normalize_pdb_id(value: Any) -> str:
    return str(value or "").strip().upper()


def normalize_ligand_base(value: Any) -> str:
    text = str(value or "").strip().replace("\\", "/")
    text = text.rsplit("/", 1)[-1]
    text = text.replace(".sanitized", "")
    text = _KNOWN_SUFFIX.sub("", text)
    text = _LEGACY_STAGE_SUFFIX.sub("", text)
    text = re.sub(r"__+", "_", text).rstrip("_")
    return text.casefold()


def _row_ligand(row: Mapping[str, Any]) -> str:
    for key in (
        "ligand_base",
        "ligand",
        "ligand_file",
        "Ligand_ID",
        "ligand_id",
        "name",
    ):
        if not _is_missing(row.get(key)):
            return str(row.get(key) or "")
    return ""


def _row_is_control(row: Mapping[str, Any]) -> bool:
    if _truthy(row.get("is_control")):
        return True
    role = str(row.get("row_role") or row.get("role") or "").strip().lower()
    return role in {"control", "positive_control", "negative_control"}


def _row_is_decoy(row: Mapping[str, Any], ligand: str) -> bool:
    explicit = str(row.get("is_decoy") or "").strip()
    if explicit:
        return _truthy(explicit)
    identifier = str(ligand or "").strip().replace("\\", "/").rsplit("/", 1)[-1]
    return bool(_DECOY_PREFIX.search(identifier))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f"CSV has no header: {path}")
        duplicates = sorted(
            name for name, count in Counter(fields).items() if count > 1
        )
        if duplicates:
            raise ValueError(f"CSV has duplicate columns {duplicates}: {path}")
        return fields, [dict(row) for row in reader]


def _artifact_pdb_id(root: Path, path: Path) -> str:
    relative = path.relative_to(root)
    if len(relative.parts) < 2:
        raise ValueError(f"score artifact is not under a PDB directory: {path}")
    pdb_id = normalize_pdb_id(relative.parts[0])
    if not pdb_id:
        raise ValueError(f"could not derive PDB ID from score artifact: {path}")
    return pdb_id


def _relative_source(root: Path, path: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _values_equal(
    left: str,
    right: str,
    field_name: str,
    *,
    numeric_score: bool = False,
) -> bool:
    if _is_missing(left) and _is_missing(right):
        return True
    if numeric_score or field_name in _NUMERIC_FIELDS:
        try:
            return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=1e-12)
        except (TypeError, ValueError):
            return False
    return str(left).strip() == str(right).strip()


def _infer_final_score_source(
    row: Mapping[str, Any], final_score: str
) -> tuple[str, str]:
    explicit = str(row.get("final_score_source") or "").strip()
    if explicit:
        return explicit, "artifact_explicit"
    matches = [
        field_name
        for field_name in (
            "z_vs_decoys_blend",
            "z_vs_decoys_consensus",
            "scorch_composite",
        )
        if _finite_text(row.get(field_name))
        and _values_equal(final_score, str(row.get(field_name) or ""), field_name)
    ]
    if len(matches) == 1:
        return f"legacy_exact_match_{matches[0]}", "inferred_exact_value_match"
    if matches:
        joined = "+".join(matches)
        return f"legacy_ambiguous_exact_match_{joined}", "inferred_ambiguous_match"
    return "legacy_unclassified_final_score", "artifact_missing_no_exact_match"


def _reference_values(
    row: Mapping[str, Any],
    *,
    reference_score_column: str,
    output_score_column: str,
    source_file: str,
    source_sha256: str,
    inventory: SourceInventory,
) -> dict[str, str] | None:
    selected_score = _finite_text(row.get(reference_score_column))
    if not selected_score:
        inventory.rows_without_reference_score += 1
        return None
    values: dict[str, str] = {output_score_column: selected_score}
    for field_name in PROVENANCE_FIELDS:
        if field_name == output_score_column:
            continue
        raw = str(row.get(field_name) or "").strip()
        if _is_missing(raw):
            values[field_name] = ""
        elif field_name in _NUMERIC_FIELDS:
            values[field_name] = _finite_text(raw)
            if not values[field_name]:
                inventory.rows_invalid_numeric_provenance += 1
        else:
            values[field_name] = raw

    source_file_field, source_hash_field, source_provenance_field = (
        _score_source_fields(output_score_column)
    )
    if output_score_column == "final_score":
        source, provenance = _infer_final_score_source(row, selected_score)
        if provenance != "artifact_explicit":
            inventory.inferred_score_sources += 1
        values["final_score_source"] = source
        values[source_provenance_field] = provenance
    else:
        values[source_provenance_field] = f"artifact_column:{reference_score_column}"
    values[source_file_field] = source_file
    values[source_hash_field] = source_sha256
    return values


def _record_conflicts(
    existing: ReferenceRecord,
    candidate: ReferenceRecord,
    *,
    comparison_fields: Sequence[str],
    output_score_column: str,
) -> list[dict[str, Any]]:
    conflicts: list[dict[str, Any]] = []
    for field_name in comparison_fields:
        left = str(existing.values.get(field_name) or "")
        right = str(candidate.values.get(field_name) or "")
        if _is_missing(left) or _is_missing(right):
            continue
        if not _values_equal(
            left,
            right,
            field_name,
            numeric_score=field_name == output_score_column,
        ):
            conflicts.append(
                {
                    "pdb_id": existing.key[0],
                    "ligand_base": existing.key[1],
                    "field": field_name,
                    "existing_value": left,
                    "candidate_value": right,
                    "existing_source_file": existing.source_file,
                    "candidate_source_file": candidate.source_file,
                    "existing_source_sha256": existing.source_sha256,
                    "candidate_source_sha256": candidate.source_sha256,
                }
            )
    return conflicts


def load_reference_scores(
    reference_root: str | Path,
    *,
    reference_score_column: str = "final_score",
    output_score_column: str | None = None,
) -> ReferenceLoadResult:
    selected_reference_column = str(reference_score_column or "").strip()
    selected_output_column = str(
        output_score_column or selected_reference_column
    ).strip()
    if not selected_reference_column:
        raise ValueError("reference_score_column cannot be blank")
    if not selected_output_column:
        raise ValueError("output_score_column cannot be blank")
    if selected_output_column in {"pdb_id", "ligand_base"}:
        raise ValueError(
            "output_score_column cannot replace a strict join column: "
            f"{selected_output_column}"
        )

    root = Path(reference_root).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"post-docked reference root does not exist: {root}")
    artifact_paths = sorted(
        path for path in root.rglob(ARTIFACT_NAME) if path.is_file()
    )
    if not artifact_paths:
        raise ValueError(f"no per-PDB {ARTIFACT_NAME} artifacts found under {root}")

    result = ReferenceLoadResult(discovered_files=len(artifact_paths))
    comparison_fields = tuple(
        dict.fromkeys((selected_output_column, *PROVENANCE_FIELDS))
    )
    for artifact_path in artifact_paths:
        path_pdb = _artifact_pdb_id(root, artifact_path)
        source_file = _relative_source(root, artifact_path)
        source_sha256 = _sha256(artifact_path)
        inventory = SourceInventory(source_file, source_sha256, path_pdb)
        fields, rows = _read_rows(artifact_path)
        if selected_reference_column not in fields:
            raise ValueError(
                "reference artifact lacks requested score column "
                f"{selected_reference_column!r}: {artifact_path}"
            )
        for row in rows:
            inventory.rows_total += 1
            row_pdb = normalize_pdb_id(row.get("pdb_id"))
            if row_pdb and row_pdb != path_pdb:
                raise ValueError(
                    f"PDB mismatch in {artifact_path}: directory={path_pdb!r}, row={row_pdb!r}"
                )
            ligand = _row_ligand(row)
            if _row_is_control(row):
                inventory.rows_control += 1
                continue
            if _row_is_decoy(row, ligand):
                inventory.rows_decoy += 1
                continue
            ligand_base = normalize_ligand_base(ligand)
            if not ligand_base:
                inventory.rows_missing_ligand += 1
                continue
            values = _reference_values(
                row,
                reference_score_column=selected_reference_column,
                output_score_column=selected_output_column,
                source_file=source_file,
                source_sha256=source_sha256,
                inventory=inventory,
            )
            if values is None:
                continue
            inventory.rows_scored += 1
            key = (path_pdb, ligand_base)
            candidate = ReferenceRecord(
                key=key,
                values=values,
                source_file=source_file,
                source_sha256=source_sha256,
            )
            existing = result.records.get(key)
            if existing is None:
                result.records[key] = candidate
                continue
            conflicts = _record_conflicts(
                existing,
                candidate,
                comparison_fields=comparison_fields,
                output_score_column=selected_output_column,
            )
            if conflicts:
                result.conflicts.extend(conflicts)
            else:
                result.duplicate_identical += 1
        result.inventory.append(inventory)
    return result


def _nonmissing_cell_hash(
    rows: Sequence[Mapping[str, str]],
    original_fields: Sequence[str],
    positions: Iterable[tuple[int, str]],
) -> str:
    digest = hashlib.sha256()
    field_set = set(original_fields)
    for row_index, field_name in positions:
        if field_name not in field_set:
            raise AssertionError(f"nonmissing hash field disappeared: {field_name}")
        payload = [row_index, field_name, str(rows[row_index].get(field_name) or "")]
        digest.update(json.dumps(payload, ensure_ascii=True).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _write_csv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=True, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _inventory_rows(
    inventory: Sequence[SourceInventory],
    *,
    reference_score_column: str,
    output_score_column: str,
) -> list[dict[str, Any]]:
    return [
        {
            "source_file": item.source_file,
            "source_sha256": item.source_sha256,
            "pdb_id": item.pdb_id,
            "reference_score_column": reference_score_column,
            "output_score_column": output_score_column,
            "rows_total": item.rows_total,
            "rows_scored": item.rows_scored,
            "rows_without_reference_score": item.rows_without_reference_score,
            "rows_decoy": item.rows_decoy,
            "rows_control": item.rows_control,
            "rows_missing_ligand": item.rows_missing_ligand,
            "rows_invalid_numeric_provenance": item.rows_invalid_numeric_provenance,
            "inferred_score_sources": item.inferred_score_sources,
        }
        for item in inventory
    ]


def _source_inventory_hash(inventory: Sequence[SourceInventory]) -> str:
    return _json_sha256(
        [
            {"source_file": item.source_file, "source_sha256": item.source_sha256}
            for item in inventory
        ]
    )


def backfill_post_docked_final_scores(
    dataset: str | Path,
    reference_root: str | Path,
    out: str | Path,
    audit_dir: str | Path,
    *,
    overwrite_existing: bool = False,
    reference_score_column: str = "final_score",
    output_score_column: str | None = None,
) -> dict[str, Any]:
    selected_reference_column = str(reference_score_column or "").strip()
    selected_output_column = str(
        output_score_column or selected_reference_column
    ).strip()
    dataset_path = Path(dataset).expanduser().resolve()
    reference_path = Path(reference_root).expanduser().resolve()
    out_path = Path(out).expanduser().resolve()
    audit_path = Path(audit_dir).expanduser().resolve()
    if not dataset_path.is_file():
        raise FileNotFoundError(f"dataset does not exist: {dataset_path}")
    if out_path == dataset_path:
        raise ValueError("--out must differ from --dataset")

    input_fields, rows = _read_rows(dataset_path)
    required = {"pdb_id", "ligand_base"}
    missing_required = sorted(required.difference(input_fields))
    if missing_required:
        raise ValueError(f"dataset missing strict join columns: {missing_required}")

    reference_fields = _reference_fields(selected_output_column)
    reference = load_reference_scores(
        reference_path,
        reference_score_column=selected_reference_column,
        output_score_column=selected_output_column,
    )
    audit_path.mkdir(parents=True, exist_ok=True)
    inventory_fields = (
        "source_file",
        "source_sha256",
        "pdb_id",
        "reference_score_column",
        "output_score_column",
        "rows_total",
        "rows_scored",
        "rows_without_reference_score",
        "rows_decoy",
        "rows_control",
        "rows_missing_ligand",
        "rows_invalid_numeric_provenance",
        "inferred_score_sources",
    )
    _write_csv(
        audit_path / "source_inventory.csv",
        inventory_fields,
        _inventory_rows(
            reference.inventory,
            reference_score_column=selected_reference_column,
            output_score_column=selected_output_column,
        ),
    )
    conflict_fields = (
        "pdb_id",
        "ligand_base",
        "field",
        "existing_value",
        "candidate_value",
        "existing_source_file",
        "candidate_source_file",
        "existing_source_sha256",
        "candidate_source_sha256",
    )
    _write_csv(
        audit_path / "reference_conflicts.csv", conflict_fields, reference.conflicts
    )
    if reference.conflicts:
        failure = {
            "status": "blocked_reference_conflicts",
            "reference_conflict_count": len(reference.conflicts),
            "reference_root": str(reference_path),
            "reference_score_column": selected_reference_column,
            "output_score_column": selected_output_column,
            "audit_dir": str(audit_path),
        }
        _write_json(audit_path / "backfill_manifest.json", failure)
        raise ValueError(
            f"reference contains {len(reference.conflicts)} duplicate-key score/provenance conflicts; "
            f"see {audit_path / 'reference_conflicts.csv'}"
        )
    if not reference.records:
        failure = {
            "status": "blocked_no_usable_reference_scores",
            "reference_artifact_count": reference.discovered_files,
            "reference_root": str(reference_path),
            "reference_score_column": selected_reference_column,
            "output_score_column": selected_output_column,
            "score_fallback_allowed": False,
            "consensus_fallback_allowed": False,
            "reference_inventory_sha256": _source_inventory_hash(reference.inventory),
            "audit_dir": str(audit_path),
        }
        _write_json(audit_path / "backfill_manifest.json", failure)
        raise ValueError(
            "reference artifacts contain no finite non-decoy values for "
            f"{selected_reference_column!r}; no alternate score will be substituted. See "
            f"{audit_path / 'source_inventory.csv'}"
        )

    output_fields = list(input_fields)
    for field_name in reference_fields:
        if field_name not in output_fields:
            output_fields.append(field_name)

    original_nonmissing_positions = [
        (row_index, field_name)
        for row_index, row in enumerate(rows)
        for field_name in input_fields
        if not _is_missing(row.get(field_name))
    ]
    before_nonmissing_hash = _nonmissing_cell_hash(
        rows, input_fields, original_nonmissing_positions
    )

    matched_keys: set[tuple[str, str]] = set()
    dataset_keys: set[tuple[str, str]] = set()
    missing_join_rows: list[dict[str, Any]] = []
    existing_value_conflicts: list[dict[str, Any]] = []
    fill_counts: Counter[str] = Counter()
    overwrite_counts: Counter[str] = Counter()
    matched_rows = 0
    scored_rows = 0
    for row_index, row in enumerate(rows):
        pdb_id = normalize_pdb_id(row.get("pdb_id"))
        ligand_base = normalize_ligand_base(row.get("ligand_base"))
        if not pdb_id or not ligand_base:
            missing_join_rows.append(
                {
                    "row_index": row_index,
                    "pdb_id": row.get("pdb_id", ""),
                    "ligand_base": row.get("ligand_base", ""),
                    "reason": "missing_normalized_join_key",
                }
            )
            continue
        key = (pdb_id, ligand_base)
        dataset_keys.add(key)
        record = reference.records.get(key)
        if record is None:
            continue
        matched_rows += 1
        matched_keys.add(key)
        if not _is_missing(record.values.get(selected_output_column)):
            scored_rows += 1
        for field_name in reference_fields:
            source_value = str(record.values.get(field_name) or "")
            if _is_missing(source_value):
                continue
            current_value = str(row.get(field_name) or "")
            if _is_missing(current_value):
                row[field_name] = source_value
                fill_counts[field_name] += 1
                continue
            if _values_equal(
                current_value,
                source_value,
                field_name,
                numeric_score=field_name == selected_output_column,
            ):
                continue
            existing_value_conflicts.append(
                {
                    "row_index": row_index,
                    "pdb_id": pdb_id,
                    "ligand_base": ligand_base,
                    "field": field_name,
                    "dataset_value": current_value,
                    "reference_value": source_value,
                    "action": "overwritten" if overwrite_existing else "preserved",
                    "reference_source_file": record.source_file,
                    "reference_source_sha256": record.source_sha256,
                }
            )
            if overwrite_existing:
                row[field_name] = source_value
                overwrite_counts[field_name] += 1

    after_nonmissing_hash = _nonmissing_cell_hash(
        rows, input_fields, original_nonmissing_positions
    )
    nonmissing_preserved = before_nonmissing_hash == after_nonmissing_hash
    if not overwrite_existing and not nonmissing_preserved:
        raise AssertionError(
            "fill-only backfill changed a pre-existing nonmissing cell"
        )

    unmatched_dataset_keys = sorted(dataset_keys.difference(reference.records))
    unmatched_reference_keys = sorted(set(reference.records).difference(dataset_keys))
    _write_csv(
        audit_path / "missing_join_key_rows.csv",
        ("row_index", "pdb_id", "ligand_base", "reason"),
        missing_join_rows,
    )
    _write_csv(
        audit_path / "unmatched_dataset_keys.csv",
        ("pdb_id", "ligand_base"),
        ({"pdb_id": key[0], "ligand_base": key[1]} for key in unmatched_dataset_keys),
    )
    _write_csv(
        audit_path / "unmatched_reference_keys.csv",
        ("pdb_id", "ligand_base"),
        ({"pdb_id": key[0], "ligand_base": key[1]} for key in unmatched_reference_keys),
    )
    _write_csv(
        audit_path / "existing_value_conflicts.csv",
        (
            "row_index",
            "pdb_id",
            "ligand_base",
            "field",
            "dataset_value",
            "reference_value",
            "action",
            "reference_source_file",
            "reference_source_sha256",
        ),
        existing_value_conflicts,
    )
    column_rows = [
        {
            "field": field_name,
            "filled_missing_cells": fill_counts[field_name],
            "overwritten_cells": overwrite_counts[field_name],
        }
        for field_name in reference_fields
    ]
    _write_csv(
        audit_path / "column_fill_summary.csv",
        ("field", "filled_missing_cells", "overwritten_cells"),
        column_rows,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(out_path, output_fields, rows)
    output_fields_check, output_rows_check = _read_rows(out_path)
    row_count_preserved = len(output_rows_check) == len(rows)
    row_order_hash_before = _json_sha256(
        [[row.get("pdb_id", ""), row.get("ligand_base", "")] for row in rows]
    )
    row_order_hash_after = _json_sha256(
        [
            [row.get("pdb_id", ""), row.get("ligand_base", "")]
            for row in output_rows_check
        ]
    )
    row_order_preserved = row_order_hash_before == row_order_hash_after
    schema_prefix_preserved = output_fields_check[: len(input_fields)] == input_fields
    if not (row_count_preserved and row_order_preserved and schema_prefix_preserved):
        raise AssertionError("output failed row-count/order or schema-prefix invariant")

    coverage_row: dict[str, Any] = {
        "dataset_rows": len(rows),
        "dataset_unique_keys": len(dataset_keys),
        "reference_scored_keys": len(reference.records),
        "matched_rows": matched_rows,
        "matched_unique_keys": len(matched_keys),
        "rows_with_reference_score": scored_rows,
        "missing_join_key_rows": len(missing_join_rows),
        "unmatched_dataset_unique_keys": len(unmatched_dataset_keys),
        "unmatched_reference_unique_keys": len(unmatched_reference_keys),
        "filled_cells": sum(fill_counts.values()),
        "overwritten_cells": sum(overwrite_counts.values()),
        "existing_value_conflicts": len(existing_value_conflicts),
    }
    if selected_reference_column == "final_score":
        coverage_row["rows_with_reference_final_score"] = scored_rows
    _write_csv(audit_path / "coverage.csv", tuple(coverage_row), [coverage_row])

    manifest: dict[str, Any] = {
        "status": "ok",
        "mode": "overwrite_existing" if overwrite_existing else "fill_only",
        "dataset": str(dataset_path),
        "reference_root": str(reference_path),
        "reference_score_column": selected_reference_column,
        "output_score_column": selected_output_column,
        "output": str(out_path),
        "audit_dir": str(audit_path),
        "artifact_name": ARTIFACT_NAME,
        "score_fallback_allowed": False,
        "consensus_fallback_allowed": False,
        "reference_artifact_count": reference.discovered_files,
        "reference_duplicate_identical_rows": reference.duplicate_identical,
        **coverage_row,
        "row_count_preserved": row_count_preserved,
        "row_order_preserved": row_order_preserved,
        "schema_prefix_preserved": schema_prefix_preserved,
        "preexisting_nonmissing_cells_preserved": nonmissing_preserved,
        "input_sha256": _sha256(dataset_path),
        "reference_inventory_sha256": _source_inventory_hash(reference.inventory),
        "output_sha256": _sha256(out_path),
        "preexisting_nonmissing_cells_sha256_before": before_nonmissing_hash,
        "preexisting_nonmissing_cells_sha256_after": after_nonmissing_hash,
        "row_order_sha256_before": row_order_hash_before,
        "row_order_sha256_after": row_order_hash_after,
    }
    _write_json(audit_path / "backfill_manifest.json", manifest)
    return manifest


__all__ = [
    "ARTIFACT_NAME",
    "PROVENANCE_FIELDS",
    "REFERENCE_FIELDS",
    "backfill_post_docked_final_scores",
    "load_reference_scores",
    "normalize_ligand_base",
    "normalize_pdb_id",
]
