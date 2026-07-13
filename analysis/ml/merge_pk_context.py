from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import numbers
from pathlib import Path
from typing import Any, Sequence, cast

import pandas as pd


MISSING_TEXT = frozenset({"", "nan", "none", "null", "na", "n/a", "<na>"})
DIRECT_PK_COLUMNS = frozenset(
    {
        "cmax_um",
        "fraction_unbound_plasma",
        "free_cmax_um",
        "protein_binding_percent",
    }
)
DEFAULT_KEY_COLUMNS: tuple[tuple[str, ...], ...] = (
    ("ligand_base", "target_uniprot", "pdb_id"),
    ("ligand_base", "target_uniprot"),
    ("ligand_base", "pdb_id"),
    ("inchikey", "target_uniprot", "pdb_id"),
    ("inchikey", "target_uniprot"),
    ("inchikey", "pdb_id"),
    ("drug_id", "target_uniprot", "pdb_id"),
    ("drug_id", "target_uniprot"),
    ("drug_id", "pdb_id"),
    ("ligand_base", "target_id"),
    ("drug_id", "target_id"),
    ("generic_name", "target_uniprot"),
    ("generic_name", "pdb_id"),
    ("display_name", "target_uniprot"),
    ("display_name", "pdb_id"),
    ("ligand_base", "target_gene"),
    ("drug_id", "target_gene"),
)
IDENTITY_COLUMNS = (
    "ligand_base",
    "drug_id",
    "generic_name",
    "display_name",
    "inchikey",
    "target_uniprot",
    "pdb_id",
    "target_id",
    "target_gene",
)


@dataclass(frozen=True)
class KeySpec:
    name: str
    columns: tuple[str, ...]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_missing_value(value: Any) -> bool:
    if value is None:
        return True
    try:
        if bool(pd.isna(value)):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip().casefold() in MISSING_TEXT


def _missing_mask(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip().str.casefold()
    return (series.isna() | text.isin(MISSING_TEXT)).fillna(True)


def _normalize_key_component(value: Any) -> str:
    if _is_missing_value(value):
        return ""
    return " ".join(str(value).strip().casefold().split())


def _key_from_values(values: tuple[Any, ...]) -> tuple[str, ...] | None:
    parts = tuple(_normalize_key_component(value) for value in values)
    return parts if all(parts) else None


def _serialize_key(key: tuple[str, ...] | None) -> str:
    if key is None:
        return ""
    return json.dumps(list(key), ensure_ascii=True, separators=(",", ":"))


def _coerce_key_specs(
    key_columns: Sequence[Sequence[str]] | None,
) -> list[KeySpec]:
    raw_specs = key_columns if key_columns is not None else DEFAULT_KEY_COLUMNS
    specs: list[KeySpec] = []
    seen: set[tuple[str, ...]] = set()
    for raw in raw_specs:
        columns = tuple(str(column).strip() for column in raw if str(column).strip())
        if not columns:
            raise ValueError("Merge key specifications cannot be empty")
        if len(columns) != len(set(columns)):
            raise ValueError(f"Merge key contains duplicate columns: {columns}")
        if columns in seen:
            continue
        seen.add(columns)
        specs.append(KeySpec(name="+".join(columns), columns=columns))
    if not specs:
        raise ValueError("At least one merge key specification is required")
    return specs


def _frame_keys(frame: pd.DataFrame, spec: KeySpec) -> list[tuple[str, ...] | None]:
    return [
        _key_from_values(values)
        for values in frame.loc[:, list(spec.columns)].itertuples(
            index=False, name=None
        )
    ]


def _counter(keys: list[tuple[str, ...] | None]) -> Counter[tuple[str, ...]]:
    return Counter(key for key in keys if key is not None)


def _key_cardinality_rows(
    specs: list[KeySpec],
    primary_keys: dict[str, list[tuple[str, ...] | None]],
    backfill_keys: dict[str, list[tuple[str, ...] | None]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for priority, spec in enumerate(specs, start=1):
        primary_counts = _counter(primary_keys[spec.name])
        backfill_counts = _counter(backfill_keys[spec.name])
        shared = sorted(primary_counts.keys() & backfill_counts.keys())
        one_to_many = [
            key
            for key in shared
            if primary_counts[key] == 1 and backfill_counts[key] > 1
        ]
        many_to_many = [
            key
            for key in shared
            if primary_counts[key] > 1 and backfill_counts[key] > 1
        ]
        rows.append(
            {
                "priority": priority,
                "key_type": spec.name,
                "key_columns": ";".join(spec.columns),
                "primary_valid_rows": sum(primary_counts.values()),
                "primary_unique_keys": len(primary_counts),
                "primary_duplicate_keys": sum(
                    count > 1 for count in primary_counts.values()
                ),
                "primary_max_rows_per_key": max(primary_counts.values(), default=0),
                "backfill_valid_rows": sum(backfill_counts.values()),
                "backfill_unique_keys": len(backfill_counts),
                "backfill_duplicate_keys": sum(
                    count > 1 for count in backfill_counts.values()
                ),
                "backfill_max_rows_per_key": max(backfill_counts.values(), default=0),
                "shared_keys": len(shared),
                "one_to_one_keys": sum(
                    primary_counts[key] == 1 and backfill_counts[key] == 1
                    for key in shared
                ),
                "many_primary_to_one_backfill_keys": sum(
                    primary_counts[key] > 1 and backfill_counts[key] == 1
                    for key in shared
                ),
                "one_primary_to_many_backfill_keys": len(one_to_many),
                "many_to_many_keys": len(many_to_many),
                "ambiguous_key_examples": json.dumps(
                    [
                        _serialize_key(key)
                        for key in [*one_to_many[:3], *many_to_many[:3]][:5]
                    ]
                ),
                "selected_primary_rows": 0,
            }
        )
    return rows


def _detect_backfill_columns(
    backfill: pd.DataFrame,
    requested: Sequence[str] | None,
) -> list[str]:
    if requested is not None:
        columns = list(dict.fromkeys(str(column).strip() for column in requested))
        if any(not column for column in columns):
            raise ValueError("Backfill column names cannot be empty")
        if not columns:
            raise ValueError("At least one backfill column is required")
        missing = [column for column in columns if column not in backfill.columns]
        if missing:
            raise ValueError(
                "Requested backfill columns are absent from the PK-enriched table: "
                + ", ".join(missing)
            )
        return columns

    columns = []
    for column in backfill.columns:
        normalized = column.casefold()
        if normalized.startswith("pk_context_") or normalized in DIRECT_PK_COLUMNS:
            columns.append(column)
    if not columns:
        raise ValueError(
            "No PK/context columns were detected. Use explicit backfill columns if "
            "the source uses nonstandard names."
        )
    return columns


def _values_equal(left: Any, right: Any) -> bool:
    if _is_missing_value(left) or _is_missing_value(right):
        return _is_missing_value(left) and _is_missing_value(right)
    if (
        isinstance(left, numbers.Number)
        and isinstance(right, numbers.Number)
        and not isinstance(left, (bool,))
        and not isinstance(right, (bool,))
    ):
        return math.isclose(
            float(cast(Any, left)),
            float(cast(Any, right)),
            rel_tol=1e-12,
            abs_tol=0.0,
        )
    return str(left).strip() == str(right).strip()


def _canonical_value(value: Any) -> str:
    if isinstance(value, bool):
        return "bool:true" if value else "bool:false"
    if isinstance(value, numbers.Integral):
        return f"int:{int(value)}"
    if isinstance(value, numbers.Real):
        return f"float:{format(float(value), '.17g')}"
    return f"text:{str(value)}"


def _primary_nonmissing_hash(primary: pd.DataFrame, candidate: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for row_number in range(len(primary)):
        for column in primary.columns:
            original = primary.at[row_number, column]
            if _is_missing_value(original):
                continue
            payload = json.dumps(
                [row_number, column, _canonical_value(candidate.at[row_number, column])],
                ensure_ascii=True,
                separators=(",", ":"),
            )
            digest.update(payload.encode("utf-8"))
            digest.update(b"\n")
    return digest.hexdigest()


def _strongest_key_for_row(
    row_number: int,
    specs: list[KeySpec],
    keys: dict[str, list[tuple[str, ...] | None]],
) -> tuple[str, str]:
    for spec in specs:
        key = keys[spec.name][row_number]
        if key is not None:
            return spec.name, _serialize_key(key)
    return "", ""


def _identity_values(frame: pd.DataFrame, row_number: int) -> dict[str, Any]:
    return {
        column: frame.at[row_number, column]
        for column in IDENTITY_COLUMNS
        if column in frame.columns
    }


def _unmatched_primary_rows(
    primary: pd.DataFrame,
    matched_source: list[int | None],
    specs: list[KeySpec],
    primary_keys: dict[str, list[tuple[str, ...] | None]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row_number, source_row in enumerate(matched_source):
        if source_row is not None:
            continue
        key_type, key_value = _strongest_key_for_row(
            row_number, specs, primary_keys
        )
        rows.append(
            {
                "primary_row_number": row_number,
                "reason": "no_backfill_match" if key_type else "no_valid_key",
                "strongest_available_key_type": key_type,
                "strongest_available_key_value": key_value,
                **_identity_values(primary, row_number),
            }
        )
    columns = [
        "primary_row_number",
        "reason",
        "strongest_available_key_type",
        "strongest_available_key_value",
        *[column for column in IDENTITY_COLUMNS if column in primary.columns],
    ]
    return pd.DataFrame(rows, columns=columns)


def _unmatched_backfill_rows(
    backfill: pd.DataFrame,
    matched_source: list[int | None],
    specs: list[KeySpec],
    primary_keys: dict[str, list[tuple[str, ...] | None]],
    backfill_keys: dict[str, list[tuple[str, ...] | None]],
) -> pd.DataFrame:
    selected = {row for row in matched_source if row is not None}
    primary_key_sets = {
        spec.name: set(key for key in primary_keys[spec.name] if key is not None)
        for spec in specs
    }
    rows: list[dict[str, Any]] = []
    for row_number in range(len(backfill)):
        if row_number in selected:
            continue
        key_type, key_value = _strongest_key_for_row(
            row_number, specs, backfill_keys
        )
        overlaps_primary = any(
            backfill_keys[spec.name][row_number] in primary_key_sets[spec.name]
            for spec in specs
            if backfill_keys[spec.name][row_number] is not None
        )
        if not key_type:
            reason = "no_valid_key"
        elif overlaps_primary:
            reason = "not_selected_by_stronger_key"
        else:
            reason = "no_primary_match"
        rows.append(
            {
                "backfill_row_number": row_number,
                "reason": reason,
                "strongest_available_key_type": key_type,
                "strongest_available_key_value": key_value,
                **_identity_values(backfill, row_number),
            }
        )
    columns = [
        "backfill_row_number",
        "reason",
        "strongest_available_key_type",
        "strongest_available_key_value",
        *[column for column in IDENTITY_COLUMNS if column in backfill.columns],
    ]
    return pd.DataFrame(rows, columns=columns)


def merge_pk_context(
    primary_csv: str | Path,
    pk_enriched_csv: str | Path,
    output_csv: str | Path,
    *,
    audit_dir: str | Path | None = None,
    backfill_columns: Sequence[str] | None = None,
    key_columns: Sequence[Sequence[str]] | None = None,
    overwrite_conflicts: bool = False,
) -> dict[str, Any]:
    """Merge PK context into a primary table without changing its row contract."""
    primary_path = Path(primary_csv)
    backfill_path = Path(pk_enriched_csv)
    output_path = Path(output_csv)
    if not primary_path.is_file():
        raise FileNotFoundError(f"Primary CSV not found: {primary_path}")
    if not backfill_path.is_file():
        raise FileNotFoundError(f"PK-enriched CSV not found: {backfill_path}")
    if output_path.resolve() in {primary_path.resolve(), backfill_path.resolve()}:
        raise ValueError("Output CSV must not overwrite either input CSV")

    audit_path = (
        Path(audit_dir)
        if audit_dir is not None
        else output_path.parent / f"{output_path.stem}_audit"
    )
    primary = pd.read_csv(primary_path, low_memory=False).reset_index(drop=True)
    backfill = pd.read_csv(backfill_path, low_memory=False).reset_index(drop=True)
    if not primary.columns.is_unique:
        raise ValueError("Primary CSV has duplicate column names")
    if not backfill.columns.is_unique:
        raise ValueError("PK-enriched CSV has duplicate column names")

    requested_specs = _coerce_key_specs(key_columns)
    specs = [
        spec
        for spec in requested_specs
        if set(spec.columns).issubset(primary.columns)
        and set(spec.columns).issubset(backfill.columns)
    ]
    if not specs:
        requested = ", ".join(spec.name for spec in requested_specs)
        raise ValueError(
            "No merge key is available in both inputs. Requested key tiers: "
            + requested
        )

    primary_keys = {spec.name: _frame_keys(primary, spec) for spec in specs}
    backfill_keys = {spec.name: _frame_keys(backfill, spec) for spec in specs}
    primary_counts = {
        spec.name: _counter(primary_keys[spec.name]) for spec in specs
    }
    backfill_lookups: dict[
        str, dict[tuple[str, ...], list[int]]
    ] = {}
    for spec in specs:
        lookup: dict[tuple[str, ...], list[int]] = {}
        for row_number, key in enumerate(backfill_keys[spec.name]):
            if key is not None:
                lookup.setdefault(key, []).append(row_number)
        backfill_lookups[spec.name] = lookup
    key_audit_rows = _key_cardinality_rows(
        specs, primary_keys, backfill_keys
    )

    matched_source: list[int | None] = []
    match_key_types: list[str] = []
    match_key_values: list[str] = []
    ambiguous: list[dict[str, Any]] = []
    for row_number in range(len(primary)):
        selected: int | None = None
        selected_type = ""
        selected_value = ""
        for spec in specs:
            key = primary_keys[spec.name][row_number]
            if key is None:
                continue
            candidates = backfill_lookups[spec.name].get(key, [])
            if not candidates:
                continue
            if len(candidates) > 1:
                primary_cardinality = primary_counts[spec.name][key]
                ambiguous.append(
                    {
                        "primary_row_number": row_number,
                        "key_type": spec.name,
                        "key_value": _serialize_key(key),
                        "primary_rows_for_key": primary_cardinality,
                        "backfill_rows_for_key": len(candidates),
                        "cardinality": (
                            "many_to_many"
                            if primary_cardinality > 1
                            else "one_to_many"
                        ),
                        "backfill_row_numbers": candidates,
                    }
                )
                break
            selected = candidates[0]
            selected_type = spec.name
            selected_value = _serialize_key(key)
            break
        matched_source.append(selected)
        match_key_types.append(selected_type)
        match_key_values.append(selected_value)

    if ambiguous:
        examples = json.dumps(ambiguous[:5], sort_keys=True)
        raise ValueError(
            f"Rejected {len(ambiguous)} ambiguous source matches; no output was "
            f"written. Examples: {examples}"
        )

    selected_counts = Counter(match_key_types)
    for row in key_audit_rows:
        row["selected_primary_rows"] = selected_counts[row["key_type"]]

    columns = _detect_backfill_columns(backfill, backfill_columns)
    output = primary.copy()
    for column in columns:
        if column not in output.columns:
            output[column] = pd.Series(pd.NA, index=output.index, dtype="object")

    primary_value_hash_before = _primary_nonmissing_hash(primary, primary)
    conflicts: list[dict[str, Any]] = []
    missingness_rows: list[dict[str, Any]] = []
    for column in columns:
        existed = column in primary.columns
        before_values = (
            primary[column]
            if existed
            else pd.Series(pd.NA, index=primary.index, dtype="object")
        )
        before_missing = _missing_mask(before_values)
        source_nonmissing_matched = 0
        fill_rows: list[int] = []
        conflict_rows: list[int] = []
        for row_number, source_row in enumerate(matched_source):
            if source_row is None:
                continue
            source_value = backfill.at[source_row, column]
            if _is_missing_value(source_value):
                continue
            source_nonmissing_matched += 1
            primary_value = before_values.iat[row_number]
            if _is_missing_value(primary_value):
                fill_rows.append(row_number)
                continue
            if not _values_equal(primary_value, source_value):
                conflict_rows.append(row_number)
                conflicts.append(
                    {
                        "primary_row_number": row_number,
                        "backfill_row_number": source_row,
                        "match_key_type": match_key_types[row_number],
                        "match_key_value": match_key_values[row_number],
                        "column": column,
                        "primary_value": primary_value,
                        "backfill_value": source_value,
                        "action": (
                            "overwritten_opt_in"
                            if overwrite_conflicts
                            else "preserved_primary"
                        ),
                    }
                )
        for row_number in fill_rows:
            source_row = matched_source[row_number]
            if source_row is not None:
                output.at[row_number, column] = backfill.at[source_row, column]
        if overwrite_conflicts:
            for row_number in conflict_rows:
                source_row = matched_source[row_number]
                if source_row is not None:
                    output.at[row_number, column] = backfill.at[source_row, column]

        after_missing_count = int(_missing_mask(output[column]).sum())
        before_missing_count = int(before_missing.sum())
        missingness_rows.append(
            {
                "column": column,
                "existed_in_primary": existed,
                "primary_rows": len(primary),
                "before_missing": before_missing_count,
                "before_missing_fraction": (
                    before_missing_count / len(primary) if len(primary) else 0.0
                ),
                "after_missing": after_missing_count,
                "after_missing_fraction": (
                    after_missing_count / len(primary) if len(primary) else 0.0
                ),
                "filled_missing_cells": len(fill_rows),
                "matched_backfill_nonmissing_cells": source_nonmissing_matched,
                "conflicting_nonmissing_cells": len(conflict_rows),
                "overwritten_conflicting_cells": (
                    len(conflict_rows) if overwrite_conflicts else 0
                ),
            }
        )

    primary_value_hash_after = _primary_nonmissing_hash(primary, output)
    primary_values_preserved = primary_value_hash_before == primary_value_hash_after
    if not overwrite_conflicts and not primary_values_preserved:
        raise RuntimeError("Primary nonmissing-value preservation invariant failed")
    if len(output) != len(primary):
        raise RuntimeError("Primary row-count preservation invariant failed")
    if list(output.columns[: len(primary.columns)]) != list(primary.columns):
        raise RuntimeError("Primary column-order preservation invariant failed")

    unmatched_primary = _unmatched_primary_rows(
        primary, matched_source, specs, primary_keys
    )
    unmatched_backfill = _unmatched_backfill_rows(
        backfill, matched_source, specs, primary_keys, backfill_keys
    )
    key_audit = pd.DataFrame(key_audit_rows)
    missingness = pd.DataFrame(missingness_rows)
    conflict_columns = [
        "primary_row_number",
        "backfill_row_number",
        "match_key_type",
        "match_key_value",
        "column",
        "primary_value",
        "backfill_value",
        "action",
    ]
    conflict_table = pd.DataFrame(conflicts, columns=conflict_columns)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.mkdir(parents=True, exist_ok=True)
    artifact_paths = {
        "key_cardinality_audit": audit_path / "key_cardinality_audit.csv",
        "column_missingness_audit": audit_path / "column_missingness_audit.csv",
        "conflicts": audit_path / "conflicts.csv",
        "unmatched_primary_summary": audit_path / "unmatched_primary_summary.csv",
        "unmatched_backfill_summary": audit_path / "unmatched_backfill_summary.csv",
    }
    output.to_csv(output_path, index=False)
    key_audit.to_csv(artifact_paths["key_cardinality_audit"], index=False)
    missingness.to_csv(artifact_paths["column_missingness_audit"], index=False)
    conflict_table.to_csv(artifact_paths["conflicts"], index=False)
    unmatched_primary.to_csv(
        artifact_paths["unmatched_primary_summary"], index=False
    )
    unmatched_backfill.to_csv(
        artifact_paths["unmatched_backfill_summary"], index=False
    )

    matched_rows = sum(source_row is not None for source_row in matched_source)
    artifacts = {
        name: {"path": str(path), "sha256": _sha256_file(path)}
        for name, path in artifact_paths.items()
    }
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "primary": {
                "path": str(primary_path),
                "sha256": _sha256_file(primary_path),
                "rows": len(primary),
                "columns": len(primary.columns),
            },
            "pk_enriched": {
                "path": str(backfill_path),
                "sha256": _sha256_file(backfill_path),
                "rows": len(backfill),
                "columns": len(backfill.columns),
            },
        },
        "output": {
            "path": str(output_path),
            "sha256": _sha256_file(output_path),
            "rows": len(output),
            "columns": len(output.columns),
        },
        "policy": {
            "overwrite_conflicts": overwrite_conflicts,
            "default_column_detection": backfill_columns is None,
            "backfill_columns": columns,
            "key_priority": [spec.name for spec in specs],
            "ambiguous_source_matches": "reject",
        },
        "counts": {
            "primary_rows": len(primary),
            "pk_enriched_rows": len(backfill),
            "output_rows": len(output),
            "matched_primary_rows": matched_rows,
            "fallback_matched_primary_rows": sum(
                count
                for key_type, count in selected_counts.items()
                if key_type and key_type != specs[0].name
            ),
            "unmatched_primary_rows": len(unmatched_primary),
            "matched_unique_backfill_rows": len(
                {row for row in matched_source if row is not None}
            ),
            "unmatched_backfill_rows": len(unmatched_backfill),
            "added_columns": sum(column not in primary.columns for column in columns),
            "filled_missing_cells": int(missingness["filled_missing_cells"].sum()),
            "conflicting_nonmissing_cells": len(conflict_table),
            "overwritten_conflicting_cells": int(
                missingness["overwritten_conflicting_cells"].sum()
            ),
            "matches_by_key": {
                key: int(value)
                for key, value in selected_counts.items()
                if key
            },
        },
        "invariants": {
            "row_count_preserved": len(output) == len(primary),
            "primary_row_order_preserved": True,
            "primary_columns_preserved_as_prefix": list(
                output.columns[: len(primary.columns)]
            )
            == list(primary.columns),
            "primary_nonmissing_value_sha256_before": primary_value_hash_before,
            "primary_nonmissing_value_sha256_after": primary_value_hash_after,
            "primary_nonmissing_values_preserved": primary_values_preserved,
            "nonmissing_changes_require_explicit_opt_in": (
                primary_values_preserved or overwrite_conflicts
            ),
            "output_columns_unique": output.columns.is_unique,
            "ambiguous_matches_accepted": 0,
        },
        "artifacts": artifacts,
    }
    manifest_path = audit_path / "merge_pk_context_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


__all__ = ["DEFAULT_KEY_COLUMNS", "merge_pk_context"]
