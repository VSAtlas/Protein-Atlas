from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence


FINAL_SCORE_CONTRACT_VERSION = "final_score_materialized_v2"
DEFAULT_EXISTING_SCORE_COLUMN = "final_score"
CANONICAL_CONSENSUS_Z_COLUMN = "z_vs_decoys_consensus"
COMPARISON_RUN_CONSENSUS_Z_COLUMN = "z_vs_compare_run_consensus"
DEFAULT_FALLBACK_COLUMN = CANONICAL_CONSENSUS_Z_COLUMN

_ALLOWED_CONSENSUS_Z_FALLBACKS = {
    CANONICAL_CONSENSUS_Z_COLUMN,
    COMPARISON_RUN_CONSENSUS_Z_COLUMN,
}
ACCEPTED_FINAL_SCORE_SOURCES = frozenset(
    {
        "z_vs_decoys_blend",
        CANONICAL_CONSENSUS_Z_COLUMN,
        COMPARISON_RUN_CONSENSUS_Z_COLUMN,
    }
)
ACCEPTED_ATLAS_SCORE_SOURCES = frozenset(
    {
        *ACCEPTED_FINAL_SCORE_SOURCES,
        "z_vs_decoys_consensus_pre",
        "consensus_z_score",
        "exact_reference_run_fda_and_dud_consensus_z",
        "borrowed_reference_run_dud_consensus_z",
    }
)

_MISSING_TEXT = {"", "na", "n/a", "nan", "none", "null"}
_RAW_PERCENTILE_TOKENS = (
    "consensus_score_fallback",
    "rank_percentile",
    "raw_consensus",
    "raw consensus",
    "consensus percentile",
    "consensus_percentile",
)
_STAGE2_SOURCE = re.compile(
    r"(?:^|[^a-z0-9])(?:stage[_ -]?2|post[_ -]?dock(?:ed|ing)?)(?:$|[^a-z0-9])",
    re.IGNORECASE,
)
_IDENTIFIER_COLUMNS = (
    "pdb_id",
    "ligand_base",
    "drug_id",
    "target_id",
    "target_gene",
    "run_id",
)


def _clean(value: Any) -> str:
    return str(value or "").strip()


def _is_missing(value: Any) -> bool:
    return _clean(value).lower() in _MISSING_TEXT


def _finite_text(value: Any) -> str:
    text = _clean(value)
    if _is_missing(text):
        return ""
    try:
        parsed = float(text)
    except (TypeError, ValueError):
        return ""
    return text if math.isfinite(parsed) else ""


def _contains_raw_percentile_provenance(value: Any) -> bool:
    text = _clean(value).lower()
    return any(token in text for token in _RAW_PERCENTILE_TOKENS)


def _indicates_stage2_or_post_docked(value: Any) -> bool:
    return bool(_STAGE2_SOURCE.search(_clean(value)))


def _indicates_direct_rescore(value: Any) -> bool:
    return canonical_final_score_source(value) == "z_vs_decoys_blend"


def canonical_final_score_source(value: Any) -> str:
    normalized = _clean(value).lower()
    if normalized in {
        CANONICAL_CONSENSUS_Z_COLUMN,
        "legacy_exact_match_z_vs_decoys_consensus",
    }:
        return CANONICAL_CONSENSUS_Z_COLUMN
    if normalized in {
        COMPARISON_RUN_CONSENSUS_Z_COLUMN,
        "comparison_run_consensus_fallback",
    }:
        return COMPARISON_RUN_CONSENSUS_Z_COLUMN
    if normalized in {
        "z_vs_decoys_blend",
        "legacy_exact_match_z_vs_decoys_blend",
    }:
        return "z_vs_decoys_blend"
    return ""


def canonical_atlas_score_source(value: Any) -> str:
    final_source = canonical_final_score_source(value)
    if final_source:
        return final_source
    normalized = _clean(value).lower()
    return normalized if normalized in ACCEPTED_ATLAS_SCORE_SOURCES else ""


def _is_decoy_context(row: Mapping[str, Any]) -> bool:
    if _clean(row.get("is_decoy")).lower() in {"1", "true", "yes"}:
        return True
    for column in (
        "run_mode",
        "library",
        "_internal_spd_addons_run_mode",
        "_internal_spd_addons_library_original",
    ):
        token = _clean(row.get(column)).lower()
        if token in {"dud", "decoy", "decoys"} or token.startswith(
            ("dud_", "decoy_", "decoys_")
        ):
            return True
    return False


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
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


def _write_csv(
    path: Path,
    fields: Sequence[str],
    rows: Sequence[Mapping[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=True, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _row_identity_payload(
    rows: Sequence[Mapping[str, Any]], fields: Sequence[str]
) -> list[list[Any]]:
    identifiers = [field for field in _IDENTIFIER_COLUMNS if field in fields]
    if not identifiers:
        return [[index] for index in range(len(rows))]
    return [
        [index, *[row.get(field, "") for field in identifiers]]
        for index, row in enumerate(rows)
    ]


def _non_score_nonmissing_hash(
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    mutable_fields: set[str],
) -> str:
    payload = [
        [row_index, field, row.get(field, "")]
        for row_index, row in enumerate(rows)
        for field in fields
        if field not in mutable_fields and not _is_missing(row.get(field, ""))
    ]
    return _json_sha256(payload)


def _fallback_provenance_columns(
    fields: Sequence[str], fallback_column: str
) -> list[str]:
    prefix = f"{fallback_column}_"
    return [field for field in fields if field.startswith(prefix)]


def _identifier_fields(row: Mapping[str, Any], fields: Sequence[str]) -> dict[str, str]:
    return {
        field: _clean(row.get(field, ""))
        for field in _IDENTIFIER_COLUMNS
        if field in fields
    }


def _raw_detection_rows(
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    *,
    existing_score_column: str,
    fallback_columns: Sequence[str],
) -> list[dict[str, Any]]:
    fallback_provenance = {
        column: _fallback_provenance_columns(fields, column)
        for column in fallback_columns
    }
    detections: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        existing = _finite_text(row.get(existing_score_column))
        relevant: list[tuple[str, str]] = []
        if existing:
            explicit_source = _clean(row.get("final_score_source"))
            if explicit_source:
                relevant.append(("final_score_source", explicit_source))
            elif _indicates_stage2_or_post_docked(row.get("z_selected_source")):
                relevant.append(
                    ("z_selected_source", _clean(row.get("z_selected_source")))
                )
        for fallback_column in fallback_columns:
            if not _finite_text(row.get(fallback_column)):
                continue
            relevant.extend(
                (field, _clean(row.get(field)))
                for field in fallback_provenance[fallback_column]
            )
        for field, value in relevant:
            if value and _contains_raw_percentile_provenance(value):
                detections.append(
                    {
                        "row_number": index + 1,
                        **_identifier_fields(row, fields),
                        "provenance_column": field,
                        "provenance_value": value,
                        "reason": "raw_percentile_provenance_token",
                    }
                )
    return detections


def _orphaned_rescore_rows(
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    *,
    existing_score_column: str,
) -> list[dict[str, Any]]:
    raw_columns = [
        column
        for column in ("scorch_composite", "SCORCH_score_used")
        if column in fields
    ]
    detections: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        if _is_decoy_context(row):
            continue
        raw_values = [
            (column, _finite_text(row.get(column)))
            for column in raw_columns
            if _finite_text(row.get(column))
        ]
        if not raw_values:
            continue
        existing = _finite_text(row.get(existing_score_column))
        explicit_source = _clean(row.get("final_score_source"))
        inferred_source = _clean(row.get("z_selected_source"))
        selected_source = explicit_source or inferred_source
        if existing and _indicates_direct_rescore(selected_source):
            continue
        detections.append(
            {
                "row_number": index + 1,
                **_identifier_fields(row, fields),
                "raw_rescore_columns": "|".join(item[0] for item in raw_values),
                "raw_rescore_values": "|".join(item[1] for item in raw_values),
                "existing_score_value": existing,
                "existing_score_source": selected_source,
                "reason": "raw_rescore_not_selected_as_direct_final_score",
            }
        )
    return detections


def _failed_manifest(
    *,
    dataset_path: Path,
    out_path: Path,
    audit_path: Path,
    existing_score_column: str,
    fallback_column: str,
    reason: str,
    counts: Mapping[str, int],
) -> dict[str, Any]:
    return {
        "status": "failed",
        "contract_version": FINAL_SCORE_CONTRACT_VERSION,
        "reason": reason,
        "dataset": str(dataset_path),
        "output": str(out_path),
        "audit_dir": str(audit_path),
        "existing_score_column": existing_score_column,
        "fallback_column": fallback_column,
        "input_sha256": _sha256_file(dataset_path),
        "counts": dict(counts),
        "raw_consensus_score_allowed": False,
        "invented_z_scores_allowed": False,
    }


def materialize_final_score(
    dataset: str | Path,
    out: str | Path,
    audit_dir: str | Path,
    *,
    existing_score_column: str = DEFAULT_EXISTING_SCORE_COLUMN,
    fallback_column: str = DEFAULT_FALLBACK_COLUMN,
    expected_input_sha256: str | None = None,
) -> dict[str, Any]:
    dataset_path = Path(dataset)
    out_path = Path(out)
    audit_path = Path(audit_dir)
    score_column = _clean(existing_score_column)
    fallback = _clean(fallback_column)
    if not score_column or not fallback:
        raise ValueError("score and fallback column names must be nonblank")
    if fallback not in _ALLOWED_CONSENSUS_Z_FALLBACKS:
        raise ValueError(
            "final_score fallback must be a DUD-standardized consensus z-score; "
            f"received {fallback!r}"
        )

    input_fields, input_rows = _read_csv(dataset_path)
    consensus_input_columns = [
        column
        for column in (
            CANONICAL_CONSENSUS_Z_COLUMN,
            COMPARISON_RUN_CONSENSUS_Z_COLUMN,
        )
        if column in input_fields
    ]
    if fallback not in input_fields:
        raise ValueError(
            f"dataset missing requested fallback column: {fallback}"
        )

    output_fields = list(input_fields)
    for field in (
        score_column,
        "final_score_source",
        "final_score_contract_version",
    ):
        if field not in output_fields:
            output_fields.append(field)

    audit_path.mkdir(parents=True, exist_ok=True)
    expected_hash = _clean(expected_input_sha256).lower()
    if expected_hash and not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
        raise ValueError("--expected-input-sha256 must be a 64-character hex digest")
    input_hash = _sha256_file(dataset_path)
    if expected_hash and input_hash != expected_hash:
        failure_manifest = _failed_manifest(
            dataset_path=dataset_path,
            out_path=out_path,
            audit_path=audit_path,
            existing_score_column=score_column,
            fallback_column=fallback,
            reason="input_sha256_mismatch",
            counts={"dataset_rows": len(input_rows)},
        )
        failure_manifest["expected_input_sha256"] = expected_hash
        _write_json(
            audit_path / "materialize_final_score_manifest.json", failure_manifest
        )
        raise ValueError(
            f"materializer input SHA mismatch: expected {expected_hash}, found {input_hash}"
        )

    raw_detections = _raw_detection_rows(
        input_rows,
        input_fields,
        existing_score_column=score_column,
        fallback_columns=consensus_input_columns,
    )
    raw_detection_fields = (
        "row_number",
        *[field for field in _IDENTIFIER_COLUMNS if field in input_fields],
        "provenance_column",
        "provenance_value",
        "reason",
    )
    _write_csv(
        audit_path / "raw_percentile_detections.csv",
        raw_detection_fields,
        raw_detections,
    )
    if raw_detections:
        failure_manifest = _failed_manifest(
            dataset_path=dataset_path,
            out_path=out_path,
            audit_path=audit_path,
            existing_score_column=score_column,
            fallback_column=fallback,
            reason="raw_percentile_fallback_detected",
            counts={
                "dataset_rows": len(input_rows),
                "raw_percentile_detections": len(raw_detections),
            },
        )
        _write_json(
            audit_path / "materialize_final_score_manifest.json", failure_manifest
        )
        raise ValueError(
            f"raw-percentile fallback provenance detected in {len(raw_detections)} rows"
        )

    orphaned_rescores = _orphaned_rescore_rows(
        input_rows,
        input_fields,
        existing_score_column=score_column,
    )
    orphaned_rescore_fields = (
        "row_number",
        *[field for field in _IDENTIFIER_COLUMNS if field in input_fields],
        "raw_rescore_columns",
        "raw_rescore_values",
        "existing_score_value",
        "existing_score_source",
        "reason",
    )
    _write_csv(
        audit_path / "orphaned_rescore_detections.csv",
        orphaned_rescore_fields,
        orphaned_rescores,
    )
    if orphaned_rescores:
        failure_manifest = _failed_manifest(
            dataset_path=dataset_path,
            out_path=out_path,
            audit_path=audit_path,
            existing_score_column=score_column,
            fallback_column=fallback,
            reason="orphaned_raw_rescore_detected",
            counts={
                "dataset_rows": len(input_rows),
                "orphaned_raw_rescore_rows": len(orphaned_rescores),
            },
        )
        failure_manifest["expected_input_sha256"] = expected_hash
        _write_json(
            audit_path / "materialize_final_score_manifest.json",
            failure_manifest,
        )
        raise ValueError(
            f"raw rescoring is orphaned from final_score in {len(orphaned_rescores)} rows"
        )

    output_rows: list[dict[str, str]] = []
    rejected_rows: list[dict[str, Any]] = []
    existing_finite = 0
    existing_explicit = 0
    existing_inferred = 0
    fallback_finite = 0
    fallback_selected = 0
    for index, source_row in enumerate(input_rows):
        row = dict(source_row)
        existing = _finite_text(source_row.get(score_column))
        fallback_value = _finite_text(source_row.get(fallback))
        explicit_source = _clean(source_row.get("final_score_source"))
        z_selected_source = _clean(source_row.get("z_selected_source"))
        canonical_explicit_source = canonical_final_score_source(explicit_source)
        canonical_inferred_source = canonical_final_score_source(z_selected_source)
        if existing:
            existing_finite += 1
        if fallback_value:
            fallback_finite += 1

        selected_score = ""
        selected_source = ""
        if existing and canonical_explicit_source:
            selected_score = existing
            selected_source = canonical_explicit_source
            existing_explicit += 1
        elif existing and canonical_inferred_source:
            selected_score = existing
            selected_source = canonical_inferred_source
            existing_inferred += 1
        elif existing:
            action = (
                "rejected_replaced_by_fallback"
                if fallback_value
                else "rejected_left_missing"
            )
            rejected_rows.append(
                {
                    "row_number": index + 1,
                    **_identifier_fields(source_row, input_fields),
                    "rejected_score_column": score_column,
                    "rejected_score_value": existing,
                    "final_score_source": explicit_source,
                    "z_selected_source": z_selected_source,
                    "fallback_column": fallback,
                    "fallback_value": fallback_value,
                    "action": action,
                    "reason": "finite_existing_score_lacks_accepted_z_provenance",
                }
            )

        if not selected_score and fallback_value:
            selected_score = fallback_value
            selected_source = fallback
            fallback_selected += 1

        row[score_column] = selected_score
        row["final_score_source"] = selected_source
        row["final_score_contract_version"] = FINAL_SCORE_CONTRACT_VERSION
        output_rows.append(row)

    selected_without_source = [
        index + 1
        for index, row in enumerate(output_rows)
        if _finite_text(row.get(score_column))
        and _is_missing(row.get("final_score_source"))
    ]
    if selected_without_source:
        failure_manifest = _failed_manifest(
            dataset_path=dataset_path,
            out_path=out_path,
            audit_path=audit_path,
            existing_score_column=score_column,
            fallback_column=fallback,
            reason="selected_score_missing_source",
            counts={
                "dataset_rows": len(input_rows),
                "selected_scores_missing_source": len(selected_without_source),
            },
        )
        _write_json(
            audit_path / "materialize_final_score_manifest.json", failure_manifest
        )
        raise AssertionError("selected final score lacks source provenance")

    mutable_fields = {
        score_column,
        "final_score_source",
        "final_score_contract_version",
    }
    before_non_score_hash = _non_score_nonmissing_hash(
        input_rows, input_fields, mutable_fields
    )
    after_non_score_hash = _non_score_nonmissing_hash(
        output_rows, input_fields, mutable_fields
    )
    row_order_hash_before = _json_sha256(
        _row_identity_payload(input_rows, input_fields)
    )
    row_order_hash_after = _json_sha256(
        _row_identity_payload(output_rows, input_fields)
    )
    invariants = {
        "row_count_preserved": len(output_rows) == len(input_rows),
        "row_order_preserved": row_order_hash_after == row_order_hash_before,
        "input_columns_preserved_as_prefix": output_fields[: len(input_fields)]
        == input_fields,
        "non_score_nonmissing_values_preserved": after_non_score_hash
        == before_non_score_hash,
        "selected_scores_have_sources": not selected_without_source,
        "input_sha256_matches_expected": not expected_hash
        or input_hash == expected_hash,
        "orphaned_raw_rescore_absent": not orphaned_rescores,
        "raw_percentile_fallback_absent": not raw_detections,
        "requested_fallback_column_used_exactly": fallback in input_fields,
        "fallback_provenance_columns_preserved": all(
            field in output_fields
            for field in _fallback_provenance_columns(input_fields, fallback)
        ),
    }
    failed_invariants = [name for name, passed in invariants.items() if not passed]
    if failed_invariants:
        raise AssertionError(
            "final-score materialization invariants failed: "
            + ", ".join(failed_invariants)
        )

    _write_csv(out_path, output_fields, output_rows)
    written_fields, written_rows = _read_csv(out_path)
    if written_fields != output_fields or len(written_rows) != len(output_rows):
        raise AssertionError("written output failed schema or row-count verification")
    written_non_score_hash = _non_score_nonmissing_hash(
        written_rows, input_fields, mutable_fields
    )
    if written_non_score_hash != before_non_score_hash:
        raise AssertionError("written output changed non-score nonmissing cells")

    rejected_fields = (
        "row_number",
        *[field for field in _IDENTIFIER_COLUMNS if field in input_fields],
        "rejected_score_column",
        "rejected_score_value",
        "final_score_source",
        "z_selected_source",
        "fallback_column",
        "fallback_value",
        "action",
        "reason",
    )
    _write_csv(
        audit_path / "rejected_final_scores.csv",
        rejected_fields,
        rejected_rows,
    )

    source_counts = Counter(
        _clean(row.get("final_score_source")) or "missing" for row in output_rows
    )
    scored_rows = sum(bool(_finite_text(row.get(score_column))) for row in output_rows)
    source_count_rows = [
        {
            "final_score_source": source,
            "n_rows": count,
            "fraction_of_all_rows": count / len(output_rows) if output_rows else 0.0,
            "fraction_of_scored_rows": (
                count / scored_rows if scored_rows and source != "missing" else 0.0
            ),
        }
        for source, count in sorted(source_counts.items())
    ]
    _write_csv(
        audit_path / "final_score_source_counts.csv",
        (
            "final_score_source",
            "n_rows",
            "fraction_of_all_rows",
            "fraction_of_scored_rows",
        ),
        source_count_rows,
    )

    coverage = {
        "dataset_rows": len(input_rows),
        "existing_finite_scores": existing_finite,
        "preserved_explicit_source_scores": existing_explicit,
        "preserved_inferred_accepted_scores": existing_inferred,
        "rejected_unprovenanced_existing_scores": len(rejected_rows),
        "finite_fallback_scores": fallback_finite,
        "selected_fallback_scores": fallback_selected,
        "materialized_finite_scores": scored_rows,
        "materialized_missing_scores": len(output_rows) - scored_rows,
        "materialized_score_coverage": (
            scored_rows / len(output_rows) if output_rows else 0.0
        ),
    }
    _write_csv(audit_path / "coverage.csv", tuple(coverage), [coverage])
    _write_csv(
        audit_path / "invariants.csv",
        ("invariant", "passed"),
        [
            {"invariant": name, "passed": bool(passed)}
            for name, passed in invariants.items()
        ],
    )

    fallback_provenance = sorted(
        _fallback_provenance_columns(input_fields, fallback)
    )
    manifest: dict[str, Any] = {
        "status": "ok",
        "contract_version": FINAL_SCORE_CONTRACT_VERSION,
        "dataset": str(dataset_path),
        "output": str(out_path),
        "audit_dir": str(audit_path),
        "existing_score_column": score_column,
        "fallback_input_column": fallback,
        "consensus_input_columns": consensus_input_columns,
        "fallback_column": fallback,
        "expected_input_sha256": expected_hash,
        "fallback_provenance_columns_carried": fallback_provenance,
        "policy": {
            "existing_score_requires_accepted_source": True,
            "accepted_final_score_sources": sorted(ACCEPTED_FINAL_SCORE_SOURCES),
            "accepted_direct_rescore_source": "z_vs_decoys_blend",
            "fallback_source": fallback,
            "same_run_and_comparison_run_consensus_are_distinct": True,
            "comparison_run_consensus_semantics": (
                "same-PDB raw consensus standardized by the referenced DUD decoy mean and SD"
            ),
            "raw_consensus_score_allowed": False,
            "raw_rescore_without_direct_final_score_allowed": False,
            "invented_z_scores_allowed": False,
            "unprovenanced_existing_score_action": "quarantine_then_use_finite_z_fallback_if_available",
        },
        "coverage": coverage,
        "source_counts": dict(sorted(source_counts.items())),
        "hashes": {
            "input_sha256": input_hash,
            "output_sha256": _sha256_file(out_path),
            "row_order_sha256_before": row_order_hash_before,
            "row_order_sha256_after": row_order_hash_after,
            "non_score_nonmissing_sha256_before": before_non_score_hash,
            "non_score_nonmissing_sha256_after": written_non_score_hash,
            "score_contract_sha256": _json_sha256(
                [
                    [
                        index,
                        row.get(score_column, ""),
                        row.get(fallback, ""),
                        row.get("final_score_source", ""),
                        row.get("final_score_contract_version", ""),
                    ]
                    for index, row in enumerate(output_rows)
                ]
            ),
        },
        "invariants": invariants,
        "artifacts": {
            "output": str(out_path),
            "coverage": str(audit_path / "coverage.csv"),
            "source_counts": str(audit_path / "final_score_source_counts.csv"),
            "rejected_scores": str(audit_path / "rejected_final_scores.csv"),
            "orphaned_rescore_detections": str(
                audit_path / "orphaned_rescore_detections.csv"
            ),
            "raw_percentile_detections": str(
                audit_path / "raw_percentile_detections.csv"
            ),
            "invariants": str(audit_path / "invariants.csv"),
        },
    }
    manifest_path = audit_path / "materialize_final_score_manifest.json"
    manifest["artifacts"]["manifest"] = str(manifest_path)
    _write_json(manifest_path, manifest)
    return manifest


__all__ = [
    "ACCEPTED_ATLAS_SCORE_SOURCES",
    "ACCEPTED_FINAL_SCORE_SOURCES",
    "CANONICAL_CONSENSUS_Z_COLUMN",
    "COMPARISON_RUN_CONSENSUS_Z_COLUMN",
    "DEFAULT_EXISTING_SCORE_COLUMN",
    "DEFAULT_FALLBACK_COLUMN",
    "FINAL_SCORE_CONTRACT_VERSION",
    "canonical_atlas_score_source",
    "canonical_final_score_source",
    "materialize_final_score",
]
