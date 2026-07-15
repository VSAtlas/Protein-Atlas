"""Deterministic, provenance-preserving release-attempt selection."""

from __future__ import annotations

import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence


AttemptKey = tuple[str, str, str, str]


def _text(value: Any) -> str:
    return str(value or "").strip()


def _normalize_ph(value: Any) -> str:
    token = _text(value)
    if token.lower() in {"base", "none", "null"} or token.lower().endswith(".csv"):
        return ""
    return token


def _selector_key(value: Mapping[str, Any]) -> AttemptKey:
    return (
        _text(value.get("pdb_id")).upper(),
        _text(value.get("variant")).upper(),
        _normalize_ph(value.get("ph_label")),
        _text(value.get("ligand_canonical_id")),
    )


def _key_label(key: AttemptKey) -> str:
    pdb_id, variant, ph_label, ligand = key
    return (
        f"pdb_id={pdb_id!r}, variant={variant!r}, ph_label={ph_label!r}, "
        f"ligand={ligand!r}"
    )


def _selectors(
    values: Sequence[Mapping[str, Any]],
) -> dict[AttemptKey, tuple[int, Mapping[str, Any]]]:
    result: dict[AttemptKey, tuple[int, Mapping[str, Any]]] = {}
    for manifest_index, value in enumerate(values, start=1):
        key = _selector_key(value)
        if key in result:
            raise ValueError(f"duplicate attempt selector for {_key_label(key)}")
        result[key] = (manifest_index, value)
    return result


def _completion_relpath(path: str, docked: Path) -> str:
    try:
        return Path(path).relative_to(docked).as_posix()
    except ValueError as exc:
        raise ValueError(
            f"completion attempt path is outside paths.docked: {path}"
        ) from exc


def _select_completion_attempts(
    connection: sqlite3.Connection,
    run_id: str,
    docked: Path,
    selectors: Mapping[AttemptKey, tuple[int, Mapping[str, Any]]],
) -> set[AttemptKey]:
    groups: dict[AttemptKey, list[dict[str, Any]]] = defaultdict(list)
    rows = connection.execute(
        """SELECT a.attempt_id, r.pdb_id, r.variant, r.ph_label,
                  l.canonical_id, c.completion_path, c.completion_sha256
           FROM docking_attempts a
           JOIN pair_cells p ON p.pair_cell_id=a.pair_cell_id
           JOIN receptor_contexts r
             ON r.receptor_context_id=p.receptor_context_id
           JOIN ligands l ON l.ligand_id=p.ligand_id
           JOIN completion_records c
             ON c.completion_record_id=a.completion_record_id
           WHERE r.run_id=? AND a.status='completed'
           ORDER BY a.attempt_id""",
        (run_id,),
    ).fetchall()
    for row in rows:
        key = (str(row[1]), str(row[2]), str(row[3]), str(row[4]))
        groups[key].append(
            {
                "attempt_id": int(row[0]),
                "completion_relpath": _completion_relpath(str(row[5]), docked),
                "completion_sha256": str(row[6]).lower(),
            }
        )

    explicitly_used: set[AttemptKey] = set()
    for key, candidates in groups.items():
        selected = selectors.get(key)
        has_explicit = bool(
            selected
            and (
                _text(selected[1].get("completion_relpath"))
                or _text(selected[1].get("completion_sha256"))
            )
        )
        if len(candidates) > 1 and not has_explicit:
            choices = ", ".join(
                f"{item['completion_relpath']} sha256={item['completion_sha256']}"
                for item in candidates
            )
            raise ValueError(
                "multiple successful completion attempts for "
                f"{_key_label(key)}; add an explicit attempt_selections entry; "
                f"candidates: {choices}"
            )
        manifest_index: int | None = None
        if has_explicit and selected is not None:
            manifest_index, selector = selected
            expected_path = _text(selector.get("completion_relpath"))
            expected_sha = _text(selector.get("completion_sha256")).lower()
            matches = [
                item
                for item in candidates
                if item["completion_relpath"] == expected_path
                and item["completion_sha256"] == expected_sha
            ]
            if len(matches) != 1:
                raise ValueError(
                    "completion selector matched "
                    f"{len(matches)} successful attempts for {_key_label(key)}"
                )
            chosen = matches[0]
            method = "manifest_explicit"
            explicitly_used.add(key)
        else:
            chosen = candidates[0]
            method = "only_successful_attempt"
        connection.execute(
            """UPDATE docking_attempts
               SET selected_for_release=1, selection_method=?,
                   selection_manifest_index=?
               WHERE attempt_id=?""",
            (method, manifest_index, chosen["attempt_id"]),
        )

    for key, (_, selector) in selectors.items():
        has_completion = bool(
            _text(selector.get("completion_relpath"))
            or _text(selector.get("completion_sha256"))
        )
        if has_completion and key not in explicitly_used:
            raise ValueError(
                "completion selector does not match a successful attempt for "
                f"{_key_label(key)}"
            )
    return explicitly_used


def _materialize_result_attempt(
    connection: sqlite3.Connection, candidate: Mapping[str, Any]
) -> None:
    connection.execute(
        """UPDATE pair_cells SET final_status=?, failure_code=?, failure_reason=?,
           has_result=1, pose_valid=?, pose_validation_method=?,
           pose_validation_scope=?, pose_validation_thresholds_json=?,
           is_control=?, is_decoy=?, atlas_score=?,
           atlas_score_source=?, selected_docking_score=?, consensus_score=?,
           final_score=?, final_score_source=?,
           final_score_source_reconstructed=?,
           final_score_source_classification=?,
           final_score_source_evidence_json=?, final_rank=?, source_csv=?,
           result_json=? WHERE pair_cell_id=?""",
        (
            candidate["final_status"],
            candidate["failure_code"],
            candidate["failure_reason"],
            candidate["pose_valid"],
            candidate["pose_validation_method"],
            candidate["pose_validation_scope"],
            candidate["pose_validation_thresholds_json"],
            candidate["is_control"],
            candidate["is_decoy"],
            candidate["atlas_score"],
            candidate["atlas_score_source"],
            candidate["selected_docking_score"],
            candidate["consensus_score"],
            candidate["final_score"],
            candidate["final_score_source"],
            candidate["final_score_source_reconstructed"],
            candidate["final_score_source_classification"],
            candidate["final_score_source_evidence_json"],
            candidate["final_rank"],
            candidate["source_csv"],
            candidate["result_json"],
            candidate["pair_cell_id"],
        ),
    )


def _select_result_attempts(
    connection: sqlite3.Connection,
    run_id: str,
    selectors: Mapping[AttemptKey, tuple[int, Mapping[str, Any]]],
) -> set[AttemptKey]:
    groups: dict[AttemptKey, list[dict[str, Any]]] = defaultdict(list)
    rows = connection.execute(
        """SELECT a.result_attempt_id, a.pair_cell_id,
                  r.pdb_id, r.variant, r.ph_label, l.canonical_id,
                  a.source_row_number, a.result_sha256, a.final_status,
                  a.failure_code, a.failure_reason, a.pose_valid,
                  a.pose_validation_method, a.pose_validation_scope,
                  a.pose_validation_thresholds_json,
                  a.is_control, a.is_decoy, a.atlas_score,
                  a.atlas_score_source, a.selected_docking_score,
                  a.consensus_score, a.final_score, a.final_score_source,
                  a.final_score_source_reconstructed,
                  a.final_score_source_classification,
                  a.final_score_source_evidence_json,
                  a.final_rank, a.source_csv, a.result_json
           FROM result_attempts a
           JOIN pair_cells p ON p.pair_cell_id=a.pair_cell_id
           JOIN receptor_contexts r
             ON r.receptor_context_id=p.receptor_context_id
           JOIN ligands l ON l.ligand_id=p.ligand_id
           WHERE r.run_id=? ORDER BY a.result_attempt_id""",
        (run_id,),
    ).fetchall()
    fields = (
        "result_attempt_id",
        "pair_cell_id",
        "pdb_id",
        "variant",
        "ph_label",
        "ligand_canonical_id",
        "source_row_number",
        "result_sha256",
        "final_status",
        "failure_code",
        "failure_reason",
        "pose_valid",
        "pose_validation_method",
        "pose_validation_scope",
        "pose_validation_thresholds_json",
        "is_control",
        "is_decoy",
        "atlas_score",
        "atlas_score_source",
        "selected_docking_score",
        "consensus_score",
        "final_score",
        "final_score_source",
        "final_score_source_reconstructed",
        "final_score_source_classification",
        "final_score_source_evidence_json",
        "final_rank",
        "source_csv",
        "result_json",
    )
    for raw in rows:
        row = dict(zip(fields, raw, strict=True))
        key = (
            str(row["pdb_id"]),
            str(row["variant"]),
            str(row["ph_label"]),
            str(row["ligand_canonical_id"]),
        )
        groups[key].append(row)

    explicitly_used: set[AttemptKey] = set()
    for key, candidates in groups.items():
        selected = selectors.get(key)
        has_explicit = bool(
            selected
            and (
                selected[1].get("result_row_number") is not None
                or _text(selected[1].get("result_sha256"))
            )
        )
        if len(candidates) > 1 and not has_explicit:
            choices = ", ".join(
                f"record={item['source_row_number']} sha256={item['result_sha256']}"
                for item in candidates
            )
            raise ValueError(
                f"multiple result attempts for {_key_label(key)}; add an "
                f"explicit attempt_selections entry; candidates: {choices}"
            )
        manifest_index: int | None = None
        if has_explicit and selected is not None:
            manifest_index, selector = selected
            expected_row = selector.get("result_row_number")
            expected_sha = _text(selector.get("result_sha256")).lower()
            matches = [
                item
                for item in candidates
                if item["source_row_number"] == expected_row
                and str(item["result_sha256"]).lower() == expected_sha
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"result selector matched {len(matches)} attempts for "
                    f"{_key_label(key)}"
                )
            chosen = matches[0]
            method = "manifest_explicit"
            explicitly_used.add(key)
        else:
            chosen = candidates[0]
            method = "only_result_attempt"
        connection.execute(
            """UPDATE result_attempts
               SET selected_for_release=1, selection_method=?,
                   selection_manifest_index=?
               WHERE result_attempt_id=?""",
            (method, manifest_index, chosen["result_attempt_id"]),
        )
        _materialize_result_attempt(connection, chosen)

    for key, (_, selector) in selectors.items():
        has_result = bool(
            selector.get("result_row_number") is not None
            or _text(selector.get("result_sha256"))
        )
        if has_result and key not in explicitly_used:
            raise ValueError(
                f"result selector does not match an attempt for {_key_label(key)}"
            )
    return explicitly_used


def apply_attempt_selections(
    connection: sqlite3.Connection,
    run_id: str,
    docked: Path,
    selector_values: Sequence[Mapping[str, Any]],
) -> None:
    """Select at most one successful completion and result row per pair.

    Unique candidates are selected automatically. Any ambiguous candidate set
    fails closed unless a manifest selector identifies exactly one attempt by
    both its location/row number and content hash.
    """
    selectors = _selectors(selector_values)
    _select_completion_attempts(connection, run_id, docked, selectors)
    _select_result_attempts(connection, run_id, selectors)
