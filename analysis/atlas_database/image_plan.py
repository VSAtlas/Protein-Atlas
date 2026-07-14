"""Plan conference pose images from an immutable Atlas release database.

The planner never renders images and never substitutes secondary scores.  It
emits explicit ``atlas screenshot`` argument vectors plus structured gaps for
prerequisites that are absent from the release.
"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


def build_release_image_plan(
    database_path: Path,
    *,
    known_pairs: Sequence[Mapping[str, Any]] | None = None,
    image_output_root: str = "release_images",
) -> dict[str, Any]:
    """Return native/top-valid/known-pair screenshot commands by receptor context."""

    connection = sqlite3.connect(Path(database_path))
    connection.row_factory = sqlite3.Row
    try:
        if known_pairs is not None:
            supplied_known = list(known_pairs)
        else:
            supplied_known = _database_known_pairs(connection)
        contexts = connection.execute(
            """SELECT r.receptor_context_id, r.run_id, r.pdb_id, r.variant,
                      r.ph_label,
                      COALESCE(a.native_redock_status, r.native_redock_status)
                          AS native_redock_status,
                      COALESCE(a.native_redock_reason, r.native_redock_reason)
                          AS native_redock_reason
               FROM receptor_contexts r
               LEFT JOIN receptor_annotations a
                 ON a.receptor_context_id=r.receptor_context_id
               ORDER BY r.run_id, r.pdb_id, r.variant, r.ph_label"""
        ).fetchall()
        safe_output_root = _safe_relative_root(image_output_root)
        plans = [
            _context_plan(connection, row, supplied_known, safe_output_root)
            for row in contexts
        ]
    finally:
        connection.close()
    return {
        "schema_version": 2,
        "database_file": Path(database_path).name,
        "selection_policy": {
            "native_control": "explicit is_control pair with a valid pose result",
            "top_valid": "top 5 non-control, non-decoy, valid poses by final_score descending",
            "known_pair": "annotation-backed pair when valid and not already selected",
            "artifact_preflight": (
                "selected rows report a gap unless a pair-linked artifact is verified"
            ),
            "best_invalid": "highest final_score among non-control, non-decoy pose_valid=0 cells",
            "score_field": "final_score",
        },
        "contexts": plans,
        "summary": {
            "context_count": len(plans),
            "selection_count": sum(len(plan["selections"]) for plan in plans),
            "gap_count": sum(len(plan["gaps"]) for plan in plans),
        },
    }


def write_release_image_plan(
    database_path: Path,
    output_path: Path,
    *,
    known_pairs: Sequence[Mapping[str, Any]] | None = None,
    image_output_root: str = "release_images",
) -> dict[str, Any]:
    """Build and write a deterministic JSON image plan."""

    plan = build_release_image_plan(
        database_path,
        known_pairs=known_pairs,
        image_output_root=image_output_root,
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return plan


def _database_known_pairs(
    connection: sqlite3.Connection,
) -> list[Mapping[str, Any]]:
    rows = connection.execute(
        """SELECT receptor_context_id, ligand_canonical_id,
        evidence_reference, selection_label
        FROM known_pair_selections ORDER BY receptor_context_id"""
    ).fetchall()
    return [
        {
            "receptor_context_id": row["receptor_context_id"],
            "canonical_id": row["ligand_canonical_id"],
            "evidence_reference": row["evidence_reference"],
            "selection_label": row["selection_label"],
        }
        for row in rows
    ]


def _context_plan(
    connection: sqlite3.Connection,
    context: sqlite3.Row,
    known_pairs: Sequence[Mapping[str, Any]],
    image_output_root: str,
) -> dict[str, Any]:
    context_id = int(context["receptor_context_id"])
    rows = connection.execute(
        """SELECT p.pair_cell_id, p.final_status, p.has_result, p.pose_valid,
                  (SELECT COUNT(*) FROM artifacts ar
                   WHERE ar.pair_cell_id=p.pair_cell_id
                     AND ar.artifact_role='docking_pose') AS artifact_count,
                  (SELECT COUNT(*) FROM artifacts ar
                   WHERE ar.pair_cell_id=p.pair_cell_id AND ar.verified=1
                     AND ar.artifact_role='docking_pose')
                      AS verified_artifact_count,
                  p.is_control, p.is_decoy, p.final_score,
                  l.canonical_id, l.display_name
           FROM pair_cells p JOIN ligands l ON l.ligand_id = p.ligand_id
           WHERE p.receptor_context_id = ?
           ORDER BY p.pair_cell_id""",
        (context_id,),
    ).fetchall()
    selections: list[dict[str, Any]] = []
    gaps: list[dict[str, str]] = []
    selected_ligands: set[str] = set()
    if _text(context["variant"]) or _text(context["ph_label"]):
        gaps.append(
            _gap(
                "screenshot_context_resolution_unverified",
                "atlas screenshot uses variant/pH for archive matching but resolves the receptor through its generic PDB lookup",
            )
        )

    control_rows = [row for row in rows if int(row["is_control"] or 0) == 1]
    native_status = _text(context["native_redock_status"])
    if not native_status:
        gaps.append(
            _gap("native_status_missing", "native redock status is not recorded")
        )
    if len(control_rows) != 1:
        gaps.append(
            _gap(
                "native_control_not_unique",
                f"expected one control pair; found {len(control_rows)}",
            )
        )
    elif not _renderable(control_rows[0]):
        gaps.append(
            _gap(
                "native_control_pose_unavailable",
                "control pair lacks a successful valid pose",
            )
        )
    elif native_status:
        _select(
            selections,
            selected_ligands,
            context,
            control_rows[0],
            "native_control",
            image_output_root,
        )
        _record_artifact_gap(gaps, control_rows[0])

    valid_scored = [
        row
        for row in rows
        if int(row["is_control"] or 0) == 0
        and int(row["is_decoy"] or 0) == 0
        and _renderable(row)
        and row["final_score"] is not None
    ]
    valid_scored.sort(
        key=lambda row: (-float(row["final_score"]), _text(row["canonical_id"]))
    )
    if not valid_scored:
        gaps.append(
            _gap(
                "top_valid_final_score_missing",
                "no valid non-control ligand has final_score",
            )
        )
    elif len(valid_scored) < 5:
        gaps.append(
            _gap(
                "top_valid_fewer_than_five",
                f"only {len(valid_scored)} valid scored ligands are available",
            )
        )
    for index, row in enumerate(valid_scored[:5], 1):
        _select(
            selections,
            selected_ligands,
            context,
            row,
            f"top_valid_{index}",
            image_output_root,
        )
        _record_artifact_gap(gaps, row)

    invalid_scored = [
        row
        for row in rows
        if int(row["is_control"] or 0) == 0
        and int(row["is_decoy"] or 0) == 0
        and int(row["has_result"] or 0) == 1
        and row["pose_valid"] == 0
        and row["final_score"] is not None
    ]
    invalid_scored.sort(
        key=lambda row: (-float(row["final_score"]), _text(row["canonical_id"]))
    )
    if not invalid_scored:
        gaps.append(
            _gap(
                "best_invalid_not_available",
                "no scored non-control invalid pose is available",
            )
        )
    else:
        _select(
            selections,
            selected_ligands,
            context,
            invalid_scored[0],
            "best_invalid",
            image_output_root,
        )
        _record_artifact_gap(gaps, invalid_scored[0])

    known = _known_for_context(context, known_pairs)
    if not known:
        gaps.append(
            _gap(
                "known_pair_not_supplied",
                "no known-pair ligand was supplied for this receptor context",
            )
        )
    elif len(known) > 1:
        gaps.append(
            _gap(
                "known_pair_ambiguous",
                f"{len(known)} known-pair ligands were supplied; choose one representative",
            )
        )
    else:
        requested = known[0]
        ligand = _text(
            requested.get("canonical_id")
            or requested.get("ligand_id")
            or requested.get("drug_id")
        )
        if not ligand:
            gaps.append(
                _gap(
                    "known_pair_ligand_missing",
                    "known-pair selection lacks canonical_id/ligand_id/drug_id",
                )
            )
        else:
            match = next(
                (row for row in rows if _text(row["canonical_id"]) == ligand), None
            )
            if match is None:
                gaps.append(
                    _gap(
                        "known_pair_not_in_matrix",
                        f"known ligand {ligand} is absent from this receptor context",
                    )
                )
            elif not _renderable(match):
                gaps.append(
                    _gap(
                        "known_pair_pose_unavailable",
                        f"known ligand {ligand} lacks a successful valid pose",
                    )
                )
            elif ligand in selected_ligands:
                pass
            else:
                _select(
                    selections,
                    selected_ligands,
                    context,
                    match,
                    "known_pair",
                    image_output_root,
                )
                _record_artifact_gap(gaps, match)

    return {
        "receptor_context_id": context_id,
        "run_id": context["run_id"],
        "pdb_id": context["pdb_id"],
        "variant": context["variant"],
        "ph_label": context["ph_label"],
        "native_redock_status": context["native_redock_status"],
        "native_redock_reason": context["native_redock_reason"],
        "selections": selections,
        "gaps": gaps,
    }


def _known_for_context(
    context: sqlite3.Row, known_pairs: Sequence[Mapping[str, Any]]
) -> list[Mapping[str, Any]]:
    result = []
    for value in known_pairs:
        explicit_context = value.get("receptor_context_id")
        if explicit_context not in (None, ""):
            try:
                if int(str(explicit_context)) == int(context["receptor_context_id"]):
                    result.append(value)
            except (TypeError, ValueError):
                continue
            continue
        if _text(value.get("run_id")) not in ("", _text(context["run_id"])):
            continue
        if (
            _text(value.get("pdb_id") or value.get("target_id")).upper()
            != _text(context["pdb_id"]).upper()
        ):
            continue
        if _text(value.get("variant")) not in ("", _text(context["variant"])):
            continue
        if _text(value.get("ph_label") or value.get("ph")) not in (
            "",
            _text(context["ph_label"]),
        ):
            continue
        result.append(value)
    return result


def _renderable(row: sqlite3.Row) -> bool:
    return (
        int(row["has_result"] or 0) == 1
        and int(row["pose_valid"] or 0) == 1
        and _text(row["final_status"]).lower() in {"valid", "success", "completed"}
    )


def _record_artifact_gap(gaps: list[dict[str, str]], row: sqlite3.Row) -> None:
    artifact_count = int(row["artifact_count"] or 0)
    verified_count = int(row["verified_artifact_count"] or 0)
    if verified_count:
        return
    ligand = _text(row["canonical_id"])
    code = (
        "pair_artifact_not_indexed"
        if artifact_count == 0
        else "pair_artifact_not_verified"
    )
    gaps.append(
        _gap(
            code,
            f"{ligand} has no verified pair-linked docking_pose artifact; "
            "preflight the screenshot command before conference use",
        )
    )


def _select(
    selections: list[dict[str, Any]],
    selected_ligands: set[str],
    context: sqlite3.Row,
    row: sqlite3.Row,
    role: str,
    image_output_root: str,
) -> None:
    ligand = _text(row["canonical_id"])
    selected_ligands.add(ligand)
    argv = [
        "atlas",
        "screenshot",
        _text(context["run_id"]),
        "--pdb",
        _text(context["pdb_id"]),
        "--ligand",
        ligand,
    ]
    if _text(context["variant"]):
        argv.extend(["--variant", _text(context["variant"])])
    if _text(context["ph_label"]):
        argv.extend(["--ph", _text(context["ph_label"])])
    output_stem = (
        f"context-{int(context['receptor_context_id'])}/{role}-{_safe_token(ligand)}"
    )
    pocket_out = f"{image_output_root}/{output_stem}/pocket"
    full_out = f"{image_output_root}/{output_stem}/full"
    pocket_argv = [
        *argv,
        "--view-context",
        "pocket",
        "--gallery-html",
        "--output-subdir",
        pocket_out,
    ]
    full_argv = [
        *argv,
        "--view-context",
        "full",
        "--output-subdir",
        full_out,
    ]
    selections.append(
        {
            "role": role,
            "pair_cell_id": int(row["pair_cell_id"]),
            "ligand_id": ligand,
            "ligand_display_name": row["display_name"],
            "final_score": row["final_score"],
            "final_status": row["final_status"],
            "pose_valid": row["pose_valid"],
            "command_argv": pocket_argv,
            "pocket_three_view_command_argv": pocket_argv,
            "pocket_expected_views": ["front", "side", "top"],
            "full_context_command_argv": full_argv,
            "full_protein_thumbnail_expected_view": "front",
            "output_directories": {"pocket": pocket_out, "full": full_out},
            "output_subdirectories": {"pocket": pocket_out, "full": full_out},
        }
    )


def _gap(code: str, message: str) -> dict[str, str]:
    return {"code": code, "message": message}


def _text(value: Any) -> str:
    return str(value).strip() if value is not None else ""


def _safe_relative_root(value: str) -> str:
    path = Path(_text(value))
    if not _text(value) or path.is_absolute() or ".." in path.parts:
        raise ValueError(
            "image_output_root must be a non-empty relative path without '..'"
        )
    return path.as_posix().rstrip("/")


def _safe_token(value: str) -> str:
    token = "".join(char if char.isalnum() or char in "._-" else "-" for char in value)
    return token.strip("-") or "ligand"
