"""Deterministic, fail-closed plans for repeating one Atlas pipeline stage.

This module deliberately plans work only.  It never imports a docking engine,
starts a subprocess, or mutates an existing run.  A future executor can consume
the persisted pair ledger after checking the recorded prerequisite decisions.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sqlite3
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Mapping, Sequence

from analysis.atlas_database.score_source_contract import (
    score_source_pose_linkage_requirement,
)
from config.output_paths import run_output_dir


PLAN_SCHEMA = "atlas_stage_repeat_plan_v1"
PAIR_SCHEMA = "atlas_stage_repeat_pairs_v1"
SELECTION_ALGORITHM = "sha256_seeded_pair_order_v1"
SUPPORTED_STAGES = ("vina", "gnina", "scorch", "mmgbsa", "pose-validation")
MAX_STAGE_CPUS = 32
MAX_IN_MEMORY_PLAN_PAIRS = 25_000
PAIR_LEDGER_SHARD_ROWS = 5_000

_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_POSE_PATH_FIELDS = (
    "selected_pose_path",
    "final_score_pose_path",
    "final_pose_path",
    "docking_pose_path",
)
_POSE_STAGE_FIELDS = (
    "selected_pose_stage",
    "final_score_pose_stage",
    "docking_pose_stage",
)


class StageRepeatPlanError(ValueError):
    """Raised when a repeat plan cannot be defined without guessing."""


def _text(value: Any) -> str:
    return str(value or "").strip()


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_json_object(value: Any) -> dict[str, Any]:
    if isinstance(value, Mapping):
        return dict(value)
    token = _text(value)
    if not token:
        return {}
    try:
        parsed = json.loads(token)
    except (TypeError, ValueError, json.JSONDecodeError):
        return {}
    return dict(parsed) if isinstance(parsed, Mapping) else {}


def _normalize_ph(value: Any) -> str:
    token = _text(value)
    return "" if token.lower() in {"", "base", "none", "null"} else token


def _normalize_variant(value: Any) -> str:
    token = _text(value)
    return (
        "" if token.lower() in {"", "base", "legacy", "none", "null"} else token.upper()
    )


def _normalize_stage(value: Any) -> str:
    return _text(value).lower()


def _normalize_locator(value: Any) -> str:
    token = _text(value).replace("\\", "/")
    if not token:
        return ""
    try:
        return PurePosixPath(token).as_posix()
    except (TypeError, ValueError):
        return token


def _dedupe(values: Iterable[str]) -> list[str]:
    return sorted({_text(value) for value in values if _text(value)})


@dataclass(frozen=True)
class ArtifactRef:
    role: str
    stage: str = ""
    mode: str = ""
    original_path: str = ""
    archive_path: str = ""
    member_name: str = ""
    sha256: str = ""
    verified: bool = False
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def locator_tokens(self) -> set[str]:
        tokens = {
            _normalize_locator(self.original_path),
            _normalize_locator(self.member_name),
        }
        if self.archive_path and self.member_name:
            tokens.add(
                f"{_normalize_locator(self.archive_path)}::{_normalize_locator(self.member_name)}"
            )
        return {token for token in tokens if token}

    def public_record(self) -> dict[str, Any]:
        return {
            "role": self.role,
            "stage": self.stage,
            "mode": self.mode,
            "original_path": self.original_path,
            "archive_path": self.archive_path,
            "member_name": self.member_name,
            "sha256": self.sha256,
            "verified": self.verified,
        }


@dataclass(frozen=True)
class SelectedResultRef:
    result_attempt_id: int | None = None
    result_sha256: str = ""
    completion_record_id: int | None = None
    completion_sha256: str = ""
    completion_link_method: str = ""
    completion_link_evidence: Mapping[str, Any] = field(default_factory=dict)
    engine: str = ""
    stage: str = ""
    score_source: str = ""
    selected_docking_attempt_count: int = 0
    result_json: Mapping[str, Any] = field(default_factory=dict)

    def public_record(self) -> dict[str, Any]:
        return {
            "result_attempt_id": self.result_attempt_id,
            "result_sha256": self.result_sha256,
            "completion_record_id": self.completion_record_id,
            "completion_sha256": self.completion_sha256,
            "completion_link_method": self.completion_link_method,
            "completion_link_evidence": dict(self.completion_link_evidence),
            "engine": self.engine,
            "stage": self.stage,
            "score_source": self.score_source,
            "selected_docking_attempt_count": self.selected_docking_attempt_count,
        }


@dataclass
class PairCandidate:
    run_id: str
    pdb_id: str
    variant: str
    ph_label: str
    ligand_id: str
    is_control: bool = False
    is_decoy: bool = False
    pair_cell_id: int | None = None
    receptor_context_id: int | None = None
    ligand_pk: int | None = None
    center: tuple[float, float, float] | None = None
    box: tuple[float, float, float] | None = None
    selected_result: SelectedResultRef = field(default_factory=SelectedResultRef)
    source_blockers: list[str] = field(default_factory=list)
    selected_result_blockers: list[str] = field(default_factory=list)
    artifacts: dict[str, list[ArtifactRef]] = field(default_factory=dict)
    source_row: Mapping[str, Any] = field(default_factory=dict)

    @property
    def key(self) -> tuple[str, str, str, str, str]:
        return (self.run_id, self.pdb_id, self.variant, self.ph_label, self.ligand_id)

    @property
    def context_key(self) -> tuple[str, str, str]:
        return (self.pdb_id, self.variant, self.ph_label)

    def key_record(self) -> dict[str, str]:
        return {
            "run_id": self.run_id,
            "pdb_id": self.pdb_id,
            "variant": self.variant,
            "ph_label": self.ph_label,
            "ligand_id": self.ligand_id,
        }


@dataclass(frozen=True)
class PlanRequest:
    run_id: str
    stage: str
    receptors: tuple[str, ...]
    all_receptors: bool
    ligands: tuple[str, ...]
    all_ligands: bool
    variants: tuple[str, ...] = ()
    ph_labels: tuple[str, ...] = ()
    all_contexts: bool = False
    fraction: float | None = 1.0
    count: int | None = None
    seed: int = 0
    cpus: int = 1
    source_stage: str = ""
    selection_strategy: str = ""


def _parse_vector(raw: Any) -> tuple[float, float, float] | None:
    if isinstance(raw, str):
        try:
            raw = json.loads(raw)
        except (TypeError, ValueError, json.JSONDecodeError):
            return None
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or len(raw) != 3:
        return None
    values: list[float] = []
    for item in raw:
        try:
            value = float(item)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(value):
            return None
        values.append(value)
    return (values[0], values[1], values[2])


def _row_vector(
    row: Mapping[str, Any], prefix: str
) -> tuple[float, float, float] | None:
    values: list[float] = []
    for axis in ("x", "y", "z"):
        try:
            value = float(row[f"{prefix}_{axis}"])
        except (KeyError, TypeError, ValueError):
            return None
        if not math.isfinite(value):
            return None
        values.append(value)
    return (values[0], values[1], values[2])


def _validate_request(request: PlanRequest) -> None:
    if not _RUN_ID.fullmatch(request.run_id) or request.run_id in {".", ".."}:
        raise StageRepeatPlanError("run_id contains unsafe path characters")
    if request.stage not in SUPPORTED_STAGES:
        raise StageRepeatPlanError(
            f"stage must be one of: {', '.join(SUPPORTED_STAGES)}"
        )
    if request.all_receptors == bool(request.receptors):
        raise StageRepeatPlanError(
            "choose exactly one of explicit receptors or all_receptors"
        )
    if request.all_ligands == bool(request.ligands):
        raise StageRepeatPlanError(
            "choose exactly one of explicit ligands or all_ligands"
        )
    if request.fraction is not None and request.count is not None:
        raise StageRepeatPlanError("fraction and count are mutually exclusive")
    if request.fraction is None and request.count is None:
        raise StageRepeatPlanError("fraction or count must be recorded")
    if request.fraction is not None and not (0.0 < request.fraction <= 1.0):
        raise StageRepeatPlanError("fraction must be greater than 0 and at most 1")
    if request.count is not None and request.count < 1:
        raise StageRepeatPlanError("count must be at least 1")
    if request.cpus < 1 or request.cpus > MAX_STAGE_CPUS:
        raise StageRepeatPlanError(f"cpus must be between 1 and {MAX_STAGE_CPUS}")
    if request.stage in {"vina", "gnina"} and not request.source_stage:
        raise StageRepeatPlanError(
            f"--source-stage is required for {request.stage} so the output stage is exact"
        )
    strategy = _text(request.selection_strategy).lower()
    if strategy not in {"", "hash", "top-score"}:
        raise StageRepeatPlanError("selection_strategy must be hash or top-score")
    if strategy == "top-score":
        raise StageRepeatPlanError(
            "top-score selection is intentionally unavailable in planner v1; "
            "score field, direction, provenance, and tie policy require an explicit contract"
        )
    partial = request.count is not None or (
        request.fraction is not None and request.fraction < 1.0
    )
    if partial and not strategy:
        raise StageRepeatPlanError(
            "partial plans require --selection-strategy; use hash only for an "
            "explicitly unbiased reproducible subset"
        )


def _candidate_from_csv(row: Mapping[str, Any], requested_run_id: str) -> PairCandidate:
    run_id = _text(row.get("run_id"))
    blockers: list[str] = []
    if not run_id:
        blockers.append("source_run_id_missing")
        run_id = requested_run_id
    elif run_id != requested_run_id:
        blockers.append("source_run_id_mismatch")
    pdb_id = _text(row.get("pdb_id")).upper()
    ligand_id = _text(
        row.get("ligand_base")
        or row.get("ligand_canonical_id")
        or row.get("ligand_id")
        or row.get("ligand")
    )
    if not pdb_id:
        blockers.append("pdb_id_missing")
    if not ligand_id:
        blockers.append("ligand_id_missing")
    result_sha = _sha256_bytes(_json(row).encode("utf-8"))
    pose_stage = next(
        (_text(row.get(key)) for key in _POSE_STAGE_FIELDS if _text(row.get(key))), ""
    )
    selected = SelectedResultRef(
        result_sha256=result_sha,
        engine=_text(row.get("selected_pose_engine") or row.get("final_score_engine")),
        stage=pose_stage,
        score_source=_text(
            row.get("final_score_source") or row.get("z_selected_source")
        ),
        result_json=dict(row),
    )
    return PairCandidate(
        run_id=run_id,
        pdb_id=pdb_id,
        variant=_normalize_variant(row.get("variant")),
        ph_label=_normalize_ph(row.get("ph_label")),
        ligand_id=ligand_id,
        is_control=_text(row.get("is_control")) == "1",
        is_decoy=_text(row.get("is_decoy")) == "1",
        center=_row_vector(row, "center"),
        box=_row_vector(row, "box"),
        selected_result=selected,
        source_blockers=blockers,
        source_row=dict(row),
    )


def load_csv_candidates(path: Path, run_id: str) -> list[PairCandidate]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise StageRepeatPlanError("master row CSV has no header")
            candidates: list[PairCandidate] = []
            for row in reader:
                candidates.append(_candidate_from_csv(row, run_id))
                if len(candidates) > MAX_IN_MEMORY_PLAN_PAIRS:
                    raise StageRepeatPlanError(
                        "stage-repeat source exceeds the bounded planner cap "
                        f"({MAX_IN_MEMORY_PLAN_PAIRS} pairs); "
                        "streaming_stage_repeat_plan_required"
                    )
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise StageRepeatPlanError(f"cannot read master row CSV: {exc}") from exc
    if not candidates:
        raise StageRepeatPlanError("master row CSV contains no pair rows")
    return candidates


def _selected_result_from_sql(row: sqlite3.Row) -> tuple[SelectedResultRef, list[str]]:
    blockers: list[str] = []
    selected_count = int(row["selected_result_count"] or 0)
    if selected_count == 0:
        blockers.append("selected_result_missing")
    elif selected_count != 1:
        blockers.append("selected_result_ambiguous")
    completion_id = row["completion_record_id"]
    if selected_count == 1 and completion_id is None:
        blockers.append("selected_result_completion_link_missing")
    link_evidence = _parse_json_object(row["completion_link_evidence_json"])
    if selected_count == 1 and completion_id is not None:
        if not _text(row["completion_link_method"]) or not link_evidence:
            blockers.append("selected_result_completion_lineage_incomplete")
        selected_docking_count = int(row["selected_docking_attempt_count"] or 0)
        if selected_docking_count == 0:
            blockers.append("matching_selected_docking_attempt_missing")
        elif selected_docking_count != 1:
            blockers.append("matching_selected_docking_attempt_ambiguous")
    return (
        SelectedResultRef(
            result_attempt_id=(
                int(row["result_attempt_id"])
                if row["result_attempt_id"] is not None
                else None
            ),
            result_sha256=_text(row["result_sha256"]),
            completion_record_id=(
                int(completion_id) if completion_id is not None else None
            ),
            completion_sha256=_text(row["completion_sha256"]),
            completion_link_method=_text(row["completion_link_method"]),
            completion_link_evidence=link_evidence,
            engine=_text(row["engine"]),
            stage=_text(row["stage"]),
            score_source=_text(
                row["final_score_source_reconstructed"]
                or row["final_score_source"]
                or row["atlas_score_source"]
            ),
            selected_docking_attempt_count=int(
                row["selected_docking_attempt_count"] or 0
            ),
            result_json=_parse_json_object(row["result_json"]),
        ),
        blockers,
    )


def load_sqlite_candidates(path: Path, run_id: str) -> list[PairCandidate]:
    try:
        connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
        connection.row_factory = sqlite3.Row
    except sqlite3.Error as exc:
        raise StageRepeatPlanError(
            f"cannot open Atlas SQLite read-only: {exc}"
        ) from exc
    try:
        version = int(connection.execute("PRAGMA user_version").fetchone()[0])
        if version < 5:
            raise StageRepeatPlanError(
                f"Atlas SQLite schema {version} lacks selected-attempt lineage; expected >= 5"
            )
        pair_count = int(
            connection.execute(
                """SELECT COUNT(*) FROM pair_cells pc
                   JOIN receptor_contexts rc
                     ON rc.receptor_context_id=pc.receptor_context_id
                   WHERE rc.run_id=?""",
                (run_id,),
            ).fetchone()[0]
        )
        if pair_count > MAX_IN_MEMORY_PLAN_PAIRS:
            raise StageRepeatPlanError(
                "stage-repeat source exceeds the bounded planner cap "
                f"({MAX_IN_MEMORY_PLAN_PAIRS} pairs); "
                "streaming_stage_repeat_plan_required"
            )
        result_columns = {
            _text(info[1])
            for info in connection.execute("PRAGMA table_info(result_attempts)")
        }
        reconstructed_projection = (
            "ra.final_score_source_reconstructed"
            if "final_score_source_reconstructed" in result_columns
            else "NULL AS final_score_source_reconstructed"
        )
        rows = connection.execute(
            f"""
            SELECT pc.pair_cell_id, pc.receptor_context_id, pc.ligand_id,
                   r.run_id, rc.pdb_id, rc.variant, rc.ph_label,
                   rc.center_json, rc.box_json, l.canonical_id,
                   pc.is_control, pc.is_decoy,
                   (SELECT COUNT(*) FROM result_attempts rx
                    WHERE rx.pair_cell_id=pc.pair_cell_id
                      AND rx.selected_for_release=1) AS selected_result_count,
                   ra.result_attempt_id, ra.result_sha256,
                   ra.completion_record_id, ra.completion_link_method,
                   ra.completion_link_evidence_json, ra.result_json,
                   ra.final_score_source, ra.atlas_score_source,
                   {reconstructed_projection},
                   cr.completion_sha256, cr.engine, cr.stage,
                   (SELECT COUNT(*) FROM docking_attempts dx
                    WHERE dx.pair_cell_id=pc.pair_cell_id
                      AND dx.selected_for_release=1
                      AND dx.completion_record_id=ra.completion_record_id)
                       AS selected_docking_attempt_count
            FROM pair_cells pc
            JOIN receptor_contexts rc
              ON rc.receptor_context_id=pc.receptor_context_id
            JOIN runs r ON r.run_id=rc.run_id
            JOIN ligands l ON l.ligand_id=pc.ligand_id
            LEFT JOIN result_attempts ra
              ON ra.result_attempt_id=(
                  SELECT MIN(ry.result_attempt_id) FROM result_attempts ry
                  WHERE ry.pair_cell_id=pc.pair_cell_id
                    AND ry.selected_for_release=1
              )
            LEFT JOIN completion_records cr
              ON cr.completion_record_id=ra.completion_record_id
            WHERE r.run_id=?
            ORDER BY rc.pdb_id, rc.variant, rc.ph_label, l.canonical_id
            """,
            (run_id,),
        )
        candidates: list[PairCandidate] = []
        for row in rows:
            selected, blockers = _selected_result_from_sql(row)
            candidates.append(
                PairCandidate(
                    run_id=_text(row["run_id"]),
                    pdb_id=_text(row["pdb_id"]).upper(),
                    variant=_normalize_variant(row["variant"]),
                    ph_label=_normalize_ph(row["ph_label"]),
                    ligand_id=_text(row["canonical_id"]),
                    is_control=bool(row["is_control"]),
                    is_decoy=bool(row["is_decoy"]),
                    pair_cell_id=int(row["pair_cell_id"]),
                    receptor_context_id=int(row["receptor_context_id"]),
                    ligand_pk=int(row["ligand_id"]),
                    center=_parse_vector(row["center_json"]),
                    box=_parse_vector(row["box_json"]),
                    selected_result=selected,
                    selected_result_blockers=blockers,
                )
            )
    except sqlite3.Error as exc:
        raise StageRepeatPlanError(f"cannot query Atlas SQLite: {exc}") from exc
    finally:
        connection.close()
    if not candidates:
        raise StageRepeatPlanError(f"Atlas SQLite has no pair cells for run {run_id!r}")
    return candidates


def _artifact_from_sql(row: sqlite3.Row) -> ArtifactRef:
    return ArtifactRef(
        role=_text(row["artifact_role"]),
        stage=_text(row["stage"]),
        mode=_text(row["mode"]),
        original_path=_text(row["original_path"]),
        archive_path=_text(row["archive_path"]),
        member_name=_text(row["member_name"]),
        sha256=_text(row["sha256"]),
        verified=bool(row["verified"]),
        metadata=_parse_json_object(row["artifact_json"]),
    )


def enrich_sqlite_artifacts(path: Path, candidates: Sequence[PairCandidate]) -> None:
    """Attach typed artifacts only for already-selected plan rows."""
    if not candidates:
        return
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        for offset in range(0, len(candidates), 300):
            chunk = candidates[offset : offset + 300]
            pair_map: dict[int, list[PairCandidate]] = {}
            context_map: dict[int, list[PairCandidate]] = {}
            ligand_map: dict[int, list[PairCandidate]] = {}
            for candidate in chunk:
                if candidate.pair_cell_id is not None:
                    pair_map.setdefault(candidate.pair_cell_id, []).append(candidate)
                if candidate.receptor_context_id is not None:
                    context_map.setdefault(candidate.receptor_context_id, []).append(
                        candidate
                    )
                if candidate.ligand_pk is not None:
                    ligand_map.setdefault(candidate.ligand_pk, []).append(candidate)
            if not pair_map:
                continue
            pair_marks = ",".join("?" for _ in pair_map)
            context_marks = ",".join("?" for _ in context_map) or "NULL"
            ligand_marks = ",".join("?" for _ in ligand_map) or "NULL"
            query = f"""
                SELECT artifact_id, pair_cell_id, receptor_context_id, ligand_id,
                       artifact_role, artifact_scope, stage, mode, original_path,
                       archive_path, member_name, sha256, verified, artifact_json
                FROM artifacts
                WHERE pair_cell_id IN ({pair_marks})
                   OR (artifact_scope='receptor' AND receptor_context_id IN ({context_marks}))
                   OR (artifact_scope='ligand' AND ligand_id IN ({ligand_marks}))
                ORDER BY artifact_id
            """
            params = [*pair_map, *context_map, *ligand_map]
            for row in connection.execute(query, params):
                scope = _text(row["artifact_scope"])
                if scope == "pair":
                    targets = pair_map.get(int(row["pair_cell_id"]), [])
                elif scope == "receptor":
                    targets = context_map.get(int(row["receptor_context_id"]), [])
                elif scope == "ligand":
                    targets = ligand_map.get(int(row["ligand_id"]), [])
                else:
                    targets = []
                artifact = _artifact_from_sql(row)
                for candidate in targets:
                    candidate.artifacts.setdefault(artifact.role, []).append(artifact)
    except sqlite3.Error as exc:
        raise StageRepeatPlanError(f"cannot load typed artifacts: {exc}") from exc
    finally:
        connection.close()


def _resolve_csv_path(repo_root: Path, raw: Any) -> Path | None:
    token = _text(raw)
    if not token:
        return None
    path = Path(token).expanduser()
    return path if path.is_absolute() else repo_root / path


def _csv_artifact(
    repo_root: Path,
    *,
    role: str,
    path_value: Any,
    stage: str = "",
    declared_sha256: Any = "",
    metadata: Mapping[str, Any] | None = None,
) -> ArtifactRef | None:
    path = _resolve_csv_path(repo_root, path_value)
    if path is None:
        return None
    digest = _text(declared_sha256).lower()
    exists = path.is_file()
    if exists:
        actual = _sha256_file(path)
        verified = (not digest or digest == actual) and bool(_SHA256.fullmatch(actual))
        digest = actual
    else:
        verified = False
    return ArtifactRef(
        role=role,
        stage=stage,
        original_path=str(path),
        sha256=digest,
        verified=verified,
        metadata=dict(metadata or {}),
    )


def enrich_csv_artifacts(repo_root: Path, candidates: Sequence[PairCandidate]) -> None:
    for candidate in candidates:
        row = candidate.source_row
        receptor = _csv_artifact(
            repo_root,
            role="prepared_receptor",
            path_value=row.get("prepared_receptor_path"),
            declared_sha256=row.get("prepared_receptor_sha256"),
        )
        ligand = _csv_artifact(
            repo_root,
            role="prepared_ligand",
            path_value=row.get("prepared_ligand_path"),
            declared_sha256=row.get("prepared_ligand_sha256"),
        )
        pose_values = _dedupe(_text(row.get(key)) for key in _POSE_PATH_FIELDS)
        pose_stage = next(
            (_text(row.get(key)) for key in _POSE_STAGE_FIELDS if _text(row.get(key))),
            "",
        )
        pose = None
        if len(pose_values) == 1:
            pose = _csv_artifact(
                repo_root,
                role="docking_pose",
                path_value=pose_values[0],
                stage=pose_stage,
                declared_sha256=row.get("selected_pose_sha256"),
                metadata={
                    "result_sha256": candidate.selected_result.result_sha256,
                    "link_method": "explicit_selected_pose_field_in_result_row",
                },
            )
        elif len(pose_values) > 1:
            candidate.source_blockers.append("selected_pose_fields_conflict")
        for artifact in (receptor, ligand, pose):
            if artifact is not None:
                candidate.artifacts.setdefault(artifact.role, []).append(artifact)


def _scope_candidates(
    candidates: Sequence[PairCandidate], request: PlanRequest
) -> list[PairCandidate]:
    receptors = {token.upper() for token in request.receptors}
    ligands = set(request.ligands)
    variants = {_normalize_variant(token) for token in request.variants}
    ph_labels = {_normalize_ph(token) for token in request.ph_labels}

    available_receptors = {item.pdb_id for item in candidates}
    available_ligands = {item.ligand_id for item in candidates}
    missing_receptors = sorted(receptors - available_receptors)
    missing_ligands = sorted(ligands - available_ligands)
    if missing_receptors:
        raise StageRepeatPlanError(
            f"requested receptors absent from source: {', '.join(missing_receptors)}"
        )
    if missing_ligands:
        raise StageRepeatPlanError(
            f"requested ligands absent from source: {', '.join(missing_ligands)}"
        )

    scoped = [
        item
        for item in candidates
        if (request.all_receptors or item.pdb_id in receptors)
        and (request.all_ligands or item.ligand_id in ligands)
        and (not variants or item.variant in variants)
        and (not ph_labels or item.ph_label in ph_labels)
    ]
    if not scoped:
        raise StageRepeatPlanError("scope filters select no receptor-ligand pairs")

    if not request.all_contexts and not variants and not ph_labels:
        contexts_by_receptor: dict[str, set[tuple[str, str]]] = {}
        for item in scoped:
            contexts_by_receptor.setdefault(item.pdb_id, set()).add(
                (item.variant, item.ph_label)
            )
        ambiguous = sorted(
            pdb_id
            for pdb_id, contexts in contexts_by_receptor.items()
            if len(contexts) > 1
        )
        if ambiguous:
            raise StageRepeatPlanError(
                "multiple receptor contexts are in scope; pass --variant/--ph or "
                f"--all-contexts explicitly (receptors: {', '.join(ambiguous[:10])})"
            )

    seen: dict[tuple[str, str, str, str, str], int] = {}
    duplicates: list[tuple[str, str, str, str, str]] = []
    for item in scoped:
        seen[item.key] = seen.get(item.key, 0) + 1
        if seen[item.key] == 2:
            duplicates.append(item.key)
    if duplicates:
        examples = ", ".join("/".join(key[1:]) for key in duplicates[:5])
        raise StageRepeatPlanError(
            f"source has duplicate rows for exact pair contexts: {examples}"
        )
    return scoped


def _selection_digest(item: PairCandidate, seed: int) -> str:
    payload = "\0".join((str(seed), *item.key))
    return _sha256_bytes(payload.encode("utf-8"))


def _select_candidates(
    candidates: Sequence[PairCandidate], request: PlanRequest
) -> list[tuple[str, PairCandidate]]:
    ordered = sorted(
        ((_selection_digest(item, request.seed), item) for item in candidates),
        key=lambda pair: (pair[0], pair[1].key),
    )
    if request.count is not None:
        if request.count > len(ordered):
            raise StageRepeatPlanError(
                f"count {request.count} exceeds scoped pair count {len(ordered)}"
            )
        count = request.count
    else:
        fraction = float(request.fraction or 0.0)
        count = max(1, int(math.ceil(len(ordered) * fraction)))
    return ordered[:count]


def _artifact_identity(artifact: ArtifactRef) -> tuple[str, ...]:
    return (
        artifact.role,
        _normalize_stage(artifact.stage),
        _normalize_locator(artifact.original_path),
        _normalize_locator(artifact.archive_path),
        _normalize_locator(artifact.member_name),
        artifact.sha256.lower(),
    )


def _usable_artifacts(candidate: PairCandidate, role: str) -> list[ArtifactRef]:
    unique: dict[tuple[str, ...], ArtifactRef] = {}
    for artifact in candidate.artifacts.get(role, []):
        if not artifact.verified or not _SHA256.fullmatch(artifact.sha256):
            continue
        if not artifact.original_path and not (
            artifact.archive_path and artifact.member_name
        ):
            continue
        unique[_artifact_identity(artifact)] = artifact
    return list(unique.values())


def _metadata_values(metadata: Mapping[str, Any], key: str) -> list[Any]:
    values: list[Any] = []
    if key in metadata:
        values.append(metadata[key])
    for container_key in ("provenance", "selected_result", "lineage"):
        nested = metadata.get(container_key)
        if isinstance(nested, Mapping) and key in nested:
            values.append(nested[key])
    return values


def _artifact_linked_to_selected_result(
    artifact: ArtifactRef, selected: SelectedResultRef
) -> bool:
    metadata = artifact.metadata
    exact_checks: tuple[tuple[str, Any], ...] = (
        ("result_attempt_id", selected.result_attempt_id),
        ("result_sha256", selected.result_sha256),
        ("selected_result_sha256", selected.result_sha256),
        ("completion_record_id", selected.completion_record_id),
        ("completion_sha256", selected.completion_sha256),
    )
    for key, expected in exact_checks:
        if expected in (None, ""):
            continue
        if any(
            _text(value) == _text(expected) for value in _metadata_values(metadata, key)
        ):
            return True

    result_paths = {
        _normalize_locator(selected.result_json.get(key))
        for key in _POSE_PATH_FIELDS
        if _text(selected.result_json.get(key))
    }
    return bool(result_paths & artifact.locator_tokens())


def _choose_single_artifact(
    candidate: PairCandidate,
    role: str,
    blockers: list[str],
) -> ArtifactRef | None:
    usable = _usable_artifacts(candidate, role)
    if not usable:
        blockers.append(f"{role}_missing_or_unverified")
        return None
    if len(usable) != 1:
        blockers.append(f"{role}_ambiguous")
        return None
    return usable[0]


def _score_source_requires_contributing_pose_lineage(
    selected: SelectedResultRef,
) -> bool | None:
    source = _text(selected.score_source).lower()
    if not source:
        source = _text(
            selected.result_json.get("final_score_source")
            or selected.result_json.get("z_selected_source")
        ).lower()
    requirement = score_source_pose_linkage_requirement(source)
    if requirement is None:
        return None
    return any(
        token in requirement
        for token in ("no_unique", "no_immutable", "component_selected")
    )


def _declared_contributing_pose_hashes(selected: SelectedResultRef) -> set[str]:
    for key in ("contributing_pose_sha256s", "selected_score_pose_sha256s"):
        raw = selected.result_json.get(key)
        if not isinstance(raw, list):
            continue
        values = {_text(value).lower() for value in raw}
        if values and all(_SHA256.fullmatch(value) for value in values):
            return values
    return set()


def _score_pose_role(artifact: ArtifactRef) -> str:
    for key in ("score_pose_role", "pose_score_role"):
        values = _metadata_values(artifact.metadata, key)
        if values:
            return _text(values[0]).lower()
    return ""


def _choose_selected_poses(
    candidate: PairCandidate,
    request: PlanRequest,
    blockers: list[str],
) -> list[ArtifactRef]:
    poses = _usable_artifacts(candidate, "docking_pose")
    if not poses:
        blockers.append("selected_docking_pose_missing_or_unverified")
        return []

    contributing_lineage = _score_source_requires_contributing_pose_lineage(
        candidate.selected_result
    )
    if contributing_lineage is None:
        blockers.append("score_source_pose_linkage_unknown")
        return []
    if contributing_lineage:
        declared = _declared_contributing_pose_hashes(candidate.selected_result)
        by_hash: dict[str, list[ArtifactRef]] = {}
        for pose in poses:
            if not pose.stage:
                continue
            if request.source_stage and _normalize_stage(
                pose.stage
            ) != _normalize_stage(request.source_stage):
                continue
            if not _artifact_linked_to_selected_result(pose, candidate.selected_result):
                continue
            if _score_pose_role(pose) != "contributing_pose":
                continue
            by_hash.setdefault(pose.sha256.lower(), []).append(pose)
        if (
            not declared
            or set(by_hash) != declared
            or any(len(items) != 1 for items in by_hash.values())
        ):
            blockers.append("score_source_pose_identity_unresolved")
            return []
        return sorted((items[0] for items in by_hash.values()), key=_artifact_identity)

    expected_stage = _normalize_stage(
        request.source_stage or candidate.selected_result.stage
    )
    if not expected_stage:
        blockers.append("selected_pose_stage_missing")
        return []
    stage_matches = [
        pose for pose in poses if _normalize_stage(pose.stage) == expected_stage
    ]
    if not stage_matches:
        blockers.append("selected_docking_pose_stage_mismatch")
        return []
    linked = [
        pose
        for pose in stage_matches
        if _artifact_linked_to_selected_result(pose, candidate.selected_result)
    ]
    if not linked:
        blockers.append("pose_artifact_not_linked_to_selected_result")
        return []
    unique = {_artifact_identity(pose): pose for pose in linked}
    if len(unique) != 1:
        blockers.append("selected_docking_pose_ambiguous")
        return []
    return list(unique.values())


def _evaluate_candidate(
    candidate: PairCandidate, request: PlanRequest, selection_digest: str
) -> dict[str, Any]:
    blockers = list(dict.fromkeys(candidate.source_blockers))
    artifacts: dict[str, Any] = {}
    selected_poses: list[ArtifactRef] = []
    receptor = _choose_single_artifact(candidate, "prepared_receptor", blockers)
    if receptor is not None:
        artifacts["prepared_receptor"] = receptor.public_record()

    needs_ligand = request.stage in {
        "vina",
        "gnina",
        "scorch",
        "mmgbsa",
        "pose-validation",
    }
    if needs_ligand:
        ligand = _choose_single_artifact(candidate, "prepared_ligand", blockers)
        if ligand is not None:
            artifacts["prepared_ligand"] = ligand.public_record()

    needs_pose = request.stage in {"scorch", "mmgbsa", "pose-validation"}
    if needs_pose:
        blockers.extend(candidate.selected_result_blockers)
        selected_poses = _choose_selected_poses(candidate, request, blockers)
        if (
            selected_poses
            and _score_source_requires_contributing_pose_lineage(
                candidate.selected_result
            )
            is True
        ):
            artifacts["contributing_docking_poses"] = [
                pose.public_record() for pose in selected_poses
            ]
        elif len(selected_poses) == 1:
            artifacts["selected_docking_pose"] = selected_poses[0].public_record()

    if request.stage in {"vina", "gnina"}:
        if candidate.center is None:
            blockers.append("search_box_center_missing")
        if candidate.box is None:
            blockers.append("search_box_size_missing")

    if request.stage in {"scorch", "mmgbsa", "pose-validation"}:
        selected = candidate.selected_result
        if not selected.result_sha256:
            blockers.append("selected_result_hash_missing")
        if candidate.pair_cell_id is not None:
            # Database-backed plans require causal completion lineage.  CSV rows
            # instead prove association through an explicit selected-pose field.
            if selected.completion_record_id is None:
                blockers.append("selected_result_completion_link_missing")
            if (
                not selected.completion_link_method
                or not selected.completion_link_evidence
            ):
                blockers.append("selected_result_completion_lineage_incomplete")
            if selected.selected_docking_attempt_count != 1:
                blockers.append("matching_selected_docking_attempt_not_unique")

    blockers = list(dict.fromkeys(blockers))
    source_stages = sorted(
        {_normalize_stage(pose.stage) for pose in selected_poses if pose.stage}
    )
    if request.source_stage:
        source_stage = request.source_stage
    elif len(source_stages) == 1:
        source_stage = source_stages[0]
    elif source_stages:
        source_stage = ""
    else:
        source_stage = candidate.selected_result.stage
    row_core = {
        **candidate.key_record(),
        "selection_digest": selection_digest,
        "stage": request.stage,
        "source_stage": source_stage,
        "source_stages": source_stages,
    }
    plan_row_id = _sha256_bytes(_json(row_core).encode("utf-8"))
    return {
        "schema": PAIR_SCHEMA,
        "plan_row_id": plan_row_id,
        **row_core,
        "is_control": candidate.is_control,
        "is_decoy": candidate.is_decoy,
        "status": "ready" if not blockers else "blocked",
        "blockers": blockers,
        "selected_result": candidate.selected_result.public_record(),
        "search_box": {
            "center": list(candidate.center) if candidate.center is not None else None,
            "size": list(candidate.box) if candidate.box is not None else None,
        },
        "artifacts": artifacts,
    }


_STAGE_ENTRYPOINTS: dict[str, dict[str, str]] = {
    "vina": {
        "kind": "python_callable",
        "entrypoint": "docking.docking_stage_runner:run_one_stage",
    },
    "gnina": {
        "kind": "python_callable",
        "entrypoint": "docking.docking_gnina:run_gnina_for_stage",
    },
    "scorch": {
        "kind": "python_module",
        "entrypoint": "post_docking.rescoring.rescoring_scorch",
    },
    "mmgbsa": {
        "kind": "python_callable",
        "entrypoint": "post_docking.mmgbsa.mmgbsa_pipeline:_maybe_run_mmgbsa_for_pdb",
    },
    "pose-validation": {
        "kind": "python_script",
        "entrypoint": "tools/pose_bust.py",
    },
}


def _executor_status(stage: str) -> dict[str, Any]:
    gap = "exact_pair_plan_executor_not_implemented"
    if stage == "pose-validation":
        gap = "pose_validator_exact_pair_allowlist_not_implemented"
    return {
        "executor_status": "planner_only",
        "can_execute": False,
        "exact_pair_allowlist_supported": False,
        "existing_stage_entrypoint": _STAGE_ENTRYPOINTS[stage],
        "adapter": None,
        "command": None,
        "gap": gap,
        "source_run_mutation_allowed": False,
    }


def _jsonl_payloads(rows: Sequence[Mapping[str, Any]]) -> list[bytes]:
    payloads: list[bytes] = []
    for offset in range(0, len(rows), PAIR_LEDGER_SHARD_ROWS):
        chunk = rows[offset : offset + PAIR_LEDGER_SHARD_ROWS]
        payloads.append("".join(f"{_json(row)}\n" for row in chunk).encode("utf-8"))
    return payloads


def _plan_core(
    request: PlanRequest,
    *,
    source_kind: str,
    source_sha256: str,
    source_path: Path,
    scoped_count: int,
    selected_count: int,
    summary: Mapping[str, int],
    pair_index_sha256: str,
    ledger_payloads: Sequence[bytes],
) -> dict[str, Any]:
    partial = request.count is not None or (
        request.fraction is not None and request.fraction < 1.0
    )
    ledger_shards = [
        {
            "index": index,
            "row_count": min(
                PAIR_LEDGER_SHARD_ROWS,
                selected_count - ((index - 1) * PAIR_LEDGER_SHARD_ROWS),
            ),
            "sha256": _sha256_bytes(payload),
        }
        for index, payload in enumerate(ledger_payloads, start=1)
    ]
    return {
        "schema": PLAN_SCHEMA,
        "run_id": request.run_id,
        "stage": request.stage,
        "source_stage": request.source_stage,
        "source": {
            "kind": source_kind,
            "path": str(source_path.resolve()),
            "sha256": source_sha256,
        },
        "scope": {
            "all_receptors": request.all_receptors,
            "receptors": list(request.receptors),
            "all_ligands": request.all_ligands,
            "ligands": list(request.ligands),
            "all_contexts": request.all_contexts,
            "variants": list(request.variants),
            "ph_labels": list(request.ph_labels),
            "scoped_pair_count": scoped_count,
        },
        "selection": {
            "strategy": _text(request.selection_strategy).lower() if partial else "all",
            "algorithm": SELECTION_ALGORITHM,
            "seed": request.seed,
            "fraction": request.fraction,
            "count": request.count,
            "unit": "global_pair",
            "selected_count": selected_count,
        },
        "resources": {
            "cpus": request.cpus,
            "hard_max_cpus": MAX_STAGE_CPUS,
            "in_memory_pair_cap": MAX_IN_MEMORY_PLAN_PAIRS,
        },
        "scalability": {
            "full_matrix_supported": False,
            "reason": "bounded_planner_v1_requires_streaming_work_database",
            "hard_pair_cap": MAX_IN_MEMORY_PLAN_PAIRS,
            "json_embeds_pair_rows": False,
            "ledger_shard_rows": PAIR_LEDGER_SHARD_ROWS,
        },
        "dry_run": True,
        "execution": {
            "status": "not_executed",
            "planner_only": True,
            "mutates_source_run": False,
            "launch_command_recorded": False,
            **_executor_status(request.stage),
        },
        "pair_index": {
            "format": "csv",
            "sha256": pair_index_sha256,
            "row_count": selected_count,
        },
        "pair_ledger": {
            "schema": PAIR_SCHEMA,
            "format": "jsonl",
            "row_count": selected_count,
            "shards": ledger_shards,
        },
        "summary": dict(summary),
    }


def _pairs_csv_bytes(rows: Sequence[Mapping[str, Any]]) -> bytes:
    fields = (
        "plan_row_id",
        "run_id",
        "pdb_id",
        "variant",
        "ph_label",
        "ligand_id",
        "stage",
        "source_stage",
        "selection_digest",
        "status",
        "blockers_json",
    )
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(
            {
                key: (
                    _json(row.get("blockers", []))
                    if key == "blockers_json"
                    else row.get(key, "")
                )
                for key in fields
            }
        )
    return buffer.getvalue().encode("utf-8")


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(payload)
    os.replace(temporary, path)


def _materialized_plan(
    core: Mapping[str, Any],
    *,
    plan_id: str,
    pairs_path: Path,
    ledger_paths: Sequence[Path],
    created_at: str,
) -> dict[str, Any]:
    pair_ledger = {
        **core["pair_ledger"],
        "shards": [
            {**record, "path": path.name}
            for record, path in zip(
                core["pair_ledger"]["shards"], ledger_paths, strict=True
            )
        ],
    }
    return {
        **core,
        "plan_id": plan_id,
        "created_at": created_at,
        "pair_index": {**core["pair_index"], "path": pairs_path.name},
        "pair_ledger": pair_ledger,
    }


def _reuse_existing_plan(
    *,
    plan_path: Path,
    pairs_path: Path,
    ledger_paths: Sequence[Path],
    pairs_payload: bytes,
    ledger_payloads: Sequence[bytes],
    core: Mapping[str, Any],
    plan_id: str,
) -> dict[str, Any] | None:
    """Return an identical persisted plan without rewriting any output bytes."""

    expected_paths = [plan_path, pairs_path, *ledger_paths]
    expected_ledger_paths = set(ledger_paths)
    unexpected_ledgers = sorted(
        path
        for path in plan_path.parent.glob(f"{plan_path.stem}.pairs-*.jsonl")
        if path not in expected_ledger_paths
    )
    if unexpected_ledgers:
        names = ", ".join(path.name for path in unexpected_ledgers[:5])
        raise StageRepeatPlanError(
            f"conflicting stage-repeat ledger shard(s) already exist: {names}"
        )

    present = [path.exists() or path.is_symlink() for path in expected_paths]
    if not any(present):
        return None
    if not all(present):
        missing = ", ".join(
            path.name for path, exists in zip(expected_paths, present, strict=True)
            if not exists
        )
        raise StageRepeatPlanError(
            "incomplete stage-repeat output set already exists; refusing to "
            f"overwrite it (missing: {missing})"
        )
    if any(path.is_symlink() or not path.is_file() for path in expected_paths):
        raise StageRepeatPlanError(
            "stage-repeat output set contains a symlink or non-file; refusing reuse"
        )

    try:
        persisted_payload = plan_path.read_bytes()
        persisted_raw = json.loads(persisted_payload)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise StageRepeatPlanError(
            "existing stage-repeat plan is unreadable; refusing overwrite"
        ) from exc
    if not isinstance(persisted_raw, Mapping):
        raise StageRepeatPlanError(
            "existing stage-repeat plan is not a JSON object; refusing overwrite"
        )
    persisted = dict(persisted_raw)
    created_at = _text(persisted.get("created_at"))
    if not created_at:
        raise StageRepeatPlanError(
            "existing stage-repeat plan lacks created_at; refusing overwrite"
        )
    expected_plan = _materialized_plan(
        core,
        plan_id=plan_id,
        pairs_path=pairs_path,
        ledger_paths=ledger_paths,
        created_at=created_at,
    )
    expected_plan_payload = (
        json.dumps(expected_plan, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    try:
        pair_matches = pairs_path.read_bytes() == pairs_payload
        ledgers_match = all(
            path.read_bytes() == payload
            for path, payload in zip(ledger_paths, ledger_payloads, strict=True)
        )
    except OSError as exc:
        raise StageRepeatPlanError(
            "existing stage-repeat output is unreadable; refusing overwrite"
        ) from exc
    if (
        persisted_payload != expected_plan_payload
        or not pair_matches
        or not ledgers_match
    ):
        raise StageRepeatPlanError(
            "conflicting stage-repeat output already exists; refusing overwrite"
        )
    return {
        "plan": expected_plan,
        "plan_path": plan_path,
        "plan_sha256": _sha256_bytes(persisted_payload),
        "pairs_path": pairs_path,
        "ledger_paths": list(ledger_paths),
        "reused_existing": True,
    }


def build_stage_repeat_plan(
    repo_root: Path,
    request: PlanRequest,
    *,
    database: Path | None = None,
    master_rows: Path | None = None,
    output: Path | None = None,
) -> dict[str, Any]:
    """Build and persist a deterministic dry-run pair plan."""
    _validate_request(request)
    if database is not None and master_rows is not None:
        raise StageRepeatPlanError("database and master_rows are mutually exclusive")
    if database is not None:
        source_path = database.expanduser()
        source_kind = "atlas_sqlite"
        loader = load_sqlite_candidates
    else:
        source_path = (
            master_rows.expanduser()
            if master_rows is not None
            else run_output_dir(repo_root, "data", request.run_id) / "master_rows.csv"
        )
        source_kind = "master_rows_csv"
        loader = load_csv_candidates
    if not source_path.is_absolute():
        source_path = repo_root / source_path
    if not source_path.is_file():
        raise StageRepeatPlanError(f"source does not exist: {source_path}")

    source_sha256 = _sha256_file(source_path)
    candidates = loader(source_path, request.run_id)
    scoped = _scope_candidates(candidates, request)
    selected = _select_candidates(scoped, request)
    selected_candidates = [item for _digest, item in selected]
    if source_kind == "atlas_sqlite":
        enrich_sqlite_artifacts(source_path, selected_candidates)
    else:
        enrich_csv_artifacts(repo_root, selected_candidates)

    rows = [
        _evaluate_candidate(candidate, request, digest)
        for digest, candidate in selected
    ]
    summary = {
        "selected": len(rows),
        "ready": sum(row["status"] == "ready" for row in rows),
        "blocked": sum(row["status"] == "blocked" for row in rows),
    }
    pairs_payload = _pairs_csv_bytes(rows)
    ledger_payloads = _jsonl_payloads(rows)
    core = _plan_core(
        request,
        source_kind=source_kind,
        source_sha256=source_sha256,
        source_path=source_path,
        scoped_count=len(scoped),
        selected_count=len(rows),
        summary=summary,
        pair_index_sha256=_sha256_bytes(pairs_payload),
        ledger_payloads=ledger_payloads,
    )
    plan_id = _sha256_bytes(_json(core).encode("utf-8"))
    if output is None:
        output_dir = (
            run_output_dir(repo_root, "data", request.run_id) / "stage_repeat_plans"
        )
        plan_path = output_dir / f"{request.stage}-{plan_id[:16]}.json"
    else:
        plan_path = output.expanduser()
        if not plan_path.is_absolute():
            plan_path = repo_root / plan_path

    pairs_path = plan_path.with_name(f"{plan_path.stem}.pairs.csv")
    ledger_paths = [
        plan_path.with_name(f"{plan_path.stem}.pairs-{index:05d}.jsonl")
        for index in range(1, len(ledger_payloads) + 1)
    ]
    reused = _reuse_existing_plan(
        plan_path=plan_path,
        pairs_path=pairs_path,
        ledger_paths=ledger_paths,
        pairs_payload=pairs_payload,
        ledger_payloads=ledger_payloads,
        core=core,
        plan_id=plan_id,
    )
    if reused is not None:
        return reused

    plan = _materialized_plan(
        core,
        plan_id=plan_id,
        pairs_path=pairs_path,
        ledger_paths=ledger_paths,
        created_at=datetime.now(timezone.utc).isoformat(),
    )
    plan_payload = (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode("utf-8")
    _atomic_write(pairs_path, pairs_payload)
    for ledger_path, payload in zip(ledger_paths, ledger_payloads, strict=True):
        _atomic_write(ledger_path, payload)
    _atomic_write(plan_path, plan_payload)
    return {
        "plan": plan,
        "plan_path": plan_path,
        "plan_sha256": _sha256_bytes(plan_payload),
        "pairs_path": pairs_path,
        "ledger_paths": ledger_paths,
        "reused_existing": False,
    }


__all__ = [
    "MAX_IN_MEMORY_PLAN_PAIRS",
    "MAX_STAGE_CPUS",
    "PAIR_LEDGER_SHARD_ROWS",
    "PLAN_SCHEMA",
    "PAIR_SCHEMA",
    "SUPPORTED_STAGES",
    "PlanRequest",
    "StageRepeatPlanError",
    "build_stage_repeat_plan",
]
