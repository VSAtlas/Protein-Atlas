"""Build a failure-complete SQLite snapshot from existing Atlas artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import yaml  # type: ignore[import-untyped]

from analysis.atlas_database.annotations import ingest_scientific_annotations
from analysis.atlas_database.manifest import iter_run_entries, load_release_manifest
from analysis.atlas_database.schema import SCHEMA_VERSION, create_schema
from config.output_paths import output_root

_LIGAND_SUFFIXES = (
    ".sanitized",
    "_ctrl_redock_vina",
    "_ctrl_redock",
)
_STAGE_SUFFIX = re.compile(r"_(?:vina_)?stage\d+$", re.IGNORECASE)


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _text(value: Any) -> str:
    return str(value or "").strip()


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _normalize_ph(value: Any) -> str:
    token = _text(value)
    if token.lower() in {"base", "none", "null"} or token.lower().endswith(".csv"):
        return ""
    return token


def _float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return (
        parsed
        if parsed == parsed and parsed not in (float("inf"), float("-inf"))
        else None
    )


def _int(value: Any) -> int | None:
    parsed = _float(value)
    return int(parsed) if parsed is not None else None


def canonical_ligand_id(value: Any) -> str:
    """Normalize a ligand artifact name for joins without discarding identity."""
    name = Path(_text(value)).name
    for suffix in (".pdbqt.gz", ".pdbqt", ".sdf.gz", ".sdf", ".mol2.gz", ".mol2"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = _STAGE_SUFFIX.sub("", name)
    for suffix in _LIGAND_SUFFIXES:
        if name.lower().endswith(suffix.lower()):
            name = name[: -len(suffix)]
    return name.strip()


def _run_paths(repo_root: Path, entry: Mapping[str, Any]) -> dict[str, Path]:
    run_id = _text(entry.get("run_id"))
    overrides = _mapping(entry.get("paths"))
    return {
        "manifest": Path(_text(overrides.get("run_manifest")))
        if overrides.get("run_manifest")
        else output_root(repo_root, "manifests") / run_id / "run_manifest.yaml",
        "master_rows": Path(_text(overrides.get("master_rows")))
        if overrides.get("master_rows")
        else output_root(repo_root, "data") / run_id / "master_rows.csv",
        "docked": Path(_text(overrides.get("docked")))
        if overrides.get("docked")
        else output_root(repo_root, "docked") / run_id,
        "archive_index": Path(_text(overrides.get("archive_index")))
        if overrides.get("archive_index")
        else output_root(repo_root, "docked")
        / run_id
        / "_artifact_archives"
        / "archive_index.json",
    }


def _completion_manifest_paths(docked: Path) -> list[Path]:
    if not docked.is_dir():
        return []
    return sorted(docked.rglob("completion_*.json"))


def _canonical_completion_ids(values: list[str], field: str) -> set[str]:
    canonical_ids: set[str] = set()
    for index, value in enumerate(values):
        canonical = canonical_ligand_id(value)
        if not canonical:
            raise ValueError(f"{field}[{index}] produces an empty canonical ligand ID")
        if canonical in canonical_ids:
            raise ValueError(
                f"{field} contains duplicate canonical ligand ID: {canonical}"
            )
        canonical_ids.add(canonical)
    return canonical_ids


def _load_completion_payload(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot parse JSON ({exc})") from exc
    if not isinstance(payload, Mapping):
        raise ValueError("top-level value must be an object")
    if not _text(payload.get("pdb_id")):
        raise ValueError("pdb_id must be non-empty")
    expected = payload.get("expected_ligands")
    if not isinstance(expected, list):
        raise ValueError("expected_ligands must be a list")
    if any(not isinstance(item, str) or not item.strip() for item in expected):
        raise ValueError("expected_ligands must contain only non-empty strings")
    expected_ids = _canonical_completion_ids(expected, "expected_ligands")
    missing = payload.get("missing_ligands_after", [])
    if not isinstance(missing, list):
        raise ValueError("missing_ligands_after must be a list when present")
    if any(not isinstance(item, str) or not item.strip() for item in missing):
        raise ValueError("missing_ligands_after must contain only non-empty strings")
    missing_ids = _canonical_completion_ids(missing, "missing_ligands_after")
    unexpected_missing = sorted(missing_ids - expected_ids)
    if unexpected_missing:
        raise ValueError(
            "missing_ligands_after contains canonical ligand IDs absent from "
            f"expected_ligands: {', '.join(unexpected_missing)}"
        )
    markers = payload.get("failure_markers", {})
    if not isinstance(markers, Mapping):
        raise ValueError("failure_markers must be an object when present")
    return dict(payload)


def audit_release_inputs(manifest_path: Path, repo_root: Path) -> dict[str, Any]:
    """Audit release inputs without creating a database."""
    manifest = load_release_manifest(manifest_path)
    runs: list[dict[str, Any]] = []
    for entry in iter_run_entries(manifest):
        paths = _run_paths(repo_root, entry)
        manifest_exists = paths["manifest"].is_file()
        master_exists = paths["master_rows"].is_file()
        completion_paths = _completion_manifest_paths(paths["docked"])
        completion_errors: list[dict[str, str]] = []
        for completion_path in completion_paths:
            try:
                _load_completion_payload(completion_path)
            except ValueError as exc:
                completion_errors.append(
                    {"path": str(completion_path), "error": str(exc)}
                )
        completion_valid_count = len(completion_paths) - len(completion_errors)
        runs.append(
            {
                "run_id": entry["run_id"],
                "paths": {name: str(path) for name, path in paths.items()},
                "run_manifest_present": manifest_exists,
                "master_rows_present": master_exists,
                "archive_index_present": paths["archive_index"].is_file(),
                "completion_manifest_count": len(completion_paths),
                "completion_manifest_valid_count": completion_valid_count,
                "completion_manifest_invalid_count": len(completion_errors),
                "completion_manifest_errors": completion_errors,
                "failure_complete_inputs_present": manifest_exists
                and master_exists
                and completion_valid_count > 0
                and not completion_errors,
            }
        )
    errors = [
        f"{run['run_id']}: missing run manifest"
        for run in runs
        if not run["run_manifest_present"]
    ]
    for run in runs:
        for issue in run["completion_manifest_errors"]:
            errors.append(
                f"{run['run_id']}: invalid completion manifest "
                f"{issue['path']}: {issue['error']}"
            )
    return {
        "release_id": manifest["release_id"],
        "manifest_path": str(manifest_path),
        "manifest_sha256": _sha256(manifest_path),
        "runs": runs,
        "run_count": len(runs),
        "failure_complete_run_count": sum(
            1 for run in runs if run["failure_complete_inputs_present"]
        ),
        "errors": errors,
    }


def _load_mapping(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    return loaded if isinstance(loaded, dict) else {}


def _context_values(key: str, entry: Mapping[str, Any]) -> tuple[str, str, str]:
    parts = key.split("|")
    pdb_id = _text(entry.get("pdb_id") or (parts[0] if parts else "")).upper()
    variant = _text(
        entry.get("variant") or (parts[1] if len(parts) > 1 else "")
    ).upper()
    ph = _normalize_ph(
        entry.get("ph") or entry.get("ph_label") or (parts[2] if len(parts) > 2 else "")
    )
    return pdb_id, variant, ph


def _pocket(entry: Mapping[str, Any]) -> tuple[str, Any, Any]:
    stages = _mapping(entry.get("stages"))
    pocket = _mapping(stages.get("pocket_detection"))
    details = _mapping(pocket.get("details"))
    return (
        _text(details.get("method") or details.get("pocket_method")),
        details.get("center"),
        details.get("box_size"),
    )


def _upsert_context(
    connection: sqlite3.Connection,
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
    **values: Any,
) -> int:
    connection.execute(
        """INSERT INTO receptor_contexts
        (run_id, pdb_id, variant, ph_label, library, status, failure_reason,
         input_receptor_path, prepared_receptor_path, pocket_method, center_json,
         box_json, manifest_entry_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(run_id, pdb_id, variant, ph_label) DO UPDATE SET
          library=COALESCE(excluded.library, library),
          status=COALESCE(excluded.status, status),
          failure_reason=COALESCE(excluded.failure_reason, failure_reason),
          input_receptor_path=COALESCE(excluded.input_receptor_path, input_receptor_path),
          prepared_receptor_path=COALESCE(excluded.prepared_receptor_path, prepared_receptor_path),
          pocket_method=COALESCE(excluded.pocket_method, pocket_method),
          center_json=COALESCE(excluded.center_json, center_json),
          box_json=COALESCE(excluded.box_json, box_json),
          manifest_entry_json=COALESCE(excluded.manifest_entry_json, manifest_entry_json)
        """,
        (
            run_id,
            pdb_id,
            variant,
            ph,
            values.get("library"),
            values.get("status"),
            values.get("failure_reason"),
            values.get("input_receptor_path"),
            values.get("prepared_receptor_path"),
            values.get("pocket_method"),
            _json(values["center"]) if values.get("center") is not None else None,
            _json(values["box"]) if values.get("box") is not None else None,
            _json(values["manifest_entry"])
            if values.get("manifest_entry") is not None
            else None,
        ),
    )
    row = connection.execute(
        "SELECT receptor_context_id FROM receptor_contexts WHERE run_id=? AND pdb_id=? AND variant=? AND ph_label=?",
        (run_id, pdb_id, variant, ph),
    ).fetchone()
    return int(row[0])


def _resolve_context(
    connection: sqlite3.Connection,
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
) -> tuple[str, str]:
    """Fill omitted legacy completion dimensions from an unambiguous context."""
    if variant and ph:
        return variant, ph
    rows = connection.execute(
        "SELECT variant, ph_label FROM receptor_contexts WHERE run_id=? AND pdb_id=?",
        (run_id, pdb_id),
    ).fetchall()
    candidates = [
        (str(row[0]), str(row[1]))
        for row in rows
        if (not variant or str(row[0]) == variant) and (not ph or str(row[1]) == ph)
    ]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise ValueError(
            "completion context is ambiguous for "
            f"run_id={run_id!r}, pdb_id={pdb_id!r}, variant={variant!r}, "
            f"ph_label={ph!r}: {candidates!r}"
        )
    return variant, ph


def _ligand(connection: sqlite3.Connection, canonical_id: str, **values: Any) -> int:
    connection.execute(
        """INSERT INTO ligands(canonical_id, display_name, source_path)
        VALUES (?, ?, ?) ON CONFLICT(canonical_id) DO UPDATE SET
        display_name=COALESCE(excluded.display_name, display_name),
        source_path=COALESCE(excluded.source_path, source_path)""",
        (canonical_id, values.get("display_name"), values.get("source_path")),
    )
    return int(
        connection.execute(
            "SELECT ligand_id FROM ligands WHERE canonical_id=?", (canonical_id,)
        ).fetchone()[0]
    )


def _pair(connection: sqlite3.Connection, context_id: int, ligand_id: int) -> int:
    connection.execute(
        "INSERT INTO pair_cells(receptor_context_id, ligand_id, final_status) VALUES (?, ?, 'expected') ON CONFLICT(receptor_context_id, ligand_id) DO NOTHING",
        (context_id, ligand_id),
    )
    return int(
        connection.execute(
            "SELECT pair_cell_id FROM pair_cells WHERE receptor_context_id=? AND ligand_id=?",
            (context_id, ligand_id),
        ).fetchone()[0]
    )


def _manifest_contexts(
    connection: sqlite3.Connection, run_id: str, manifest: Mapping[str, Any]
) -> None:
    proteins = _mapping(manifest.get("proteins"))
    for key, raw in proteins.items():
        if not isinstance(raw, Mapping):
            continue
        pdb_id, variant, ph = _context_values(str(key), raw)
        stages = _mapping(raw.get("stages"))
        prep = _mapping(stages.get("prep"))
        details = _mapping(prep.get("details"))
        pocket_method, center, box = _pocket(raw)
        _upsert_context(
            connection,
            run_id,
            pdb_id,
            variant,
            ph,
            library=_text(raw.get("library")) or None,
            status=_text(raw.get("status")) or None,
            failure_reason=_text(raw.get("error")) or None,
            input_receptor_path=_text(details.get("input_pdb")) or None,
            prepared_receptor_path=_text(
                details.get("receptor_pdbqt") or details.get("cleaned_pdb")
            )
            or None,
            pocket_method=pocket_method or None,
            center=center,
            box=box,
            manifest_entry=raw,
        )


def _completion_files(
    connection: sqlite3.Connection, run_id: str, docked: Path
) -> None:
    if not docked.is_dir():
        return
    for path in _completion_manifest_paths(docked):
        payload = _load_completion_payload(path)
        pdb_id = _text(payload.get("pdb_id")).upper()
        variant = _text(payload.get("variant")).upper()
        ph = _normalize_ph(payload.get("ph_label"))
        variant, ph = _resolve_context(connection, run_id, pdb_id, variant, ph)
        context_id = _upsert_context(connection, run_id, pdb_id, variant, ph)
        missing = {
            canonical_ligand_id(item)
            for item in payload.get("missing_ligands_after", [])
        }
        markers = _mapping(payload.get("failure_markers"))
        engine, stage = _text(payload.get("engine")), _text(payload.get("stage"))
        chunk = _text(payload.get("chunk_id"))
        connection.execute(
            """INSERT INTO completion_records
            (run_id, receptor_context_id, engine, stage, chunk_id,
             completion_path, completion_sha256, completion_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                run_id,
                context_id,
                engine,
                stage,
                chunk,
                str(path),
                _sha256(path),
                _json(payload),
            ),
        )
        completion_record_id = int(
            connection.execute(
                """SELECT completion_record_id FROM completion_records
            WHERE run_id=? AND completion_path=?""",
                (run_id, str(path)),
            ).fetchone()[0]
        )
        for ligand_path in payload.get("expected_ligands", []):
            canonical = canonical_ligand_id(ligand_path)
            if not canonical:
                continue
            ligand_id = _ligand(connection, canonical, source_path=_text(ligand_path))
            pair_id = _pair(connection, context_id, ligand_id)
            is_missing = canonical in missing
            reason = _text(
                markers.get(canonical) or markers.get(Path(_text(ligand_path)).name)
            )
            status = "failed" if is_missing or reason else "completed"
            connection.execute(
                "UPDATE pair_cells SET expected=1 WHERE pair_cell_id=?", (pair_id,)
            )
            connection.execute(
                """INSERT INTO docking_attempts
                (pair_cell_id, completion_record_id, engine, stage, chunk_id,
                 status, failure_code, failure_reason)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(pair_cell_id, completion_record_id, engine, stage,
                            chunk_id) DO UPDATE SET
                status=excluded.status, failure_code=excluded.failure_code,
                failure_reason=excluded.failure_reason""",
                (
                    pair_id,
                    completion_record_id,
                    engine,
                    stage,
                    chunk,
                    status,
                    "missing_engine_output"
                    if is_missing
                    else ("engine_failure" if reason else None),
                    reason
                    or ("expected ligand has no engine output" if is_missing else None),
                ),
            )


def _master_rows(connection: sqlite3.Connection, run_id: str, path: Path) -> None:
    if not path.is_file():
        return
    seen_master_keys: dict[tuple[int, str], int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row_number, row in enumerate(csv.DictReader(handle), start=2):
            pdb_id = _text(row.get("pdb_id")).upper()
            variant = _text(row.get("variant")).upper()
            ph = _normalize_ph(row.get("ph_label"))
            context_id = _upsert_context(
                connection,
                run_id,
                pdb_id,
                variant,
                ph,
                library=_text(row.get("library")) or None,
                pocket_method=_text(row.get("pocket_method")) or None,
                center=[
                    _float(row.get(name))
                    for name in ("center_x", "center_y", "center_z")
                ],
                box=[_float(row.get(name)) for name in ("box_x", "box_y", "box_z")],
            )
            canonical = canonical_ligand_id(
                row.get("ligand_base") or row.get("ligand_file") or row.get("ligand")
            )
            if not canonical:
                continue
            key = (context_id, canonical)
            previous_row = seen_master_keys.get(key)
            if previous_row is not None:
                raise ValueError(
                    f"{path}: duplicate master result key for run={run_id!r}, "
                    f"context_id={context_id}, ligand={canonical!r} in CSV rows "
                    f"{previous_row} and {row_number}; explicit attempt selection "
                    "is required"
                )
            seen_master_keys[key] = row_number
            ligand_id = _ligand(
                connection,
                canonical,
                display_name=_text(row.get("ligand_display")) or None,
            )
            pair_id = _pair(connection, context_id, ligand_id)
            valid_token = _text(row.get("pose_valid_any"))
            pose_valid = (
                1 if valid_token == "1" else (0 if valid_token == "0" else None)
            )
            score = _float(row.get("z_selected"))
            has_numeric = score is not None or any(
                _float(row.get(name)) is not None
                for name in ("selected_docking_score", "consensus_score", "final_score")
            )
            if pose_valid == 0:
                status = "invalid"
            elif pose_valid == 1:
                status = "valid"
            elif has_numeric:
                status = "calculated_unvalidated"
            else:
                status = "unsuccessful"
            reason = _text(
                row.get("pose_invalid_reason_top")
                or row.get("stage_fallback_reason")
                or row.get("dud_eval_status_reason")
            )
            connection.execute(
                """UPDATE pair_cells SET final_status=?, failure_code=?, failure_reason=?,
                has_result=1, pose_valid=?, is_control=?, is_decoy=?, atlas_score=?,
                atlas_score_source=?, selected_docking_score=?, consensus_score=?,
                final_score=?, final_score_source=?, final_rank=?, source_csv=?,
                result_json=? WHERE pair_cell_id=?""",
                (
                    status,
                    "pose_invalid"
                    if pose_valid == 0
                    else ("no_numeric_result" if not has_numeric else None),
                    reason or None,
                    pose_valid,
                    int(_text(row.get("is_control")) == "1"),
                    int(_text(row.get("is_decoy")) == "1"),
                    score,
                    _text(row.get("z_selected_source")) or None,
                    _float(row.get("selected_docking_score")),
                    _float(row.get("consensus_score")),
                    _float(row.get("final_score")),
                    _text(row.get("final_score_source")) or None,
                    _int(row.get("final_rank")),
                    _text(row.get("source_csv")) or None,
                    _json(row),
                    pair_id,
                ),
            )


def _finalize_missing(connection: sqlite3.Connection, run_id: str) -> None:
    connection.execute(
        """UPDATE pair_cells SET final_status='missing_result',
        failure_code=COALESCE(failure_code, 'missing_final_result'),
        failure_reason=COALESCE(failure_reason, 'scheduled pair has no master result row')
        WHERE has_result=0 AND expected=1 AND receptor_context_id IN
        (SELECT receptor_context_id FROM receptor_contexts WHERE run_id=?)""",
        (run_id,),
    )
    contexts = connection.execute(
        "SELECT receptor_context_id FROM receptor_contexts WHERE run_id=?", (run_id,)
    ).fetchall()
    for (context_id,) in contexts:
        controls = connection.execute(
            "SELECT final_status, failure_reason FROM pair_cells "
            "WHERE receptor_context_id=? AND is_control=1",
            (context_id,),
        ).fetchall()
        reason: str | None
        if not controls:
            status, reason = "missing", "no native-control result row"
        elif any(row[0] == "valid" for row in controls):
            status = "observed_valid_control_result_unqualified"
            reason = "control row observed; RMSD qualification not ingested"
        else:
            status = "observed_invalid_or_unsuccessful_control_result"
            reason = next((str(row[1]) for row in controls if row[1]), None)
        connection.execute(
            "UPDATE receptor_contexts SET native_redock_status=?, "
            "native_redock_reason=? WHERE receptor_context_id=?",
            (status, reason, context_id),
        )


def _artifacts(connection: sqlite3.Connection, run_id: str, path: Path) -> None:
    if not path.is_file():
        return
    try:
        index = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return
    for group in index.get("groups", []):
        if not isinstance(group, Mapping):
            continue
        group_key = _mapping(group.get("group"))
        pdb_id = _text(group_key.get("pdb_id")).upper()
        variant = _text(group_key.get("variant")).upper()
        ph = _normalize_ph(group_key.get("ph"))
        context_id = _upsert_context(connection, run_id, pdb_id, variant, ph)
        for entry in group.get("entries", []):
            if not isinstance(entry, Mapping):
                continue
            candidate = canonical_ligand_id(
                entry.get("member_name") or entry.get("original_path")
            )
            ligand_row = connection.execute(
                "SELECT ligand_id FROM ligands WHERE canonical_id=?", (candidate,)
            ).fetchone()
            pair_row = (
                connection.execute(
                    "SELECT pair_cell_id FROM pair_cells "
                    "WHERE receptor_context_id=? AND ligand_id=?",
                    (context_id, int(ligand_row[0])),
                ).fetchone()
                if ligand_row
                else None
            )
            pair_id = int(pair_row[0]) if pair_row else None
            connection.execute(
                """INSERT OR REPLACE INTO artifacts
                (run_id, receptor_context_id, pair_cell_id, stage, mode, original_path,
                 archive_path, member_name, sha256, size_bytes, file_type, verified, artifact_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    run_id,
                    context_id,
                    pair_id,
                    _text(entry.get("stage_dir")) or None,
                    _text(entry.get("mode")) or None,
                    _text(entry.get("original_path")) or None,
                    _text(entry.get("archive_path")) or None,
                    _text(entry.get("member_name")) or None,
                    _text(entry.get("sha256")) or None,
                    _int(entry.get("size_bytes")),
                    _text(entry.get("file_type")) or None,
                    int(bool(group.get("verified"))),
                    _json(entry),
                ),
            )


def _summary(
    connection: sqlite3.Connection, audit: Mapping[str, Any], database_path: Path
) -> dict[str, Any]:
    counts = {
        table: int(connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
        for table in (
            "runs",
            "receptor_contexts",
            "ligands",
            "pair_cells",
            "completion_records",
            "docking_attempts",
            "artifacts",
            "protein_identities",
            "receptor_annotations",
            "known_pair_selections",
        )
    }
    statuses = {
        row[0]: int(row[1])
        for row in connection.execute(
            "SELECT final_status, COUNT(*) FROM pair_cells GROUP BY final_status ORDER BY final_status"
        )
    }
    audited_completion_count = sum(
        int(run["completion_manifest_valid_count"]) for run in audit["runs"]
    )
    completion_import_complete = (
        counts["completion_records"] == audited_completion_count
    )
    return {
        "release_id": audit["release_id"],
        "database_path": str(database_path),
        "database_schema_version": SCHEMA_VERSION,
        "input_audit": audit,
        "counts": counts,
        "pair_status_counts": statuses,
        "completion_import_complete": completion_import_complete,
        "failure_complete": (
            audit["failure_complete_run_count"] == audit["run_count"]
            and completion_import_complete
            and statuses.get("expected", 0) == 0
        ),
    }


def build_release_database(
    manifest_path: Path,
    database_path: Path,
    repo_root: Path,
    *,
    overwrite: bool = False,
    summary_path: Path | None = None,
) -> dict[str, Any]:
    """Validate, import, and summarize one immutable Atlas release snapshot."""
    manifest = load_release_manifest(manifest_path)
    audit = audit_release_inputs(manifest_path, repo_root)
    if audit["errors"]:
        raise ValueError("release input audit failed: " + "; ".join(audit["errors"]))
    if database_path.exists():
        if not overwrite:
            raise FileExistsError(f"database already exists: {database_path}")
        database_path.unlink()
    database_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(database_path)
    try:
        create_schema(connection)
        now = datetime.now(timezone.utc).isoformat()
        connection.execute(
            "INSERT INTO releases VALUES (?, ?, ?, ?, ?, ?, ?)",
            (
                _text(manifest["release_id"]),
                _text(manifest["schema_version"]),
                _text(manifest.get("title")) or None,
                _text(manifest.get("created_at")) or None,
                _sha256(manifest_path),
                _json(manifest),
                now,
            ),
        )
        for name, value in manifest["scientific_policies"].items():
            connection.execute(
                "INSERT INTO scientific_policies VALUES (?, ?, ?)",
                (manifest["release_id"], name, _json(value)),
            )
        for entry in iter_run_entries(manifest):
            run_id = _text(entry["run_id"])
            paths = _run_paths(repo_root, entry)
            run_manifest = _load_mapping(paths["manifest"])
            connection.execute(
                """INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    run_id,
                    manifest["release_id"],
                    _text(run_manifest.get("status")) or None,
                    str(paths["manifest"]),
                    _sha256(paths["manifest"]),
                    _json(run_manifest.get("command")),
                    _json(run_manifest.get("git")),
                    _json(run_manifest.get("paths")),
                    _json(run_manifest.get("timing")),
                    _json(run_manifest.get("resources")),
                    _json(run_manifest),
                ),
            )
            _manifest_contexts(connection, run_id, run_manifest)
            _completion_files(connection, run_id, paths["docked"])
            _master_rows(connection, run_id, paths["master_rows"])
            _finalize_missing(connection, run_id)
            _artifacts(connection, run_id, paths["archive_index"])
        annotation_counts = ingest_scientific_annotations(
            connection, manifest, manifest_path
        )
        connection.commit()
        summary = _summary(connection, audit, database_path)
        summary["annotation_counts"] = annotation_counts
    finally:
        connection.close()
    destination = summary_path or database_path.with_suffix(".summary.json")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary
