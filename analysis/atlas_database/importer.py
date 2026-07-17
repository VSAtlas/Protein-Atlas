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
from analysis.atlas_database.attempt_selection import (
    AttemptKey,
    _key_label as _attempt_key_label,
    _selector_key as _attempt_selector_key,
    _selectors as _attempt_selectors,
    apply_attempt_selections,
)
from analysis.atlas_database.artifact_contract import (
    ArtifactContractError,
    artifact_role_counts,
    artifact_scope,
    load_artifact_index,
)
from analysis.atlas_database.manifest import iter_run_entries, load_release_manifest
from analysis.atlas_database.pose_validation_contract import (
    LEGACY_POSE_VALIDATION_METHOD,
    LEGACY_POSE_VALIDATION_SCOPE,
    legacy_pose_validation_thresholds_json,
)
from analysis.atlas_database.schema import SCHEMA_VERSION, create_schema
from analysis.atlas_database.score_source_contract import (
    build_score_source_materialization,
)
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


def _record_sha256(value: Mapping[str, Any]) -> str:
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _artifact_index_audit(path: Path) -> tuple[dict[str, int], list[str]]:
    if not path.is_file():
        return {}, []
    try:
        return artifact_role_counts(load_artifact_index(path)), []
    except ArtifactContractError as exc:
        return {}, [str(exc)]


def _group_attempt_candidate(
    first: dict[AttemptKey, Any],
    duplicates: dict[AttemptKey, list[Any]],
    key: AttemptKey,
    candidate: Any,
) -> None:
    if key in duplicates:
        duplicates[key].append(candidate)
    elif key in first:
        duplicates[key] = [first.pop(key), candidate]
    else:
        first[key] = candidate


def _completion_selector_present(selector: Mapping[str, Any]) -> bool:
    return bool(
        _text(selector.get("completion_relpath"))
        or _text(selector.get("completion_sha256"))
    )


def _result_selector_present(selector: Mapping[str, Any]) -> bool:
    return bool(
        selector.get("result_row_number") is not None
        or _text(selector.get("result_sha256"))
    )


def _pooled_attempt_key(
    value: Mapping[str, Any], token_pool: dict[str, str]
) -> AttemptKey:
    key = _attempt_selector_key(value)
    return (
        token_pool.setdefault(key[0], key[0]),
        token_pool.setdefault(key[1], key[1]),
        token_pool.setdefault(key[2], key[2]),
        token_pool.setdefault(key[3], key[3]),
    )


def _manifest_context_keys(manifest: Mapping[str, Any]) -> set[tuple[str, str, str]]:
    contexts: set[tuple[str, str, str]] = set()
    for key, raw in _mapping(manifest.get("proteins")).items():
        if isinstance(raw, Mapping):
            contexts.add(_context_values(str(key), raw))
    return contexts


def _completion_attempt_audit(
    *,
    run_id: str,
    docked: Path,
    completion_paths: list[Path],
    run_manifest: Mapping[str, Any],
    selectors: Mapping[AttemptKey, tuple[int, Mapping[str, Any]]],
) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    first: dict[AttemptKey, tuple[str, str]] = {}
    duplicates: dict[AttemptKey, list[tuple[str, str]]] = {}
    contexts = _manifest_context_keys(run_manifest)
    token_pool: dict[str, str] = {}
    valid_count = 0

    for path in completion_paths:
        try:
            payload = _load_completion_payload(path)
            pdb_id = _text(payload.get("pdb_id")).upper()
            variant = _text(payload.get("variant")).upper()
            ph = _normalize_ph(payload.get("ph_label"))
            variant, ph = _resolve_context_candidates(
                [
                    (candidate_variant, candidate_ph)
                    for candidate_pdb, candidate_variant, candidate_ph in contexts
                    if candidate_pdb == pdb_id
                ],
                run_id,
                pdb_id,
                variant,
                ph,
            )
        except ValueError as exc:
            errors.append({"path": str(path), "error": str(exc)})
            continue

        valid_count += 1
        contexts.add((pdb_id, variant, ph))
        missing = {
            canonical_ligand_id(item)
            for item in payload.get("missing_ligands_after", [])
        }
        markers = _mapping(payload.get("failure_markers"))
        candidate = (path.relative_to(docked).as_posix(), _sha256(path))
        for ligand_path in payload.get("expected_ligands", []):
            canonical = canonical_ligand_id(ligand_path)
            reason = _text(
                markers.get(canonical) or markers.get(Path(_text(ligand_path)).name)
            )
            if canonical in missing or reason:
                continue
            key = _pooled_attempt_key(
                {
                    "pdb_id": pdb_id,
                    "variant": variant,
                    "ph_label": ph,
                    "ligand_canonical_id": canonical,
                },
                token_pool,
            )
            _group_attempt_candidate(first, duplicates, key, candidate)

    selection_errors: list[str] = []
    for key in sorted(duplicates):
        candidates = duplicates[key]
        selected = selectors.get(key)
        has_explicit = bool(selected and _completion_selector_present(selected[1]))
        if not has_explicit:
            choices = ", ".join(
                f"{path} sha256={digest}" for path, digest in candidates
            )
            selection_errors.append(
                "multiple successful completion attempts for "
                f"{_attempt_key_label(key)}; add an explicit attempt_selections "
                f"entry; candidates: {choices}"
            )

    for key, (_, selector) in selectors.items():
        if not _completion_selector_present(selector):
            continue
        if key not in first and key not in duplicates:
            selection_errors.append(
                "completion selector does not match a successful attempt for "
                f"{_attempt_key_label(key)}"
            )
            continue
        candidates = duplicates.get(key) or [first[key]]
        expected = (
            _text(selector.get("completion_relpath")),
            _text(selector.get("completion_sha256")).lower(),
        )
        match_count = sum(candidate == expected for candidate in candidates)
        if match_count != 1:
            selection_errors.append(
                f"completion selector matched {match_count} successful attempts "
                f"for {_attempt_key_label(key)}"
            )

    return {
        "errors": errors,
        "valid_count": valid_count,
        "successful_attempt_count": len(first)
        + sum(len(values) for values in duplicates.values()),
        "successful_pair_count": len(first) + len(duplicates),
        "duplicate_successful_pair_count": len(duplicates),
        "selection_errors": selection_errors,
    }


def _result_attempt_audit(
    path: Path,
    selectors: Mapping[AttemptKey, tuple[int, Mapping[str, Any]]],
) -> dict[str, Any]:
    if not path.is_file():
        return {
            "errors": [],
            "row_count": 0,
            "pair_count": 0,
            "duplicate_pair_count": 0,
            "selection_errors": [
                "result selector does not match an attempt for "
                f"{_attempt_key_label(key)}"
                for key, (_, selector) in selectors.items()
                if _result_selector_present(selector)
            ],
        }

    first_rows: dict[AttemptKey, int] = {}
    duplicate_keys: set[AttemptKey] = set()
    selector_match_counts: dict[AttemptKey, int] = {}
    token_pool: dict[str, str] = {}
    errors: list[dict[str, str]] = []
    row_count = 0
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                pdb_id = _text(row.get("pdb_id")).upper()
                variant = _text(row.get("variant")).upper()
                ph = _normalize_ph(row.get("ph_label"))
                canonical = canonical_ligand_id(
                    row.get("ligand_base")
                    or row.get("ligand_file")
                    or row.get("ligand")
                )
                if not canonical:
                    continue
                row_count += 1
                key = _pooled_attempt_key(
                    {
                        "pdb_id": pdb_id,
                        "variant": variant,
                        "ph_label": ph,
                        "ligand_canonical_id": canonical,
                    },
                    token_pool,
                )
                if key in first_rows:
                    duplicate_keys.add(key)
                else:
                    first_rows[key] = row_number
                selected = selectors.get(key)
                if selected and _result_selector_present(selected[1]):
                    expected = (
                        selected[1].get("result_row_number"),
                        _text(selected[1].get("result_sha256")).lower(),
                    )
                    candidate = (row_number, _record_sha256(row))
                    if candidate == expected:
                        selector_match_counts[key] = (
                            selector_match_counts.get(key, 0) + 1
                        )
    except (OSError, UnicodeDecodeError, csv.Error, TypeError, ValueError) as exc:
        errors.append({"path": str(path), "error": str(exc)})

    choice_keys = {
        key
        for key in duplicate_keys
        if not (selectors.get(key) and _result_selector_present(selectors[key][1]))
    }
    choices: dict[AttemptKey, list[tuple[int, str]]] = {key: [] for key in choice_keys}
    if choice_keys and not errors:
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                for row_number, row in enumerate(csv.DictReader(handle), start=2):
                    canonical = canonical_ligand_id(
                        row.get("ligand_base")
                        or row.get("ligand_file")
                        or row.get("ligand")
                    )
                    if not canonical:
                        continue
                    key = _pooled_attempt_key(
                        {
                            "pdb_id": _text(row.get("pdb_id")).upper(),
                            "variant": _text(row.get("variant")).upper(),
                            "ph_label": _normalize_ph(row.get("ph_label")),
                            "ligand_canonical_id": canonical,
                        },
                        token_pool,
                    )
                    if key in choices:
                        choices[key].append((row_number, _record_sha256(row)))
        except (OSError, UnicodeDecodeError, csv.Error, TypeError, ValueError) as exc:
            errors.append({"path": str(path), "error": str(exc)})

    selection_errors: list[str] = []
    for key in sorted(choice_keys):
        candidate_text = ", ".join(
            f"record={row_number} sha256={digest}"
            for row_number, digest in choices[key]
        )
        selection_errors.append(
            f"multiple result attempts for {_attempt_key_label(key)}; add an "
            f"explicit attempt_selections entry; candidates: {candidate_text}"
        )

    for key, (_, selector) in selectors.items():
        if not _result_selector_present(selector):
            continue
        if key not in first_rows:
            selection_errors.append(
                f"result selector does not match an attempt for {_attempt_key_label(key)}"
            )
            continue
        match_count = selector_match_counts.get(key, 0)
        if match_count != 1:
            selection_errors.append(
                f"result selector matched {match_count} attempts for "
                f"{_attempt_key_label(key)}"
            )

    return {
        "errors": errors,
        "row_count": row_count,
        "pair_count": len(first_rows),
        "duplicate_pair_count": len(duplicate_keys),
        "selection_errors": selection_errors,
    }


def audit_release_inputs(manifest_path: Path, repo_root: Path) -> dict[str, Any]:
    """Audit release inputs without creating a database."""
    manifest = load_release_manifest(manifest_path)
    runs: list[dict[str, Any]] = []
    for entry in iter_run_entries(manifest):
        paths = _run_paths(repo_root, entry)
        manifest_exists = paths["manifest"].is_file()
        master_exists = paths["master_rows"].is_file()
        completion_paths = _completion_manifest_paths(paths["docked"])
        run_manifest: dict[str, Any] = {}
        run_manifest_errors: list[dict[str, str]] = []
        if manifest_exists:
            try:
                run_manifest = _load_mapping(paths["manifest"])
            except (OSError, UnicodeDecodeError, yaml.YAMLError) as exc:
                run_manifest_errors.append(
                    {"path": str(paths["manifest"]), "error": str(exc)}
                )
        selectors = _attempt_selectors(entry.get("attempt_selections", []))
        completion_audit = _completion_attempt_audit(
            run_id=_text(entry["run_id"]),
            docked=paths["docked"],
            completion_paths=completion_paths,
            run_manifest=run_manifest,
            selectors=selectors,
        )
        result_audit = _result_attempt_audit(paths["master_rows"], selectors)
        attempt_selection_errors = [
            *completion_audit["selection_errors"],
            *result_audit["selection_errors"],
        ]
        artifact_roles, artifact_errors = _artifact_index_audit(paths["archive_index"])
        completion_errors = completion_audit["errors"]
        completion_valid_count = completion_audit["valid_count"]
        runs.append(
            {
                "run_id": entry["run_id"],
                "paths": {name: str(path) for name, path in paths.items()},
                "run_manifest_present": manifest_exists,
                "run_manifest_errors": run_manifest_errors,
                "master_rows_present": master_exists,
                "master_rows_valid": master_exists and not result_audit["errors"],
                "master_rows_errors": result_audit["errors"],
                "master_result_row_count": result_audit["row_count"],
                "master_result_pair_count": result_audit["pair_count"],
                "duplicate_result_pair_count": result_audit["duplicate_pair_count"],
                "archive_index_present": paths["archive_index"].is_file(),
                "archive_index_valid": not artifact_errors,
                "archive_index_errors": artifact_errors,
                "artifact_role_counts": artifact_roles,
                "completion_manifest_count": len(completion_paths),
                "completion_manifest_valid_count": completion_valid_count,
                "completion_manifest_invalid_count": len(completion_errors),
                "completion_manifest_errors": completion_errors,
                "successful_completion_attempt_count": completion_audit[
                    "successful_attempt_count"
                ],
                "successful_completion_pair_count": completion_audit[
                    "successful_pair_count"
                ],
                "duplicate_successful_completion_pair_count": completion_audit[
                    "duplicate_successful_pair_count"
                ],
                "attempt_selection_valid": not attempt_selection_errors,
                "attempt_selection_errors": attempt_selection_errors,
                "failure_complete_inputs_present": manifest_exists
                and master_exists
                and not run_manifest_errors
                and not result_audit["errors"]
                and completion_valid_count > 0
                and not completion_errors
                and not attempt_selection_errors,
            }
        )
    errors = [
        f"{run['run_id']}: missing run manifest"
        for run in runs
        if not run["run_manifest_present"]
    ]
    for run in runs:
        for issue in run["run_manifest_errors"]:
            errors.append(
                f"{run['run_id']}: invalid run manifest "
                f"{issue['path']}: {issue['error']}"
            )
        for issue in run["completion_manifest_errors"]:
            errors.append(
                f"{run['run_id']}: invalid completion manifest "
                f"{issue['path']}: {issue['error']}"
            )
        for issue in run["master_rows_errors"]:
            errors.append(
                f"{run['run_id']}: invalid master result table "
                f"{issue['path']}: {issue['error']}"
            )
        for issue in run["attempt_selection_errors"]:
            errors.append(f"{run['run_id']}: {issue}")
        for issue in run["archive_index_errors"]:
            errors.append(
                f"{run['run_id']}: invalid artifact index "
                f"{run['paths']['archive_index']}: {issue}"
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


def _resolve_context_candidates(
    context_values: list[tuple[str, str]],
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
) -> tuple[str, str]:
    if variant and ph:
        return variant, ph
    candidates = [
        (candidate_variant, candidate_ph)
        for candidate_variant, candidate_ph in context_values
        if (not variant or candidate_variant == variant)
        and (not ph or candidate_ph == ph)
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


def _resolve_context(
    connection: sqlite3.Connection,
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
) -> tuple[str, str]:
    """Fill omitted legacy completion dimensions from an unambiguous context."""
    rows = connection.execute(
        "SELECT variant, ph_label FROM receptor_contexts WHERE run_id=? AND pdb_id=?",
        (run_id, pdb_id),
    ).fetchall()
    return _resolve_context_candidates(
        [(str(row[0]), str(row[1])) for row in rows],
        run_id,
        pdb_id,
        variant,
        ph,
    )


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
    input_csv_sha256 = _sha256(path)
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
            final_score = _float(row.get("final_score"))
            result_sha256 = _record_sha256(row)
            (
                final_score_source_reconstructed,
                final_score_source_classification,
                final_score_source_evidence_json,
            ) = build_score_source_materialization(
                row,
                input_csv_sha256=input_csv_sha256,
                source_row_number=row_number,
                result_sha256=result_sha256,
            )
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
                """INSERT INTO result_attempts
                (pair_cell_id, input_csv_path, input_csv_sha256, source_row_number,
                 result_sha256, final_status, failure_code, failure_reason,
                 pose_valid, pose_validation_method, pose_validation_scope,
                 pose_validation_thresholds_json, is_control, is_decoy,
                 atlas_score, atlas_score_source,
                 selected_docking_score, consensus_score, final_score,
                 final_score_source, final_score_source_reconstructed,
                 final_score_source_classification,
                 final_score_source_evidence_json,
                 final_rank, source_csv, result_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    pair_id,
                    str(path),
                    input_csv_sha256,
                    row_number,
                    result_sha256,
                    status,
                    "pose_invalid"
                    if pose_valid == 0
                    else ("no_numeric_result" if not has_numeric else None),
                    reason or None,
                    pose_valid,
                    (LEGACY_POSE_VALIDATION_METHOD if pose_valid is not None else None),
                    LEGACY_POSE_VALIDATION_SCOPE if pose_valid is not None else None,
                    (
                        legacy_pose_validation_thresholds_json()
                        if pose_valid is not None
                        else None
                    ),
                    int(_text(row.get("is_control")) == "1"),
                    int(_text(row.get("is_decoy")) == "1"),
                    score,
                    _text(row.get("z_selected_source")) or None,
                    _float(row.get("selected_docking_score")),
                    _float(row.get("consensus_score")),
                    final_score,
                    _text(row.get("final_score_source")) or None,
                    final_score_source_reconstructed,
                    final_score_source_classification,
                    final_score_source_evidence_json,
                    _int(row.get("final_rank")),
                    _text(row.get("source_csv")) or None,
                    _json(row),
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
    index = load_artifact_index(path)
    for group in index.get("groups", []):
        group_key = _mapping(group.get("group"))
        pdb_id = _text(group_key.get("pdb_id")).upper()
        variant = _text(group_key.get("variant")).upper()
        ph = _normalize_ph(group_key.get("ph"))
        variant, ph = _resolve_context(connection, run_id, pdb_id, variant, ph)
        context_row = connection.execute(
            """SELECT receptor_context_id FROM receptor_contexts
               WHERE run_id=? AND pdb_id=? AND variant=? AND ph_label=?""",
            (run_id, pdb_id, variant, ph),
        ).fetchone()
        if context_row is None:
            raise ValueError(
                "artifact group does not match an imported receptor context: "
                f"run_id={run_id!r}, pdb_id={pdb_id!r}, variant={variant!r}, "
                f"ph_label={ph!r}"
            )
        context_id = int(context_row[0])
        for entry in group.get("entries", []):
            role = _text(entry.get("artifact_role"))
            scope = artifact_scope(role)
            associated_context_id = (
                context_id if scope in {"receptor", "pair"} else None
            )
            ligand_id: int | None = None
            if scope in {"ligand", "pair"}:
                declared_ligand = _text(entry.get("ligand_canonical_id"))
                canonical = canonical_ligand_id(declared_ligand)
                if not canonical or canonical != declared_ligand:
                    raise ValueError(
                        "artifact ligand_canonical_id must already be canonical: "
                        f"{declared_ligand!r}"
                    )
                ligand_id = _ligand(connection, canonical)
            pair_row = None
            if scope == "pair" and ligand_id is not None:
                pair_row = connection.execute(
                    "SELECT pair_cell_id FROM pair_cells "
                    "WHERE receptor_context_id=? AND ligand_id=?",
                    (context_id, ligand_id),
                ).fetchone()
            pair_id = int(pair_row[0]) if pair_row else None
            connection.execute(
                """INSERT INTO artifacts
                (run_id, receptor_context_id, ligand_id, pair_cell_id,
                 artifact_role, artifact_scope, stage, mode, original_path,
                 archive_path, member_name, sha256, size_bytes, file_type,
                 verified, artifact_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    run_id,
                    associated_context_id,
                    ligand_id,
                    pair_id,
                    role,
                    scope,
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
            "ligand_stereo_evidence",
            "receptor_audits",
            "pair_cells",
            "completion_records",
            "docking_attempts",
            "result_attempts",
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
            apply_attempt_selections(
                connection, run_id, paths["docked"], entry.get("attempt_selections", [])
            )
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
