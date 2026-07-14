"""CLI for auditing or building an Atlas release database."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sqlite3
import time
from pathlib import Path
from typing import Any

from analysis.atlas_database import (
    audit_release_inputs,
    build_release_database,
    load_release_manifest,
)
from config.output_paths import output_root


_PUBLIC_PATH_COLUMNS = {
    "runs": ("manifest_path",),
    "receptor_contexts": ("input_receptor_path", "prepared_receptor_path"),
    "ligands": ("source_path",),
    "pair_cells": ("source_csv",),
    "completion_records": ("completion_path",),
    "result_attempts": ("input_csv_path", "source_csv"),
    "artifacts": ("original_path", "archive_path"),
    "protein_identities": ("source_path",),
    "receptor_annotations": ("source_path",),
    "receptor_audits": ("source_path",),
    "known_pair_selections": ("source_path",),
}
_PUBLIC_JSON_COLUMNS = {
    "releases": ("manifest_json",),
    "scientific_policies": ("policy_json",),
    "runs": (
        "command_json",
        "git_json",
        "paths_json",
        "timing_json",
        "resources_json",
        "manifest_json",
    ),
    "receptor_contexts": ("manifest_entry_json",),
    "ligands": ("prepared_state_json",),
    "pair_cells": ("result_json", "pose_validation_thresholds_json"),
    "completion_records": ("completion_json",),
    "result_attempts": (
        "result_json", "pose_validation_thresholds_json",
        "completion_link_evidence_json",
    ),
    "artifacts": ("artifact_json",),
    "protein_identities": ("source_record_json", "provenance_json"),
    "receptor_annotations": ("source_record_json", "provenance_json"),
    "receptor_audits": ("audit_json",),
    "known_pair_selections": ("source_record_json", "provenance_json"),
}
_ABSOLUTE_PUBLIC_PATH = re.compile(
    r"(?<![A-Za-z0-9._:/-])/(?!/)[^\s\"'\\,\]\}]+"
)
_FORBIDDEN_PUBLIC_PREFIXES = (
    "/stor/",
    "/home/",
    "/tmp/",
    "/scratch/",
    "/work/",
    "/mnt/",
    "/var/",
    "/opt/",
    "/root/",
    "/Users/",
)
_OMITTED_COMPLETION_JSON = json.dumps(
    {
        "public_projection": "omitted",
        "reconstruct_with": "completion_sha256 and the private evidence ledger",
    },
    sort_keys=True,
    separators=(",", ":"),
)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _default_out_dir(repo_root: Path, manifest_path: Path) -> Path:
    manifest = load_release_manifest(manifest_path)
    return output_root(repo_root, "data") / str(manifest["release_id"])


def _safe_public_path(value: str, repo_root: Path) -> str:
    path = Path(value)
    if not path.is_absolute():
        return value
    try:
        return path.resolve().relative_to(repo_root).as_posix()
    except ValueError:
        path_hash = hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]
        return f"external/{path_hash}/{path.name or 'root'}"


def _sanitize_public_string(value: str, repo_root: Path) -> str:
    if Path(value).is_absolute():
        return _safe_public_path(value, repo_root)
    return _ABSOLUTE_PUBLIC_PATH.sub(
        lambda match: _safe_public_path(match.group(0), repo_root), value
    )


def _sanitize_json_value(value: Any, repo_root: Path, counts: dict[str, int]) -> Any:
    if isinstance(value, dict):
        return {
            str(key): _sanitize_json_value(item, repo_root, counts)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_sanitize_json_value(item, repo_root, counts) for item in value]
    if isinstance(value, str):
        sanitized = _sanitize_public_string(value, repo_root)
        if sanitized != value:
            counts["json_path_values_redacted"] += 1
        return sanitized
    return value


def _sanitize_public_database(database_path: Path, repo_root: Path) -> dict[str, int]:
    """Create a compact public projection while retaining the private ledger."""
    counts = {
        "path_columns_redacted": 0,
        "json_path_values_redacted": 0,
        "score_sources_materialized": 0,
        "pair_result_json_omitted": 0,
        "result_attempt_json_omitted": 0,
        "completion_json_omitted": 0,
    }
    with sqlite3.connect(database_path) as connection:
        rows = connection.execute(
            """SELECT pair_cell_id, result_json FROM pair_cells
            WHERE final_score IS NOT NULL
              AND COALESCE(TRIM(final_score_source), '') = ''
              AND result_json IS NOT NULL"""
        ).fetchall()
        for pair_cell_id, value in rows:
            try:
                payload = json.loads(str(value))
            except (TypeError, ValueError):
                continue
            source = str(payload.get("final_score_source") or "").strip()
            if source:
                connection.execute(
                    "UPDATE pair_cells SET final_score_source=? WHERE pair_cell_id=?",
                    (source, pair_cell_id),
                )
                counts["score_sources_materialized"] += 1
        counts["pair_result_json_omitted"] = int(
            connection.execute(
                "SELECT COUNT(*) FROM pair_cells WHERE result_json IS NOT NULL"
            ).fetchone()[0]
        )
        connection.execute(
            "UPDATE pair_cells SET result_json=NULL WHERE result_json IS NOT NULL"
        )
        counts["result_attempt_json_omitted"] = int(
            connection.execute(
                "SELECT COUNT(*) FROM result_attempts WHERE result_json IS NOT NULL"
            ).fetchone()[0]
        )
        connection.execute(
            "UPDATE result_attempts SET result_json=NULL WHERE result_json IS NOT NULL"
        )
        counts["completion_json_omitted"] = int(
            connection.execute(
                "SELECT COUNT(*) FROM completion_records WHERE completion_json != ?",
                (_OMITTED_COMPLETION_JSON,),
            ).fetchone()[0]
        )
        connection.execute(
            "UPDATE completion_records SET completion_json=?",
            (_OMITTED_COMPLETION_JSON,),
        )
        for table, columns in _PUBLIC_PATH_COLUMNS.items():
            for column in columns:
                rows = connection.execute(
                    f"SELECT rowid, {column} FROM {table} WHERE {column} IS NOT NULL"
                ).fetchall()
                for rowid, value in rows:
                    sanitized = _sanitize_public_string(str(value), repo_root)
                    if sanitized != value:
                        counts["path_columns_redacted"] += 1
                        connection.execute(
                            f"UPDATE {table} SET {column}=? WHERE rowid=?",
                            (sanitized, rowid),
                        )
        for table, columns in _PUBLIC_JSON_COLUMNS.items():
            for column in columns:
                rows = connection.execute(
                    f"SELECT rowid, {column} FROM {table} WHERE {column} IS NOT NULL"
                ).fetchall()
                for rowid, value in rows:
                    try:
                        decoded = json.loads(str(value))
                    except (TypeError, ValueError):
                        continue
                    sanitized = json.dumps(
                        _sanitize_json_value(decoded, repo_root, counts),
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    connection.execute(
                        f"UPDATE {table} SET {column}=? WHERE rowid=?",
                        (sanitized, rowid),
                    )
        connection.commit()
        connection.execute("VACUUM")
    return counts


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_entry(path: Path) -> dict[str, str]:
    filename_labels = {
        "release_readiness.json": "Release readiness audit",
        "image_plan.json": "Conference image plan",
        "database_redaction.json": "Database public-projection metadata",
    }
    labels = {
        ".sqlite": "SQLite release database",
        ".csv": "Complete pair table (CSV)",
        ".parquet": "Complete pair table (Parquet)",
        ".json": "Database public-projection metadata",
    }
    return {
        "label": filename_labels.get(
            path.name, labels.get(path.suffix.lower(), path.name)
        ),
        "url": f"downloads/{path.name}",
        "content_hash": f"sha256:{_sha256(path)}",
    }


def _assert_public_database_safe(database_path: Path) -> None:
    with sqlite3.connect(database_path) as connection:
        tables = [
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
            )
        ]
        for table in tables:
            text_columns = [
                row[1]
                for row in connection.execute(f"PRAGMA table_info({table})")
                if "TEXT" in str(row[2]).upper()
            ]
            for column in text_columns:
                for (value,) in connection.execute(
                    f"SELECT {column} FROM {table} WHERE {column} IS NOT NULL"
                ):
                    text = str(value)
                    if (
                        Path(text).is_absolute()
                        or _ABSOLUTE_PUBLIC_PATH.search(text)
                        or any(
                            prefix in text for prefix in _FORBIDDEN_PUBLIC_PREFIXES
                        )
                    ):
                        raise ValueError(
                            f"public database contains a machine-local path in {table}.{column}"
                        )


def _assert_public_text_safe(site_dir: Path) -> None:
    text_suffixes = {".csv", ".css", ".html", ".js", ".json"}
    for path in site_dir.rglob("*"):
        if path.is_file() and path.suffix.lower() in text_suffixes:
            text = path.read_text(encoding="utf-8")
            if any(prefix in text for prefix in _FORBIDDEN_PUBLIC_PREFIXES):
                raise ValueError(
                    f"public file contains a machine-local path: {path.name}"
                )


def _remove_publication_site(site_dir: Path) -> None:
    """Remove only the generated site, tolerating brief shared-filesystem lag."""
    for attempt in range(3):
        try:
            shutil.rmtree(site_dir)
            return
        except OSError:
            if attempt == 2:
                raise
            time.sleep(0.05 * (attempt + 1))


def _build_publication_site(
    manifest_path: Path,
    repo_root: Path,
    out_dir: Path,
    *,
    overwrite: bool,
    include_parquet: bool,
) -> dict[str, Any]:
    from analysis.atlas_database.exports import export_release_database
    from analysis.atlas_database.image_plan import write_release_image_plan
    from analysis.atlas_database.readiness import audit_release_readiness
    from analysis.reporting.docking_atlas_static import generate_static_explorer

    site_dir = out_dir / "site"
    if site_dir.is_symlink():
        raise ValueError(f"refusing to replace symlinked publication site: {site_dir}")
    if site_dir.exists():
        if not overwrite:
            raise FileExistsError(
                f"publication site already exists: {site_dir}; use --overwrite"
            )
        _remove_publication_site(site_dir)

    downloads_dir = site_dir / "downloads"
    downloads_dir.mkdir(parents=True)
    private_database_path = out_dir / "docking_atlas.private.sqlite"
    build_release_database(
        manifest_path,
        private_database_path,
        repo_root,
        overwrite=overwrite,
        summary_path=out_dir / "build_summary.private.json",
    )
    database_path = downloads_dir / "docking_atlas.sqlite"
    shutil.copy2(private_database_path, database_path)
    source_database_sha256 = _sha256(private_database_path)
    redaction_counts = _sanitize_public_database(database_path, repo_root)
    public_database_sha256 = _sha256(database_path)
    redaction_path = downloads_dir / "database_redaction.json"
    _write_json(
        redaction_path,
        {
            "projection": "public_compact_path_redacted_copy",
            "source_database_sha256": f"sha256:{source_database_sha256}",
            "public_database_sha256": f"sha256:{public_database_sha256}",
            "redaction_counts": redaction_counts,
            "manifest_sha256_semantics": (
                "The releases.manifest_sha256 field hashes the original release "
                "manifest. The private database retains raw result/completion JSON; "
                "the public database is a compact, path-redacted projection."
            ),
        },
    )
    _assert_public_database_safe(database_path)
    export_summary = export_release_database(
        database_path,
        downloads_dir,
        include_parquet=include_parquet,
    )
    readiness_path = downloads_dir / "release_readiness.json"
    readiness = audit_release_readiness(database_path, output_path=readiness_path)
    image_plan_path = downloads_dir / "image_plan.json"
    image_plan = write_release_image_plan(
        database_path,
        image_plan_path,
        image_output_root="release_images",
    )

    download_paths = [
        database_path,
        downloads_dir / "pairs.csv",
        readiness_path,
        image_plan_path,
        redaction_path,
    ]
    parquet_path = downloads_dir / "pairs.parquet"
    if parquet_path.is_file():
        download_paths.append(parquet_path)
    downloads = [_download_entry(path) for path in download_paths]

    browser_path = downloads_dir / "release_browser.json"
    payload = json.loads(browser_path.read_text(encoding="utf-8"))
    payload["downloads"] = downloads
    browser_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    site_manifest = generate_static_explorer(payload, site_dir)
    site_manifest["downloads"] = downloads
    _write_json(site_dir / "site_manifest.json", site_manifest)
    _assert_public_text_safe(site_dir)
    _assert_public_database_safe(database_path)
    from analysis.atlas_database.release_bundle import prepare_release_bundle

    integrity = prepare_release_bundle(site_dir)
    if integrity["status"] != "passed":
        raise ValueError(
            "generated release bundle failed integrity verification: "
            f"{integrity['error_count']} errors"
        )

    return {
        "release_id": site_manifest["release_id"],
        "site_dir": "site",
        "entrypoint": "site/index.html",
        "downloads": downloads,
        "pair_count": export_summary["pair_count"],
        "primary_score_count": export_summary["primary_score_count"],
        "rank_eligible_count": export_summary["rank_eligible_count"],
        "readiness_blocker_count": readiness["blocker_count"],
        "readiness_explicitly_qualified_receptor_count": readiness["metrics"][
            "native_redocking"
        ]["explicitly_qualified_count"],
        "image_plan_context_count": image_plan["summary"]["context_count"],
        "image_plan_selection_count": image_plan["summary"]["selection_count"],
        "image_plan_gap_count": image_plan["summary"]["gap_count"],
        "pairs_parquet": export_summary["pairs_parquet"],
        "parquet_error": export_summary["parquet_error"],
        "bundle_integrity_status": integrity["status"],
        "bundle_inventory": "site/release_inventory.json",
        "bundle_checksums": "site/release_checksums.sha256",
        "bundle_verification_report": "site/release_verification.json",
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit or build a failure-complete Atlas release database."
    )
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("audit", "build"):
        child = sub.add_parser(command)
        child.add_argument("--manifest", required=True, type=Path)
        child.add_argument("--repo-root", type=Path, default=Path("."))
        child.add_argument("--out-dir", type=Path)
    sub.choices["audit"].add_argument("--strict", action="store_true")
    sub.choices["build"].add_argument("--overwrite", action="store_true")
    sub.choices["build"].add_argument(
        "--no-parquet", action="store_false", dest="include_parquet", default=True
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    repo_root = args.repo_root.resolve()
    manifest_path = args.manifest.resolve()
    out_dir = (
        args.out_dir.resolve()
        if args.out_dir is not None
        else _default_out_dir(repo_root, manifest_path)
    )

    if args.command == "audit":
        summary = audit_release_inputs(manifest_path, repo_root)
        _write_json(out_dir / "audit.json", summary)
        strict_failure = bool(args.strict) and (
            bool(summary.get("errors"))
            or summary.get("failure_complete_run_count") != summary.get("run_count")
        )
        return_code = 2 if strict_failure or summary.get("errors") else 0
    else:
        summary = _build_publication_site(
            manifest_path,
            repo_root,
            out_dir,
            overwrite=bool(args.overwrite),
            include_parquet=bool(args.include_parquet),
        )
        _write_json(out_dir / "build_summary.json", summary)
        return_code = 0

    print(json.dumps(summary, indent=2, sort_keys=True))
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
