from __future__ import annotations

import json
import sqlite3
from copy import deepcopy
from pathlib import Path

import pytest

from analysis.atlas_database.artifact_contract import (
    ARTIFACT_PUBLICATION_POLICIES,
    ARTIFACT_ROLE_SCOPES,
    ArtifactContractError,
    artifact_public_metadata_allowed,
    artifact_public_record_allowed,
    artifact_publication_policy,
    artifact_scope,
    validate_artifact_index,
)
from analysis.atlas_database.manifest import validate_release_manifest
from analysis.atlas_database.schema import create_schema
from analysis.cli.build_atlas_database import _sanitize_public_database


def _artifact_index() -> dict[str, object]:
    return {
        "groups": [
            {
                "group": {"pdb_id": "1ABC", "variant": "HOLO", "ph": "base"},
                "verified": True,
                "entries": [
                    {
                        "artifact_role": "docking_pose",
                        "ligand_canonical_id": "drug-a",
                        "archive_path": "artifacts.tar.zst",
                        "member_name": "docked/run/1ABC/drug-a.pdbqt",
                        "sha256": "a" * 64,
                        "size_bytes": 12,
                        "file_type": "pdbqt",
                    }
                ],
            }
        ]
    }


def _release_manifest() -> dict[str, object]:
    return {
        "release_id": "atlas-contract-test",
        "schema_version": 1,
        "runs": [{"run_id": "run-a"}],
        "scientific_policies": {
            "receptor_selection": {"policy": "declared"},
            "native_redocking": {"policy": "declared"},
            "failure_handling": {"policy": "declared"},
            "normalization": {"policy": "declared"},
        },
    }


def test_typed_artifact_index_accepts_explicit_pair_association() -> None:
    assert validate_artifact_index(_artifact_index()) == []


def test_typed_artifact_index_rejects_unknown_role_and_missing_ligand() -> None:
    index = _artifact_index()
    entry = index["groups"][0]["entries"][0]  # type: ignore[index]
    entry["artifact_role"] = "generic_file"
    entry.pop("ligand_canonical_id")

    errors = validate_artifact_index(index)

    assert any("artifact_role must be one of" in error for error in errors)

    entry["artifact_role"] = "docking_pose"
    errors = validate_artifact_index(index)
    assert any("ligand_canonical_id is required" in error for error in errors)


def test_raw_logs_are_pair_scoped_and_require_explicit_ligand() -> None:
    index = _artifact_index()
    entry = index["groups"][0]["entries"][0]  # type: ignore[index]
    entry["artifact_role"] = "raw_stdout"
    entry.pop("ligand_canonical_id")

    errors = validate_artifact_index(index)

    assert any("ligand_canonical_id is required" in error for error in errors)


def test_approved_artifact_publication_policy_is_explicit_and_fail_closed() -> None:
    assert set(ARTIFACT_PUBLICATION_POLICIES) == set(ARTIFACT_ROLE_SCOPES)
    assert artifact_scope("native_ligand") == "receptor"
    assert artifact_scope("native_redock_pose") == "receptor"
    assert artifact_scope("release_manifest") == "run"
    assert artifact_publication_policy("prepared_receptor") == "public_approved"
    assert (
        artifact_publication_policy("docking_pose") == "public_if_selected_for_release"
    )
    assert artifact_publication_policy("raw_stdout") == "private"
    assert not artifact_public_metadata_allowed("raw_stdout", {})
    assert not artifact_public_metadata_allowed("docking_pose", {})
    assert artifact_public_metadata_allowed(
        "docking_pose", {"selected_for_release": True}
    )
    assert artifact_public_record_allowed(
        "docking_pose",
        {"selected_for_release": True},
        verified=1,
        sha256="a" * 64,
    )
    assert not artifact_public_record_allowed(
        "prepared_receptor", {}, verified=0, sha256="a" * 64
    )
    assert not artifact_public_record_allowed(
        "prepared_receptor", {}, verified=1, sha256="not-a-sha256"
    )
    with pytest.raises(ArtifactContractError, match="unknown artifact_role"):
        artifact_publication_policy("all_poses")


def test_typed_artifact_index_rejects_duplicate_archive_member() -> None:
    index = _artifact_index()
    group = index["groups"][0]  # type: ignore[index]
    group["entries"].append(deepcopy(group["entries"][0]))

    errors = validate_artifact_index(index)

    assert any("duplicates archive_path/member_name" in error for error in errors)


def test_typed_artifact_index_rejects_string_selection_marker() -> None:
    index = _artifact_index()
    entry = index["groups"][0]["entries"][0]  # type: ignore[index]
    entry["selected_for_release"] = "true"

    errors = validate_artifact_index(index)

    assert any("selected_for_release must be a boolean" in error for error in errors)


def test_attempt_selector_rejects_unsafe_or_incomplete_identity() -> None:
    manifest = _release_manifest()
    manifest["runs"] = [
        {
            "run_id": "run-a",
            "attempt_selections": [
                {
                    "pdb_id": "1ABC",
                    "variant": "HOLO",
                    "ph_label": "base",
                    "ligand_canonical_id": "drug-a",
                    "completion_relpath": "../retry/completion_stage3.json",
                    "completion_sha256": "not-a-hash",
                }
            ],
        }
    ]

    errors = validate_release_manifest(manifest)

    assert any(
        "completion_relpath must be a safe POSIX path" in error for error in errors
    )
    assert any(
        "completion_sha256 must be 64 hexadecimal digits" in error for error in errors
    )


def test_public_path_redaction_handles_provider_neutral_mounts() -> None:
    from analysis.cli.build_atlas_database import _sanitize_public_string

    repo_root = Path("/repo")
    redacted = _sanitize_public_string(
        "failure at /scratch/project/raw.log and /mnt/archive/input.sdf",
        repo_root,
    )

    assert "/scratch/" not in redacted
    assert "/mnt/" not in redacted
    assert redacted.count("external/") == 2
    assert (
        _sanitize_public_string("https://example.org/api/releases/v1", repo_root)
        == "https://example.org/api/releases/v1"
    )


def test_public_sqlite_omits_disallowed_artifacts_before_sanitization(
    tmp_path: Path,
) -> None:
    database = tmp_path / "public.sqlite"
    with sqlite3.connect(database) as connection:
        create_schema(connection)
        connection.execute(
            """INSERT INTO releases(
            release_id, schema_version, manifest_sha256, manifest_json, imported_at
            ) VALUES ('release', '1', ?, '{}', '2026-07-14T00:00:00Z')""",
            ("a" * 64,),
        )
        connection.execute(
            """INSERT INTO runs(
            run_id, release_id, manifest_path, manifest_sha256, manifest_json
            ) VALUES ('run', 'release', '/stor/private/manifest.yaml', ?, '{}')""",
            ("b" * 64,),
        )
        artifacts = [
            (
                "prepared_receptor",
                "receptor",
                "/stor/private/receptor.pdbqt",
                "/mnt/private/archive.tar.zst",
                "prepared/receptor.pdbqt",
                {"source_path": "/scratch/private/receptor.pdb"},
                1,
                "1" * 64,
            ),
            (
                "docking_pose",
                "pair",
                "/stor/private/selected_pose.pdbqt",
                "/mnt/private/archive.tar.zst",
                "poses/selected.pdbqt",
                {
                    "selected_for_release": True,
                    "source_path": "/scratch/private/selected_pose.pdbqt",
                },
                1,
                "2" * 64,
            ),
            (
                "docking_pose",
                "pair",
                "/stor/private/unselected_pose.pdbqt",
                "/mnt/private/archive.tar.zst",
                "poses/unselected.pdbqt",
                {"selected_for_release": False},
                1,
                "3" * 64,
            ),
            (
                "raw_stdout",
                "pair",
                "/stor/private/stdout.txt",
                "/mnt/private/archive.tar.zst",
                "logs/stdout.txt",
                {},
                1,
                "4" * 64,
            ),
            (
                "all_docking_poses",
                "pair",
                "/stor/private/all_poses.pdbqt",
                "/mnt/private/archive.tar.zst",
                "poses/all.pdbqt",
                {},
                1,
                "5" * 64,
            ),
            (
                "all_pose_scores",
                "pair",
                "/stor/private/all_scores.csv",
                "/mnt/private/archive.tar.zst",
                "scores/all.csv",
                {},
                1,
                "6" * 64,
            ),
            (
                "pose_image",
                "pair",
                "/stor/private/unverified.png",
                "/mnt/private/archive.tar.zst",
                "images/unverified.png",
                {},
                0,
                "7" * 64,
            ),
            (
                "pose_validation_report",
                "pair",
                "/stor/private/invalid-hash.json",
                "/mnt/private/archive.tar.zst",
                "validation/invalid-hash.json",
                {},
                1,
                "not-a-sha256",
            ),
        ]
        connection.executemany(
            """INSERT INTO artifacts(
            run_id, artifact_role, artifact_scope, original_path, archive_path,
            member_name, artifact_json, verified, sha256
            ) VALUES ('run', ?, ?, ?, ?, ?, ?, ?, ?)""",
            [
                (
                    *row[:5],
                    json.dumps(row[5], sort_keys=True),
                    row[6],
                    row[7],
                )
                for row in artifacts
            ],
        )
        connection.commit()

    counts = _sanitize_public_database(database, tmp_path)

    assert counts["artifact_rows_examined"] == 8
    assert counts["artifact_rows_retained"] == 2
    assert counts["artifact_rows_omitted"] == 6
    assert counts["artifact_private_rows_omitted"] == 3
    assert counts["artifact_conditional_rows_omitted"] == 1
    assert counts["artifact_unknown_role_rows_omitted"] == 0
    assert counts["artifact_unverified_rows_omitted"] == 1
    assert counts["artifact_invalid_sha256_rows_omitted"] == 1
    with sqlite3.connect(database) as connection:
        rows = connection.execute(
            """SELECT artifact_role, original_path, archive_path, artifact_json,
            verified, sha256
            FROM artifacts ORDER BY artifact_id"""
        ).fetchall()
    assert [row[0] for row in rows] == ["prepared_receptor", "docking_pose"]
    assert [row[4] for row in rows] == [1, 1]
    assert [row[5] for row in rows] == ["1" * 64, "2" * 64]
    serialized = json.dumps(rows)
    assert "/stor/" not in serialized
    assert "/mnt/" not in serialized
    assert "/scratch/" not in serialized
    assert json.loads(rows[1][3])["selected_for_release"] is True
