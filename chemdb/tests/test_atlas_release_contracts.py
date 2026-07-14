from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from analysis.atlas_database.artifact_contract import validate_artifact_index
from analysis.atlas_database.manifest import validate_release_manifest


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


def test_typed_artifact_index_rejects_duplicate_archive_member() -> None:
    index = _artifact_index()
    group = index["groups"][0]  # type: ignore[index]
    group["entries"].append(deepcopy(group["entries"][0]))

    errors = validate_artifact_index(index)

    assert any("duplicates archive_path/member_name" in error for error in errors)


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
