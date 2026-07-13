from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

import pytest
import yaml

from analysis.atlas_database import (
    ReleaseManifestError,
    build_release_database,
    load_release_manifest,
    validate_release_manifest,
)
from analysis.cli import build_atlas_database as build_database_cli
from cli.qol.publish import _cmd_publish


def _policies() -> dict[str, dict[str, str]]:
    return {
        "receptor_selection": {"policy": "experimentally determined holo"},
        "native_redocking": {"policy": "mandatory"},
        "failure_handling": {"policy": "retain every scheduled pair"},
        "normalization": {"policy": "within-receptor z score"},
    }


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_release_manifest_validation_reports_all_structural_errors(tmp_path: Path) -> None:
    invalid = {
        "release_id": "unsafe release/id",
        "schema_version": 2,
        "runs": ["run-a", {"run_id": "run-a"}, {}],
        "scientific_policies": {
            "receptor_selection": {},
            "native_redocking": {"policy": "mandatory"},
        },
    }

    errors = validate_release_manifest(invalid)

    assert "release_id must use only letters, numbers, '.', '_', or '-'" in errors
    assert "schema_version must be 1" in errors
    assert "duplicate run_id: run-a" in errors
    assert "runs[2] must be a run ID or mapping with run_id" in errors
    assert "scientific_policies.receptor_selection must be a non-empty mapping" in errors
    assert "scientific_policies.failure_handling must be a non-empty mapping" in errors
    assert "scientific_policies.normalization must be a non-empty mapping" in errors

    path = tmp_path / "invalid_release.json"
    path.write_text(json.dumps(invalid), encoding="utf-8")
    with pytest.raises(ReleaseManifestError, match="schema_version must be 1"):
        load_release_manifest(path)


def test_legacy_ph_and_failure_complete_sqlite_build(tmp_path: Path) -> None:
    run_id = "legacy-ph-run"
    run_manifest = tmp_path / "run_manifest.yaml"
    master_rows = tmp_path / "master_rows.csv"
    docked = tmp_path / "docked"
    completion = docked / "1ABC" / "completion_stage3.json"
    completion.parent.mkdir(parents=True)

    run_manifest.write_text(
        yaml.safe_dump(
            {
                "run_id": run_id,
                "status": "completed",
                "proteins": {
                    "1ABC|HOLO|base": {
                        "pdb_id": "1ABC",
                        "variant": "HOLO",
                        "ph_label": "base",
                        "status": "completed",
                        "library": "FDA",
                        "stages": {
                            "prep": {
                                "details": {
                                    "input_pdb": "/stor/private/shared.pdb",
                                    "receptor_pdbqt": "/home/private/shared.pdb",
                                }
                            }
                        },
                    }
                },
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    expected = [
        "native_ctrl_redock_vina.pdbqt",
        "invalid_pose.pdbqt",
        "unsuccessful.pdbqt",
        "missing_output.pdbqt",
    ]
    completion.write_text(
        json.dumps(
            {
                "pdb_id": "1ABC",
                "variant": "HOLO",
                "ph_label": "completion_stage3.csv",
                "engine": "vina",
                "stage": "stage3",
                "chunk_id": "0",
                "expected_ligands": expected,
                "missing_ligands_after": ["missing_output.pdbqt"],
                "failure_markers": {
                    "missing_output": "vina output was not produced"
                },
            }
        ),
        encoding="utf-8",
    )
    _write_csv(
        master_rows,
        [
            {
                "pdb_id": "1ABC",
                "variant": "HOLO",
                "ph_label": "base",
                "ligand_base": "native_ctrl_redock_vina",
                "ligand_display": "Native control",
                "pose_valid_any": "1",
                "is_control": "1",
                "is_decoy": "0",
                "z_selected": "2.5",
                "z_selected_source": "decoy_standardized",
                "pose_invalid_reason_top": "",
            },
            {
                "pdb_id": "1ABC",
                "variant": "HOLO",
                "ph_label": "base",
                "ligand_base": "invalid_pose",
                "ligand_display": "Invalid pose",
                "pose_valid_any": "0",
                "is_control": "0",
                "is_decoy": "0",
                "z_selected": "1.2",
                "z_selected_source": "decoy_standardized",
                "pose_invalid_reason_top": "steric clash",
            },
            {
                "pdb_id": "1ABC",
                "variant": "HOLO",
                "ph_label": "base",
                "ligand_base": "unsuccessful",
                "ligand_display": "Unsuccessful",
                "pose_valid_any": "",
                "is_control": "0",
                "is_decoy": "0",
                "z_selected": "",
                "z_selected_source": "missing_consensus_decoy_null",
                "pose_invalid_reason_top": "engine returned no score",
            },
        ],
    )
    release_manifest = tmp_path / "release.yaml"
    release_manifest.write_text(
        yaml.safe_dump(
            {
                "release_id": "atlas-v0.1-test",
                "schema_version": 1,
                "runs": [
                    {
                        "run_id": run_id,
                        "paths": {
                            "run_manifest": str(run_manifest),
                            "master_rows": str(master_rows),
                            "docked": str(docked),
                        },
                    }
                ],
                "scientific_policies": _policies(),
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    database = tmp_path / "docking_atlas.sqlite"

    summary = build_release_database(release_manifest, database, tmp_path)

    assert summary["failure_complete"] is True
    assert summary["pair_status_counts"] == {
        "invalid": 1,
        "missing_result": 1,
        "unsuccessful": 1,
        "valid": 1,
    }
    with sqlite3.connect(database) as connection:
        contexts = connection.execute(
            "SELECT pdb_id, variant, ph_label, native_redock_status "
            "FROM receptor_contexts"
        ).fetchall()
        assert contexts == [
            ("1ABC", "HOLO", "", "observed_valid_control_result_unqualified")
        ]
        rows = connection.execute(
            "SELECT ligands.canonical_id, pair_cells.final_status, "
            "pair_cells.failure_code, pair_cells.failure_reason "
            "FROM pair_cells JOIN ligands USING (ligand_id) "
            "ORDER BY ligands.canonical_id"
        ).fetchall()
    assert rows == [
        ("invalid_pose", "invalid", "pose_invalid", "steric clash"),
        (
            "missing_output",
            "missing_result",
            "missing_final_result",
            "scheduled pair has no master result row",
        ),
        ("native", "valid", None, None),
        ("unsuccessful", "unsuccessful", "no_numeric_result", "engine returned no score"),
    ]

    publication_dir = tmp_path / "publication"
    publication_dir.mkdir()
    private_sentinel = publication_dir / "keep-private.txt"
    private_sentinel.write_text("keep", encoding="utf-8")
    build_args = [
        "build",
        "--manifest",
        str(release_manifest),
        "--repo-root",
        str(tmp_path),
        "--out-dir",
        str(publication_dir),
        "--no-parquet",
    ]

    assert build_database_cli.main(build_args) == 0

    site_dir = publication_dir / "site"
    downloads_dir = site_dir / "downloads"
    public_build_summary = json.loads(
        (publication_dir / "build_summary.json").read_text(encoding="utf-8")
    )
    assert public_build_summary["site_dir"] == "site"
    assert public_build_summary["rank_eligible_count"] == 0
    assert public_build_summary["readiness_blocker_count"] > 0
    assert public_build_summary["readiness_explicitly_qualified_receptor_count"] == 0
    assert public_build_summary["image_plan_context_count"] == 1
    assert public_build_summary["image_plan_selection_count"] == 1
    assert public_build_summary["image_plan_gap_count"] > 0
    assert str(tmp_path) not in json.dumps(public_build_summary)
    for relative_path in (
        "index.html",
        "analysis.html",
        "site_manifest.json",
        "assets/atlas.css",
        "assets/atlas.js",
        "assets/analysis_data.json",
        "downloads/docking_atlas.sqlite",
        "downloads/pairs.csv",
        "downloads/release_browser.json",
        "downloads/export_summary.json",
        "downloads/database_redaction.json",
        "downloads/release_readiness.json",
        "downloads/image_plan.json",
    ):
        assert (site_dir / relative_path).is_file()
    assert not (downloads_dir / "pairs.parquet").exists()

    browser_payload = json.loads(
        (downloads_dir / "release_browser.json").read_text(encoding="utf-8")
    )
    assert [item["url"] for item in browser_payload["downloads"]] == [
        "downloads/docking_atlas.sqlite",
        "downloads/pairs.csv",
        "downloads/release_readiness.json",
        "downloads/image_plan.json",
        "downloads/database_redaction.json",
    ]
    assert all(
        item["content_hash"].startswith("sha256:")
        for item in browser_payload["downloads"]
    )
    site_manifest = json.loads(
        (site_dir / "site_manifest.json").read_text(encoding="utf-8")
    )
    assert site_manifest["downloads"] == browser_payload["downloads"]
    for public_text_path in site_dir.rglob("*"):
        if public_text_path.suffix in {".html", ".json", ".csv", ".js", ".css"}:
            assert str(tmp_path) not in public_text_path.read_text(encoding="utf-8")
    with sqlite3.connect(downloads_dir / "docking_atlas.sqlite") as connection:
        public_manifest_path = connection.execute(
            "SELECT manifest_path FROM runs WHERE run_id=?", (run_id,)
        ).fetchone()[0]
        assert public_manifest_path == "run_manifest.yaml"
        for (value,) in connection.execute(
            "SELECT manifest_json FROM releases UNION ALL "
            "SELECT manifest_json FROM runs"
        ):
            assert str(tmp_path) not in value
        input_path, prepared_path = connection.execute(
            "SELECT input_receptor_path, prepared_receptor_path "
            "FROM receptor_contexts"
        ).fetchone()
        assert input_path.startswith("external/")
        assert prepared_path.startswith("external/")
        assert input_path != prepared_path

    private_database = publication_dir / "docking_atlas.private.sqlite"
    assert private_database.is_file()
    with sqlite3.connect(private_database) as connection:
        assert connection.execute(
            "SELECT input_receptor_path FROM receptor_contexts"
        ).fetchone()[0] == "/stor/private/shared.pdb"
    redaction = json.loads(
        (downloads_dir / "database_redaction.json").read_text(encoding="utf-8")
    )
    assert redaction["projection"] == "public_path_redacted_copy"
    assert redaction["source_database_sha256"].startswith("sha256:")
    assert redaction["public_database_sha256"].startswith("sha256:")
    assert redaction["source_database_sha256"] != redaction["public_database_sha256"]
    assert redaction["redaction_counts"]["path_columns_redacted"] >= 3
    assert "original release manifest" in redaction["manifest_sha256_semantics"]
    readiness = json.loads(
        (downloads_dir / "release_readiness.json").read_text(encoding="utf-8")
    )
    native_metrics = readiness["metrics"]["native_redocking"]
    assert native_metrics["status_counts"] == {
        "observed_valid_control_result_unqualified": 1
    }
    assert native_metrics["explicitly_qualified_count"] == 0
    assert native_metrics["not_explicitly_qualified_count"] == 1
    image_plan = json.loads(
        (downloads_dir / "image_plan.json").read_text(encoding="utf-8")
    )
    assert image_plan["database_file"] == "docking_atlas.sqlite"
    assert image_plan["summary"] == {
        "context_count": 1,
        "selection_count": 1,
        "gap_count": len(image_plan["contexts"][0]["gaps"]),
    }
    selection = image_plan["contexts"][0]["selections"][0]
    assert selection["role"] == "native_control"
    assert selection["command_argv"][:2] == ["atlas", "screenshot"]
    assert all(
        not Path(output_dir).is_absolute()
        for output_dir in selection["output_directories"].values()
    )

    with pytest.raises(FileExistsError, match="publication site already exists"):
        build_database_cli.main(build_args)
    stale = site_dir / "stale.txt"
    stale.write_text("remove", encoding="utf-8")
    assert build_database_cli.main([*build_args, "--overwrite"]) == 0
    assert not stale.exists()
    assert private_sentinel.read_text(encoding="utf-8") == "keep"


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (
            ["audit", "--manifest", "release.yaml", "--out-dir", "audit", "--strict"],
            ["audit", "--manifest", "release.yaml", "--out-dir", "audit", "--strict"],
        ),
        (
            [
                "build",
                "--manifest",
                "release.yaml",
                "--out-dir",
                "site",
                "--overwrite",
                "--no-parquet",
            ],
            [
                "build",
                "--manifest",
                "release.yaml",
                "--out-dir",
                "site",
                "--overwrite",
                "--no-parquet",
            ],
        ),
    ],
)
def test_publish_cli_forwards_manifest_commands(
    monkeypatch: pytest.MonkeyPatch,
    argv: list[str],
    expected: list[str],
) -> None:
    from analysis.cli import build_atlas_database

    captured: list[str] = []

    def fake_main(forwarded: list[str]) -> int:
        captured.extend(forwarded)
        return 7

    monkeypatch.setattr(build_atlas_database, "main", fake_main)

    assert _cmd_publish(argv) == 7
    assert captured == expected
