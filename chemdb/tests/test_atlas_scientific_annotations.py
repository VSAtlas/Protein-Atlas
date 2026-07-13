from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

import pytest
import yaml

from analysis.atlas_database.annotations import (
    ScientificAnnotationError,
    ingest_scientific_annotations,
)
from analysis.atlas_database.image_plan import build_release_image_plan
from analysis.atlas_database.manifest import validate_release_manifest
from analysis.atlas_database.exports import export_release_database
from analysis.atlas_database.readiness import audit_release_readiness
from analysis.atlas_database.schema import create_schema


def _seed_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    create_schema(connection)
    connection.execute(
        "INSERT INTO releases VALUES (?, ?, ?, ?, ?, ?, ?)",
        ("release-annotations", "1", None, None, "hash", "{}", "now"),
    )
    connection.execute(
        "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (
            "run-1",
            "release-annotations",
            "completed",
            "manifests/run-1/run_manifest.yaml",
            "hash",
            "{}",
            "{}",
            "{}",
            "{}",
            "{}",
            "{}",
        ),
    )
    connection.execute(
        """INSERT INTO receptor_contexts
        (receptor_context_id, run_id, pdb_id, variant, ph_label)
        VALUES (1, 'run-1', '1ABC', 'APO', '7.4')"""
    )
    for ligand_id in range(1, 7):
        connection.execute(
            "INSERT INTO ligands(ligand_id, canonical_id) VALUES (?, ?)",
            (ligand_id, f"drug-{ligand_id}"),
        )
        connection.execute(
            """INSERT INTO pair_cells
            (receptor_context_id, ligand_id, final_status, has_result,
             pose_valid, final_score)
            VALUES (1, ?, 'valid', 1, 1, ?)""",
            (ligand_id, float(7 - ligand_id)),
        )
    connection.commit()
    return connection


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _base_manifest() -> dict[str, object]:
    return {
        "release_id": "release-annotations",
        "schema_version": 1,
        "runs": ["run-1"],
        "scientific_policies": {
            "receptor_selection": {"policy": "explicit"},
            "native_redocking": {"policy": "explicit"},
            "failure_handling": {"policy": "complete"},
            "normalization": {"policy": "declared"},
        },
    }


def test_ingests_explicit_annotations_with_provenance_and_no_inference(
    tmp_path: Path,
) -> None:
    database = tmp_path / "annotations.sqlite"
    connection = _seed_database(database)
    connection.execute(
        "INSERT INTO ligands(ligand_id, canonical_id) VALUES (7, 'invalid-high')"
    )
    connection.execute(
        """INSERT INTO pair_cells
        (receptor_context_id, ligand_id, final_status, has_result, pose_valid,
         final_score) VALUES (1, 7, 'invalid', 1, 0, 999.0)"""
    )
    proteins = tmp_path / "proteins.yaml"
    receptors = tmp_path / "receptors.csv"
    known_pairs = tmp_path / "known_pairs.yaml"
    manifest_path = tmp_path / "release.yaml"
    proteins.write_text(
        yaml.safe_dump(
            {
                "records": [
                    {
                        "protein_key": "protein-alpha",
                        "uniprot_id": "P12345",
                        "gene_symbol": "GENE1",
                        "display_name": "Protein Alpha",
                        "provenance": {
                            "source_record": "curation-row-7",
                            "reviewed_by": "scientist",
                        },
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    _write_csv(
        receptors,
        [
            {
                "run_id": "run-1",
                "pdb_id": "1ABC",
                "variant": "APO",
                "ph_label": "7.4",
                "protein_key": "protein-alpha",
                "receptor_classification": "experimentally_holo",
                "qualification_status": "selected_under_policy",
                "qualification_reason": "explicit curator decision",
                "native_redock_status": "qualification_passed",
                "native_redock_reason": "RMSD met the frozen redock threshold",
                "native_redock_rmsd": 1.25,
            }
        ],
    )
    known_pairs.write_text(
        yaml.safe_dump(
            [
                {
                    "run_id": "run-1",
                    "pdb_id": "1ABC",
                    "variant": "APO",
                    "ph_label": "7.4",
                    "ligand_canonical_id": "drug-6",
                    "selection_label": "representative literature pair",
                    "evidence_reference": "PMID:123456",
                    "provenance": {"curation_date": "2026-07-13"},
                }
            ]
        ),
        encoding="utf-8",
    )
    manifest = _base_manifest()
    manifest["annotations"] = {
        "proteins": proteins.name,
        "receptors": {"path": receptors.name},
        "known_pairs": known_pairs.name,
    }
    manifest_path.write_text(yaml.safe_dump(manifest), encoding="utf-8")

    counts = ingest_scientific_annotations(connection, manifest, manifest_path)
    connection.commit()

    assert counts == {"proteins": 1, "receptors": 1, "known_pairs": 1}
    protein = connection.execute(
        """SELECT protein_key, uniprot_id, gene_symbol, display_name,
        source_path, source_sha256, source_record_index, provenance_json
        FROM protein_identities"""
    ).fetchone()
    assert protein[:4] == ("protein-alpha", "P12345", "GENE1", "Protein Alpha")
    assert protein[4] == str(proteins.resolve())
    assert len(protein[5]) == 64
    assert protein[6] == 1
    assert json.loads(protein[7])["reviewed_by"] == "scientist"

    context_variant, context_redock = connection.execute(
        "SELECT variant, native_redock_status FROM receptor_contexts"
    ).fetchone()
    assert (context_variant, context_redock) == ("APO", None)
    receptor = connection.execute(
        """SELECT p.protein_key, a.receptor_classification,
        a.qualification_status, a.qualification_reason, a.native_redock_status,
        a.native_redock_reason, a.native_redock_rmsd
        FROM receptor_annotations a JOIN protein_identities p
        ON p.protein_identity_id=a.protein_identity_id"""
    ).fetchone()
    assert receptor == (
        "protein-alpha",
        "experimentally_holo",
        "selected_under_policy",
        "explicit curator decision",
        "qualification_passed",
        "RMSD met the frozen redock threshold",
        1.25,
    )
    known = connection.execute(
        """SELECT ligand_canonical_id, pair_cell_id, evidence_reference,
        source_path, provenance_json FROM known_pair_selections"""
    ).fetchone()
    assert known[0] == "drug-6"
    assert known[1] is not None
    assert known[2] == "PMID:123456"
    assert known[3] == str(known_pairs.resolve())
    assert json.loads(known[4]) == {"curation_date": "2026-07-13"}
    connection.close()

    plan = build_release_image_plan(database)
    roles = [item["role"] for item in plan["contexts"][0]["selections"]]
    assert roles == [
        "top_valid_1",
        "top_valid_2",
        "top_valid_3",
        "top_valid_4",
        "top_valid_5",
        "known_pair",
    ]
    assert "best_invalid" not in roles
    assert "invalid-high" not in {
        item["ligand_id"] for item in plan["contexts"][0]["selections"]
    }
    assert "never selected" in plan["selection_policy"]["invalid_poses"]
    assert "pair_artifact_not_indexed" in {
        gap["code"] for gap in plan["contexts"][0]["gaps"]
    }
    readiness = audit_release_readiness(database)
    assert readiness["metrics"]["native_redocking"]["explicitly_qualified_count"] == 1
    export_dir = tmp_path / "exports"
    export_release_database(database, export_dir, include_parquet=False)
    browser = json.loads(
        (export_dir / "release_browser.json").read_text(encoding="utf-8")
    )
    target = browser["targets"][0]
    assert target["protein_key"] == "protein-alpha"
    assert target["uniprot_id"] == "P12345"
    assert target["gene_symbol"] == "GENE1"
    assert target["display_name"] == "Protein Alpha"
    assert target["receptor_classification"] == "experimentally_holo"
    assert target["qualification_status"] == "selected_under_policy"
    assert target["native_redock_status"] == "qualification_passed"
    assert target["native_redock_reason"] == "RMSD met the frozen redock threshold"


def test_duplicate_or_incomplete_context_annotations_fail_without_fallback(
    tmp_path: Path,
) -> None:
    database = tmp_path / "duplicates.sqlite"
    connection = _seed_database(database)
    manifest_path = tmp_path / "release.yaml"
    duplicates = tmp_path / "duplicate_receptors.yaml"
    record = {
        "run_id": "run-1",
        "pdb_id": "1ABC",
        "variant": "APO",
        "ph_label": "7.4",
        "qualification_status": "reviewed",
    }
    duplicates.write_text(yaml.safe_dump([record, dict(record)]), encoding="utf-8")
    manifest = _base_manifest()
    manifest["annotations"] = {"receptors": duplicates.name}

    with pytest.raises(
        ScientificAnnotationError, match="duplicate annotations.receptors"
    ):
        ingest_scientific_annotations(connection, manifest, manifest_path)
    assert (
        connection.execute("SELECT COUNT(*) FROM receptor_annotations").fetchone()[0]
        == 0
    )

    incomplete = tmp_path / "incomplete.csv"
    _write_csv(
        incomplete,
        [
            {
                "run_id": "run-1",
                "pdb_id": "1ABC",
                "qualification_status": "reviewed",
            }
        ],
    )
    manifest["annotations"] = {"receptors": incomplete.name}
    with pytest.raises(ScientificAnnotationError, match="requires explicit variant"):
        ingest_scientific_annotations(connection, manifest, manifest_path)
    negative_rmsd = tmp_path / "negative_rmsd.yaml"
    negative_rmsd.write_text(
        yaml.safe_dump(
            [
                {
                    **record,
                    "native_redock_rmsd": -0.1,
                }
            ]
        ),
        encoding="utf-8",
    )
    manifest["annotations"] = {"receptors": negative_rmsd.name}
    with pytest.raises(
        ScientificAnnotationError, match="non-negative native_redock_rmsd"
    ):
        ingest_scientific_annotations(connection, manifest, manifest_path)
    connection.close()


def test_release_manifest_rejects_unknown_or_pathless_annotation_sources() -> None:
    manifest = _base_manifest()
    manifest["annotations"] = {
        "proteins": {"path": ""},
        "automatically_inferred_targets": "guesses.csv",
    }

    errors = validate_release_manifest(manifest)

    assert "annotations.proteins must be a path string or mapping with path" in errors
    assert "unknown annotations source: automatically_inferred_targets" in errors


def test_schema_rejects_older_version_before_mutating_database(
    tmp_path: Path,
) -> None:
    database = tmp_path / "legacy.sqlite"
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE legacy_table(value TEXT)")
        connection.execute("PRAGMA user_version = 1")

        with pytest.raises(
            ValueError, match="unsupported Atlas database schema version 1"
        ):
            create_schema(connection)

        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
    assert tables == {"legacy_table"}


def test_frozen_image_plan_ignores_inline_known_pair_without_annotation_record(
    tmp_path: Path,
) -> None:
    database = tmp_path / "inline-known-pair.sqlite"
    connection = _seed_database(database)
    connection.execute(
        "UPDATE releases SET manifest_json=?",
        (
            json.dumps(
                {
                    "image_plan": {
                        "known_pairs": [
                            {
                                "pdb_id": "1ABC",
                                "canonical_id": "drug-6",
                            }
                        ]
                    }
                }
            ),
        ),
    )
    connection.commit()
    connection.close()

    plan = build_release_image_plan(database)

    assert "known_pair" not in {
        item["role"] for item in plan["contexts"][0]["selections"]
    }
