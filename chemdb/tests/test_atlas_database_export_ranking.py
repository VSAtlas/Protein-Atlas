from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

from analysis.atlas_database.exports import export_release_database
from analysis.atlas_database.readiness import audit_release_readiness
from analysis.atlas_database.pose_validation_contract import (
    QUALIFYING_POSE_VALIDATION_SCOPE,
)
from analysis.atlas_database.schema import create_schema


def _insert_pair(
    connection: sqlite3.Connection,
    context_id: int,
    ligand_id: int,
    *,
    final_score: float | None,
    pose_valid: int | None,
    atlas_score: float | None = None,
) -> None:
    result = (
        json.dumps({"final_score_source": "z_vs_decoys_blend"})
        if final_score is not None
        else "{}"
    )
    connection.execute(
        """INSERT INTO pair_cells
        (receptor_context_id, ligand_id, final_status, has_result, pose_valid,
         pose_validation_scope, final_score, atlas_score, atlas_score_source,
         result_json)
        VALUES (?, ?, 'valid', 1, ?, ?, ?, ?, ?, ?)""",
        (
            context_id,
            ligand_id,
            pose_valid,
            QUALIFYING_POSE_VALIDATION_SCOPE if pose_valid is not None else None,
            final_score,
            atlas_score,
            "decoy_standardized" if atlas_score is not None else None,
            result,
        ),
    )


def test_export_ranks_only_pose_valid_final_scores_without_fallback(
    tmp_path: Path,
) -> None:
    database = tmp_path / "release.sqlite"
    with sqlite3.connect(database) as connection:
        create_schema(connection)
        connection.execute(
            "INSERT INTO releases VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("release-test", "1", None, None, "manifest-hash", "{}", "now"),
        )
        connection.execute(
            "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                "run-1",
                "release-test",
                "completed",
                "manifests/run-1/run_manifest.yaml",
                "run-hash",
                "{}",
                "{}",
                "{}",
                "{}",
                "{}",
                "{}",
            ),
        )
        connection.executemany(
            """INSERT INTO receptor_contexts
            (receptor_context_id, run_id, pdb_id, variant, ph_label)
            VALUES (?, 'run-1', ?, 'HOLO', '7.4')""",
            [(1, "P1"), (2, "P2")],
        )
        connection.executemany(
            """INSERT INTO receptor_annotations
            (receptor_context_id, receptor_classification, qualification_status,
             native_redock_status, source_path, source_sha256, source_record_index,
             source_record_json) VALUES (?, 'HOLO', 'qualification_passed',
                    'qualification_passed', 'annotations.csv', 'source-hash', ?, '{}')""",
            [(1, 1), (2, 2)],
        )
        connection.executemany(
            "INSERT INTO ligands(ligand_id, canonical_id) VALUES (?, ?)",
            [(1, "ranked"), (2, "invalid"), (3, "unvalidated"), (4, "atlas-only")],
        )
        _insert_pair(connection, 1, 1, final_score=2.5, pose_valid=1)
        _insert_pair(connection, 2, 1, final_score=1.0, pose_valid=1)
        _insert_pair(connection, 1, 2, final_score=9.0, pose_valid=0)
        _insert_pair(connection, 1, 3, final_score=8.0, pose_valid=None)
        _insert_pair(
            connection,
            1,
            4,
            final_score=None,
            pose_valid=1,
            atlas_score=99.0,
        )
        connection.commit()

    out_dir = tmp_path / "exports"
    summary = export_release_database(database, out_dir, include_parquet=False)

    assert summary["pair_count"] == 5
    assert summary["primary_score_count"] == 4
    assert summary["rank_eligible_count"] == 2
    payload = json.loads((out_dir / "release_browser.json").read_text(encoding="utf-8"))
    assert payload["score_contract"] == {
        "primary_field": "final_score",
        "source_field": "final_score_source",
        "direction": "higher_is_better",
        "no_fallback": True,
        "normalized_comparison_field": "atlas_score",
        "normalized_comparison_source_field": "atlas_score_source",
    }
    rows = {(row["pdb_id"], row["drug_id"]): row for row in payload["pairs"]}
    high = rows[("P1", "ranked")]
    low = rows[("P2", "ranked")]
    assert high["final_score"] == 2.5
    assert high["final_score_source"] == "z_vs_decoys_blend"
    assert high["rank_eligible"] == 1
    assert high["rank_within_receptor"] == 1
    assert high["rank_across_receptors"] == 1
    assert high["score_component_coverage"] == "final_score"
    assert high["atlas_score"] is None
    assert high["atlas_score_source"] is None
    assert low["rank_within_receptor"] == 1
    assert low["rank_across_receptors"] == 2

    assert rows[("P1", "invalid")]["ranking_eligibility_reason"] == "pose_invalid"
    assert rows[("P1", "invalid")]["rank_within_receptor"] is None
    assert (
        rows[("P1", "unvalidated")]["ranking_eligibility_reason"]
        == "pose_validation_missing"
    )
    atlas_only = rows[("P1", "atlas-only")]
    assert atlas_only["ranking_eligibility_reason"] == "missing_final_score"
    assert atlas_only["primary_score_present"] == 0
    assert atlas_only["rank_within_receptor"] is None
    assert atlas_only["final_score"] is None
    assert atlas_only["atlas_score"] == 99.0
    assert atlas_only["score_component_coverage"] == "atlas_score"

    with (out_dir / "pairs.csv").open(encoding="utf-8", newline="") as handle:
        csv_rows = list(csv.DictReader(handle))
    assert len(csv_rows) == 5
    assert {row["ligand_canonical_id"] for row in csv_rows} == {
        "ranked",
        "invalid",
        "unvalidated",
        "atlas-only",
    }


def test_readiness_requires_selected_result_completion_lineage(
    tmp_path: Path,
) -> None:
    database = tmp_path / "attempt-lineage.sqlite"
    with sqlite3.connect(database) as connection:
        create_schema(connection)
        connection.execute(
            "INSERT INTO releases VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("lineage-release", "1", None, None, "hash", "{}", "now"),
        )
        connection.execute(
            "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                "run-lineage",
                "lineage-release",
                "completed",
                "manifest.yaml",
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
            VALUES (1, 'run-lineage', 'P1', 'HOLO', '7.4')"""
        )
        connection.execute(
            "INSERT INTO ligands(ligand_id, canonical_id) VALUES (1, 'drug-a')"
        )
        _insert_pair(connection, 1, 1, final_score=1.0, pose_valid=1)
        pair_cell_id = int(
            connection.execute("SELECT pair_cell_id FROM pair_cells").fetchone()[0]
        )
        connection.execute(
            """INSERT INTO completion_records
            (completion_record_id, run_id, receptor_context_id, completion_path,
             completion_sha256, completion_json)
            VALUES (1, 'run-lineage', 1, 'completion.json', ?, '{}')""",
            ("a" * 64,),
        )
        connection.execute(
            """INSERT INTO docking_attempts
            (pair_cell_id, completion_record_id, status, selected_for_release)
            VALUES (?, 1, 'completed', 1)""",
            (pair_cell_id,),
        )
        connection.execute(
            """INSERT INTO result_attempts
            (pair_cell_id, input_csv_path, input_csv_sha256, source_row_number,
             result_sha256, final_status, selected_for_release)
            VALUES (?, 'scores.csv', ?, 2, ?, 'valid', 1)""",
            (pair_cell_id, "b" * 64, "c" * 64),
        )
        connection.commit()

    readiness = audit_release_readiness(database)
    checks = {item["check"]: item for item in readiness["checklist"]}
    assert checks["selected_result_completion_lineage"]["missing_count"] == 1

    with sqlite3.connect(database) as connection:
        connection.execute(
            """UPDATE result_attempts
            SET completion_record_id=1,
                completion_link_method='completion_manifest_pair_record',
                completion_link_evidence_json=?
            WHERE selected_for_release=1""",
            (json.dumps({"record_key": "P1|HOLO|7.4|drug-a"}),),
        )
        connection.commit()

    readiness = audit_release_readiness(database)
    checks = {item["check"]: item for item in readiness["checklist"]}
    assert checks["selected_result_completion_lineage"]["missing_count"] == 0
    assert readiness["metrics"]["attempt_lineage"]["completion_linked_count"] == 1


def test_headline_ranking_requires_quality_redock_and_non_apo_context(
    tmp_path: Path,
) -> None:
    database = tmp_path / "qualified-release.sqlite"
    with sqlite3.connect(database) as connection:
        create_schema(connection)
        connection.execute(
            "INSERT INTO releases VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("qualified-release", "1", None, None, "hash", "{}", "now"),
        )
        connection.execute(
            "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                "run-qualified",
                "qualified-release",
                "completed",
                "manifest.yaml",
                "hash",
                "{}",
                "{}",
                "{}",
                "{}",
                "{}",
                "{}",
            ),
        )
        connection.executemany(
            """INSERT INTO receptor_contexts
            (receptor_context_id, run_id, pdb_id, variant, ph_label,
             native_redock_status)
            VALUES (?, 'run-qualified', ?, ?, '7.4', ?)""",
            [
                (1, "P1", "HOLO", None),
                (2, "P2", "APO", None),
                (3, "P3", "HOLO", None),
                (4, "P4", "HOLO", "observed_valid_control_result_unqualified"),
                (5, "P5", "HOLO", None),
            ],
        )
        connection.execute(
            "INSERT INTO ligands(ligand_id, canonical_id) VALUES (1, 'drug-a')"
        )
        connection.executemany(
            """INSERT INTO receptor_annotations
            (receptor_context_id, receptor_classification,
             qualification_status, native_redock_status, source_path,
             source_sha256, source_record_index, source_record_json)
            VALUES (?, ?, ?, ?, 'annotations.csv', 'hash', ?, '{}')""",
            [
                (1, "HOLO", "qualification_passed", "qualified", 1),
                (2, "APO", "qualification_passed", "qualified", 2),
                (3, "HOLO", None, "qualified", 3),
                (4, "HOLO", "qualified", None, 4),
                (5, None, "qualified", "qualified", 5),
            ],
        )
        for context_id in range(1, 6):
            _insert_pair(
                connection,
                context_id,
                1,
                final_score=float(10 - context_id),
                pose_valid=1,
            )
        connection.commit()

    out_dir = tmp_path / "qualified-exports"
    summary = export_release_database(database, out_dir, include_parquet=False)
    payload = json.loads((out_dir / "release_browser.json").read_text(encoding="utf-8"))
    rows = {row["pdb_id"]: row for row in payload["pairs"]}

    assert summary["pair_count"] == 5
    assert summary["rank_eligible_count"] == 1
    assert rows["P1"]["rank_eligible"] == 1
    assert rows["P2"]["ranking_eligibility_reason"] == "apo_receptor"
    assert rows["P3"]["ranking_eligibility_reason"] == "receptor_quality_not_qualified"
    assert rows["P4"]["ranking_eligibility_reason"] == "native_redock_not_qualified"
    assert rows["P5"]["ranking_eligibility_reason"] == "receptor_classification_missing"
    assert all(
        row["rank_within_receptor"] is None
        for row in rows.values()
        if row["pdb_id"] != "P1"
    )

    readiness = audit_release_readiness(database)
    assert readiness["metrics"]["scores"]["ranking_eligibility_reason_counts"] == {
        "apo_receptor": 1,
        "eligible": 1,
        "native_redock_not_qualified": 1,
        "receptor_classification_missing": 1,
        "receptor_quality_not_qualified": 1,
    }
    checks = {item["check"]: item for item in readiness["checklist"]}
    assert checks["receptor_classification_controlled_coverage"]["missing_count"] == 1
    assert checks["receptor_quality_explicit_qualification"]["missing_count"] == 1
    assert checks["native_redock_explicit_qualification"]["missing_count"] == 1
