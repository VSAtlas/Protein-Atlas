from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

from analysis.atlas_database.exports import export_release_database
from analysis.atlas_database.receptor_evidence import CHEMISTRY_CLASSIFICATION_METHOD
from analysis.atlas_database.readiness import audit_release_readiness
from analysis.atlas_database.pose_validation_contract import (
    QUALIFYING_POSE_VALIDATION_SCOPE,
)
from analysis.atlas_database.score_source_contract import (
    build_score_source_materialization,
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
    with_result_attempt: bool = True,
) -> None:
    source_row_number = context_id * 1_000 + ligand_id
    result_sha256 = f"{source_row_number:064x}"
    source = "z_vs_decoys_consensus" if final_score is not None else None
    materialized = (
        build_score_source_materialization(
            {
                "final_score": str(final_score),
                "final_score_source": source,
                "z_vs_decoys_consensus": str(final_score),
            },
            input_csv_sha256="a" * 64,
            source_row_number=source_row_number,
            result_sha256=result_sha256,
        )
        if final_score is not None
        else (None, None, None)
    )
    _, classification, evidence = materialized
    result = json.dumps({"final_score_source": source}) if source else "{}"
    pair_id = connection.execute(
        """INSERT INTO pair_cells
        (receptor_context_id, ligand_id, final_status, has_result, pose_valid,
         pose_validation_scope, final_score, atlas_score, atlas_score_source,
         final_score_source, final_score_source_classification,
         final_score_source_evidence_json, result_json)
        VALUES (?, ?, 'valid', 1, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            context_id,
            ligand_id,
            pose_valid,
            QUALIFYING_POSE_VALIDATION_SCOPE if pose_valid is not None else None,
            final_score,
            atlas_score,
            "decoy_standardized" if atlas_score is not None else None,
            source,
            classification,
            evidence,
            result,
        ),
    ).lastrowid
    if final_score is not None and with_result_attempt:
        connection.execute(
            """INSERT INTO result_attempts
            (pair_cell_id, input_csv_path, input_csv_sha256, source_row_number,
             result_sha256, final_status, final_score, final_score_source,
             final_score_source_classification, final_score_source_evidence_json,
             selected_for_release) VALUES (?, 'scores.csv', ?, ?, ?, 'valid', ?, ?, ?, ?, 1)""",
            (
                pair_id,
                "a" * 64,
                source_row_number,
                result_sha256,
                final_score,
                source,
                classification,
                evidence,
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
             native_redock_status, classification_method, chemistry_evidence_status,
             prepared_receptor_sha256, canonical_chemistry_policy_sha256,
             source_path, source_sha256, source_record_index, source_record_json)
            VALUES (?, 'HOLO', 'qualification_passed', 'qualification_passed',
                    ?, 'observed', ?, ?, 'annotations.csv', 'source-hash', ?, '{}')""",
            [
                (1, CHEMISTRY_CLASSIFICATION_METHOD, "c" * 64, "d" * 64, 1),
                (2, CHEMISTRY_CLASSIFICATION_METHOD, "e" * 64, "f" * 64, 2),
            ],
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
        public_pair_id = int(
            connection.execute(
                """SELECT pair_cell_id FROM pair_cells
                WHERE receptor_context_id=1 AND ligand_id=1"""
            ).fetchone()[0]
        )
        connection.executemany(
            """INSERT INTO artifacts(
            run_id, receptor_context_id, pair_cell_id, ligand_id,
            artifact_role, artifact_scope, member_name, sha256, verified,
            artifact_json)
            VALUES ('run-1', 1, ?, 1, ?, 'pair', ?, ?, ?, ?)""",
            [
                (
                    public_pair_id,
                    "docking_pose",
                    "selected.sdf",
                    "1" * 64,
                    1,
                    '{"selected_for_release":true}',
                ),
                (
                    public_pair_id,
                    "docking_pose",
                    "unselected.sdf",
                    "2" * 64,
                    1,
                    '{"selected_for_release":false}',
                ),
                (public_pair_id, "raw_stdout", "stdout.txt", "3" * 64, 1, "{}"),
                (public_pair_id, "pose_image", "image.png", "4" * 64, 0, "{}"),
                (public_pair_id, "pose_validation_report", "report.json", "bad", 1, "{}"),
            ],
        )
        connection.commit()

    out_dir = tmp_path / "exports"
    summary = export_release_database(database, out_dir, include_parquet=False)

    assert summary["pair_count"] == 5
    assert summary["primary_score_count"] == 4
    assert summary["rank_eligible_count"] == 0
    payload = json.loads((out_dir / "release_browser.json").read_text(encoding="utf-8"))
    score_contract = payload["score_contract"]
    assert score_contract["primary_field"] == "final_score"
    assert score_contract["declared_source_field"] == "final_score_source"
    assert score_contract["reconstructed_source_field"] == (
        "final_score_source_reconstructed"
    )
    assert score_contract["effective_source_field"] == ("final_score_source_effective")
    assert score_contract["direction"] == "higher_is_better"
    assert score_contract["no_fallback"] is True
    assert score_contract["normalized_comparison_field"] == "atlas_score"
    assert score_contract["normalized_comparison_source_field"] == (
        "atlas_score_source"
    )
    assert score_contract["source_classification"]["contract_version"] == (
        "atlas_final_score_source_v1"
    )
    rows = {(row["pdb_id"], row["drug_id"]): row for row in payload["pairs"]}
    high = rows[("P1", "ranked")]
    low = rows[("P2", "ranked")]
    assert high["final_score"] == 2.5
    assert high["final_score_source"] == "z_vs_decoys_consensus"
    assert high["final_score_source_effective"] == ("declared:consensus_vs_decoy_z")
    assert high["final_score_source_classification"] == "declared_classified"
    assert high["final_score_source_eligible"] == 1
    assert high["rank_eligible"] == 0
    assert high["ranking_eligibility_reason"] == "final_score_pose_linkage_unresolved"
    assert high["rank_within_receptor"] is None
    assert high["rank_across_receptors"] is None
    assert high["score_component_coverage"] == "final_score"
    assert high["atlas_score"] is None
    assert high["atlas_score_source"] is None
    assert high["artifact_count"] == 1
    assert high["verified_artifact_count"] == 1
    assert [row["artifact_role"] for row in payload["artifacts"]] == ["docking_pose"]
    assert low["rank_within_receptor"] is None
    assert low["rank_across_receptors"] is None

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
        _insert_pair(
            connection, 1, 1, final_score=1.0, pose_valid=1, with_result_attempt=False
        )
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
             qualification_status, native_redock_status, classification_method,
             chemistry_evidence_status, prepared_receptor_sha256,
             canonical_chemistry_policy_sha256, source_path, source_sha256,
             source_record_index, source_record_json)
            VALUES (?, ?, ?, ?, ?, 'observed', ?, ?,
                    'annotations.csv', 'hash', ?, '{}')""",
            [
                (
                    1,
                    "HOLO",
                    "qualification_passed",
                    "qualified",
                    CHEMISTRY_CLASSIFICATION_METHOD,
                    "a" * 64,
                    "b" * 64,
                    1,
                ),
                (
                    2,
                    "APO",
                    "qualification_passed",
                    "qualified",
                    CHEMISTRY_CLASSIFICATION_METHOD,
                    "c" * 64,
                    "d" * 64,
                    2,
                ),
                (
                    3,
                    "HOLO",
                    None,
                    "qualified",
                    CHEMISTRY_CLASSIFICATION_METHOD,
                    "e" * 64,
                    "f" * 64,
                    3,
                ),
                (
                    4,
                    "HOLO",
                    "qualified",
                    None,
                    CHEMISTRY_CLASSIFICATION_METHOD,
                    "1" * 64,
                    "2" * 64,
                    4,
                ),
                (
                    5,
                    None,
                    "qualified",
                    "qualified",
                    CHEMISTRY_CLASSIFICATION_METHOD,
                    "3" * 64,
                    "4" * 64,
                    5,
                ),
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
    assert summary["rank_eligible_count"] == 0
    assert summary["apo_exploratory_rank_eligible_count"] == 0
    assert rows["P1"]["rank_eligible"] == 0
    assert (
        rows["P1"]["ranking_eligibility_reason"]
        == "final_score_pose_linkage_unresolved"
    )
    assert rows["P1"]["apo_exploratory_ranking_eligibility_reason"] == ("holo_receptor")
    assert rows["P2"]["ranking_eligibility_reason"] == "apo_receptor"
    assert rows["P2"]["apo_exploratory_rank_eligible"] == 0
    assert (
        rows["P2"]["apo_exploratory_ranking_eligibility_reason"]
        == "final_score_pose_linkage_unresolved"
    )
    assert rows["P2"]["rank_within_receptor"] is None
    assert rows["P2"]["apo_rank_within_receptor"] is None
    assert rows["P2"]["apo_rank_across_receptors"] is None
    assert rows["P3"]["ranking_eligibility_reason"] == "receptor_quality_not_qualified"
    assert rows["P4"]["ranking_eligibility_reason"] == "native_redock_not_qualified"
    assert rows["P5"]["ranking_eligibility_reason"] == "receptor_classification_missing"
    assert all(
        row["rank_within_receptor"] is None
        for row in rows.values()
        if row["pdb_id"] != "P1"
    )

    readiness = audit_release_readiness(database)
    score_metrics = readiness["metrics"]["scores"]
    assert score_metrics["ranking_eligibility_reason_counts"] == {
        "apo_receptor": 1,
        "final_score_pose_linkage_unresolved": 1,
        "native_redock_not_qualified": 1,
        "receptor_classification_missing": 1,
        "receptor_quality_not_qualified": 1,
    }
    assert score_metrics["apo_exploratory_rank_eligible_count"] == 0
    assert score_metrics["apo_exploratory_ranking_eligibility_reason_counts"] == {
        "final_score_pose_linkage_unresolved": 1,
        "holo_receptor": 3,
        "receptor_classification_missing": 1,
    }
    assert score_metrics["final_score_pose_linkage_unresolved_count"] == 5
    checks = {item["check"]: item for item in readiness["checklist"]}
    assert checks["receptor_classification_controlled_coverage"]["missing_count"] == 1
    assert checks["receptor_quality_explicit_qualification"]["missing_count"] == 1
    assert checks["native_redock_explicit_qualification"]["missing_count"] == 1
