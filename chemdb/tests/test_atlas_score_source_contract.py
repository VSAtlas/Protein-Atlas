from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from analysis.atlas_database.exports import export_release_database
from analysis.atlas_database.readiness import audit_release_readiness
from analysis.atlas_database.score_source_contract import (
    AVAILABLE_COMPONENT_BLEND_VS_DECOY_Z,
    CONSENSUS_VS_DECOY_Z,
    FINAL_SCORE_SOURCE_CONTRACT_SHA256,
    FINAL_SCORE_SOURCE_CONTRACT_VERSION,
    LEGACY_NEUTRAL_IMPUTED_BLEND_VS_DECOY_Z,
    LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z,
    build_score_source_materialization,
    classify_materialized_source,
    score_source_pose_linkage_requirement,
)
from analysis.atlas_database.schema import create_schema


FILE_HASH = "a" * 64
ROW_HASH = "b" * 64


def _materialized(row: dict[str, str]) -> dict[str, object]:
    reconstructed, status, evidence = build_score_source_materialization(
        row,
        input_csv_sha256=FILE_HASH,
        source_row_number=17,
        result_sha256=ROW_HASH,
    )
    return {
        "final_score": float(row["final_score"]),
        "final_score_source": row.get("final_score_source") or None,
        "final_score_source_reconstructed": reconstructed,
        "final_score_source_classification": status,
        "final_score_source_evidence_json": evidence,
        "selected_result_input_csv_sha256": FILE_HASH,
        "selected_result_source_row_number": 17,
        "selected_result_sha256": ROW_HASH,
    }


def test_unique_documented_normalized_column_is_reconstructed_with_hashes() -> None:
    materialized = _materialized(
        {
            "final_score": "1.2500",
            "final_score_source": "",
            "z_selected_source": "stage2",
            "z_vs_decoys_blend": "1.25",
            "t_vs_decoys_blend": "1.250",
            "z_vs_decoys_consensus": "0.7",
            "scorch_composite": "0.91",
            "scorch_pct": "0.75",
            "cnn_pct": "0.5",
            "ml_blend_score": "0.6625",
            "blend_mu_decoy": "0.5",
            "blend_sigma_decoy": "0.13",
            "blend_n_decoys": "20",
        }
    )

    assert materialized["final_score_source"] is None
    assert materialized["final_score_source_reconstructed"] == (
        f"reconstructed:{LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z}"
    )
    classification = classify_materialized_source(materialized)
    assert classification.eligible is True
    assert classification.classification == "reconstructed_verified"
    assert classification.canonical_source == LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z
    assert classification.pose_linkage_requirement == (
        "scorch_rescored_stage_lineage_no_immutable_pose_identifier"
    )
    evidence = json.loads(str(materialized["final_score_source_evidence_json"]))
    fingerprint = evidence.pop("row_semantic_fingerprint_sha256")
    assert len(fingerprint) == 64
    assert evidence == {
        "ambiguity_state": "unique",
        "canonical_source": LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z,
        "contract_sha256": FINAL_SCORE_SOURCE_CONTRACT_SHA256,
        "contract_version": FINAL_SCORE_SOURCE_CONTRACT_VERSION,
        "final_score": "1.25",
        "input_csv_sha256": FILE_HASH,
        "matched_columns": ["z_vs_decoys_blend", "t_vs_decoys_blend"],
        "reason": (
            "unique_documented_normalized_source_exact_value_match_with_legacy_alias"
        ),
        "result_sha256": ROW_HASH,
        "source_row_number": 17,
        "status": "reconstructed_verified",
        "row_semantic_fingerprint_version": ("atlas_score_source_row_semantics_v1"),
    }
    assert "/" not in str(materialized["final_score_source_evidence_json"])
    assert "scorch_composite" not in str(
        materialized["final_score_source_evidence_json"]
    )


def test_pose_linkage_requirement_accepts_effective_source_prefixes() -> None:
    legacy = LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z
    expected = "scorch_rescored_stage_lineage_no_immutable_pose_identifier"

    assert score_source_pose_linkage_requirement(legacy) == expected
    assert score_source_pose_linkage_requirement(f"declared:{legacy}") == expected
    assert score_source_pose_linkage_requirement(f"reconstructed:{legacy}") == expected
    assert score_source_pose_linkage_requirement("unknown_score") is None


def test_ambiguous_equal_normalized_columns_and_raw_scorch_fail_closed() -> None:
    ambiguous = _materialized(
        {
            "final_score": "2",
            "final_score_source": "",
            "z_vs_decoys_blend": "2.0",
            "z_vs_decoys_consensus": "2.00",
        }
    )
    evidence = json.loads(str(ambiguous["final_score_source_evidence_json"]))
    assert ambiguous["final_score_source_reconstructed"] is None
    assert ambiguous["final_score_source_classification"] == (
        "reconstruction_ambiguous"
    )
    assert evidence["matched_columns"] == [
        "z_vs_decoys_consensus",
        "z_vs_decoys_blend",
    ]
    assert classify_materialized_source(ambiguous).eligible is False

    raw_scorch = _materialized(
        {
            "final_score": "0.9",
            "final_score_source": "scorch_composite",
            "scorch_composite": "0.9",
        }
    )
    classification = classify_materialized_source(raw_scorch)
    assert classification.classification == "declared_unclassified"
    assert classification.eligible is False


def test_reconstruction_requires_selected_file_and_row_provenance() -> None:
    materialized = _materialized(
        {
            "final_score": "-0.5",
            "final_score_source": "",
            "z_vs_decoys_consensus": "-0.50",
        }
    )
    assert classify_materialized_source(materialized).canonical_source == (
        CONSENSUS_VS_DECOY_Z
    )
    assert classify_materialized_source(materialized).pose_linkage_requirement == (
        "aggregate_no_unique_pose"
    )
    materialized["selected_result_sha256"] = "c" * 64
    classification = classify_materialized_source(materialized)
    assert classification.eligible is False
    assert classification.classification == "reconstruction_evidence_invalid"
    assert "row hash does not match" in classification.reason


def test_controlled_declared_source_requires_bound_exact_value_evidence() -> None:
    materialized = _materialized(
        {
            "final_score": "1.5",
            "final_score_source": "z_vs_decoys_consensus",
            "z_vs_decoys_consensus": "1.5",
        }
    )
    assert materialized["final_score_source_reconstructed"] is None
    classification = classify_materialized_source(materialized)
    assert classification.effective_source == f"declared:{CONSENSUS_VS_DECOY_Z}"
    assert classification.classification == "declared_classified"
    assert classification.eligible is True

    materialized["final_score"] = 0.0
    classification = classify_materialized_source(materialized)
    assert classification.eligible is False
    assert classification.classification == "declared_evidence_invalid"
    assert "final score does not match" in classification.reason

    mismatch = _materialized(
        {
            "final_score": "1.5",
            "final_score_source": "z_vs_decoys_consensus",
            "z_vs_decoys_consensus": "0.25",
            "z_vs_decoys_blend": "1.5",
        }
    )
    assert mismatch["final_score_source_classification"] == ("declared_value_mismatch")
    classification = classify_materialized_source(mismatch)
    assert classification.eligible is False
    assert classification.classification == "declared_value_mismatch"


def test_declared_blend_revision_requires_matching_row_semantics() -> None:
    source_row = {
        "final_score": "1.25",
        "final_score_source": AVAILABLE_COMPONENT_BLEND_VS_DECOY_Z,
        "z_vs_decoys_blend": "1.25",
        "t_vs_decoys_blend": "1.250",
        "ml_blend_mode": "scorch_only",
        "scorch_pct": "0.8",
        "ml_blend_score": "0.8",
        "ml_blend_scorch_weight_effective": "1.0",
        "ml_blend_cnn_weight_effective": "0.0",
        "blend_mu_decoy": "0.55",
        "blend_sigma_decoy": "0.2",
        "blend_n_decoys": "20",
    }
    materialized = _materialized(source_row)
    assert materialized["final_score_source_classification"] == ("declared_classified")
    classification = classify_materialized_source(materialized)
    assert classification.eligible is True
    assert classification.canonical_source == AVAILABLE_COMPONENT_BLEND_VS_DECOY_Z
    evidence = json.loads(str(materialized["final_score_source_evidence_json"]))
    assert evidence["canonical_source"] == AVAILABLE_COMPONENT_BLEND_VS_DECOY_Z
    assert evidence["matched_columns"] == [
        "z_vs_decoys_blend",
        "t_vs_decoys_blend",
    ]
    assert len(evidence["row_semantic_fingerprint_sha256"]) == 64

    tampered = dict(materialized)
    evidence["row_semantic_fingerprint_sha256"] = "c" * 64
    tampered["final_score_source_evidence_json"] = json.dumps(evidence)
    classification = classify_materialized_source(tampered)
    assert classification.eligible is False
    assert "does not verify the exact source revision" in classification.reason

    wrong_revision_row = dict(source_row)
    wrong_revision_row["final_score_source"] = LEGACY_NEUTRAL_IMPUTED_BLEND_VS_DECOY_Z
    wrong_revision = _materialized(wrong_revision_row)
    assert wrong_revision["final_score_source_classification"] == (
        "declared_semantics_unverified"
    )
    assert classify_materialized_source(wrong_revision).eligible is False


def test_reconstruction_evidence_binds_exact_canonical_revision() -> None:
    materialized = _materialized(
        {
            "final_score": "1.25",
            "final_score_source": "",
            "z_vs_decoys_blend": "1.25",
            "scorch_composite": "0.91",
            "scorch_pct": "0.75",
            "cnn_pct": "0.5",
            "ml_blend_score": "0.6625",
            "blend_mu_decoy": "0.5",
            "blend_sigma_decoy": "0.13",
            "blend_n_decoys": "20",
        }
    )
    assert classify_materialized_source(materialized).canonical_source == (
        LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z
    )

    materialized["final_score_source_reconstructed"] = (
        f"reconstructed:{LEGACY_NEUTRAL_IMPUTED_BLEND_VS_DECOY_Z}"
    )
    classification = classify_materialized_source(materialized)
    assert classification.eligible is False
    assert classification.classification == "reconstruction_evidence_invalid"
    assert "exact canonical source" in classification.reason


def test_export_keeps_declared_source_blank_and_blocks_unlinked_reconstruction(
    tmp_path: Path,
) -> None:
    database = tmp_path / "release.sqlite"
    source_row = {
        "final_score": "2.25",
        "final_score_source": "",
        "z_vs_decoys_blend": "2.250",
        "z_vs_decoys_consensus": "0.5",
        "scorch_composite": "0.9",
        "scorch_pct": "0.9",
        "cnn_pct": "0.5",
        "ml_blend_score": "0.76",
        "blend_mu_decoy": "0.4",
        "blend_sigma_decoy": "0.16",
        "blend_n_decoys": "20",
    }
    materialized = _materialized(source_row)
    with sqlite3.connect(database) as connection:
        create_schema(connection)
        connection.execute(
            "INSERT INTO releases VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("release-source", "1", None, None, "hash", "{}", "now"),
        )
        connection.execute(
            "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                "run-source",
                "release-source",
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
            VALUES (1, 'run-source', '1ABC', 'HOLO', '7.4')"""
        )
        connection.execute(
            """INSERT INTO receptor_annotations
            (receptor_context_id, receptor_classification, qualification_status,
             native_redock_status, classification_method,
             chemistry_evidence_status, prepared_receptor_sha256,
             canonical_chemistry_policy_sha256, source_path, source_sha256,
             source_record_index, source_record_json)
            VALUES (1, 'HOLO', 'qualification_passed', 'qualification_passed',
                    'atlas_prepared_receptor_retained_chemistry_v2', 'observed',
                    ?, ?, 'annotations.csv', 'hash', 1, '{}')""",
            ("c" * 64, "d" * 64),
        )
        connection.execute(
            "INSERT INTO ligands(ligand_id, canonical_id) VALUES (1, 'drug-a')"
        )
        connection.execute(
            """INSERT INTO pair_cells
            (pair_cell_id, receptor_context_id, ligand_id, final_status,
             has_result, pose_valid, pose_validation_scope, final_score,
             final_score_source, final_score_source_reconstructed,
             final_score_source_classification,
             final_score_source_evidence_json)
            VALUES (1, 1, 1, 'valid', 1, 1, 'selected_final_score_pose',
                    ?, NULL, ?, ?, ?)""",
            (
                materialized["final_score"],
                materialized["final_score_source_reconstructed"],
                materialized["final_score_source_classification"],
                materialized["final_score_source_evidence_json"],
            ),
        )
        connection.execute(
            """INSERT INTO result_attempts
            (pair_cell_id, input_csv_path, input_csv_sha256, source_row_number,
             result_sha256, final_status, final_score,
             final_score_source_reconstructed,
             final_score_source_classification,
             final_score_source_evidence_json, selected_for_release)
            VALUES (1, 'master_rows.csv', ?, 17, ?, 'valid', ?, ?, ?, ?, 1)""",
            (
                FILE_HASH,
                ROW_HASH,
                materialized["final_score"],
                materialized["final_score_source_reconstructed"],
                materialized["final_score_source_classification"],
                materialized["final_score_source_evidence_json"],
            ),
        )
        connection.commit()

    summary = export_release_database(
        database, tmp_path / "exports", include_parquet=False
    )
    payload = json.loads(
        (tmp_path / "exports" / "release_browser.json").read_text(encoding="utf-8")
    )
    pair = payload["pairs"][0]
    assert pair["final_score_source"] is None
    assert pair["final_score_source_effective"] == (
        f"reconstructed:{LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z}"
    )
    assert pair["final_score_source_classification"] == "reconstructed_verified"
    assert pair["final_score_pose_linkage_requirement"] == (
        "scorch_rescored_stage_lineage_no_immutable_pose_identifier"
    )
    assert pair["final_score_source_eligible"] == 1
    assert pair["rank_eligible"] == 0
    assert pair["ranking_eligibility_reason"] == ("final_score_pose_linkage_unresolved")
    assert pair["rank_within_receptor"] is None
    assert summary["classified_final_score_source_count"] == 1
    readiness = audit_release_readiness(database)
    assert readiness["readiness_schema_version"] == 6
    score_metrics = readiness["metrics"]["scores"]
    assert score_metrics["final_score_source_counts"] == {"missing": 1}
    assert score_metrics["effective_final_score_source_counts"] == {
        f"reconstructed:{LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z}": 1
    }
    assert score_metrics["final_score_source_classification_counts"] == {
        "reconstructed_verified": 1
    }
    assert score_metrics["classified_final_score_source_count"] == 1
    checks = {item["check"]: item for item in readiness["checklist"]}
    assert checks["final_score_source_classification"]["status"] == "complete"
