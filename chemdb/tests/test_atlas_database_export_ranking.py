from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

from analysis.atlas_database.exports import export_release_database
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
         final_score, atlas_score, atlas_score_source, result_json)
        VALUES (?, ?, 'valid', 1, ?, ?, ?, ?, ?)""",
        (
            context_id,
            ligand_id,
            pose_valid,
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
    payload = json.loads(
        (out_dir / "release_browser.json").read_text(encoding="utf-8")
    )
    assert payload["score_contract"] == {
        "primary_field": "final_score",
        "source_field": "final_score_source",
        "direction": "higher_is_better",
        "no_fallback": True,
        "normalized_comparison_field": "atlas_score",
        "normalized_comparison_source_field": "atlas_score_source",
    }
    rows = {
        (row["pdb_id"], row["drug_id"]): row for row in payload["pairs"]
    }
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
