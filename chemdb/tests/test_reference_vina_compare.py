from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import pytest

from analysis.ml.reference_vina_compare import (
    _scorch_null_manifest_row,
    _scorch_null_source_paths,
)
from docking.docking_consensus_score import _external_decoy_consensus_scores
from post_docking.rescoring.rescore_reranker import _combo_scope_from_relative


def test_scorch_guard_covers_decoy_and_score_references(tmp_path: Path) -> None:
    decoy = tmp_path / "dud_consensus.csv"
    score = tmp_path / "consensus.csv"
    decoy.write_text("score\n1\n", encoding="utf-8")
    score.write_text("score\n2\n", encoding="utf-8")
    null = {
        "decoy_reference_file": json.dumps([str(decoy)]),
        "score_reference_file": json.dumps([str(score)]),
    }

    paths = _scorch_null_source_paths({"1ABC": null})

    assert paths == sorted([decoy.resolve(), score.resolve()])


def test_scorch_manifest_row_preserves_both_reference_kinds() -> None:
    null = {
        "pdb_id": "1ABC",
        "n_decoys": 200,
        "n_unique": 190,
        "mu": 0.2,
        "sigma": 0.1,
        "null_sha256": "a" * 64,
        "decoy_reference_file": json.dumps(["dud.csv"]),
        "decoy_reference_sha256": json.dumps(["b" * 64]),
        "score_reference_file": json.dumps(["fda.csv"]),
        "score_reference_sha256": json.dumps(["c" * 64]),
    }

    row = _scorch_null_manifest_row(null, stage1_decoy_n=210)

    assert row["decoy_reference_file"] == json.dumps(["dud.csv"])
    assert row["score_reference_file"] == json.dumps(["fda.csv"])
    assert row["score_reference_sha256"] == json.dumps(["c" * 64])
    assert row["stage1_coverage_fraction"] == pytest.approx(200 / 210)


def test_external_consensus_null_records_exact_source(tmp_path: Path) -> None:
    dud_path = tmp_path / "dud_consensus_docking_scores.csv"
    rows = [
        {
            "run_id": "reference",
            "is_decoy": "1",
            "consensus_score": "0.2",
        },
        {
            "run_id": "reference",
            "is_decoy": "1",
            "consensus_score": "0.8",
        },
        {
            "run_id": "reference",
            "is_decoy": "0",
            "consensus_score": "9.9",
        },
        {
            "run_id": "other",
            "is_decoy": "1",
            "consensus_score": "7.7",
        },
    ]
    with dud_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    scores, source = _external_decoy_consensus_scores(
        tmp_path,
        current_prefix="",
        expected_run_id="reference",
        logger=logging.getLogger("test"),
    )

    assert scores == [0.2, 0.8]
    assert source == str(dud_path)


@pytest.mark.parametrize(
    ("relative", "expected"),
    [
        (Path("1abc"), ("1ABC", "", "")),
        (Path("1abc/HOLO"), ("1ABC", "HOLO", "")),
        (Path("1abc/pH7_0"), ("1ABC", "", "pH7_0")),
        (Path("1abc/HOLO/pH7_0"), ("1ABC", "HOLO", "pH7_0")),
    ],
)
def test_reranker_scope_supports_all_atlas_layouts(
    relative: Path,
    expected: tuple[str, str, str],
) -> None:
    assert _combo_scope_from_relative(relative) == expected


def test_reranker_scope_rejects_unstructured_depth() -> None:
    with pytest.raises(ValueError, match="PDB.*variant.*pH"):
        _combo_scope_from_relative(Path("1ABC/HOLO/pH7_0/unexpected"))
