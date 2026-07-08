from __future__ import annotations

import pandas as pd

from analysis.dud_eval_labels import (
    _collapse_best_scores,
    _make_placeholder_row,
    parse_name_and_label,
)


def test_parse_name_and_label_handles_pose_and_microstate_suffixes() -> None:
    lig_id_a, is_active_a = parse_name_and_label(
        "x/actives_final_00012__ms_0123abcd4567ef89_pose2.pdbqt"
    )
    assert lig_id_a == "actives_final_00012"
    assert is_active_a == 1

    lig_id_d, is_active_d = parse_name_and_label("decoys_final_00044_pose3.mol2")
    assert lig_id_d == "decoys_final_00044"
    assert is_active_d == 0

    lig_id_plain, is_active_plain = parse_name_and_label("ligand_without_label.sdf")
    assert lig_id_plain == "ligand_without_label"
    assert is_active_plain is None


def test_collapse_best_scores_keeps_best_and_drops_nan() -> None:
    df = pd.DataFrame(
        {
            "lig_id": ["a", "a", "b", "b", "c"],
            "score": [3.0, 1.0, 2.0, float("nan"), float("nan")],
            "is_active": [1, 1, 0, 0, 0],
        }
    )

    best, dropped = _collapse_best_scores(df, "score", best_is_min=True)

    assert dropped == 1
    assert list(best["lig_id"]) == ["a", "b"]
    assert list(best["best_score"]) == [1.0, 2.0]
    assert list(best["is_active"]) == [1, 0]


def test_make_placeholder_row_keeps_legacy_schema() -> None:
    row = _make_placeholder_row(
        pdb_id="TEST",
        variant="HOLO",
        ph_tag="pH7_0",
        run_id="RUN1",
        bedroc_alpha=20.0,
        status_reason="missing_csv",
    )

    assert row["run_id"] == "RUN1"
    assert row["pdb_id"] == "TEST"
    assert row["variant"] == "HOLO"
    assert row["pH"] == "pH7_0"
    assert row["status_reason"] == "missing_csv"
    assert "BEDROC_alpha_20" in row.index
