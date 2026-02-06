from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.dud_eval_schema import (
    filter_df_by_run_id,
    guess_consensus_ligfile_col,
    guess_consensus_score_col,
    guess_ligfile_col,
    guess_reranked_scorch_score_col,
    guess_score_col,
    parse_valid_mask,
)


def test_guess_columns_for_docking_and_consensus_and_reranked() -> None:
    dock_df = pd.DataFrame(
        {
            "Ligand_Path": ["a.pdbqt"],
            "Docking_Score": [-7.0],
        }
    )
    assert guess_ligfile_col(dock_df, None) == "Ligand_Path"
    assert guess_score_col(dock_df, None) == "Docking_Score"

    consensus_df = pd.DataFrame(
        {
            "ligand": ["a.pdbqt"],
            "consensus_score": [0.8],
        }
    )
    assert guess_consensus_ligfile_col(consensus_df, None) == "ligand"
    assert guess_consensus_score_col(consensus_df, None) == "consensus_score"

    reranked_df = pd.DataFrame(
        {
            "ligand": ["a.pdbqt"],
            "final_score": [0.9],
            "final_rank": [1],
        }
    )
    assert guess_reranked_scorch_score_col(reranked_df, None) == "final_score"


def test_parse_valid_mask_for_numeric_and_tokens() -> None:
    numeric = pd.Series([1, 0, -2, 3, None])
    assert parse_valid_mask(numeric).tolist() == [True, False, True, True, False]

    text = pd.Series(["true", "1", "False", "yes", "0", "", None])
    assert parse_valid_mask(text).tolist() == [True, True, False, True, False, False, False]


def test_filter_df_by_run_id_handles_match_and_fallback_modes() -> None:
    df = pd.DataFrame(
        {
            "run_id": ["RUN-1", "RUN_1", "RUN2"],
            "score": [1.0, 2.0, 3.0],
        }
    )

    matched, status = filter_df_by_run_id(
        df,
        "RUN1",
        pdb_id="TEST",
        csv_path=Path("docking_score_long.csv"),
        strict=False,
    )
    assert status == "ok"
    assert len(matched) == 2

    fallback, fallback_status = filter_df_by_run_id(
        df,
        "RUN3",
        pdb_id="TEST",
        csv_path=Path("docking_score_long.csv"),
        strict=False,
    )
    assert fallback_status == "run_id_mismatch_legacy_fallback"
    assert len(fallback) == 2

    strict_df, strict_status = filter_df_by_run_id(
        df,
        "RUN3",
        pdb_id="TEST",
        csv_path=Path("docking_score_long.csv"),
        strict=True,
    )
    assert strict_status == "run_id_mismatch_new_layout"
    assert strict_df.empty

    no_run_df = pd.DataFrame({"score": [1.0, 2.0]})
    no_col_out, no_col_status = filter_df_by_run_id(
        no_run_df,
        "RUN1",
        pdb_id="TEST",
        csv_path=Path("docking_score_long.csv"),
        strict=False,
    )
    assert no_col_status == "ok_no_run_id_col"
    assert len(no_col_out) == len(no_run_df)
