from __future__ import annotations

import pandas as pd
import pytest

from analysis.ml.build_ml_dataset import add_standard_ml_columns
from analysis.ml.source_benchmark_tables import _spd_table


def test_final_score_requires_and_preserves_accepted_row_source() -> None:
    frame = pd.DataFrame(
        {
            "final_score": [1.0, 2.0, 3.0, 4.0, 5.0],
            "final_score_source": [
                "z_vs_decoys_consensus",
                "z_vs_compare_run_consensus",
                "z_vs_decoys_blend",
                "final_score",
                "",
            ],
        }
    )

    result = add_standard_ml_columns(frame)

    assert result.loc[:2, "atlas_score"].tolist() == pytest.approx([1.0, 2.0, 3.0])
    assert result.loc[:2, "atlas_score_source_for_ml"].tolist() == [
        "z_vs_decoys_consensus",
        "z_vs_compare_run_consensus",
        "z_vs_decoys_blend",
    ]
    assert result.loc[3:, "atlas_score"].isna().all()
    assert result.loc[3:, "atlas_score_source_for_ml"].tolist() == ["final_score", ""]


def test_spd_table_keeps_standardized_final_score_out_of_raw_consensus() -> None:
    frame = pd.DataFrame(
        {
            "drug_id": ["drug-a"],
            "target_id": ["target-a"],
            "final_score": [2.5],
            "final_score_source": ["z_vs_compare_run_consensus"],
            "spd_exposure_relevant": [1],
        }
    )

    staged = _spd_table(frame)

    assert "consensus_score" not in staged.columns
    result = add_standard_ml_columns(staged)
    assert result.loc[0, "atlas_score"] == pytest.approx(2.5)
    assert result.loc[0, "atlas_score_source_for_ml"] == "z_vs_compare_run_consensus"
