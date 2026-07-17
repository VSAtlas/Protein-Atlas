from __future__ import annotations

import pandas as pd

from analysis.ml.score_scale_audit import audit_score_scale_frame


def test_scale_audit_fails_unsafe_atlas_score_derivation() -> None:
    frame = pd.DataFrame(
        {
            "atlas_score": [1.0, 2.0],
            "atlas_score_source_for_ml": [
                "z_vs_decoys_consensus",
                "final_score",
            ],
            "consensus_z_decoy_n": [20, 20],
            "consensus_z_decoy_sigma": [1.0, 1.0],
        }
    )

    result, sources, _ = audit_score_scale_frame(
        frame, feature_names=["atlas_score"]
    )

    assert result["status"] == "failed"
    assert result["publication_blocker"] is True
    assert result["n_unsafe_atlas_score"] == 1
    atlas_sources = sources.loc[
        sources["source_column"].eq("atlas_score_source_for_ml")
    ]
    assert set(atlas_sources["source_value"]) == {"z_vs_decoys_consensus", "final_score"}


def test_cross_run_final_score_source_remains_distinct_and_accepted() -> None:
    frame = pd.DataFrame(
        {
            "final_score": [2.5],
            "final_score_source": ["z_vs_compare_run_consensus"],
            "z_vs_compare_run_consensus": [2.5],
            "z_vs_compare_run_consensus_decoy_n": [30],
            "z_vs_compare_run_consensus_decoy_sigma": [1.2],
        }
    )

    result, sources, _ = audit_score_scale_frame(
        frame, feature_names=["final_score"]
    )

    assert result["n_unsafe_final_score"] == 0
    assert result["status"] == "passed"
    assert "z_vs_compare_run_consensus" in set(sources["source_value"])
