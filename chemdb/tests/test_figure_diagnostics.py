from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.ml import figure_diagnostics


def test_small_dataset_top_ks_are_clipped_then_deduped(
    tmp_path: Path, monkeypatch
) -> None:
    assert figure_diagnostics._normalize_top_ks([10, 20, 50, 100], 3) == [1, 3]
    monkeypatch.setattr(
        figure_diagnostics,
        "_plot_topk",
        lambda table, fig_path, *, title_prefix: None,
    )
    frame = pd.DataFrame(
        {
            "label": [1, 0, 1],
            "score": [0.9, 0.8, 0.1],
        }
    )

    result = figure_diagnostics.write_topk_outputs(
        frame,
        label_col="label",
        score_col="score",
        top_ks=[2, 3, 10, 20],
        out_dir=tmp_path,
        title_prefix="test",
    )

    table = result.tables["topk"]
    assert table["effective_k"].tolist() == [2, 3]
    assert table["top_k"].tolist() == [2, 3]
    assert table["effective_k"].is_unique


def test_failure_outputs_are_empty_when_threshold_has_no_errors(
    tmp_path: Path,
) -> None:
    frame = pd.DataFrame(
        {
            "drug_id": ["a", "b", "c", "d"],
            "label": [0, 0, 1, 1],
            "score": [0.1, 0.2, 0.8, 0.9],
        }
    )

    result = figure_diagnostics.write_failure_analysis_outputs(
        frame,
        label_col="label",
        score_col="score",
        out_dir=tmp_path,
        threshold=0.5,
        top_n=10,
        title_prefix="test",
    )

    assert result.tables["top_false_positives"].empty
    assert result.tables["top_false_negatives"].empty
    assert pd.read_csv(result.outputs["top_false_positives_csv"]).empty
    assert pd.read_csv(result.outputs["top_false_negatives_csv"]).empty


def test_failure_candidates_filter_by_threshold_error_before_ranking(
    tmp_path: Path,
) -> None:
    frame = pd.DataFrame(
        {
            "drug_id": ["fp", "tn-near", "tn", "fn-low", "fn-near", "tp"],
            "label": [0, 0, 0, 1, 1, 1],
            "score": [0.95, 0.75, 0.2, 0.05, 0.65, 0.9],
        }
    )

    result = figure_diagnostics.write_failure_analysis_outputs(
        frame,
        label_col="label",
        score_col="score",
        out_dir=tmp_path,
        threshold=0.8,
        top_n=10,
        title_prefix="test",
    )

    false_positives = result.tables["top_false_positives"]
    false_negatives = result.tables["top_false_negatives"]
    assert false_positives["drug_id"].tolist() == ["fp"]
    assert false_positives["error_type"].tolist() == ["false_positive"]
    assert false_negatives["drug_id"].tolist() == ["fn-low", "fn-near"]
    assert false_negatives["score"].tolist() == [0.05, 0.65]
    assert false_negatives["error_type"].tolist() == [
        "false_negative",
        "false_negative",
    ]
