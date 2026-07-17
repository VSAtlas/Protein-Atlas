from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from analysis.ml.materialize_final_score import (
    COMPARISON_RUN_CONSENSUS_Z_COLUMN,
    materialize_final_score,
)


def _materialize(
    tmp_path: Path, frame: pd.DataFrame, *, fallback: str | None = None
) -> tuple[pd.DataFrame, dict[str, object]]:
    source = tmp_path / f"input-{fallback or 'default'}.csv"
    output = tmp_path / f"output-{fallback or 'default'}.csv"
    audit = tmp_path / f"audit-{fallback or 'default'}"
    frame.to_csv(source, index=False)
    kwargs = {"fallback_column": fallback} if fallback else {}
    manifest = materialize_final_score(source, output, audit, **kwargs)
    return pd.read_csv(output), manifest


def test_materializer_uses_only_the_requested_fallback_scale(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "pdb_id": ["1ABC"],
            "z_vs_decoys_consensus": [1.25],
            "z_vs_compare_run_consensus": [3.5],
        }
    )

    same_run, same_manifest = _materialize(tmp_path, frame)
    assert same_run.loc[0, "final_score"] == pytest.approx(1.25)
    assert same_run.loc[0, "final_score_source"] == "z_vs_decoys_consensus"
    assert same_run.loc[0, "z_vs_compare_run_consensus"] == pytest.approx(3.5)
    assert same_manifest["fallback_input_column"] == "z_vs_decoys_consensus"

    cross_run, cross_manifest = _materialize(
        tmp_path,
        frame,
        fallback=COMPARISON_RUN_CONSENSUS_Z_COLUMN,
    )
    assert cross_run.loc[0, "final_score"] == pytest.approx(3.5)
    assert (
        cross_run.loc[0, "final_score_source"]
        == COMPARISON_RUN_CONSENSUS_Z_COLUMN
    )
    assert cross_run.loc[0, "z_vs_decoys_consensus"] == pytest.approx(1.25)
    assert cross_manifest["fallback_input_column"] == COMPARISON_RUN_CONSENSUS_Z_COLUMN
    assert cross_manifest["policy"]["same_run_and_comparison_run_consensus_are_distinct"]


def test_materializer_requires_the_explicit_fallback_column(tmp_path: Path) -> None:
    source = tmp_path / "input.csv"
    pd.DataFrame({"z_vs_decoys_consensus": [1.0]}).to_csv(source, index=False)

    with pytest.raises(ValueError, match="missing requested fallback column"):
        materialize_final_score(
            source,
            tmp_path / "output.csv",
            tmp_path / "audit",
            fallback_column=COMPARISON_RUN_CONSENSUS_Z_COLUMN,
        )
