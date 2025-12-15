from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path

import pytest

from chemdb.target_difficulty import evaluate_difficulty_for_target


@pytest.mark.parametrize(
    "roc_auc,expected",
    [
        (0.85, "easy"),
        (0.70, "medium"),
        (0.50, "hard"),
        (float("nan"), "degenerate"),
    ],
)
def test_evaluate_difficulty_for_target_maps_buckets(tmp_path: Path, monkeypatch, caplog, roc_auc: float, expected: str) -> None:
    # Stub dud_eval.evaluate_target so we don't depend on full DUD-E inputs.
    def _fake_evaluate_target(**_kwargs):
        return SimpleNamespace(metrics={"ROC_AUC": roc_auc, "N": 10, "n_actives": 2})

    import dud_eval

    monkeypatch.setattr(dud_eval, "evaluate_target", _fake_evaluate_target, raising=True)

    csv_path = tmp_path / "docking_score_long.csv"
    csv_path.write_text("ligand,score,label\nL1,-8.0,active\nL2,-7.0,decoy\n")
    analysis_root = tmp_path / "analysis"

    with caplog.at_level("INFO"):
        td = evaluate_difficulty_for_target(
            pdb_id="TEST",
            csv_path=csv_path,
            analysis_root=analysis_root,
        )

    assert td.difficulty == expected
    assert any("[difficulty]" in r.message for r in caplog.records)

