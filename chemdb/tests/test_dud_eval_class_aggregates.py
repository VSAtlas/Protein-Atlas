from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.dud_eval_core.orchestrate import _emit_class_aggregate_plots
from analysis.dud_eval_core.types import TargetEvaluation, TargetSpec


def _target_eval(
    *,
    n_total: int,
    n_actives: int,
    roc_auc: float,
    ef1: float,
    ef2: float,
    ef5: float,
    ef10: float,
    bedroc: float,
    pr_auc: float,
    logauc: float,
    logauc_adj: float,
) -> TargetEvaluation:
    frame = pd.DataFrame(
        {
            "is_active": [1] * n_actives + [0] * (n_total - n_actives),
            "best_score": list(range(n_total)),
        }
    )
    return TargetEvaluation(
        metrics=pd.Series(
            {
                "run_id": "R1",
                "pdb_id": "TEST",
                "N": float(n_total),
                "n_actives": float(n_actives),
                "actives_fraction": float(n_actives) / float(n_total),
                "ROC_AUC": roc_auc,
                "PR_AUC": pr_auc,
                "logAUC": logauc,
                "logAUC_adj": logauc_adj,
                "BEDROC_alpha_20": bedroc,
                "EF@1%": ef1,
                "EF@2%": ef2,
                "EF@5%": ef5,
                "EF@10%": ef10,
                "status_reason": "ok",
            }
        ),
        ligand_basenames=set(),
        has_run_id_column=False,
        run_ids=set(),
        status_reason="ok",
        z_selected=frame,
    )


@pytest.mark.parametrize(
    ("mode_slug", "score_high_is_better", "hist_xlabel", "title_suffix", "mode_label"),
    [
        ("docking", False, "Best docking score (lower = better)", "", "docking"),
        (
            "consensus",
            True,
            "Consensus score (higher = better; negated for plot)",
            " (consensus)",
            "consensus",
        ),
        (
            "post_docked_scorch",
            True,
            "Reranked SCORCH score (higher = better)",
            " (reranked_scorch)",
            "post_docked_scorch",
        ),
    ],
)
def test_class_aggregate_uses_macro_target_metrics_and_writes_boxplot(
    tmp_path: Path,
    mode_slug: str,
    score_high_is_better: bool,
    hist_xlabel: str,
    title_suffix: str,
    mode_label: str,
) -> None:
    target_specs_by_key = {
        "t1": TargetSpec(
            target_key="t1",
            pdb_id="1SQT",
            variant=None,
            ph_tag=None,
            csv_path=None,
            source="test",
        ),
        "t2": TargetSpec(
            target_key="t2",
            pdb_id="1W7X",
            variant=None,
            ph_tag=None,
            csv_path=None,
            source="test",
        ),
        "t3": TargetSpec(
            target_key="t3",
            pdb_id="1MV9",
            variant=None,
            ph_tag=None,
            csv_path=None,
            source="test",
        ),
    }
    eval_results = {
        "t1": _target_eval(
            n_total=20,
            n_actives=2,
            roc_auc=0.8,
            pr_auc=0.7,
            logauc=0.5,
            logauc_adj=0.4,
            bedroc=0.6,
            ef1=10.0,
            ef2=9.0,
            ef5=8.0,
            ef10=7.0,
        ),
        "t2": _target_eval(
            n_total=40,
            n_actives=4,
            roc_auc=0.6,
            pr_auc=0.5,
            logauc=0.3,
            logauc_adj=0.2,
            bedroc=0.4,
            ef1=5.0,
            ef2=4.0,
            ef5=3.0,
            ef10=2.0,
        ),
        "t3": _target_eval(
            n_total=30,
            n_actives=3,
            roc_auc=0.9,
            pr_auc=0.8,
            logauc=0.7,
            logauc_adj=0.6,
            bedroc=0.75,
            ef1=12.0,
            ef2=11.0,
            ef5=10.0,
            ef10=9.0,
        ),
    }

    class_results, generation_rows = _emit_class_aggregate_plots(
        analysis_root=tmp_path,
        run_id="R1",
        target_specs_by_key=target_specs_by_key,
        eval_results=eval_results,
        included_target_keys=set(target_specs_by_key),
        mode_slug=mode_slug,
        score_high_is_better=score_high_is_better,
        hist_xlabel=hist_xlabel,
        title_suffix=title_suffix,
        mode_label=mode_label,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        class_min_targets=1,
        write_placeholder_images=True,
    )

    protease = class_results["Protease"]
    assert protease.metrics is not None
    assert protease.metrics["EF@1%"] == pytest.approx(7.5)
    assert protease.metrics["EF@1%_median"] == pytest.approx(7.5)
    assert protease.metrics["EF@1%_iqr"] == pytest.approx(2.5)
    assert protease.metrics["EF@1%_std"] == pytest.approx(3.5355339059)
    assert protease.metrics["ROC_AUC"] == pytest.approx(0.7)
    assert protease.metrics["ROC_AUC_median"] == pytest.approx(0.7)
    assert protease.metrics["aggregation_method"] == "macro_per_target"
    assert protease.metrics["pooled_rows"] == 60
    assert protease.metrics["pooled_actives"] == 6

    report_row = next(row for row in generation_rows if row["class_name"] == "Protease")
    assert report_row["aggregation_method"] == "macro_per_target"
    assert report_row["EF@1%"] == pytest.approx(7.5)
    assert report_row["EF@1%_median"] == pytest.approx(7.5)
    assert report_row["EF@1%_iqr"] == pytest.approx(2.5)
    assert report_row["EF@1%_std"] == pytest.approx(3.5355339059)

    metrics_path = tmp_path / "class_aggregates" / mode_slug / "protease" / "metrics.tsv"
    plot_path = (
        tmp_path / "class_aggregates" / mode_slug / "protease" / "metric_boxplots.png"
    )
    placeholder_plot = (
        tmp_path
        / "class_aggregates"
        / mode_slug
        / "nuclear_receptor"
        / "metric_boxplots.png"
    )
    overall_ef_plot = (
        tmp_path
        / "class_aggregates"
        / mode_slug
        / "overall_by_metric"
        / "ef_at_1pct_by_class.png"
    )
    overall_roc_plot = (
        tmp_path
        / "class_aggregates"
        / mode_slug
        / "overall_by_metric"
        / "roc_auc_by_class.png"
    )
    assert metrics_path.exists()
    assert plot_path.exists()
    assert placeholder_plot.exists()
    assert overall_ef_plot.exists()
    assert overall_roc_plot.exists()
