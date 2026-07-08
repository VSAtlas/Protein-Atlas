from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from analysis.dud_eval_engine import (
    emit_consensus_summary,
    emit_reranked_scorch_summary,
    evaluate_target,
    evaluate_target_consensus,
    evaluate_target_post_docked_reranked_scorch,
)


def _write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _assert_metric_snapshot(metrics: pd.Series) -> None:
    assert metrics["N"] == 4
    assert metrics["n_actives"] == 2
    assert metrics["ROC_AUC"] == pytest.approx(1.0)
    assert metrics["PR_AUC"] == pytest.approx(1.0)
    assert metrics["logAUC"] == pytest.approx(1.0)
    assert metrics["logAUC_adj"] == pytest.approx(0.85538)
    assert metrics["BEDROC_alpha_20"] == pytest.approx(1.0000454019910097)
    assert metrics["EF@1%"] == pytest.approx(2.0)
    assert metrics["EF@2%"] == pytest.approx(2.0)
    assert metrics["EF@5%"] == pytest.approx(2.0)
    assert metrics["EF@10%"] == pytest.approx(2.0)


def test_modular_engine_smoke_with_snapshot_outputs(tmp_path: Path) -> None:
    long_csv = tmp_path / "docking_score_long.csv"
    _write_lines(
        long_csv,
        [
            "ligand_file,score",
            "actives_final_00001__ms_a.pdbqt,-9.0",
            "actives_final_00001__ms_b.pdbqt,-8.0",
            "actives_final_00002_pose1.pdbqt,-7.0",
            "decoys_final_00001__ms_c.pdbqt,-4.0",
            "decoys_final_00002_pose2.pdbqt,-3.0",
        ],
    )

    consensus_csv = tmp_path / "consensus_docking_scores.csv"
    _write_lines(
        consensus_csv,
        [
            "ligand,consensus_score,n_engines_with_data",
            "actives_final_00001__ms_a.pdbqt,0.9,3",
            "actives_final_00001__ms_b.pdbqt,0.7,3",
            "actives_final_00002_pose1.pdbqt,0.6,2",
            "decoys_final_00001__ms_c.pdbqt,0.4,3",
            "decoys_final_00002_pose2.pdbqt,0.1,1",
        ],
    )

    reranked_csv = tmp_path / "consensus_reranked_scorch.csv"
    _write_lines(
        reranked_csv,
        [
            "ligand,final_score,rescored_flag,run_id,pdb_id,variant,ph_label",
            "actives_final_00001__ms_a.pdbqt,0.95,1,R1,TEST,HOLO,pH7_0",
            "actives_final_00001__ms_b.pdbqt,0.85,1,R1,TEST,HOLO,pH7_0",
            "actives_final_00002_pose1.pdbqt,0.75,1,R1,TEST,HOLO,pH7_0",
            "decoys_final_00001__ms_c.pdbqt,0.25,1,R1,TEST,HOLO,pH7_0",
            "decoys_final_00002_pose2.pdbqt,0.15,1,R1,TEST,HOLO,pH7_0",
        ],
    )

    docking_out = tmp_path / "out" / "docking"
    consensus_out = tmp_path / "out" / "consensus"
    reranked_out = tmp_path / "out" / "reranked"

    docking_eval = evaluate_target(
        pdb_id="TEST",
        csv_path=long_csv,
        out_dir=docking_out,
        lig_col_cli=None,
        score_col_cli=None,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        run_id=None,
    )
    consensus_eval = evaluate_target_consensus(
        pdb_id="TEST",
        csv_path=consensus_csv,
        out_dir=consensus_out,
        lig_col_cli=None,
        score_col_cli=None,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        run_id=None,
    )
    reranked_eval = evaluate_target_post_docked_reranked_scorch(
        pdb_id="TEST",
        csv_path=reranked_csv,
        out_dir=reranked_out,
        lig_col_cli=None,
        score_col_cli=None,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        run_id="R1",
        variant_hint="HOLO",
        ph_tag="pH7_0",
        new_layout_active=True,
        layout_label="post_docked_new",
    )

    assert docking_eval is not None and docking_eval.metrics is not None
    assert consensus_eval is not None and consensus_eval.metrics is not None
    assert reranked_eval is not None and reranked_eval.metrics is not None

    _assert_metric_snapshot(docking_eval.metrics)
    _assert_metric_snapshot(consensus_eval.metrics)
    _assert_metric_snapshot(reranked_eval.metrics)

    for out_dir in (docking_out, consensus_out, reranked_out):
        metrics_path = out_dir / "metrics.tsv"
        assert metrics_path.exists()
        metrics_df = pd.read_csv(metrics_path, sep="\t")
        assert len(metrics_df) == 1
        assert {"pdb_id", "N", "ROC_AUC", "PR_AUC"}.issubset(metrics_df.columns)

    analysis_root = tmp_path / "analysis" / "dud_eval" / "smoke"
    consensus_df = pd.DataFrame(
        [
            {
                **consensus_eval.metrics.to_dict(),
                "pdb_id": "TEST",
                "variant": "",
                "pH": "",
                "target_name": "Target",
                "library_name": "lib",
            }
        ]
    )
    reranked_df = pd.DataFrame(
        [
            {
                **reranked_eval.metrics.to_dict(),
                "pdb_id": "TEST",
                "variant": "",
                "pH": "",
                "target_name": "Target",
                "library_name": "lib",
            }
        ]
    )

    emit_consensus_summary(consensus_df, analysis_root, run_label="smoke", pretty_enabled=True)
    emit_reranked_scorch_summary(
        reranked_df,
        analysis_root,
        run_label="smoke",
        pretty_enabled=True,
    )

    consensus_summary = analysis_root / "consensus_summary_smoke.tsv"
    reranked_summary = (
        analysis_root / "post_docked" / "consensus_reranked_scorch_summary_smoke.tsv"
    )
    assert consensus_summary.exists()
    assert reranked_summary.exists()
    assert (analysis_root / "consensus_summary_smoke_pretty.txt").exists()
    assert (
        analysis_root
        / "post_docked"
        / "consensus_reranked_scorch_summary_smoke_pretty.txt"
    ).exists()

    consensus_table = pd.read_csv(consensus_summary, sep="\t")
    reranked_table = pd.read_csv(reranked_summary, sep="\t")
    assert len(consensus_table) == 1
    assert len(reranked_table) == 1
    assert list(consensus_table.columns[:3]) == ["variant", "pH", "run_id"]
    assert list(reranked_table.columns[:3]) == ["variant", "pH", "run_id"]
