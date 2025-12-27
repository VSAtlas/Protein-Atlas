from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import dud_eval


def _write_csv(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_reranked_pretty_header_read_and_eval(tmp_path: Path) -> None:
    reranked_path = tmp_path / "consensus_reranked_scorch.csv"
    _write_csv(
        reranked_path,
        [
            "# run_id=RUN1",
            "# pdb_id=TEST, variant=HOLO, ph=pH7_0",
            "ligand,scorch_composite,SCORCH_score_used,SCORCH_certainty_used,final_rank,consensus_score,run_id,pdb_id,variant,ph_label",
            "actives_final_00001.pdbqt,0.9,1.0,0.9,1,0.5,RUN1,TEST,HOLO,pH7_0",
            "decoys_final_00001.pdbqt,0.1,0.2,0.5,2,0.2,RUN1,TEST,HOLO,pH7_0",
        ],
    )

    df = dud_eval.read_reranked_scorch_csv(reranked_path)
    assert "ligand" in df.columns
    assert len(df) == 2

    out_dir = tmp_path / "out"
    result = dud_eval.evaluate_target_post_docked_reranked_scorch(
        pdb_id="TEST",
        csv_path=reranked_path,
        out_dir=out_dir,
        lig_col_cli=None,
        score_col_cli=None,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        run_id="RUN1",
        variant_hint="HOLO",
        ph_tag="pH7_0",
        new_layout_active=True,
        layout_label="post_docked_new",
    )

    assert result is not None
    assert result.metrics is not None
    assert result.metrics["N"] == 2


def test_reranked_summary_files_written(tmp_path: Path) -> None:
    analysis_root = tmp_path / "analysis" / "dud_eval" / "demo_run"
    df = pd.DataFrame(
        [
            {
                "run_id": "demo_run",
                "pdb_id": "TEST",
                "variant": "HOLO",
                "pH": "pH7_0",
                "target_name": "Target",
                "library_name": "lib",
                "N": 2,
                "n_actives": 1,
                "actives_fraction": 0.5,
                "ROC_AUC": 1.0,
                "PR_AUC": 1.0,
                "logAUC": 0.2,
                "logAUC_adj": 0.1,
                "BEDROC_alpha_20": 0.5,
                "EF@1%": 1.0,
                "status_reason": "ok",
            }
        ]
    )

    dud_eval.emit_reranked_scorch_summary(df, analysis_root, run_label="demo_run", pretty_enabled=True)

    summary_path = analysis_root / "post_docked" / "consensus_reranked_scorch_summary_demo_run.tsv"
    pretty_path = analysis_root / "post_docked" / "consensus_reranked_scorch_summary_demo_run_pretty.txt"
    assert summary_path.exists()
    assert pretty_path.exists()
