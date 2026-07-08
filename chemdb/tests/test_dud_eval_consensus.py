from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis import dud_eval


def _write_csv(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_consensus_dedup_uses_max(tmp_path: Path) -> None:
    csv_path = tmp_path / "consensus_docking_scores.csv"
    _write_csv(
        csv_path,
        [
            "ligand,consensus_score,n_engines_with_data",
            "actives_final_00001__ms_a.pdbqt,0.3,3",
            "actives_final_00001__ms_b.pdbqt,0.9,3",
            "decoys_final_00002__ms_c.pdbqt,0.4,3",
        ],
    )

    out_dir = tmp_path / "out"
    result = dud_eval.evaluate_target_consensus(
        pdb_id="TEST",
        csv_path=csv_path,
        out_dir=out_dir,
        lig_col_cli=None,
        score_col_cli=None,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        run_id=None,
    )

    assert result is not None
    assert result.metrics is not None
    assert result.metrics["N"] == 2
    assert result.metrics["ROC_AUC"] == 1.0


def test_consensus_filters_no_data_rows(tmp_path: Path) -> None:
    csv_path = tmp_path / "consensus_docking_scores.csv"
    _write_csv(
        csv_path,
        [
            "ligand,consensus_score,n_engines_with_data",
            "actives_final_00001__ms_a.pdbqt,1.0,1",
            "decoys_final_00001__ms_a.pdbqt,0.2,1",
            "decoys_final_00002__ms_b.pdbqt,0.0,0",
        ],
    )

    out_dir = tmp_path / "out"
    result = dud_eval.evaluate_target_consensus(
        pdb_id="TEST",
        csv_path=csv_path,
        out_dir=out_dir,
        lig_col_cli=None,
        score_col_cli=None,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        run_id=None,
    )

    assert result is not None
    assert result.metrics is not None
    assert result.metrics["N"] == 2


def test_consensus_summary_files_written(tmp_path: Path) -> None:
    analysis_root = tmp_path / "analysis" / "dud_eval" / "demo_run"
    df = pd.DataFrame(
        [
            {
                "run_id": "demo_run",
                "pdb_id": "TEST",
                "variant": "",
                "pH": "",
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
            }
        ]
    )

    dud_eval.emit_consensus_summary(
        df, analysis_root, run_label="demo_run", pretty_enabled=True
    )

    summary_path = analysis_root / "consensus_summary_demo_run.tsv"
    pretty_path = analysis_root / "consensus_summary_demo_run_pretty.txt"
    assert summary_path.exists()
    assert pretty_path.exists()
