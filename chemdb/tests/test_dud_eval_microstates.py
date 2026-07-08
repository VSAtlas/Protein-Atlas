from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis import dud_eval


def test_parse_name_and_label_collapses_microstates() -> None:
    lig_id_decoy, is_active_decoy = dud_eval.parse_name_and_label(
        "decoys_final_00044__ms_e5d9fb2ba4a3a8c4.pdbqt"
    )
    assert lig_id_decoy == "decoys_final_00044"
    assert is_active_decoy == 0

    lig_id_active, is_active_active = dud_eval.parse_name_and_label(
        "actives_final_00012__ms_0123abcd4567ef89.pdbqt"
    )
    assert lig_id_active == "actives_final_00012"
    assert is_active_active == 1

    lig_id_plain, _ = dud_eval.parse_name_and_label("decoys_final_00044.pdbqt")
    assert lig_id_plain == "decoys_final_00044"


def test_evaluate_target_collapses_microstates(tmp_path: Path) -> None:
    csv_path = tmp_path / "docking_score_long.csv"
    csv_lines = [
        "ligand_file,score",
        "actives_final_00001__ms_aaaaaaaaaaaaaaaa.pdbqt,-6.0",
        "actives_final_00001__ms_bbbbbbbbbbbbbbbb.pdbqt,-8.0",
        "decoys_final_00002__ms_cccccccccccccccc.pdbqt,-3.0",
        "decoys_final_00002__ms_dddddddddddddddd.pdbqt,-5.0",
        "decoys_final_00003__ms_eeeeeeeeeeeeeeee.pdbqt,-4.0",
    ]
    csv_path.write_text("\n".join(csv_lines) + "\n", encoding="utf-8")

    out_dir = tmp_path / "out"
    result = dud_eval.evaluate_target(
        pdb_id="TEST",
        csv_path=csv_path,
        out_dir=out_dir,
        lig_col_cli=None,
        score_col_cli=None,
        bedroc_alpha=20.0,
        logauc_lambda=1e-3,
        run_id=None,
        valid_only=False,
    )

    assert result is not None
    assert result.metrics is not None
    assert result.metrics["N"] == 3  # base ligands: 1 active + 2 decoys
    assert result.metrics["n_actives"] == 1
