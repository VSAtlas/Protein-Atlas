from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from record_data import write_scores_csv

RUN_ID = "test_run"


def _base_cfg(root: Path) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    cfg = {
        "OVERALL_DIR": str(root),
        "INPUT_DIR": str(root / "input_pdbs"),
        "OUTPUT_DIR": str(root / "processed_pdbs"),
        "DOCKED_DIR": str(root / "docked"),
        "OUTPUT_LIGANDS_DIR": str(root / "prepped_ligands"),
        "USE_GNINA": "false",
        "RUN_ID": RUN_ID,
    }
    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "OUTPUT_LIGANDS_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)
    return cfg


def _score_history() -> dict:
    return {
        "stage1": {
            "ligand1.pdbqt": {
                "score": -7.5,
                "valid": True,
                "reason": "",
                "heavy_atoms": 12,
                "le": None,
                "pains_flag": False,
            }
        }
    }


def test_write_scores_csv_with_prefix(tmp_path):
    cfg = _base_cfg(tmp_path / "prefixed")
    csv_path = write_scores_csv(cfg, "TEST1", _score_history(), csv_prefix="dud_")
    csv_file = Path(csv_path)

    assert csv_file.name == "dud_docking_score_summary.csv"
    assert csv_file.exists()

    long_path = csv_file.with_name("dud_docking_score_long.csv")
    assert long_path.exists()
    expected_parent = Path(cfg["DOCKED_DIR"]) / RUN_ID / "TEST1"
    assert csv_file.parent.resolve() == expected_parent.resolve()


def test_write_scores_csv_default_prefix(tmp_path):
    cfg = _base_cfg(tmp_path / "default")
    csv_path = write_scores_csv(cfg, "TEST2", _score_history())
    csv_file = Path(csv_path)

    assert csv_file.name == "docking_score_summary.csv"
    assert csv_file.exists()

    long_path = csv_file.with_name("docking_score_long.csv")
    assert long_path.exists()
    expected_parent = Path(cfg["DOCKED_DIR"]) / RUN_ID / "TEST2"
    assert csv_file.parent.resolve() == expected_parent.resolve()
