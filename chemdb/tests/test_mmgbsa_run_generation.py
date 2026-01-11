import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from post_docking.mmgbsa.run_mmgbsa import run_mmgbsa


def test_mmgbsa_run_generation_input_only(tmp_path: Path) -> None:
    work_dir = tmp_path / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    complex_prmtop = work_dir / "complex.prmtop"
    receptor_prmtop = work_dir / "receptor.prmtop"
    ligand_prmtop = work_dir / "ligand.prmtop"
    traj = work_dir / "mdcrd"

    for path in (complex_prmtop, receptor_prmtop, ligand_prmtop, traj):
        path.write_text("X", encoding="utf-8")

    cfg = {
        "MMGBSA_MMPBSA_ENABLED": True,
        "MMGBSA_MMPBSA_RUN": False,
        "MMGBSA_MMPBSA_STARTFRAME": 1,
        "MMGBSA_MMPBSA_ENDFRAME": 1,
        "MMGBSA_MMPBSA_INTERVAL": 1,
        "MMGBSA_MMPBSA_VERBOSE": 2,
        "MMGBSA_GB_IGB": 5,
        "MMGBSA_GB_SALTCON": 0.150,
        "MMGBSA_MMPBSA_INPUT_NAME": "mmpbsa.in",
        "MMGBSA_MMPBSA_LOG_NAME": "mmpbsa.log",
        "MMGBSA_MMPBSA_OUT_DAT": "FINAL_RESULTS_MMPBSA.dat",
        "MMGBSA_MMPBSA_OUT_CSV": "FINAL_RESULTS_MMPBSA.csv",
        "MMGBSA_DEFAULT_TRAJ_NAME": "mdcrd",
    }

    result = run_mmgbsa(
        complex_prmtop=str(complex_prmtop),
        receptor_prmtop=str(receptor_prmtop),
        ligand_prmtop=str(ligand_prmtop),
        trajectory_path=str(traj),
        work_dir=str(work_dir),
        cfg=cfg,
        force=False,
        run=False,
    )

    input_path = work_dir / "mmpbsa.in"
    assert input_path.exists()
    text = input_path.read_text(encoding="utf-8")
    assert "startframe=1" in text
    assert "endframe=1" in text
    assert "interval=1" in text
    assert "verbose=2" in text
    assert "igb=5" in text
    assert "saltcon=0.150" in text

    cmd_preview = result["cmd_preview"]
    assert "MMPBSA.py" in cmd_preview
    assert "-cp complex.prmtop" in cmd_preview
    assert "-rp receptor.prmtop" in cmd_preview
    assert "-lp ligand.prmtop" in cmd_preview
    assert "-y mdcrd" in cmd_preview
    assert "-o FINAL_RESULTS_MMPBSA.dat" in cmd_preview
    assert "-eo FINAL_RESULTS_MMPBSA.csv" in cmd_preview

    assert result["out_dat"] == str(work_dir / "FINAL_RESULTS_MMPBSA.dat")
    assert result["out_csv"] == str(work_dir / "FINAL_RESULTS_MMPBSA.csv")


# Manual integration example (not executed in CI):
# python run_mmgbsa.py \
#   --work-dir /path/to/mmgbsa_topologies/stage1/ligand \
#   --complex-prmtop /path/to/complex.prmtop \
#   --receptor-prmtop /path/to/receptor.prmtop \
#   --ligand-prmtop /path/to/ligand.prmtop \
#   --traj /path/to/mdcrd
