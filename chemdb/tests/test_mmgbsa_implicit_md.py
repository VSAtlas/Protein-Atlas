import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import post_docking.mmgbsa.mmgbsa_trajectory as mmgbsa_trajectory


def test_run_implicit_md_writes_inputs_and_calls_sander(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prmtop = tmp_path / "complex.prmtop"
    inpcrd = tmp_path / "complex.inpcrd"
    prmtop.write_text("X", encoding="utf-8")
    inpcrd.write_text("X", encoding="utf-8")

    calls = []

    def fake_select_sander_runner(logger, cfg):
        return [], "sander", "TEST"

    def fake_run(cmd, cwd=None, stdout=None, stderr=None):
        calls.append(list(cmd))
        cmd_list = list(cmd)
        if "-r" in cmd_list:
            out_path = Path(cwd) / cmd_list[cmd_list.index("-r") + 1]
            out_path.write_text("rst7", encoding="utf-8")
        if "-x" in cmd_list:
            out_path = Path(cwd) / cmd_list[cmd_list.index("-x") + 1]
            out_path.write_text("traj", encoding="utf-8")

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(
        mmgbsa_trajectory, "_select_sander_runner", fake_select_sander_runner
    )
    monkeypatch.setattr(mmgbsa_trajectory.subprocess, "run", fake_run)

    cfg = {
        "MMGBSA_MD_ENGINE": "sander",
        "MMGBSA_MD_DT_PS": 0.002,
        "MMGBSA_MD_HEAT_PS": 10,
        "MMGBSA_MD_EQUIL_PS": 20,
        "MMGBSA_MD_PROD_PS": 10,
        "MMGBSA_MD_FRAME_STRIDE_PS": 2,
        "MMGBSA_MD_RESTRAIN_PROTEIN_HEAVY": True,
        "MMGBSA_MD_IGB": 5,
        "MMGBSA_MD_SALTCON": 0.150,
        "MMGBSA_MD_TRAJ_NAME": "prod.nc",
        "MMGBSA_MD_TRAJ_FORMAT": "nc",
    }

    result = mmgbsa_trajectory.run_implicit_md(
        complex_prmtop=str(prmtop),
        complex_inpcrd=str(inpcrd),
        out_dir=str(tmp_path),
        cfg=cfg,
        replicate_index=1,
        seed=123,
        force=True,
        run=True,
    )

    md_dir = tmp_path / "rep1" / "md"
    assert (md_dir / "min.in").exists()
    assert (md_dir / "heat.in").exists()
    assert (md_dir / "equil.in").exists()
    assert (md_dir / "prod.in").exists()
    assert "ig=123" in (md_dir / "heat.in").read_text(encoding="utf-8")
    assert "ntwx=" in (md_dir / "prod.in").read_text(encoding="utf-8")
    assert len(calls) == 4

    traj_path = Path(result["traj_path"])
    assert traj_path.exists()
    assert traj_path.stat().st_size > 0
    assert result["n_frames_est"] > 1
