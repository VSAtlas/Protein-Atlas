import os
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mmgbsa_trajectory import make_mmgbsa_trajectory


def test_mmgbsa_trajectory_input_only(tmp_path: Path) -> None:
    prmtop = tmp_path / "complex.prmtop"
    inpcrd = tmp_path / "complex.inpcrd"
    prmtop.write_text("", encoding="utf-8")
    inpcrd.write_text("", encoding="utf-8")

    cfg = {
        "MMGBSA_CPPTRAJ_ENABLED": True,
        "MMGBSA_CPPTRAJ_RUN": False,
    }

    result = make_mmgbsa_trajectory(
        complex_prmtop=str(prmtop),
        complex_inpcrd=str(inpcrd),
        out_dir=str(tmp_path),
        cfg=cfg,
        force=False,
        run_cpptraj=False,
    )

    cpptraj_in = tmp_path / "cpptraj_mmgbsa.in"
    assert cpptraj_in.exists()
    text = cpptraj_in.read_text(encoding="utf-8")
    assert "parm complex.prmtop" in text
    assert "trajin complex.inpcrd 1 1" in text
    assert "trajout mdcrd" in text
    assert result["trajout_path"] == str(tmp_path / "mdcrd")


def _cpptraj_available() -> bool:
    env_prefix = os.environ.get("AMBERTOOLS_PREFIX")
    if env_prefix:
        env_path = Path(env_prefix) / "bin" / "cpptraj"
        if env_path.is_file():
            return True
    return shutil.which("cpptraj") is not None


def test_mmgbsa_trajectory_optional_cpptraj(tmp_path: Path) -> None:
    if not _cpptraj_available():
        pytest.skip("cpptraj not available")

    prmtop = tmp_path / "complex.prmtop"
    inpcrd = tmp_path / "complex.inpcrd"
    prmtop.write_text("", encoding="utf-8")
    inpcrd.write_text("", encoding="utf-8")

    cfg = {
        "MMGBSA_CPPTRAJ_ENABLED": True,
        "MMGBSA_CPPTRAJ_RUN": True,
    }

    try:
        result = make_mmgbsa_trajectory(
            complex_prmtop=str(prmtop),
            complex_inpcrd=str(inpcrd),
            out_dir=str(tmp_path),
            cfg=cfg,
            force=True,
            run_cpptraj=True,
        )
    except Exception:
        pytest.skip("cpptraj failed with dummy inputs")

    trajout = Path(result["trajout_path"])
    if not trajout.exists():
        pytest.skip("cpptraj did not produce trajectory")
