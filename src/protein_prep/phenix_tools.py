"""Phenix and MolProbity command helpers used by receptor cleaning."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union

from protein_prep.os_utils import _powershell, _run_and_log, _win_path
from protein_prep.prep_utils import _cfg


PHENIX_DIR = _cfg("PHENIX_DIR", "", "phenix_dir")
PHENIX_LIB_PATH = _cfg("PHENIX_LIB_PATH", "", "phenix_lib_path")
PHENIX_CLEAN_SCRIPT = _cfg("PHENIX_CLEAN_SCRIPT", "", "phenix_clean_script")
PDBTOOLS_BAT = str(Path(PHENIX_DIR) / "phenix.pdbtools.bat") if PHENIX_DIR else ""
MOLPROBITY_BAT = str(Path(PHENIX_DIR) / "phenix.molprobity.bat") if PHENIX_DIR else ""
PHENIX_PYTHON_BAT = str(Path(PHENIX_DIR) / "phenix.python.bat") if PHENIX_DIR else ""


def _phenix_detect() -> tuple[str, list[str] | None, dict]:
    """
    Return (mode, cmd_list, env) for phenix.pdbtools detection.
    mode ∈ {'linux-path','linux-explicit','wsl-bat','unavailable'}
    """
    env = dict(os.environ)
    if PHENIX_LIB_PATH:
        env["PYTHONPATH"] = PHENIX_LIB_PATH

    # 1) On PATH (native Linux installs)
    if shutil.which("phenix.pdbtools"):
        return ("linux-path", ["phenix.pdbtools"], env)

    # 2) Explicit Linux path via PHENIX_DIR (your case)
    if PHENIX_DIR:
        explicit = Path(PHENIX_DIR) / "phenix.pdbtools"
        if explicit.exists() and os.access(str(explicit), os.X_OK):
            return ("linux-explicit", [str(explicit)], env)

    # 3) Windows .bat (usable from WSL via PowerShell)
    if PDBTOOLS_BAT and Path(PDBTOOLS_BAT).exists():
        return ("wsl-bat", [PDBTOOLS_BAT], env)

    return ("unavailable", None, env)


def run_phenix_pdbtools(
    input_pdb: Union[str, Path],
    output_pdb: Union[str, Path],
    remove_waters: bool = True,
) -> bool:
    """
    Run Phenix pdbtools if available (Linux, or Windows via PowerShell under WSL).
    Return True if successful. Always log branch as: [phenix] mode=...
    """
    inp = os.path.abspath(str(input_pdb)).replace("\\", "/")
    outp = os.path.abspath(str(output_pdb)).replace("\\", "/")

    mode, base_cmd, env = _phenix_detect()
    logging.info("[phenix] mode=%s", mode)

    if mode == "unavailable" or not base_cmd:
        logging.warning("Phenix not available; skipping pdbtools polish")
        return False

    cmd = base_cmd + [inp, f"output.file_name={outp}"]
    if remove_waters:
        cmd.append('remove="resname HOH"')

    try:
        _run_and_log(cmd, env=env, check=True)
        if os.path.exists(outp) and os.path.getsize(outp) > 0:
            logging.info("phenix.pdbtools wrote %s", outp)
            return True
        logging.warning(
            "phenix.pdbtools finished but did not create expected output: %s", outp
        )
        return False
    except Exception as e:
        logging.error("phenix.pdbtools failed: %s", e)
        # Diagnostic: phenix.python import iotbx (same env)
        py_cmd = None
        if PHENIX_DIR and (Path(PHENIX_DIR) / "phenix.python").exists():
            py_cmd = [str(Path(PHENIX_DIR) / "phenix.python")]
        elif shutil.which("phenix.python"):
            py_cmd = ["phenix.python"]
        elif PHENIX_PYTHON_BAT and Path(PHENIX_PYTHON_BAT).exists():
            py_cmd = [PHENIX_PYTHON_BAT]
        if py_cmd:
            try:
                _run_and_log(
                    py_cmd + ["-c", "from iotbx import pdb; print('iotbx_ok')"],
                    env=env,
                    check=False,
                )
            except Exception:
                pass
        return False


def run_windows_phenix_clean_script(
    loop_fixed_pdb: Union[str, Path], nolig_dir: Union[str, Path]
) -> int:
    """Call your PHENIX_CLEAN_SCRIPT using Windows Phenix Python from WSL. Returns returncode."""
    if (
        not PHENIX_PYTHON_BAT
        or not Path(PHENIX_PYTHON_BAT).exists()
        or not PHENIX_CLEAN_SCRIPT
    ):
        return 127
    win_script = _win_path(PHENIX_CLEAN_SCRIPT)
    win_in = _win_path(loop_fixed_pdb)
    win_outdir = _win_path(nolig_dir)
    ps = (
        "$ErrorActionPreference='Stop';"
        f"if (!(Test-Path '{win_outdir}')) {{ New-Item -ItemType Directory -Force -Path '{win_outdir}' | Out-Null }};"
        f"Push-Location '{win_outdir}';"
        f"& '{_win_path(PHENIX_PYTHON_BAT)}' '{win_script}' '{win_in}' '{win_outdir}';"
        "Pop-Location;"
    )
    r = _powershell(ps)
    return r.returncode


def run_molprobity_validate(
    pdb_path: Union[str, Path], work_dir: Optional[Union[str, Path]] = None
) -> int:
    """
    Try running MolProbity validation via 'phenix.molprobity' if available.
    Returns the process return code (0 means success). Non-fatal if missing.
    """
    exe = shutil.which("phenix.molprobity")
    if not exe:
        logging.info("MolProbity not available; skipping validation.")
        return 127
    work = Path(work_dir or Path(pdb_path).with_suffix(".molprobity")).resolve()
    work.mkdir(parents=True, exist_ok=True)
    cmd = [exe, str(Path(pdb_path).resolve())]
    logging.info("Running MolProbity: %s", " ".join(cmd))
    r = subprocess.run(
        cmd, cwd=str(work), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if r.returncode != 0:
        logging.warning(
            "MolProbity failed (rc=%s)\nSTDERR:\n%s", r.returncode, (r.stderr or "")
        )
    else:
        logging.info("MolProbity finished; results under %s", work)
    return r.returncode
