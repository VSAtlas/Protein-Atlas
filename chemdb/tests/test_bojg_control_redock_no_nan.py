from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.slow
def test_bojg_control_redock_scores_are_finite() -> None:
    root = Path(__file__).resolve().parents[2]
    run_id = "bojg_ctrl_no_nan"
    cmd = [
        sys.executable,
        str(root / "main.py"),
        "-bojg",
        "-fast",
        "--no-docking",
        "--run-id",
        run_id,
    ]

    res = subprocess.run(
        cmd,
        cwd=root,
        env=os.environ.copy(),
        text=True,
        capture_output=True,
    )
    combined = (res.stdout or "") + "\n" + (res.stderr or "")
    assert res.returncode == 0, (
        "main.py BOJG no-docking run failed\n"
        f"stdout:\n{res.stdout}\n\nstderr:\n{res.stderr}"
    )

    control_lines = [ln for ln in combined.splitlines() if "[control-redock]" in ln]
    assert control_lines, "No [control-redock] lines were found in run output."

    nan_lines = [ln for ln in control_lines if re.search(r"\bscore=nan\b", ln)]
    assert not nan_lines, "Control redock produced NaN score lines:\n" + "\n".join(
        nan_lines[:20]
    )
