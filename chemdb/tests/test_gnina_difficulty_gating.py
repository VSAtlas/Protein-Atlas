from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"
DOCKED_ROOT = REPO_ROOT / "docked"
CONFIGS_ROOT = REPO_ROOT / "configs"

pytestmark = pytest.mark.slow


def _extract_log_path(output: str) -> Path | None:
    for line in output.splitlines():
        if "Log ->" in line:
            _, _, tail = line.partition("Log ->")
            path_str = tail.strip()
            if path_str:
                return Path(path_str)
    return None


def _find_recent_path(root: Path, pattern: str, start_ts: float) -> Path | None:
    candidates: list[tuple[float, Path]] = []
    for p in root.rglob(pattern):
        try:
            st = p.stat()
        except FileNotFoundError:
            continue
        if st.st_mtime >= start_ts:
            candidates.append((st.st_mtime, p))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


@pytest.mark.parametrize(
    "forced_bucket,expect_gnina",
    [
        ("hard", True),
        ("easy", False),
    ],
)
def test_gnina_runs_only_for_hard_or_degenerate(forced_bucket: str, expect_gnina: bool) -> None:
    if not MAIN_PY.exists():
        pytest.skip("main.py missing")
    if shutil.which("gnina") is None:
        pytest.skip("gnina executable not available on PATH")

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    env["ENABLE_GNINA"] = "1"
    env["USE_GNINA"] = "1"
    env["FORCE_TARGET_DIFFICULTY"] = forced_bucket

    cmd = [
        sys.executable,
        str(MAIN_PY),
        "-test",
        "-fast",
        "--single",
        "dexamethasone",
    ]

    start_ts = time.time()
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        text=True,
        capture_output=True,
    )
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    log_path = _extract_log_path(combined)

    assert proc.returncode == 0, combined

    vina_out = _find_recent_path(DOCKED_ROOT, "stage1/*.pdbqt", start_ts)
    assert vina_out is not None, "Expected Vina output from stage1"

    gnina_manifest = _find_recent_path(CONFIGS_ROOT, "gnina_stage1/gnina.json", start_ts)
    # GNINA may skip overwriting existing outputs in resume-like paths; do not rely on mtime.
    gnina_out_any = next(DOCKED_ROOT.rglob("gnina_stage1/*.pdbqt"), None)

    if expect_gnina:
        assert gnina_manifest is not None, "Expected gnina.json manifest for hard target"
        assert gnina_out_any is not None or "[gnina.call]" in combined, (
            "Expected GNINA to be invoked for hard target"
        )
    else:
        assert gnina_manifest is None, "Did not expect gnina.json for easy target"
        assert "[gnina.call]" not in combined, "Did not expect GNINA to be invoked for easy target"

    if expect_gnina and log_path is not None and log_path.exists():
        content = log_path.read_text(encoding="utf-8", errors="replace")
        vina_idx = content.find("[vina.call]")
        gnina_idx = content.find("[gnina.stage]")
        assert vina_idx != -1 and gnina_idx != -1, "Expected vina.call and gnina.stage in logs"
        assert vina_idx < gnina_idx, "Expected Vina to run before GNINA"
