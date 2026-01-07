from __future__ import annotations

import csv
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

# Paths/constants reused from other slow integration tests
REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"
TEST_PDB_ID = "TEST"
TEST_PDB_PATH = REPO_ROOT / "input_pdbs" / f"{TEST_PDB_ID}.pdb"
TEST_LIBRARY_SUBDIR = "test_library_10"
PREPPED_TEST_LIB = REPO_ROOT / "prepped_ligands" / TEST_LIBRARY_SUBDIR
TEST_RUN_ID = "test_run"
DOCKED_ROOT = REPO_ROOT / "docked" / TEST_RUN_ID / TEST_PDB_ID
SUCCESS_SENTINEL = "No proteins recorded as failed."

# Keep integration runs lightweight by disabling MMGBSA in this test.
MMGBSA_DISABLE_ENV = {
    "MMGBSA_ENABLED": "false",
    "MMGBSA_MD_ENABLED": "false",
    "MMGBSA_MMPBSA_ENABLED": "false",
    "MMGBSA_MMPBSA_RUN": "false",
    "MD_FIVE_REPLICATE": "false",
}

pytestmark = pytest.mark.slow


def _have_test_inputs() -> bool:
    if not MAIN_PY.exists():
        return False
    if not TEST_PDB_PATH.exists():
        return False
    if not PREPPED_TEST_LIB.exists():
        return False
    return any(PREPPED_TEST_LIB.glob("*.pdbqt"))


def _extract_log_path(output: str) -> Path | None:
    for line in output.splitlines():
        if "Log ->" in line:
            _, _, tail = line.partition("Log ->")
            path_str = tail.strip()
            if path_str:
                return Path(path_str)
    return None


def _find_recent_long_csv(start_time: float) -> Path | None:
    """
    Return the newest docking_score_long.csv under docked/<RUN_ID>/TEST/... whose
    mtime is >= start_time. We accept any variant/pH subdirectory.
    """
    if not DOCKED_ROOT.exists():
        return None

    candidates: list[tuple[float, Path]] = []
    for csv_path in DOCKED_ROOT.rglob("docking_score_long.csv"):
        try:
            stat = csv_path.stat()
        except FileNotFoundError:
            continue
        if stat.st_mtime < start_time or stat.st_size <= 0:
            continue
        candidates.append((stat.st_mtime, csv_path))

    if not candidates:
        return None
    candidates.sort(key=lambda pair: pair[0], reverse=True)
    return candidates[0][1]


def test_decoy_t_scores_present_in_fast_test_mode() -> None:
    """
    Run the lightweight FDA+DUD flow via:

        python main.py -test -test-fda -fast

    and assert that the resulting docking_score_long.csv includes a
    non-zero t_vs_decoys column.
    """
    if not _have_test_inputs():
        pytest.skip("Skipping decoy T-score test: TEST inputs or test library missing.")

    assert MAIN_PY.exists(), f"main.py not found at {MAIN_PY}"

    env = os.environ.copy()
    # Ensure the code-path that mixes FDA + DUD is active.
    env["TEST_MODE_ENABLE"] = "fda+dud"
    env["ATLAS_RUN_ID"] = TEST_RUN_ID
    env["PYTHONUNBUFFERED"] = "1"
    env.update(MMGBSA_DISABLE_ENV)

    cmd = [sys.executable, str(MAIN_PY), "-test", "-test-fda", "-fast", "--run-id", TEST_RUN_ID]
    start_time = time.time()
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        text=True,
        capture_output=True,
    )

    combined_output = (proc.stdout or "") + "\n" + (proc.stderr or "")
    log_path = _extract_log_path(combined_output)
    print("=== decoy T-score test ===")
    if log_path is not None:
        try:
            rel_log = log_path.relative_to(REPO_ROOT)
        except ValueError:
            rel_log = log_path
        print(f"LOG: {rel_log}")
    print("Last 20 lines of output:")
    lines = combined_output.splitlines()
    for line in lines[-20:]:
        print(line)

    assert proc.returncode == 0, (
        f"main.py exited with {proc.returncode}\n"
        f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    assert SUCCESS_SENTINEL in combined_output, "Expected success sentinel in output."

    long_csv = _find_recent_long_csv(start_time)
    assert long_csv is not None, "Unable to locate docking_score_long.csv from this run."

    with long_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames and "t_vs_decoys" in reader.fieldnames, (
            f"t_vs_decoys column missing in {long_csv}"
        )
        t_values = []
        for row in reader:
            value_str = (row.get("t_vs_decoys") or "").strip()
            if not value_str:
                continue
            try:
                t_values.append(float(value_str))
            except ValueError:
                continue

    assert t_values, f"No T-score values recorded in {long_csv}"
    assert any(abs(val) > 1e-9 for val in t_values), (
        f"All T-scores were zero in {long_csv}; expected decoy comparison to produce non-zero."
    )
