from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

# Root of the repo (chemdb/tests/ -> chemdb -> repo root)
REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"

TEST_PDB_ID = "TEST"
TEST_PDB_PATH = REPO_ROOT / "input_pdbs" / f"{TEST_PDB_ID}.pdb"
TEST_RUN_ID = "test_run"

# Small prepped library used when TEST_MODE_ENABLE=dud
# (config.txt should map TEST/dud -> prepped_ligands/test_library_10)
TEST_LIBRARY_SUBDIR = "test_library_10"
PREPPED_TEST_LIB = REPO_ROOT / "prepped_ligands" / TEST_LIBRARY_SUBDIR

# Root under which we expect docked outputs for TEST to land:
# - no variant: docked/<RUN_ID>/TEST/...
# - variant:    docked/<RUN_ID>/TEST/APO/... or docked/<RUN_ID>/TEST/HOLO/...
# - variant+pH: docked/<RUN_ID>/TEST/APO/pH*/... or docked/<RUN_ID>/TEST/HOLO/pH*/...
DOCKED_ROOT = REPO_ROOT / "docked" / TEST_RUN_ID / TEST_PDB_ID

SUCCESS_SENTINEL = "No proteins recorded as failed."

# Keep integration runs lightweight by disabling MMGBSA in these tests.
MMGBSA_DISABLE_ENV = {
    "MMGBSA_ENABLED": "false",
    "MMGBSA_MD_ENABLED": "false",
    "MMGBSA_MMPBSA_ENABLED": "false",
    "MMGBSA_MMPBSA_RUN": "false",
    "MD_FIVE_REPLICATE": "false",
}

# Mark the whole module as slow so it only runs when explicitly requested
pytestmark = pytest.mark.slow

SCENARIOS = [
    (
        "legacy_ph_on",
        {
            "APO_HOLO_MODE": "legacy",
            "PH_ENSEMBLE": "1",
        },
    ),
    (
        "apo_ph_on",
        {
            "APO_HOLO_MODE": "APO",
            "PH_ENSEMBLE": "1",
        },
    ),
    (
        "holo_ph_on",
        {
            "APO_HOLO_MODE": "HOLO",
            "PH_ENSEMBLE": "1",
        },
    ),
]


def _have_test_inputs() -> bool:
    """Return True if TEST.pdb and the small prepped library exist."""
    if not MAIN_PY.exists():
        return False
    if not TEST_PDB_PATH.exists():
        return False
    if not PREPPED_TEST_LIB.exists():
        return False
    # Require at least one ligand file so we know the library is populated
    has_ligand = any(PREPPED_TEST_LIB.glob("*.pdbqt"))
    return has_ligand


def _extract_log_path(output: str) -> Path | None:
    """Parse the 'Log -> ...' line from main.py output if present."""
    for line in output.splitlines():
        if "Log ->" in line:
            _, _, tail = line.partition("Log ->")
            path_str = tail.strip()
            if path_str:
                return Path(path_str)
    return None


def _find_recent_summary(
    start_time: float, expected_variant: str | None
) -> Path | None:
    """
    Find the newest docking_score_summary.csv under docked/<RUN_ID>/TEST that was
    written or updated after start_time.

    If expected_variant is provided ('APO' or 'HOLO'), we also require the
    summary path to live under docked/<RUN_ID>/TEST/<expected_variant>/...
    """
    if not DOCKED_ROOT.exists():
        return None

    candidates: list[tuple[float, Path]] = []
    for summary in DOCKED_ROOT.rglob("docking_score_summary.csv"):
        try:
            stat = summary.stat()
        except FileNotFoundError:
            continue
        if stat.st_mtime < start_time or stat.st_size <= 0:
            continue

        if expected_variant is not None:
            try:
                rel = summary.relative_to(DOCKED_ROOT)
            except ValueError:
                # Not under DOCKED_ROOT for some reason; skip
                continue
            # Require first path component to match the variant
            # e.g. APO/pH7_6/stage3/docking_score_summary.csv
            if not rel.parts:
                continue
            if rel.parts[0] != expected_variant:
                continue

        candidates.append((stat.st_mtime, summary))

    if not candidates:
        return None

    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


@pytest.mark.parametrize("label,env_overrides", SCENARIOS)
def test_full_run_produces_summary_without_failures(
    label: str, env_overrides: dict[str, str]
) -> None:
    """
    Slow end-to-end test that runs:

        python main.py -pdb TEST -fast

    under different APO_HOLO_MODE / PH_ENSEMBLE combinations and asserts:

      - main.py exits with code 0
      - the success sentinel "No proteins recorded as failed." is printed
      - at least one docking_score_summary.csv under docked/<RUN_ID>/TEST is
        updated by this run and contains at least one data row.
      - for APO/HOLO modes, the summary is located under docked/<RUN_ID>/TEST/APO/...
        or docked/<RUN_ID>/TEST/HOLO/... respectively.
    """
    if not _have_test_inputs():
        pytest.skip("Skipping full-run test: TEST.pdb or prepped test library missing.")

    assert MAIN_PY.exists(), f"main.py not found at {MAIN_PY}"

    env = os.environ.copy()
    env.update(env_overrides)
    env.update(MMGBSA_DISABLE_ENV)
    env["ATLAS_RUN_ID"] = TEST_RUN_ID
    # Force small test-mode library
    env["TEST_MODE_ENABLE"] = "dud"
    # Keep Python unbuffered so output ordering is sane in CI
    env["PYTHONUNBUFFERED"] = "1"

    cmd = [
        sys.executable,
        str(MAIN_PY),
        "-pdb",
        TEST_PDB_ID,
        "-fast",
        "--run-id",
        TEST_RUN_ID,
    ]

    start_time = time.time()
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        text=True,
        capture_output=True,
    )

    combined_output = (proc.stdout or "") + "\n" + (proc.stderr or "")

    # Helpful context if something goes wrong
    log_path = _extract_log_path(combined_output)
    print(f"=== {label} ===")
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
        f"{label}: non-zero exit {proc.returncode}\n"
        f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    assert (
        SUCCESS_SENTINEL in combined_output
    ), f"{label}: expected success sentinel in output."

    # Determine expected variant layout for this scenario
    mode = env_overrides.get("APO_HOLO_MODE", "").upper()
    expected_variant: str | None
    if mode in {"APO", "HOLO"}:
        expected_variant = mode
    else:
        expected_variant = None  # legacy: no variant in the docked path root

    summary_path = _find_recent_summary(start_time, expected_variant=expected_variant)
    assert summary_path is not None, (
        f"{label}: expected at least one docking_score_summary.csv under {DOCKED_ROOT} "
        f"updated by this run"
        + (f" (and under variant {expected_variant})" if expected_variant else "")
    )

    summary_stat = summary_path.stat()
    assert summary_stat.st_size > 0, f"{label}: {summary_path} is empty"
    assert summary_stat.st_mtime >= start_time, (
        f"{label}: {summary_path} was not updated by this run "
        f"(mtime={summary_stat.st_mtime}, start={start_time})"
    )

    # Basic content sanity: at least header + one data row
    non_empty_lines = [ln for ln in summary_path.read_text().splitlines() if ln.strip()]
    assert len(non_empty_lines) > 1, (
        f"{label}: {summary_path} has no data rows "
        f"(only {len(non_empty_lines)} non-empty lines)"
    )
