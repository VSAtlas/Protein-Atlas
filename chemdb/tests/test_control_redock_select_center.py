import os
import re
import sys
import time
import subprocess
from pathlib import Path

import pytest
import yaml

PDBS_OF_INTEREST = [
    "1SYN",
    "1UYG",
    "1XL2",
    "2I78",
    "2OF2",
    "2OWB",
    "TEST",
]


def _run_id_from_stdout(stdout: str) -> str | None:
    """
    Extract the run_id from the standard output of main.py to avoid collisions
    with other concurrently running processes that may create manifests.
    """
    match = re.search(r"run_id=([\w\-.]+)", stdout)
    if match:
        return match.group(1)
    return None


def _latest_manifest_after(manifests_dir: Path, ts: float) -> str | None:
    """
    Choose the most recently modified manifest directory created after the
    provided timestamp. This is a fallback when stdout parsing fails.
    """
    candidates: list[tuple[float, str]] = []
    for d in manifests_dir.iterdir():
        if not d.is_dir():
            continue
        try:
            mtime = d.stat().st_mtime
        except OSError:
            continue
        if mtime >= ts:
            candidates.append((mtime, d.name))
    if not candidates:
        return None
    candidates.sort(reverse=True)
    return candidates[0][1]


@pytest.mark.slow
def test_control_redock_sets_method_control_for_dud_targets(tmp_path):
    root = Path(__file__).resolve().parents[2]
    manifests_dir = root / "manifests"
    manifests_dir.mkdir(exist_ok=True)

    start_ts = time.time()
    pdbs_arg = " ".join(PDBS_OF_INTEREST)
    cmd = [
        sys.executable,
        os.fspath(root / "main.py"),
        "-fast",
        "-pdbs",
        pdbs_arg,
        "--single",
        "dexamethasone",
    ]

    result = subprocess.run(
        cmd,
        cwd=root,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, (
        "main.py failed with non-zero exit code\n"
        f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
    )

    run_id = _run_id_from_stdout(result.stdout)
    if not run_id:
        run_id = _latest_manifest_after(manifests_dir, start_ts)
    assert run_id, "Failed to resolve run_id for control-redock test run"

    manifest_path = manifests_dir / run_id / "run_manifest.yaml"
    assert manifest_path.is_file(), f"Manifest not found at {manifest_path}"

    with manifest_path.open("r") as f:
        data = yaml.safe_load(f) or {}
    proteins = data.get("proteins", {}) or {}

    required = set(PDBS_OF_INTEREST)
    seen_control = {pdb_id: False for pdb_id in required}

    for entry in proteins.values():
        pdb_id = entry.get("pdb_id")
        if pdb_id not in required:
            continue

        status = entry.get("status")
        stages = entry.get("stages", {}) or {}
        pocket = stages.get("pocket_detection", {}) or {}
        details = pocket.get("details", {}) or {}
        method = details.get("method")

        if status == "completed" and method == "control":
            seen_control[pdb_id] = True

    missing = [p for p, ok in seen_control.items() if not ok]

    assert not missing, (
        "Expected control redocking to be used for pocket_detection "
        "for all DUD targets, but these PDBs had no completed entry "
        "with method='control': "
        + ", ".join(sorted(missing))
        + f"\nRun ID: {run_id}"
    )
