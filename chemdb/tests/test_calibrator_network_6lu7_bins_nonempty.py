from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]


def _count_nonempty(path: Path) -> int:
    return sum(1 for line in path.open() if line.strip())


@pytest.mark.network
@pytest.mark.slow
def test_calibrator_bins_nonempty_for_6lu7_network(tmp_path):  # noqa: ARG001
    if os.environ.get("RUN_NETWORK_TESTS") != "1":
        pytest.skip("set RUN_NETWORK_TESTS=1 to run")

    try:
        resp = requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/target.json"
            "?target_components__accession=P0DTD1&limit=1",
            timeout=10,
        )
        resp.raise_for_status()
    except Exception as exc:  # pragma: no cover - network preflight guard
        pytest.skip(f"network preflight failed: {exc}")

    run_tag = "pytest_network_6lu7"
    cmd = [
        sys.executable,
        "calibrator/chemdbl_calibrator.py",
        "--pdbs",
        "6LU7",
        "--timeout",
        "30",
        "--retries",
        "3",
        "--run-tag",
        run_tag,
    ]
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")
    subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )

    out_dir = REPO_ROOT / "extracted_ligands" / "6LU7_calibrator"
    strong_path = out_dir / "strong_binders.smi"
    weak_path = out_dir / "weak_binders.smi"
    non_path = out_dir / "non_binders.smi"
    for path in (strong_path, weak_path, non_path):
        assert path.is_file(), f"missing output file: {path}"
        assert _count_nonempty(path) > 0, f"expected non-empty file: {path}"

    log_path = REPO_ROOT / "calibrator" / "calibrator_logs" / run_tag / "6LU7.log"
    assert log_path.is_file()
