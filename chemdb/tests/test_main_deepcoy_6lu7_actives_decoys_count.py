from __future__ import annotations

import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_DIR = REPO_ROOT / "extracted_ligands" / "deepcoy" / "6LU7_P0DTD1"
ACTIVES_PATH = TARGET_DIR / "6LU7_P0DTD1_actives.smi"
DECOYS_PATH = TARGET_DIR / "deepcoy_decoys.smi"
RAW_ACTIVES_PATH = TARGET_DIR / "6LU7_P0DTD1_actives.raw.smi"


def _count_nonempty_lines(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                count += 1
    return count


def _tail(text: str, limit: int = 2000) -> str:
    return text[-limit:] if text else ""


@pytest.mark.network
@pytest.mark.integration
def test_main_force_deepcoy_6lu7_produces_38_actives_and_38_decoys(
    tmp_path: Path,
) -> None:  # noqa: ARG001
    if os.environ.get("RUN_NETWORK_TESTS") != "1":
        pytest.skip("set RUN_NETWORK_TESTS=1 to run")

    try:
        resp = requests.get("https://data.rcsb.org/rest/v1/core/entry/6LU7", timeout=10)
        resp.raise_for_status()
    except Exception as exc:  # pragma: no cover - network guard path
        pytest.skip(f"network preflight failed: {exc}")

    run_id = f"{uuid.uuid4().hex[:8]}_deepcoy_actives_per_decoy"
    env = os.environ.copy()
    env["DEEPCOY_DECOYS_PER_ACTIVE"] = "1"
    env.setdefault("PYTHONUNBUFFERED", "1")

    cmd = [
        sys.executable,
        "main.py",
        "-6lu7",
        "-fast",
        "--test-fda",
        "--run-id",
        run_id,
        "--force-deepcoy",
    ]
    try:
        result = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            env=env,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        pytest.fail(
            "main.py failed for --force-deepcoy\n"
            f"returncode: {exc.returncode}\n"
            f"stdout tail:\n{_tail(exc.stdout)}\n"
            f"stderr tail:\n{_tail(exc.stderr)}"
        )

    stdout = result.stdout or ""
    stderr = result.stderr or ""

    if not TARGET_DIR.is_dir():
        pytest.fail(
            f"DeepCoy output dir missing: {TARGET_DIR}\n"
            f"stdout tail:\n{_tail(stdout)}\n"
            f"stderr tail:\n{_tail(stderr)}"
        )

    missing = [p for p in (ACTIVES_PATH, DECOYS_PATH) if not p.exists()]
    if missing:
        missing_str = ", ".join(str(p) for p in missing)
        pytest.fail(
            f"DeepCoy output files missing: {missing_str}\n"
            f"stdout tail:\n{_tail(stdout)}\n"
            f"stderr tail:\n{_tail(stderr)}"
        )

    n_actives = _count_nonempty_lines(ACTIVES_PATH)
    n_decoys = _count_nonempty_lines(DECOYS_PATH)
    n_raw = (
        _count_nonempty_lines(RAW_ACTIVES_PATH) if RAW_ACTIVES_PATH.exists() else None
    )

    if not (n_actives == 38 and n_decoys == 38):
        msg = f"DeepCoy counts unexpected: actives={n_actives}, decoys={n_decoys}"
        if n_raw is not None:
            msg += f", raw_actives={n_raw}"
        msg += (
            f"\nactives path: {ACTIVES_PATH}\n"
            f"decoys path: {DECOYS_PATH}\n"
            f"stdout tail:\n{_tail(stdout)}\n"
            f"stderr tail:\n{_tail(stderr)}"
        )
        pytest.fail(msg)
