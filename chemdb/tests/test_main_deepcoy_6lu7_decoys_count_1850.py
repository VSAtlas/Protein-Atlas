from __future__ import annotations

import contextlib
import os
import re
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import pytest
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"
PDB_ID = "6LU7"
PDB_PATH = REPO_ROOT / "input_pdbs" / f"{PDB_ID}.pdb"
TARGET_DIR = REPO_ROOT / "extracted_ligands" / "deepcoy" / "6LU7_P0DTD1"
DECOYS_PATH = TARGET_DIR / "deepcoy_decoys.smi"


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


@contextlib.contextmanager
def _temporary_decoys_per_active(value: int) -> None:
    candidates = [
        REPO_ROOT / "config.txt",
        REPO_ROOT.parent / "config.txt",
    ]
    config_paths = [p for p in candidates if p.exists()]
    originals: dict[Path, str] = {}
    try:
        for path in config_paths:
            text = path.read_text()
            originals[path] = text
            lines = text.splitlines()
            replaced = False
            new_lines: list[str] = []
            for line in lines:
                if re.match(r"^\\s*DEEPCOY_DECOYS_PER_ACTIVE\\s*=", line):
                    new_lines.append(f"DEEPCOY_DECOYS_PER_ACTIVE={value}")
                    replaced = True
                else:
                    new_lines.append(line)
            if not replaced:
                new_lines.append(f"DEEPCOY_DECOYS_PER_ACTIVE={value}")
            path.write_text(
                "\n".join(new_lines) + ("\n" if text.endswith("\n") else "")
            )
        yield
    finally:
        for path, content in originals.items():
            path.write_text(content)


@pytest.mark.network
@pytest.mark.integration
def test_main_6lu7_fast_test_fda_produces_1850_decoys(tmp_path: Path) -> None:  # noqa: ARG001
    if os.environ.get("RUN_NETWORK_TESTS") != "1":
        pytest.skip("set RUN_NETWORK_TESTS=1 to run")
    if not MAIN_PY.exists():
        pytest.skip(f"main.py not found at {MAIN_PY}")
    if not PDB_PATH.exists():
        pytest.skip(f"{PDB_ID} input pdb missing at {PDB_PATH}")

    try:
        resp = requests.get("https://data.rcsb.org/rest/v1/core/entry/6LU7", timeout=10)
        resp.raise_for_status()
    except Exception as exc:  # pragma: no cover - network guard path
        pytest.skip(f"network preflight failed: {exc}")

    run_id = f"{uuid.uuid4().hex[:8]}_deepcoy_1850"
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.pop("DEEPCOY_DECOYS_PER_ACTIVE", None)
    env.pop("DEEPCOY_NUM_SAMPLES", None)

    deepcoy_work_dir = REPO_ROOT / "DeepCoy_duds" / "deepcoy_work" / PDB_ID
    deepcoy_chunks_dir = TARGET_DIR.parent / ".deepcoy_chunks" / TARGET_DIR.name
    docked_dir = REPO_ROOT / "docked" / run_id
    cleanup_targets = [TARGET_DIR, deepcoy_work_dir, deepcoy_chunks_dir, docked_dir]
    pre_existing = {path: path.exists() for path in cleanup_targets}

    cmd = [
        sys.executable,
        "main.py",
        "-6lu7",
        "-fast",
        "-test-fda",
        "--run-id",
        run_id,
    ]
    try:
        with _temporary_decoys_per_active(50):
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
                    "main.py failed for -6lu7 -fast -test-fda\n"
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
        if not DECOYS_PATH.exists():
            pytest.fail(
                f"DeepCoy decoys missing: {DECOYS_PATH}\n"
                f"stdout tail:\n{_tail(stdout)}\n"
                f"stderr tail:\n{_tail(stderr)}"
            )

        decoy_count = _count_nonempty_lines(DECOYS_PATH)
        assert (
            decoy_count == 1850
        ), f"Expected 1850 merged decoys, got {decoy_count}"
    finally:
        for path, existed in pre_existing.items():
            if not existed and path.exists():
                shutil.rmtree(path, ignore_errors=True)
