from __future__ import annotations

import contextlib
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

# Root of the repo (chemdb/tests/ -> chemdb -> repo root)
REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"
PDB_ID = "6LU7"
PDB_PATH = REPO_ROOT / "input_pdbs" / f"{PDB_ID}.pdb"
MIN_DECOYS = 5

LOG_PATH_RE = re.compile(r"log_file=(?P<path>\\S+)")
DEEPCOY_DIR_RE = re.compile(r"Decoys should be generated in:\\s*(?P<dir>.+)")
DEEPCOY_WARN_RE = re.compile(r"\\[deepcoy\\]\\s+WARNING: DeepCoy generated\\s+(?P<count>\\d+)\\s+decoys")
ACTIVES_PATTERNS = [
    re.compile(r"-> Found\\s+(?P<count>\\d+)\\s+unique SMILES strings\\."),
    re.compile(r"\\[deepcoy\\.actives\\.sdf\\].*\\bn=(?P<count>\\d+)\\b"),
]


def _resolve_log_path(output: str, run_id: str) -> Path:
    match = LOG_PATH_RE.search(output)
    if match:
        candidate = Path(match.group("path"))
        if candidate.exists():
            return candidate
    return REPO_ROOT / "logs" / f"main_{run_id}.log"


def _extract_decoy_dir(log_text: str) -> Path | None:
    match = DEEPCOY_DIR_RE.search(log_text)
    if match:
        candidate = Path(match.group("dir").strip())
        if candidate.exists():
            return candidate
    return None


def _find_recent_deepcoy_dir(start_time: float) -> Path | None:
    base = REPO_ROOT / "extracted_ligands" / "deepcoy"
    if not base.exists():
        return None

    candidates: list[tuple[float, Path]] = []
    for child in base.iterdir():
        if not child.is_dir():
            continue
        try:
            stat = child.stat()
        except FileNotFoundError:
            continue
        if stat.st_mtime >= start_time:
            candidates.append((stat.st_mtime, child))

    if not candidates:
        return None

    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1]


def _count_nonempty_lines(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                count += 1
    return count


@contextlib.contextmanager
def _temporary_decoys_per_active_one():
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
                    new_lines.append("DEEPCOY_DECOYS_PER_ACTIVE=1")
                    replaced = True
                else:
                    new_lines.append(line)
            if not replaced:
                new_lines.append("DEEPCOY_DECOYS_PER_ACTIVE=1")
            path.write_text("\n".join(new_lines) + ("\n" if text.endswith("\n") else ""))
        yield
    finally:
        for path, content in originals.items():
            path.write_text(content)


@pytest.mark.slow
def test_deepcoy_generation_for_6lu7_produces_decoys() -> None:
    """
    Regression harness for:
        python main.py -pdb 6LU7 --run-id deepcoyBIG --force-deepcoy -fast
    Uses a unique run id and adds -test-fda/-fast to keep inputs bounded while
    exercising DeepCoy autogen + DUD merge. Forces decoys-per-active=1 for speed.
    """
    if not MAIN_PY.exists():
        pytest.skip(f"main.py not found at {MAIN_PY}")
    if not PDB_PATH.exists():
        pytest.skip(f"{PDB_ID} input pdb missing at {PDB_PATH}")

    run_id = f"deepcoy_repro_{PDB_ID}_testfda_dpa1_{int(time.time())}"
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.update(
        {
            "MMGBSA_ENABLED": "false",
            "MMGBSA_MD_ENABLED": "false",
            "MMGBSA_MMPBSA_ENABLED": "false",
            "MMGBSA_MMPBSA_RUN": "false",
            "MD_FIVE_REPLICATE": "false",
        }
    )

    cmd = [
        sys.executable,
        str(MAIN_PY),
        "-pdb",
        PDB_ID,
        "-fast",
        "--run-id",
        run_id,
        "--force-deepcoy",
        "-test-fda",
    ]

    start_time = time.time()
    with _temporary_decoys_per_active_one():
        proc = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            env=env,
            text=True,
            capture_output=True,
        )
    output = (proc.stdout or "") + "\n" + (proc.stderr or "")

    assert proc.returncode == 0, (
        f"main.py exited with {proc.returncode} for --force-deepcoy\n"
        f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )

    log_path = _resolve_log_path(output, run_id)
    assert log_path.exists(), f"Log file not found (expected {log_path})"
    log_text = log_path.read_text()

    assert "## 5. Executing DeepCoy workflow..." in log_text, "DeepCoy workflow was not invoked"
    assert "simplified fallback" not in log_text.lower(), "DeepCoy fell back to simplified preprocessing path"
    assert (
        "[deepcoy] WARNING: DeepCoy generated" not in log_text
    ), "DeepCoy emitted underproduction warning"

    warn_match = DEEPCOY_WARN_RE.search(log_text)
    decoy_dir = _extract_decoy_dir(log_text)
    if decoy_dir is None or not decoy_dir.exists():
        decoy_dir = _find_recent_deepcoy_dir(start_time)

    assert decoy_dir is not None and decoy_dir.exists(), "DeepCoy decoy output dir not found"

    decoys_smi = decoy_dir / "deepcoy_decoys.smi"
    assert decoys_smi.exists(), f"deepcoy_decoys.smi missing in {decoy_dir}"

    decoy_lines = _count_nonempty_lines(decoys_smi)
    reported_count = warn_match.group("count") if warn_match else None

    n_actives = None
    for pattern in ACTIVES_PATTERNS:
        match = pattern.search(log_text)
        if match:
            try:
                n_actives = int(match.group("count"))
            except (KeyError, ValueError):
                n_actives = None
            break

    expected_min = MIN_DECOYS
    if n_actives is not None:
        expected_min = max(expected_min, int(0.5 * n_actives))

    failure_msg = f"DeepCoy generated {decoy_lines} decoys; expected at least {expected_min}"
    if reported_count is not None:
        failure_msg += f" (warning reported {reported_count})"
    assert decoy_lines >= expected_min, failure_msg
