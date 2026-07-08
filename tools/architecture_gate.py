#!/usr/bin/env python3
"""Run architecture-oriented quality gates with repo-local baselines."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
XENON_BASELINE = REPO_ROOT / "tools" / "xenon_baseline.txt"
XENON_IGNORES = [
    ".cache",
    ".cache_test",
    ".direnv",
    ".git",
    ".home",
    "__pycache__",
    "src/ml/cache",
    "src/ml/data",
    "src/ml/outputs",
]


def _read_baseline(path: Path) -> list[str]:
    entries: list[str] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            entries.append(stripped)
    return entries


def _run(command: list[str], display: str | None = None) -> int:
    shown = display or " ".join(command)
    print(f"[architecture] {shown}", flush=True)
    return subprocess.run(command, cwd=REPO_ROOT).returncode


def _run_tach() -> int:
    return _run(["tach", "check"])


def _run_xenon() -> int:
    baseline = _read_baseline(XENON_BASELINE)
    command = [
        "xenon",
        "--max-absolute",
        "B",
        "--max-modules",
        "B",
        "--max-average",
        "A",
        "--ignore",
        ",".join(XENON_IGNORES),
    ]
    if baseline:
        command.extend(["--exclude", ",".join(baseline)])
    command.append("src")
    return _run(
        command,
        display=(
            "xenon --max-absolute B --max-modules B --max-average A "
            f"--ignore {','.join(XENON_IGNORES)} "
            f"--exclude <{len(baseline)} baseline paths> src"
        ),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-tach", action="store_true")
    parser.add_argument("--skip-xenon", action="store_true")
    args = parser.parse_args(argv)

    exit_code = 0
    if not args.skip_tach:
        exit_code = _run_tach() or exit_code
    if not args.skip_xenon:
        exit_code = _run_xenon() or exit_code
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
