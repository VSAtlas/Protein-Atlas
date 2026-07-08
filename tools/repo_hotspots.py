#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path


EXCLUDE_DIRS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    ".home",
    ".probe_regression",
    ".cache_test",
    "node_modules",
    "docked",
    "post_docked",
    "input_pdbs",
}


@dataclass(frozen=True)
class Hotspot:
    path: Path
    bytes_size: int
    lines: int


def _iter_py_files(root: Path) -> list[Path]:
    results: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [name for name in dirnames if name not in EXCLUDE_DIRS]
        for name in filenames:
            if not name.endswith(".py"):
                continue
            results.append(Path(dirpath) / name)
    return results


def _line_count(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for _ in handle)


def _collect(root: Path) -> list[Hotspot]:
    hotspots: list[Hotspot] = []
    for path in _iter_py_files(root):
        try:
            size = path.stat().st_size
            lines = _line_count(path)
        except OSError:
            continue
        hotspots.append(Hotspot(path=path.relative_to(root), bytes_size=size, lines=lines))
    hotspots.sort(key=lambda x: (x.bytes_size, x.lines), reverse=True)
    return hotspots


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="List large Python files to prioritize simplification.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of rows to print.",
    )
    parser.add_argument(
        "--min-bytes",
        type=int,
        default=50_000,
        help="Minimum file size (bytes) to include.",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    rows = [h for h in _collect(root) if h.bytes_size >= int(args.min_bytes)]
    rows = rows[: max(1, int(args.top))]

    if not rows:
        print("No files matched the requested threshold.")
        return 0

    print("bytes\tlines\tpath")
    for item in rows:
        print(f"{item.bytes_size}\t{item.lines}\t{item.path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
