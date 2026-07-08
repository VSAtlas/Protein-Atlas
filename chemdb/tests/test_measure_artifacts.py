import csv
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "measure_artifacts.py"


@pytest.fixture(scope="module")
def measure_module():
    spec = importlib.util.spec_from_file_location("measure_artifacts", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_bytes(path: Path, size: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * size)


def test_analyze_root_totals_extensions_and_dirs(measure_module, tmp_path):
    root = tmp_path / "docked"
    _write_bytes(root / "file0", 4)
    _write_bytes(root / "a" / "x.txt", 10)
    _write_bytes(root / "a" / "y.bin", 6)
    _write_bytes(root / "a" / "b" / "z.TXT", 8)

    result = measure_module.analyze_root(
        root,
        dir_depth=1,
        top_n=10,
        follow_symlinks=False,
        exclude_patterns=[],
        quiet=True,
        verbose=False,
        dir_depth_full=False,
    )

    assert result["total_files"] == 4
    assert result["total_bytes"] == 28

    ext = result["extension_stats"]
    assert ext[".txt"]["count"] == 2
    assert ext[".txt"]["bytes"] == 18
    assert ext[".bin"]["count"] == 1
    assert ext[".bin"]["bytes"] == 6
    assert ext["(none)"]["count"] == 1
    assert ext["(none)"]["bytes"] == 4

    dirs = result["directory_stats"]
    assert dirs["."]["count"] == 4
    assert dirs["."]["bytes"] == 28
    assert dirs["a"]["count"] == 3
    assert dirs["a"]["bytes"] == 24


def test_exclude_pattern_skips_subtree(measure_module, tmp_path):
    root = tmp_path / "docked"
    _write_bytes(root / "run" / "microstates" / "a.dat", 5)
    _write_bytes(root / "run" / "keep" / "b.dat", 7)

    result = measure_module.analyze_root(
        root,
        dir_depth=3,
        top_n=10,
        follow_symlinks=False,
        exclude_patterns=["**/microstates/**"],
        quiet=True,
        verbose=False,
        dir_depth_full=False,
    )

    assert result["total_files"] == 1
    assert result["total_bytes"] == 7


def test_walk_tree_follow_symlinks_toggle(measure_module, tmp_path):
    root = tmp_path / "docked"
    _write_bytes(root / "real.dat", 9)
    target = root / "real.dat"
    link = root / "real_link.dat"

    try:
        os.symlink(target, link)
    except (NotImplementedError, OSError):
        pytest.skip("Symlink creation unsupported in this environment")

    nofollow = list(measure_module.walk_tree(root, follow_symlinks=False, exclude_patterns=[], quiet=True))
    follow = list(measure_module.walk_tree(root, follow_symlinks=True, exclude_patterns=[], quiet=True))

    nofollow_files = [r for r in nofollow if r[3]]
    follow_files = [r for r in follow if r[3]]

    assert len(nofollow_files) == 1
    assert len(follow_files) == 2


def test_cli_writes_single_root_reports(tmp_path):
    root = tmp_path / "docked"
    _write_bytes(root / "a" / "x.txt", 3)
    _write_bytes(root / "a" / "b" / "y.bin", 11)
    out_dir = tmp_path / "reports"

    cp = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--root",
            str(root),
            "--dir-depth",
            "2",
            "--top",
            "5",
            "--out",
            str(out_dir),
        ],
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )

    assert "Total files" in cp.stdout
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "by_extension.csv").exists()
    assert (out_dir / "by_dir.csv").exists()
    assert (out_dir / "top_files.csv").exists()

    summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["totals"]["files"] == 2
    assert summary["totals"]["bytes"] == 14


def test_cli_compare_mode_outputs_delta(tmp_path):
    root_a = tmp_path / "a"
    root_b = tmp_path / "b"
    _write_bytes(root_a / "keep.txt", 3)

    _write_bytes(root_b / "keep.txt", 5)
    _write_bytes(root_b / "new.bin", 2)

    out_dir = tmp_path / "compare_reports"

    cp = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--compare",
            str(root_a),
            str(root_b),
            "--out",
            str(out_dir),
            "--dir-depth",
            "3",
            "--top",
            "10",
        ],
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )

    assert "Compare (B - A)" in cp.stdout
    assert (out_dir / "compare.json").exists()
    assert (out_dir / "compare_by_extension.csv").exists()
    assert (out_dir / "compare_by_dir.csv").exists()

    compare = json.loads((out_dir / "compare.json").read_text(encoding="utf-8"))
    assert compare["delta"]["files"] == 1
    assert compare["delta"]["bytes"] == 4

    with (out_dir / "compare_by_extension.csv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    row_by_ext = {row["extension"]: row for row in rows}
    assert ".bin" in row_by_ext
    assert int(row_by_ext[".bin"]["delta_count"]) == 1
    assert int(row_by_ext[".bin"]["delta_bytes"]) == 2
