import csv
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "measure_artifacts.py"


def _write_bytes(path: Path, size: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * size)


def test_run_id_capture_and_diff_reports(tmp_path):
    work = tmp_path / "repo"
    work.mkdir(parents=True)

    docked = work / "docked"
    _write_bytes(docked / "a.txt", 10)
    _write_bytes(docked / "drop.tmp", 8)

    run_id = "TEST_RUN_001"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--run-id",
            run_id,
            "--capture",
            "baseline",
            "--root",
            str(docked),
            "--exclude",
            "**/*.tmp",
        ],
        cwd=str(work),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    _write_bytes(docked / "a.txt", 20)  # modified
    _write_bytes(docked / "b.pdbqt", 50)  # created

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--run-id",
            run_id,
            "--capture",
            "final",
            "--root",
            str(docked),
            "--exclude",
            "**/*.tmp",
        ],
        cwd=str(work),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    diff_cp = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--run-id",
            run_id,
            "--diff",
            "--root",
            str(docked),
            "--dir-depth",
            "3",
            "--top",
            "50",
            "--exclude",
            "**/*.tmp",
        ],
        cwd=str(work),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert "Created files: 1" in diff_cp.stdout
    assert "Modified files: 1" in diff_cp.stdout

    run_dir = work / "size_artifacts" / run_id
    reports = run_dir / "reports"

    assert (run_dir / "baseline.jsonl").exists()
    assert (run_dir / "final.jsonl").exists()
    assert (reports / "summary.json").exists()
    assert (reports / "by_extension.csv").exists()
    assert (reports / "by_dir.csv").exists()
    assert (reports / "top_files.csv").exists()
    assert (reports / "created_files.csv").exists()
    assert (reports / "modified_files.csv").exists()
    assert (reports / "deleted_files.csv").exists()

    summary = json.loads((reports / "summary.json").read_text(encoding="utf-8"))
    # footprint is created + modified (new size for modified)
    assert summary["totals"]["files"] == 2
    assert summary["totals"]["bytes"] == 70

    with (reports / "created_files.csv").open("r", encoding="utf-8", newline="") as handle:
        created_rows = list(csv.DictReader(handle))
    assert len(created_rows) == 1
    assert created_rows[0]["path"] == "b.pdbqt"
    assert int(created_rows[0]["size"]) == 50

    with (reports / "modified_files.csv").open("r", encoding="utf-8", newline="") as handle:
        modified_rows = list(csv.DictReader(handle))
    assert len(modified_rows) == 1
    assert modified_rows[0]["path"] == "a.txt"
    assert int(modified_rows[0]["old_size"]) == 10
    assert int(modified_rows[0]["new_size"]) == 20
