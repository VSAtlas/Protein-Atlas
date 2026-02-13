from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

RUN_ID = "deepcoy_repro_6LU7_testfda_1768000719"


def _create_fake_logs(tmp_path: Path) -> None:
    failed_dir = tmp_path / "failed"
    failed_dir.mkdir()
    failed_log = failed_dir / f"{RUN_ID}__6LU7__holo.log"
    failed_log.write_text(
        (
            "[FAILED PDB]\n"
            "  pdb_id        = 6LU7\n"
            "  variant_label = holo\n"
            "  mode          = holo\n"
            f"  pdb_file      = 6LU7.pdb run_id={RUN_ID}\n"
            "\n"
            "2026-01-09 12:00:00,000 - ERROR - root - DeepCoy generation failed: non-zero exit status\n"
            "Traceback (most recent call last):\n"
            '  File "deepcoy_driver.py", line 1, in <module>\n'
            '    raise RuntimeError("deepcoy crashed")\n'
            "RuntimeError: deepcoy crashed\n"
        ),
        encoding="utf-8",
    )

    deepcoy_dir = tmp_path / "DeepCoy_duds" / "deepcoy_work" / "6LU7"
    deepcoy_dir.mkdir(parents=True)
    (deepcoy_dir / "deepcoy.stderr.txt").write_text(
        (
            "Traceback (most recent call last):\n"
            '  File "run.py", line 10, in <module>\n'
            "    1/0\n"
            "ZeroDivisionError: division by zero\n"
        ),
        encoding="utf-8",
    )

    docked_dir = tmp_path / "docked" / RUN_ID / "6LU7"
    docked_dir.mkdir(parents=True)
    (docked_dir / "protein.log").write_text(
        (
            "2026-01-09 12:00:01,000 - INFO - root - start\n"
            "2026-01-09 12:00:02,000 - ERROR - root - Something bad happened in docking\n"
        ),
        encoding="utf-8",
    )


def _posix_path(value: str) -> str:
    return Path(value).as_posix()


def test_summarize_run_errors_hermetic(tmp_path: Path) -> None:
    _create_fake_logs(tmp_path)

    repo_root = Path(__file__).resolve().parents[2]
    script_path = repo_root / "analysis" / "summarize_run_errors.py"
    json_out = tmp_path / "out.json"

    cmd = [
        sys.executable,
        str(script_path),
        RUN_ID,
        "--deep",
        "--root",
        str(tmp_path),
        "--recent-hours",
        "0",
        "--top",
        "10",
        "--max-examples",
        "1",
        "--tail-kb",
        "64",
        "--json-out",
        str(json_out),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
        cwd=repo_root,
    )

    stdout = result.stdout
    assert "failed/" in stdout
    assert "DeepCoy_duds/deepcoy_work/6LU7/" in stdout
    assert "Traceback (most recent call last):" in stdout
    assert "unique_signatures:" in stdout or "raw_hits:" in stdout

    data = json.loads(json_out.read_text(encoding="utf-8"))
    examples = [example for sig in data["signatures"] for example in sig["examples"]]

    assert any("/failed/" in _posix_path(example["path"]) for example in examples)
    assert any(
        "/DeepCoy_duds/deepcoy_work/6LU7/" in _posix_path(example["path"])
        for example in examples
    )
    assert any(example.get("category") == "traceback" for example in examples)
