from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from cli import qol_cli


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_atlas_help_sections_and_key_options(tmp_path) -> None:
    input_dir = tmp_path / "input_pdbs"
    overall_dir = tmp_path / "overall"
    input_dir.mkdir()
    overall_dir.mkdir()
    env = os.environ.copy()
    env["INPUT_DIR"] = str(input_dir)
    env["OVERALL_DIR"] = str(overall_dir)
    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "main.py"), "--help"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    assert proc.returncode == 0
    output = proc.stdout
    for section in (
        "Inputs:",
        "Run control:",
        "Docking modes:",
        "Ligand/library modes:",
        "pH/APO/HOLO options:",
        "Reporting:",
        "Debug/developer options:",
    ):
        assert section in output
    for option in (
        "--pdb",
        "--run-id",
        "--single",
        "--ph-ligand-mode",
        "--verify-tools",
        "--doctor",
        "--print-effective-config",
    ):
        assert option in output


def test_atlas_dev_verify_help_text() -> None:
    dev_proc = subprocess.run(
        [sys.executable, "-m", "cli.atlas_main_cli", "dev", "--help"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert dev_proc.returncode == 0
    dev_output = dev_proc.stdout
    assert "Developer maintenance workflows." in dev_output
    assert "verify" in dev_output
    assert "Run the repo quality gate." in dev_output

    verify_proc = subprocess.run(
        [sys.executable, "-m", "cli.atlas_main_cli", "dev", "verify", "--help"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert verify_proc.returncode == 0
    verify_output = verify_proc.stdout
    for option in ("--fix", "--smoke", "--full"):
        assert option in verify_output


def test_atlas_dev_verify_invokes_quality_gate_full(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(qol_cli.subprocess, "run", fake_run)

    rc = qol_cli._cmd_dev(["verify", "--full"])

    assert rc == 0
    assert captured["cmd"] == [
        "bash",
        str(REPO_ROOT / "tools" / "quality_gate.sh"),
        "--full",
    ]
    assert captured["kwargs"] == {"cwd": REPO_ROOT, "check": False}
