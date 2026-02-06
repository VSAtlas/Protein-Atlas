import subprocess
from pathlib import Path

import pytest

import protein_prep.os_utils as os_utils


def test_win_path_uses_wslpath_when_running_in_wsl(monkeypatch):
    calls = {}

    monkeypatch.setattr(os_utils, "_is_wsl", lambda: True)

    def fake_check_output(cmd, text):
        calls["cmd"] = cmd
        calls["text"] = text
        return "C:\\tmp\\input.pdb\n"

    monkeypatch.setattr(os_utils.subprocess, "check_output", fake_check_output)

    converted = os_utils._win_path("/tmp/input.pdb")

    assert converted == "C:\\tmp\\input.pdb"
    assert calls["cmd"] == ["wslpath", "-w", "/tmp/input.pdb"]
    assert calls["text"] is True


def test_powershell_builds_expected_command(monkeypatch):
    captured = {}

    def fake_run(cmd, stdout, stderr, text):
        captured["cmd"] = cmd
        captured["stdout"] = stdout
        captured["stderr"] = stderr
        captured["text"] = text
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(os_utils.subprocess, "run", fake_run)

    cp = os_utils._powershell("Write-Output hi")

    assert cp.returncode == 0
    assert captured["cmd"] == [
        "powershell.exe",
        "-NoProfile",
        "-NonInteractive",
        "-Command",
        "Write-Output hi",
    ]
    assert captured["text"] is True


def test_run_and_log_preserves_command_and_env(monkeypatch):
    captured = {}

    def fake_run(cmd, stdout, stderr, text, env):
        captured["cmd"] = cmd
        captured["stdout"] = stdout
        captured["stderr"] = stderr
        captured["text"] = text
        captured["env"] = env
        return subprocess.CompletedProcess(cmd, 0, "ok", "")

    monkeypatch.setattr(os_utils.subprocess, "run", fake_run)

    cp = os_utils._run_and_log(["phenix.pdbtools", "in.pdb"], env={"A": "1"}, check=True)

    assert cp.returncode == 0
    assert captured["cmd"] == ["phenix.pdbtools", "in.pdb"]
    assert captured["env"] == {"A": "1"}
    assert captured["text"] is True


def test_run_and_log_raises_on_nonzero_when_check_true(monkeypatch):
    def fake_run(cmd, stdout, stderr, text, env):
        return subprocess.CompletedProcess(cmd, 2, "out", "err")

    monkeypatch.setattr(os_utils.subprocess, "run", fake_run)

    with pytest.raises(subprocess.CalledProcessError):
        os_utils._run_and_log(["phenix.pdbtools"], env={}, check=True)


def test_as_path_returns_path_instance():
    as_path = os_utils._as_path("/tmp/example")
    assert isinstance(as_path, Path)
    assert as_path == Path("/tmp/example")
