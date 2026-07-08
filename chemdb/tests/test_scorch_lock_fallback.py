from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

from post_docking.rescoring.scorch_lock_fallback import (
    build_env_python_retry_cmd,
    is_mamba_lock_error,
    run_timed_subprocess,
    run_with_lock_fallback,
)


def _cp(cmd: list[str], rc: int, *, stdout: str = "", stderr: str = "") -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(cmd, rc, stdout=stdout, stderr=stderr)


def test_is_mamba_lock_error_detects_known_signature() -> None:
    stderr = "error    libmamba Could not set lock (Resource temporarily unavailable)"
    assert is_mamba_lock_error(stderr, "") is True


def test_is_mamba_lock_error_ignores_non_lock_failure() -> None:
    stderr = "SCORCH failed: model load error"
    assert is_mamba_lock_error(stderr, "") is False


def test_build_env_python_retry_cmd_from_micromamba_command() -> None:
    script = Path("/opt/scorch/scorch.py")
    cmd = [
        "micromamba",
        "run",
        "-n",
        "scorch-env",
        "python",
        str(script),
        "--receptor",
        "r.pdbqt",
        "--ligand",
        "l.sdf",
        "--out",
        "o.csv",
    ]
    retry = build_env_python_retry_cmd(
        cmd,
        scorch_python=Path("/envs/scorch/bin/python"),
        scorch_script=script,
    )
    assert retry is not None
    assert retry[:2] == ["/envs/scorch/bin/python", str(script)]
    assert "--receptor" in retry
    assert "r.pdbqt" in retry


def test_build_env_python_retry_cmd_derives_prefix_python() -> None:
    script = Path("/opt/scorch/scorch.py")
    cmd = [
        "micromamba",
        "run",
        "-p",
        "/work/envs/scorch-env",
        "python",
        str(script),
        "--receptor",
        "r.pdbqt",
    ]
    retry = build_env_python_retry_cmd(
        cmd,
        scorch_python=None,
        scorch_script=script,
    )
    assert retry is not None
    assert retry[:2] == ["/work/envs/scorch-env/bin/python", str(script)]
    assert retry[-2:] == ["--receptor", "r.pdbqt"]


def test_run_with_lock_fallback_retries_once_on_lock_error(monkeypatch) -> None:
    calls: list[list[str]] = []
    script = Path("/opt/scorch/scorch.py")
    mamba_cmd = [
        "micromamba",
        "run",
        "-n",
        "scorch-env",
        "python",
        str(script),
        "--receptor",
        "r.pdbqt",
        "--ligand",
        "l.sdf",
        "--out",
        "o.csv",
    ]
    lock_stderr = "error libmamba Could not set lock (Resource temporarily unavailable)"
    responses = [
        _cp(mamba_cmd, 1, stderr=lock_stderr),
        _cp(["/envs/scorch/bin/python", str(script)], 0),
    ]

    def _fake_run(cmd, **kwargs):
        del kwargs
        cmd_list = [str(x) for x in cmd]
        calls.append(cmd_list)
        return responses[len(calls) - 1]

    monkeypatch.setattr("post_docking.rescoring.scorch_lock_fallback.subprocess.run", _fake_run)
    monkeypatch.delenv("ATLAS_SCORCH_DISABLE_LOCK_FALLBACK", raising=False)

    proc, retried = run_with_lock_fallback(
        cmd=mamba_cmd,
        run_kwargs={"capture_output": True, "text": True},
        logger=logging.getLogger("test.scorch.lock"),
        component="[scorch-rescore]",
        source="vina",
        stage="vina_best",
        scorch_python=Path("/envs/scorch/bin/python"),
        scorch_script=script,
    )
    assert retried is True
    assert proc.returncode == 0
    assert len(calls) == 2
    assert calls[0][0] == "micromamba"
    assert calls[1][0] == "/envs/scorch/bin/python"


def test_run_with_lock_fallback_does_not_retry_non_lock_error(monkeypatch) -> None:
    calls: list[list[str]] = []
    script = Path("/opt/scorch/scorch.py")
    mamba_cmd = ["micromamba", "run", "-n", "scorch-env", "python", str(script)]

    def _fake_run(cmd, **kwargs):
        del kwargs
        cmd_list = [str(x) for x in cmd]
        calls.append(cmd_list)
        return _cp(cmd_list, 1, stderr="generic failure")

    monkeypatch.setattr("post_docking.rescoring.scorch_lock_fallback.subprocess.run", _fake_run)
    monkeypatch.delenv("ATLAS_SCORCH_DISABLE_LOCK_FALLBACK", raising=False)

    proc, retried = run_with_lock_fallback(
        cmd=mamba_cmd,
        run_kwargs={"capture_output": True, "text": True},
        logger=logging.getLogger("test.scorch.lock"),
        component="[scorch-rescore]",
        source="vina",
        stage="vina_best",
        scorch_python=Path("/envs/scorch/bin/python"),
        scorch_script=script,
    )
    assert retried is False
    assert proc.returncode == 1
    assert len(calls) == 1


def test_run_with_lock_fallback_can_be_disabled(monkeypatch) -> None:
    calls: list[list[str]] = []
    script = Path("/opt/scorch/scorch.py")
    mamba_cmd = ["micromamba", "run", "-n", "scorch-env", "python", str(script)]
    lock_stderr = "error libmamba Could not set lock (Resource temporarily unavailable)"

    def _fake_run(cmd, **kwargs):
        del kwargs
        cmd_list = [str(x) for x in cmd]
        calls.append(cmd_list)
        return _cp(cmd_list, 1, stderr=lock_stderr)

    monkeypatch.setattr("post_docking.rescoring.scorch_lock_fallback.subprocess.run", _fake_run)
    monkeypatch.setenv("ATLAS_SCORCH_DISABLE_LOCK_FALLBACK", "1")

    proc, retried = run_with_lock_fallback(
        cmd=mamba_cmd,
        run_kwargs={"capture_output": True, "text": True},
        logger=logging.getLogger("test.scorch.lock"),
        component="[scorch-rescore]",
        source="vina",
        stage="vina_best",
        scorch_python=Path("/envs/scorch/bin/python"),
        scorch_script=script,
    )
    assert retried is False
    assert proc.returncode == 1
    assert len(calls) == 1


def test_run_with_lock_fallback_runtime_retry_uses_threads_one(monkeypatch) -> None:
    calls: list[list[str]] = []
    script = Path("/opt/scorch/scorch.py")
    cmd = [
        "micromamba",
        "run",
        "-n",
        "scorch-env",
        "python",
        str(script),
        "--threads",
        "4",
        "--receptor",
        "r.pdbqt",
        "--ligand",
        "l.sdf",
        "--out",
        "o.csv",
    ]
    runtime_stderr = "joblib.externals.loky.process_executor.BrokenProcessPool: A task has failed to un-serialize."
    responses = [
        _cp(cmd, 1, stderr=runtime_stderr),
        _cp(cmd, 0, stderr=""),
    ]

    def _fake_run(run_cmd, **kwargs):
        del kwargs
        cmd_list = [str(x) for x in run_cmd]
        calls.append(cmd_list)
        return responses[len(calls) - 1]

    monkeypatch.setattr("post_docking.rescoring.scorch_lock_fallback.subprocess.run", _fake_run)
    monkeypatch.delenv("ATLAS_SCORCH_DISABLE_RUNTIME_RETRY", raising=False)
    monkeypatch.delenv("ATLAS_SCORCH_DISABLE_LOCK_FALLBACK", raising=False)

    proc, retried = run_with_lock_fallback(
        cmd=cmd,
        run_kwargs={"capture_output": True, "text": True},
        logger=logging.getLogger("test.scorch.lock"),
        component="[scorch-rescore]",
        source="vina",
        stage="vina_best",
        scorch_python=Path("/envs/scorch/bin/python"),
        scorch_script=script,
    )
    assert retried is True
    assert proc.returncode == 0
    assert len(calls) == 2
    assert "--threads" in calls[1]
    assert calls[1][calls[1].index("--threads") + 1] == "1"


def test_run_with_lock_fallback_runtime_retry_can_be_disabled(monkeypatch) -> None:
    calls: list[list[str]] = []
    script = Path("/opt/scorch/scorch.py")
    cmd = [
        "micromamba",
        "run",
        "-n",
        "scorch-env",
        "python",
        str(script),
        "--threads",
        "4",
    ]
    runtime_stderr = "Fatal Python error: Cannot recover from MemoryErrors while normalizing exceptions."

    def _fake_run(run_cmd, **kwargs):
        del kwargs
        cmd_list = [str(x) for x in run_cmd]
        calls.append(cmd_list)
        return _cp(cmd_list, 1, stderr=runtime_stderr)

    monkeypatch.setattr("post_docking.rescoring.scorch_lock_fallback.subprocess.run", _fake_run)
    monkeypatch.setenv("ATLAS_SCORCH_DISABLE_RUNTIME_RETRY", "1")
    monkeypatch.delenv("ATLAS_SCORCH_DISABLE_LOCK_FALLBACK", raising=False)

    proc, retried = run_with_lock_fallback(
        cmd=cmd,
        run_kwargs={"capture_output": True, "text": True},
        logger=logging.getLogger("test.scorch.lock"),
        component="[scorch-rescore]",
        source="vina",
        stage="vina_best",
        scorch_python=Path("/envs/scorch/bin/python"),
        scorch_script=script,
    )
    assert retried is False
    assert proc.returncode == 1
    assert len(calls) == 1


def test_run_timed_subprocess_returns_timeout_code() -> None:
    proc = run_timed_subprocess(
        [sys.executable, "-c", "import time; time.sleep(5)"],
        run_kwargs={"capture_output": True, "text": True},
        timeout_sec=0.1,
    )
    assert proc.returncode == 124
    assert "timeout" in (proc.stderr or "").lower()
