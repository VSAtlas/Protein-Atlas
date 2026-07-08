from __future__ import annotations

import logging
from pathlib import Path

import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def test_run_pose_bust_uses_chunked_posebusters(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    class _Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    def _fake_run(cmd, capture_output, text):  # type: ignore[no-untyped-def]
        captured["cmd"] = list(cmd)
        captured["capture_output"] = capture_output
        captured["text"] = text
        return _Proc()

    monkeypatch.setattr(rescoring_scorch.subprocess, "run", _fake_run)

    ok = rescoring_scorch._run_pose_bust(
        run_id="rid",
        repo_root=tmp_path,
        docked_root=tmp_path / "docked",
        post_docked_root=tmp_path / "post_docked",
        overwrite=False,
        max_workers=2,
        logger=logging.getLogger("test-scorch-pose-bust"),
    )

    assert ok is True
    cmd = captured.get("cmd")
    assert isinstance(cmd, list)
    assert "--posebusters-chunk" in cmd
    idx = cmd.index("--posebusters-chunk")
    assert cmd[idx + 1] == "96"
    assert "--posebusters-max-workers" in cmd


def test_launch_pose_bust_async_uses_popen(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    class _Proc:
        pid = 4321

    def _fake_popen(*args, **kwargs):  # type: ignore[no-untyped-def]
        captured["args"] = args
        captured["kwargs"] = kwargs
        return _Proc()

    monkeypatch.setattr(rescoring_scorch.subprocess, "Popen", _fake_popen)

    pid = rescoring_scorch._launch_pose_bust_async(
        run_id="rid_async",
        repo_root=tmp_path,
        docked_root=tmp_path / "docked",
        post_docked_root=tmp_path / "post_docked",
        overwrite=True,
        max_workers=9,
        logger=logging.getLogger("test-scorch-pose-bust"),
    )

    assert pid == 4321
    popen_args = captured.get("args")
    assert isinstance(popen_args, tuple) and popen_args
    cmd = list(popen_args[0])
    assert "--posebusters-max-workers" in cmd
    max_idx = cmd.index("--posebusters-max-workers")
    assert cmd[max_idx + 1] == "2"
    assert "--overwrite" in cmd
