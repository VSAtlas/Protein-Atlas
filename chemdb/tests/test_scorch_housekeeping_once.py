from __future__ import annotations

import logging
import threading
import time
from pathlib import Path

from post_docking.rescoring.scorch_housekeeping_once import run_housekeeping_once


def test_run_housekeeping_once_single_executor(tmp_path: Path) -> None:
    repo_root = tmp_path
    run_id = "rid_once"
    step = "pose_bust_dud"
    counter = {"n": 0}

    def action() -> None:
        counter["n"] += 1
        time.sleep(0.2)

    t = threading.Thread(
        target=run_housekeeping_once,
        kwargs={
            "repo_root": repo_root,
            "run_id": run_id,
            "step": step,
            "action": action,
            "logger": logging.getLogger("test-scorch-once-thread"),
            "wait_timeout_sec": 5.0,
            "poll_sec": 0.05,
        },
    )
    t.start()
    time.sleep(0.05)

    ok = run_housekeeping_once(
        repo_root=repo_root,
        run_id=run_id,
        step=step,
        action=action,
        logger=logging.getLogger("test-scorch-once-main"),
        wait_timeout_sec=5.0,
        poll_sec=0.05,
    )
    t.join(timeout=5.0)

    assert ok is True
    assert counter["n"] == 1
    done = (
        repo_root
        / "manifests"
        / run_id
        / "scorch_housekeeping"
        / f"{step}.done.json"
    )
    assert done.exists()


def test_run_housekeeping_once_wait_timeout(tmp_path: Path) -> None:
    repo_root = tmp_path
    run_id = "rid_timeout"
    step = "pose_bust_dud"
    marker_dir = repo_root / "manifests" / run_id / "scorch_housekeeping"
    marker_dir.mkdir(parents=True, exist_ok=True)
    (marker_dir / f"{step}.claim").write_text("busy\n", encoding="utf-8")
    ran = {"value": False}

    def action() -> None:
        ran["value"] = True

    ok = run_housekeeping_once(
        repo_root=repo_root,
        run_id=run_id,
        step=step,
        action=action,
        logger=logging.getLogger("test-scorch-once-timeout"),
        wait_timeout_sec=0.3,
        poll_sec=0.05,
    )
    assert ok is False
    assert ran["value"] is False
