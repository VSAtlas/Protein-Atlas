from __future__ import annotations

import contextlib
import threading
from types import SimpleNamespace

from cli import run_lifecycle


class _FakeScorchQueue:
    def close_and_wait(self) -> dict[str, int]:
        return {"failed": 0}


class _FastClock:
    def __init__(self, start: float = 1000.0, step: float = 10.0) -> None:
        self._now = float(start)
        self._step = float(step)

    def time(self) -> float:
        self._now += self._step
        return self._now


def _dist_ctx() -> SimpleNamespace:
    return SimpleNamespace(
        enabled=True,
        task_id=0,
        task_count=4,
        task_min_id=0,
        leader_task_id=0,
        barrier_timeout_sec=600.0,
        barrier_poll_sec=1.0,
    )


def test_scorch_barrier_timeout_fail_open_runs_retention(
    monkeypatch,
    tmp_path,
) -> None:
    clock = _FastClock()
    monkeypatch.setattr(run_lifecycle.time, "time", clock.time)
    monkeypatch.setattr(run_lifecycle, "write_phase_marker", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        run_lifecycle,
        "wait_for_all_phase_markers",
        lambda *args, **kwargs: (False, [1, 2]),
    )
    monkeypatch.setattr(
        run_lifecycle,
        "_reconcile_missing_post_scorch_outputs",
        lambda *args, **kwargs: {"attempted": 0, "succeeded": 0, "failed": 0, "skipped_owner": 0},
    )
    monkeypatch.setattr(run_lifecycle, "acquire_global_slot", lambda *args, **kwargs: contextlib.nullcontext())
    monkeypatch.setenv("ATLAS_SCORCH_BARRIER_STRICT", "0")

    calls: list[list[tuple[str, str, str]]] = []

    def _record_retention(*args, **kwargs) -> None:
        calls.append(list(kwargs.get("pending_combo_list", []) or args[2]))

    monkeypatch.setattr(run_lifecycle, "_maybe_run_artifact_retention_for_combos", _record_retention)

    run_lifecycle.finalize_scheduler_and_retention(
        cfg={"OVERALL_DIR": str(tmp_path)},
        run_id="rid",
        global_start=0.0,
        run_scope_completion_times={},
        run_chunk_rebalance_count=0,
        dist_ctx=_dist_ctx(),
        scorch_queue_service=_FakeScorchQueue(),
        pending_coverage_refresh={},
        pending_retention_combos={("1ABC", "HOLO", "base")},
        retention_lock=threading.Lock(),
    )

    assert len(calls) == 1
    assert calls[0][0][0] == "1ABC"


def test_scorch_barrier_timeout_strict_skips_retention(
    monkeypatch,
    tmp_path,
) -> None:
    clock = _FastClock()
    monkeypatch.setattr(run_lifecycle.time, "time", clock.time)
    monkeypatch.setattr(run_lifecycle, "write_phase_marker", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        run_lifecycle,
        "wait_for_all_phase_markers",
        lambda *args, **kwargs: (False, [1, 2]),
    )
    monkeypatch.setattr(
        run_lifecycle,
        "_reconcile_missing_post_scorch_outputs",
        lambda *args, **kwargs: {"attempted": 0, "succeeded": 0, "failed": 0, "skipped_owner": 0},
    )
    monkeypatch.setattr(run_lifecycle, "acquire_global_slot", lambda *args, **kwargs: contextlib.nullcontext())
    monkeypatch.setenv("ATLAS_SCORCH_BARRIER_STRICT", "1")

    calls: list[list[tuple[str, str, str]]] = []

    def _record_retention(*args, **kwargs) -> None:
        calls.append(list(kwargs.get("pending_combo_list", []) or args[2]))

    monkeypatch.setattr(run_lifecycle, "_maybe_run_artifact_retention_for_combos", _record_retention)

    run_lifecycle.finalize_scheduler_and_retention(
        cfg={"OVERALL_DIR": str(tmp_path)},
        run_id="rid",
        global_start=0.0,
        run_scope_completion_times={},
        run_chunk_rebalance_count=0,
        dist_ctx=_dist_ctx(),
        scorch_queue_service=_FakeScorchQueue(),
        pending_coverage_refresh={},
        pending_retention_combos={("1ABC", "HOLO", "base")},
        retention_lock=threading.Lock(),
    )

    assert calls == []
