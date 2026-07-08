from __future__ import annotations

from concurrent.futures import Future
import logging

from cli.postrun_hooks_support import _resolve_scorch_decoy_prefix
from cli.postrun_hooks_scorch_queue import ScorchMixedQueueService


def test_scorch_decoy_prefix_honors_slurm_env(monkeypatch) -> None:
    monkeypatch.setenv("TEST_MODE_ENABLE", "dud")
    monkeypatch.setenv("DUD_PREFIX", "dud")
    monkeypatch.setenv("DECOY_PREFIX", "dud")

    assert _resolve_scorch_decoy_prefix({}) == "dud"


def test_sharded_queue_budgets_against_actual_pending_work() -> None:
    queue = ScorchMixedQueueService(
        run_id="rid",
        worker_cap=8,
        total_cpu=112,
        execution_mode="subprocess",
        logger=logging.getLogger("test.scorch.queue"),
        verbose=False,
    )
    cfg = {"ATLAS_SCORCH_SHARDS_ENABLE": "1"}
    try:
        with queue._lock:
            cpu_budget, free_budget, active_width = queue._cpu_budget_hints_locked(
                cfg=cfg,
                free_now=112,
            )
            assert (cpu_budget, free_budget, active_width) == (112, 112, 1)

            pending: Future[None] = Future()
            queue._futures[pending] = {"phase": "final"}
            cpu_budget, free_budget, active_width = queue._cpu_budget_hints_locked(
                cfg=cfg,
                free_now=112,
            )
            assert (cpu_budget, free_budget, active_width) == (56, 56, 2)
    finally:
        queue._pool.shutdown(wait=False, cancel_futures=True)


def test_queue_skips_missing_consensus_without_expected_hint(
    tmp_path,
    monkeypatch,
) -> None:
    queue = ScorchMixedQueueService(
        run_id="rid",
        worker_cap=2,
        total_cpu=16,
        execution_mode="subprocess",
        logger=logging.getLogger("test.scorch.queue"),
        verbose=False,
    )
    calls: list[str] = []
    monkeypatch.setattr(
        "cli.postrun_hooks_scorch_queue._maybe_run_scorch_rescore_for_pdb",
        lambda *args, **kwargs: calls.append("called") or 0,
    )
    cfg = {
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
        "TEST_MODE_ENABLE": "dud",
    }
    try:
        assert queue.submit(cfg=cfg, pdb_id="1ABC") is True
        assert calls == []
        with queue._lock:
            assert queue._submitted == 0
            assert queue._futures == {}
    finally:
        queue._pool.shutdown(wait=False, cancel_futures=True)


def test_queue_submits_missing_consensus_when_expected_hint_present(
    tmp_path,
    monkeypatch,
) -> None:
    queue = ScorchMixedQueueService(
        run_id="rid",
        worker_cap=2,
        total_cpu=16,
        execution_mode="subprocess",
        logger=logging.getLogger("test.scorch.queue"),
        verbose=False,
    )
    calls: list[str] = []
    monkeypatch.setattr(
        "cli.postrun_hooks_scorch_queue._maybe_run_scorch_rescore_for_pdb",
        lambda *args, **kwargs: calls.append("called") or 0,
    )
    cfg = {
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
        "TEST_MODE_ENABLE": "dud",
    }
    try:
        assert queue.submit(cfg=cfg, pdb_id="1ABC", allowed_count_hint=4) is True
        summary = queue.close_and_wait()
        assert summary["submitted"] == 1
        assert summary["failed"] == 0
        assert calls == ["called"]
    finally:
        try:
            queue._pool.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
