from __future__ import annotations

import threading
import time

import main
from docking.global_scheduler import (
    GlobalAdmissionScheduler,
    acquire_global_cores,
    acquire_global_slot,
    sem_free_slots,
)



def test_resolve_global_scheduler_plan_adaptive_defaults() -> None:
    cfg = {
        "ENABLE_GLOBAL_SCHEDULER": True,
        "GLOBAL_ADMISSION_POLICY": "adaptive",
        "GLOBAL_SCHEDULER_CPUS": 12,
        "GLOBAL_MIN_PROTEINS": 1,
        "GLOBAL_MAX_PROTEINS": 0,
    }
    plan = main._resolve_global_scheduler_plan(cfg, pdb_count=20, cpu=8)

    assert plan["enabled"] is True
    assert plan["policy"] == "adaptive"
    assert plan["scheduler_cpus"] == 12
    assert plan["max_parallel"] == 12
    assert plan["initial_target"] == 6



def test_resolve_global_scheduler_plan_fixed_with_caps() -> None:
    cfg = {
        "ENABLE_GLOBAL_SCHEDULER": True,
        "GLOBAL_ADMISSION_POLICY": "fixed",
        "GLOBAL_SCHEDULER_CPUS": 16,
        "GLOBAL_MIN_PROTEINS": 2,
        "GLOBAL_MAX_PROTEINS": 4,
    }
    plan = main._resolve_global_scheduler_plan(cfg, pdb_count=10, cpu=32)

    assert plan["policy"] == "fixed"
    assert plan["min_parallel"] == 2
    assert plan["max_parallel"] == 4
    assert plan["initial_target"] == 4



def test_acquire_global_slot_uses_runtime_semaphore() -> None:
    sem = threading.BoundedSemaphore(2)
    cfg = {"GLOBAL_DOCKING_SEM": sem}

    start_free = sem_free_slots(sem, default=2)
    assert start_free == 2

    with acquire_global_slot(cfg):
        mid_free = sem_free_slots(sem, default=2)
        assert mid_free == 1

    end_free = sem_free_slots(sem, default=2)
    assert end_free == 2


def test_global_admission_backfills_smaller_request() -> None:
    scheduler = GlobalAdmissionScheduler(total_cores=4)
    cfg = {"GLOBAL_DOCKING_SCHEDULER": scheduler}

    held = scheduler.acquire(cores=3, task_id="held")
    order: list[str] = []

    def _run(name: str, cores: int, hold_sec: float) -> None:
        with acquire_global_slot(cfg, cores=cores, task_id=name):
            order.append(name)
            time.sleep(hold_sec)

    t_big = threading.Thread(target=_run, args=("big", 3, 0.02))
    t_big.start()
    time.sleep(0.02)

    t_small = threading.Thread(target=_run, args=("small", 1, 0.02))
    t_small.start()
    t_small.join(timeout=1.0)
    assert not t_small.is_alive()
    assert order == ["small"]

    scheduler.release(held)
    t_big.join(timeout=1.0)
    assert not t_big.is_alive()
    assert order == ["small", "big"]


def test_global_admission_elastic_grants_fit_available() -> None:
    scheduler = GlobalAdmissionScheduler(total_cores=4)
    cfg = {"GLOBAL_DOCKING_SCHEDULER": scheduler}

    held = scheduler.acquire(cores=3, min_cores=3, task_id="held")
    try:
        with acquire_global_cores(
            cfg,
            cores=4,
            min_cores=1,
            priority=0,
            task_id="elastic",
            task_type="scorch_score_chunk",
        ) as granted:
            assert granted == 1
    finally:
        scheduler.release(held)


def test_global_admission_tail_prefers_short_fitting_request() -> None:
    scheduler = GlobalAdmissionScheduler(total_cores=4)
    cfg = {"GLOBAL_DOCKING_SCHEDULER": scheduler}

    held = scheduler.acquire(cores=4, min_cores=4, task_id="held")
    order: list[str] = []
    gate = threading.Barrier(3)

    def _run(name: str, est_sec: float) -> None:
        gate.wait(timeout=1.0)
        with acquire_global_cores(
            cfg,
            cores=2,
            min_cores=1,
            est_duration_sec=est_sec,
            task_id=name,
            task_type="dock_stage_ligand",
        ):
            order.append(name)
            time.sleep(0.01)

    t_slow = threading.Thread(target=_run, args=("slow", 9.0))
    t_fast = threading.Thread(target=_run, args=("fast", 0.2))
    t_slow.start()
    t_fast.start()
    gate.wait(timeout=1.0)
    time.sleep(0.02)
    scheduler.release(held)

    t_slow.join(timeout=1.0)
    t_fast.join(timeout=1.0)
    assert not t_slow.is_alive()
    assert not t_fast.is_alive()
    assert order and order[0] == "fast"


def test_global_admission_tail_caps_large_request() -> None:
    scheduler = GlobalAdmissionScheduler(total_cores=8)
    cfg = {"GLOBAL_DOCKING_SCHEDULER": scheduler}

    held = scheduler.acquire(cores=4, min_cores=4, task_id="held")
    try:
        with acquire_global_cores(
            cfg,
            cores=8,
            min_cores=1,
            est_duration_sec=2.0,
            task_id="large",
            task_type="scorch_score_chunk",
        ) as granted:
            assert granted <= 2
    finally:
        scheduler.release(held)
