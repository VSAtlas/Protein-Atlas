from __future__ import annotations

import math

import post_docking.rescoring.rescoring_scorch as rescoring_scorch
from post_docking.rescoring import scorch_main
from post_docking.rescoring import scorch_parallel


def test_allocate_scorch_parallelism_aggressive_low_pressure() -> None:
    jobs, threads, util = rescoring_scorch._allocate_scorch_parallelism(
        cpu_budget=112,
        free_cores=112,
        backlog=16,
        allowed_count=128,
        scheduler_queue_depth=0,
        local_queue_depth=0,
        profile="aggressive",
    )
    assert util == 0.97
    assert jobs == 16
    assert threads == 4
    assert jobs * threads <= math.floor(util * 112)


def test_allocate_scorch_parallelism_high_pressure_reduces_utilization() -> None:
    jobs_aggr, threads_aggr, util_aggr = rescoring_scorch._allocate_scorch_parallelism(
        cpu_budget=112,
        free_cores=112,
        backlog=64,
        allowed_count=128,
        scheduler_queue_depth=60,
        local_queue_depth=20,
        profile="aggressive",
    )
    jobs_cons, threads_cons, util_cons = rescoring_scorch._allocate_scorch_parallelism(
        cpu_budget=112,
        free_cores=112,
        backlog=64,
        allowed_count=128,
        scheduler_queue_depth=60,
        local_queue_depth=20,
        profile="conservative",
    )
    assert util_aggr == 0.88
    assert util_cons == 0.80
    assert jobs_cons <= jobs_aggr
    assert threads_cons <= 4
    assert threads_aggr <= 6


def test_allocate_scorch_parallelism_invariant() -> None:
    cases = [
        (16, 16, 2, 64, 0, 0, "aggressive"),
        (64, 48, 32, 512, 15, 10, "aggressive"),
        (128, 120, 64, 2048, 80, 20, "conservative"),
        (32, 12, 8, 150, 2, 1, "conservative"),
    ]
    for cpu_budget, free_cores, backlog, allowed, sq, lq, profile in cases:
        jobs, threads, util = rescoring_scorch._allocate_scorch_parallelism(
            cpu_budget=cpu_budget,
            free_cores=free_cores,
            backlog=backlog,
            allowed_count=allowed,
            scheduler_queue_depth=sq,
            local_queue_depth=lq,
            profile=profile,
        )
        util_cores = math.floor(util * min(cpu_budget, free_cores))
        assert jobs >= 1
        assert threads >= 1
        assert jobs * threads <= max(1, util_cores)


def test_elastic_scorch_threads_respects_base_thread_cap(monkeypatch) -> None:
    monkeypatch.delenv("ATLAS_SCORCH_PARALLEL_PROFILE", raising=False)

    threads = scorch_parallel.elastic_scorch_threads(
        base_threads=1,
        allowed_count=512,
        free_cores=112,
        scheduler_queue_depth=0,
        local_queue_depth=0,
    )

    assert threads == 1


def test_normalize_parallel_profile_defaults_to_aggressive() -> None:
    assert rescoring_scorch._normalize_parallel_profile("aggressive") == "aggressive"
    assert rescoring_scorch._normalize_parallel_profile("conservative") == "conservative"
    assert rescoring_scorch._normalize_parallel_profile("unknown-value") == "aggressive"


def test_adaptive_chunk_size_creates_enough_bench2_shards(monkeypatch) -> None:
    monkeypatch.delenv("ATLAS_SCORCH_CHUNK_MIN", raising=False)
    size = scorch_parallel.adaptive_chunk_size(
        allowed_count=128,
        base_chunk_size=128,
        free_cores=112,
        scheduler_queue_depth=0,
        local_queue_depth=0,
        jobs=8,
        want_threads=14,
    )
    assert size == 16


def test_standalone_cli_jobs_respected_up_to_cpu_cap() -> None:
    jobs, threads, _util = scorch_main._resolve_standalone_cli_parallelism(
        requested_jobs=52,
        requested_threads=1,
        cpu_threads=4,
        profile="aggressive",
    )
    assert jobs == 4
    assert threads == 1

    jobs, threads, _util = scorch_main._resolve_standalone_cli_parallelism(
        requested_jobs=52,
        requested_threads=1,
        cpu_threads=56,
        profile="aggressive",
    )
    assert jobs == 52
    assert threads == 1


def test_engine_source_filter_skips_disabled_optional_engines(monkeypatch) -> None:
    monkeypatch.delenv("USE_GNINA", raising=False)
    monkeypatch.delenv("USE_LEDOCK", raising=False)
    monkeypatch.delenv("USE_DOCK6", raising=False)
    cfg = {
        "USE_GNINA": "false",
        "USE_LEDOCK": "false",
        "USE_DOCK6": "false",
    }
    assert scorch_main._engine_source_enabled(cfg, "vina") is True
    assert scorch_main._engine_source_enabled(cfg, "gnina") is False
    assert scorch_main._engine_source_enabled(cfg, "ledock") is False
    assert scorch_main._engine_source_enabled(cfg, "dock6") is False


def test_scheduler_snapshot_uses_cpu_when_no_scheduler(monkeypatch) -> None:
    monkeypatch.delenv("CPU", raising=False)
    monkeypatch.delenv("SLURM_CPUS_ON_NODE", raising=False)
    assert rescoring_scorch._scheduler_runtime_snapshot({"CPU": 112}) == (112, 0)
