from __future__ import annotations

import logging

from cli.distributed_context import detect_allocated_cpus
from cli import planner_chunking


def test_distributed_mode_uses_allocated_cpus(monkeypatch) -> None:
    cfg = {"CPU": 4, "GLOBAL_SCHEDULER_CPUS": 2}
    monkeypatch.setattr(planner_chunking, "detect_allocated_cpus", lambda _fallback: 16)

    cfg_cpu, alloc_cpu, effective_cpu, effective_scheduler = (
        planner_chunking._resolve_runtime_cpu_settings(cfg, distributed_enabled=True)
    )

    assert cfg_cpu == 4
    assert alloc_cpu == 16
    assert effective_cpu == 16
    assert effective_scheduler == 16


def test_non_distributed_mode_keeps_config_cap(monkeypatch) -> None:
    cfg = {"CPU": 4, "GLOBAL_SCHEDULER_CPUS": 12}
    monkeypatch.setattr(planner_chunking, "detect_allocated_cpus", lambda _fallback: 32)

    cfg_cpu, alloc_cpu, effective_cpu, effective_scheduler = (
        planner_chunking._resolve_runtime_cpu_settings(cfg, distributed_enabled=False)
    )

    assert cfg_cpu == 4
    assert alloc_cpu == 32
    assert effective_cpu == 4
    assert effective_scheduler == 4


def test_non_distributed_mode_caps_by_alloc_and_scheduler(monkeypatch) -> None:
    cfg = {"CPU": 32, "GLOBAL_SCHEDULER_CPUS": 6}
    monkeypatch.setattr(planner_chunking, "detect_allocated_cpus", lambda _fallback: 8)

    cfg_cpu, alloc_cpu, effective_cpu, effective_scheduler = (
        planner_chunking._resolve_runtime_cpu_settings(cfg, distributed_enabled=False)
    )

    assert cfg_cpu == 32
    assert alloc_cpu == 8
    assert effective_cpu == 8
    assert effective_scheduler == 6


def test_distributed_cpu_node_log_includes_slurm_cpu_env(monkeypatch, caplog) -> None:
    class _Ctx:
        task_id = 3
        task_count = 8
        task_min_id = 1
        leader_task_id = 1

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "64")
    monkeypatch.setenv("SLURM_JOB_CPUS_PER_NODE", "64(x2)")
    monkeypatch.setattr(planner_chunking.socket, "gethostname", lambda: "node-a")

    caplog.set_level(logging.INFO)
    planner_chunking._log_distributed_cpu_node(
        run_id="rid_cpu",
        dist_ctx=_Ctx(),
        cfg_cpu_raw=4,
        alloc_cpu=16,
        effective_cpu=16,
        effective_scheduler_cpus=16,
    )
    text = caplog.text
    assert "[distributed.cpu.node]" in text
    assert "run_id=rid_cpu" in text
    assert "task_id=3" in text
    assert "task_count=8" in text
    assert "hostname=node-a" in text
    assert "effective_cpu=16" in text
    assert "slurm_cpus_per_task=16" in text
    assert "slurm_cpus_on_node=64" in text
    assert "slurm_job_cpus_per_node=64(x2)" in text


def test_detect_allocated_cpus_can_prefer_whole_node(monkeypatch) -> None:
    monkeypatch.setenv("ATLAS_USE_NODE_CPUS", "1")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "112")
    monkeypatch.setenv("SLURM_JOB_CPUS_PER_NODE", "112")

    assert detect_allocated_cpus(4) == 112


def test_detect_allocated_cpus_uses_job_node_count_when_srun_reports_one(
    monkeypatch,
) -> None:
    monkeypatch.setenv("ATLAS_USE_NODE_CPUS", "1")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "1")
    monkeypatch.setenv("SLURM_JOB_CPUS_PER_NODE", "48(x3)")

    assert detect_allocated_cpus(4) == 48


def test_detect_allocated_cpus_prefers_task_without_whole_node(monkeypatch) -> None:
    monkeypatch.delenv("ATLAS_USE_NODE_CPUS", raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "112")

    assert detect_allocated_cpus(4) == 8
