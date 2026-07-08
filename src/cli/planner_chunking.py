from __future__ import annotations

import math
import os
import socket
import logging
from typing import Any, Mapping

from cli.distributed_context import detect_allocated_cpus

_CHUNK_TARGET_SIZE = 96
_CHUNK_MIN_SIZE = 8
_CHUNK_PRODUCTION_MIN_SIZE = 32
_CHUNK_MAX_SIZE = 256
_CHUNK_OVERSHARD_FACTOR = 6
_CHUNK_LOCAL_WORKER_CPU_DIVISOR = 1
_CHUNK_LOCAL_WORKER_MAX = 512


def _env_int(name: str, default: int, *, min_value: int, max_value: int) -> int:
    raw = os.environ.get(name)
    if raw is None or not str(raw).strip():
        return int(default)
    try:
        parsed = int(str(raw).strip())
    except Exception:
        return int(default)
    return max(int(min_value), min(int(max_value), int(parsed)))


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return bool(default)
    token = str(raw).strip().lower()
    if token in {"1", "true", "yes", "on"}:
        return True
    if token in {"0", "false", "no", "off"}:
        return False
    return bool(default)


def _adaptive_chunk_min_size(total_ligands: int, workers_total: int) -> int:
    explicit = os.environ.get("ATLAS_CHUNK_MIN_SIZE")
    if explicit is not None and str(explicit).strip():
        return _env_int(
            "ATLAS_CHUNK_MIN_SIZE",
            _CHUNK_MIN_SIZE,
            min_value=1,
            max_value=_CHUNK_MAX_SIZE,
        )
    if not _env_bool("ATLAS_ADAPTIVE_CHUNK_MIN_SIZE", True):
        return int(_CHUNK_MIN_SIZE)
    total = max(1, int(total_ligands))
    workers = max(1, int(workers_total))
    ligands_per_worker = float(total) / float(workers)
    if workers >= 128 and ligands_per_worker <= 8.0:
        return 4
    if workers >= 128 and ligands_per_worker <= 16.0:
        return 6
    return int(_CHUNK_MIN_SIZE)


def _planner_chunk_min_size(
    total_ligands: int,
    workers_total: int,
    *,
    production: bool = False,
) -> int:
    explicit = os.environ.get("ATLAS_CHUNK_MIN_SIZE")
    if explicit is not None and str(explicit).strip():
        return _adaptive_chunk_min_size(total_ligands, workers_total)
    if production:
        return int(_CHUNK_PRODUCTION_MIN_SIZE)
    return _adaptive_chunk_min_size(total_ligands, workers_total)


def _cfg_cpu_total(cfg: Mapping[str, Any]) -> int:
    try:
        return max(1, int(cfg.get("CPU", os.cpu_count() or 1) or 1))
    except Exception:
        return max(1, int(os.cpu_count() or 1))


def _distributed_chunk_worker_plan(cpu_total: int) -> tuple[int, int, int]:
    cpu_total = max(1, int(cpu_total))
    workers = max(1, min(int(_CHUNK_LOCAL_WORKER_MAX), cpu_total))
    cpu_per_worker = 1
    return workers, cpu_per_worker, workers * cpu_per_worker


def _split_ligands_into_chunks(
    ligand_bases: list[str],
    *,
    task_count: int,
    local_workers_per_task: int = 1,
    ligand_weights: Mapping[str, float] | None = None,
    min_chunk_size: int | None = None,
    scale_with_workers: bool = True,
) -> list[list[str]]:
    total = len(ligand_bases)
    if total <= 0:
        return [[]]

    ordered = sorted(ligand_bases)
    workers_total = max(1, int(task_count) * max(1, int(local_workers_per_task)))
    if min_chunk_size is None:
        min_chunk_size = _adaptive_chunk_min_size(total, workers_total)
    min_chunk_size = max(1, min(int(min_chunk_size), int(_CHUNK_MAX_SIZE)))
    base_chunks = max(1, int(math.ceil(total / float(_CHUNK_TARGET_SIZE))))
    max_chunks_by_min = max(1, int(total // max(1, min_chunk_size)))
    if scale_with_workers:
        desired_by_workers = max(
            1,
            int(workers_total) * int(_CHUNK_OVERSHARD_FACTOR),
        )
        desired_chunks = min(max_chunks_by_min, max(base_chunks, desired_by_workers))
    else:
        desired_chunks = min(max_chunks_by_min, base_chunks)
    desired_chunks = max(1, desired_chunks)
    chunk_size = int(math.ceil(total / float(desired_chunks)))
    chunk_size = max(1, min(_CHUNK_MAX_SIZE, chunk_size))
    if ligand_weights:
        bins: list[list[str]] = [[] for _ in range(desired_chunks)]
        loads = [0.0 for _ in range(desired_chunks)]
        counts = [0 for _ in range(desired_chunks)]

        def _weight(token: str) -> float:
            try:
                return max(0.01, float(ligand_weights.get(str(token), 1.0)))
            except Exception:
                return 1.0

        for ligand in sorted(ordered, key=lambda item: (-_weight(item), item)):
            target_idx = min(
                range(desired_chunks),
                key=lambda idx: (
                    counts[idx] >= chunk_size,
                    loads[idx],
                    counts[idx],
                    idx,
                ),
            )
            bins[target_idx].append(ligand)
            loads[target_idx] += _weight(ligand)
            counts[target_idx] += 1
        return [sorted(chunk) for chunk in bins if chunk]

    chunks: list[list[str]] = []
    for i in range(0, total, chunk_size):
        chunks.append(ordered[i : i + chunk_size])
    return chunks or [[]]


def _distributed_chunk_local_workers(cfg: Mapping[str, Any]) -> int:
    workers, _, _ = _distributed_chunk_worker_plan(_cfg_cpu_total(cfg))
    return workers


def _distributed_chunk_cpu_per_worker(
    cfg: Mapping[str, Any], *, local_workers: int
) -> int:
    workers = max(1, int(local_workers))
    cpu_total = _cfg_cpu_total(cfg)
    return max(1, int(cpu_total // workers))


def _resolve_runtime_cpu_settings(
    cfg: Mapping[str, Any], *, distributed_enabled: bool
) -> tuple[int, int, int, int]:
    try:
        cfg_cpu_raw = int(cfg.get("CPU", os.cpu_count() or 1) or 1)
    except Exception:
        cfg_cpu_raw = int(os.cpu_count() or 1)
    cfg_cpu_raw = max(1, cfg_cpu_raw)
    alloc_cpu = max(1, int(detect_allocated_cpus(cfg_cpu_raw)))

    if distributed_enabled:
        effective_cpu = alloc_cpu
        effective_scheduler_cpus = alloc_cpu
        return cfg_cpu_raw, alloc_cpu, effective_cpu, effective_scheduler_cpus

    effective_cpu = max(1, min(cfg_cpu_raw, alloc_cpu))
    try:
        sched_cpu_raw = int(
            cfg.get("GLOBAL_SCHEDULER_CPUS", effective_cpu) or effective_cpu
        )
    except Exception:
        sched_cpu_raw = effective_cpu
    effective_scheduler_cpus = max(1, min(sched_cpu_raw, effective_cpu))
    return cfg_cpu_raw, alloc_cpu, effective_cpu, effective_scheduler_cpus


def _log_distributed_cpu_node(
    *,
    run_id: str,
    dist_ctx: Any,
    cfg_cpu_raw: int,
    alloc_cpu: int,
    effective_cpu: int,
    effective_scheduler_cpus: int,
) -> None:
    logging.info(
        "[distributed.cpu.node] run_id=%s task_id=%d task_count=%d task_min_id=%d leader_task_id=%d hostname=%s "
        "cfg_cpu=%d alloc_cpu=%d effective_cpu=%d effective_scheduler_cpus=%d "
        "slurm_cpus_per_task=%s slurm_cpus_on_node=%s slurm_job_cpus_per_node=%s",
        str(run_id),
        int(getattr(dist_ctx, "task_id", 0)),
        int(getattr(dist_ctx, "task_count", 1)),
        int(getattr(dist_ctx, "task_min_id", 0)),
        int(getattr(dist_ctx, "leader_task_id", 0)),
        socket.gethostname(),
        int(cfg_cpu_raw),
        int(alloc_cpu),
        int(effective_cpu),
        int(effective_scheduler_cpus),
        str(os.environ.get("SLURM_CPUS_PER_TASK", "")),
        str(os.environ.get("SLURM_CPUS_ON_NODE", "")),
        str(os.environ.get("SLURM_JOB_CPUS_PER_NODE", "")),
    )
