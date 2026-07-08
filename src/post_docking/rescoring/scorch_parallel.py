from __future__ import annotations

import math
import os
from typing import Any, List, Set

from post_docking.rescoring.scorch_types import ScorchTask

SCORCH_CHUNK_MIN = 16
SCORCH_TAIL_SPLIT_TRIGGER = 64
SCORCH_HEAVY_ALLOWED_THRESHOLD = 512
SCORCH_HEAVY_CHUNK_CAP = 96
SCORCH_PARALLEL_PROFILE_DEFAULT = "aggressive"
SCORCH_PARALLEL_PROFILE_ENV = "ATLAS_SCORCH_PARALLEL_PROFILE"


def _chunk_min_from_env() -> int:
    for name in ("ATLAS_SCORCH_CHUNK_MIN", "SCORCH_CHUNK_MIN"):
        raw = os.environ.get(name)
        if raw is None:
            continue
        try:
            value = int(str(raw).strip())
        except Exception:
            continue
        if value > 0:
            return int(value)
    return int(SCORCH_CHUNK_MIN)


def normalize_parallel_profile(raw: Any) -> str:
    token = str(raw or "").strip().lower()
    if token in {"conservative", "safe", "balanced"}:
        return "conservative"
    return SCORCH_PARALLEL_PROFILE_DEFAULT


def util_target_for_profile(
    *,
    profile: str,
    free_cores: int,
    pressure: int,
    backlog: int,
) -> float:
    profile_norm = normalize_parallel_profile(profile)
    free = max(1, int(free_cores))
    queue_pressure = max(0, int(pressure))
    pending = max(1, int(backlog))
    if profile_norm == "conservative":
        u_low, u_mid, u_high = 0.80, 0.85, 0.90
        pressure_trigger = max(2, free // 2)
    else:
        u_low, u_mid, u_high = 0.88, 0.93, 0.97
        pressure_trigger = max(2, free // 3)
    if queue_pressure >= pressure_trigger:
        return u_low
    if pending <= max(2, free // 4):
        return u_high
    return u_mid


def allocate_scorch_parallelism(
    *,
    cpu_budget: int,
    free_cores: int,
    backlog: int,
    allowed_count: int,
    scheduler_queue_depth: int = 0,
    local_queue_depth: int = 0,
    profile: str = SCORCH_PARALLEL_PROFILE_DEFAULT,
) -> tuple[int, int, float]:
    budget = max(1, int(cpu_budget))
    free = max(1, min(int(free_cores), budget))
    pending = max(1, int(backlog))
    allowed_n = max(1, int(allowed_count))
    pressure = max(int(scheduler_queue_depth), int(local_queue_depth))
    profile_norm = normalize_parallel_profile(profile)
    util_target = util_target_for_profile(
        profile=profile_norm, free_cores=free, pressure=pressure, backlog=pending
    )
    if allowed_n <= 48:
        demand_from_allowed = 1
    elif allowed_n <= 128:
        demand_from_allowed = 2
    elif allowed_n <= 512:
        demand_from_allowed = int(math.ceil(float(allowed_n) / 64.0))
    else:
        demand_from_allowed = int(math.ceil(float(allowed_n) / 96.0))
    pending_effective = max(pending, demand_from_allowed)

    if allowed_n < 96:
        base_threads = 1
    elif allowed_n < 384:
        base_threads = 2
    elif allowed_n < 1024:
        base_threads = 3
    else:
        base_threads = 4

    bonus = 0
    if free >= 32:
        bonus += 1
    if free >= 96:
        bonus += 1

    if profile_norm == "conservative":
        threads_max = 4
        pressure_trigger = max(2, free // 2)
    else:
        threads_max = 6
        pressure_trigger = max(2, free // 3)
    penalty = 1 if pressure >= pressure_trigger else 0

    threads = base_threads + bonus - penalty
    threads = max(1, min(int(threads_max), int(threads)))

    util_cores = max(1, int(math.floor(float(util_target) * float(free))))
    jobs = max(1, min(pending_effective, util_cores // max(1, threads)))

    if profile_norm == "aggressive" and free >= 32 and pending <= 3 and allowed_n >= 64:
        floor_cores = 12
        threads_hard_cap = max(threads_max, 6)
        if free >= 64:
            floor_cores = 18
            threads_hard_cap = max(threads_hard_cap, 8)
        if free >= 96:
            floor_cores = 24
            threads_hard_cap = max(threads_hard_cap, 10)
        floor_cores = max(1, min(util_cores, int(floor_cores)))
        used_cores = max(1, int(jobs) * int(threads))
        if used_cores < floor_cores:
            min_threads = int(math.ceil(float(floor_cores) / float(max(1, jobs))))
            threads = max(threads, min_threads)
            threads = max(1, min(int(threads_hard_cap), int(threads), int(util_cores)))
            max_jobs = max(1, int(util_cores // max(1, threads)))
            min_jobs = int(math.ceil(float(floor_cores) / float(max(1, threads))))
            jobs = max(jobs, min(max_jobs, max(min_jobs, pending_effective)))

    while jobs > 1 and (jobs * threads) > util_cores:
        jobs -= 1
    if jobs * threads > util_cores:
        threads = max(1, min(threads, util_cores))
    return int(jobs), int(threads), float(util_target)


def adaptive_chunk_size(
    *,
    allowed_count: int,
    base_chunk_size: int,
    free_cores: int,
    scheduler_queue_depth: int,
    local_queue_depth: int,
    jobs: int,
    want_threads: int = 1,
) -> int:
    total = max(0, int(allowed_count))
    if total <= 0:
        return 0
    chunk_min = max(1, _chunk_min_from_env())
    if total <= chunk_min:
        return total

    size = max(chunk_min, int(base_chunk_size))
    if total >= SCORCH_HEAVY_ALLOWED_THRESHOLD:
        size = min(size, SCORCH_HEAVY_CHUNK_CAP)
    if int(want_threads) >= 5:
        size = min(size, 64)

    if free_cores > 0:
        if free_cores <= 2:
            size = min(size, 32)
        elif free_cores <= 4:
            size = min(size, 48)
        elif free_cores <= 8:
            size = min(size, 64)

    pressure = max(int(scheduler_queue_depth), int(local_queue_depth))
    if pressure <= max(2, int(jobs)):
        size = min(size, 64 if total < SCORCH_HEAVY_ALLOWED_THRESHOLD else 80)

    if free_cores >= 16 and pressure <= max(2, int(jobs)):
        max_productive_shards = max(1, int(math.ceil(float(total) / float(chunk_min))))
        core_target_shards = max(
            1,
            int(math.ceil(float(max(1, int(free_cores))) / float(max(1, int(want_threads))))),
        )
        desired_shards = min(
            max_productive_shards,
            max(max(1, int(jobs) * 2), core_target_shards),
        )
        if desired_shards > 1:
            size = min(
                size,
                max(chunk_min, int(math.ceil(float(total) / float(desired_shards)))),
            )

    if total <= size * 2:
        size = min(size, max(chunk_min, total // 2))

    return max(1, min(size, total))


def tail_split_allowed_chunks(
    chunks: List[Set[str]],
    *,
    tail_mode: bool,
    min_split: int,
) -> List[Set[str]]:
    if not tail_mode or len(chunks) > 2:
        return chunks

    split_threshold = max(2, int(min_split) * 2)
    out: List[Set[str]] = []
    for chunk in chunks:
        if len(chunk) < split_threshold:
            out.append(chunk)
            continue
        ordered = sorted(chunk)
        half = max(1, len(ordered) // 2)
        left = set(ordered[:half])
        right = set(ordered[half:])
        if left:
            out.append(left)
        if right:
            out.append(right)
    return out or chunks


def chunk_allowed_bases(bases: Set[str], chunk_size: int) -> List[Set[str]]:
    ordered = sorted(set(bases))
    if not ordered:
        return []
    size = max(1, int(chunk_size))
    chunks: List[Set[str]] = []
    for i in range(0, len(ordered), size):
        chunks.append(set(ordered[i : i + size]))
    return chunks


def task_allowed_count(task: ScorchTask) -> int:
    return len(task[3])


def task_estimate_seconds(task: ScorchTask, thread_budget: int) -> float:
    allowed_n = max(1, task_allowed_count(task))
    threads = max(1, int(thread_budget))
    return max(0.5, (allowed_n * 0.08) / float(threads))


def elastic_scorch_threads(
    *,
    base_threads: int,
    allowed_count: int,
    free_cores: int,
    scheduler_queue_depth: int,
    local_queue_depth: int,
) -> int:
    base = max(1, int(base_threads))
    free = max(1, int(free_cores))
    profile = normalize_parallel_profile(
        os.environ.get(SCORCH_PARALLEL_PROFILE_ENV, SCORCH_PARALLEL_PROFILE_DEFAULT)
    )
    _jobs_alloc, threads_alloc, _util = allocate_scorch_parallelism(
        cpu_budget=max(base, free),
        free_cores=free,
        backlog=max(1, int(local_queue_depth) + 1),
        allowed_count=max(1, int(allowed_count)),
        scheduler_queue_depth=max(0, int(scheduler_queue_depth)),
        local_queue_depth=max(0, int(local_queue_depth)),
        profile=profile,
    )
    return max(1, min(base, int(threads_alloc)))


def split_task_for_rechunk(task: ScorchTask, suffix_seed: str) -> List[ScorchTask]:
    (
        spec,
        combo,
        receptor,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
    ) = task
    if len(allowed_bases) < SCORCH_TAIL_SPLIT_TRIGGER:
        return [task]

    allowed_n = len(allowed_bases)
    if allowed_n >= SCORCH_HEAVY_ALLOWED_THRESHOLD:
        parts = 4
    elif allowed_n >= (SCORCH_TAIL_SPLIT_TRIGGER * 2):
        parts = 3
    else:
        parts = 2

    split_size = max(SCORCH_CHUNK_MIN, math.ceil(allowed_n / float(parts)))
    chunks = chunk_allowed_bases(allowed_bases, split_size)
    chunks = tail_split_allowed_chunks(
        chunks,
        tail_mode=True,
        min_split=max(SCORCH_CHUNK_MIN, split_size // 3),
    )
    if len(chunks) <= 1:
        return [task]

    root = chunk_tag or f"rt{suffix_seed}"
    out: List[ScorchTask] = []
    for idx, bases in enumerate(chunks):
        out.append(
            (
                spec,
                combo,
                receptor,
                bases,
                control_bases,
                run_mode,
                stage_dirs,
                score_csv,
                {
                    base: stage
                    for base, stage in selected_stage_by_base.items()
                    if base in bases
                },
                {
                    base: score
                    for base, score in selected_score_by_base.items()
                    if base in bases
                },
                f"{root}r{idx:02d}",
            )
        )
    return out
