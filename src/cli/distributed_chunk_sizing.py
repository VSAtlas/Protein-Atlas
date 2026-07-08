from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class ChunkSizingPolicy:
    target_size: int
    min_size: int
    max_size: int
    parallelism_floor_mult: float


def adaptive_chunk_sizing_enabled(*, distributed_enabled: bool) -> bool:
    return bool(distributed_enabled)


def resolve_chunk_sizing_policy(
    *,
    ligand_count: int,
    distributed_enabled: bool,
) -> ChunkSizingPolicy:
    total = max(0, int(ligand_count))
    if not adaptive_chunk_sizing_enabled(distributed_enabled=distributed_enabled):
        return ChunkSizingPolicy(
            target_size=96,
            min_size=32,
            max_size=256,
            parallelism_floor_mult=6.0,
        )

    # Production-focused sizing for large runs while preserving enough parallelism.
    if total >= 20000:
        return ChunkSizingPolicy(
            target_size=512,
            min_size=192,
            max_size=1024,
            parallelism_floor_mult=1.0,
        )
    if total >= 5000:
        return ChunkSizingPolicy(
            target_size=320,
            min_size=128,
            max_size=640,
            parallelism_floor_mult=1.25,
        )
    if total >= 1500:
        return ChunkSizingPolicy(
            target_size=192,
            min_size=64,
            max_size=384,
            parallelism_floor_mult=1.75,
        )
    return ChunkSizingPolicy(
        target_size=128,
        min_size=64,
        max_size=320,
        parallelism_floor_mult=1.25,
    )


def split_ligands_into_chunks_adaptive(
    ligand_bases: Sequence[str],
    *,
    task_count: int,
    local_workers_per_task: int,
    distributed_enabled: bool,
) -> list[list[str]]:
    total = len(ligand_bases)
    if total <= 0:
        return [[]]

    policy = resolve_chunk_sizing_policy(
        ligand_count=total,
        distributed_enabled=distributed_enabled,
    )

    ordered = sorted(str(x) for x in ligand_bases if str(x).strip())
    if not ordered:
        return [[]]

    workers_total = max(1, int(task_count)) * max(1, int(local_workers_per_task))
    target_chunks_by_size = max(1, int(math.ceil(total / float(policy.target_size))))
    parallel_floor = max(
        1,
        int(math.ceil(float(workers_total) * float(policy.parallelism_floor_mult))),
    )
    max_chunks_by_min = max(1, int(total // max(1, int(policy.min_size))))
    desired_chunks = max(target_chunks_by_size, min(parallel_floor, max_chunks_by_min))
    desired_chunks = max(1, min(desired_chunks, max_chunks_by_min))

    chunk_size = int(math.ceil(float(total) / float(desired_chunks)))
    chunk_size = max(1, min(int(policy.max_size), chunk_size))

    chunks: list[list[str]] = []
    for i in range(0, total, chunk_size):
        chunks.append(ordered[i : i + chunk_size])
    return chunks or [[]]
