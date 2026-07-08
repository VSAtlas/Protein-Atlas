from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Set


@dataclass(frozen=True)
class ChunkPartitionStats:
    chunk_count: int
    chunk_total: int
    chunk_unique: int
    duplicate_count: int
    missing_count: int
    extra_count: int
    preserved: bool


def summarize_chunk_partition(
    selected_bases: Set[str],
    chunks: Sequence[Set[str]],
) -> ChunkPartitionStats:
    chunk_total = 0
    chunk_union: Set[str] = set()
    for chunk in chunks:
        chunk_total += len(chunk)
        chunk_union.update(chunk)
    chunk_unique = len(chunk_union)
    duplicate_count = max(0, int(chunk_total) - int(chunk_unique))
    missing_count = len(set(selected_bases) - chunk_union)
    extra_count = len(chunk_union - set(selected_bases))
    preserved = bool(missing_count == 0 and extra_count == 0)
    return ChunkPartitionStats(
        chunk_count=int(len(chunks)),
        chunk_total=int(chunk_total),
        chunk_unique=int(chunk_unique),
        duplicate_count=int(duplicate_count),
        missing_count=int(missing_count),
        extra_count=int(extra_count),
        preserved=preserved,
    )
