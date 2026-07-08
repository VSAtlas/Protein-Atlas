from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

from post_docking.rescoring.scorch_types import ScorchTask, SelectionResult

ComboKey = Tuple[str, str, str]


@dataclass(frozen=True)
class OrchestrationLimits:
    chunk_size: int
    chunk_min: int
    tail_split_trigger: int
    runtime_rechunk_stale_sec: float
    hedge_timeout_sec: float
    hedge_max_inflight: int
    hedge_min_allowed: int


@dataclass(frozen=True)
class OrchestrationDeps:
    scheduler_runtime_snapshot: Callable[..., Tuple[int, int]]
    discover_mode_dirs: Callable[..., Dict[str, List[Path]]]
    load_control_bases: Callable[..., Set[str]]
    stage_dir_candidates: Callable[..., Sequence[str]]
    score_csv_for_spec: Callable[..., Tuple[Optional[Path], Sequence[str], bool]]
    select_top_bases_from_score_csv: Callable[..., SelectionResult]
    adaptive_chunk_size: Callable[..., int]
    chunk_allowed_bases: Callable[..., List[Set[str]]]
    tail_split_allowed_chunks: Callable[..., List[Set[str]]]
    score_stage: Callable[..., Tuple[bool, Optional[Path]]]
    task_estimate_seconds: Callable[..., float]
    task_allowed_count: Callable[..., int]
    split_task_for_rechunk: Callable[..., List[ScorchTask]]
    elastic_scorch_threads: Callable[..., int]
    aggregate_combo: Callable[..., Optional[Path]]
    annotate_scorch_z_scores: Callable[[Path, Path, logging.Logger], None]
    emit_task_event: Callable[..., None]
    emit_bench_event: Callable[..., None]
    summarize_chunk_partition: Callable[..., Any]
    find_consensus_csv: Optional[Callable[[Path], Optional[Path]]] = None
    rerank_consensus_with_scorch: Optional[Callable[..., bool]] = None


@dataclass
class OrchestrationResult:
    failed_jobs: int
    combos_with_tasks: Set[ComboKey]
    combo_failed_local: Dict[ComboKey, bool]
    local_completed_combos: Set[ComboKey]


@dataclass(frozen=True)
class DoneSentinelDeps:
    emit_task_event: Callable[..., None]
    mark_combo_done: Callable[..., None]


def receptor_path_for_combo(processed_root: Path, combo: ComboKey) -> Path:
    pdb_id, variant, ph = combo
    receptor_root = processed_root / pdb_id / variant / "receptor"
    if str(ph).strip():
        return receptor_root / "ph_ensemble" / f"{pdb_id}_{ph}.pdbqt"
    return receptor_root / f"{pdb_id}.pdbqt"
