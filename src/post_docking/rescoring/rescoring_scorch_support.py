from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from post_docking.rescoring import scorch_parallel as _scorch_parallel_mod
from post_docking.rescoring import scorch_prefix_and_stages as _scorch_prefix_mod
from post_docking.rescoring import scorch_selection as _scorch_selection_mod
from post_docking.rescoring import scorch_stage_io as _scorch_stage_io_mod
from post_docking.rescoring.scorch_types import ScorchTask, StageSpec

DECOY_PREFIX_KEY = "DECOY_PREFIX"
DUD_PREFIX_KEY = "DUD_PREFIX"
DECOY_PREFIX_DEFAULT = "dud"
DECOY_PREFIX_VALUE = DECOY_PREFIX_DEFAULT

SCORCH_CHUNK_SIZE = 128
SCORCH_CHUNK_MIN = 16
SCORCH_TAIL_SPLIT_TRIGGER = 64
SCORCH_RUNTIME_RECHUNK_STALE_SEC = 3.0
SCORCH_HEDGE_TIMEOUT_SEC = 30.0
SCORCH_HEDGE_MAX_INFLIGHT = 1
SCORCH_HEDGE_MIN_ALLOWED = 64
SCORCH_PARALLEL_PROFILE_DEFAULT = "aggressive"
SCORCH_PARALLEL_PROFILE_ENV = "ATLAS_SCORCH_PARALLEL_PROFILE"


def _normalize_parallel_profile(raw: Any) -> str:
    return _scorch_parallel_mod.normalize_parallel_profile(raw)


def _util_target_for_profile(
    *,
    profile: str,
    free_cores: int,
    pressure: int,
    backlog: int,
) -> float:
    return _scorch_parallel_mod.util_target_for_profile(
        profile=profile,
        free_cores=free_cores,
        pressure=pressure,
        backlog=backlog,
    )


def _allocate_scorch_parallelism(
    *,
    cpu_budget: int,
    free_cores: int,
    backlog: int,
    allowed_count: int,
    scheduler_queue_depth: int = 0,
    local_queue_depth: int = 0,
    profile: str = SCORCH_PARALLEL_PROFILE_DEFAULT,
) -> tuple[int, int, float]:
    return _scorch_parallel_mod.allocate_scorch_parallelism(
        cpu_budget=cpu_budget,
        free_cores=free_cores,
        backlog=backlog,
        allowed_count=allowed_count,
        scheduler_queue_depth=scheduler_queue_depth,
        local_queue_depth=local_queue_depth,
        profile=profile,
    )


def _decoy_prefix_value() -> str:
    return DECOY_PREFIX_VALUE or DECOY_PREFIX_DEFAULT


def _normalize_decoy_prefix(value: Optional[object]) -> str:
    return _scorch_prefix_mod.normalize_decoy_prefix(
        value,
        default=DECOY_PREFIX_DEFAULT,
    )


def _resolve_decoy_prefix_from_config(cfg: Dict[str, object]) -> str:
    return _scorch_prefix_mod.resolve_decoy_prefix_from_config(
        cfg,
        decoy_prefix_key=DECOY_PREFIX_KEY,
        dud_prefix_key=DUD_PREFIX_KEY,
        default_prefix=DECOY_PREFIX_DEFAULT,
    )


def _set_decoy_prefix(
    value: Optional[object], logger: Optional[logging.Logger] = None
) -> str:
    global DECOY_PREFIX_VALUE
    DECOY_PREFIX_VALUE = _normalize_decoy_prefix(value)
    _scorch_stage_io_mod.set_decoy_prefix(DECOY_PREFIX_VALUE)
    if logger:
        logger.info("[scorch.preflight] decoy_prefix=%s", DECOY_PREFIX_VALUE)
    return DECOY_PREFIX_VALUE


def _parse_test_mode_value(raw: object) -> list[str]:
    return _scorch_prefix_mod.parse_test_mode_value(raw)


def _parse_test_mode_tokens(cfg: Dict[str, object]) -> list[str]:
    return _scorch_prefix_mod.parse_test_mode_tokens(cfg, prefer_env=True)


def _resolve_decoy_prefix_override(args: argparse.Namespace) -> Optional[str]:
    return _scorch_prefix_mod.resolve_decoy_prefix_override_with_env(
        args,
        env_keys=(DECOY_PREFIX_KEY, DUD_PREFIX_KEY),
        default_prefix=DECOY_PREFIX_DEFAULT,
    )


def _decoy_prefixes_from_test_mode(
    cfg: Dict[str, object], override: Optional[str]
) -> list[str]:
    return _scorch_prefix_mod.decoy_prefixes_from_test_mode(
        cfg,
        override,
        default_prefix=DECOY_PREFIX_DEFAULT,
        prefer_env=True,
    )


def _decoy_stage_dirs(order: Sequence[int]) -> Tuple[str, ...]:
    return _scorch_prefix_mod.decoy_stage_dirs(
        order,
        decoy_prefix=_decoy_prefix_value(),
    )


def _decoy_engine_stage_dirs(engine: str, order: Sequence[int]) -> Tuple[str, ...]:
    return _scorch_prefix_mod.decoy_engine_stage_dirs(
        engine,
        order,
        decoy_prefix=_decoy_prefix_value(),
    )


def _decoy_engine_stage_dirs_legacy(
    engine: str, order: Sequence[int]
) -> Tuple[str, ...]:
    return _scorch_prefix_mod.decoy_engine_stage_dirs_legacy(
        engine,
        order,
        decoy_prefix=_decoy_prefix_value(),
    )


def _decoy_engine_stage_dirs_all(engine: str) -> Tuple[str, ...]:
    return _scorch_prefix_mod.decoy_engine_stage_dirs_all(
        engine,
        decoy_prefix=_decoy_prefix_value(),
    )


def _decoy_post_stage_dirs() -> Tuple[str, ...]:
    return _scorch_prefix_mod.decoy_post_stage_dirs(
        decoy_prefix=_decoy_prefix_value(),
    )


def _vina_stage_dirs() -> Tuple[str, ...]:
    return _scorch_prefix_mod.vina_stage_dirs(decoy_prefix=_decoy_prefix_value())


def _gnina_stage_dirs() -> Tuple[str, ...]:
    return _scorch_prefix_mod.gnina_stage_dirs(decoy_prefix=_decoy_prefix_value())


def _dock6_stage_dirs() -> Tuple[str, ...]:
    return _scorch_prefix_mod.dock6_stage_dirs(decoy_prefix=_decoy_prefix_value())


def _ledock_stage_dirs() -> Tuple[str, ...]:
    return _scorch_prefix_mod.ledock_stage_dirs(decoy_prefix=_decoy_prefix_value())


def _post_stage_dirs() -> Tuple[str, ...]:
    return _scorch_prefix_mod.post_stage_dirs(decoy_prefix=_decoy_prefix_value())


def discover_stage3_roots(variant_root: Path) -> Dict[str, Path]:
    return _scorch_prefix_mod.discover_stage3_roots(
        variant_root,
        decoy_prefix=_decoy_prefix_value(),
    )


def _discover_mode_dirs(
    combo: Tuple[str, str, str],
    run_root: Path,
    post_root: Path,
    specs: Sequence[StageSpec],
    logger: logging.Logger,
) -> Dict[str, List[Path]]:
    return _scorch_selection_mod.discover_mode_dirs(
        combo,
        run_root,
        post_root,
        specs,
        logger,
        decoy_prefix=_decoy_prefix_value(),
    )


def _prefer_existing(
    root: Optional[Path], primary: Sequence[str], legacy: Sequence[str]
) -> Tuple[str, ...]:
    return _scorch_prefix_mod._prefer_existing(root, primary, legacy)


def stage_dir_candidates(
    source: str, mode: str, root: Optional[Path] = None
) -> Tuple[str, ...]:
    return _scorch_prefix_mod.stage_dir_candidates(
        source,
        mode,
        root,
        decoy_prefix=_decoy_prefix_value(),
    )


def _scheduler_runtime_snapshot(cfg: Dict[str, object]) -> Tuple[int, int]:
    scheduler = cfg.get("GLOBAL_DOCKING_SCHEDULER") if isinstance(cfg, dict) else None
    free_cores = 0
    queued = 0
    if scheduler is None:
        for raw in (
            cfg.get("CPU") if isinstance(cfg, dict) else None,
            os.environ.get("CPU"),
            os.environ.get("SLURM_CPUS_ON_NODE"),
            os.environ.get("SLURM_CPUS_PER_TASK"),
        ):
            try:
                free_cores = max(1, int(str(raw).strip()))
                break
            except Exception:
                continue
        return free_cores, queued
    if hasattr(scheduler, "available_cores"):
        try:
            free_cores = max(0, int(scheduler.available_cores()))
        except Exception:
            free_cores = 0
    if hasattr(scheduler, "queued_requests"):
        try:
            queued = max(0, int(scheduler.queued_requests()))
        except Exception:
            queued = 0
    return free_cores, queued


def _adaptive_chunk_size(
    *,
    allowed_count: int,
    base_chunk_size: int,
    free_cores: int,
    scheduler_queue_depth: int,
    local_queue_depth: int,
    jobs: int,
    want_threads: int = 1,
) -> int:
    return _scorch_parallel_mod.adaptive_chunk_size(
        allowed_count=allowed_count,
        base_chunk_size=base_chunk_size,
        free_cores=free_cores,
        scheduler_queue_depth=scheduler_queue_depth,
        local_queue_depth=local_queue_depth,
        jobs=jobs,
        want_threads=want_threads,
    )


def _tail_split_allowed_chunks(
    chunks: List[Set[str]],
    *,
    tail_mode: bool,
    min_split: int,
) -> List[Set[str]]:
    return _scorch_parallel_mod.tail_split_allowed_chunks(
        chunks,
        tail_mode=tail_mode,
        min_split=min_split,
    )


def _chunk_allowed_bases(bases: Set[str], chunk_size: int) -> List[Set[str]]:
    return _scorch_parallel_mod.chunk_allowed_bases(bases, chunk_size)


def _task_allowed_count(task: ScorchTask) -> int:
    return _scorch_parallel_mod.task_allowed_count(task)


def _task_estimate_seconds(task: ScorchTask, thread_budget: int) -> float:
    return _scorch_parallel_mod.task_estimate_seconds(task, thread_budget)


def _elastic_scorch_threads(
    *,
    base_threads: int,
    allowed_count: int,
    free_cores: int,
    scheduler_queue_depth: int,
    local_queue_depth: int,
) -> int:
    return _scorch_parallel_mod.elastic_scorch_threads(
        base_threads=base_threads,
        allowed_count=allowed_count,
        free_cores=free_cores,
        scheduler_queue_depth=scheduler_queue_depth,
        local_queue_depth=local_queue_depth,
    )


def _split_task_for_rechunk(task: ScorchTask, suffix_seed: str) -> List[ScorchTask]:
    return _scorch_parallel_mod.split_task_for_rechunk(task, suffix_seed)
