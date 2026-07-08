from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

from post_docking.rescoring.scorch_orchestration_types import (
    ComboKey,
    OrchestrationDeps,
    OrchestrationLimits,
    receptor_path_for_combo,
)
from post_docking.rescoring.scorch_types import ScorchTask, SelectionResult, StageSpec


@dataclass
class TaskBuildState:
    tasks: List[ScorchTask]
    combos_with_tasks: Set[ComboKey]
    combo_modes: Dict[ComboKey, Set[str]]
    combo_failed_local: Dict[ComboKey, bool]
    skipped_missing_receptor: int
    skipped_missing_consensus: int


def build_orchestration_tasks(
    *,
    cfg: dict[str, object],
    include_fda: bool,
    combos: Set[ComboKey],
    specs: List[StageSpec],
    run_root: Path,
    post_run_root: Path,
    processed_root: Path,
    source_post_run_root: Path,
    top_fraction: float,
    jobs: int,
    threads: int,
    component: str,
    logger: logging.Logger,
    limits: OrchestrationLimits,
    deps: OrchestrationDeps,
) -> TaskBuildState:
    base_chunk_size = int(limits.chunk_size)
    tasks: List[ScorchTask] = []
    skipped_missing_receptor = 0
    skipped_missing_consensus = 0
    combos_with_tasks: Set[ComboKey] = set()
    combo_modes: Dict[ComboKey, Set[str]] = {}
    combo_failed_local: Dict[ComboKey, bool] = {}
    control_cache: Dict[str, Set[str]] = {}
    source_post_run_root = source_post_run_root or post_run_root

    for combo in sorted(combos):
        before_task_count = len(tasks)
        free_cores, scheduler_queue_depth = deps.scheduler_runtime_snapshot(cfg)
        tail_mode = len(tasks) <= max(2, int(jobs) * 2)
        pdb_id, variant, ph = combo
        receptor = receptor_path_for_combo(processed_root, combo)
        if not receptor.exists():
            logger.error(
                "%s action=score status=skip reason=missing_receptor pdb_id=%s variant=%s ph=%s receptor=%s",
                component,
                pdb_id,
                variant,
                ph,
                receptor,
            )
            skipped_missing_receptor += 1
            continue

        dock_combo_dir = run_root / pdb_id / variant / ph
        consensus_csv = (
            deps.find_consensus_csv(dock_combo_dir)
            if deps.find_consensus_csv is not None
            else None
        )
        if consensus_csv is None:
            consensus_csv = dock_combo_dir / "consensus_docking_scores.csv"

        if pdb_id not in control_cache:
            control_cache[pdb_id] = deps.load_control_bases(processed_root, pdb_id, logger)

        if not consensus_csv.exists() or consensus_csv.stat().st_size == 0:
            logger.warning(
                "%s action=select status=skip reason=missing_consensus pdb_id=%s variant=%s ph=%s path=%s",
                component,
                pdb_id,
                variant,
                ph,
                consensus_csv,
            )
            skipped_missing_consensus += 1
            continue

        mode_dirs = deps.discover_mode_dirs(
            combo,
            run_root,
            source_post_run_root,
            specs,
            logger,
        )
        if mode_dirs["dud"]:
            logger.info(
                "[scorch.run] mode=dud stage_root=%s n_allowed=%d",
                mode_dirs["dud"][0],
                len(control_cache[pdb_id]),
            )
            for spec in specs:
                if spec.source in {"vina", "gnina"}:
                    stage_root = run_root / pdb_id / variant / ph
                else:
                    stage_root = post_run_root / pdb_id / variant / ph
                stage_dirs_override: Optional[Sequence[str]] = deps.stage_dir_candidates(
                    spec.source,
                    "dud",
                    stage_root,
                )
                generic_dud_stage_dirs = tuple(
                    stage_dir
                    for stage_dir in deps.stage_dir_candidates(
                        spec.source,
                        "fda",
                        stage_root,
                    )
                    if (stage_root / stage_dir) in mode_dirs["dud"]
                )
                if generic_dud_stage_dirs:
                    stage_dirs_override = generic_dud_stage_dirs
                score_csv, score_cols, higher_is_better = deps.score_csv_for_spec(
                    run_root,
                    combo,
                    spec,
                    "dud",
                )
                selection = (
                    deps.select_top_bases_from_score_csv(
                        score_csv if score_csv is not None else Path(""),
                        top_fraction,
                        control_cache[pdb_id],
                        higher_is_better=higher_is_better,
                        logger=logger,
                        score_cols=score_cols,
                    )
                    if score_csv is not None
                    else SelectionResult(
                        allowed_bases=set(control_cache[pdb_id]),
                        n_pool=0,
                        k=0,
                        controls_total=len(control_cache[pdb_id]),
                        selected_stage_by_base={},
                        selected_score_by_base={},
                    )
                )
                allowed_bases = selection.allowed_bases
                n_pool = selection.n_pool
                k = selection.k
                controls_total = selection.controls_total
                logger.info(
                    "[select.pool] mode=dud engine=%s n_controls=%d n_pool=%d k=%d",
                    spec.source,
                    controls_total,
                    n_pool,
                    k,
                )
                logger.info(
                    "[select.allowed] mode=dud engine=%s score_csv=%s n_candidates=%d n_allowed=%d",
                    spec.source,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    len(allowed_bases),
                )
                logger.info(
                    "%s action=select source=%s run_mode=dud pdb_id=%s variant=%s ph=%s score_csv=%s n_candidates=%d k=%d controls=%d allowed_total=%d",
                    component,
                    spec.source,
                    pdb_id,
                    variant,
                    ph,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    k,
                    controls_total,
                    len(allowed_bases),
                )
                adaptive_chunk_size = deps.adaptive_chunk_size(
                    allowed_count=len(allowed_bases),
                    base_chunk_size=base_chunk_size,
                    free_cores=free_cores,
                    scheduler_queue_depth=scheduler_queue_depth,
                    local_queue_depth=len(tasks),
                    jobs=max(1, int(jobs)),
                    want_threads=max(1, int(threads)),
                )
                dud_chunks = deps.chunk_allowed_bases(allowed_bases, adaptive_chunk_size)
                dud_chunks = deps.tail_split_allowed_chunks(
                    dud_chunks,
                    tail_mode=tail_mode,
                    min_split=max(int(limits.chunk_min), adaptive_chunk_size // 2),
                )
                if not dud_chunks:
                    dud_chunks = [set(allowed_bases)]
                dud_partition = deps.summarize_chunk_partition(allowed_bases, dud_chunks)
                logger.info(
                    "[select.allowed.partition] mode=dud engine=%s n_allowed=%d chunk_count=%d chunk_total=%d chunk_unique=%d duplicate=%d preserved=%s",
                    spec.source,
                    len(allowed_bases),
                    dud_partition.chunk_count,
                    dud_partition.chunk_total,
                    dud_partition.chunk_unique,
                    dud_partition.duplicate_count,
                    str(dud_partition.preserved).lower(),
                )
                if not dud_partition.preserved:
                    logger.warning(
                        "[select.allowed.partition] mode=dud engine=%s status=warn missing=%d extra=%d",
                        spec.source,
                        dud_partition.missing_count,
                        dud_partition.extra_count,
                    )
                for chunk_idx, chunk_bases in enumerate(dud_chunks):
                    chunk_tag = f"part{chunk_idx:03d}" if len(dud_chunks) > 1 else None
                    tasks.append(
                        (
                            spec,
                            combo,
                            receptor,
                            chunk_bases,
                            control_cache[pdb_id],
                            "dud",
                            stage_dirs_override,
                            score_csv,
                            {
                                base: stage
                                for base, stage in selection.selected_stage_by_base.items()
                                if base in chunk_bases
                            },
                            {
                                base: score
                                for base, score in selection.selected_score_by_base.items()
                                if base in chunk_bases
                            },
                            chunk_tag,
                        )
                    )
            combo_modes.setdefault(combo, set()).add("dud")
            combo_failed_local.setdefault(combo, False)

        if include_fda and mode_dirs["fda"]:
            logger.info(
                "[scorch.run] mode=fda stage_root=%s n_allowed=%d",
                mode_dirs["fda"][0],
                len(control_cache[pdb_id]),
            )
            for spec in specs:
                score_csv, score_cols, higher_is_better = deps.score_csv_for_spec(
                    run_root,
                    combo,
                    spec,
                    "fda",
                )
                selection = (
                    deps.select_top_bases_from_score_csv(
                        score_csv if score_csv is not None else Path(""),
                        top_fraction,
                        control_cache[pdb_id],
                        higher_is_better=higher_is_better,
                        logger=logger,
                        score_cols=score_cols,
                    )
                    if score_csv is not None
                    else SelectionResult(
                        allowed_bases=set(control_cache[pdb_id]),
                        n_pool=0,
                        k=0,
                        controls_total=len(control_cache[pdb_id]),
                        selected_stage_by_base={},
                        selected_score_by_base={},
                    )
                )
                allowed_bases = selection.allowed_bases
                n_pool = selection.n_pool
                k = selection.k
                controls_total = selection.controls_total
                logger.info(
                    "[select.pool] mode=fda engine=%s n_controls=%d n_pool=%d k=%d",
                    spec.source,
                    controls_total,
                    n_pool,
                    k,
                )
                logger.info(
                    "[select.allowed] mode=fda engine=%s score_csv=%s n_candidates=%d n_allowed=%d",
                    spec.source,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    len(allowed_bases),
                )
                logger.info(
                    "%s action=select source=%s run_mode=fda pdb_id=%s variant=%s ph=%s score_csv=%s n_candidates=%d k=%d controls=%d allowed_total=%d",
                    component,
                    spec.source,
                    pdb_id,
                    variant,
                    ph,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    k,
                    controls_total,
                    len(allowed_bases),
                )
                adaptive_chunk_size = deps.adaptive_chunk_size(
                    allowed_count=len(allowed_bases),
                    base_chunk_size=base_chunk_size,
                    free_cores=free_cores,
                    scheduler_queue_depth=scheduler_queue_depth,
                    local_queue_depth=len(tasks),
                    jobs=max(1, int(jobs)),
                    want_threads=max(1, int(threads)),
                )
                fda_chunks = deps.chunk_allowed_bases(allowed_bases, adaptive_chunk_size)
                fda_chunks = deps.tail_split_allowed_chunks(
                    fda_chunks,
                    tail_mode=tail_mode,
                    min_split=max(int(limits.chunk_min), adaptive_chunk_size // 2),
                )
                if not fda_chunks:
                    fda_chunks = [set(allowed_bases)]
                fda_partition = deps.summarize_chunk_partition(allowed_bases, fda_chunks)
                logger.info(
                    "[select.allowed.partition] mode=fda engine=%s n_allowed=%d chunk_count=%d chunk_total=%d chunk_unique=%d duplicate=%d preserved=%s",
                    spec.source,
                    len(allowed_bases),
                    fda_partition.chunk_count,
                    fda_partition.chunk_total,
                    fda_partition.chunk_unique,
                    fda_partition.duplicate_count,
                    str(fda_partition.preserved).lower(),
                )
                if not fda_partition.preserved:
                    logger.warning(
                        "[select.allowed.partition] mode=fda engine=%s status=warn missing=%d extra=%d",
                        spec.source,
                        fda_partition.missing_count,
                        fda_partition.extra_count,
                    )
                for chunk_idx, chunk_bases in enumerate(fda_chunks):
                    chunk_tag = f"part{chunk_idx:03d}" if len(fda_chunks) > 1 else None
                    tasks.append(
                        (
                            spec,
                            combo,
                            receptor,
                            chunk_bases,
                            control_cache[pdb_id],
                            "fda",
                            None,
                            score_csv,
                            {
                                base: stage
                                for base, stage in selection.selected_stage_by_base.items()
                                if base in chunk_bases
                            },
                            {
                                base: score
                                for base, score in selection.selected_score_by_base.items()
                                if base in chunk_bases
                            },
                            chunk_tag,
                        )
                    )
            combo_modes.setdefault(combo, set()).add("fda")
            combo_failed_local.setdefault(combo, False)

        if len(tasks) > before_task_count:
            combos_with_tasks.add(combo)

    if len(tasks) <= max(2, int(jobs) * 2):
        tail_split_tasks: List[ScorchTask] = []
        for (
            spec,
            combo,
            receptor,
            allowed_bases,
            control_bases,
            run_mode,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
            chunk_tag,
        ) in tasks:
            if len(allowed_bases) < int(limits.tail_split_trigger):
                tail_split_tasks.append(
                    (
                        spec,
                        combo,
                        receptor,
                        allowed_bases,
                        control_bases,
                        run_mode,
                        stage_dirs_override,
                        score_csv,
                        selected_stage_by_base,
                        selected_score_by_base,
                        chunk_tag,
                    )
                )
                continue
            split_size = max(int(limits.chunk_min), len(allowed_bases) // 2)
            split_chunks = deps.tail_split_allowed_chunks(
                deps.chunk_allowed_bases(allowed_bases, split_size),
                tail_mode=True,
                min_split=max(int(limits.chunk_min), split_size // 2),
            )
            if len(split_chunks) <= 1:
                tail_split_tasks.append(
                    (
                        spec,
                        combo,
                        receptor,
                        allowed_bases,
                        control_bases,
                        run_mode,
                        stage_dirs_override,
                        score_csv,
                        selected_stage_by_base,
                        selected_score_by_base,
                        chunk_tag,
                    )
                )
                continue
            tag_root = chunk_tag or "tail"
            for part_idx, part_bases in enumerate(split_chunks):
                tail_split_tasks.append(
                    (
                        spec,
                        combo,
                        receptor,
                        part_bases,
                        control_bases,
                        run_mode,
                        stage_dirs_override,
                        score_csv,
                        {
                            base: stage
                            for base, stage in selected_stage_by_base.items()
                            if base in part_bases
                        },
                        {
                            base: score
                            for base, score in selected_score_by_base.items()
                            if base in part_bases
                        },
                        f"{tag_root}t{part_idx:02d}",
                    )
                )
        tasks = tail_split_tasks
    return TaskBuildState(
        tasks=tasks,
        combos_with_tasks=combos_with_tasks,
        combo_modes=combo_modes,
        combo_failed_local=combo_failed_local,
        skipped_missing_receptor=skipped_missing_receptor,
        skipped_missing_consensus=skipped_missing_consensus,
    )
