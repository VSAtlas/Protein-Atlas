from __future__ import annotations

import logging
import shutil
import time
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, cast

from post_docking.rescoring.scorch_orchestration_types import ComboKey, OrchestrationDeps, OrchestrationLimits
from post_docking.rescoring.scorch_types import ScorchTask


@dataclass
class TaskExecuteResult:
    failed_jobs: int
    completed: int
    combo_failed_local: Dict[ComboKey, bool]


def execute_orchestration_tasks(
    *,
    cfg: dict[str, object],
    tasks: List[ScorchTask],
    run_root: Path,
    source_post_run_root: Path,
    jobs: int,
    threads: int,
    overwrite: bool,
    component: str,
    logger: logging.Logger,
    limits: OrchestrationLimits,
    deps: OrchestrationDeps,
    combo_failed_local: Dict[ComboKey, bool],
) -> TaskExecuteResult:
    completed = 0
    failed_jobs = 0

    if jobs <= 1:
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
            ok, _ = deps.score_stage(
                cfg,
                spec,
                combo,
                run_root,
                source_post_run_root,
                receptor,
                threads,
                overwrite,
                logger,
                allowed_bases,
                control_bases=control_bases,
                run_mode=run_mode,
                stage_dirs_override=stage_dirs_override,
                score_csv=score_csv,
                selected_stage_by_base=selected_stage_by_base,
                selected_score_by_base=selected_score_by_base,
                chunk_tag=chunk_tag,
            )
            if ok:
                completed += 1
            else:
                failed_jobs += 1
                combo_failed_local[combo] = True
    else:
        def _task_key(task: ScorchTask) -> str:
            spec, combo, _receptor, _allowed, _ctrl, run_mode, _dirs, _score, _selected_stage, _selected_score, chunk_tag = task
            return f"{combo[0]}:{combo[1]}:{combo[2]}:{spec.source}:{run_mode}:{chunk_tag or 'all'}"

        def _task_sort_key(task: ScorchTask) -> Tuple[float, int]:
            return (
                deps.task_estimate_seconds(task, max(1, int(threads))),
                deps.task_allowed_count(task),
            )

        pending: List[ScorchTask] = sorted(list(tasks), key=_task_sort_key)
        running: Dict[Future[Tuple[bool, Optional[Path]]], Dict[str, Any]] = {}
        hedged_groups: Set[str] = set()
        group_success: Dict[str, bool] = {}
        group_winner: Dict[str, Dict[str, Any]] = {}
        group_failures: Dict[str, int] = defaultdict(int)
        stranded_since: Optional[float] = None
        rechunk_generation = 0

        with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
            while pending or running:
                while pending and len(running) < max(1, int(jobs)):
                    task = pending.pop(0)
                    group_id = _task_key(task)
                    if group_success.get(group_id, False):
                        continue
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
                    ) = task
                    free_cores, scheduler_queue_depth = deps.scheduler_runtime_snapshot(cfg)
                    elastic_threads = deps.elastic_scorch_threads(
                        base_threads=max(1, int(threads)),
                        allowed_count=len(allowed_bases),
                        free_cores=free_cores,
                        scheduler_queue_depth=scheduler_queue_depth,
                        local_queue_depth=len(pending),
                    )
                    fut = pool.submit(
                        deps.score_stage,
                        cfg,
                        spec,
                        combo,
                        run_root,
                        source_post_run_root,
                        receptor,
                        elastic_threads,
                        overwrite,
                        logger,
                        allowed_bases,
                        control_bases,
                        run_mode,
                        stage_dirs_override,
                        score_csv,
                        selected_stage_by_base,
                        selected_score_by_base,
                        chunk_tag,
                        group_id,
                        "primary",
                        rechunk_generation,
                    )
                    running[fut] = {
                        "task": task,
                        "combo": combo,
                        "run_mode": run_mode,
                        "group_id": group_id,
                        "role": "primary",
                        "start_ts": time.monotonic(),
                        "chunk_tag": chunk_tag,
                        "canonical_chunk_tag": chunk_tag,
                    }

                if pending:
                    free_cores, _scheduler_queue_depth = deps.scheduler_runtime_snapshot(cfg)
                    has_big_pending = any(
                        deps.task_allowed_count(t) >= int(limits.tail_split_trigger)
                        for t in pending
                    )
                    big_pending_n = sum(
                        1
                        for t in pending
                        if deps.task_allowed_count(t) >= int(limits.tail_split_trigger)
                    )
                    pending_n = max(1, len(pending))
                    mostly_big_pending = big_pending_n >= max(1, (pending_n * 2) // 3)
                    if (
                        free_cores > 0
                        and has_big_pending
                        and (mostly_big_pending or len(running) <= max(1, int(jobs) // 2))
                    ):
                        now = time.monotonic()
                        if stranded_since is None:
                            stranded_since = now
                        elif now - stranded_since >= float(limits.runtime_rechunk_stale_sec):
                            largest_idx = max(
                                range(len(pending)),
                                key=lambda i: deps.task_allowed_count(pending[i]),
                            )
                            original = pending.pop(largest_idx)
                            rechunk_generation += 1
                            split = deps.split_task_for_rechunk(
                                original,
                                suffix_seed=f"{rechunk_generation:03d}",
                            )
                            if len(split) > 1:
                                pending.extend(split)
                                pending.sort(key=_task_sort_key)
                                logger.info(
                                    "%s action=rechunk status=ok original_allowed=%d split_into=%d",
                                    component,
                                    deps.task_allowed_count(original),
                                    len(split),
                                )
                            else:
                                pending.insert(largest_idx, original)
                            stranded_since = now
                    else:
                        stranded_since = None

                if running and (
                    len(hedged_groups) < int(limits.hedge_max_inflight)
                    and (len(running) + len(pending)) <= max(2, int(jobs))
                ):
                    now = time.monotonic()
                    hedge_candidate: Optional[
                        Tuple[Future[Tuple[bool, Optional[Path]]], Dict[str, Any]]
                    ] = None
                    for fut, meta in running.items():
                        if meta.get("role") != "primary":
                            continue
                        group_id = str(meta.get("group_id", ""))
                        if not group_id or group_id in hedged_groups:
                            continue
                        task_any = meta.get("task")
                        if not isinstance(task_any, tuple):
                            continue
                        task = cast(ScorchTask, task_any)
                        if deps.task_allowed_count(task) < int(limits.hedge_min_allowed):
                            continue
                        elapsed = now - float(meta.get("start_ts", now))
                        if elapsed < float(limits.hedge_timeout_sec):
                            continue
                        hedge_candidate = (fut, meta)
                        break

                    if hedge_candidate is not None and len(running) < max(1, int(jobs)):
                        _fut, meta = hedge_candidate
                        task = meta["task"]
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
                        ) = task
                        group_id = str(meta["group_id"])
                        hedge_chunk = f"hedge_{int(now * 1000) % 100000}"
                        hedge_task: ScorchTask = (
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
                            hedge_chunk,
                        )
                        free_cores, scheduler_queue_depth = deps.scheduler_runtime_snapshot(cfg)
                        hedge_threads = deps.elastic_scorch_threads(
                            base_threads=max(1, int(threads)),
                            allowed_count=len(allowed_bases),
                            free_cores=free_cores,
                            scheduler_queue_depth=scheduler_queue_depth,
                            local_queue_depth=len(pending),
                        )
                        hedge_fut = pool.submit(
                            deps.score_stage,
                            cfg,
                            spec,
                            combo,
                            run_root,
                            source_post_run_root,
                            receptor,
                            hedge_threads,
                            overwrite,
                            logger,
                            allowed_bases,
                            control_bases,
                            run_mode,
                            stage_dirs_override,
                            score_csv,
                            selected_stage_by_base,
                            selected_score_by_base,
                            hedge_chunk,
                            group_id,
                            "hedge",
                            rechunk_generation,
                        )
                        running[hedge_fut] = {
                            "task": hedge_task,
                            "combo": combo,
                            "run_mode": run_mode,
                            "group_id": group_id,
                            "role": "hedge",
                            "start_ts": now,
                            "chunk_tag": hedge_chunk,
                            "canonical_chunk_tag": task[10],
                        }
                        hedged_groups.add(group_id)
                        logger.info(
                            "%s action=hedge status=launched group=%s combo=%s",
                            component,
                            group_id,
                            combo,
                        )

                if not running:
                    continue

                done, _ = wait(
                    list(running.keys()),
                    timeout=0.5,
                    return_when=FIRST_COMPLETED,
                )
                if not done:
                    continue

                for future in done:
                    meta = running.pop(future, {})
                    if not meta:
                        continue
                    task = meta["task"]
                    combo = meta["combo"]
                    group_id = str(meta.get("group_id", _task_key(task)))
                    role = str(meta.get("role", "primary"))
                    chunk_tag = meta.get("chunk_tag")
                    canonical_chunk_tag = meta.get("canonical_chunk_tag")
                    try:
                        ok, out_path = future.result()
                    except Exception as exc:  # defensive
                        ok, out_path = False, None
                        spec = task[0]
                        run_mode = task[5]
                        logger.error(
                            "%s action=score status=failed reason=worker_exception source=%s stage=%s combo=%s run_mode=%s error=%s",
                            component,
                            spec.source,
                            spec.stage_dir,
                            combo,
                            run_mode,
                            exc,
                        )

                    if ok:
                        if not group_success.get(group_id, False):
                            group_success[group_id] = True
                            group_winner[group_id] = {
                                "role": role,
                                "out_path": str(out_path) if out_path is not None else "",
                                "chunk_tag": chunk_tag,
                                "canonical_chunk_tag": canonical_chunk_tag,
                            }
                            completed += 1
                        else:
                            deps.emit_bench_event(
                                "HEDGE_CANCELLED",
                                hedge_group_id=group_id,
                                hedge_role=role,
                            )
                            if out_path is not None:
                                try:
                                    Path(out_path).unlink(missing_ok=True)
                                except Exception:
                                    pass
                    else:
                        group_failures[group_id] += 1
                        has_peer = any(
                            str(v.get("group_id", "")) == group_id
                            for v in running.values()
                        )
                        if not has_peer and not group_success.get(group_id, False):
                            failed_jobs += 1
                            combo_failed_local[combo] = True

                    if not any(
                        str(v.get("group_id", "")) == group_id
                        for v in running.values()
                    ):
                        hedged_groups.discard(group_id)
                        winner = group_winner.get(group_id, {})
                        if winner.get("role") == "hedge":
                            raw_src = str(winner.get("out_path", "")).strip()
                            if raw_src:
                                src = Path(raw_src)
                                if src.exists():
                                    hedge_tag = str(winner.get("chunk_tag", "") or "")
                                    canonical_tag = winner.get("canonical_chunk_tag")
                                    dst = src
                                    if hedge_tag:
                                        suffix = f".{hedge_tag}.csv"
                                        if str(src).endswith(suffix):
                                            if canonical_tag:
                                                dst = Path(
                                                    str(src)[: -len(suffix)]
                                                    + f".{canonical_tag}.csv"
                                                )
                                            else:
                                                dst = Path(
                                                    str(src)[: -len(suffix)] + ".csv"
                                                )
                                    if dst != src:
                                        try:
                                            if dst.exists():
                                                dst.unlink(missing_ok=True)
                                        except Exception:
                                            pass
                                        try:
                                            shutil.copy2(src, dst)
                                            deps.emit_bench_event(
                                                "HEDGE_PROMOTED",
                                                hedge_group_id=group_id,
                                                src=str(src),
                                                dst=str(dst),
                                            )
                                        except Exception:
                                            pass
    return TaskExecuteResult(
        failed_jobs=failed_jobs,
        completed=completed,
        combo_failed_local=combo_failed_local,
    )
