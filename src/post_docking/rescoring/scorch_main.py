from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from config.normalize import _to_bool
from config.runtime_config import load_config
from docking.global_scheduler import acquire_global_cores
from post_docking.rescoring import scorch_orchestration as _scorch_orchestration_mod
from post_docking.rescoring import scorch_provisional_cache as _scorch_cache_mod
from post_docking.rescoring import rescoring_scorch_support as _scorch_support_mod
from post_docking.rescoring.scorch_device import (
    apply_scorch_device_plan_to_cfg,
    build_scorch_device_plan,
    cap_scorch_jobs_for_device,
    scorch_chunk_size_for_cfg,
)
from post_docking.rescoring.scorch_events import emit_bench_event, emit_task_event
from post_docking.rescoring.scorch_housekeeping_once import run_housekeeping_once
from post_docking.rescoring.scorch_coverage import (
    accept_quarantined_coverage,
    evaluate_scorch_coverage,
    remove_done_sentinel,
    resolve_top_fraction,
    write_coverage_summary,
)
from post_docking.rescoring.scorch_runtime_state import (
    COMPONENT,
    SCORCH_TOP_FRACTION_DEFAULT,
    SCORCH_TOP_FRACTION_KEY,
    _apply_scorch_config_defaults,
    _preflight,
    _resolve_roots,
    configure_logging,
    gpu_adaptive_chunk_size,
    parse_args,
    preserve_allowed_chunks,
    preserve_task_for_rechunk,
)
from post_docking.rescoring.scorch_scoring_bridge import (
    _aggregate_combo,
    _decoy_prefixes_from_test_mode,
    _discover_mode_dirs,
    _filter_combos,
    _filter_done_combos_with_stats,
    _launch_pose_bust_async,
    _load_control_bases,
    _log_no_combos_after_filters,
    _mark_combo_done,
    _parse_test_mode_tokens,
    _posebusters_missing,
    _prep_missing,
    _resolve_decoy_prefix_override,
    _run_prep_for_scorch,
    _score_csv_for_spec,
    _score_stage,
    _select_top_bases_from_score_csv,
    _set_decoy_prefix,
    annotate_scorch_z_scores,
    discover_combos,
    stage_dir_candidates,
)
from post_docking.rescoring.scorch_selection_invariants import summarize_chunk_partition
from post_docking.rescoring.scorch_types import StageSpec

SCORCH_CHUNK_SIZE = _scorch_support_mod.SCORCH_CHUNK_SIZE


def _restore_undercovered_done_combos(
    *,
    cfg: Dict[str, object],
    candidate_combos: Set[Tuple[str, str, str]],
    pending_combos: Set[Tuple[str, str, str]],
    run_root: Path,
    post_run_root: Path,
    overwrite: bool,
    decoy_prefix: str,
    logger,
) -> Set[Tuple[str, str, str]]:
    if overwrite:
        return pending_combos
    restored = set(pending_combos)
    skipped = set(candidate_combos) - set(pending_combos)
    if not skipped:
        return restored
    for combo in sorted(skipped):
        summary = evaluate_scorch_coverage(
            docked_run_root=run_root,
            post_run_root=post_run_root,
            pdb_id=combo[0],
            variant=combo[1],
            ph=combo[2],
            top_fraction=resolve_top_fraction(cfg),
            decoy_prefix=decoy_prefix,
            accept_quarantined=accept_quarantined_coverage(cfg),
        )
        try:
            write_coverage_summary(post_run_root, summary)
        except Exception:
            logger.debug(
                "%s action=coverage_summary status=failed combo=%s",
                COMPONENT,
                combo,
                exc_info=True,
            )
        if summary.all_complete:
            continue
        try:
            remove_done_sentinel(post_run_root, summary)
        except Exception:
            logger.debug(
                "%s action=clear_done status=failed combo=%s",
                COMPONENT,
                combo,
                exc_info=True,
            )
        restored.add(combo)
        logger.warning(
            "%s action=idempotency status=retry reason=scorch_coverage_incomplete combo=%s details=%s",
            COMPONENT,
            combo,
            ",".join(summary.reasons) or "unknown",
        )
    return restored
SCORCH_CHUNK_MIN = _scorch_support_mod.SCORCH_CHUNK_MIN
SCORCH_TAIL_SPLIT_TRIGGER = _scorch_support_mod.SCORCH_TAIL_SPLIT_TRIGGER
SCORCH_RUNTIME_RECHUNK_STALE_SEC = _scorch_support_mod.SCORCH_RUNTIME_RECHUNK_STALE_SEC
SCORCH_HEDGE_TIMEOUT_SEC = _scorch_support_mod.SCORCH_HEDGE_TIMEOUT_SEC
SCORCH_HEDGE_MAX_INFLIGHT = _scorch_support_mod.SCORCH_HEDGE_MAX_INFLIGHT
SCORCH_HEDGE_MIN_ALLOWED = _scorch_support_mod.SCORCH_HEDGE_MIN_ALLOWED
SCORCH_PARALLEL_PROFILE_DEFAULT = _scorch_support_mod.SCORCH_PARALLEL_PROFILE_DEFAULT
SCORCH_PARALLEL_PROFILE_ENV = _scorch_support_mod.SCORCH_PARALLEL_PROFILE_ENV
_normalize_parallel_profile = _scorch_support_mod._normalize_parallel_profile
_allocate_scorch_parallelism = _scorch_support_mod._allocate_scorch_parallelism
_adaptive_chunk_size = _scorch_support_mod._adaptive_chunk_size
_tail_split_allowed_chunks = _scorch_support_mod._tail_split_allowed_chunks
_chunk_allowed_bases = _scorch_support_mod._chunk_allowed_bases
_task_allowed_count = _scorch_support_mod._task_allowed_count
_task_estimate_seconds = _scorch_support_mod._task_estimate_seconds
_elastic_scorch_threads = _scorch_support_mod._elastic_scorch_threads
_split_task_for_rechunk = _scorch_support_mod._split_task_for_rechunk
_decoy_prefix_value = _scorch_support_mod._decoy_prefix_value
_scheduler_runtime_snapshot = _scorch_support_mod._scheduler_runtime_snapshot

_find_consensus_csv_import = None
_rerank_consensus_with_scorch_import = None


def _engine_source_enabled(cfg: Dict[str, object], source: str) -> bool:
    source_norm = str(source or "").strip().lower()
    if source_norm == "vina":
        return True
    flag_by_source = {
        "gnina": "USE_GNINA",
        "ledock": "USE_LEDOCK",
        "dock6": "USE_DOCK6",
    }
    flag = flag_by_source.get(source_norm)
    if not flag:
        return True
    raw = os.environ.get(flag)
    if raw is None:
        raw = cfg.get(flag)
    return _to_bool(raw, False)


def _resolve_standalone_cli_parallelism(
    *,
    requested_jobs: int,
    requested_threads: int,
    cpu_threads: int,
    profile: str,
) -> tuple[int, int, float]:
    req_jobs = max(1, int(requested_jobs))
    req_threads = max(1, int(requested_threads))
    cpu_cap = max(1, int(cpu_threads))
    _alloc_jobs, alloc_threads, alloc_util = _allocate_scorch_parallelism(
        cpu_budget=cpu_cap,
        free_cores=cpu_cap,
        backlog=req_jobs,
        allowed_count=SCORCH_TAIL_SPLIT_TRIGGER,
        scheduler_queue_depth=0,
        local_queue_depth=0,
        profile=profile,
    )
    effective_threads = max(1, min(req_threads, int(alloc_threads)))
    jobs_cpu_cap = max(1, int(cpu_cap) // int(effective_threads))
    effective_jobs = max(1, min(req_jobs, int(jobs_cpu_cap)))
    return int(effective_jobs), int(effective_threads), float(alloc_util)


try:
    from post_docking.rescoring.rescore_reranker import (
        find_consensus_csv as _find_consensus_csv_loaded,
        rerank_consensus_with_scorch as _rerank_consensus_with_scorch_loaded,
    )

    _find_consensus_csv_import = _find_consensus_csv_loaded
    _rerank_consensus_with_scorch_import = _rerank_consensus_with_scorch_loaded
except Exception:  # pragma: no cover - optional dependency
    pass

find_consensus_csv_fn = _find_consensus_csv_import
rerank_consensus_with_scorch_fn = _rerank_consensus_with_scorch_import


def main(
    argv: Optional[Sequence[str]] = None,
    *,
    cfg_override: Optional[Dict[str, object]] = None,
) -> int:
    args = parse_args(argv)
    logger = configure_logging(args.verbose)

    if cfg_override is None:
        try:
            cfg = load_config()
        except Exception as exc:
            logger.error(
                "%s action=preflight status=failed reason=config_load_error error=%s",
                COMPONENT,
                exc,
            )
            return 1
    else:
        cfg = dict(cfg_override)

    repo_root_for_defaults = (
        Path(args.repo_root).expanduser().resolve()
        if args.repo_root
        else Path(__file__).resolve().parents[3]
    )
    _apply_scorch_config_defaults(cfg, repo_root=repo_root_for_defaults)

    raw_use = cfg.get("USE_SCORCH", False)
    use_scorch = False
    if isinstance(raw_use, bool):
        use_scorch = raw_use
    else:
        s = str(raw_use).strip().lower()
        use_scorch = s in ("1", "true", "yes", "on")

    top_fraction = SCORCH_TOP_FRACTION_DEFAULT
    if SCORCH_TOP_FRACTION_KEY in cfg:
        try:
            candidate = float(
                cfg.get(SCORCH_TOP_FRACTION_KEY, SCORCH_TOP_FRACTION_DEFAULT)
            )
            if 0 < candidate <= 1:
                top_fraction = candidate
            else:
                logger.warning(
                    "%s action=preflight status=warn reason=invalid_top_fraction value=%s using_default=%.2f",
                    COMPONENT,
                    candidate,
                    SCORCH_TOP_FRACTION_DEFAULT,
                )
        except Exception as exc:
            logger.warning(
                "%s action=preflight status=warn reason=parse_top_fraction_failed error=%s using_default=%.2f",
                COMPONENT,
                exc,
                SCORCH_TOP_FRACTION_DEFAULT,
            )

    if not use_scorch:
        logger.info(
            "%s action=skip status=ok reason=USE_SCORCH_false run_id=%s",
            COMPONENT,
            args.run_id,
        )
        return 0

    device_plan = build_scorch_device_plan(cfg, logger)
    cfg = apply_scorch_device_plan_to_cfg(cfg, device_plan)

    scheduler_present = bool(cfg.get("GLOBAL_DOCKING_SCHEDULER"))
    try:
        cpu_threads = max(1, int(cfg.get("CPU", os.cpu_count() or 1)))
    except Exception:
        cpu_threads = max(1, os.cpu_count() or 1)
    if not scheduler_present:
        profile = _normalize_parallel_profile(
            os.environ.get(
                SCORCH_PARALLEL_PROFILE_ENV,
                cfg.get("SCORCH_PARALLEL_PROFILE", SCORCH_PARALLEL_PROFILE_DEFAULT),
            )
        )
        alloc_jobs, alloc_threads, alloc_util = _resolve_standalone_cli_parallelism(
            requested_jobs=int(args.jobs),
            requested_threads=int(args.threads),
            cpu_threads=int(cpu_threads),
            profile=profile,
        )
        if int(args.threads) != int(alloc_threads):
            logger.info(
                "%s action=threads source=allocator profile=%s util_target=%.2f requested=%d effective=%d",
                COMPONENT,
                profile,
                alloc_util,
                args.threads,
                alloc_threads,
            )
        args.threads = max(1, int(alloc_threads))
        if int(args.jobs) != int(alloc_jobs):
            logger.info(
                "%s action=jobs source=allocator cpu=%d profile=%s util_target=%.2f requested=%d effective=%d",
                COMPONENT,
                cpu_threads,
                profile,
                alloc_util,
                args.jobs,
                alloc_jobs,
            )
        args.jobs = max(1, int(alloc_jobs))

    capped_jobs = cap_scorch_jobs_for_device(int(args.jobs), device_plan)
    if capped_jobs != int(args.jobs):
        logger.info(
            "%s action=jobs source=device requested=%d effective=%d device=%s backend=%s gpu_count=%d",
            COMPONENT,
            args.jobs,
            capped_jobs,
            device_plan.effective,
            device_plan.backend,
            device_plan.gpu_count,
        )
        args.jobs = capped_jobs

    if not _preflight(cfg, logger):
        return 1

    repo_root, docked_root, post_root, processed_root = _resolve_roots(args, cfg)
    run_root = docked_root / args.run_id
    public_post_run_root = post_root / args.run_id
    provisional_run = _scorch_cache_mod.provisional_run_enabled(cfg)
    if provisional_run:
        args.skip_autofix = True
        post_root = _scorch_cache_mod.artifact_post_root(cfg, args.run_id)
    post_run_root = post_root / args.run_id
    if not run_root.exists():
        logger.error(
            "%s action=preflight status=failed reason=missing_run_root path=%s",
            COMPONENT,
            run_root,
        )
        return 1

    all_specs: List[StageSpec] = [
        StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv"),
        StageSpec("gnina", "gnina_best", "scorch_scores_gnina_best.csv"),
        StageSpec("ledock", "ledock_pdbqt", "scorch_scores_ledock.csv"),
        StageSpec("dock6", "dock6_pdbqt", "scorch_scores_dock6.csv"),
    ]
    specs = [spec for spec in all_specs if _engine_source_enabled(cfg, spec.source)]
    disabled_sources = [spec.source for spec in all_specs if spec not in specs]
    if disabled_sources:
        logger.info(
            "%s action=engine_filter disabled=%s enabled=%s",
            COMPONENT,
            ",".join(disabled_sources),
            ",".join(spec.source for spec in specs),
        )

    override_prefix = _resolve_decoy_prefix_override(args)
    decoy_prefixes = _decoy_prefixes_from_test_mode(cfg, override_prefix)
    if override_prefix:
        logger.info(
            "%s action=decoy_prefixes source=override prefixes=%s",
            COMPONENT,
            ",".join(decoy_prefixes),
        )
    else:
        tokens = _parse_test_mode_tokens(cfg)
        logger.info(
            "%s action=decoy_prefixes source=test_mode tokens=%s prefixes=%s",
            COMPONENT,
            ",".join(tokens),
            ",".join(decoy_prefixes),
        )

    ran_pose_bust = False
    any_failed = False
    combo_had_tasks: Set[Tuple[str, str, str]] = set()
    combo_failed: Dict[Tuple[str, str, str], bool] = {}

    for idx, decoy_prefix in enumerate(decoy_prefixes):
        _set_decoy_prefix(decoy_prefix, logger)
        discovered = discover_combos(
            run_root,
            public_post_run_root,
            pdb_id=args.pdb_id,
            variant=args.variant,
            ph=args.ph,
        )
        candidate_combos = _filter_combos(
            discovered,
            pdb_id=args.pdb_id,
            variant=args.variant,
            ph=args.ph,
        )
        combos, done_stats = _filter_done_combos_with_stats(
            candidate_combos,
            post_run_root,
            args.overwrite,
            logger,
        )
        combos = _restore_undercovered_done_combos(
            cfg=cfg,
            candidate_combos=candidate_combos,
            pending_combos=combos,
            run_root=run_root,
            post_run_root=post_run_root,
            overwrite=bool(args.overwrite),
            decoy_prefix=str(decoy_prefix),
            logger=logger,
        )
        logger.info(
            "%s action=discover status=ok combos_discovered=%d combos_filtered=%d decoy_prefix=%s filters[pdb_id=%s variant=%s ph=%s]",
            COMPONENT,
            len(discovered),
            len(combos),
            decoy_prefix,
            args.pdb_id or "*",
            args.variant or "*",
            args.ph or "*",
        )
        if not combos:
            _log_no_combos_after_filters(
                logger,
                run_id=str(args.run_id),
                decoy_prefix=str(decoy_prefix),
                pdb_id_filter=str(args.pdb_id or "*"),
                variant_filter=str(args.variant or "*"),
                ph_filter=str(args.ph or "*"),
                idempotent_only=(
                    int(done_stats.get("combos_input", 0)) > 0
                    and int(done_stats.get("combos_skipped", 0))
                    >= int(done_stats.get("combos_input", 0))
                ),
                combos_input=int(done_stats.get("combos_input", 0)),
                combos_skipped=int(done_stats.get("combos_skipped", 0)),
                filtered_scope=bool(args.pdb_id or args.variant or args.ph),
            )
            continue

        if not args.skip_autofix:
            if not ran_pose_bust and _posebusters_missing(post_run_root, combos):
                pose_task_id = f"scorch:housekeeping:pose_bust:{args.run_id}:{decoy_prefix}"
                pose_want_cores = max(1, min(2, int(args.threads)))
                emit_task_event(
                    "TASK_SUBMITTED",
                    task_id=pose_task_id,
                    task_type="housekeeping",
                    want_cores=pose_want_cores,
                    min_cores=1,
                    est_duration_sec=2.0,
                )
                emit_task_event(
                    "TASK_START",
                    task_id=pose_task_id,
                    task_type="housekeeping",
                    want_cores=pose_want_cores,
                    min_cores=1,
                    granted_cores=0,
                    est_duration_sec=2.0,
                )
                pose_step = f"pose_bust_{decoy_prefix}"

                def _pose_bust_action() -> None:
                    _launch_pose_bust_async(
                        args.run_id,
                        repo_root,
                        docked_root,
                        post_root,
                        args.overwrite,
                        pose_want_cores,
                        logger,
                    )

                run_housekeeping_once(
                    repo_root=repo_root,
                    run_id=str(args.run_id),
                    step=pose_step,
                    action=_pose_bust_action,
                    logger=logger,
                )
                emit_task_event(
                    "TASK_END",
                    task_id=pose_task_id,
                    task_type="housekeeping",
                    want_cores=pose_want_cores,
                    min_cores=1,
                    granted_cores=0,
                    est_duration_sec=2.0,
                    status="launched",
                )
                ran_pose_bust = True
            if _prep_missing(run_root, post_run_root, combos):
                prep_workers = max(1, int(args.jobs))
                prep_task_id = (
                    f"scorch:materialize:prep:{args.run_id}:{args.pdb_id or '*'}:{args.variant or '*'}:{args.ph or '*'}:{decoy_prefix}"
                )
                emit_task_event(
                    "TASK_SUBMITTED",
                    task_id=prep_task_id,
                    task_type="scorch_materialize",
                    want_cores=int(prep_workers),
                    min_cores=1,
                    est_duration_sec=2.0,
                )
                with acquire_global_cores(
                    cfg,
                    cores=int(prep_workers),
                    min_cores=1,
                    priority=1,
                    est_duration_sec=2.0,
                    task_id=prep_task_id,
                    task_type="scorch_materialize",
                ) as granted:
                    emit_task_event(
                        "TASK_START",
                        task_id=prep_task_id,
                        task_type="scorch_materialize",
                        want_cores=int(prep_workers),
                        min_cores=1,
                        granted_cores=int(granted),
                        est_duration_sec=2.0,
                    )
                    _run_prep_for_scorch(
                        args.run_id,
                        repo_root,
                        docked_root,
                        post_root,
                        args.overwrite,
                        logger,
                        decoy_prefix=decoy_prefix,
                        pdb_id=args.pdb_id,
                        variant=args.variant,
                        ph=args.ph,
                        max_workers=int(max(1, granted)),
                    )
                    emit_task_event(
                        "TASK_END",
                        task_id=prep_task_id,
                        task_type="scorch_materialize",
                        want_cores=int(prep_workers),
                        min_cores=1,
                        granted_cores=int(granted),
                        est_duration_sec=2.0,
                        status="ok",
                    )

        discovered = discover_combos(
            run_root,
            public_post_run_root,
            pdb_id=args.pdb_id,
            variant=args.variant,
            ph=args.ph,
        )
        candidate_combos = _filter_combos(
            discovered,
            pdb_id=args.pdb_id,
            variant=args.variant,
            ph=args.ph,
        )
        combos, done_stats = _filter_done_combos_with_stats(
            candidate_combos,
            post_run_root,
            args.overwrite,
            logger,
        )
        combos = _restore_undercovered_done_combos(
            cfg=cfg,
            candidate_combos=candidate_combos,
            pending_combos=combos,
            run_root=run_root,
            post_run_root=post_run_root,
            overwrite=bool(args.overwrite),
            decoy_prefix=str(decoy_prefix),
            logger=logger,
        )
        if not combos:
            _log_no_combos_after_filters(
                logger,
                run_id=str(args.run_id),
                decoy_prefix=str(decoy_prefix),
                pdb_id_filter=str(args.pdb_id or "*"),
                variant_filter=str(args.variant or "*"),
                ph_filter=str(args.ph or "*"),
                idempotent_only=(
                    int(done_stats.get("combos_input", 0)) > 0
                    and int(done_stats.get("combos_skipped", 0))
                    >= int(done_stats.get("combos_input", 0))
                ),
                combos_input=int(done_stats.get("combos_input", 0)),
                combos_skipped=int(done_stats.get("combos_skipped", 0)),
                filtered_scope=bool(args.pdb_id or args.variant or args.ph),
            )
            continue

        include_fda = override_prefix is not None or idx == 0
        scorch_chunk_size = scorch_chunk_size_for_cfg(cfg, SCORCH_CHUNK_SIZE)
        gpu_chunking = (
            str(cfg.get("SCORCH_DEVICE_EFFECTIVE", "cpu")).strip().lower() == "gpu"
            and int(scorch_chunk_size) > int(SCORCH_CHUNK_SIZE)
        )
        tail_split_chunks = (
            preserve_allowed_chunks if gpu_chunking else _tail_split_allowed_chunks
        )
        adaptive_chunk_size = gpu_adaptive_chunk_size if gpu_chunking else _adaptive_chunk_size
        split_for_rechunk = (
            preserve_task_for_rechunk if gpu_chunking else _split_task_for_rechunk
        )
        tail_split_trigger = (
            max(int(SCORCH_TAIL_SPLIT_TRIGGER), int(scorch_chunk_size) + 1)
            if gpu_chunking
            else int(SCORCH_TAIL_SPLIT_TRIGGER)
        )
        logger.info(
            "%s action=chunking device=%s chunk_size=%d gpu_chunking=%s min_gpu_ligands=%s tail_split_trigger=%d",
            COMPONENT,
            str(cfg.get("SCORCH_DEVICE_EFFECTIVE", "cpu")),
            int(scorch_chunk_size),
            bool(gpu_chunking),
            str(cfg.get("SCORCH_GPU_MIN_LIGANDS_EFFECTIVE", "")),
            int(tail_split_trigger),
        )
        orchestration_limits = _scorch_orchestration_mod.OrchestrationLimits(
            chunk_size=scorch_chunk_size,
            chunk_min=SCORCH_CHUNK_MIN,
            tail_split_trigger=tail_split_trigger,
            runtime_rechunk_stale_sec=SCORCH_RUNTIME_RECHUNK_STALE_SEC,
            hedge_timeout_sec=SCORCH_HEDGE_TIMEOUT_SEC,
            hedge_max_inflight=SCORCH_HEDGE_MAX_INFLIGHT,
            hedge_min_allowed=SCORCH_HEDGE_MIN_ALLOWED,
        )
        orchestration_deps = _scorch_orchestration_mod.OrchestrationDeps(
            scheduler_runtime_snapshot=_scheduler_runtime_snapshot,
            discover_mode_dirs=_discover_mode_dirs,
            load_control_bases=_load_control_bases,
            stage_dir_candidates=stage_dir_candidates,
            score_csv_for_spec=_score_csv_for_spec,
            select_top_bases_from_score_csv=_select_top_bases_from_score_csv,
            adaptive_chunk_size=adaptive_chunk_size,
            chunk_allowed_bases=_chunk_allowed_bases,
            tail_split_allowed_chunks=tail_split_chunks,
            score_stage=_score_stage,
            task_estimate_seconds=_task_estimate_seconds,
            task_allowed_count=_task_allowed_count,
            split_task_for_rechunk=split_for_rechunk,
            elastic_scorch_threads=_elastic_scorch_threads,
            aggregate_combo=_aggregate_combo,
            annotate_scorch_z_scores=annotate_scorch_z_scores,
            emit_task_event=emit_task_event,
            emit_bench_event=emit_bench_event,
            summarize_chunk_partition=summarize_chunk_partition,
            find_consensus_csv=find_consensus_csv_fn,
            rerank_consensus_with_scorch=rerank_consensus_with_scorch_fn,
        )
        orchestration_result = _scorch_orchestration_mod.run_combo_task_orchestration(
            cfg=cfg,
            include_fda=include_fda,
            combos=combos,
            specs=specs,
            run_root=run_root,
            post_run_root=post_run_root,
            source_post_run_root=public_post_run_root,
            processed_root=processed_root,
            top_fraction=top_fraction,
            jobs=int(args.jobs),
            threads=int(args.threads),
            overwrite=bool(args.overwrite),
            decoy_prefix=_decoy_prefix_value(),
            component=COMPONENT,
            logger=logger,
            limits=orchestration_limits,
            deps=orchestration_deps,
        )
        if orchestration_result.failed_jobs != 0:
            any_failed = True
        for combo in orchestration_result.combos_with_tasks:
            combo_had_tasks.add(combo)
            if combo not in orchestration_result.local_completed_combos:
                combo_failed[combo] = True
            elif combo not in combo_failed:
                combo_failed[combo] = False
            if orchestration_result.combo_failed_local.get(combo, False):
                combo_failed[combo] = True

    done_sentinel_deps = _scorch_orchestration_mod.DoneSentinelDeps(
        emit_task_event=emit_task_event,
        mark_combo_done=_mark_combo_done,
    )
    coverage_incomplete_combos = _scorch_orchestration_mod.mark_done_sentinels(
        cfg=cfg,
        run_root=run_root,
        post_run_root=post_run_root,
        combos_with_tasks=combo_had_tasks,
        combo_failed=combo_failed,
        decoy_prefix=_decoy_prefix_value(),
        component=COMPONENT,
        logger=logger,
        deps=done_sentinel_deps,
    )
    if coverage_incomplete_combos:
        any_failed = True

    return 1 if any_failed else 0
