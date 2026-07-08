# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import math
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any, Mapping, Optional

from cli.distributed_context import (
    renew_scorch_scope_claim,
    resolve_distributed_context,
    scorch_scope_key,
    try_claim_scorch_scope,
    write_scorch_scope_state,
)
from cli.postrun_hooks_common import _is_no_library_docking
from cli.postrun_hooks_support import (
    _allocate_scorch_jobs_threads,
    _as_int,
    _is_scorch_enabled,
    _normalize_scorch_parallel_profile,
    _resolve_hook_roots,
    _resolve_scorch_decoy_prefix,
    _resolve_scorch_parallel_profile,
    _scheduler_queue_depth,
    _scorch_cmd_prefix,
    _scorch_subprocess_env,
    _with_repo_src_on_pythonpath,
)
from config.output_paths import runtime_root
from post_docking.rescoring.scorch_device import (
    apply_scorch_device_plan_to_cfg,
    build_scorch_device_plan,
    cap_scorch_jobs_for_device,
    scorch_child_env,
)
from post_docking.rescoring.scorch_selection import load_control_bases
from post_docking.rescoring.scorch_coverage import (
    combo_done_sentinel,
    evaluate_scorch_coverage,
    remove_done_sentinel,
    resolve_top_fraction,
    write_coverage_summary,
)
from post_docking.rescoring.scorch_lock_fallback import run_timed_subprocess
from post_docking.rescoring.scorch_shards import shard_mode_enabled
from post_docking.rescoring import scorch_provisional_cache as _scorch_cache_mod


def _resolve_sem_free_slots():
    from docking.global_scheduler import sem_free_slots as impl

    return impl


def _resolve_acquire_global_cores():
    from docking.global_scheduler import acquire_global_cores as impl

    return impl


sem_free_slots = _resolve_sem_free_slots()
acquire_global_cores = _resolve_acquire_global_cores()


def _load_scorch_control_bases(
    cfg: Mapping[str, Any],
    pdb_id: str,
    logger: logging.Logger,
) -> set[str]:
    try:
        processed_root = runtime_root(
            cfg,
            "OUTPUT_DIR",
            "processed_pdbs",
            prefer_existing=True,
        )
        return set(
            load_control_bases(
                processed_root,
                str(pdb_id).upper(),
                logger,
                component="[scorch-coverage]",
            )
        )
    except Exception:
        logger.debug(
            "[scorch-coverage] action=controls status=failed pdb_id=%s",
            pdb_id,
            exc_info=True,
        )
        return set()


def _update_manifest_for_postprocessing(*args, **kwargs):
    from cli import postrun_hooks_runtime as _runtime

    return _runtime.update_manifest_for_postprocessing(*args, **kwargs)


def _float_cfg_or_env(
    cfg: Mapping[str, Any],
    *names: str,
    default: float,
) -> float:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        try:
            value = float(str(raw).strip())
        except Exception:
            continue
        if value >= 0.0:
            return float(value)
    return float(default)


def _int_cfg_or_env(
    cfg: Mapping[str, Any],
    *names: str,
    default: int,
) -> int:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        try:
            value = int(str(raw).strip())
        except Exception:
            continue
        if value >= 0:
            return int(value)
    return int(default)


def _cfg_or_env_is_set(cfg: Mapping[str, Any], *names: str) -> bool:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is not None and str(raw).strip():
            return True
    return False


def _run_scorch_hook_command(
    *,
    cmd: list[str],
    cwd: Path,
    env: Mapping[str, str],
    cfg: Mapping[str, Any],
    logger: logging.Logger,
    heartbeat_label: str,
    heartbeat_fn: Optional[Any] = None,
) -> subprocess.CompletedProcess[str]:
    timeout_sec = _float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_HOOK_TIMEOUT_SEC",
        "SCORCH_HOOK_TIMEOUT_SEC",
        default=0.0,
    )
    heartbeat_sec = _float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_HOOK_HEARTBEAT_SEC",
        "SCORCH_HOOK_HEARTBEAT_SEC",
        default=120.0,
    )
    started = time.time()

    def _heartbeat() -> None:
        elapsed = max(0.0, time.time() - started)
        logger.info(
            "[scorch-rescore.invoke] heartbeat=true scope=%s elapsed_s=%.1f timeout_s=%.1f",
            heartbeat_label,
            elapsed,
            float(timeout_sec),
        )
        if heartbeat_fn is not None:
            try:
                heartbeat_fn(elapsed)
            except Exception:
                pass

    run_kwargs: dict[str, Any] = {
        "cwd": str(cwd),
        "check": False,
        "env": dict(env),
    }
    heartbeat_enabled = heartbeat_fn is not None and not os.environ.get(
        "PYTEST_CURRENT_TEST"
    )
    if timeout_sec > 0 or heartbeat_enabled:
        return run_timed_subprocess(
            cmd,
            run_kwargs=run_kwargs,
            timeout_sec=float(timeout_sec) if timeout_sec > 0 else None,
            heartbeat_fn=_heartbeat if heartbeat_enabled else None,
            heartbeat_sec=float(max(1.0, heartbeat_sec)),
        )
    return subprocess.run(cmd, **run_kwargs)


def _scorch_scope_done_sentinel(
    post_docked_root: Path,
    *,
    run_id: str,
    pdb_id: str,
    variant: Optional[str],
    ph: Optional[str],
) -> Path:
    return combo_done_sentinel(
        Path(post_docked_root) / str(run_id),
        (str(pdb_id).upper(), str(variant or ""), str(ph or "")),
    )


def _validate_scorch_scope_coverage(
    *,
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant: Optional[str],
    ph: Optional[str],
    docked_root: Path,
    post_docked_root: Path,
    logger: logging.Logger,
    clear_stale_done: bool,
) -> bool:
    decoy_prefix = _resolve_scorch_decoy_prefix(cfg)
    summary = evaluate_scorch_coverage(
        docked_run_root=Path(docked_root) / str(run_id),
        post_run_root=Path(post_docked_root) / str(run_id),
        pdb_id=str(pdb_id),
        variant=variant or "",
        ph=ph or "",
        top_fraction=resolve_top_fraction(cfg),
        decoy_prefix=decoy_prefix,
        control_bases=_load_scorch_control_bases(cfg, str(pdb_id), logger),
        manifest_run_root=runtime_root(cfg, "MANIFESTS_DIR", "manifests") / str(run_id),
    )
    try:
        write_coverage_summary(Path(post_docked_root) / str(run_id), summary)
    except Exception:
        logger.debug(
            "[scorch-coverage] action=write_summary status=failed run_id=%s pdb_id=%s",
            run_id,
            pdb_id,
            exc_info=True,
        )
    if summary.all_complete:
        logger.info(
            "[scorch-coverage] status=complete run_id=%s pdb_id=%s variant=%s ph=%s fda=%d/%d dud=%d/%d final_fda=%d final_dud=%d",
            run_id,
            pdb_id,
            variant or "*",
            ph or "*",
            summary.fda.raw_rows,
            summary.fda.expected_rows,
            summary.dud.raw_rows,
            summary.dud.expected_rows,
            summary.fda.final_rescored_rows,
            summary.dud.final_rescored_rows,
        )
        return True

    cleared = False
    if clear_stale_done:
        try:
            cleared = remove_done_sentinel(Path(post_docked_root) / str(run_id), summary)
        except Exception:
            logger.debug(
                "[scorch-coverage] action=clear_done status=failed run_id=%s pdb_id=%s",
                run_id,
                pdb_id,
                exc_info=True,
            )
    logger.warning(
        "[scorch-coverage] status=incomplete run_id=%s pdb_id=%s variant=%s ph=%s cleared_done=%s reasons=%s fda_raw=%d/%d fda_final=%d/%d dud_raw=%d/%d dud_final=%d/%d",
        run_id,
        pdb_id,
        variant or "*",
        ph or "*",
        str(bool(cleared)).lower(),
        ",".join(summary.reasons) or "unknown",
        summary.fda.raw_rows,
        summary.fda.expected_rows,
        summary.fda.final_rescored_rows,
        summary.fda.expected_rows,
        summary.dud.raw_rows,
        summary.dud.expected_rows,
        summary.dud.final_rescored_rows,
        summary.dud.expected_rows,
    )
    return False

def _maybe_run_scorch_rescore(
    cfg: Mapping[str, Any], run_id: str, verbose: bool = False
) -> None:
    """
    Best-effort post-run SCORCH rescoring. Never raises.
    """
    logger = logging.getLogger("scorch-rescore-hook")
    if not run_id:
        logger.info("[scorch-rescore.skip] reason=missing_run_id")
        return
    if _is_no_library_docking(cfg):
        logger.info("[scorch-rescore.skip] reason=no_library_docking run_id=%s", run_id)
        return

    roots = _resolve_hook_roots(cfg)
    repo_root = roots.code_root
    docked_root = roots.docked_root
    post_docked_root = roots.post_docked_root
    script_path = repo_root / "src/post_docking/rescoring/rescoring_scorch.py"
    if not script_path.exists():
        logger.warning(
            "[scorch-rescore.skip] reason=missing_script path=%s", script_path
        )
        return

    total_cpu = _as_int(cfg.get("CPU") or (os.cpu_count() or 1), os.cpu_count() or 1)
    total_cpu = max(1, total_cpu)

    device_plan = build_scorch_device_plan(cfg, logger)
    cfg_device = apply_scorch_device_plan_to_cfg(cfg, device_plan)

    profile = _resolve_scorch_parallel_profile(cfg_device)
    jobs, threads, util_target = _allocate_scorch_jobs_threads(
        cpu_budget=total_cpu,
        free_cores=total_cpu,
        backlog=max(1, total_cpu),
        allowed_count=128,
        scheduler_queue_depth=0,
        local_queue_depth=max(0, total_cpu - 1),
        profile=profile,
    )
    jobs = cap_scorch_jobs_for_device(jobs, device_plan)
    decoy_prefix = _resolve_scorch_decoy_prefix(cfg_device)
    scorch_env = _with_repo_src_on_pythonpath(
        _scorch_subprocess_env(cfg_device), repo_root
    )
    scorch_env = scorch_child_env(cfg_device, base_env=scorch_env)

    try:
        cmd = _scorch_cmd_prefix(script_path, cfg_device) + [
            "--run-id",
            run_id,
            "--repo-root",
            str(repo_root),
            "--docked-root",
            str(docked_root),
            "--post-docked-root",
            str(post_docked_root),
            "--threads",
            str(threads),
            "--jobs",
            str(jobs),
            "--decoy-prefix",
            str(decoy_prefix),
        ]
    except Exception as exc:
        logger.warning(
            "[scorch-rescore.invoke] action=skip reason=scorch_env_unresolved run_id=%s error=%s",
            run_id,
            exc,
        )
        return
    if verbose or logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        cmd.append("--verbose")

    logger.info(
        "[scorch-rescore.invoke] run_id=%s jobs=%d threads=%d total_cpu=%d profile=%s util_target=%.2f device=%s backend=%s gpu_count=%d docked_root=%s post_docked_root=%s cmd=%s",
        run_id,
        jobs,
        threads,
        total_cpu,
        profile,
        util_target,
        device_plan.effective,
        device_plan.backend,
        device_plan.gpu_count,
        docked_root,
        post_docked_root,
        shlex.join(str(c) for c in cmd),
    )
    try:
        proc = _run_scorch_hook_command(
            cmd=cmd,
            cwd=repo_root,
            env=scorch_env,
            cfg=cfg_device,
            logger=logger,
            heartbeat_label=f"run:{run_id}",
        )
    except Exception:
        logger.warning(
            "[scorch-rescore.invoke] action=skip reason=execution_failed", exc_info=True
        )
        return

    if proc.returncode != 0:
        logger.warning(
            "[scorch-rescore.result] status=failed run_id=%s returncode=%d",
            run_id,
            proc.returncode,
        )
    else:
        logger.info("[scorch-rescore.result] status=ok run_id=%s", run_id)


def _maybe_run_scorch_rescore_for_pdb(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    *,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
    backlog_hint: Optional[int] = None,
    free_cores_hint: Optional[int] = None,
    scheduler_queue_depth_hint: Optional[int] = None,
    allowed_count_hint: Optional[int] = None,
    cpu_budget_hint: Optional[int] = None,
    profile_hint: Optional[str] = None,
    execution_mode: str = "subprocess",
    overwrite: bool = False,
    verbose: bool = False,
    provisional: bool = False,
) -> Optional[int]:
    """
    Best-effort per-protein SCORCH rescoring. Never raises.
    Returns subprocess-like return code on execution path, else None for skipped/error.
    """
    logger = logging.getLogger("scorch-rescore-hook")
    if not run_id:
        logger.info("[scorch-rescore.skip] reason=missing_run_id")
        return None
    if not pdb_id:
        logger.info("[scorch-rescore.skip] reason=missing_pdb_id run_id=%s", run_id)
        return None
    if _is_no_library_docking(cfg):
        logger.info(
            "[scorch-rescore.skip] reason=no_library_docking run_id=%s pdb_id=%s",
            run_id,
            pdb_id,
        )
        return None

    if not _is_scorch_enabled(cfg):
        logger.info(
            "[scorch-rescore.skip] reason=use_scorch_false run_id=%s pdb_id=%s",
            run_id,
            pdb_id,
        )
        return None

    roots = _resolve_hook_roots(cfg)
    repo_root = roots.code_root
    docked_root = roots.docked_root
    post_docked_root = roots.post_docked_root
    script_path = repo_root / "src/post_docking/rescoring/rescoring_scorch.py"
    if not script_path.exists():
        logger.warning(
            "[scorch-rescore.skip] reason=missing_script path=%s", script_path
        )
        return None

    dist_ctx = None
    final_scope_key: Optional[str] = None
    scorch_scope_claimed = False
    distributed_shard_scope = False
    if not provisional:
        done_sentinel = _scorch_scope_done_sentinel(
            post_docked_root,
            run_id=str(run_id),
            pdb_id=str(pdb_id),
            variant=variant,
            ph=ph,
        )
        if not overwrite and done_sentinel.exists():
            if _validate_scorch_scope_coverage(
                cfg=cfg,
                run_id=str(run_id),
                pdb_id=str(pdb_id),
                variant=variant,
                ph=ph,
                docked_root=docked_root,
                post_docked_root=post_docked_root,
                logger=logger,
                clear_stale_done=True,
            ):
                logger.info(
                    "[scorch-rescore.skip] reason=scope_done_coverage_ok run_id=%s pdb_id=%s variant=%s ph=%s done=%s",
                    run_id,
                    pdb_id,
                    variant or "*",
                    ph or "*",
                    done_sentinel,
                )
                return None
            logger.warning(
                "[scorch-rescore.retry] reason=scope_done_coverage_failed run_id=%s pdb_id=%s variant=%s ph=%s done=%s",
                run_id,
                pdb_id,
                variant or "*",
                ph or "*",
                done_sentinel,
            )
        try:
            dist_ctx = resolve_distributed_context(cfg, str(run_id))
        except Exception:
            dist_ctx = None
        if dist_ctx is not None and bool(getattr(dist_ctx, "enabled", False)):
            distributed_shard_scope = bool(shard_mode_enabled(cfg))
            if distributed_shard_scope:
                logger.info(
                    "[scorch-rescore.scope-claim] action=skip reason=distributed_shard_mode run_id=%s pdb_id=%s variant=%s ph=%s",
                    run_id,
                    pdb_id,
                    variant or "*",
                    ph or "*",
                )
            else:
                final_scope_key = scorch_scope_key(
                    run_id=str(run_id),
                    pdb_id=str(pdb_id),
                    variant=variant,
                    ph=ph,
                    phase="final",
                )
                lease_sec = float(
                    max(
                        30,
                        _int_cfg_or_env(
                            cfg,
                            "ATLAS_SCORCH_SCOPE_CLAIM_LEASE_SEC",
                            "SCORCH_SCOPE_CLAIM_LEASE_SEC",
                            default=21600,
                        ),
                    )
                )
                try:
                    scorch_scope_claimed = bool(
                        try_claim_scorch_scope(
                            dist_ctx,
                            cfg,
                            scope_key=str(final_scope_key),
                            pdb_id=str(pdb_id),
                            variant=variant,
                            ph=ph,
                            phase="final",
                            lease_sec=float(lease_sec),
                            force=bool(overwrite),
                        )
                    )
                except Exception:
                    scorch_scope_claimed = False
                    logger.warning(
                        "[scorch-rescore.scope-claim] action=error run_id=%s pdb_id=%s variant=%s ph=%s",
                        run_id,
                        pdb_id,
                        variant or "*",
                        ph or "*",
                        exc_info=True,
                    )
                if not scorch_scope_claimed:
                    logger.info(
                        "[scorch-rescore.skip] reason=scope_claim_unavailable run_id=%s pdb_id=%s variant=%s ph=%s scope_key=%s",
                        run_id,
                        pdb_id,
                        variant or "*",
                        ph or "*",
                        final_scope_key,
                    )
                    return None

    cfg_cpu = max(
        1,
        _as_int(cfg.get("CPU") or (os.cpu_count() or 1), os.cpu_count() or 1),
    )
    hinted_cpu = cfg_cpu
    if cpu_budget_hint is not None:
        try:
            hinted_cpu = min(cfg_cpu, max(1, int(cpu_budget_hint)))
        except Exception:
            hinted_cpu = cfg_cpu
    total_cpu = max(1, int(hinted_cpu))
    jobs = total_cpu
    threads = 1
    device_plan = build_scorch_device_plan(cfg, logger)
    cfg_device = apply_scorch_device_plan_to_cfg(cfg, device_plan)
    if provisional:
        cfg_device = dict(cfg_device)
        cfg_device["SCORCH_PROVISIONAL_ENABLE"] = True
        cfg_device["SCORCH_PROVISIONAL_CACHE_WRITE"] = True
        cfg_device["SCORCH_PROVISIONAL_REUSE"] = False
        provisional_cpu_cap = max(
            1,
            _as_int(
                cfg_device.get("SCORCH_PROVISIONAL_CPU_BUDGET"),
                min(2, max(1, int(total_cpu))),
            ),
        )
        total_cpu = max(1, min(int(total_cpu), int(provisional_cpu_cap)))
        cfg_device["CPU"] = int(total_cpu)
        final_fraction = 0.10
        try:
            final_fraction = float(cfg_device.get("SCORCH_TOP_FRACTION", 0.10) or 0.10)
        except Exception:
            final_fraction = 0.10
        provisional_fraction_raw = cfg_device.get("SCORCH_PROVISIONAL_TOP_FRACTION")
        if provisional_fraction_raw is None:
            provisional_fraction = max(final_fraction, min(1.0, final_fraction * 2.0))
        else:
            try:
                provisional_fraction = float(provisional_fraction_raw)
            except Exception:
                provisional_fraction = max(final_fraction, min(1.0, final_fraction * 2.0))
        if not (0.0 < float(provisional_fraction) <= 1.0):
            provisional_fraction = max(final_fraction, min(1.0, final_fraction * 2.0))
        cfg_device["SCORCH_TOP_FRACTION"] = float(provisional_fraction)
        _scorch_cache_mod.artifact_post_root(cfg_device, run_id).mkdir(
            parents=True,
            exist_ok=True,
        )
    decoy_prefix = _resolve_scorch_decoy_prefix(cfg_device)
    scorch_env = _with_repo_src_on_pythonpath(
        _scorch_subprocess_env(cfg_device), repo_root
    )
    scorch_env = scorch_child_env(cfg_device, base_env=scorch_env)
    scheduler = None
    if isinstance(cfg_device, Mapping):
        scheduler = cfg_device.get("GLOBAL_DOCKING_SCHEDULER")
    free_now = (
        max(1, int(free_cores_hint))
        if free_cores_hint is not None
        else (
            max(1, int(sem_free_slots(scheduler, total_cpu)))
            if scheduler is not None
            else total_cpu
        )
    )
    backlog_raw = max(
        1,
        int(backlog_hint) if backlog_hint is not None else total_cpu,
    )
    scheduler_queue_depth = (
        max(0, int(scheduler_queue_depth_hint))
        if scheduler_queue_depth_hint is not None
        else _scheduler_queue_depth(scheduler)
    )
    local_queue_depth = max(0, backlog_raw - 1)
    if backlog_hint is None:
        local_queue_depth = 0
    profile = (
        _normalize_scorch_parallel_profile(profile_hint)
        if profile_hint is not None
        else _resolve_scorch_parallel_profile(cfg_device)
    )
    allowed_count = max(
        1,
        int(allowed_count_hint)
        if allowed_count_hint is not None
        else int(cfg_device.get("_SCORCH_ALLOWED_COUNT", 128) or 128),
    )
    # Use ligand cardinality to keep per-combo SCORCH work saturated even when
    # queue backlog is small during the distributed tail.
    demand_from_allowed = max(1, int(math.ceil(float(allowed_count) / 64.0)))
    backlog = max(int(backlog_raw), int(demand_from_allowed))
    jobs, threads, util_target = _allocate_scorch_jobs_threads(
        cpu_budget=total_cpu,
        free_cores=free_now,
        backlog=backlog,
        allowed_count=allowed_count,
        scheduler_queue_depth=scheduler_queue_depth,
        local_queue_depth=local_queue_depth,
        profile=profile,
    )
    jobs = cap_scorch_jobs_for_device(jobs, device_plan)
    if shard_mode_enabled(cfg_device):
        chunk_min = max(
            1,
            _int_cfg_or_env(
                cfg_device,
                "ATLAS_SCORCH_CHUNK_MIN",
                "SCORCH_CHUNK_MIN",
                default=16,
            ),
        )
        max_productive_shards = max(
            1,
            int(math.ceil(float(max(1, int(allowed_count))) / float(chunk_min))),
        )
        desired_util_cores_local = max(
            1,
            int(math.floor(float(util_target) * float(max(1, int(free_now))))),
        )
        user_workers_cap = _cfg_or_env_is_set(
            cfg_device,
            "ATLAS_SCORCH_SHARD_LOCAL_WORKERS",
            "SCORCH_SHARD_LOCAL_WORKERS",
        )
        user_threads_cap = _cfg_or_env_is_set(
            cfg_device,
            "ATLAS_SCORCH_SHARD_THREADS",
            "SCORCH_SHARD_THREADS",
        )
        if device_plan.effective != "gpu":
            default_thread_cap = 1
            productive_cpu_cap = max(1, min(max_productive_shards, int(free_now)))
            # CPU SCORCH subprocesses are memory/feature-prep heavy.  Treat the
            # queue-provided CPU budget as a hard lane budget and use at most
            # about half of it for simultaneous SCORCH subprocesses; this keeps
            # distributed resume runs from pinning a full node with Python/TF
            # workers that then time out in feature preparation.
            budget_process_cap = max(1, int(total_cpu) // 2)
            node_floor = max(8, max(1, int(total_cpu) // 4))
            target_worker_cap = max(int(desired_util_cores_local), int(node_floor))
            default_worker_cap = max(
                1,
                min(
                    max_productive_shards,
                    productive_cpu_cap,
                    target_worker_cap,
                    budget_process_cap,
                ),
            )
        elif profile == "conservative" or provisional:
            default_thread_cap = min(6, max(1, int(threads)))
            default_worker_cap = max(
                1,
                min(max_productive_shards, max(2, total_cpu // 16 + 1)),
            )
        else:
            default_thread_cap = min(16, max(4, max(1, int(free_now)) // 8))
            default_worker_cap = max(1, min(max_productive_shards, max(4, total_cpu // 4)))
        shard_workers = _int_cfg_or_env(
            cfg_device,
            "ATLAS_SCORCH_SHARD_LOCAL_WORKERS",
            "SCORCH_SHARD_LOCAL_WORKERS",
            default=default_worker_cap,
        )
        shard_threads = _int_cfg_or_env(
            cfg_device,
            "ATLAS_SCORCH_SHARD_THREADS",
            "SCORCH_SHARD_THREADS",
            default=default_thread_cap,
        )
        thread_cap = max(1, int(shard_threads))
        worker_cap = max(1, int(shard_workers))
        if user_threads_cap:
            threads = max(1, min(int(threads), thread_cap))
        else:
            threads = max(1, min(int(threads), thread_cap))
        if user_workers_cap:
            jobs = max(1, min(int(jobs), worker_cap, max_productive_shards))
        else:
            jobs = max(1, min(int(jobs), worker_cap, max_productive_shards))
            if device_plan.effective != "gpu":
                jobs = max(
                    int(jobs),
                    min(int(worker_cap), int(max_productive_shards)),
                )
        jobs = max(1, int(jobs))
        threads = max(1, int(threads))
        if not _cfg_or_env_is_set(
            cfg_device,
            "ATLAS_SCORCH_HOOK_TIMEOUT_SEC",
            "SCORCH_HOOK_TIMEOUT_SEC",
        ):
            cfg_device = dict(cfg_device)
            cfg_device["ATLAS_SCORCH_HOOK_TIMEOUT_SEC"] = _float_cfg_or_env(
                cfg_device,
                "ATLAS_SCORCH_DISTRIBUTED_HOOK_TIMEOUT_SEC",
                "SCORCH_DISTRIBUTED_HOOK_TIMEOUT_SEC",
                default=5400.0,
            )

    module_args = [
        "--run-id",
        str(run_id),
        "--repo-root",
        str(repo_root),
        "--docked-root",
        str(docked_root),
        "--post-docked-root",
        str(post_docked_root),
        "--pdb-id",
        str(pdb_id),
        "--threads",
        str(threads),
        "--jobs",
        str(jobs),
        "--decoy-prefix",
        str(decoy_prefix),
    ]
    if variant:
        module_args.extend(["--variant", str(variant)])
    if ph:
        module_args.extend(["--ph", str(ph)])
    if overwrite:
        module_args.append("--overwrite")
    if provisional:
        module_args.append("--skip-autofix")
    if verbose or logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        module_args.append("--verbose")
    try:
        cmd = _scorch_cmd_prefix(script_path, cfg_device) + module_args
    except Exception as exc:
        logger.warning(
            "[scorch-rescore.invoke] action=skip reason=scorch_env_unresolved run_id=%s pdb_id=%s error=%s",
            run_id,
            pdb_id,
            exc,
        )
        if scorch_scope_claimed and dist_ctx is not None and final_scope_key is not None:
            try:
                write_scorch_scope_state(
                    dist_ctx,
                    cfg_device,
                    scope_key=str(final_scope_key),
                    status="failed",
                    payload={
                        "pdb_id": str(pdb_id).upper(),
                        "variant": str(variant or "*").upper(),
                        "ph": str(ph or "base"),
                        "phase": "final",
                        "error": "scorch_env_unresolved",
                    },
                    release_claim=True,
                )
            except Exception:
                pass
        if not provisional and not distributed_shard_scope:
            _update_manifest_for_postprocessing(
                cfg_device,
                run_id,
                str(pdb_id),
                variant,
                ph_tag=ph,
                status="failed",
                error="scorch_env_unresolved",
                details={"scorch_phase": "final", "scorch_env_error": str(exc)},
            )
        return None
    argv = list(module_args)
    desired_util_cores = max(1, int(math.floor(float(util_target) * float(free_now))))
    allocated_util_cores = max(1, int(jobs) * int(threads))

    logger.info(
        "[scorch-rescore.invoke] scope=per_pdb phase=%s run_id=%s pdb_id=%s variant=%s ph=%s jobs=%d threads=%d total_cpu=%d cfg_cpu=%d free_cores=%d backlog=%d backlog_effective=%d allowed_count=%d scheduler_q=%d profile=%s util_target=%.2f desired_util_cores=%d allocated_util_cores=%d scheduler_shared=%s mode=%s device=%s backend=%s gpu_count=%d docked_root=%s post_docked_root=%s cmd=%s",
        "provisional" if provisional else "final",
        run_id,
        pdb_id,
        variant or "*",
        ph or "*",
        jobs,
        threads,
        total_cpu,
        cfg_cpu,
        free_now,
        backlog_raw,
        backlog,
        allowed_count,
        scheduler_queue_depth,
        profile,
        util_target,
        desired_util_cores,
        allocated_util_cores,
        bool(scheduler is not None),
        str(execution_mode),
        device_plan.effective,
        device_plan.backend,
        device_plan.gpu_count,
        docked_root,
        post_docked_root,
        shlex.join(str(c) for c in cmd),
    )
    scorch_started = time.time()
    if not provisional and not distributed_shard_scope:
        _update_manifest_for_postprocessing(
            cfg_device,
            run_id,
            str(pdb_id),
            variant,
            ph_tag=ph,
            status="running",
            details={
                "scorch_phase": "final",
                "scorch_jobs": int(jobs),
                "scorch_threads": int(threads),
                "scorch_profile": str(profile),
                "scorch_util_target": float(util_target),
                "scorch_cpu_budget": int(total_cpu),
                "scorch_scheduler_shared": bool(scheduler is not None),
                **device_plan.manifest_details(),
            },
        )
    try:
        mode_token = str(execution_mode or "subprocess").strip().lower()
        if mode_token == "inprocess":
            from post_docking.rescoring import rescoring_scorch as _scorch_mod

            try:
                rc = int(_scorch_mod.main(argv=argv, cfg_override=dict(cfg_device)))
            except SystemExit as exc:
                rc = int(exc.code or 0)
        else:
            def _manifest_heartbeat(elapsed: float) -> None:
                if (
                    scorch_scope_claimed
                    and dist_ctx is not None
                    and final_scope_key is not None
                ):
                    try:
                        renew_scorch_scope_claim(
                            dist_ctx,
                            cfg_device,
                            scope_key=str(final_scope_key),
                        )
                    except Exception:
                        pass
                if provisional or distributed_shard_scope:
                    return
                _update_manifest_for_postprocessing(
                    cfg_device,
                    run_id,
                    str(pdb_id),
                    variant,
                    ph_tag=ph,
                    status="running",
                    details={
                        "scorch_phase": "final",
                        "scorch_mode": mode_token,
                        "scorch_last_heartbeat_elapsed_sec": float(elapsed),
                    },
                )

            proc = _run_scorch_hook_command(
                cmd=cmd,
                cwd=repo_root,
                env=scorch_env,
                cfg=cfg_device,
                logger=logger,
                heartbeat_label=f"{run_id}:{pdb_id}:{variant or '*'}:{ph or '*'}",
                heartbeat_fn=_manifest_heartbeat,
            )
            rc = int(proc.returncode)
    except Exception:
        logger.warning(
            "[scorch-rescore.invoke] action=skip reason=execution_failed run_id=%s pdb_id=%s",
            run_id,
            pdb_id,
            exc_info=True,
        )
        if scorch_scope_claimed and dist_ctx is not None and final_scope_key is not None:
            try:
                write_scorch_scope_state(
                    dist_ctx,
                    cfg_device,
                    scope_key=str(final_scope_key),
                    status="failed",
                    payload={
                        "pdb_id": str(pdb_id).upper(),
                        "variant": str(variant or "*").upper(),
                        "ph": str(ph or "base"),
                        "phase": "final",
                        "error": f"{mode_token}_execution_failed",
                    },
                    release_claim=True,
                )
            except Exception:
                pass
        if not provisional and not distributed_shard_scope:
            _update_manifest_for_postprocessing(
                cfg_device,
                run_id,
                str(pdb_id),
                variant,
                ph_tag=ph,
                status="failed",
                elapsed_sec=time.time() - scorch_started,
                error=f"{mode_token}_execution_failed",
                details={"scorch_mode": mode_token},
            )
        return None

    logger.info(
        "[scorch-rescore.result] scope=per_pdb phase=%s run_id=%s pdb_id=%s variant=%s ph=%s returncode=%d mode=%s",
        "provisional" if provisional else "final",
        run_id,
        pdb_id,
        variant or "*",
        ph or "*",
        rc,
        str(mode_token),
    )
    coverage_ok = True
    if not provisional and int(rc) == 0:
        coverage_ok = _validate_scorch_scope_coverage(
            cfg=cfg_device,
            run_id=str(run_id),
            pdb_id=str(pdb_id),
            variant=variant,
            ph=ph,
            docked_root=docked_root,
            post_docked_root=post_docked_root,
            logger=logger,
            clear_stale_done=True,
        )
        if not coverage_ok:
            rc = 2
    if not provisional:
        if scorch_scope_claimed and dist_ctx is not None and final_scope_key is not None:
            try:
                write_scorch_scope_state(
                    dist_ctx,
                    cfg_device,
                    scope_key=str(final_scope_key),
                    status="completed" if int(rc) == 0 else "failed",
                    payload={
                        "pdb_id": str(pdb_id).upper(),
                        "variant": str(variant or "*").upper(),
                        "ph": str(ph or "base"),
                        "phase": "final",
                        "returncode": int(rc),
                        "elapsed_sec": time.time() - scorch_started,
                    },
                    release_claim=True,
                )
            except Exception:
                pass
        if not distributed_shard_scope:
            _update_manifest_for_postprocessing(
                cfg_device,
                run_id,
                str(pdb_id),
                variant,
                ph_tag=ph,
                status="completed" if int(rc) == 0 else "failed",
                elapsed_sec=time.time() - scorch_started,
                error=None
                if int(rc) == 0
                else (
                    "scorch_coverage_incomplete"
                    if not coverage_ok
                    else f"returncode={int(rc)}"
                ),
                details={
                    "scorch_phase": "final",
                    "scorch_mode": str(mode_token),
                    "scorch_returncode": int(rc),
                    "scorch_coverage_ok": bool(coverage_ok),
                    **device_plan.manifest_details(),
                },
            )
    return int(rc)


def _run_provisional_scorch_rescore_for_pdb(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    *,
    provisional_cpu: int,
    variant: Optional[str],
    ph: Optional[str],
    allowed_count_hint: Optional[int],
    execution_mode: str,
    verbose: bool,
) -> Optional[int]:
    task_id = f"scorch:provisional:{run_id}:{pdb_id}:{variant or '*'}:{ph or '*'}"
    with acquire_global_cores(
        dict(cfg),
        cores=max(1, int(provisional_cpu)),
        min_cores=1,
        priority=-1,
        est_duration_sec=10.0,
        task_id=task_id,
        task_type="scorch_provisional",
    ) as granted:
        return _maybe_run_scorch_rescore_for_pdb(
            cfg,
            run_id,
            pdb_id,
            variant=variant,
            ph=ph,
            backlog_hint=1,
            free_cores_hint=int(granted),
            scheduler_queue_depth_hint=0,
            allowed_count_hint=allowed_count_hint,
            cpu_budget_hint=int(granted),
            profile_hint="conservative",
            execution_mode=execution_mode,
            overwrite=False,
            verbose=verbose,
            provisional=True,
        )


def _log_rescore_verification(
    run_id: str, cfg: Optional[Mapping[str, Any]] = None
) -> None:
    """
    Lightweight best-effort check for consensus and SCORCH outputs.
    """
    logger = logging.getLogger("post-run-checks")
    if not run_id:
        logger.info("[post-check.skip] reason=missing_run_id")
        return
    if _is_no_library_docking(cfg or {}):
        logger.info("[post-check.skip] reason=no_library_docking run_id=%s", run_id)
        return

    roots = _resolve_hook_roots(cfg or {})
    docked_base = roots.docked_root
    post_base = roots.post_docked_root

    docked_root = docked_base / run_id
    post_root = post_base / run_id

    def _count_and_sample(root: Path, pattern: str) -> tuple[int, Optional[Path]]:
        count = 0
        sample: Optional[Path] = None
        if not root.exists():
            return count, sample
        for path in root.rglob(pattern):
            if not path.is_file():
                continue
            count += 1
            if sample is None:
                sample = path
        return count, sample

    consensus_count, consensus_sample = _count_and_sample(
        docked_root, "consensus_docking_scores.csv"
    )
    scorch_count, scorch_sample = _count_and_sample(post_root, "scorch_scores_all.csv")
    reranked_count, reranked_sample = _count_and_sample(
        post_root, "consensus_reranked_scorch.csv"
    )

    if consensus_sample is not None:
        logger.info(
            "[post-check.consensus] found=%d sample=%s",
            consensus_count,
            consensus_sample,
        )
    else:
        logger.warning("[post-check.consensus] reason=missing_files run_id=%s", run_id)

    if scorch_sample is not None:
        logger.info(
            "[post-check.scorch] found=%d sample=%s", scorch_count, scorch_sample
        )
    else:
        logger.warning("[post-check.scorch] reason=missing_files run_id=%s", run_id)

    if reranked_sample is not None:
        logger.info(
            "[post-check.scorch-reranked] found=%d sample=%s",
            reranked_count,
            reranked_sample,
        )
    else:
        logger.warning(
            "[post-check.scorch-reranked] reason=missing_files run_id=%s", run_id
        )

    try:
        decoy_prefix = _resolve_scorch_decoy_prefix(cfg or {})
        top_fraction = resolve_top_fraction(cfg or {})
        incomplete: list[str] = []
        audited = 0
        for consensus_csv in sorted(docked_root.rglob("consensus_docking_scores.csv")):
            rel_parts = consensus_csv.parent.relative_to(docked_root).parts
            if not rel_parts:
                continue
            pdb_id = rel_parts[0]
            variant = rel_parts[1] if len(rel_parts) >= 2 else ""
            ph = rel_parts[2] if len(rel_parts) >= 3 else ""
            audited += 1
            summary = evaluate_scorch_coverage(
                docked_run_root=docked_root,
                post_run_root=post_root,
                pdb_id=pdb_id,
                variant=variant,
                ph=ph,
                top_fraction=top_fraction,
                decoy_prefix=decoy_prefix,
                control_bases=_load_scorch_control_bases(cfg or {}, pdb_id, logger),
                manifest_run_root=runtime_root(
                    cfg or {}, "MANIFESTS_DIR", "manifests"
                )
                / str(run_id),
            )
            if not summary.all_complete:
                incomplete.append(f"{pdb_id}:{'|'.join(summary.reasons) or 'unknown'}")
        if incomplete:
            logger.warning(
                "[post-check.scorch-coverage] status=incomplete run_id=%s audited=%d incomplete=%d examples=%s",
                run_id,
                audited,
                len(incomplete),
                ",".join(incomplete[:5]),
            )
        elif audited:
            logger.info(
                "[post-check.scorch-coverage] status=complete run_id=%s audited=%d",
                run_id,
                audited,
            )
    except Exception:
        logger.warning(
            "[post-check.scorch-coverage] status=error run_id=%s",
            run_id,
            exc_info=True,
        )

    try:
        from cli.run_manifest_runtime import load_run_manifest as _load_run_manifest

        manifest = _load_run_manifest(dict(cfg or {}), run_id) or {}
    except Exception:
        logger.debug(
            "[post-check.scorch-failures] reason=manifest_unavailable run_id=%s",
            run_id,
            exc_info=True,
        )
        manifest = {}
    proteins = manifest.get("proteins") if isinstance(manifest, Mapping) else None
    failed_examples: list[str] = []
    if isinstance(proteins, Mapping):
        for key, entry in proteins.items():
            if not isinstance(entry, Mapping):
                continue
            stages = entry.get("stages")
            post = stages.get("postprocessing") if isinstance(stages, Mapping) else None
            if not isinstance(post, Mapping):
                continue
            if str(post.get("status") or "").lower() != "failed":
                continue
            details = post.get("details")
            details = details if isinstance(details, Mapping) else {}
            failed_examples.append(str(key))
            if len(failed_examples) <= 5:
                logger.warning(
                    "[post-check.scorch-failure] combo=%s error=%s returncode=%s mode=%s phase=%s device=%s",
                    key,
                    post.get("error"),
                    details.get("scorch_returncode"),
                    details.get("scorch_mode"),
                    details.get("scorch_phase"),
                    details.get("scorch_device_effective"),
                )
    if failed_examples:
        logger.warning(
            "[post-check.scorch-failures] failed=%d examples=%s",
            len(failed_examples),
            ",".join(failed_examples[:5]),
        )
