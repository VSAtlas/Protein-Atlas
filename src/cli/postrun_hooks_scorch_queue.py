# -*- coding: utf-8 -*-
from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
import logging
import os
from pathlib import Path
import threading
import time
from typing import Any, Dict, Mapping, Optional

from cli.postrun_hooks_common import _is_no_library_docking
from cli.postrun_hooks_scorch import (
    _maybe_run_scorch_rescore_for_pdb,
    _run_provisional_scorch_rescore_for_pdb,
)
from cli.postrun_hooks_support import (
    _as_int,
    _is_scorch_enabled,
    _resolve_hook_roots,
    _resolve_scorch_decoy_prefix,
    _resolve_scorch_execution_mode,
    _resolve_scorch_parallel_profile,
    _scheduler_queue_depth,
)
from post_docking.rescoring.scorch_device import (
    build_scorch_device_plan,
    cap_scorch_queue_workers,
)
from post_docking.rescoring.scorch_shards import shard_mode_enabled
from post_docking.rescoring import scorch_provisional_cache as _scorch_cache_mod
from post_docking.rescoring.scorch_coverage import find_combo_dir, normalize_combo


def _resolve_sem_free_slots():
    from docking.global_scheduler import sem_free_slots as impl

    return impl


sem_free_slots = _resolve_sem_free_slots()


def _csv_has_data_row(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size <= 0:
        return False
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            # Header plus at least one data row.
            next(handle, None)
            for line in handle:
                if line.strip():
                    return True
    except OSError:
        return False
    return False


def _consensus_ready_for_scorch(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    *,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
) -> tuple[bool, Optional[str]]:
    roots = _resolve_hook_roots(cfg)
    docked_run_root = roots.docked_root / str(run_id)
    decoy_prefix = _resolve_scorch_decoy_prefix(cfg)
    consensus_names = (
        "consensus_docking_scores.csv",
        f"{decoy_prefix}_consensus_docking_scores.csv",
        "dud_consensus_docking_scores.csv",
        "decoys_consensus_docking_scores.csv",
    )
    dock_dir = find_combo_dir(
        docked_run_root,
        normalize_combo(pdb_id, variant, ph),
        consensus_names,
    )
    if dock_dir is None:
        return False, None
    for name in consensus_names:
        candidate = dock_dir / name
        if _csv_has_data_row(candidate):
            return True, str(candidate)
    for candidate in sorted(dock_dir.glob("*_consensus_docking_scores.csv")):
        if candidate.name == "consensus_docking_scores.csv":
            continue
        if _csv_has_data_row(candidate):
            return True, str(candidate)
    return False, str(dock_dir)


class ScorchMixedQueueService:
    """
    Persistent in-process SCORCH queue to avoid per-protein process launch gaps.
    """

    def __init__(
        self,
        *,
        run_id: str,
        worker_cap: int,
        total_cpu: int,
        execution_mode: str,
        logger: logging.Logger,
        verbose: bool,
        drain_timeout_sec: float = 0.0,
    ) -> None:
        self._run_id = str(run_id)
        self._total_cpu = max(1, int(total_cpu))
        self._worker_cap = max(1, int(worker_cap))
        self._execution_mode = str(execution_mode or "subprocess")
        self._drain_timeout_sec = max(0.0, float(drain_timeout_sec or 0.0))
        self._target_workers = 1
        self._peak_target_workers = 1
        self._logger = logger
        self._verbose = bool(verbose)
        self._lock = threading.Lock()
        self._closed = False
        self._pool = ThreadPoolExecutor(
            max_workers=self._worker_cap,
            thread_name_prefix="scorch-mixed-q",
        )
        self._futures: Dict[Future[Optional[int]], Dict[str, Any]] = {}
        self._inflight_keys: set[tuple[str, str, str, str]] = set()
        self._submitted = 0
        self._failed = 0
        self._resized = 0
        self._started_at = time.time()
        self._busy_seconds = 0.0

    def _release_scope_key_locked(self, meta: Mapping[str, Any]) -> None:
        scope_key = meta.get("scope_key")
        if (
            isinstance(scope_key, tuple)
            and len(scope_key) == 4
            and all(isinstance(part, str) for part in scope_key)
        ):
            self._inflight_keys.discard(
                (
                    str(scope_key[0]),
                    str(scope_key[1]),
                    str(scope_key[2]),
                    str(scope_key[3]),
                )
            )

    def _record_done_future_locked(
        self,
        future: Future[Optional[int]],
        meta: Mapping[str, Any],
    ) -> None:
        pdb_id = str(meta.get("pdb_id", ""))
        variant = str(meta.get("variant", "*"))
        ph = str(meta.get("ph", "*"))
        self._release_scope_key_locked(meta)
        try:
            rc = future.result()
        except Exception:
            self._failed += 1
            self._logger.warning(
                "[scorch-rescore.queue.result] run_id=%s pdb_id=%s variant=%s ph=%s returncode=exception",
                self._run_id,
                pdb_id,
                variant,
                ph,
                exc_info=True,
            )
            return
        try:
            submitted_at = float(meta.get("submitted_at") or 0.0)
        except Exception:
            submitted_at = 0.0
        if submitted_at > 0.0:
            self._busy_seconds += max(0.0, time.time() - submitted_at)
        if rc is None or int(rc) != 0:
            self._failed += 1
            self._logger.warning(
                "[scorch-rescore.queue.result] run_id=%s pdb_id=%s variant=%s ph=%s returncode=%s",
                self._run_id,
                pdb_id,
                variant,
                ph,
                rc,
            )

    def _reap_done_locked(self) -> int:
        done = [future for future in self._futures if future.done()]
        for future in done:
            meta = self._futures.pop(future, {})
            self._record_done_future_locked(future, meta)
        return int(len(done))

    def _pending_count_locked(self) -> int:
        self._reap_done_locked()
        return sum(1 for fut in self._futures if not fut.done())

    def _provisional_pending_count_locked(self) -> int:
        self._reap_done_locked()
        return sum(
            1
            for fut, meta in self._futures.items()
            if not fut.done() and str(meta.get("phase", "")) == "provisional"
        )

    def _scope_active_locked(
        self,
        *,
        pdb_id: str,
        variant: Optional[str],
        ph: Optional[str],
        phase_prefix: Optional[str] = None,
    ) -> bool:
        scope_tail = (
            str(pdb_id or "").strip().upper(),
            str(variant or "*").strip().upper() or "*",
            str(ph or "base").strip() or "base",
        )
        for key in self._inflight_keys:
            if len(key) != 4 or tuple(key[1:]) != scope_tail:
                continue
            if phase_prefix is None or str(key[0]).startswith(phase_prefix):
                return True
        return False

    def _min_cpu_per_worker(self, cfg: Mapping[str, Any]) -> int:
        default_min = 24
        if shard_mode_enabled(cfg):
            default_min = max(8, min(24, max(1, int(self._total_cpu) // 8)))
        raw = os.environ.get("ATLAS_SCORCH_QUEUE_MIN_CPU_PER_WORKER")
        if raw is None:
            raw = cfg.get("SCORCH_QUEUE_MIN_CPU_PER_WORKER")
        return max(1, min(self._total_cpu, _as_int(raw, default_min)))

    def _desired_workers_locked(self, cfg: Mapping[str, Any], backlog: int) -> int:
        scheduler = cfg.get("GLOBAL_DOCKING_SCHEDULER")
        free_cores = (
            sem_free_slots(scheduler, self._total_cpu)
            if scheduler is not None
            else self._total_cpu
        )
        demand = max(1, int(backlog))
        profile = _resolve_scorch_parallel_profile(cfg)
        min_cpu_per_worker = self._min_cpu_per_worker(cfg)
        max_workers_by_cpu = max(1, int(free_cores) // int(min_cpu_per_worker))
        target = 1
        if free_cores >= max(4, self._total_cpu // 6):
            target += 1
        if free_cores >= max(8, self._total_cpu // 3):
            target += 1
        if free_cores >= max(12, self._total_cpu // 2):
            target += 1
        if free_cores >= max(16, (self._total_cpu * 3) // 4):
            target += 1
        if demand >= 4:
            target += 1
        if demand >= 8:
            target += 1
        if demand >= 12:
            target += 1
        if demand >= 24:
            target += 1
        if profile != "conservative" and free_cores >= max(16, self._total_cpu // 2):
            target = max(target, min(self._worker_cap, max(2, int(free_cores) // 24)))
        if profile != "conservative" and demand > 1:
            target = max(target, min(self._worker_cap, demand))
        return max(1, min(self._worker_cap, int(max_workers_by_cpu), int(target)))

    def _cpu_budget_hints_locked(
        self,
        *,
        cfg: Mapping[str, Any],
        free_now: int,
    ) -> tuple[int, int, int]:
        active_width = max(
            1,
            min(
                int(self._worker_cap),
                int(self._target_workers),
            ),
        )
        if shard_mode_enabled(cfg):
            # Sharded SCORCH wrappers each fan out to local subprocess workers.
            # Budget against the work that is actually pending on this process,
            # not the theoretical queue cap.  Resume tails often have only one
            # or two real targets per node; dividing by the cap leaves most of
            # an exclusive node idle even though each target can consume many
            # shard workers.
            active_width = max(
                1,
                min(
                    int(self._worker_cap),
                    int(self._pending_count_locked()) + 1,
                ),
            )
        cpu_budget = max(1, int(self._total_cpu) // int(active_width))
        free_budget = max(1, int(max(1, free_now)) // int(active_width))
        return int(cpu_budget), int(free_budget), int(active_width)

    def _refresh_worker_target_locked(self, cfg: Mapping[str, Any]) -> None:
        self._reap_done_locked()
        backlog = self._pending_count_locked() + 1
        desired = self._desired_workers_locked(cfg, backlog)
        if desired == int(self._target_workers):
            return
        self._target_workers = int(desired)
        self._peak_target_workers = max(self._peak_target_workers, int(desired))
        self._resized += 1
        # ThreadPoolExecutor does not expose a public resize API; this runtime tune
        # controls how many workers may run concurrently as backlog/fit pressure shifts.
        try:
            setattr(self._pool, "_max_workers", int(desired))
        except Exception:
            pass
        self._logger.info(
            "[scorch-rescore.queue.scale] run_id=%s workers=%d cap=%d total_cpu=%d backlog=%d",
            self._run_id,
            desired,
            self._worker_cap,
            self._total_cpu,
            backlog,
        )

    def submit(
        self,
        *,
        cfg: Mapping[str, Any],
        pdb_id: str,
        variant: Optional[str] = None,
        ph: Optional[str] = None,
        allowed_count_hint: Optional[int] = None,
    ) -> bool:
        with self._lock:
            self._reap_done_locked()
            if self._closed:
                return False
            scope_key = (
                "final",
                str(pdb_id or "").strip().upper(),
                str(variant or "*").strip().upper() or "*",
                str(ph or "base").strip() or "base",
            )
            if scope_key in self._inflight_keys:
                self._logger.debug(
                    "[scorch-rescore.queue] action=dedupe_inflight run_id=%s pdb_id=%s variant=%s ph=%s",
                    self._run_id,
                    scope_key[1],
                    scope_key[2],
                    scope_key[3],
                )
                return True
            if self._scope_active_locked(
                pdb_id=str(pdb_id),
                variant=variant,
                ph=ph,
                phase_prefix="provisional:",
            ):
                self._logger.info(
                    "[scorch-rescore.queue] action=final_suppresses_provisional run_id=%s pdb_id=%s variant=%s ph=%s",
                    self._run_id,
                    str(pdb_id),
                    str(variant or "*"),
                    str(ph or "*"),
                )
            consensus_ready, consensus_source = _consensus_ready_for_scorch(
                cfg,
                self._run_id,
                str(pdb_id),
                variant=variant,
                ph=ph,
            )
            if not consensus_ready and allowed_count_hint is None:
                self._logger.info(
                    "[scorch-rescore.queue.skip] reason=missing_docking_consensus run_id=%s pdb_id=%s variant=%s ph=%s source=%s expected_hint=none",
                    self._run_id,
                    str(pdb_id),
                    str(variant or "*"),
                    str(ph or "*"),
                    str(consensus_source or ""),
                )
                return True
            backlog_hint = self._pending_count_locked() + 1
            scheduler = cfg.get("GLOBAL_DOCKING_SCHEDULER")
            free_now = (
                sem_free_slots(scheduler, self._total_cpu)
                if scheduler is not None
                else self._total_cpu
            )
            scheduler_queue_depth = _scheduler_queue_depth(scheduler)
            profile = _resolve_scorch_parallel_profile(cfg)
            self._refresh_worker_target_locked(cfg)
            cpu_budget_hint, free_cores_hint, active_width = (
                self._cpu_budget_hints_locked(cfg=cfg, free_now=int(free_now))
            )
            self._logger.info(
                "[scorch-rescore.queue.submit] run_id=%s pdb_id=%s variant=%s ph=%s active_width=%d cpu_budget=%d free_cores=%d total_cpu=%d free_now=%d",
                self._run_id,
                str(pdb_id),
                str(variant or "*"),
                str(ph or "*"),
                int(active_width),
                int(cpu_budget_hint),
                int(free_cores_hint),
                int(self._total_cpu),
                int(free_now),
            )
            future = self._pool.submit(
                _maybe_run_scorch_rescore_for_pdb,
                cfg,
                self._run_id,
                str(pdb_id),
                variant=variant,
                ph=ph,
                backlog_hint=int(backlog_hint),
                free_cores_hint=int(free_cores_hint),
                scheduler_queue_depth_hint=int(scheduler_queue_depth),
                allowed_count_hint=(
                    max(1, int(allowed_count_hint))
                    if allowed_count_hint is not None
                    else None
                ),
                cpu_budget_hint=int(cpu_budget_hint),
                profile_hint=profile,
                execution_mode=str(self._execution_mode),
                verbose=self._verbose,
            )
            self._inflight_keys.add(scope_key)
            self._futures[future] = {
                "pdb_id": str(pdb_id),
                "variant": str(variant or "*"),
                "ph": str(ph or "*"),
                "scope_key": scope_key,
                "submitted_at": time.time(),
            }
            self._submitted += 1
            return True

    def submit_provisional(
        self,
        *,
        cfg: Mapping[str, Any],
        pdb_id: str,
        variant: Optional[str] = None,
        ph: Optional[str] = None,
        phase_key: str = "reservoir",
        allowed_count_hint: Optional[int] = None,
        idle_worker_count: int = 0,
        active_worker_count: int = 0,
    ) -> bool:
        if not _scorch_cache_mod.provisional_run_enabled(cfg):
            return False
        reserved_workers = max(
            1,
            _as_int(cfg.get("SCORCH_PROVISIONAL_RESERVED_WORKERS"), 1),
        )
        if int(idle_worker_count) < int(reserved_workers):
            return False
        with self._lock:
            self._reap_done_locked()
            if self._closed:
                return False
            scope_key = (
                f"provisional:{str(phase_key or 'reservoir')}",
                str(pdb_id or "").strip().upper(),
                str(variant or "*").strip().upper() or "*",
                str(ph or "base").strip() or "base",
            )
            if scope_key in self._inflight_keys:
                return True
            if self._scope_active_locked(
                pdb_id=str(pdb_id),
                variant=variant,
                ph=ph,
                phase_prefix="final",
            ):
                return False
            if self._scope_active_locked(
                pdb_id=str(pdb_id),
                variant=variant,
                ph=ph,
                phase_prefix="provisional:",
            ):
                return False
            provisional_cap = max(
                1,
                _as_int(cfg.get("SCORCH_PROVISIONAL_MAX_INFLIGHT"), 1),
            )
            if self._provisional_pending_count_locked() >= int(provisional_cap):
                return False
            scheduler = cfg.get("GLOBAL_DOCKING_SCHEDULER")
            free_now = (
                sem_free_slots(scheduler, self._total_cpu)
                if scheduler is not None
                else int(idle_worker_count)
            )
            if int(free_now) < int(reserved_workers):
                return False
            provisional_cpu = max(
                1,
                _as_int(
                    cfg.get("SCORCH_PROVISIONAL_CPU_BUDGET"),
                    min(2, max(1, int(self._total_cpu))),
                ),
            )
            provisional_cpu = min(
                int(provisional_cpu),
                max(1, int(idle_worker_count)),
                max(1, int(free_now)),
            )
            future = self._pool.submit(
                _run_provisional_scorch_rescore_for_pdb,
                cfg,
                self._run_id,
                str(pdb_id),
                provisional_cpu=int(provisional_cpu),
                variant=variant,
                ph=ph,
                allowed_count_hint=(
                    max(1, int(allowed_count_hint))
                    if allowed_count_hint is not None
                    else None
                ),
                execution_mode=str(self._execution_mode),
                verbose=self._verbose,
            )
            self._inflight_keys.add(scope_key)
            self._futures[future] = {
                "pdb_id": str(pdb_id),
                "variant": str(variant or "*"),
                "ph": str(ph or "*"),
                "scope_key": scope_key,
                "submitted_at": time.time(),
                "phase": "provisional",
                "active_worker_count": int(active_worker_count),
                "idle_worker_count": int(idle_worker_count),
            }
            self._submitted += 1
            self._logger.info(
                "[scorch-rescore.queue] action=submit_provisional run_id=%s pdb_id=%s variant=%s ph=%s phase_key=%s cpu_budget=%d idle_workers=%d active_workers=%d",
                self._run_id,
                str(pdb_id),
                str(variant or "*"),
                str(ph or "*"),
                str(phase_key or "reservoir"),
                int(provisional_cpu),
                int(idle_worker_count),
                int(active_worker_count),
            )
            return True

    def close_and_wait(self) -> Dict[str, int]:
        with self._lock:
            self._closed = True
            futures = list(self._futures.keys())
        pending = set(futures)
        wait_started = time.time()
        wait_deadline = (
            wait_started + float(self._drain_timeout_sec)
            if self._drain_timeout_sec > 0.0
            else None
        )
        last_wait_log = wait_started
        while pending:
            now = time.time()
            if wait_deadline is not None and now >= wait_deadline:
                with self._lock:
                    still_pending = {
                        future: self._futures.pop(future, {})
                        for future in list(pending)
                        if future in self._futures and not future.done()
                    }
                    for meta in still_pending.values():
                        self._release_scope_key_locked(meta)
                    self._failed += len(still_pending)
                for meta in still_pending.values():
                    self._logger.warning(
                        "[scorch-rescore.queue.timeout] run_id=%s pdb_id=%s variant=%s ph=%s phase=%s elapsed_s=%.1f timeout_s=%.1f",
                        self._run_id,
                        str(meta.get("pdb_id", "")),
                        str(meta.get("variant", "*")),
                        str(meta.get("ph", "*")),
                        str(meta.get("phase", "final")),
                        max(0.0, now - wait_started),
                        float(self._drain_timeout_sec),
                    )
                try:
                    self._pool.shutdown(wait=False, cancel_futures=True)
                except TypeError:
                    self._pool.shutdown(wait=False)
                pending = set()
                break
            done, pending = wait(
                pending,
                timeout=15.0,
                return_when=FIRST_COMPLETED,
            )
            if not done:
                now = time.time()
                if now - last_wait_log >= 15.0:
                    self._logger.info(
                        "[scorch-rescore.queue.wait] run_id=%s pending=%d submitted=%d elapsed_s=%.1f",
                        self._run_id,
                        len(pending),
                        int(self._submitted),
                        max(0.0, now - wait_started),
                    )
                    last_wait_log = now
                continue
            with self._lock:
                for future in done:
                    meta = self._futures.pop(future, {})
                    self._record_done_future_locked(future, meta)
                pending = {
                    future
                    for future in pending
                    if future in self._futures and not future.done()
                }
            now = time.time()
            if now - last_wait_log >= 15.0:
                self._logger.info(
                    "[scorch-rescore.queue.wait] run_id=%s pending=%d submitted=%d completed=%d failed=%d elapsed_s=%.1f",
                    self._run_id,
                    len(pending),
                    int(self._submitted),
                    int(self._submitted - len(pending)),
                    int(self._failed),
                    max(0.0, now - wait_started),
                )
                last_wait_log = now
        self._pool.shutdown(wait=True)
        elapsed = max(0.001, time.time() - float(self._started_at))
        busy_ratio = max(
            0.0,
            min(
                1.0,
                float(self._busy_seconds)
                / max(0.001, float(elapsed) * max(1, int(self._peak_target_workers))),
            ),
        )
        self._logger.info(
            "[scorch.utilization] run_id=%s mode=%s worker_busy_ratio=%.3f elapsed_s=%.2f busy_s=%.2f workers_peak=%d",
            self._run_id,
            str(self._execution_mode),
            float(busy_ratio),
            float(elapsed),
            float(self._busy_seconds),
            int(self._peak_target_workers),
        )
        self._logger.info(
            "[scorch-rescore.queue.summary] run_id=%s submitted=%d failed=%d workers_current=%d workers_peak=%d cap=%d resized=%d mode=%s",
            self._run_id,
            self._submitted,
            self._failed,
            self._target_workers,
            self._peak_target_workers,
            self._worker_cap,
            self._resized,
            str(self._execution_mode),
        )
        return {
            "submitted": int(self._submitted),
            "failed": int(self._failed),
            "workers": int(self._target_workers),
            "workers_peak": int(self._peak_target_workers),
            "worker_cap": int(self._worker_cap),
            "resized": int(self._resized),
            "busy_ratio_milli": int(round(float(busy_ratio) * 1000.0)),
        }


def _start_scorch_mixed_queue_service(
    cfg: Mapping[str, Any], run_id: str, *, verbose: bool = False
) -> Optional[ScorchMixedQueueService]:
    logger = logging.getLogger("scorch-rescore-hook")
    if not run_id or not _is_scorch_enabled(cfg):
        return None
    if _is_no_library_docking(cfg):
        logger.info("[scorch-rescore.queue.skip] reason=no_library_docking run_id=%s", run_id)
        return None

    total_cpu = max(
        1,
        _as_int(cfg.get("CPU") or (os.cpu_count() or 1), os.cpu_count() or 1),
    )
    execution_mode = _resolve_scorch_execution_mode(cfg)
    device_plan = build_scorch_device_plan(cfg, logger)
    worker_cap = max(1, min(16, max(2, total_cpu // 8 + 1)))
    worker_cap = cap_scorch_queue_workers(worker_cap, device_plan)
    cap_raw = os.environ.get("ATLAS_SCORCH_QUEUE_MAX_WORKERS")
    if cap_raw is None:
        cap_raw = cfg.get("SCORCH_QUEUE_MAX_WORKERS")
    if (
        cap_raw is None
        and device_plan.effective != "gpu"
        and shard_mode_enabled(cfg)
    ):
        shard_outer_cap = max(2, min(8, int(total_cpu) // 14))
        cpu_shard_cap = max(2, min(int(worker_cap), int(shard_outer_cap)))
        worker_cap = max(1, int(cpu_shard_cap))
    if cap_raw is not None and str(cap_raw).strip():
        try:
            worker_cap = max(1, min(int(worker_cap), int(cap_raw)))
        except Exception:
            logger.warning(
                "[scorch-rescore.queue] action=cap_parse_failed raw=%s",
                cap_raw,
            )
    drain_timeout_raw = os.environ.get("ATLAS_SCORCH_QUEUE_DRAIN_TIMEOUT_SEC")
    if drain_timeout_raw is None:
        drain_timeout_raw = cfg.get("SCORCH_QUEUE_DRAIN_TIMEOUT_SEC")
    try:
        drain_timeout_sec = float(str(drain_timeout_raw).strip()) if drain_timeout_raw is not None else 0.0
    except Exception:
        drain_timeout_sec = 0.0
    service = ScorchMixedQueueService(
        run_id=run_id,
        worker_cap=worker_cap,
        total_cpu=total_cpu,
        execution_mode=execution_mode,
        logger=logger,
        verbose=verbose,
        drain_timeout_sec=drain_timeout_sec,
    )
    logger.info(
        "[scorch-rescore.queue.start] run_id=%s workers_initial=%d workers_cap=%d total_cpu=%d mode=%s device=%s backend=%s gpu_count=%d",
        run_id,
        1,
        worker_cap,
        total_cpu,
        execution_mode,
        device_plan.effective,
        device_plan.backend,
        device_plan.gpu_count,
    )
    return service
