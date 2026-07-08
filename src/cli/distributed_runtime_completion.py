from __future__ import annotations

import logging
import threading
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Optional

from cli.run_context import ConfigDict

ComboScope = tuple[str, str, str]


@dataclass
class DeferredScorchCompletionBookkeeper:
    cfg_v: ConfigDict
    run_id: str
    label: str
    dist_task_id: int
    combo_total_chunks: Mapping[ComboScope, int]
    combo_total_ligands: Mapping[ComboScope, int]
    combo_owner_task_id: Mapping[ComboScope, int]
    scorch_queue_service: Any
    variant_param_for_scorch: Optional[str]
    hybrid_start_ts: float
    acquire_global_slot: Callable[..., Any]
    maybe_run_scorch_rescore_for_pdb: Callable[..., Any]
    normalize_ph_tag_token: Callable[[Optional[str]], str]
    combo_progress_reader: Optional[
        Callable[[], Mapping[ComboScope, tuple[int, ...]]]
    ] = None
    before_queue_combo_scorch: Optional[Callable[[ComboScope], None]] = None
    global_reconcile_interval_sec: float = 3.0
    combo_terminal_chunks: dict[ComboScope, int] = field(default_factory=dict)
    combo_completed_chunks_live: dict[ComboScope, int] = field(default_factory=dict)
    combo_failed_chunks_live: dict[ComboScope, int] = field(default_factory=dict)
    combo_scorch_submitted: set[ComboScope] = field(default_factory=set)
    combo_scorch_submitting: set[ComboScope] = field(default_factory=set)
    combo_provisional_submitted: set[tuple[ComboScope, str]] = field(default_factory=set)
    combo_provisional_submitting: set[tuple[ComboScope, str]] = field(default_factory=set)
    hybrid_completed_scope_times: dict[ComboScope, float] = field(default_factory=dict)
    hybrid_first_complete_elapsed: Optional[float] = None
    _last_global_reconcile_ts: float = 0.0
    _global_reconcile_inflight: bool = False
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def snapshot_for_rebalance(
        self,
    ) -> tuple[dict[ComboScope, int], dict[ComboScope, int], Optional[float]]:
        with self._lock:
            return (
                dict(self.combo_terminal_chunks),
                dict(self.combo_completed_chunks_live),
                (
                    None
                    if self.hybrid_first_complete_elapsed is None
                    else float(self.hybrid_first_complete_elapsed)
                ),
            )

    def _record_combo_scope_completion_locked(self, combo_scope: ComboScope) -> None:
        if combo_scope in self.hybrid_completed_scope_times:
            return
        total_chunks = int(self.combo_total_chunks.get(combo_scope, 0))
        done_chunks = int(self.combo_terminal_chunks.get(combo_scope, 0))
        completed_chunks = int(self.combo_completed_chunks_live.get(combo_scope, 0))
        failed_chunks = int(self.combo_failed_chunks_live.get(combo_scope, 0))
        if (
            total_chunks <= 0
            or done_chunks < total_chunks
            or completed_chunks < total_chunks
            or failed_chunks > 0
        ):
            return
        elapsed = float(max(0.0, time.time() - self.hybrid_start_ts))
        self.hybrid_completed_scope_times[combo_scope] = float(elapsed)
        if (
            self.hybrid_first_complete_elapsed is None
            or elapsed < float(self.hybrid_first_complete_elapsed)
        ):
            self.hybrid_first_complete_elapsed = float(elapsed)
        logging.info(
            "[hybrid.progress] scope=%s ratio=1.000 eta_s=0.0",
            "|".join(str(x) for x in combo_scope),
        )

    def _reserve_combo_scorch_submission_locked(
        self, combo_scope: ComboScope
    ) -> bool:
        if combo_scope in self.combo_scorch_submitted:
            return False
        if combo_scope in self.combo_scorch_submitting:
            return False
        self.combo_scorch_submitting.add(combo_scope)
        return True

    def _finish_combo_scorch_submission(
        self, combo_scope: ComboScope, *, success: bool
    ) -> None:
        with self._lock:
            self.combo_scorch_submitting.discard(combo_scope)
            if success:
                self.combo_scorch_submitted.add(combo_scope)

    def _queue_deferred_combo_scorch(
        self, combo_scope: ComboScope, *, reason: str
    ) -> bool:
        if self.scorch_queue_service is None:
            return False
        pdb_scope, _variant_scope, ph_scope = combo_scope
        ph_value = (
            None if str(ph_scope).strip().lower() == "base" else str(ph_scope)
        )
        allowed_hint = int(max(1, self.combo_total_ligands.get(combo_scope, 1)))
        if self.before_queue_combo_scorch is not None:
            try:
                self.before_queue_combo_scorch(combo_scope)
            except Exception:
                logging.warning(
                    "[scorch-rescore.queue] action=pre_submit_aggregate_failed pdb_id=%s variant=%s ph=%s reason=%s",
                    str(pdb_scope),
                    self.label,
                    str(ph_scope),
                    str(reason),
                    exc_info=True,
                )
        try:
            queued = bool(
                self.scorch_queue_service.submit(
                    cfg=self.cfg_v,
                    pdb_id=str(pdb_scope),
                    variant=self.variant_param_for_scorch,
                    ph=ph_value,
                    allowed_count_hint=allowed_hint,
                )
            )
        except Exception:
            queued = False
            logging.warning(
                "[scorch-rescore.queue] action=submit_failed_deferred pdb_id=%s variant=%s ph=%s reason=%s",
                str(pdb_scope),
                self.label,
                str(ph_scope),
                str(reason),
                exc_info=True,
            )
        return bool(queued)

    def _provisional_enabled(self) -> bool:
        value = self.cfg_v.get("SCORCH_PROVISIONAL_ENABLE")
        if isinstance(value, bool):
            return bool(value)
        return str(value or "").strip().lower() in {"1", "true", "yes", "on"}

    def _scorch_shard_mode_enabled(self) -> bool:
        try:
            from post_docking.rescoring.scorch_shards import shard_mode_enabled

            return bool(shard_mode_enabled(self.cfg_v))
        except Exception:
            return False

    def _provisional_min_progress(self) -> float:
        try:
            return max(
                0.0,
                min(1.0, float(self.cfg_v.get("SCORCH_PROVISIONAL_MIN_PROGRESS", 0.10))),
            )
        except Exception:
            return 0.10

    def maybe_queue_provisional_combo_scorch(
        self,
        combo_scope: ComboScope,
        *,
        reason: str,
        active_workers: int,
        idle_workers: int,
        unresolved_chunks: int,
    ) -> bool:
        if not self._provisional_enabled() or self.scorch_queue_service is None:
            return False
        if not hasattr(self.scorch_queue_service, "submit_provisional"):
            return False
        if int(idle_workers) <= 0:
            return False
        with self._lock:
            owner_task = int(self.combo_owner_task_id.get(combo_scope, -1))
            if owner_task != int(self.dist_task_id):
                return False
            total_chunks = int(self.combo_total_chunks.get(combo_scope, 0))
            completed_chunks = int(
                self.combo_completed_chunks_live.get(combo_scope, 0)
            )
            terminal_chunks = int(self.combo_terminal_chunks.get(combo_scope, 0))
            if total_chunks <= 0 or completed_chunks <= 0:
                return False
            if terminal_chunks >= total_chunks:
                return False
            progress = float(completed_chunks) / float(max(1, total_chunks))
            if progress < self._provisional_min_progress():
                return False
            try:
                max_per_combo = max(
                    1,
                    int(self.cfg_v.get("SCORCH_PROVISIONAL_MAX_PER_COMBO", 1)),
                )
            except Exception:
                max_per_combo = 1
            submitted_for_combo = sum(
                1
                for submitted_scope, _phase in self.combo_provisional_submitted
                if submitted_scope == combo_scope
            )
            inflight_for_combo = sum(
                1
                for submitted_scope, _phase in self.combo_provisional_submitting
                if submitted_scope == combo_scope
            )
            if submitted_for_combo + inflight_for_combo >= int(max_per_combo):
                return False
            progress_bucket = max(1, min(9, int(progress * 10.0)))
            phase_key = f"p{progress_bucket * 10:03d}"
            reserve_key = (combo_scope, phase_key)
            if (
                reserve_key in self.combo_provisional_submitted
                or reserve_key in self.combo_provisional_submitting
            ):
                return False
            self.combo_provisional_submitting.add(reserve_key)
            allowed_hint = int(
                max(1, self.combo_total_ligands.get(combo_scope, completed_chunks))
            )
        pdb_scope, _variant_scope, ph_scope = combo_scope
        ph_value = None if str(ph_scope).strip().lower() == "base" else str(ph_scope)
        queued = False
        try:
            queued = bool(
                self.scorch_queue_service.submit_provisional(
                    cfg=self.cfg_v,
                    pdb_id=str(pdb_scope),
                    variant=self.variant_param_for_scorch,
                    ph=ph_value,
                    phase_key=phase_key,
                    allowed_count_hint=allowed_hint,
                    idle_worker_count=int(idle_workers),
                    active_worker_count=int(active_workers),
                )
            )
        except Exception:
            queued = False
            logging.warning(
                "[scorch-rescore.provisional] action=submit_failed pdb_id=%s variant=%s ph=%s reason=%s",
                str(pdb_scope),
                self.label,
                str(ph_scope),
                str(reason),
                exc_info=True,
            )
        with self._lock:
            reserve_key = (combo_scope, phase_key)
            self.combo_provisional_submitting.discard(reserve_key)
            if queued:
                self.combo_provisional_submitted.add(reserve_key)
        if queued:
            logging.info(
                "[scorch-rescore.provisional] action=queued scope=%s progress=%s active_workers=%d idle_workers=%d unresolved_chunks=%d reason=%s",
                "|".join(str(x) for x in combo_scope),
                phase_key,
                int(active_workers),
                int(idle_workers),
                int(unresolved_chunks),
                str(reason),
            )
        return bool(queued)

    def reconcile_provisional_owned_submissions(
        self,
        *,
        reason: str,
        active_workers: int,
        idle_workers: int,
        unresolved_chunks: int,
    ) -> int:
        if not self._provisional_enabled() or int(idle_workers) <= 0:
            return 0
        queued = 0
        with self._lock:
            scopes = [
                scope
                for scope, owner_task in sorted(self.combo_owner_task_id.items())
                if int(owner_task) == int(self.dist_task_id)
                and int(self.combo_completed_chunks_live.get(scope, 0)) > 0
                and int(self.combo_failed_chunks_live.get(scope, 0)) <= 0
                and int(self.combo_terminal_chunks.get(scope, 0))
                < int(self.combo_total_chunks.get(scope, 0))
            ]
        for scope in scopes:
            if self.maybe_queue_provisional_combo_scorch(
                scope,
                reason=reason,
                active_workers=active_workers,
                idle_workers=idle_workers,
                unresolved_chunks=unresolved_chunks,
            ):
                queued += 1
                break
        return int(queued)

    def _read_global_combo_progress(
        self,
    ) -> Mapping[ComboScope, tuple[int, ...]]:
        if self.combo_progress_reader is None:
            return {}
        try:
            progress = self.combo_progress_reader()
        except Exception:
            logging.warning(
                "[scorch-rescore.defer] action=global_progress_read_failed variant=%s task_id=%d",
                self.label.upper(),
                int(self.dist_task_id),
                exc_info=True,
            )
            return {}
        return progress if hasattr(progress, "items") else {}

    @staticmethod
    def _coerce_progress_counts(raw_counts: object) -> tuple[int, int, int]:
        try:
            done_chunks = int(raw_counts[0])  # type: ignore[index]
            completed_chunks = int(raw_counts[1])  # type: ignore[index]
        except Exception:
            raise ValueError("invalid combo progress counts") from None
        try:
            failed_chunks = int(raw_counts[2])  # type: ignore[index]
        except Exception:
            failed_chunks = max(0, int(done_chunks) - int(completed_chunks))
        return (
            max(0, int(done_chunks)),
            max(0, int(completed_chunks)),
            max(0, int(failed_chunks)),
        )

    def _reserve_ready_owned_scopes_from_progress(
        self,
        progress: Mapping[ComboScope, tuple[int, ...]],
        *,
        reason: str,
    ) -> list[ComboScope]:
        ready_scopes: list[ComboScope] = []
        with self._lock:
            for combo_scope, raw_counts in progress.items():
                try:
                    done_chunks, completed_chunks, failed_chunks = (
                        self._coerce_progress_counts(raw_counts)
                    )
                except Exception:
                    continue
                total_chunks = int(self.combo_total_chunks.get(combo_scope, 0))
                owner_task = int(self.combo_owner_task_id.get(combo_scope, -1))
                if owner_task != int(self.dist_task_id):
                    continue
                if total_chunks <= 0 or done_chunks < total_chunks:
                    continue
                if completed_chunks <= 0:
                    continue
                if (
                    combo_scope in self.combo_scorch_submitted
                    or combo_scope in self.combo_scorch_submitting
                ):
                    continue
                self.combo_terminal_chunks[combo_scope] = max(
                    int(self.combo_terminal_chunks.get(combo_scope, 0)),
                    int(done_chunks),
                )
                self.combo_completed_chunks_live[combo_scope] = max(
                    int(self.combo_completed_chunks_live.get(combo_scope, 0)),
                    int(completed_chunks),
                )
                self.combo_failed_chunks_live[combo_scope] = max(
                    int(self.combo_failed_chunks_live.get(combo_scope, 0)),
                    int(failed_chunks),
                )
                if failed_chunks > 0:
                    logging.warning(
                        "[scorch-rescore.defer] action=skip_final_failed_chunks_present variant=%s task_id=%d scope=%s terminal_chunks=%d completed_chunks=%d failed_chunks=%d total_chunks=%d reason=%s",
                        self.label.upper(),
                        int(self.dist_task_id),
                        "|".join(str(x) for x in combo_scope),
                        int(done_chunks),
                        int(completed_chunks),
                        int(failed_chunks),
                        int(total_chunks),
                        str(reason),
                    )
                    continue
                if completed_chunks < total_chunks:
                    continue
                self._record_combo_scope_completion_locked(combo_scope)
                if not self._reserve_combo_scorch_submission_locked(combo_scope):
                    continue
                logging.info(
                    "[scorch-rescore.defer] action=scope_ready variant=%s task_id=%d scope=%s terminal_chunks=%d completed_chunks=%d total_chunks=%d reason=%s",
                    self.label.upper(),
                    int(self.dist_task_id),
                    "|".join(str(x) for x in combo_scope),
                    int(done_chunks),
                    int(completed_chunks),
                    int(total_chunks),
                    str(reason),
                )
                ready_scopes.append(combo_scope)
        return ready_scopes

    def reconcile_ready_owned_deferred_submissions(
        self, *, reason: str, force: bool = False
    ) -> int:
        if self.scorch_queue_service is None or self.combo_progress_reader is None:
            return 0
        now = time.time()
        with self._lock:
            if self._global_reconcile_inflight:
                return 0
            if (
                not force
                and self._last_global_reconcile_ts > 0.0
                and now - float(self._last_global_reconcile_ts)
                < float(self.global_reconcile_interval_sec)
            ):
                return 0
            self._global_reconcile_inflight = True
            self._last_global_reconcile_ts = float(now)
        queued_count = 0
        try:
            progress = self._read_global_combo_progress()
            ready_scopes = self._reserve_ready_owned_scopes_from_progress(
                progress,
                reason=reason,
            )
            for combo_scope in ready_scopes:
                queued = False
                try:
                    queued = self._queue_deferred_combo_scorch(
                        combo_scope,
                        reason=reason,
                    )
                finally:
                    self._finish_combo_scorch_submission(
                        combo_scope,
                        success=bool(queued),
                    )
                if queued:
                    queued_count += 1
        finally:
            with self._lock:
                self._global_reconcile_inflight = False
        return int(queued_count)

    def maybe_queue_deferred_combo_scorch(
        self, combo_scope: ComboScope, *, reason: str
    ) -> bool:
        with self._lock:
            total_chunks = int(self.combo_total_chunks.get(combo_scope, 0))
            done_chunks = int(self.combo_terminal_chunks.get(combo_scope, 0))
            completed_chunks = int(self.combo_completed_chunks_live.get(combo_scope, 0))
            failed_chunks = int(self.combo_failed_chunks_live.get(combo_scope, 0))
            owner_task = int(self.combo_owner_task_id.get(combo_scope, -1))
            already_submitted = combo_scope in self.combo_scorch_submitted
            submission_inflight = combo_scope in self.combo_scorch_submitting
            if already_submitted or submission_inflight:
                return True
            if owner_task != int(self.dist_task_id):
                return False
            if total_chunks <= 0 or done_chunks < total_chunks:
                return False
            if failed_chunks > 0:
                logging.warning(
                    "[scorch-rescore.defer] action=skip_final_failed_chunks_present variant=%s task_id=%d scope=%s terminal_chunks=%d completed_chunks=%d failed_chunks=%d total_chunks=%d reason=%s",
                    self.label.upper(),
                    int(self.dist_task_id),
                    "|".join(str(x) for x in combo_scope),
                    int(done_chunks),
                    int(completed_chunks),
                    int(failed_chunks),
                    int(total_chunks),
                    str(reason),
                )
                return False
            if completed_chunks < total_chunks:
                return False
            if not self._reserve_combo_scorch_submission_locked(combo_scope):
                return True
        queued = False
        try:
            queued = self._queue_deferred_combo_scorch(combo_scope, reason=reason)
        finally:
            self._finish_combo_scorch_submission(combo_scope, success=bool(queued))
        return bool(queued)

    def mark_chunk_completed(self, combo_scope: ComboScope, *, reason: str) -> None:
        with self._lock:
            self.combo_terminal_chunks[combo_scope] = (
                int(self.combo_terminal_chunks.get(combo_scope, 0)) + 1
            )
            self.combo_completed_chunks_live[combo_scope] = (
                int(self.combo_completed_chunks_live.get(combo_scope, 0)) + 1
            )
            self._record_combo_scope_completion_locked(combo_scope)
        self.maybe_queue_deferred_combo_scorch(combo_scope, reason=reason)

    def mark_chunk_terminal_failure(
        self, combo_scope: ComboScope, *, reason: str
    ) -> None:
        with self._lock:
            self.combo_terminal_chunks[combo_scope] = (
                int(self.combo_terminal_chunks.get(combo_scope, 0)) + 1
            )
            self.combo_failed_chunks_live[combo_scope] = (
                int(self.combo_failed_chunks_live.get(combo_scope, 0)) + 1
            )
            self._record_combo_scope_completion_locked(combo_scope)
        logging.warning(
            "[scorch-rescore.defer] action=record_terminal_failure_skip_final scope=%s reason=%s",
            "|".join(str(x) for x in combo_scope),
            str(reason),
        )

    def flush_owned_deferred_submissions(
        self,
        *,
        combo_completed_chunks: Mapping[ComboScope, int],
        combo_failed_chunks: Optional[Mapping[ComboScope, int]] = None,
    ) -> None:
        if self.scorch_queue_service is None or not self.combo_owner_task_id:
            return
        owned_combo_scopes = [
            scope
            for scope, owner_task in sorted(self.combo_owner_task_id.items())
            if int(owner_task) == int(self.dist_task_id)
            and int(combo_completed_chunks.get(scope, 0)) > 0
        ]
        with self._lock:
            already_submitted = {
                scope for scope in owned_combo_scopes if scope in self.combo_scorch_submitted
            }
            already_submitting = {
                scope for scope in owned_combo_scopes if scope in self.combo_scorch_submitting
            }
        logging.info(
            "[scorch-rescore.defer] variant=%s task_id=%d combos_owned=%d combos_with_completed_chunks=%d combos_already_submitted=%d combos_submission_inflight=%d",
            self.label.upper(),
            int(self.dist_task_id),
            len(owned_combo_scopes),
            len(combo_completed_chunks),
            len(already_submitted),
            len(already_submitting),
        )
        for combo_scope in owned_combo_scopes:
            if combo_scope in already_submitted or combo_scope in already_submitting:
                continue
            total_chunks = int(self.combo_total_chunks.get(combo_scope, 0))
            completed_chunks = int(combo_completed_chunks.get(combo_scope, 0))
            failed_chunks = int(
                (combo_failed_chunks or {}).get(
                    combo_scope, self.combo_failed_chunks_live.get(combo_scope, 0)
                )
            )
            with self._lock:
                self.combo_completed_chunks_live[combo_scope] = max(
                    int(self.combo_completed_chunks_live.get(combo_scope, 0)),
                    int(completed_chunks),
                )
                self.combo_failed_chunks_live[combo_scope] = max(
                    int(self.combo_failed_chunks_live.get(combo_scope, 0)),
                    int(failed_chunks),
                )
                self.combo_terminal_chunks[combo_scope] = max(
                    int(self.combo_terminal_chunks.get(combo_scope, 0)),
                    int(completed_chunks) + int(failed_chunks),
                )
            if total_chunks <= 0 or completed_chunks < total_chunks:
                logging.info(
                    "[scorch-rescore.defer] action=skip_final_incomplete_scope variant=%s task_id=%d scope=%s completed_chunks=%d failed_chunks=%d total_chunks=%d",
                    self.label.upper(),
                    int(self.dist_task_id),
                    "|".join(str(x) for x in combo_scope),
                    int(completed_chunks),
                    int(failed_chunks),
                    int(total_chunks),
                )
                continue
            if failed_chunks > 0:
                logging.warning(
                    "[scorch-rescore.defer] action=skip_final_failed_chunks_present variant=%s task_id=%d scope=%s completed_chunks=%d failed_chunks=%d total_chunks=%d reason=final_flush",
                    self.label.upper(),
                    int(self.dist_task_id),
                    "|".join(str(x) for x in combo_scope),
                    int(completed_chunks),
                    int(failed_chunks),
                    int(total_chunks),
                )
                continue
            with self._lock:
                if not self._reserve_combo_scorch_submission_locked(combo_scope):
                    continue
            pdb_scope, _variant_scope, ph_scope = combo_scope
            ph_value = (
                None if str(ph_scope).strip().lower() == "base" else str(ph_scope)
            )
            allowed_hint = int(max(1, self.combo_total_ligands.get(combo_scope, 1)))
            queued_scorch = False
            direct_success = False
            try:
                queued_scorch = bool(
                    self.scorch_queue_service.submit(
                        cfg=self.cfg_v,
                        pdb_id=str(pdb_scope),
                        variant=self.variant_param_for_scorch,
                        ph=ph_value,
                        allowed_count_hint=allowed_hint,
                    )
                )
            except Exception:
                queued_scorch = False
                logging.warning(
                    "[scorch-rescore.queue] action=submit_failed_deferred pdb_id=%s variant=%s ph=%s",
                    str(pdb_scope),
                    self.label,
                    str(ph_scope),
                    exc_info=True,
                )
            if queued_scorch:
                self._finish_combo_scorch_submission(combo_scope, success=True)
                continue
            try:
                with self.acquire_global_slot(
                    self.cfg_v,
                    cores=1,
                    min_cores=1,
                    priority=2,
                    task_id=f"main:scorch:deferred:{pdb_scope}:{self.label}:{ph_scope}",
                    task_type="scorch_score_chunk",
                ):
                    rc = self.maybe_run_scorch_rescore_for_pdb(
                        self.cfg_v,
                        self.run_id,
                        str(pdb_scope),
                        variant=self.variant_param_for_scorch,
                        ph=ph_value,
                        allowed_count_hint=allowed_hint,
                        verbose=False,
                    )
                    direct_success = bool(rc == 0)
            except Exception:
                logging.warning(
                    "[scorch-rescore.hook] action=skip reason=unexpected_exception_deferred pdb_id=%s variant=%s ph=%s",
                    str(pdb_scope),
                    self.label,
                    str(ph_scope),
                    exc_info=True,
                )
            finally:
                self._finish_combo_scorch_submission(
                    combo_scope, success=bool(direct_success)
                )
