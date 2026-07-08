from __future__ import annotations

import logging
import os
import threading
import time
from typing import Optional

from cli.distributed_chunk_planner import (
    _combo_pdb_id_from_file,
    _combo_prep_key,
    _combo_receptor_ready,
    _normalize_ph_tag_token,
    _scope_prep_key,
)
from cli.distributed_chunk_runtime_state import (
    add_assigned_pdb_id,
    apply_combo_prep_reuse_payload,
    build_combo_prep_ready_payload,
    combo_prep_has_reuse_data,
    combo_prep_status,
)
from cli.distributed_context import (
    combo_prep_failure_retryable,
    combo_prep_claim_path,
    combo_prep_state_path,
    read_combo_prep_state,
    renew_chunk_claim,
    try_begin_combo_prep,
    write_combo_prep_state,
)
from cli.distributed_runtime_claim_chunk_result import write_chunk_result_and_cache
from cli.distributed_runtime_claim_prep import (
    acquire_scope_prep_lock,
    combo_receptor_ready_after_settle,
    release_scope_prep_lock,
)
from cli.distributed_runtime_claim_types import ClaimLoopSession
from cli.distributed_runtime_completion import ComboScope
from cli.run_process_one import (
    _CHUNK_FAILURE_REASON_KEY,
    _CHUNK_FAILURE_RETRYABLE_KEY,
)


def process_claimed_chunk(
    session: ClaimLoopSession,
    claimed_id: str,
    claimed_prev_attempt: int,
) -> None:
    chunk = session.context.by_chunk_id.get(claimed_id)
    attempt_no = int(claimed_prev_attempt) + 1
    if not isinstance(chunk, dict):
        fail_status = (
            "terminal_failed"
            if int(attempt_no) >= int(session.cfg.chunk_max_attempts)
            else "failed"
        )
        write_chunk_result_and_cache(session, 
            chunk_id=claimed_id,
            status=fail_status,
            payload={
                "reason": "missing_chunk_payload",
                "attempt": int(attempt_no),
                "max_attempts": int(session.cfg.chunk_max_attempts),
            },
        )
        with session.progress_lock:
            session.context.bar.update(1)
        return

    pdb_file = str(chunk.get("pdb_file") or "")
    ph_tag_raw = chunk.get("ph_tag")
    ph_tag = None if ph_tag_raw is None else str(ph_tag_raw)
    chunk_ligands_raw = chunk.get("ligand_bases") or []
    chunk_ligands = [str(x) for x in chunk_ligands_raw if str(x).strip()]
    chunk_run_mode = str(chunk.get("run_mode") or "").strip()
    chunk_library_name = str(chunk.get("library_name") or "").strip()
    chunk_library_root = str(chunk.get("library_root") or "").strip()
    combo_scope: ComboScope = (
        _combo_pdb_id_from_file(pdb_file),
        str(session.context.label or "").strip().upper() or "BASE",
        _normalize_ph_tag_token(ph_tag),
    )
    heartbeat_stop = threading.Event()
    heartbeat_thread: Optional[threading.Thread] = None

    def _start_chunk_claim_heartbeat() -> None:
        nonlocal heartbeat_thread
        interval = max(
            5.0,
            min(
                float(session.cfg.chunk_claim_heartbeat_sec),
                float(session.cfg.chunk_claim_lease_sec) / 3.0,
            ),
        )

        def _heartbeat_loop() -> None:
            while not heartbeat_stop.wait(interval):
                if not renew_chunk_claim(
                    session.context.dist_ctx,
                    session.context.cfg_v,
                    chunk_id=claimed_id,
                    lease_sec=session.cfg.chunk_claim_lease_sec,
                ):
                    return

        heartbeat_thread = threading.Thread(
            target=_heartbeat_loop,
            name=f"chunk-claim-hb-{claimed_id[:8]}",
            daemon=True,
        )
        heartbeat_thread.start()

    _start_chunk_claim_heartbeat()
    combo_key = _combo_prep_key(
        run_id=session.context.run_id,
        variant_label=session.context.label,
        pdb_file=pdb_file,
        ph_tag=ph_tag,
    )
    scope_prep_key = _scope_prep_key(
        run_id=session.context.run_id,
        variant_label=session.context.label,
        pdb_file=pdb_file,
        ph_tag=ph_tag,
    )
    combo_cfg = session.context.cfg_v
    prep_owner = False
    prep_attempt_no = 1
    scope_owner_task = int(session.context.combo_owner_task_id.get(combo_scope, -1))
    is_designated_prep_owner = scope_owner_task < 0 or int(scope_owner_task) == int(
        session.context.dist_ctx.task_id
    )
    prep_status, prep_payload = combo_prep_status(
        dist_ctx=session.context.dist_ctx,
        cfg_v=session.context.cfg_v,
        combo_key=scope_prep_key,
        read_combo_prep_state=read_combo_prep_state,
    )
    try:
        if prep_status == "failed":
            prep_reason = (
                str((prep_payload or {}).get("error") or "").strip()
                or "owner_marked_failed"
            )
            if combo_prep_failure_retryable(session.context.cfg_v, prep_payload):
                logging.warning(
                    "[prep.retry] scope=%s combo=%s pdb=%s reason=%s",
                    scope_prep_key,
                    combo_key,
                    _combo_pdb_id_from_file(pdb_file),
                    prep_reason,
                )
                prep_status = ""
                prep_payload = None
            else:
                raise RuntimeError(
                    f"combo_prep_failed:scope={scope_prep_key}:reason={prep_reason}"
                )

        preparing_since = None
        try:
            preparing_since = float((prep_payload or {}).get("updated_at") or 0.0)
        except Exception:
            preparing_since = None

        if prep_status == "ready":
            combo_cfg = session.context.cfg_v.copy()
            combo_cfg["FORCE_REPROCESS"] = False
            reused = apply_combo_prep_reuse_payload(combo_cfg, prep_payload)
            if reused:
                combo_cfg["_DISTRIBUTED_PREP_REUSE_ONLY"] = True
                combo_cfg["_SKIP_CONTROL_REDOCK"] = True
                combo_cfg["_SKIP_CONTROL_DOCKING"] = True
            logging.info(
                "[prep.ready] scope=%s combo=%s pdb=%s variant=%s ph=%s source=state_ready reused=%s",
                scope_prep_key,
                combo_key,
                _combo_pdb_id_from_file(pdb_file),
                session.context.label.upper(),
                _normalize_ph_tag_token(ph_tag),
                str(reused).lower(),
            )
        elif (
            prep_status == "preparing"
            and combo_prep_has_reuse_data(prep_payload)
            and _combo_receptor_ready(
                session.context.cfg_v,
                pdb_file=pdb_file,
                variant_token=session.context.variant,
                ph_tag=ph_tag,
                min_mtime=preparing_since,
            )
        ):
            combo_cfg = session.context.cfg_v.copy()
            combo_cfg["FORCE_REPROCESS"] = False
            reused = apply_combo_prep_reuse_payload(combo_cfg, prep_payload)
            if reused:
                combo_cfg["_DISTRIBUTED_PREP_REUSE_ONLY"] = True
                combo_cfg["_SKIP_CONTROL_REDOCK"] = True
                combo_cfg["_SKIP_CONTROL_DOCKING"] = True
            logging.info(
                "[prep.ready] scope=%s combo=%s pdb=%s variant=%s ph=%s source=owner_artifact reused=%s",
                scope_prep_key,
                combo_key,
                _combo_pdb_id_from_file(pdb_file),
                session.context.label.upper(),
                _normalize_ph_tag_token(ph_tag),
                str(reused).lower(),
            )
        elif is_designated_prep_owner and try_begin_combo_prep(
            session.context.dist_ctx,
            session.context.cfg_v,
            combo_key=scope_prep_key,
            lease_sec=session.cfg.combo_prep_lease_sec,
        ):
            prep_owner = True
            _, owner_payload = combo_prep_status(
                dist_ctx=session.context.dist_ctx,
                cfg_v=session.context.cfg_v,
                combo_key=scope_prep_key,
                read_combo_prep_state=read_combo_prep_state,
            )
            try:
                prep_attempt_no = max(
                    1, int((owner_payload or {}).get("prep_attempt") or 1)
                )
            except Exception:
                prep_attempt_no = 1
            logging.info(
                "[prep.owner] scope=%s combo=%s owner_task=%d state=%s",
                scope_prep_key,
                combo_key,
                session.context.dist_ctx.task_id,
                combo_prep_state_path(
                    session.context.dist_ctx,
                    session.context.cfg_v,
                    combo_key=scope_prep_key,
                ),
            )
        else:
            wait_deadline = time.time() + session.cfg.combo_prep_wait_timeout_sec
            waiting_logged = False
            while time.time() < wait_deadline:
                status_now, payload_now = combo_prep_status(
                    dist_ctx=session.context.dist_ctx,
                    cfg_v=session.context.cfg_v,
                    combo_key=scope_prep_key,
                    read_combo_prep_state=read_combo_prep_state,
                )
                if status_now == "failed":
                    prep_reason = (
                        str((payload_now or {}).get("error") or "").strip()
                        or "owner_marked_failed"
                    )
                    if combo_prep_failure_retryable(session.context.cfg_v, payload_now):
                        if not waiting_logged:
                            waiting_logged = True
                            logging.warning(
                                "[prep.wait] scope=%s combo=%s waiter_task=%d state=failed_retryable reason=%s",
                                scope_prep_key,
                                combo_key,
                                session.context.dist_ctx.task_id,
                                prep_reason,
                            )
                        time.sleep(session.cfg.combo_prep_wait_poll_sec)
                        continue
                    raise RuntimeError(
                        f"combo_prep_failed:scope={scope_prep_key}:reason={prep_reason}"
                    )
                payload_mtime = None
                try:
                    payload_mtime = float((payload_now or {}).get("updated_at") or 0.0)
                except Exception:
                    payload_mtime = None
                if status_now == "ready":
                    combo_cfg = session.context.cfg_v.copy()
                    combo_cfg["FORCE_REPROCESS"] = False
                    reused = apply_combo_prep_reuse_payload(combo_cfg, payload_now)
                    if reused:
                        combo_cfg["_DISTRIBUTED_PREP_REUSE_ONLY"] = True
                        combo_cfg["_SKIP_CONTROL_REDOCK"] = True
                        combo_cfg["_SKIP_CONTROL_DOCKING"] = True
                    logging.info(
                        "[prep.ready] scope=%s combo=%s pdb=%s variant=%s ph=%s source=wait_state_ready reused=%s",
                        scope_prep_key,
                        combo_key,
                        _combo_pdb_id_from_file(pdb_file),
                        session.context.label.upper(),
                        _normalize_ph_tag_token(ph_tag),
                        str(reused).lower(),
                    )
                    break
                if (
                    status_now == "preparing"
                    and combo_prep_has_reuse_data(payload_now)
                    and _combo_receptor_ready(
                        session.context.cfg_v,
                        pdb_file=pdb_file,
                        variant_token=session.context.variant,
                        ph_tag=ph_tag,
                        min_mtime=payload_mtime,
                    )
                ):
                    combo_cfg = session.context.cfg_v.copy()
                    combo_cfg["FORCE_REPROCESS"] = False
                    reused = apply_combo_prep_reuse_payload(combo_cfg, payload_now)
                    if reused:
                        combo_cfg["_DISTRIBUTED_PREP_REUSE_ONLY"] = True
                        combo_cfg["_SKIP_CONTROL_REDOCK"] = True
                        combo_cfg["_SKIP_CONTROL_DOCKING"] = True
                    logging.info(
                        "[prep.ready] scope=%s combo=%s pdb=%s variant=%s ph=%s source=wait_owner_artifact reused=%s",
                        scope_prep_key,
                        combo_key,
                        _combo_pdb_id_from_file(pdb_file),
                        session.context.label.upper(),
                        _normalize_ph_tag_token(ph_tag),
                        str(reused).lower(),
                    )
                    break
                if not waiting_logged:
                    waiting_logged = True
                    logging.info(
                        "[prep.wait] scope=%s combo=%s waiter_task=%d state=%s",
                        scope_prep_key,
                        combo_key,
                        session.context.dist_ctx.task_id,
                        status_now or "none",
                    )
                time.sleep(session.cfg.combo_prep_wait_poll_sec)
            else:
                if try_begin_combo_prep(
                    session.context.dist_ctx,
                    session.context.cfg_v,
                    combo_key=scope_prep_key,
                    lease_sec=session.cfg.combo_prep_lease_sec,
                ):
                    prep_owner = True
                    _, owner_payload = combo_prep_status(
                        dist_ctx=session.context.dist_ctx,
                        cfg_v=session.context.cfg_v,
                        combo_key=scope_prep_key,
                        read_combo_prep_state=read_combo_prep_state,
                    )
                    try:
                        prep_attempt_no = max(
                            1, int((owner_payload or {}).get("prep_attempt") or 1)
                        )
                    except Exception:
                        prep_attempt_no = 1
                    logging.warning(
                        "[prep.owner.promote] scope=%s combo=%s owner_task=%d reason=wait_timeout",
                        scope_prep_key,
                        combo_key,
                        session.context.dist_ctx.task_id,
                    )
                else:
                    raise RuntimeError(
                        f"combo_prep_wait_timeout:scope={scope_prep_key}"
                    )

        combo_cfg = combo_cfg.copy()
        combo_cfg["CPU"] = int(session.context.cpu_per_chunk_worker)
        combo_cfg["THREADS_PER_VINA"] = 1
        if prep_owner:
            combo_cfg["_DISTRIBUTED_PREP_READY_SIGNAL"] = {
                "state_path": str(
                    combo_prep_state_path(
                        session.context.dist_ctx,
                        session.context.cfg_v,
                        combo_key=scope_prep_key,
                    )
                ),
                "claim_path": str(
                    combo_prep_claim_path(
                        session.context.dist_ctx,
                        session.context.cfg_v,
                        combo_key=scope_prep_key,
                    )
                ),
                "run_id": str(session.context.run_id),
                "combo_key": str(scope_prep_key),
                "task_id": int(session.context.dist_ctx.task_id),
                "task_count": int(session.context.dist_ctx.task_count),
                "owner_task_id": int(session.context.dist_ctx.task_id),
                "prep_attempt": int(prep_attempt_no),
                "pdb_id": _combo_pdb_id_from_file(pdb_file),
                "pdb_file": str(pdb_file),
                "variant_label": str(session.context.label).upper(),
                "ph_tag": _normalize_ph_tag_token(ph_tag),
            }
        else:
            combo_cfg.pop("_DISTRIBUTED_PREP_READY_SIGNAL", None)
        scope_prep_lock = None
        try:
            if prep_owner:
                scope_prep_lock = acquire_scope_prep_lock(
                    context=session.context,
                    cfg=session.cfg,
                    pdb_file=pdb_file,
                )
            ok = bool(
                session.context.process_one(
                    pdb_file,
                    combo_cfg,
                    ph_tag,
                    chunk_ligand_bases=chunk_ligands,
                    chunk_id=claimed_id,
                    chunk_run_mode=chunk_run_mode or None,
                    chunk_library_name=chunk_library_name or None,
                    chunk_library_root=chunk_library_root or None,
                )
            )
        finally:
            release_scope_prep_lock(scope_prep_lock)
        if prep_owner:
            if combo_receptor_ready_after_settle(
                combo_cfg,
                pdb_file=pdb_file,
                variant_token=session.context.variant,
                ph_tag=ph_tag,
                poll_sec=session.cfg.combo_prep_wait_poll_sec,
            ):
                write_combo_prep_state(
                    session.context.dist_ctx,
                    session.context.cfg_v,
                    combo_key=scope_prep_key,
                    status="ready",
                    payload=build_combo_prep_ready_payload(
                        cfg_source=combo_cfg,
                        owner_task_id=int(session.context.dist_ctx.task_id),
                        pdb_file=pdb_file,
                        variant_label=session.context.label,
                        ph_tag=ph_tag,
                    ),
                )
            else:
                write_combo_prep_state(
                    session.context.dist_ctx,
                    session.context.cfg_v,
                    combo_key=scope_prep_key,
                    status="failed",
                    payload={
                        "owner_task_id": int(session.context.dist_ctx.task_id),
                        "pdb_id": _combo_pdb_id_from_file(pdb_file),
                        "variant": str(session.context.label).upper(),
                        "ph": _normalize_ph_tag_token(ph_tag),
                        "error": "owner_completed_without_receptor_artifact",
                        "prep_attempt": int(prep_attempt_no),
                    },
                )
        with session.failure_lock:
            add_assigned_pdb_id(
                session.context.distributed_assigned_pdb_ids,
                pdb_file=pdb_file,
            )
        if ok:
            try:
                softfail_count = int(combo_cfg.get("_SOFTFAIL_LIGAND_COUNT", 0) or 0)
            except Exception:
                softfail_count = 0
            softfail_reasons_raw = combo_cfg.get("_SOFTFAIL_LIGAND_REASONS")
            softfail_reasons: list[str] = (
                [str(x) for x in softfail_reasons_raw if str(x).strip()]
                if isinstance(softfail_reasons_raw, list)
                else []
            )
            degraded = bool(softfail_count > 0)
            if degraded:
                logging.warning(
                    "[distributed.chunk.degraded] variant=%s chunk_id=%s pdb=%s ph=%s softfail_ligands=%d reasons=%s",
                    session.context.label.upper(),
                    claimed_id,
                    pdb_file,
                    _normalize_ph_tag_token(ph_tag),
                    softfail_count,
                    ",".join(softfail_reasons) or "none",
                )
            write_chunk_result_and_cache(session, 
                chunk_id=claimed_id,
                status="completed",
                payload={
                    "pdb_file": pdb_file,
                    "ph_tag": ph_tag,
                    "run_mode": chunk_run_mode,
                    "library_name": chunk_library_name,
                    "library_root": chunk_library_root,
                    "ligand_count": len(chunk_ligands),
                    "attempt": int(attempt_no),
                    "max_attempts": int(session.cfg.chunk_max_attempts),
                    "degraded": degraded,
                    "degraded_reason": (
                        "ligand_soft_failures_present" if degraded else ""
                    ),
                    "softfail_ligand_count": int(softfail_count),
                    "softfail_reasons": softfail_reasons,
                },
            )
            session.completion_bookkeeper.mark_chunk_completed(
                combo_scope,
                reason="chunk_completed",
            )
        else:
            chunk_reason = str(
                combo_cfg.get(_CHUNK_FAILURE_REASON_KEY) or ""
            ).strip()
            chunk_retryable = bool(combo_cfg.get(_CHUNK_FAILURE_RETRYABLE_KEY))
            chunk_error = (
                "chunk_output_verification_failed"
                if chunk_reason.startswith("chunk_output_verification_failed:")
                else "process_one_failed"
            )
            fail_status = (
                "terminal_failed"
                if int(attempt_no) >= int(session.cfg.chunk_max_attempts)
                else "failed"
            )
            write_chunk_result_and_cache(session, 
                chunk_id=claimed_id,
                status=fail_status,
                payload={
                    "pdb_file": pdb_file,
                    "ph_tag": ph_tag,
                    "run_mode": chunk_run_mode,
                    "library_name": chunk_library_name,
                    "library_root": chunk_library_root,
                    "ligand_count": len(chunk_ligands),
                    "error": chunk_error,
                    "chunk_reason": chunk_reason,
                    "chunk_retryable": bool(chunk_retryable),
                    "attempt": int(attempt_no),
                    "max_attempts": int(session.cfg.chunk_max_attempts),
                },
            )
            if fail_status == "terminal_failed":
                session.completion_bookkeeper.mark_chunk_terminal_failure(
                    combo_scope,
                    reason="chunk_terminal_failed",
                )
                with session.failure_lock:
                    if claimed_id not in session.context.recorded_terminal_chunk_ids:
                        session.context.recorded_terminal_chunk_ids.add(claimed_id)
                        session.context.failed_entries.append(
                            (
                                os.path.splitext(os.path.basename(pdb_file))[0].upper(),
                                session.context.label,
                                "-",
                                "ChunkFailed",
                                (
                                    f"chunk_id={claimed_id} attempts={attempt_no} "
                                    f"reason={chunk_reason or chunk_error}"
                                ),
                            )
                        )
    except Exception as exc:
        if prep_owner:
            try:
                if combo_receptor_ready_after_settle(
                    session.context.cfg_v,
                    pdb_file=pdb_file,
                    variant_token=session.context.variant,
                    ph_tag=ph_tag,
                    poll_sec=session.cfg.combo_prep_wait_poll_sec,
                ):
                    write_combo_prep_state(
                        session.context.dist_ctx,
                        session.context.cfg_v,
                        combo_key=scope_prep_key,
                        status="ready",
                        payload=build_combo_prep_ready_payload(
                            cfg_source=combo_cfg,
                            owner_task_id=int(session.context.dist_ctx.task_id),
                            pdb_file=pdb_file,
                            variant_label=session.context.label,
                            ph_tag=ph_tag,
                            note="owner_exception_receptor_ready",
                            error=str(exc),
                        ),
                    )
                else:
                    write_combo_prep_state(
                        session.context.dist_ctx,
                        session.context.cfg_v,
                        combo_key=scope_prep_key,
                        status="failed",
                        payload={
                            "owner_task_id": int(session.context.dist_ctx.task_id),
                            "pdb_id": _combo_pdb_id_from_file(pdb_file),
                            "variant": str(session.context.label).upper(),
                            "ph": _normalize_ph_tag_token(ph_tag),
                            "error": str(exc),
                        },
                    )
            except Exception:
                logging.warning(
                    "[prep.state.write] scope=%s combo=%s action=skip reason=unexpected_exception",
                    scope_prep_key,
                    combo_key,
                    exc_info=True,
                )
        exc_text = str(exc)
        exc_text_l = exc_text.lower()
        soft_tokens = (
            "config_guard_violation",
            "postprocess_exception",
            "retry_validation_exception",
            "pose_file_missing_transient",
        )
        is_soft_exception = any(token in exc_text_l for token in soft_tokens)
        if is_soft_exception:
            logging.warning(
                "[distributed.chunk.degraded] variant=%s chunk_id=%s pdb=%s ph=%s reason=ligand_soft_exception error=%s",
                session.context.label.upper(),
                claimed_id,
                pdb_file,
                _normalize_ph_tag_token(ph_tag),
                exc_text,
            )
            write_chunk_result_and_cache(session, 
                chunk_id=claimed_id,
                status="completed",
                payload={
                    "pdb_file": pdb_file,
                    "ph_tag": ph_tag,
                    "ligand_count": len(chunk_ligands),
                    "attempt": int(attempt_no),
                    "max_attempts": int(session.cfg.chunk_max_attempts),
                    "degraded": True,
                    "degraded_reason": "ligand_soft_exception",
                    "error": exc_text,
                },
            )
            session.completion_bookkeeper.mark_chunk_completed(
                combo_scope,
                reason="chunk_completed_degraded_exception",
            )
        else:
            fail_status = (
                "terminal_failed"
                if int(attempt_no) >= int(session.cfg.chunk_max_attempts)
                else "failed"
            )
            write_chunk_result_and_cache(session, 
                chunk_id=claimed_id,
                status=fail_status,
                payload={
                    "pdb_file": pdb_file,
                    "ph_tag": ph_tag,
                    "ligand_count": len(chunk_ligands),
                    "error": str(exc),
                    "attempt": int(attempt_no),
                    "max_attempts": int(session.cfg.chunk_max_attempts),
                },
            )
            if fail_status == "terminal_failed":
                session.completion_bookkeeper.mark_chunk_terminal_failure(
                    combo_scope,
                    reason="chunk_exception_terminal",
                )
                with session.failure_lock:
                    if claimed_id not in session.context.recorded_terminal_chunk_ids:
                        session.context.recorded_terminal_chunk_ids.add(claimed_id)
                        session.context.failed_entries.append(
                            (
                                os.path.splitext(os.path.basename(pdb_file))[0].upper(),
                                session.context.label,
                                "-",
                                "ChunkException",
                                f"chunk_id={claimed_id} attempts={attempt_no} error={exc}",
                            )
                        )
            logging.exception(
                "[distributed.combo-chunk.error] variant=%s chunk_id=%s pdb=%s ph=%s",
                session.context.label.upper(),
                claimed_id,
                pdb_file,
                _normalize_ph_tag_token(ph_tag),
            )
    finally:
        heartbeat_stop.set()
        if heartbeat_thread is not None:
            heartbeat_thread.join(timeout=1.0)
        with session.progress_lock:
            session.context.bar.update(1)
