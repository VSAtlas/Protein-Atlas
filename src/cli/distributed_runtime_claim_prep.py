from __future__ import annotations

import hashlib
import logging
import os
import time
from typing import Any, Mapping, Optional

from cli.distributed_chunk_planner import (
    _combo_pdb_id_from_file,
    _combo_receptor_ready,
    _normalize_ph_tag_token,
    _scope_prep_key,
)
from cli.distributed_chunk_runtime_state import (
    combo_prep_has_reuse_data,
    combo_prep_status,
)
from cli.distributed_context import combo_prep_claim_dir, read_combo_prep_state
from cli.distributed_runtime_claim_types import (
    ClaimLoopConfig,
    ClaimLoopContext,
    ClaimLoopSession,
)
from cli.distributed_runtime_completion import ComboScope


def combo_receptor_ready_after_settle(
    cfg: Mapping[str, Any],
    *,
    pdb_file: str,
    variant_token: Optional[str],
    ph_tag: Optional[str],
    poll_sec: float,
    timeout_sec: float = 20.0,
) -> bool:
    if _combo_receptor_ready(
        cfg,
        pdb_file=pdb_file,
        variant_token=variant_token,
        ph_tag=ph_tag,
    ):
        return True
    deadline = time.time() + max(0.0, float(timeout_sec))
    sleep_for = max(0.25, min(2.0, float(poll_sec)))
    while time.time() < deadline:
        time.sleep(sleep_for)
        if _combo_receptor_ready(
            cfg,
            pdb_file=pdb_file,
            variant_token=variant_token,
            ph_tag=ph_tag,
        ):
            return True
    return False


def scope_prep_lock_token(*, run_id: str, variant_label: str, pdb_file: str) -> str:
    raw = "|".join(
        (
            str(run_id),
            str(variant_label or "").strip().upper() or "BASE",
            _combo_pdb_id_from_file(pdb_file),
        )
    )
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:24]


def acquire_scope_prep_lock(
    *,
    context: ClaimLoopContext,
    cfg: ClaimLoopConfig,
    pdb_file: str,
) -> Any:
    if not bool(getattr(context.dist_ctx, "enabled", False)):
        return None
    token = scope_prep_lock_token(
        run_id=context.run_id,
        variant_label=context.label,
        pdb_file=pdb_file,
    )
    lock_path = combo_prep_claim_dir(context.dist_ctx, context.cfg_v) / f"scope_{token}.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        "{"
        f"\"run_id\":\"{context.run_id}\","
        f"\"task_id\":{int(context.dist_ctx.task_id)},"
        f"\"pdb_id\":\"{_combo_pdb_id_from_file(pdb_file)}\","
        f"\"variant\":\"{str(context.label).upper()}\","
        f"\"claimed_at\":{time.time():.6f}"
        "}\n"
    )
    wait_logged = False
    while True:
        try:
            fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            try:
                age = max(0.0, time.time() - float(lock_path.stat().st_mtime))
                if age > float(max(30.0, cfg.combo_prep_lease_sec)):
                    lock_path.unlink()
                    continue
            except FileNotFoundError:
                continue
            except Exception:
                pass
            if not wait_logged:
                wait_logged = True
                logging.info(
                    "[prep.scope.wait] pdb=%s variant=%s task=%d lock=%s",
                    _combo_pdb_id_from_file(pdb_file),
                    str(context.label).upper(),
                    int(context.dist_ctx.task_id),
                    lock_path,
                )
            time.sleep(max(0.25, min(2.0, float(cfg.combo_prep_wait_poll_sec))))
            continue
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        logging.info(
            "[prep.scope.owner] pdb=%s variant=%s task=%d lock=%s",
            _combo_pdb_id_from_file(pdb_file),
            str(context.label).upper(),
            int(context.dist_ctx.task_id),
            lock_path,
        )
        return lock_path


def release_scope_prep_lock(lock_path: Any) -> None:
    if lock_path is None:
        return
    try:
        lock_path.unlink()
    except FileNotFoundError:
        pass
    except Exception:
        logging.warning("[prep.scope.release] action=skip reason=unlink_failed", exc_info=True)


def collect_prep_ready_scopes(session: ClaimLoopSession) -> set[ComboScope]:
    context = session.context
    ready: set[ComboScope] = set()
    seen: set[ComboScope] = set()
    for chunk in context.by_chunk_id.values():
        if not isinstance(chunk, dict):
            continue
        pdb_file = str(chunk.get("pdb_file") or "")
        ph_tag_raw = chunk.get("ph_tag")
        ph_tag = None if ph_tag_raw is None else str(ph_tag_raw)
        combo_scope: ComboScope = (
            _combo_pdb_id_from_file(pdb_file),
            str(context.label or "").strip().upper() or "BASE",
            _normalize_ph_tag_token(ph_tag),
        )
        if combo_scope in seen:
            continue
        seen.add(combo_scope)
        scope_prep_key = _scope_prep_key(
            run_id=context.run_id,
            variant_label=context.label,
            pdb_file=pdb_file,
            ph_tag=ph_tag,
        )
        prep_status, prep_payload = combo_prep_status(
            dist_ctx=context.dist_ctx,
            cfg_v=context.cfg_v,
            combo_key=scope_prep_key,
            read_combo_prep_state=read_combo_prep_state,
        )
        status = str(prep_status or "").strip().lower()
        if status == "ready":
            ready.add(combo_scope)
            continue
        preparing_since = None
        try:
            preparing_since = float((prep_payload or {}).get("updated_at") or 0.0)
        except Exception:
            preparing_since = None
        if (
            status == "preparing"
            and combo_prep_has_reuse_data(prep_payload)
            and _combo_receptor_ready(
                context.cfg_v,
                pdb_file=pdb_file,
                variant_token=context.variant,
                ph_tag=ph_tag,
                min_mtime=preparing_since,
            )
        ):
            ready.add(combo_scope)
    return ready


def chunk_prep_claimable(session: ClaimLoopSession, candidate_id: str) -> bool:
    context = session.context
    cfg = session.cfg
    chunk = context.by_chunk_id.get(str(candidate_id))
    if not isinstance(chunk, dict):
        return True
    pdb_file = str(chunk.get("pdb_file") or "")
    ph_tag_raw = chunk.get("ph_tag")
    ph_tag = None if ph_tag_raw is None else str(ph_tag_raw)
    combo_scope: ComboScope = (
        _combo_pdb_id_from_file(pdb_file),
        str(context.label or "").strip().upper() or "BASE",
        _normalize_ph_tag_token(ph_tag),
    )
    owner_task = int(
        context.combo_owner_task_id.get(combo_scope, context.dist_ctx.task_id)
    )
    is_owner = int(owner_task) == int(context.dist_ctx.task_id)
    scope_prep_key = _scope_prep_key(
        run_id=context.run_id,
        variant_label=context.label,
        pdb_file=pdb_file,
        ph_tag=ph_tag,
    )
    prep_status, prep_payload = combo_prep_status(
        dist_ctx=context.dist_ctx,
        cfg_v=context.cfg_v,
        combo_key=scope_prep_key,
        read_combo_prep_state=read_combo_prep_state,
    )
    status = str(prep_status or "").strip().lower()
    if status in {"ready", "failed"}:
        return True

    preparing_since = None
    try:
        preparing_since = float((prep_payload or {}).get("updated_at") or 0.0)
    except Exception:
        preparing_since = None
    if (
        status == "preparing"
        and combo_prep_has_reuse_data(prep_payload)
        and _combo_receptor_ready(
            context.cfg_v,
            pdb_file=pdb_file,
            variant_token=context.variant,
            ph_tag=ph_tag,
            min_mtime=preparing_since,
        )
    ):
        return True
    if status == "preparing":
        if preparing_since and (
            time.time() - float(preparing_since)
        ) >= float(cfg.combo_prep_wait_timeout_sec):
            return bool(is_owner)
        return False
    if not is_owner:
        return False
    return True
