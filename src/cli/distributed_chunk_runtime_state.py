from __future__ import annotations

import os
import threading
import time
from typing import Any, Callable, Mapping, Optional, Sequence

from cli.distributed_chunk_claim_cache import DistributedChunkClaimCache
from cli.distributed_chunk_planner import (
    _CHUNK_MAX_ATTEMPTS,
    _CHUNK_STATE_RECONCILE_SEC,
    _combo_pdb_id_from_file,
    _normalize_ph_tag_token,
)


def chunk_terminal(status: str, attempt: int, *, max_attempts: int) -> bool:
    if status in {"completed", "terminal_failed"}:
        return True
    if status == "failed" and int(attempt) >= int(max_attempts):
        return True
    return False


class ChunkStateTracker:
    def __init__(
        self,
        *,
        chunk_ids: Sequence[str],
        dist_ctx: Any,
        cfg_v: Mapping[str, Any],
        read_chunk_result: Callable[..., Any],
        max_attempts: int = _CHUNK_MAX_ATTEMPTS,
        reconcile_sec: float = _CHUNK_STATE_RECONCILE_SEC,
        use_claim_cache: bool,
        claim_cache: Optional[DistributedChunkClaimCache] = None,
    ) -> None:
        self._chunk_ids = [str(x) for x in chunk_ids if str(x).strip()]
        self._dist_ctx = dist_ctx
        self._cfg_v = cfg_v
        self._read_chunk_result = read_chunk_result
        self._max_attempts = int(max_attempts)
        self._reconcile_sec = float(reconcile_sec)
        self._use_claim_cache = bool(use_claim_cache)
        self._claim_cache = claim_cache
        self._lock = threading.Lock()
        self._last_refresh = 0.0
        self._unresolved = max(0, len(self._chunk_ids))
        self._completed = 0
        self._terminal_failed = 0

    def set_claim_cache(self, claim_cache: Optional[DistributedChunkClaimCache]) -> None:
        self._claim_cache = claim_cache

    def record_local_result(self, *, chunk_id: str, status: str, attempt: int) -> None:
        if self._use_claim_cache and self._claim_cache is not None:
            self._claim_cache.update_local_result(
                chunk_id=str(chunk_id),
                status=str(status),
                attempt=int(attempt),
            )

    def refresh(self, *, force: bool) -> tuple[int, int, int]:
        now = time.time()
        with self._lock:
            should_refresh = bool(
                force
                or self._last_refresh <= 0.0
                or (now - float(self._last_refresh) >= float(self._reconcile_sec))
            )
            if not should_refresh:
                return (
                    int(self._unresolved),
                    int(self._completed),
                    int(self._terminal_failed),
                )

        unresolved = 0
        completed = 0
        terminal_failed = 0
        if self._use_claim_cache and self._claim_cache is not None:
            unresolved, completed, terminal_failed = self._claim_cache.reconcile(
                read_result=lambda cid: self._read_chunk_result(
                    self._dist_ctx,
                    self._cfg_v,
                    chunk_id=cid,
                ),
                now=now,
            )
        else:
            for chunk_id in self._chunk_ids:
                payload = self._read_chunk_result(
                    self._dist_ctx,
                    self._cfg_v,
                    chunk_id=chunk_id,
                )
                if not isinstance(payload, dict):
                    unresolved += 1
                    continue
                status = str(payload.get("status") or "").strip().lower()
                try:
                    attempt = int(payload.get("attempt") or 1)
                except Exception:
                    attempt = 1
                if not chunk_terminal(
                    status,
                    attempt,
                    max_attempts=self._max_attempts,
                ):
                    unresolved += 1
                    continue
                if status == "completed":
                    completed += 1
                else:
                    terminal_failed += 1

        with self._lock:
            self._last_refresh = now
            self._unresolved = int(unresolved)
            self._completed = int(completed)
            self._terminal_failed = int(terminal_failed)
            return (
                int(self._unresolved),
                int(self._completed),
                int(self._terminal_failed),
            )

    def all_done(self) -> bool:
        unresolved, _, _ = self.refresh(force=False)
        return int(unresolved) <= 0

    def unresolved_chunk_count(self) -> int:
        unresolved, _, _ = self.refresh(force=False)
        return int(unresolved)


def combo_prep_status(
    *,
    dist_ctx: Any,
    cfg_v: Mapping[str, Any],
    combo_key: str,
    read_combo_prep_state: Callable[..., Any],
) -> tuple[str, dict[str, Any] | None]:
    payload = read_combo_prep_state(dist_ctx, cfg_v, combo_key=combo_key)
    if not isinstance(payload, dict):
        return "", None
    return str(payload.get("status") or "").strip().lower(), payload


def combo_prep_has_reuse_data(payload: Mapping[str, Any] | None) -> bool:
    if not isinstance(payload, dict):
        return False
    center_map = payload.get("center_by_ph")
    box_map = payload.get("box_by_ph")
    return bool(
        isinstance(center_map, dict)
        and center_map
        and isinstance(box_map, dict)
        and box_map
    )


def apply_combo_prep_reuse_payload(
    cfg_target: dict[str, Any],
    payload: Mapping[str, Any] | None,
) -> bool:
    if not combo_prep_has_reuse_data(payload):
        return False
    assert isinstance(payload, dict)
    cfg_target["_COMBO_PREP_CENTER_BY_PH"] = payload.get("center_by_ph")
    cfg_target["_COMBO_PREP_BOX_BY_PH"] = payload.get("box_by_ph")
    source_map = payload.get("center_source_by_ph")
    if isinstance(source_map, dict) and source_map:
        cfg_target["_COMBO_PREP_SOURCE_BY_PH"] = source_map
    control_stems = payload.get("control_stems")
    if isinstance(control_stems, list) and control_stems:
        cfg_target["_COMBO_PREP_CONTROL_STEMS"] = [
            str(x) for x in control_stems if str(x).strip()
        ]
    control_lookup = payload.get("control_lookup")
    if isinstance(control_lookup, dict) and control_lookup:
        cfg_target["_COMBO_PREP_CONTROL_LOOKUP"] = {
            str(k): str(v)
            for k, v in control_lookup.items()
            if str(k).strip() and str(v).strip()
        }
    return True


def build_combo_prep_ready_payload(
    *,
    cfg_source: Mapping[str, Any],
    owner_task_id: int,
    pdb_file: str,
    variant_label: str,
    ph_tag: Optional[str],
    note: Optional[str] = None,
    error: Optional[str] = None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "owner_task_id": int(owner_task_id),
        "pdb_id": _combo_pdb_id_from_file(pdb_file),
        "variant": str(variant_label).upper(),
        "ph": _normalize_ph_tag_token(ph_tag),
    }
    if note:
        payload["note"] = str(note)
    if error:
        payload["error"] = str(error)
    center_map = cfg_source.get("_COMBO_PREP_CENTER_BY_PH")
    box_map = cfg_source.get("_COMBO_PREP_BOX_BY_PH")
    if isinstance(center_map, dict) and isinstance(box_map, dict):
        if center_map and box_map:
            payload["center_by_ph"] = dict(center_map)
            payload["box_by_ph"] = dict(box_map)
    source_map = cfg_source.get("_COMBO_PREP_SOURCE_BY_PH")
    if isinstance(source_map, dict) and source_map:
        payload["center_source_by_ph"] = dict(source_map)
    control_stems = cfg_source.get("_COMBO_PREP_CONTROL_STEMS")
    if isinstance(control_stems, list) and control_stems:
        payload["control_stems"] = [str(x) for x in control_stems if str(x).strip()]
    control_lookup = cfg_source.get("_COMBO_PREP_CONTROL_LOOKUP")
    if isinstance(control_lookup, dict) and control_lookup:
        payload["control_lookup"] = {
            str(k): str(v)
            for k, v in control_lookup.items()
            if str(k).strip() and str(v).strip()
        }
    return payload


def add_assigned_pdb_id(
    assigned_ids: set[str],
    *,
    pdb_file: str,
) -> None:
    assigned_ids.add(os.path.splitext(os.path.basename(pdb_file))[0].upper())
