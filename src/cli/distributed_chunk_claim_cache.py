from __future__ import annotations

import time
from dataclasses import dataclass
from threading import RLock
from typing import Any, Callable, Mapping, Sequence


@dataclass
class ChunkResultSnapshot:
    status: str
    attempt: int


class DistributedChunkClaimCache:
    def __init__(
        self,
        *,
        candidate_ids: Sequence[str],
        max_attempts: int,
    ) -> None:
        ordered = [str(x).strip() for x in candidate_ids if str(x).strip()]
        self._candidate_ids = ordered
        self._cursor = 0
        self._max_attempts = max(1, int(max_attempts))
        self._snapshot: dict[str, ChunkResultSnapshot] = {
            chunk_id: ChunkResultSnapshot(status="", attempt=0)
            for chunk_id in ordered
        }
        self._claimable_ids: list[str] = list(ordered)
        self._claimable_cursor = 0
        self._claimable_dirty = False
        self._last_reconcile = 0.0
        self._lock = RLock()

    @property
    def candidate_count(self) -> int:
        return int(len(self._candidate_ids))

    def _is_terminal(self, status: str, attempt: int) -> bool:
        st = str(status or "").strip().lower()
        if st in {"completed", "terminal_failed"}:
            return True
        if st == "failed" and int(attempt) >= int(self._max_attempts):
            return True
        return False

    def unresolved_count(self) -> int:
        with self._lock:
            unresolved = 0
            for snap in self._snapshot.values():
                if not self._is_terminal(snap.status, snap.attempt):
                    unresolved += 1
            return int(unresolved)

    def terminal_counts(self) -> tuple[int, int]:
        with self._lock:
            completed = 0
            terminal_failed = 0
            for snap in self._snapshot.values():
                if not self._is_terminal(snap.status, snap.attempt):
                    continue
                st = str(snap.status or "").strip().lower()
                if st == "completed":
                    completed += 1
                else:
                    terminal_failed += 1
            return int(completed), int(terminal_failed)

    def update_from_payload(self, *, chunk_id: str, payload: Mapping[str, Any] | None) -> None:
        cid = str(chunk_id or "").strip()
        if not cid:
            return
        with self._lock:
            if cid not in self._snapshot:
                self._snapshot[cid] = ChunkResultSnapshot(status="", attempt=0)
            if not isinstance(payload, Mapping):
                self._snapshot[cid] = ChunkResultSnapshot(status="", attempt=0)
                self._claimable_dirty = True
                return
            status = str(payload.get("status") or "").strip().lower()
            try:
                attempt = int(payload.get("attempt") or 1)
            except Exception:
                attempt = 1
            self._snapshot[cid] = ChunkResultSnapshot(
                status=status,
                attempt=max(0, int(attempt)),
            )
            self._claimable_dirty = True

    def update_local_result(self, *, chunk_id: str, status: str, attempt: int) -> None:
        cid = str(chunk_id or "").strip()
        if not cid:
            return
        with self._lock:
            self._snapshot[cid] = ChunkResultSnapshot(
                status=str(status or "").strip().lower(),
                attempt=max(0, int(attempt)),
            )
            self._claimable_dirty = True

    def _refresh_claimable_ids(self) -> None:
        if not self._claimable_dirty:
            return
        self._claimable_ids = [
            cid
            for cid, snap in self._snapshot.items()
            if not self._is_terminal(snap.status, snap.attempt)
        ]
        if self._claimable_cursor >= len(self._claimable_ids):
            self._claimable_cursor = 0
        self._claimable_dirty = False

    def next_claim_candidate(self) -> tuple[str | None, int]:
        with self._lock:
            self._refresh_claimable_ids()
            total = len(self._claimable_ids)
            if total <= 0:
                return None, 0
            for _ in range(total):
                cid = self._claimable_ids[self._claimable_cursor % total]
                self._claimable_cursor += 1
                snap = self._snapshot.get(cid) or ChunkResultSnapshot(status="", attempt=0)
                if self._is_terminal(snap.status, snap.attempt):
                    continue
                return cid, int(max(0, snap.attempt))
            return None, 0

    def should_reconcile(self, *, now: float, interval_sec: float, force: bool = False) -> bool:
        if force:
            return True
        return (float(now) - float(self._last_reconcile)) >= float(max(0.1, interval_sec))

    def reconcile(
        self,
        *,
        read_result: Callable[[str], Mapping[str, Any] | None],
        now: float | None = None,
    ) -> tuple[int, int, int]:
        with self._lock:
            scan_ids = [
                cid
                for cid, snap in self._snapshot.items()
                if not self._is_terminal(snap.status, snap.attempt)
            ]
        for cid in scan_ids:
            payload = read_result(str(cid))
            self.update_from_payload(chunk_id=str(cid), payload=payload)
        with self._lock:
            self._refresh_claimable_ids()
            stamp = float(now if now is not None else time.time())
            self._last_reconcile = stamp
            unresolved = self.unresolved_count()
            completed, terminal_failed = self.terminal_counts()
            return int(unresolved), int(completed), int(terminal_failed)
