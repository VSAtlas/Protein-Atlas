from __future__ import annotations

import json
import threading
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator


class EventLogger:
    """Low-overhead line-buffered JSONL event logger."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = self.path.open("a", encoding="utf-8", buffering=1)
        self._lock = threading.Lock()

    def emit(self, event: str, **fields: Any) -> None:
        payload = {
            "ts_mono_ns": time.monotonic_ns(),
            "ts_wall_ns": time.time_ns(),
            "event": str(event),
        }
        payload.update(fields)
        line = json.dumps(payload, separators=(",", ":"))
        with self._lock:
            self._fh.write(line + "\n")

    def close(self) -> None:
        with self._lock:
            if not self._fh.closed:
                self._fh.flush()
                self._fh.close()


_EVENT_LOGGER: EventLogger | None = None


def set_event_logger(logger: EventLogger | None) -> None:
    global _EVENT_LOGGER
    _EVENT_LOGGER = logger


def get_event_logger() -> EventLogger | None:
    return _EVENT_LOGGER


def emit_event(event: str, **fields: Any) -> None:
    logger = get_event_logger()
    if logger is None:
        return
    logger.emit(event, **fields)


def close_event_log() -> None:
    logger = get_event_logger()
    if logger is None:
        return
    try:
        logger.close()
    finally:
        set_event_logger(None)


@contextmanager
def event_run(path: str | Path, run_meta: dict[str, Any] | None = None) -> Iterator[EventLogger]:
    logger = EventLogger(path)
    set_event_logger(logger)
    meta = dict(run_meta or {})
    try:
        logger.emit("RUN_START", **meta)
        yield logger
        logger.emit("RUN_END", **meta)
    finally:
        close_event_log()
