from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

LOG_LEVELS = {"DEBUG": 10, "INFO": 20, "WARN": 30, "ERROR": 40}

LEVEL_ABBREV = {"DEBUG": "DBG", "INFO": "INF", "WARN": "WRN", "ERROR": "ERR"}

_current_log_level = LOG_LEVELS["INFO"]

_BACKEND_LOGGED = False

def set_log_level(level_name: Optional[str]) -> None:
    global _current_log_level
    key = (level_name or "INFO").strip().upper()
    if key not in LOG_LEVELS:
        key = "INFO"
    _current_log_level = LOG_LEVELS[key]

def dbg(level: str, tag: str, message: str) -> None:
    upper = (level or "INFO").strip().upper()
    if upper not in LOG_LEVELS:
        upper = "INFO"
    if LOG_LEVELS[upper] < _current_log_level:
        return
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    label = LEVEL_ABBREV.get(upper, upper[:3])
    print(f"[{stamp}][{label}][{tag}] {message}")
