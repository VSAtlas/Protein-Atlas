"""Deferred queue helpers for pose-capture rendering."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from docking import capture_pose_common as common

def set_defer_mode(enable: bool, queue_path: Optional[str] = None):
    """Toggle deferral globally (preferred entry point for callers)."""
    common._DEFER_MODE = bool(enable)
    if queue_path:
        common._QUEUE_PATH = str(queue_path)
    return
    

def _enqueue_job(kind: str, payload: dict):
    Path(common._QUEUE_PATH).parent.mkdir(parents=True, exist_ok=True)
    rec = {"kind": kind, **payload}
    line = json.dumps(rec, ensure_ascii=False)
    with common._Q_LOCK:
        with open(common._QUEUE_PATH, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")

def _iter_jobs(path: str):
    p = Path(path)
    if not p.is_file():
        return
    with open(p, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            try:
                yield json.loads(ln)
            except Exception:
                continue
