from __future__ import annotations

import datetime as _dt
import html
import json
import time
from pathlib import Path
from typing import Any, Mapping

from config.output_paths import output_root


STAGE_ORDER = ("prep", "pocket_detection", "docking", "postprocessing")


def mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def esc(value: str) -> str:
    return html.escape(value, quote=True)


def as_int(raw: Any) -> int | None:
    try:
        return int(raw)
    except Exception:
        return None


def as_float(raw: Any) -> float | None:
    try:
        return float(raw)
    except Exception:
        return None


def format_seconds(seconds: float) -> str:
    seconds_i = max(0, int(round(seconds)))
    hours, rem = divmod(seconds_i, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h{minutes:02d}m"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def median(values: list[float]) -> float | None:
    if not values:
        return None
    vals = sorted(float(v) for v in values)
    mid = len(vals) // 2
    if len(vals) % 2:
        return round(vals[mid], 3)
    return round((vals[mid - 1] + vals[mid]) / 2.0, 3)


def parse_time(raw: Any) -> float | None:
    if not raw:
        return None
    text = str(raw).strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        return _dt.datetime.fromisoformat(text).timestamp()
    except Exception:
        return None


def elapsed_seconds(timing: Mapping[str, Any]) -> float | None:
    wall = as_float(timing.get("wall_time_sec"))
    if wall is not None:
        return float(wall)
    started = parse_time(timing.get("started_at") or timing.get("created_at"))
    finished = parse_time(timing.get("finished_at"))
    if started is None:
        return None
    end = finished or time.time()
    return max(0.0, float(end - started))


def read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def chunk_ids_from_json(path: Path) -> set[str]:
    payload = read_json(path)
    ids: set[str] = set()

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            chunk_id = value.get("chunk_id")
            if chunk_id:
                ids.add(str(chunk_id))
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(payload)
    return ids


def run_candidates(root: Path) -> list[tuple[float, str]]:
    candidates: dict[str, float] = {}
    for base in (
        output_root(root, "manifests"),
        output_root(root, "data"),
        output_root(root, "post_docked"),
        output_root(root, "docked"),
        root / "manifests",
        root / "data",
        root / "post_docked",
        root / "docked",
    ):
        if not base.exists():
            continue
        for child in base.iterdir():
            if not child.is_dir():
                continue
            try:
                mtime = child.stat().st_mtime
            except OSError:
                continue
            candidates[child.name] = max(candidates.get(child.name, 0.0), float(mtime))
    return sorted(((mtime, run_id) for run_id, mtime in candidates.items()), reverse=True)


def path_age_seconds(path: Path) -> float:
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return 0.0
    return max(0.0, time.time() - float(mtime))


def eta_label(eta: Mapping[str, Any]) -> str:
    eta_sec = eta.get("eta_sec")
    if isinstance(eta_sec, (int, float)):
        return format_seconds(float(eta_sec))
    return "unknown"
