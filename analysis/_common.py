from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from path_router.path_router import _ensure_router_roots


TRUTHY = {"1", "true", "yes", "y", "on", "supported", "positive", "hit"}
FALSY = {"0", "false", "no", "n", "off", "unsupported", "negative"}


def repo_root() -> Path:
    return _ensure_router_roots().overall


def clean_text(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def first_value(row: Mapping[str, Any], names: Sequence[str]) -> str:
    for name in names:
        value = clean_text(row.get(name))
        if value:
            return value
    return ""


def parse_float(value: Any) -> float | None:
    text = clean_text(value)
    if not text:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def format_float(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.6g}"


def parse_binary(value: Any) -> int | None:
    text = clean_text(value).lower()
    if not text:
        return None
    if text in TRUTHY:
        return 1
    if text in FALSY:
        return 0
    numeric = parse_float(text)
    if numeric is None:
        return None
    return 1 if numeric > 0 else 0


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def write_csv_rows(path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames or [])
    if not names:
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    names.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in names})


def entity_key(row: Mapping[str, Any], *names: str) -> str:
    for name in names:
        value = clean_text(row.get(name))
        if value:
            return value
    return ""


def minmax(values: Iterable[float]) -> tuple[float, float] | None:
    vals = [v for v in values if math.isfinite(v)]
    if not vals:
        return None
    lo = min(vals)
    hi = max(vals)
    return lo, hi


def scale01(value: float | None, bounds: tuple[float, float] | None) -> float | None:
    if value is None or bounds is None:
        return None
    lo, hi = bounds
    if hi == lo:
        return 0.5
    return max(0.0, min(1.0, (value - lo) / (hi - lo)))

