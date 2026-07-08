from __future__ import annotations

from typing import Any, Iterable

from analysis._common import parse_binary, parse_float


def numeric_value(value: Any) -> float | None:
    binary = parse_binary(value)
    if binary is not None:
        return float(binary)
    return parse_float(value)


def bounds(values: Iterable[float | None]) -> tuple[float, float] | None:
    vals = [float(v) for v in values if v is not None]
    if not vals:
        return None
    return min(vals), max(vals)


def normalize(value: float | None, value_bounds: tuple[float, float] | None) -> float | None:
    if value is None or value_bounds is None:
        return None
    lo, hi = value_bounds
    if hi == lo:
        return 0.5
    return max(0.0, min(1.0, (value - lo) / (hi - lo)))

