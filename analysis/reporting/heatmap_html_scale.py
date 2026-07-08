from __future__ import annotations

import logging
import math
from typing import Any, Dict, Iterable, List, Optional, Tuple

from analysis.reporting.value_utils import as_float

_LOG = logging.getLogger("heatmap-html")
_HEATMAP_SCALE_FIXED = "fixed"
_HEATMAP_SCALE_QUANTILE = "quantile"
_HEATMAP_SCALE_MODES = {_HEATMAP_SCALE_FIXED, _HEATMAP_SCALE_QUANTILE}

def _parse_num(value: Any, default_value: float) -> float:
    """Heatmap numeric config: finite floats only; non-finite / unparseable → default."""
    parsed = as_float(value)
    return default_value if parsed is None else float(parsed)


def _parse_color_value(value: str) -> Optional[Tuple[int, int, int]]:
    if not value:
        return None
    candidate = value.strip()
    if not candidate:
        return None
    try:
        from PIL import ImageColor

        rgb = ImageColor.getrgb(candidate)
        if len(rgb) == 4:
            rgb = rgb[:3]
        return (rgb[0], rgb[1], rgb[2])
    except Exception:
        pass
    if candidate.startswith("#"):
        hex_str = candidate[1:]
        if len(hex_str) == 3 and all(ch in "0123456789abcdefABCDEF" for ch in hex_str):
            return tuple(int(ch * 2, 16) for ch in hex_str)  # type: ignore[return-value]
        if len(hex_str) == 6 and all(ch in "0123456789abcdefABCDEF" for ch in hex_str):
            return (
                int(hex_str[0:2], 16),
                int(hex_str[2:4], 16),
                int(hex_str[4:6], 16),
            )
    return None


def _resolve_color(value: Optional[str], default_value: str) -> Tuple[int, int, int]:
    candidate = (value or "").strip()
    rgb = _parse_color_value(candidate) if candidate else None
    if rgb is None:
        rgb = _parse_color_value(default_value)
    if rgb is None:
        raise ValueError(f"Unable to parse heatmap color value: {value}")
    return rgb


def _rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
    r, g, b = rgb
    return f"#{r:02x}{g:02x}{b:02x}"


def _color_gradient(
    start_rgb: Tuple[int, int, int], end_rgb: Tuple[int, int, int], n: int
) -> List[str]:
    if n <= 0:
        return []
    if n == 1:
        return [_rgb_to_hex(start_rgb)]
    colors: List[str] = []
    for i in range(n):
        t = i / (n - 1)
        r = int(round(start_rgb[0] + t * (end_rgb[0] - start_rgb[0])))
        g = int(round(start_rgb[1] + t * (end_rgb[1] - start_rgb[1])))
        b = int(round(start_rgb[2] + t * (end_rgb[2] - start_rgb[2])))
        colors.append(_rgb_to_hex((r, g, b)))
    return colors


def _seq(start: float, end: float, length: int) -> List[float]:
    if length <= 1:
        return [start]
    step = (end - start) / (length - 1)
    return [start + step * i for i in range(length)]


def _build_palette(
    scale_min: float,
    scale_mid: float,
    scale_mid2: float,
    scale_max: float,
    color_low: Tuple[int, int, int],
    color_mid: Tuple[int, int, int],
    color_mid2: Tuple[int, int, int],
    color_high: Tuple[int, int, int],
) -> Tuple[List[float], List[str]]:
    total_steps = 100
    span_low = max(scale_mid - scale_min, 0)
    span_mid = max(scale_mid2 - scale_mid, 0)
    span_high = max(scale_max - scale_mid2, 0)
    total_span = span_low + span_mid + span_high
    if total_span <= 0:
        span_low = 1
        span_mid = 1
        span_high = 1
        total_span = 3
    n_low = max(1, round(total_steps * span_low / total_span))
    n_mid = max(1, round(total_steps * span_mid / total_span))
    n_high = total_steps - n_low - n_mid
    if n_high < 1:
        n_high = 1
        if n_mid > 1:
            n_mid -= 1
        elif n_low > 1:
            n_low -= 1

    breaks_low = _seq(scale_min, scale_mid, n_low + 1)
    breaks_mid = _seq(scale_mid, scale_mid2, n_mid + 1)
    breaks_high = _seq(scale_mid2, scale_max, n_high + 1)
    heatmap_breaks = breaks_low + breaks_mid[1:] + breaks_high[1:]
    heatmap_colors = (
        _color_gradient(color_low, color_mid, n_low)
        + _color_gradient(color_mid, color_mid2, n_mid)
        + _color_gradient(color_mid2, color_high, n_high)
    )
    return heatmap_breaks, heatmap_colors


def _format_breaks_for_log(breaks: Dict[str, float]) -> str:
    return (
        f"min={breaks['scale_min']:.6g} "
        f"mid={breaks['scale_mid']:.6g} "
        f"mid2={breaks['scale_mid2']:.6g} "
        f"max={breaks['scale_max']:.6g}"
    )


def _breaks_increasing(breaks: Dict[str, float]) -> bool:
    return (
        breaks["scale_min"] < breaks["scale_mid"] < breaks["scale_mid2"] < breaks["scale_max"]
    )


def _ensure_increasing_breaks_with_epsilon(
    breaks: Dict[str, float],
) -> Optional[Dict[str, float]]:
    if _breaks_increasing(breaks):
        return dict(breaks)
    start = float(breaks["scale_min"])
    end = float(breaks["scale_max"])
    if not math.isfinite(start) or not math.isfinite(end) or end <= start:
        return None
    span = abs(end - start)
    eps = max(span * 1e-6, 1e-9)
    widened = {
        "scale_min": start,
        "scale_mid": float(breaks["scale_mid"]),
        "scale_mid2": float(breaks["scale_mid2"]),
        "scale_max": end,
    }
    if widened["scale_mid"] <= widened["scale_min"]:
        widened["scale_mid"] = widened["scale_min"] + eps
    if widened["scale_mid2"] <= widened["scale_mid"]:
        widened["scale_mid2"] = widened["scale_mid"] + eps
    if widened["scale_max"] <= widened["scale_mid2"]:
        widened["scale_max"] = widened["scale_mid2"] + eps
    if not _breaks_increasing(widened):
        return None
    return widened


def _quantile_from_sorted(sorted_values: List[float], q: float) -> float:
    if not sorted_values:
        raise ValueError("No values available for quantile")
    if q <= 0:
        return sorted_values[0]
    if q >= 1:
        return sorted_values[-1]
    pos = (len(sorted_values) - 1) * q
    low_idx = int(math.floor(pos))
    high_idx = int(math.ceil(pos))
    if low_idx == high_idx:
        return sorted_values[low_idx]
    fraction = pos - low_idx
    low_val = sorted_values[low_idx]
    high_val = sorted_values[high_idx]
    return low_val + (high_val - low_val) * fraction


def _collect_finite_values(values: Iterable[Any]) -> List[float]:
    finite_values: List[float] = []
    for value in values:
        try:
            parsed = float(value)
        except Exception:
            continue
        if math.isfinite(parsed):
            finite_values.append(parsed)
    return finite_values


def _parse_scale_mode(value: Optional[str]) -> str:
    candidate = (value or "").strip().lower()
    if candidate in _HEATMAP_SCALE_MODES:
        return candidate
    return _HEATMAP_SCALE_FIXED


def _parse_quantile_value(value: Optional[str], default_value: float) -> float:
    parsed = _parse_num(value, default_value)
    if not math.isfinite(parsed):
        return default_value
    if parsed < 0.0:
        return 0.0
    if parsed > 1.0:
        return 1.0
    return parsed


def _parse_quantile_params(config: Dict[str, Optional[str]]) -> Dict[str, float]:
    q_low_default = 0.02
    q_high_default = 0.98
    q_mid2_default = 0.85
    q_low = _parse_quantile_value(config.get("HEATMAP_QUANTILE_LOW"), q_low_default)
    q_mid2 = _parse_quantile_value(config.get("HEATMAP_QUANTILE_MID2"), q_mid2_default)
    q_high = _parse_quantile_value(config.get("HEATMAP_QUANTILE_HIGH"), q_high_default)
    if not (0.0 <= q_low < q_mid2 < q_high <= 1.0):
        q_low, q_mid2, q_high = q_low_default, q_mid2_default, q_high_default
    return {"q_low": q_low, "q_mid2": q_mid2, "q_high": q_high}


def compute_fixed_breaks(config: Dict[str, Optional[str]]) -> Dict[str, float]:
    scale_min = _parse_num(config.get("HEATMAP_SCALE_MIN"), -2)
    scale_mid = _parse_num(config.get("HEATMAP_SCALE_MID"), 0)
    scale_max = _parse_num(config.get("HEATMAP_SCALE_MAX"), 2)
    if not (scale_min < scale_mid < scale_max):
        scale_min, scale_mid, scale_max = -2, 0, 2
    scale_mid2 = _parse_num(
        config.get("HEATMAP_SCALE_MID2"),
        (scale_mid + scale_max) / 2,
    )
    if not (scale_mid < scale_mid2 < scale_max):
        scale_mid2 = (scale_mid + scale_max) / 2
    return {
        "scale_min": float(scale_min),
        "scale_mid": float(scale_mid),
        "scale_mid2": float(scale_mid2),
        "scale_max": float(scale_max),
    }


def compute_quantile_breaks(
    values: Iterable[float], config: Dict[str, Optional[str]]
) -> Dict[str, float]:
    finite_values = sorted(_collect_finite_values(values))
    if not finite_values:
        raise ValueError("No finite values for quantile heatmap scale")
    q = _parse_quantile_params(config)
    return {
        "scale_min": float(_quantile_from_sorted(finite_values, q["q_low"])),
        "scale_mid": float(_quantile_from_sorted(finite_values, 0.5)),
        "scale_mid2": float(_quantile_from_sorted(finite_values, q["q_mid2"])),
        "scale_max": float(_quantile_from_sorted(finite_values, q["q_high"])),
    }


def compute_heatmap_breaks(
    values: Iterable[float], config: Dict[str, Optional[str]]
) -> Dict[str, Any]:
    fixed_breaks = compute_fixed_breaks(config)
    requested_mode_raw = config.get("HEATMAP_SCALE_MODE")
    requested_mode = _parse_scale_mode(requested_mode_raw)
    quantile_params = _parse_quantile_params(config)
    finite_values = sorted(_collect_finite_values(values))
    n_finite = len(finite_values)

    quantile_breaks: Optional[Dict[str, float]] = None
    chosen_breaks = dict(fixed_breaks)
    chosen_mode = _HEATMAP_SCALE_FIXED
    fallback_reason: Optional[str] = None
    quantile_widened = False
    quantile_unavailable_reason: Optional[str] = None

    if n_finite < 10:
        quantile_unavailable_reason = "n_finite_lt_10"
    else:
        try:
            quantile_breaks = compute_quantile_breaks(finite_values, config)
        except Exception as exc:
            quantile_unavailable_reason = f"quantile_compute_error:{exc}"
        if quantile_breaks is not None:
            if quantile_breaks["scale_min"] == quantile_breaks["scale_max"]:
                quantile_unavailable_reason = "quantile_min_eq_max"
                quantile_breaks = None
            else:
                widened = _ensure_increasing_breaks_with_epsilon(quantile_breaks)
                if widened is None:
                    quantile_unavailable_reason = "quantile_non_increasing_breaks"
                    quantile_breaks = None
                else:
                    quantile_widened = widened != quantile_breaks
                    quantile_breaks = widened

    if requested_mode == _HEATMAP_SCALE_QUANTILE:
        if quantile_breaks is None:
            fallback_reason = quantile_unavailable_reason or "quantile_unavailable"
        else:
            chosen_breaks = dict(quantile_breaks)
            chosen_mode = _HEATMAP_SCALE_QUANTILE

    result = {
        "mode": chosen_mode,
        "requested_mode": requested_mode,
        "requested_mode_raw": requested_mode_raw,
        "chosen_breaks": chosen_breaks,
        "fixed_breaks": fixed_breaks,
        "quantile_breaks": quantile_breaks,
        "stats": {
            "n_finite": n_finite,
            "q_low": quantile_params["q_low"],
            "q_mid2": quantile_params["q_mid2"],
            "q_high": quantile_params["q_high"],
            "quantile_widened": quantile_widened,
            "fallback_reason": fallback_reason,
        },
    }

    if requested_mode == _HEATMAP_SCALE_QUANTILE and fallback_reason is not None:
        _LOG.warning(
            "heatmap action=scale_breaks_fallback requested_mode=quantile reason=%s n=%d q_low=%.4f q_mid2=%.4f q_high=%.4f using=fixed breaks=%s",
            fallback_reason,
            n_finite,
            quantile_params["q_low"],
            quantile_params["q_mid2"],
            quantile_params["q_high"],
            _format_breaks_for_log(fixed_breaks),
        )
    else:
        _LOG.info(
            "heatmap action=scale_breaks mode=%s n=%d q_low=%.4f q_mid2=%.4f q_high=%.4f breaks=%s",
            chosen_mode,
            n_finite,
            quantile_params["q_low"],
            quantile_params["q_mid2"],
            quantile_params["q_high"],
            _format_breaks_for_log(chosen_breaks),
        )

    return result
def _is_finite(value: Optional[float]) -> bool:
    return value is not None and math.isfinite(value)
def _format_breaks_line(breaks: Dict[str, float]) -> str:
    return (
        f"min={breaks['scale_min']:.6g}, "
        f"mid={breaks['scale_mid']:.6g}, "
        f"mid2={breaks['scale_mid2']:.6g}, "
        f"max={breaks['scale_max']:.6g}"
    )
def _scale_mode_text(mode: str) -> str:
    if mode == _HEATMAP_SCALE_QUANTILE:
        return "Scale mode: QUANTILE (within-run; colors are relative)"
    return "Scale mode: FIXED (absolute)"
