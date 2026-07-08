from __future__ import annotations

import html
import logging
from bisect import bisect_right
from typing import Dict, List, Optional, Tuple, cast

from analysis.reporting.heatmap_html_scale import _is_finite

_LOG = logging.getLogger("heatmap-html")

def _cluster_order(
    matrix: List[List[Optional[float]]], use_dendrogram: bool
) -> Tuple[Optional[List[int]], Optional[List[int]]]:
    if not use_dendrogram:
        return None, None
    if not matrix or not matrix[0]:
        return None, None
    try:
        import numpy as np
        from scipy.cluster.hierarchy import leaves_list, linkage  # type: ignore[import-untyped]
        from scipy.spatial.distance import pdist  # type: ignore[import-untyped]
    except Exception:
        return None, None

    mat = np.array(
        [[float("nan") if v is None else float(v) for v in row] for row in matrix],
        dtype=float,
    )
    if not np.isfinite(mat).any():
        return None, None
    min_val: float = float(np.nanmin(mat))
    mat_for_cluster = np.where(np.isfinite(mat), mat, min_val)

    row_order = None
    col_order = None
    if mat_for_cluster.shape[0] > 1:
        row_order = leaves_list(linkage(pdist(mat_for_cluster), method="complete")).tolist()
    if mat_for_cluster.shape[1] > 1:
        col_order = leaves_list(linkage(pdist(mat_for_cluster.T), method="complete")).tolist()
    return row_order, col_order


def _clamp_value(value: float, min_val: float, max_val: float) -> float:
    return min(max(value, min_val), max_val)


def _value_to_color(
    value: Optional[float],
    heatmap_breaks: List[float],
    heatmap_colors: List[str],
    scale_min: float,
    scale_max: float,
) -> str:
    if not _is_finite(value):
        return "white"
    assert value is not None
    clamped = _clamp_value(value, scale_min, scale_max)
    idx = bisect_right(heatmap_breaks, clamped) - 1
    if idx < 0:
        idx = 0
    if idx >= len(heatmap_colors):
        idx = len(heatmap_colors) - 1
    return heatmap_colors[idx]


def _value_for_breaks(
    value: Optional[float], breaks: Optional[Dict[str, float]]
) -> Optional[float]:
    if not _is_finite(value):
        return None
    if breaks is None:
        return value
    return _clamp_value(
        cast(float, value),
        float(breaks["scale_min"]),
        float(breaks["scale_max"]),
    )


def _format_tooltip(lines: List[str]) -> str:
    escaped = html.escape("\n".join(lines), quote=True)
    return escaped.replace("\n", "&#10;")


def _format_value(value: Optional[float]) -> str:
    if not _is_finite(value):
        return ""
    return f"{value:.6g}"


def _format_percent(value: Optional[float], decimals: int = 2) -> str:
    if not _is_finite(value):
        return ""
    return f"{cast(float, value) * 100.0:.{max(0, int(decimals))}f}%"
