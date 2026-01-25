from __future__ import annotations

import csv
import html
import math
import re
from bisect import bisect_right
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, cast

_TRUTHY = {"1", "true", "yes", "on"}
_FALSY = {"0", "false", "no", "off"}


def _read_config_value(path: Path, key: str) -> Optional[str]:
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                line = re.sub(r"\s+#.*$", "", line).strip()
                if not line or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                if k.strip().upper() == key.upper():
                    return v.strip()
    except Exception:
        return None
    return None


def _parse_bool(value: Optional[str], default_value: bool) -> bool:
    if value is None:
        return default_value
    v = str(value).strip().lower()
    if v in _TRUTHY:
        return True
    if v in _FALSY:
        return False
    return default_value


def _parse_num(value: Optional[str], default_value: float) -> float:
    if value is None:
        return default_value
    try:
        parsed = float(str(value).strip())
    except Exception:
        return default_value
    if math.isnan(parsed):
        return default_value
    return parsed


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


def _is_finite(value: Optional[float]) -> bool:
    return value is not None and math.isfinite(value)


def _parse_t_selected(value: Any) -> Optional[float]:
    try:
        parsed = float(value)
    except Exception:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def _normalize_text(value: Any) -> str:
    return str(value or "").strip()


def _load_heatmap_rows(
    input_csv: Path, include_decoys: bool
) -> List[Dict[str, Any]]:
    with input_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "t_selected" not in fieldnames:
            raise ValueError("Missing required column: t_selected")
        has_target_id = "target_id" in fieldnames
        if not has_target_id:
            required_cols = {"pdb_id", "variant", "ph_label"}
            missing = required_cols - set(fieldnames)
            if missing:
                raise ValueError(
                    f"Missing required columns: {', '.join(sorted(missing))}"
                )
        has_ligand_display = "ligand_display" in fieldnames
        has_ligand_base = "ligand_base" in fieldnames
        if not has_ligand_display and not has_ligand_base:
            raise ValueError("Missing ligand_display/ligand_base columns in input CSV")

        rows: List[Dict[str, Any]] = []
        for row in reader:
            t_val = _parse_t_selected(row.get("t_selected"))
            if t_val is None:
                continue
            if not include_decoys and "is_decoy" in fieldnames:
                decoy_flag = (
                    str(row.get("is_decoy") or "").strip().lower() in _TRUTHY
                )
                if decoy_flag:
                    continue

            ligand_display = _normalize_text(row.get("ligand_display")) if has_ligand_display else ""
            ligand_base = _normalize_text(row.get("ligand_base")) if has_ligand_base else ""
            if not ligand_display or ligand_display.lower() == "nan":
                ligand_display = ligand_base
            ligand_name = ligand_display

            if has_target_id:
                target_id = _normalize_text(row.get("target_id"))
            else:
                target_id = "|".join(
                    [
                        _normalize_text(row.get("pdb_id")),
                        _normalize_text(row.get("variant")),
                        _normalize_text(row.get("ph_label")),
                    ]
                )

            row["t_selected"] = t_val
            row["ligand_name"] = ligand_name
            row["target_id"] = target_id
            rows.append(row)

    if not rows:
        raise ValueError("No usable t_selected values found in input CSV")
    return rows


def _aggregate_rows(
    rows: Iterable[Dict[str, Any]]
) -> Tuple[Dict[Tuple[str, str], Dict[str, Any]], Dict[str, float]]:
    agg: Dict[Tuple[str, str], Dict[str, Any]] = {}
    ligand_max: Dict[str, float] = {}
    for row in rows:
        ligand_name = row.get("ligand_name") or ""
        target_id = row.get("target_id") or ""
        t_val = row.get("t_selected")
        if not _is_finite(t_val):
            continue
        t_val = cast(float, t_val)
        key = (ligand_name, target_id)
        existing = agg.get(key)
        if existing is None or t_val > existing["t_selected"]:
            agg[key] = {
                "t_selected": t_val,
                "ligand_name": ligand_name,
                "target_id": target_id,
                "rank": _normalize_text(row.get("rank")),
                "pct_rank": _normalize_text(row.get("pct_rank")),
                "pose_valid_any": _normalize_text(row.get("pose_valid_any")),
                "pose_invalid_reason_top": _normalize_text(row.get("pose_invalid_reason_top")),
                "library": _normalize_text(row.get("library")),
            }
        current_max = ligand_max.get(ligand_name)
        if current_max is None or t_val > current_max:
            ligand_max[ligand_name] = t_val
    return agg, ligand_max


def _build_matrix(
    agg: Dict[Tuple[str, str], Dict[str, Any]],
    ligand_order: List[str],
    target_order: List[str],
) -> List[List[Optional[float]]]:
    matrix: List[List[Optional[float]]] = []
    for ligand in ligand_order:
        row = []
        for target in target_order:
            entry = agg.get((ligand, target))
            row.append(entry["t_selected"] if entry else None)
        matrix.append(row)
    return matrix


def _filter_empty(
    matrix: List[List[Optional[float]]],
    row_labels: List[str],
    col_labels: List[str],
) -> Tuple[List[List[Optional[float]]], List[str], List[str]]:
    row_keep = [any(_is_finite(v) for v in row) for row in matrix]
    col_keep = []
    for col_idx in range(len(col_labels)):
        col_keep.append(
            any(_is_finite(matrix[row_idx][col_idx]) for row_idx in range(len(row_labels)))
        )

    filtered_rows = [row for keep, row in zip(row_keep, matrix) if keep]
    filtered_row_labels = [label for keep, label in zip(row_keep, row_labels) if keep]
    if not filtered_row_labels:
        return [], [], []

    filtered_matrix: List[List[Optional[float]]] = []
    for row in filtered_rows:
        filtered_matrix.append([v for keep, v in zip(col_keep, row) if keep])
    filtered_col_labels = [label for keep, label in zip(col_keep, col_labels) if keep]
    return filtered_matrix, filtered_row_labels, filtered_col_labels


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
    min_val = np.nanmin(mat)
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


def _format_tooltip(lines: List[str]) -> str:
    escaped = html.escape("\n".join(lines), quote=True)
    return escaped.replace("\n", "&#10;")


def _format_value(value: Optional[float]) -> str:
    if not _is_finite(value):
        return ""
    return f"{value:.6g}"


def render_interactive_heatmap_html(
    repo_root: Path,
    run_id: str,
    input_csv: Path,
    top_k: int = 100,
    include_decoys: bool = False,
) -> str:
    config_path = repo_root / "config.txt"
    use_dendrogram = _parse_bool(
        _read_config_value(config_path, "USE_DENDROGRAM"), True
    )
    scale_min = _parse_num(
        _read_config_value(config_path, "HEATMAP_SCALE_MIN"), -3
    )
    scale_mid = _parse_num(
        _read_config_value(config_path, "HEATMAP_SCALE_MID"), 0
    )
    scale_max = _parse_num(
        _read_config_value(config_path, "HEATMAP_SCALE_MAX"), 3
    )
    if not (scale_min < scale_mid < scale_max):
        scale_min, scale_mid, scale_max = -3, 0, 3
    scale_mid2 = _parse_num(
        _read_config_value(config_path, "HEATMAP_SCALE_MID2"),
        (scale_mid + scale_max) / 2,
    )
    if not (scale_mid < scale_mid2 < scale_max):
        scale_mid2 = (scale_mid + scale_max) / 2

    color_low = _resolve_color(
        _read_config_value(config_path, "HEATMAP_COLOR_MIN"), "green"
    )
    color_mid = _resolve_color(
        _read_config_value(config_path, "HEATMAP_COLOR_MID"), "black"
    )
    color_mid2 = _resolve_color(
        _read_config_value(config_path, "HEATMAP_COLOR_MID2"), "orange"
    )
    color_high = _resolve_color(
        _read_config_value(config_path, "HEATMAP_COLOR_MAX"), "red"
    )

    resolved_input = input_csv
    if not resolved_input.exists():
        heatmap_path = repo_root / "data" / run_id / "heatmap_input.csv"
        master_path = repo_root / "data" / run_id / "master_rows.csv"
        if heatmap_path.exists():
            resolved_input = heatmap_path
        else:
            resolved_input = master_path
    if not resolved_input.exists():
        raise FileNotFoundError(f"Input CSV not found: {resolved_input}")

    rows = _load_heatmap_rows(resolved_input, include_decoys)
    agg, ligand_max = _aggregate_rows(rows)
    if not ligand_max:
        raise ValueError("No ligands available after filtering")

    top_k = int(top_k) if top_k and top_k > 0 else 100
    ranked_ligands = sorted(ligand_max.items(), key=lambda item: (-item[1], item[0]))
    top_ligands = sorted([name for name, _ in ranked_ligands[: min(top_k, len(ranked_ligands))]])

    agg = {
        key: value
        for key, value in agg.items()
        if key[0] in top_ligands
    }
    target_levels = sorted({target for _, target in agg.keys()})
    if not target_levels:
        raise ValueError("Heatmap matrix is empty after filtering")

    matrix = _build_matrix(agg, top_ligands, target_levels)
    matrix, row_labels, col_labels = _filter_empty(
        matrix, top_ligands, target_levels
    )
    if not row_labels or not col_labels:
        raise ValueError("Heatmap matrix is empty after filtering")

    row_order, col_order = _cluster_order(matrix, use_dendrogram)
    if row_order is not None:
        row_labels = [row_labels[idx] for idx in row_order]
    if col_order is not None:
        col_labels = [col_labels[idx] for idx in col_order]

    heatmap_breaks, heatmap_colors = _build_palette(
        scale_min,
        scale_mid,
        scale_mid2,
        scale_max,
        color_low,
        color_mid,
        color_mid2,
        color_high,
    )

    header_cells = "".join(
        f"<th class=\"hm-col-label\"><div>{html.escape(str(label))}</div></th>"
        for label in col_labels
    )
    body_rows: List[str] = []
    for ligand in row_labels:
        row_cells: List[str] = []
        for target in col_labels:
            cell = agg.get((ligand, target))
            t_val = cell["t_selected"] if cell else None
            color = _value_to_color(t_val, heatmap_breaks, heatmap_colors, scale_min, scale_max)
            tooltip_lines = [
                f"ligand: {ligand}",
                f"target: {target}",
                f"t_selected: {_format_value(t_val)}",
            ]
            data_attrs = {
                "ligand": ligand,
                "target": target,
                "t_selected": _format_value(t_val),
            }
            if cell:
                for key in (
                    "rank",
                    "pct_rank",
                    "pose_valid_any",
                    "pose_invalid_reason_top",
                    "library",
                ):
                    value = cell.get(key) or ""
                    if value != "":
                        tooltip_lines.append(f"{key}: {value}")
                        data_attrs[key] = value
            tooltip = _format_tooltip(tooltip_lines)
            attrs = " ".join(
                f'data-{k.replace("_", "-")}="{html.escape(str(v), quote=True)}"'
                for k, v in data_attrs.items()
            )
            row_cells.append(
                "<td>"
                f"<button type=\"button\" class=\"hm-cell\" style=\"background-color: {color};\" "
                f"title=\"{tooltip}\" aria-label=\"{tooltip}\" {attrs}></button>"
                "</td>"
            )
        body_rows.append(
            "<tr>"
            f"<th class=\"hm-row-label\">{html.escape(str(ligand))}</th>"
            f"{''.join(row_cells)}"
            "</tr>"
        )

    safe_run_id = re.sub(r"[^A-Za-z0-9_-]+", "-", str(run_id)).strip("-") or "run"
    container_id = f"heatmap-interactive-{safe_run_id}"
    detail_id = f"heatmap-detail-{safe_run_id}"
    style_block = """
<style>
  .heatmap-wrap { display: flex; gap: 16px; align-items: flex-start; flex-wrap: wrap; }
  .heatmap-scroll { max-height: 720px; overflow: auto; border: 1px solid #ddd; }
  .heatmap-table { border-collapse: collapse; }
  .heatmap-table th, .heatmap-table td { border: 1px solid #eee; padding: 0; }
  .heatmap-table th { background: #f7f7f7; font-weight: 600; }
  .heatmap-table .hm-row-label { position: sticky; left: 0; background: #fff; z-index: 1; padding: 4px 6px; }
  .heatmap-table thead th { position: sticky; top: 0; z-index: 2; }
  .heatmap-table .hm-col-label { height: 140px; vertical-align: bottom; padding: 0 6px; }
  .heatmap-table .hm-col-label > div { transform: rotate(-60deg); transform-origin: left bottom; white-space: nowrap; }
  .hm-cell { width: 18px; height: 18px; border: none; cursor: pointer; }
  .hm-cell:hover { outline: 1px solid #555; }
  .heatmap-detail { min-width: 220px; border: 1px solid #ddd; padding: 8px 10px; background: #fafafa; }
  .heatmap-detail-title { font-weight: 600; margin-bottom: 6px; }
  .heatmap-detail pre { margin: 0; white-space: pre-wrap; font-family: Consolas, monospace; font-size: 0.9em; }
</style>
""".strip()
    table_block = f"""
<div class="heatmap-wrap" id="{container_id}">
  <div class="heatmap-scroll">
    <table class="heatmap-table">
      <thead>
        <tr>
          <th class="hm-row-label">ligand</th>
          {header_cells}
        </tr>
      </thead>
      <tbody>
        {''.join(body_rows)}
      </tbody>
    </table>
  </div>
  <div class="heatmap-detail">
    <div class="heatmap-detail-title">Cell details</div>
    <pre id="{detail_id}">Click a cell to see details.</pre>
  </div>
</div>
""".strip()
    script_block = (
        "<script>\n"
        "  (function() {\n"
        f"    var container = document.getElementById(\"{container_id}\");\n"
        f"    var detail = document.getElementById(\"{detail_id}\");\n"
        "    if (!container || !detail) return;\n"
        "    var cells = container.querySelectorAll(\".hm-cell\");\n"
        "    cells.forEach(function(cell) {\n"
        "      cell.addEventListener(\"click\", function() {\n"
        "        var lines = [\n"
        "          \"ligand: \" + (cell.dataset.ligand || \"\"),\n"
        "          \"target: \" + (cell.dataset.target || \"\"),\n"
        "          \"t_selected: \" + (cell.dataset.tSelected || \"\")\n"
        "        ];\n"
        "        if (cell.dataset.rank) lines.push(\"rank: \" + cell.dataset.rank);\n"
        "        if (cell.dataset.pctRank) lines.push(\"pct_rank: \" + cell.dataset.pctRank);\n"
        "        if (cell.dataset.poseValidAny) lines.push(\"pose_valid_any: \" + cell.dataset.poseValidAny);\n"
        "        if (cell.dataset.poseInvalidReasonTop) {\n"
        "          lines.push(\"pose_invalid_reason_top: \" + cell.dataset.poseInvalidReasonTop);\n"
        "        }\n"
        "        if (cell.dataset.library) lines.push(\"library: \" + cell.dataset.library);\n"
        "        detail.textContent = lines.join(\"\\n\");\n"
        "      });\n"
        "    });\n"
        "  })();\n"
        "</script>"
    )
    html_out = "\n".join([style_block, table_block, script_block])
    return html_out
