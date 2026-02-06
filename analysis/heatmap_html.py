from __future__ import annotations

import base64
import csv
import html
import json
import logging
import math
import os
import re
import time
from bisect import bisect_right
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, cast

_TRUTHY = {"1", "true", "yes", "on"}
_FALSY = {"0", "false", "no", "off"}
_LOG = logging.getLogger("heatmap-html")
_CLUSTERGRAMMER_WARNED = False
_FONT_MIME_BY_EXT = {
    ".woff2": "font/woff2",
    ".woff": "font/woff",
    ".ttf": "font/ttf",
    ".eot": "application/vnd.ms-fontobject",
    ".svg": "image/svg+xml",
}


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


def _is_truthy_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        try:
            numeric = float(value)
        except Exception:
            numeric = 0.0
        if math.isnan(numeric):
            return False
        return numeric != 0.0
    return _normalize_text(value).lower() in _TRUTHY


def _build_target_id(pdb_id: Any, variant: Any, ph_label: Any) -> str:
    return "|".join(
        [
            _normalize_text(pdb_id),
            _normalize_text(variant),
            _normalize_text(ph_label),
        ]
    )


def _extract_pdb_id_from_target_id(target_id: Any) -> str:
    token = _normalize_text(target_id).split("|", 1)[0].strip().upper()
    return token


def _resolve_target_display_label(row: Dict[str, Any], target_id: str) -> str:
    pdb_id = _normalize_text(row.get("pdb_id")).upper()
    if not pdb_id:
        pdb_id = _extract_pdb_id_from_target_id(target_id)
    target_name = _normalize_text(row.get("target_name"))
    variant = _normalize_text(row.get("variant"))
    ph_label = _normalize_text(row.get("ph_label"))
    if (not variant or not ph_label) and "|" in target_id:
        parts = target_id.split("|")
        if len(parts) >= 3:
            if not variant:
                variant = _normalize_text(parts[1])
            if not ph_label:
                ph_label = _normalize_text(parts[2])
    variant_ph = " | ".join([part for part in (variant, ph_label) if part])
    if target_name and pdb_id:
        lines = [pdb_id, target_name]
        if variant_ph:
            lines.append(variant_ph)
        return "\n".join(lines)
    if pdb_id and variant_ph:
        return f"{pdb_id}\n{variant_ph}"
    if pdb_id:
        return pdb_id
    if target_name:
        return target_name
    return target_id


def _build_target_display_map(rows: List[Dict[str, Any]]) -> Dict[str, str]:
    labels: Dict[str, str] = {}
    has_name: Dict[str, bool] = {}
    for row in rows:
        target_id = _normalize_text(row.get("target_id"))
        if not target_id:
            continue
        candidate = _resolve_target_display_label(row, target_id)
        candidate_has_name = bool(_normalize_text(row.get("target_name")))
        if target_id not in labels:
            labels[target_id] = candidate
            has_name[target_id] = candidate_has_name
            continue
        # Prefer labels that include target_name metadata when available.
        if candidate_has_name and not has_name.get(target_id, False):
            labels[target_id] = candidate
            has_name[target_id] = True
    return labels


def _sample_names(values: List[str], max_items: int = 5, max_len: int = 40) -> List[str]:
    sample: List[str] = []
    for value in values[:max_items]:
        text = str(value)
        if len(text) > max_len:
            text = text[: max_len - 3] + "..."
        sample.append(text)
    return sample


def _load_heatmap_rows_csv(
    input_csv: Path, include_decoys: bool, allowed_pdb_ids: Optional[Set[str]] = None
) -> List[Dict[str, Any]]:
    normalized_allowed_pdb_ids: Optional[Set[str]] = None
    if allowed_pdb_ids is not None:
        normalized_allowed_pdb_ids = {
            str(pdb_id).strip().upper() for pdb_id in allowed_pdb_ids if str(pdb_id).strip()
        }

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
                decoy_flag = _is_truthy_value(row.get("is_decoy"))
                if decoy_flag:
                    continue

            ligand_display = _normalize_text(row.get("ligand_display")) if has_ligand_display else ""
            ligand_base = _normalize_text(row.get("ligand_base")) if has_ligand_base else ""
            if not ligand_display or ligand_display.lower() == "nan":
                ligand_display = ligand_base
            ligand_name = ligand_display
            if not ligand_name:
                continue

            if has_target_id:
                target_id = _normalize_text(row.get("target_id"))
                if not target_id and {"pdb_id", "variant", "ph_label"}.issubset(fieldnames):
                    target_id = _build_target_id(
                        row.get("pdb_id"),
                        row.get("variant"),
                        row.get("ph_label"),
                    )
            else:
                target_id = _build_target_id(
                    row.get("pdb_id"),
                    row.get("variant"),
                    row.get("ph_label"),
                )
            if not target_id:
                continue

            pdb_id = _normalize_text(row.get("pdb_id")).upper() if "pdb_id" in fieldnames else ""
            if not pdb_id:
                pdb_id = _extract_pdb_id_from_target_id(target_id)
            if (
                normalized_allowed_pdb_ids is not None
                and pdb_id not in normalized_allowed_pdb_ids
            ):
                continue

            row["t_selected"] = t_val
            row["ligand_name"] = ligand_name
            row["target_id"] = target_id
            row["pdb_id"] = pdb_id
            rows.append(row)

    if not rows:
        raise ValueError("No usable t_selected values found in input CSV")
    return rows


def _load_heatmap_rows_parquet(
    input_path: Path, include_decoys: bool, allowed_pdb_ids: Optional[Set[str]] = None
) -> List[Dict[str, Any]]:
    normalized_allowed_pdb_ids: Optional[Set[str]] = None
    if allowed_pdb_ids is not None:
        normalized_allowed_pdb_ids = {
            str(pdb_id).strip().upper() for pdb_id in allowed_pdb_ids if str(pdb_id).strip()
        }

    try:
        import pyarrow.dataset as ds  # type: ignore[import-untyped]
    except Exception as exc:  # pragma: no cover - import error path
        raise RuntimeError(
            "pyarrow is required to read parquet heatmap inputs"
        ) from exc

    dataset = ds.dataset(str(input_path), format="parquet", partitioning="hive")
    schema_names = set(dataset.schema.names)
    if "t_selected" not in schema_names:
        raise ValueError("Missing required column: t_selected")

    has_target_id = "target_id" in schema_names
    has_pdb_id = "pdb_id" in schema_names
    has_target_components = {"pdb_id", "variant", "ph_label"}.issubset(schema_names)
    if not has_target_id and not has_target_components:
        raise ValueError(
            "Missing required columns: target_id or pdb_id/variant/ph_label"
        )

    has_ligand_display = "ligand_display" in schema_names
    has_ligand_base = "ligand_base" in schema_names
    if not has_ligand_display and not has_ligand_base:
        raise ValueError("Missing ligand_display/ligand_base columns in input parquet")

    scan_columns: List[str] = ["t_selected"]
    for col in ("target_id", "pdb_id", "variant", "ph_label"):
        if col in schema_names and col not in scan_columns:
            scan_columns.append(col)
    if "target_name" in schema_names and "target_name" not in scan_columns:
        scan_columns.append("target_name")
    for col in (
        "ligand_display",
        "ligand_base",
        "is_decoy",
        "rank",
        "pct_rank",
        "pose_valid_any",
        "pose_invalid_reason_top",
        "library",
    ):
        if col in schema_names and col not in scan_columns:
            scan_columns.append(col)

    rows: List[Dict[str, Any]] = []
    scanner = dataset.scanner(columns=scan_columns, use_threads=True)
    for batch in scanner.to_batches():
        payload = batch.to_pydict()
        n_rows = batch.num_rows
        for idx in range(n_rows):
            t_val = _parse_t_selected(payload["t_selected"][idx])
            if t_val is None:
                continue

            if not include_decoys and "is_decoy" in payload:
                if _is_truthy_value(payload["is_decoy"][idx]):
                    continue

            ligand_display = (
                _normalize_text(payload["ligand_display"][idx])
                if has_ligand_display
                else ""
            )
            ligand_base = (
                _normalize_text(payload["ligand_base"][idx]) if has_ligand_base else ""
            )
            if not ligand_display or ligand_display.lower() == "nan":
                ligand_display = ligand_base
            ligand_name = ligand_display
            if not ligand_name:
                continue

            if has_target_id:
                target_id = _normalize_text(payload["target_id"][idx])
                if not target_id and has_target_components:
                    target_id = _build_target_id(
                        payload["pdb_id"][idx],
                        payload["variant"][idx],
                        payload["ph_label"][idx],
                    )
            else:
                target_id = _build_target_id(
                    payload["pdb_id"][idx],
                    payload["variant"][idx],
                    payload["ph_label"][idx],
                )
            if not target_id:
                continue

            pdb_id = _normalize_text(payload["pdb_id"][idx]).upper() if has_pdb_id else ""
            if not pdb_id:
                pdb_id = _extract_pdb_id_from_target_id(target_id)
            if (
                normalized_allowed_pdb_ids is not None
                and pdb_id not in normalized_allowed_pdb_ids
            ):
                continue

            row: Dict[str, Any] = {
                "t_selected": t_val,
                "ligand_name": ligand_name,
                "target_id": target_id,
                "pdb_id": pdb_id,
            }
            for col in (
                "rank",
                "pct_rank",
                "pose_valid_any",
                "pose_invalid_reason_top",
                "library",
                "target_name",
            ):
                if col in payload:
                    row[col] = payload[col][idx]
            rows.append(row)

    if not rows:
        raise ValueError("No usable t_selected values found in input parquet")
    return rows


def _load_heatmap_rows(
    input_path: Path,
    include_decoys: bool,
    allowed_pdb_ids: Optional[Set[str]] = None,
) -> List[Dict[str, Any]]:
    is_parquet_source = input_path.is_dir() or input_path.suffix.lower() == ".parquet"
    if is_parquet_source:
        return _load_heatmap_rows_parquet(
            input_path, include_decoys, allowed_pdb_ids=allowed_pdb_ids
        )
    return _load_heatmap_rows_csv(
        input_path, include_decoys, allowed_pdb_ids=allowed_pdb_ids
    )


def _aggregate_rows(
    rows: Iterable[Dict[str, Any]]
) -> Tuple[Dict[Tuple[str, str], Dict[str, Any]], Dict[str, float]]:
    agg: Dict[Tuple[str, str], Dict[str, Any]] = {}
    ligand_max: Dict[str, float] = {}
    for row in rows:
        ligand_name = row.get("ligand_name") or ""
        target_id = row.get("target_id") or ""
        target_label = row.get("target_label") or target_id
        t_val = row.get("t_selected")
        if not _is_finite(t_val):
            continue
        t_val = cast(float, t_val)
        key = (ligand_name, target_label)
        existing = agg.get(key)
        if existing is None or t_val > existing["t_selected"]:
            agg[key] = {
                "t_selected": t_val,
                "ligand_name": ligand_name,
                "target_id": target_label,
                "target_raw": _normalize_text(target_id),
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


def _warn_clustergrammer_failure(exc: Exception) -> None:
    global _CLUSTERGRAMMER_WARNED
    if _CLUSTERGRAMMER_WARNED:
        return
    _CLUSTERGRAMMER_WARNED = True
    _LOG.warning("heatmap action=clustergrammer_fallback error=%s", exc)


def _ensure_pandas_ix_compat(pd_module: Any) -> None:
    if hasattr(pd_module.DataFrame, "ix"):
        return
    pd_module.DataFrame.ix = property(lambda self: self.loc)  # type: ignore[attr-defined]


def _load_clustergrammer() -> Optional[Tuple[Any, Any]]:
    try:
        import pandas as pd  # type: ignore[import-untyped]
        from clustergrammer import Network  # type: ignore[import-untyped]
    except Exception:
        return None
    _ensure_pandas_ix_compat(pd)
    return Network, pd


def _clustergrammer_key(name: str) -> str:
    return str(name).replace("_", " ")


def _escape_inline_text(text: str) -> str:
    return re.sub(r"</(script|style)", r"<\\/\1", text, flags=re.IGNORECASE)


def _read_asset_text(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        return handle.read()


def _inline_font_urls(css_text: str, fonts_dir: Path) -> str:
    def replacer(match: re.Match[str]) -> str:
        raw_url = match.group(1).strip().strip("'\"")
        if not raw_url or raw_url.startswith("data:") or raw_url.startswith("http"):
            return match.group(0)
        cleaned = raw_url.split("#", 1)[0].split("?", 1)[0]
        filename = Path(cleaned).name
        if not filename:
            return match.group(0)
        font_path = fonts_dir / filename
        if not font_path.exists():
            return match.group(0)
        mime = _FONT_MIME_BY_EXT.get(font_path.suffix.lower())
        if not mime:
            return match.group(0)
        encoded = base64.b64encode(font_path.read_bytes()).decode("ascii")
        return f"url('data:{mime};base64,{encoded}')"

    return re.sub(r"url\(([^)]+)\)", replacer, css_text)


def _clustergrammer_inline_assets(repo_root: Path) -> str:
    assets_root = repo_root / "report_assets" / "clustergrammer"
    fonts_dir = assets_root / "lib" / "fonts"
    css_paths = [
        assets_root / "lib" / "css" / "bootstrap.css",
        assets_root / "lib" / "css" / "font-awesome.min.css",
        assets_root / "css" / "custom.css",
    ]
    js_paths = [
        assets_root / "lib" / "js" / "d3.js",
        assets_root / "lib" / "js" / "jquery-1.11.2.min.js",
        assets_root / "lib" / "js" / "underscore-min.js",
        assets_root / "lib" / "js" / "bootstrap.min.js",
        assets_root / "clustergrammer.js",
    ]
    css_parts = []
    for path in css_paths:
        text = _read_asset_text(path)
        text = _inline_font_urls(text, fonts_dir)
        css_parts.append(text)
    css_blob = _escape_inline_text("\n".join(css_parts))
    blocks = [f"<style>\n{css_blob}\n</style>"]
    for path in js_paths:
        js_text = _escape_inline_text(_read_asset_text(path))
        blocks.append(f"<script>\n{js_text}\n</script>")
    return "\n".join(blocks)


def _clustergrammer_asset_urls(
    repo_root: Path, run_id: str, offline_assets: bool
) -> Dict[str, str]:
    if offline_assets:
        report_dir = repo_root / "data" / run_id
        assets_root = repo_root / "report_assets" / "clustergrammer"
        rel_root = Path(os.path.relpath(assets_root, report_dir)).as_posix()
        return {
            "d3": f"{rel_root}/lib/js/d3.js",
            "jquery": f"{rel_root}/lib/js/jquery-1.11.2.min.js",
            "underscore": f"{rel_root}/lib/js/underscore-min.js",
            "bootstrap_js": f"{rel_root}/lib/js/bootstrap.min.js",
            "bootstrap_css": f"{rel_root}/lib/css/bootstrap.css",
            "font_awesome_css": f"{rel_root}/lib/css/font-awesome.min.css",
            "custom_css": f"{rel_root}/css/custom.css",
            "clustergrammer_js": f"{rel_root}/clustergrammer.js",
        }
    return {
        "d3": "https://d3js.org/d3.v3.min.js",
        "jquery": "https://code.jquery.com/jquery-1.11.2.min.js",
        "underscore": "https://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.8.3/underscore-min.js",
        "bootstrap_js": "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
        "bootstrap_css": "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css",
        "font_awesome_css": "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
        "custom_css": "https://unpkg.com/clustergrammer@1.19.5/css/custom.css",
        "clustergrammer_js": "https://unpkg.com/clustergrammer@1.19.5/clustergrammer.js",
    }


def _clustergrammer_asset_block(
    repo_root: Path,
    run_id: str,
    offline_assets: bool,
    inline_assets: bool,
) -> str:
    if inline_assets:
        try:
            return _clustergrammer_inline_assets(repo_root)
        except Exception as exc:
            _LOG.warning("heatmap action=inline_assets_failed error=%s", exc)
    assets = _clustergrammer_asset_urls(repo_root, run_id, offline_assets)
    return "\n".join(
        [
            f"<link rel=\"stylesheet\" href=\"{html.escape(assets['bootstrap_css'], quote=True)}\">",
            f"<link rel=\"stylesheet\" href=\"{html.escape(assets['font_awesome_css'], quote=True)}\">",
            f"<link rel=\"stylesheet\" href=\"{html.escape(assets['custom_css'], quote=True)}\">",
            f"<script src=\"{html.escape(assets['d3'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['jquery'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['underscore'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['bootstrap_js'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['clustergrammer_js'], quote=True)}\"></script>",
        ]
    )


def _build_clustergrammer_viz_json(
    matrix: List[List[Optional[float]]],
    row_labels: List[str],
    col_labels: List[str],
    use_dendrogram: bool,
) -> Dict[str, Any]:
    loaded = _load_clustergrammer()
    if loaded is None:
        raise ImportError("clustergrammer not available")
    Network, pd = loaded
    df = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
    df_filled = df
    if df.isna().any().any():
        min_val = df.min().min()
        fill_value = float(min_val) if pd.notna(min_val) else 0.0
        df_filled = df.fillna(fill_value)
    net = Network()
    net.load_df(df_filled)
    net.cluster(run_clustering=use_dendrogram, dendro=use_dendrogram)
    viz_json = net.export_net_json("viz")
    if isinstance(viz_json, str):
        return json.loads(viz_json)
    return cast(Dict[str, Any], viz_json)


def _render_clustergrammer_html(
    repo_root: Path,
    run_id: str,
    row_labels: List[str],
    col_labels: List[str],
    matrix: List[List[Optional[float]]],
    agg: Dict[Tuple[str, str], Dict[str, Any]],
    use_dendrogram: bool,
    offline_assets: bool,
    inline_assets: bool,
    ligand_page_prefix: str,
    payload_start: Optional[float] = None,
) -> Optional[str]:
    build_start = payload_start if payload_start is not None else time.perf_counter()
    try:
        viz_json = _build_clustergrammer_viz_json(
            matrix, row_labels, col_labels, use_dendrogram
        )
    except Exception as exc:
        _warn_clustergrammer_failure(exc)
        return None

    row_set = set(row_labels)
    col_set = set(col_labels)
    cell_meta: Dict[str, Dict[str, str]] = {}
    for (ligand, target), cell in agg.items():
        if ligand not in row_set or target not in col_set:
            continue
        key = f"{_clustergrammer_key(ligand)}||{_clustergrammer_key(target)}"
        meta_entry: Dict[str, str] = {
            "t_selected": _format_value(cast(Optional[float], cell.get("t_selected"))),
            "ligand_raw": ligand,
            "target_raw": _normalize_text(cell.get("target_raw")) or target,
        }
        for meta_key in (
            "rank",
            "pct_rank",
            "pose_valid_any",
            "pose_invalid_reason_top",
            "library",
        ):
            value = cell.get(meta_key) or ""
            if value != "":
                meta_entry[meta_key] = str(value)
        cell_meta[key] = meta_entry

    safe_run_id = re.sub(r"[^A-Za-z0-9_-]+", "-", str(run_id)).strip("-") or "run"
    container_id = f"cg-heatmap-{safe_run_id}"
    detail_id = f"cg-heatmap-detail-{safe_run_id}"
    viz_json_id = f"cg-heatmap-data-{safe_run_id}"
    meta_json_id = f"cg-heatmap-meta-{safe_run_id}"

    assets_block = _clustergrammer_asset_block(
        repo_root=repo_root,
        run_id=run_id,
        offline_assets=offline_assets,
        inline_assets=inline_assets,
    )

    style_block = """
<style>
  .cg-heatmap-wrap { display: flex; gap: 16px; align-items: flex-start; flex-wrap: wrap; }
  .cg-heatmap-viz { min-width: 640px; min-height: 720px; border: 1px solid #ddd; }
  .cg-heatmap-detail { min-width: 220px; border: 1px solid #ddd; padding: 8px 10px; background: #fafafa; }
  .cg-heatmap-detail-title { font-weight: 600; margin-bottom: 6px; }
  .cg-heatmap-detail pre { margin: 0; white-space: pre-wrap; font-family: Consolas, monospace; font-size: 0.9em; }
  .heatmap-debug { margin-top: 12px; }
  .heatmap-debug summary { cursor: pointer; font-weight: 600; }
  .heatmap-debug pre { margin: 6px 0 0; padding: 8px; border: 1px solid #ddd; background: #f5f5f5; max-height: 240px; overflow: auto; white-space: pre-wrap; font-family: Consolas, monospace; font-size: 0.85em; }
</style>
""".strip()

    viz_json_text = json.dumps(viz_json, ensure_ascii=True).replace("</", "<\\/")
    meta_json_text = json.dumps(cell_meta, ensure_ascii=True).replace("</", "<\\/")
    mat_rows = len(matrix)
    mat_cols = len(matrix[0]) if matrix else 0
    elapsed_s = time.perf_counter() - build_start
    viz_bytes = len(viz_json_text.encode("utf-8"))
    meta_bytes = len(meta_json_text.encode("utf-8"))
    _LOG.info(
        "[heatmap.payload.build.done] rows=%d cols=%d mat_shape=%dx%d bytes=%d meta_bytes=%d elapsed_s=%.3f",
        len(row_labels),
        len(col_labels),
        mat_rows,
        mat_cols,
        viz_bytes,
        meta_bytes,
        elapsed_s,
    )

    container_block = f"""
<div class="cg-heatmap-wrap" id="{container_id}-wrap">
  <div class="cg-heatmap-viz" id="{container_id}">
    <div class="wait_message">Loading heatmap...</div>
  </div>
  <div class="cg-heatmap-detail">
    <div class="cg-heatmap-detail-title">Cell details</div>
    <pre id="{detail_id}">Click a cell to see details.</pre>
  </div>
</div>
<details class="heatmap-debug" id="heatmap-debug">
  <summary>Heatmap debug</summary>
  <pre id="heatmap-debug-log"></pre>
</details>
<script type="application/json" id="{viz_json_id}">{viz_json_text}</script>
<script type="application/json" id="{meta_json_id}">{meta_json_text}</script>
""".strip()

    ligand_prefix_json = json.dumps(ligand_page_prefix or "", ensure_ascii=True)
    row_order_js = (
        "    args.row_order = \"alpha\";\n    args.col_order = \"alpha\";\n"
        if not use_dendrogram
        else ""
    )
    script_block = (
        "<script>\n"
        "  (function() {\n"
        f"    var container = document.getElementById(\"{container_id}\");\n"
        f"    var detail = document.getElementById(\"{detail_id}\");\n"
        f"    var vizEl = document.getElementById(\"{viz_json_id}\");\n"
        f"    var metaEl = document.getElementById(\"{meta_json_id}\");\n"
        "    var debugEl = document.getElementById(\"heatmap-debug-log\");\n"
        "    function debugLine(message) {\n"
        "      if (!debugEl) return;\n"
        "      debugEl.textContent += String(message) + \"\\n\";\n"
        "    }\n"
        "    function nowMs() {\n"
        "      if (window.performance && typeof window.performance.now === \"function\") {\n"
        "        return window.performance.now();\n"
        "      }\n"
        "      return Date.now();\n"
        "    }\n"
        "    function scanGlobals() {\n"
        "      var matches = [];\n"
        "      for (var key in window) {\n"
        "        if (!Object.prototype.hasOwnProperty.call(window, key)) continue;\n"
        "        if (/cluster|gram/i.test(key)) {\n"
        "          matches.push(key);\n"
        "          if (matches.length >= 20) break;\n"
        "        }\n"
        "      }\n"
        "      return matches;\n"
        "    }\n"
        "    function dumpError(label, err) {\n"
        "      var msg = String(label || \"[heatmap.error]\");\n"
        "      if (err) {\n"
        "        if (err.message) msg += \" message=\" + err.message;\n"
        "        if (err.stack) msg += \"\\n\" + err.stack;\n"
        "      }\n"
        "      debugLine(msg);\n"
        "    }\n"
        "    var prevOnError = window.onerror;\n"
        "    window.onerror = function(message, source, lineno, colno, error) {\n"
        "      var parts = [\"[heatmap.error]\", String(message || \"\"), String(source || \"\"), String(lineno || \"\"), String(colno || \"\")];\n"
        "      debugLine(parts.join(\" \"));\n"
        "      if (error && error.stack) debugLine(error.stack);\n"
        "      if (window.console && console.error) console.error(\"Heatmap error\", error || message);\n"
        "      if (typeof prevOnError === \"function\") {\n"
        "        return prevOnError(message, source, lineno, colno, error);\n"
        "      }\n"
        "      return false;\n"
        "    };\n"
        "    window.addEventListener(\"unhandledrejection\", function(evt) {\n"
        "      var reason = evt && evt.reason ? evt.reason : \"unknown\";\n"
        "      debugLine(\"[heatmap.unhandledrejection] \" + String(reason));\n"
        "      if (reason && reason.stack) debugLine(reason.stack);\n"
        "      if (window.console && console.error) console.error(\"Heatmap unhandled rejection\", reason);\n"
        "    });\n"
        "    if (!container || !vizEl || !metaEl) {\n"
        "      debugLine(\"[heatmap.debug] missing_elements container=\" + !!container + \" viz=\" + !!vizEl + \" meta=\" + !!metaEl);\n"
        "      return;\n"
        "    }\n"
        f"    var ligandPrefix = {ligand_prefix_json};\n"
        "    if (ligandPrefix && ligandPrefix.slice(-1) !== \"/\") ligandPrefix += \"/\";\n"
        "    function escapeHtml(value) {\n"
        "      var text = String(value == null ? \"\" : value);\n"
        "      return text.replace(/[&<>\\\"]/g, function(ch) {\n"
        "        return {\"&\": \"&amp;\", \"<\": \"&lt;\", \">\": \"&gt;\", \"\\\"\": \"&quot;\"}[ch] || ch;\n"
        "      });\n"
        "    }\n"
        "    function setWaitMessage(message) {\n"
        "      if (!container) return;\n"
        "      var text = message || \"Heatmap failed to render.\";\n"
        "      var wait = container.querySelector('.wait_message');\n"
        "      if (wait) {\n"
        "        wait.textContent = text;\n"
        "        return;\n"
        "      }\n"
        "      var div = document.createElement(\"div\");\n"
        "      div.className = \"wait_message\";\n"
        "      div.textContent = text;\n"
        "      container.appendChild(div);\n"
        "    }\n"
        "    function clearWaitMessage() {\n"
        "      if (!container) return;\n"
        "      var wait = container.querySelector('.wait_message');\n"
        "      if (!wait) return;\n"
        "      if (wait.remove) {\n"
        "        wait.remove();\n"
        "      } else if (wait.parentNode) {\n"
        "        wait.parentNode.removeChild(wait);\n"
        "      }\n"
        "    }\n"
        "    debugLine(\"[heatmap.debug] Clustergrammer_type=\" + (typeof Clustergrammer));\n"
        "    var readStart = nowMs();\n"
        "    var vizText = vizEl.textContent || \"\";\n"
        "    var metaText = metaEl.textContent || \"\";\n"
        "    var readEnd = nowMs();\n"
        "    debugLine(\"[heatmap.debug] json_text_ms=\" + (readEnd - readStart).toFixed(1) + \" viz_chars=\" + vizText.length + \" meta_chars=\" + metaText.length);\n"
        "    if (typeof Clustergrammer === \"undefined\") {\n"
        "      setWaitMessage(\"Clustergrammer not available.\");\n"
        "      var globals = scanGlobals();\n"
        "      debugLine(\"[heatmap.debug] clustergrammer_missing globals=\" + globals.join(\", \"));\n"
        "      if (window.console && console.error) console.error(\"Clustergrammer not available\", globals);\n"
        "      return;\n"
        "    }\n"
        "    var vizData = {};\n"
        "    var cellMeta = {};\n"
        "    var parseStart = nowMs();\n"
        "    try { vizData = JSON.parse(vizText || \"{}\"); }\n"
        "    catch (err) { setWaitMessage(\"Heatmap data parse failed.\"); dumpError(\"[heatmap.debug] viz_parse_failed\", err); if (window.console && console.error) console.error(\"Heatmap data parse failed\", err); return; }\n"
        "    var parseEnd = nowMs();\n"
        "    debugLine(\"[heatmap.debug] viz_parse_ms=\" + (parseEnd - parseStart).toFixed(1));\n"
        "    var metaParseStart = nowMs();\n"
        "    try { cellMeta = JSON.parse(metaText || \"{}\"); }\n"
        "    catch (err) { setWaitMessage(\"Heatmap metadata parse failed.\"); dumpError(\"[heatmap.debug] meta_parse_failed\", err); if (window.console && console.error) console.error(\"Heatmap metadata parse failed\", err); return; }\n"
        "    var metaParseEnd = nowMs();\n"
        "    debugLine(\"[heatmap.debug] meta_parse_ms=\" + (metaParseEnd - metaParseStart).toFixed(1));\n"
        "    var rowNodes = vizData.row_nodes || [];\n"
        "    var colNodes = vizData.col_nodes || [];\n"
        "    var mat = vizData.mat || [];\n"
        "    var matRows = Array.isArray(mat) ? mat.length : 0;\n"
        "    var matCols = matRows && Array.isArray(mat[0]) ? mat[0].length : 0;\n"
        "    debugLine(\"[heatmap.debug] row_nodes=\" + rowNodes.length + \" col_nodes=\" + colNodes.length + \" mat=\" + matRows + \"x\" + matCols);\n"
        "    function buildLines(rowName, colName, meta) {\n"
        "      var lines = [\n"
        "        \"ligand: \" + rowName,\n"
        "        \"target: \" + colName,\n"
        "        \"t_selected: \" + (meta && meta.t_selected ? meta.t_selected : \"NA\")\n"
        "      ];\n"
        "      if (meta && meta.rank) lines.push(\"rank: \" + meta.rank);\n"
        "      if (meta && meta.pct_rank) lines.push(\"pct_rank: \" + meta.pct_rank);\n"
        "      if (meta && meta.pose_invalid_reason_top) {\n"
        "        lines.push(\"pose_invalid_reason_top: \" + meta.pose_invalid_reason_top);\n"
        "      }\n"
        "      return lines;\n"
        "    }\n"
        "    function tileKey(rowName, colName) {\n"
        "      return String(rowName || \"\") + \"||\" + String(colName || \"\");\n"
        "    }\n"
        "    function lookupMeta(rowName, colName) {\n"
        "      return cellMeta[tileKey(rowName, colName)];\n"
        "    }\n"
        "    function updateDetails(rowName, colName, meta) {\n"
        "      if (!detail) return;\n"
        "      if (!meta) {\n"
        "        detail.textContent = \"No data for this cell.\";\n"
        "        return;\n"
        "      }\n"
        "      detail.textContent = buildLines(rowName, colName, meta).join(\"\\n\");\n"
        "    }\n"
        "    function tooltipHtml(rowName, colName, meta) {\n"
        "      return buildLines(rowName, colName, meta)\n"
        "        .map(function(line) { return \"<div>\" + escapeHtml(line) + \"</div>\"; })\n"
        "        .join(\"\");\n"
        "    }\n"
        "    function handleTileClick(evt) {\n"
        "      var data = this.__data__ || {};\n"
        "      var rowName = data.row_name || data.row || \"\";\n"
        "      var colName = data.col_name || data.col || \"\";\n"
        "      var meta = lookupMeta(rowName, colName);\n"
        "      updateDetails(rowName, colName, meta);\n"
        "      if (!meta || !ligandPrefix) return;\n"
        "      var ligand = meta.ligand_raw || rowName;\n"
        "      var target = meta.target_raw || colName;\n"
        "      var url = ligandPrefix + encodeURIComponent(ligand) + \".html?target=\" + encodeURIComponent(target);\n"
        "      if (evt && (evt.metaKey || evt.ctrlKey)) {\n"
        "        window.open(url, \"_blank\");\n"
        "      } else {\n"
        "        window.location.href = url;\n"
        "      }\n"
        "    }\n"
        "    function attachTileHandlers() {\n"
        "      var tiles = container.querySelectorAll(\"rect.tile\");\n"
        "      for (var i = 0; i < tiles.length; i++) {\n"
        "        var tile = tiles[i];\n"
        "        var data = tile.__data__ || {};\n"
        "        var rowName = data.row_name || data.row || \"\";\n"
        "        var colName = data.col_name || data.col || \"\";\n"
        "        if (!lookupMeta(rowName, colName)) {\n"
        "          tile.style.fill = \"#ffffff\";\n"
        "        }\n"
        "        if (tile.getAttribute(\"data-cg-bound\") === \"1\") continue;\n"
        "        tile.setAttribute(\"data-cg-bound\", \"1\");\n"
        "        tile.addEventListener(\"click\", handleTileClick);\n"
        "      }\n"
        "    }\n"
        "    var args = {\n"
        f"      root: \"#{container_id}\",\n"
        "      network_data: vizData,\n"
        "      tile_tip_callback: function(tileInfo) {\n"
        "        var rowName = tileInfo.row_name || tileInfo.row || \"\";\n"
        "        var colName = tileInfo.col_name || tileInfo.col || \"\";\n"
        "        var meta = lookupMeta(rowName, colName);\n"
        f"        var tip = document.querySelector('.{container_id}_d3-tip_tile_tip');\n"
        "        if (!tip) return;\n"
        "        tip.innerHTML = tooltipHtml(rowName, colName, meta);\n"
        "      },\n"
        "      matrix_update_callback: function() { setTimeout(attachTileHandlers, 0); }\n"
        "    };\n"
        f"{row_order_js}"
        "    var cgm = null;\n"
        "    var initStart = nowMs();\n"
        "    try { cgm = Clustergrammer(args); }\n"
        "    catch (err) {\n"
        "      var message = err && err.message ? err.message : String(err || \"unknown error\");\n"
        "      setWaitMessage(\"Heatmap failed: \" + message);\n"
        "      dumpError(\"[heatmap.debug] clustergrammer_init_failed\", err);\n"
        "      if (window.console && console.error) console.error(\"Heatmap init failed\", err);\n"
        "      return;\n"
        "    }\n"
        "    var initEnd = nowMs();\n"
        "    debugLine(\"[heatmap.debug] clustergrammer_init_ms=\" + (initEnd - initStart).toFixed(1));\n"
        "    clearWaitMessage();\n"
        "    attachTileHandlers();\n"
        "    window.addEventListener(\"resize\", function() {\n"
        "      if (cgm && typeof cgm.resize_viz === \"function\") {\n"
        "        cgm.resize_viz();\n"
        "      }\n"
        "    });\n"
        "  })();\n"
        "</script>"
    )

    html_out = "\n".join([assets_block, style_block, container_block, script_block])
    return html_out


def render_interactive_heatmap_html(
    repo_root: Path,
    run_id: str,
    input_csv: Path,
    top_k: int = 100,
    include_decoys: bool = False,
    allowed_pdb_ids: Optional[Set[str]] = None,
) -> str:
    config_path = repo_root / "config.txt"
    use_dendrogram = _parse_bool(
        _read_config_value(config_path, "USE_DENDROGRAM"), True
    )
    report_offline_assets = _parse_bool(
        _read_config_value(config_path, "REPORT_OFFLINE_ASSETS"), True
    )
    report_inline_assets = _parse_bool(
        _read_config_value(config_path, "REPORT_INLINE_ASSETS"), True
    )
    ligand_page_prefix = _read_config_value(config_path, "REPORT_LIGAND_PAGE_PREFIX")
    if ligand_page_prefix is None:
        ligand_page_prefix = "ligands/"
    ligand_page_prefix = ligand_page_prefix.strip()
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
        interactions_parquet = repo_root / "data" / run_id / "dataset" / "interactions"
        heatmap_path = repo_root / "data" / run_id / "heatmap_input.csv"
        master_path = repo_root / "data" / run_id / "master_rows.csv"
        if interactions_parquet.exists():
            resolved_input = interactions_parquet
        elif heatmap_path.exists():
            resolved_input = heatmap_path
        else:
            resolved_input = master_path
    if not resolved_input.exists():
        raise FileNotFoundError(f"Input heatmap data not found: {resolved_input}")

    rows = _load_heatmap_rows(
        resolved_input, include_decoys, allowed_pdb_ids=allowed_pdb_ids
    )
    target_display_map = _build_target_display_map(rows)
    for row in rows:
        target_id = _normalize_text(row.get("target_id"))
        if target_id:
            row["target_label"] = target_display_map.get(target_id, target_id)
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

    mat_rows = len(matrix)
    mat_cols = len(matrix[0]) if matrix else 0
    payload_start = time.perf_counter()
    _LOG.info(
        "[heatmap.payload.build.start] input=%s rows=%d cols=%d mat_shape=%dx%d top_rows=%s top_cols=%s",
        resolved_input,
        len(row_labels),
        len(col_labels),
        mat_rows,
        mat_cols,
        _sample_names(row_labels),
        _sample_names(col_labels),
    )

    cluster_html = _render_clustergrammer_html(
        repo_root=repo_root,
        run_id=run_id,
        row_labels=row_labels,
        col_labels=col_labels,
        matrix=matrix,
        agg=agg,
        use_dendrogram=use_dendrogram,
        offline_assets=report_offline_assets,
        inline_assets=report_inline_assets,
        ligand_page_prefix=ligand_page_prefix,
        payload_start=payload_start,
    )
    if cluster_html is not None:
        return cluster_html

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
                    "pose_invalid_reason_top",
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
  .heatmap-table .hm-col-label > div { transform: rotate(-60deg); transform-origin: left bottom; white-space: pre-line; line-height: 1.15; text-align: left; }
  .hm-cell { width: 18px; height: 18px; border: none; cursor: pointer; }
  .hm-cell:hover { outline: 1px solid #555; }
  .heatmap-detail { min-width: 220px; border: 1px solid #ddd; padding: 8px 10px; background: #fafafa; }
  .heatmap-detail-title { font-weight: 600; margin-bottom: 6px; }
  .heatmap-detail pre { margin: 0; white-space: pre-wrap; font-family: Consolas, monospace; font-size: 0.9em; }
  .heatmap-debug { margin-top: 12px; }
  .heatmap-debug summary { cursor: pointer; font-weight: 600; }
  .heatmap-debug pre { margin: 6px 0 0; padding: 8px; border: 1px solid #ddd; background: #f5f5f5; max-height: 240px; overflow: auto; white-space: pre-wrap; font-family: Consolas, monospace; font-size: 0.85em; }
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
<details class="heatmap-debug" id="heatmap-debug">
  <summary>Heatmap debug</summary>
  <pre id="heatmap-debug-log"></pre>
</details>
""".strip()
    script_block = (
        "<script>\n"
        "  (function() {\n"
        f"    var container = document.getElementById(\"{container_id}\");\n"
        f"    var detail = document.getElementById(\"{detail_id}\");\n"
        "    var debugEl = document.getElementById(\"heatmap-debug-log\");\n"
        "    function debugLine(message) {\n"
        "      if (!debugEl) return;\n"
        "      debugEl.textContent += String(message) + \"\\n\";\n"
        "    }\n"
        "    function nowMs() {\n"
        "      if (window.performance && typeof window.performance.now === \"function\") {\n"
        "        return window.performance.now();\n"
        "      }\n"
        "      return Date.now();\n"
        "    }\n"
        "    var prevOnError = window.onerror;\n"
        "    window.onerror = function(message, source, lineno, colno, error) {\n"
        "      var parts = [\"[heatmap.error]\", String(message || \"\"), String(source || \"\"), String(lineno || \"\"), String(colno || \"\")];\n"
        "      debugLine(parts.join(\" \"));\n"
        "      if (error && error.stack) debugLine(error.stack);\n"
        "      if (window.console && console.error) console.error(\"Heatmap error\", error || message);\n"
        "      if (typeof prevOnError === \"function\") {\n"
        "        return prevOnError(message, source, lineno, colno, error);\n"
        "      }\n"
        "      return false;\n"
        "    };\n"
        "    window.addEventListener(\"unhandledrejection\", function(evt) {\n"
        "      var reason = evt && evt.reason ? evt.reason : \"unknown\";\n"
        "      debugLine(\"[heatmap.unhandledrejection] \" + String(reason));\n"
        "      if (reason && reason.stack) debugLine(reason.stack);\n"
        "      if (window.console && console.error) console.error(\"Heatmap unhandled rejection\", reason);\n"
        "    });\n"
        "    if (!container || !detail) {\n"
        "      debugLine(\"[heatmap.debug] missing_elements container=\" + !!container + \" detail=\" + !!detail);\n"
        "      return;\n"
        "    }\n"
        f"    debugLine(\"[heatmap.debug] mode=fallback rows={len(row_labels)} cols={len(col_labels)} mat={mat_rows}x{mat_cols}\");\n"
        "    var initStart = nowMs();\n"
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
        "        if (cell.dataset.poseInvalidReasonTop) {\n"
        "          lines.push(\"pose_invalid_reason_top: \" + cell.dataset.poseInvalidReasonTop);\n"
        "        }\n"
        "        detail.textContent = lines.join(\"\\n\");\n"
        "      });\n"
        "    });\n"
        "    var initEnd = nowMs();\n"
        "    debugLine(\"[heatmap.debug] fallback_init_ms=\" + (initEnd - initStart).toFixed(1));\n"
        "  })();\n"
        "</script>"
    )
    html_out = "\n".join([style_block, table_block, script_block])
    elapsed_s = time.perf_counter() - payload_start
    _LOG.info(
        "[heatmap.payload.build.done] mode=fallback rows=%d cols=%d mat_shape=%dx%d bytes=%d elapsed_s=%.3f",
        len(row_labels),
        len(col_labels),
        mat_rows,
        mat_cols,
        len(html_out.encode("utf-8")),
        elapsed_s,
    )
    return html_out
