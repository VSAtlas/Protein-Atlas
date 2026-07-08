from __future__ import annotations

import html
import json
import logging
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, cast

from analysis.reporting.heatmap_html_assets import (
    _clustergrammer_asset_plan,  # noqa: F401
    _clustergrammer_asset_urls,  # noqa: F401
    _clustergrammer_inline_assets,  # noqa: F401
    _clustergrammer_key,
    _load_clustergrammer,  # noqa: F401
    _normalize_asset_mode,
)
from analysis.reporting.heatmap_html_clustergrammer_payload import (
    _build_clustergrammer_viz_json,  # noqa: F401
    _render_clustergrammer_html,
)
from analysis.reporting.heatmap_html_data import (
    _aggregate_rows,
    _build_matrix,
    _build_promiscuity_meta,
    _filter_empty,
    _load_heatmap_rows,
    _parse_fraction_config,
    _sample_names,
)
from analysis.reporting.heatmap_html_matrix import (
    _cluster_order,
    _format_tooltip,
    _format_value,
    _value_for_breaks,
    _value_to_color,
)
from analysis.reporting.heatmap_html_scale import (
    _build_palette,
    _format_breaks_line,
    _HEATMAP_SCALE_QUANTILE,
    _is_finite,
    _resolve_color,
    _rgb_to_hex,
    _scale_mode_text,
    compute_heatmap_breaks,
)
from analysis.reporting.heatmap_html_target_meta import (
    _build_target_display_map,
    _build_target_grouping_meta,
    _fetch_pathway_memberships_by_pdb,  # noqa: F401
    _fetch_target_uniprots_by_pdb,  # noqa: F401
)
from analysis.reporting.heatmap_report_meta import build_heatmap_report_meta
from analysis.reporting.ligand_annotation_groups import build_ligand_grouping_meta
from analysis.reporting.report_config import read_config_keys
from analysis.reporting.value_utils import normalize_text, parse_bool_or
from config.output_paths import run_output_dir

_HEATMAP_REPO_CONFIG_KEYS = (
    "USE_DENDROGRAM",
    "REPORT_OFFLINE_ASSETS",
    "REPORT_INLINE_ASSETS",
    "REPORT_ASSET_MODE",
    "BREADTH_PCT_THRESHOLD",
    "REPORT_LIGAND_PAGE_PREFIX",
    "HEATMAP_SCALE_MODE",
    "HEATMAP_SCALE_MIN",
    "HEATMAP_SCALE_MID",
    "HEATMAP_SCALE_MID2",
    "HEATMAP_SCALE_MAX",
    "HEATMAP_QUANTILE_LOW",
    "HEATMAP_QUANTILE_MID2",
    "HEATMAP_QUANTILE_HIGH",
    "HEATMAP_COLOR_MIN",
    "HEATMAP_COLOR_MID",
    "HEATMAP_COLOR_MID2",
    "HEATMAP_COLOR_MAX",
)

_LOG = logging.getLogger("heatmap-html")

# Backward-compatible re-exports for tests and monkeypatch targets.
from analysis.reporting.heatmap_html_scale import _HEATMAP_SCALE_FIXED  # noqa: E402

def render_interactive_heatmap_html(
    repo_root: Path,
    run_id: str,
    input_csv: Path,
    top_k: int = 100,
    include_decoys: bool = False,
    allowed_pdb_ids: Optional[Set[str]] = None,
    report_asset_mode: Optional[str] = None,
) -> str:
    config_path = repo_root / "config.txt"
    cfg = read_config_keys(
        config_path,
        _HEATMAP_REPO_CONFIG_KEYS,
        strip_inline_comments=True,
        keys_case_insensitive=True,
    )
    use_dendrogram = parse_bool_or(cfg.get("USE_DENDROGRAM"), default=True)
    report_offline_assets = parse_bool_or(cfg.get("REPORT_OFFLINE_ASSETS"), default=True)
    report_inline_assets = parse_bool_or(cfg.get("REPORT_INLINE_ASSETS"), default=True)
    resolved_asset_mode = _normalize_asset_mode(
        report_asset_mode or cfg.get("REPORT_ASSET_MODE")
    )
    breadth_pct_threshold = _parse_fraction_config(
        cfg.get("BREADTH_PCT_THRESHOLD"), 0.01
    )
    ligand_page_prefix = cfg.get("REPORT_LIGAND_PAGE_PREFIX")
    if ligand_page_prefix is None:
        ligand_page_prefix = "ligands/"
    ligand_page_prefix = ligand_page_prefix.strip()
    scale_config: Dict[str, Optional[str]] = {
        "HEATMAP_SCALE_MODE": cfg.get("HEATMAP_SCALE_MODE"),
        "HEATMAP_SCALE_MIN": cfg.get("HEATMAP_SCALE_MIN"),
        "HEATMAP_SCALE_MID": cfg.get("HEATMAP_SCALE_MID"),
        "HEATMAP_SCALE_MID2": cfg.get("HEATMAP_SCALE_MID2"),
        "HEATMAP_SCALE_MAX": cfg.get("HEATMAP_SCALE_MAX"),
        "HEATMAP_QUANTILE_LOW": cfg.get("HEATMAP_QUANTILE_LOW"),
        "HEATMAP_QUANTILE_MID2": cfg.get("HEATMAP_QUANTILE_MID2"),
        "HEATMAP_QUANTILE_HIGH": cfg.get("HEATMAP_QUANTILE_HIGH"),
    }

    color_low = _resolve_color(cfg.get("HEATMAP_COLOR_MIN"), "#2166ac")
    color_mid = _resolve_color(cfg.get("HEATMAP_COLOR_MID"), "#f7f7f7")
    color_mid2 = _resolve_color(cfg.get("HEATMAP_COLOR_MID2"), "#f4a582")
    color_high = _resolve_color(cfg.get("HEATMAP_COLOR_MAX"), "#b2182b")
    color_range_hex = [
        _rgb_to_hex(color_low),
        _rgb_to_hex(color_mid),
        _rgb_to_hex(color_mid2),
        _rgb_to_hex(color_high),
    ]

    resolved_input = input_csv
    if not resolved_input.exists():
        data_dir = run_output_dir(repo_root, "data", run_id)
        interactions_parquet = data_dir / "dataset" / "interactions"
        heatmap_path = data_dir / "heatmap_input.csv"
        master_path = data_dir / "master_rows.csv"
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
        target_id = normalize_text(row.get("target_id"))
        if target_id:
            row["target_label"] = target_display_map.get(target_id, target_id)
    full_agg, ligand_max = _aggregate_rows(rows)
    if not ligand_max:
        raise ValueError("No ligands available after filtering")
    row_promiscuity_meta, target_promiscuity_meta = _build_promiscuity_meta(
        full_agg, pct_threshold=breadth_pct_threshold
    )
    all_target_levels = sorted({target for _, target in full_agg.keys()})
    if not all_target_levels:
        raise ValueError("Heatmap matrix is empty after filtering")

    top_k = int(top_k) if top_k and top_k > 0 else 100
    ranked_ligands = sorted(ligand_max.items(), key=lambda item: (-item[1], item[0]))
    ranked_ligand_labels = [name for name, _ in ranked_ligands]
    display_top_k = min(top_k, len(ranked_ligand_labels))
    top_ligands = sorted(ranked_ligand_labels[:display_top_k])

    display_agg = {
        key: value
        for key, value in full_agg.items()
        if key[0] in top_ligands
    }
    full_matrix = _build_matrix(full_agg, ranked_ligand_labels, all_target_levels)
    full_matrix, full_row_labels, full_col_labels = _filter_empty(
        full_matrix, ranked_ligand_labels, all_target_levels, drop_empty_cols=False
    )
    if not full_row_labels or not full_col_labels:
        raise ValueError("Heatmap matrix is empty after filtering")
    full_row_label_set = set(full_row_labels)
    ranked_ligand_labels = [label for label in ranked_ligand_labels if label in full_row_label_set]

    target_group_meta, target_group_summary = _build_target_grouping_meta(
        repo_root,
        rows,
        full_col_labels,
    )
    ligand_group_meta, ligand_group_summary = build_ligand_grouping_meta(
        repo_root,
        rows,
        full_row_labels,
    )
    analysis_payload: Dict[str, Any] = {
        "default_threshold_pct_rank": breadth_pct_threshold,
        "default_top_k": display_top_k,
        "ranked_rows": [_clustergrammer_key(label) for label in ranked_ligand_labels],
        "default_pose_valid_internal": True,
        "target_organization": target_group_summary,
        "ligand_organization": ligand_group_summary,
    }
    analysis_payload.update(
        build_heatmap_report_meta(
            repo_root=repo_root,
            run_id=run_id,
            asset_mode=resolved_asset_mode,
            rows=rows,
            full_row_labels=full_row_labels,
            full_col_labels=full_col_labels,
            display_top_k=display_top_k,
            target_group_meta=target_group_meta,
        )
    )

    matrix = _build_matrix(display_agg, top_ligands, all_target_levels)
    matrix, row_labels, col_labels = _filter_empty(
        matrix, top_ligands, all_target_levels, drop_empty_cols=False
    )
    if not row_labels or not col_labels:
        raise ValueError("Heatmap matrix is empty after filtering")

    finite_values = [
        cast(float, value)
        for row in matrix
        for value in row
        if _is_finite(cast(Optional[float], value))
    ]
    scale_result = compute_heatmap_breaks(finite_values, scale_config)
    chosen_breaks = cast(Dict[str, float], scale_result["chosen_breaks"])
    fixed_breaks = cast(Dict[str, float], scale_result["fixed_breaks"])
    quantile_breaks = cast(Optional[Dict[str, float]], scale_result["quantile_breaks"])
    scale_stats = cast(Dict[str, Any], scale_result["stats"])
    scale_mode = str(scale_result["mode"])
    scale_min = chosen_breaks["scale_min"]
    scale_mid = chosen_breaks["scale_mid"]
    scale_mid2 = chosen_breaks["scale_mid2"]
    scale_max = chosen_breaks["scale_max"]

    scale_payload: Dict[str, Any] = {
        "scale_mode": scale_mode,
        "requested_scale_mode": scale_result["requested_mode"],
        "chosen_breaks": chosen_breaks,
        "fixed_breaks": fixed_breaks,
        "quantile_breaks": quantile_breaks,
        "n_finite": scale_stats.get("n_finite"),
        "q_low": scale_stats.get("q_low"),
        "q_mid2": scale_stats.get("q_mid2"),
        "q_high": scale_stats.get("q_high"),
        "fallback_reason": scale_stats.get("fallback_reason"),
        "quantile_widened": scale_stats.get("quantile_widened"),
        "color_range": color_range_hex,
    }

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
        row_labels=full_row_labels,
        col_labels=full_col_labels,
        matrix=full_matrix,
        agg=full_agg,
        row_promiscuity_meta=row_promiscuity_meta,
        target_promiscuity_meta=target_promiscuity_meta,
        target_group_meta=target_group_meta,
        ligand_group_meta=ligand_group_meta,
        analysis_payload=analysis_payload,
        use_dendrogram=use_dendrogram,
        offline_assets=report_offline_assets,
        inline_assets=report_inline_assets,
        asset_mode=resolved_asset_mode,
        ligand_page_prefix=ligand_page_prefix,
        scale_payload=scale_payload,
        color_scale_range=color_range_hex,
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
            cell = display_agg.get((ligand, target))
            target_pdb = str(target).split("\n", 1)[0]
            z_val_raw = cell["z_selected"] if cell else None
            z_val_fixed = _value_for_breaks(z_val_raw, fixed_breaks)
            z_val_quantile = _value_for_breaks(z_val_raw, quantile_breaks)
            z_val_display = z_val_fixed
            if scale_mode == _HEATMAP_SCALE_QUANTILE and z_val_quantile is not None:
                z_val_display = z_val_quantile
            if z_val_display is None:
                z_val_display = z_val_raw
            color = _value_to_color(
                z_val_display, heatmap_breaks, heatmap_colors, scale_min, scale_max
            )
            tooltip_lines = [
                f"ligand: {ligand}",
                f"target: {target_pdb}",
                f"Atlas score: {_format_value(z_val_display)}",
            ]
            if _format_value(z_val_raw) != _format_value(z_val_display):
                tooltip_lines.append(f"Raw Atlas score: {_format_value(z_val_raw)}")
            data_attrs = {
                "ligand": ligand,
                "target": target,
                "target_pdb": target_pdb,
                "z_selected": _format_value(z_val_display),
                "z_selected_raw": _format_value(z_val_raw),
                "z_selected_fixed": _format_value(z_val_fixed),
            }
            if z_val_quantile is not None:
                data_attrs["z_selected_quantile"] = _format_value(z_val_quantile)
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
    scale_json_id = f"heatmap-scale-meta-{safe_run_id}"
    scale_mode_label_id = f"heatmap-scale-mode-{safe_run_id}"
    scale_breaks_label_id = f"heatmap-scale-breaks-{safe_run_id}"
    scale_fixed_label_id = f"heatmap-scale-fixed-{safe_run_id}"
    scale_toggle_fixed_id = f"heatmap-scale-toggle-fixed-{safe_run_id}"
    scale_toggle_quantile_id = f"heatmap-scale-toggle-quantile-{safe_run_id}"
    mode_label = _scale_mode_text(scale_mode)
    chosen_breaks_line = _format_breaks_line(chosen_breaks)
    fixed_breaks_line = _format_breaks_line(fixed_breaks)
    requested_scale_mode = str(scale_payload.get("requested_scale_mode") or _HEATMAP_SCALE_FIXED)
    fallback_reason = str(scale_payload.get("fallback_reason") or "")
    quantile_available = quantile_breaks is not None
    quantile_checked = (
        requested_scale_mode == _HEATMAP_SCALE_QUANTILE and quantile_available
    )
    fixed_checked = not quantile_checked
    fixed_reference_html = (
        f"<div class=\"heatmap-scale-fixed\" id=\"{scale_fixed_label_id}\">Fixed reference breaks: {html.escape(fixed_breaks_line)}</div>"
        if scale_mode == _HEATMAP_SCALE_QUANTILE
        else f"<div class=\"heatmap-scale-fixed\" id=\"{scale_fixed_label_id}\"></div>"
    )
    fallback_notice_html = (
        f"<div class=\"heatmap-scale-fixed\">Quantile fallback: {html.escape(fallback_reason)}</div>"
        if requested_scale_mode == _HEATMAP_SCALE_QUANTILE
        and scale_mode != _HEATMAP_SCALE_QUANTILE
        and fallback_reason
        else ""
    )
    toggle_line = (
        "<div class=\"heatmap-scale-toggle\">"
        "<span>Scale display:</span>"
        f"<label><input type=\"radio\" name=\"hm-scale-mode-{safe_run_id}\" id=\"{scale_toggle_fixed_id}\" value=\"fixed\"{' checked' if fixed_checked else ''}> Fixed (absolute)</label>"
        f"<label><input type=\"radio\" name=\"hm-scale-mode-{safe_run_id}\" id=\"{scale_toggle_quantile_id}\" value=\"quantile\"{' checked' if quantile_checked else ''}{' disabled' if not quantile_available else ''}> Quantile (within-run){' (unavailable)' if not quantile_available else ''}</label>"
        "</div>"
    )
    scale_json_text = json.dumps(scale_payload, ensure_ascii=True).replace("</", "<\\/")
    style_block = """
<style>
  .heatmap-scale-meta { margin: 0 0 10px 0; border: 1px solid #d9d9d9; background: #fffde8; padding: 8px 10px; }
  .heatmap-scale-mode { font-weight: 700; }
  .heatmap-scale-breaks { margin-top: 4px; font-family: Consolas, monospace; font-size: 0.9em; }
  .heatmap-scale-toggle { margin-top: 6px; font-size: 0.9em; display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }
  .heatmap-scale-toggle label { display: inline-flex; align-items: center; gap: 4px; }
  .heatmap-scale-fixed { margin-top: 4px; font-family: Consolas, monospace; font-size: 0.85em; color: #444; }
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
<div class="heatmap-scale-meta heatmap-controls-card" data-heatmap-controls="card">
  <div class="heatmap-scale-mode" id="{scale_mode_label_id}">{html.escape(mode_label)}</div>
  <div class="heatmap-scale-breaks" id="{scale_breaks_label_id}">Breaks used: {html.escape(chosen_breaks_line)}</div>
  {toggle_line}
  {fixed_reference_html}
  {fallback_notice_html}
</div>
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
<details class="heatmap-debug heatmap-debug-panel" id="heatmap-debug">
  <summary>Heatmap debug</summary>
  <pre id="heatmap-debug-log"></pre>
</details>
<script type="application/json" id="{scale_json_id}">{scale_json_text}</script>
""".strip()
    script_block = (
        "<script>\n"
        "  (function() {\n"
        f"    var container = document.getElementById(\"{container_id}\");\n"
        f"    var detail = document.getElementById(\"{detail_id}\");\n"
        f"    var scaleEl = document.getElementById(\"{scale_json_id}\");\n"
        f"    var modeLabelEl = document.getElementById(\"{scale_mode_label_id}\");\n"
        f"    var breaksLabelEl = document.getElementById(\"{scale_breaks_label_id}\");\n"
        f"    var fixedLabelEl = document.getElementById(\"{scale_fixed_label_id}\");\n"
        f"    var toggleFixedEl = document.getElementById(\"{scale_toggle_fixed_id}\");\n"
        f"    var toggleQuantileEl = document.getElementById(\"{scale_toggle_quantile_id}\");\n"
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
        "    if (!container || !detail || !scaleEl) {\n"
        "      debugLine(\"[heatmap.debug] missing_elements container=\" + !!container + \" detail=\" + !!detail + \" scale=\" + !!scaleEl);\n"
        "      return;\n"
        "    }\n"
        f"    debugLine(\"[heatmap.debug] mode=fallback rows={len(row_labels)} cols={len(col_labels)} mat={mat_rows}x{mat_cols}\");\n"
        "    function asNumber(value) {\n"
        "      var parsed = Number(value);\n"
        "      return Number.isFinite(parsed) ? parsed : null;\n"
        "    }\n"
        "    function normalizeBreaks(raw) {\n"
        "      if (!raw || typeof raw !== \"object\") return null;\n"
        "      var b = {\n"
        "        scale_min: asNumber(raw.scale_min),\n"
        "        scale_mid: asNumber(raw.scale_mid),\n"
        "        scale_mid2: asNumber(raw.scale_mid2),\n"
        "        scale_max: asNumber(raw.scale_max)\n"
        "      };\n"
        "      if (b.scale_min == null || b.scale_mid == null || b.scale_mid2 == null || b.scale_max == null) {\n"
        "        return null;\n"
        "      }\n"
        "      if (!(b.scale_min < b.scale_mid && b.scale_mid < b.scale_mid2 && b.scale_mid2 < b.scale_max)) {\n"
        "        return null;\n"
        "      }\n"
        "      return b;\n"
        "    }\n"
        "    function parseHexColor(value) {\n"
        "      if (!value || typeof value !== \"string\") return null;\n"
        "      var hex = value.trim();\n"
        "      if (hex.charAt(0) === \"#\") hex = hex.slice(1);\n"
        "      if (hex.length === 3) {\n"
        "        hex = hex.charAt(0) + hex.charAt(0) + hex.charAt(1) + hex.charAt(1) + hex.charAt(2) + hex.charAt(2);\n"
        "      }\n"
        "      if (!/^[0-9a-fA-F]{6}$/.test(hex)) return null;\n"
        "      return {\n"
        "        r: parseInt(hex.slice(0, 2), 16),\n"
        "        g: parseInt(hex.slice(2, 4), 16),\n"
        "        b: parseInt(hex.slice(4, 6), 16)\n"
        "      };\n"
        "    }\n"
        "    function rgbToHex(rgb) {\n"
        "      function clampChannel(v) {\n"
        "        var n = Math.max(0, Math.min(255, Math.round(v)));\n"
        "        var s = n.toString(16);\n"
        "        return s.length === 1 ? \"0\" + s : s;\n"
        "      }\n"
        "      return \"#\" + clampChannel(rgb.r) + clampChannel(rgb.g) + clampChannel(rgb.b);\n"
        "    }\n"
        "    function lerpColor(a, b, t) {\n"
        "      var clamped = Math.max(0, Math.min(1, t));\n"
        "      return {\n"
        "        r: a.r + (b.r - a.r) * clamped,\n"
        "        g: a.g + (b.g - a.g) * clamped,\n"
        "        b: a.b + (b.b - a.b) * clamped\n"
        "      };\n"
        "    }\n"
        "    function modeLabelText(mode) {\n"
        "      if (mode === \"quantile\") {\n"
        "        return \"Scale mode: QUANTILE (within-run; colors are relative)\";\n"
        "      }\n"
        "      return \"Scale mode: FIXED (absolute)\";\n"
        "    }\n"
        "    function breaksLineText(breaks) {\n"
        "      if (!breaks) return \"Breaks used: unavailable\";\n"
        "      return \"Breaks used: min=\" + breaks.scale_min.toPrecision(6)\n"
        "        + \", mid=\" + breaks.scale_mid.toPrecision(6)\n"
        "        + \", mid2=\" + breaks.scale_mid2.toPrecision(6)\n"
        "        + \", max=\" + breaks.scale_max.toPrecision(6);\n"
        "    }\n"
        "    function fixedReferenceText(mode, fixedBreaks) {\n"
        "      if (mode !== \"quantile\" || !fixedBreaks) return \"\";\n"
        "      return \"Fixed reference breaks: min=\" + fixedBreaks.scale_min.toPrecision(6)\n"
        "        + \", mid=\" + fixedBreaks.scale_mid.toPrecision(6)\n"
        "        + \", mid2=\" + fixedBreaks.scale_mid2.toPrecision(6)\n"
        "        + \", max=\" + fixedBreaks.scale_max.toPrecision(6);\n"
        "    }\n"
        "    var scaleMeta = {};\n"
        "    try { scaleMeta = JSON.parse(scaleEl.textContent || \"{}\"); }\n"
        "    catch (err) {\n"
        "      debugLine(\"[heatmap.debug] scale_parse_failed \" + String(err && err.message ? err.message : err));\n"
        "      return;\n"
        "    }\n"
        "    var fixedBreaks = normalizeBreaks(scaleMeta.fixed_breaks);\n"
        "    var quantileBreaks = normalizeBreaks(scaleMeta.quantile_breaks);\n"
        "    var chosenBreaks = normalizeBreaks(scaleMeta.chosen_breaks) || fixedBreaks;\n"
        "    var requestedMode = scaleMeta.requested_scale_mode === \"quantile\" ? \"quantile\" : \"fixed\";\n"
        "    var quantileAvailable = !!quantileBreaks;\n"
        "    var activeMode = requestedMode === \"quantile\" && quantileAvailable ? \"quantile\" : \"fixed\";\n"
        "    var activeBreaks = null;\n"
        "    var colorRange = Array.isArray(scaleMeta.color_range) ? scaleMeta.color_range : [];\n"
        "    var colorStops = colorRange.map(parseHexColor);\n"
        "    function resolveActiveBreaks() {\n"
        "      if (activeMode === \"quantile\") {\n"
        "        return quantileBreaks || chosenBreaks || fixedBreaks;\n"
        "      }\n"
        "      return fixedBreaks || chosenBreaks;\n"
        "    }\n"
        "    function updateScaleBanner() {\n"
        "      if (modeLabelEl) modeLabelEl.textContent = modeLabelText(activeMode);\n"
        "      if (breaksLabelEl) breaksLabelEl.textContent = breaksLineText(activeBreaks);\n"
        "      if (fixedLabelEl) fixedLabelEl.textContent = fixedReferenceText(activeMode, fixedBreaks);\n"
        "    }\n"
        "    function colorForValue(value) {\n"
        "      var numeric = asNumber(value);\n"
        "      if (numeric == null || !activeBreaks) return \"white\";\n"
        "      if (!colorStops || colorStops.length < 4 || !colorStops[0] || !colorStops[1] || !colorStops[2] || !colorStops[3]) {\n"
        "        return \"white\";\n"
        "      }\n"
        "      var x0 = activeBreaks.scale_min;\n"
        "      var x1 = activeBreaks.scale_mid;\n"
        "      var x2 = activeBreaks.scale_mid2;\n"
        "      var x3 = activeBreaks.scale_max;\n"
        "      if (numeric <= x0) return colorRange[0] || \"white\";\n"
        "      if (numeric >= x3) return colorRange[3] || \"white\";\n"
        "      if (numeric <= x1) {\n"
        "        return rgbToHex(lerpColor(colorStops[0], colorStops[1], (numeric - x0) / (x1 - x0)));\n"
        "      }\n"
        "      if (numeric <= x2) {\n"
        "        return rgbToHex(lerpColor(colorStops[1], colorStops[2], (numeric - x1) / (x2 - x1)));\n"
        "      }\n"
        "      return rgbToHex(lerpColor(colorStops[2], colorStops[3], (numeric - x2) / (x3 - x2)));\n"
        "    }\n"
        "    function activeCellValueText(cell) {\n"
        "      if (!cell || !cell.dataset) return \"\";\n"
        "      if (activeMode === \"quantile\" && cell.dataset.zSelectedQuantile) {\n"
        "        return cell.dataset.zSelectedQuantile;\n"
        "      }\n"
        "      if (activeMode === \"fixed\" && cell.dataset.zSelectedFixed) {\n"
        "        return cell.dataset.zSelectedFixed;\n"
        "      }\n"
        "      return cell.dataset.zSelected || cell.dataset.tSelected || \"\";\n"
        "    }\n"
        "    function scalePosition(value, breaks) {\n"
        "      var numeric = asNumber(value);\n"
        "      if (numeric == null || !breaks) return null;\n"
        "      var x0 = breaks.scale_min;\n"
        "      var x1 = breaks.scale_mid;\n"
        "      var x2 = breaks.scale_mid2;\n"
        "      var x3 = breaks.scale_max;\n"
        "      if (!(x0 < x1 && x1 < x2 && x2 < x3)) return null;\n"
        "      if (numeric <= x0) return 0;\n"
        "      if (numeric >= x3) return 3;\n"
        "      if (numeric <= x1) return (numeric - x0) / (x1 - x0);\n"
        "      if (numeric <= x2) return 1 + (numeric - x1) / (x2 - x1);\n"
        "      return 2 + (numeric - x2) / (x3 - x2);\n"
        "    }\n"
        "    var initStart = nowMs();\n"
        "    var cells = container.querySelectorAll(\".hm-cell\");\n"
        "    var selectedCell = null;\n"
        "    function applyScaleMode() {\n"
        "      activeBreaks = resolveActiveBreaks();\n"
        "      updateScaleBanner();\n"
        "      cells.forEach(function(cell) {\n"
        "        var value = activeCellValueText(cell);\n"
        "        cell.style.backgroundColor = colorForValue(value);\n"
        "      });\n"
        "      if (selectedCell) {\n"
        "        selectedCell.click();\n"
        "      }\n"
        "      debugLine(\"[heatmap.debug] scale_mode=\" + activeMode + \" breaks=\" + breaksLineText(activeBreaks));\n"
        "    }\n"
        "    if (toggleFixedEl) {\n"
        "      toggleFixedEl.checked = activeMode === \"fixed\";\n"
        "      toggleFixedEl.addEventListener(\"change\", function() {\n"
        "        if (!toggleFixedEl.checked) return;\n"
        "        activeMode = \"fixed\";\n"
        "        applyScaleMode();\n"
        "      });\n"
        "    }\n"
        "    if (toggleQuantileEl) {\n"
        "      toggleQuantileEl.disabled = !quantileAvailable;\n"
        "      toggleQuantileEl.checked = activeMode === \"quantile\";\n"
        "      toggleQuantileEl.addEventListener(\"change\", function() {\n"
        "        if (!toggleQuantileEl.checked || !quantileAvailable) return;\n"
        "        activeMode = \"quantile\";\n"
        "        applyScaleMode();\n"
        "      });\n"
        "    }\n"
        "    cells.forEach(function(cell) {\n"
        "      cell.addEventListener(\"click\", function() {\n"
        "        selectedCell = cell;\n"
        "        var displayValue = activeCellValueText(cell);\n"
        "        var lines = [\n"
        "          \"ligand: \" + (cell.dataset.ligand || \"\"),\n"
        "          \"target: \" + (cell.dataset.targetPdb || cell.dataset.target || \"\"),\n"
        "          \"Atlas score: \" + (displayValue || \"\")\n"
        "        ];\n"
        "        if (cell.dataset.zSelectedRaw && displayValue && cell.dataset.zSelectedRaw !== displayValue) {\n"
        "          lines.push(\"Raw Atlas score: \" + cell.dataset.zSelectedRaw);\n"
        "        }\n"
        "        if (cell.dataset.pctRank) lines.push(\"pct_rank: \" + cell.dataset.pctRank);\n"
        "        if (cell.dataset.poseInvalidReasonTop) {\n"
        "          lines.push(\"pose_invalid_reason_top: \" + cell.dataset.poseInvalidReasonTop);\n"
        "        }\n"
        "        detail.textContent = lines.join(\"\\n\");\n"
        "      });\n"
        "    });\n"
        "    applyScaleMode();\n"
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
