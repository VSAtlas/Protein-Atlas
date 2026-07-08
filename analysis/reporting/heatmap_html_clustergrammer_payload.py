from __future__ import annotations

import copy
import html
import json
import logging
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, cast

from analysis.reporting.heatmap_html_runtime import heatmap_call
from analysis.reporting.heatmap_html_assets import (
    _clustergrammer_key,
    _warn_clustergrammer_failure,
)
from analysis.reporting.heatmap_html_matrix import (
    _format_percent,
    _format_value,
    _value_for_breaks,
)
from analysis.reporting.heatmap_html_scale import (
    _HEATMAP_SCALE_FIXED,
    _HEATMAP_SCALE_QUANTILE,
    _format_breaks_line,
    _is_finite,
    _scale_mode_text,
)
from analysis.reporting.heatmap_html_clustergrammer_template import (
    render_clustergrammer_client_html,
)
from analysis.reporting.heatmap_html_target_meta import _target_expression_report_fields
from analysis.reporting.side_effect_overlap_meta import (
    ligand_side_effect_report_fields,
    target_side_effect_report_fields,
)
from analysis.reporting.value_utils import (
    normalize_side_effect_label,
    normalize_text,
)
from analysis.reporting.target_annotation_groups import ADME_DEFAULT, SAFETY_DEFAULT
from analysis.target_ids import pdb_id_from_target_id

_LOG = logging.getLogger("heatmap-html")

def _build_clustergrammer_viz_json(
    matrix: List[List[Optional[float]]],
    row_labels: List[str],
    col_labels: List[str],
    use_dendrogram: bool,
) -> Dict[str, Any]:
    loaded = heatmap_call("_load_clustergrammer")
    if loaded is None:
        raise ImportError("clustergrammer not available")
    Network, pd = loaded
    df = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
    df_filled = df
    if df.isna().any().any():
        # Keep missing cells at the neutral midpoint instead of painting them
        # with the most negative observed value in the matrix.
        df_filled = df.fillna(0.0)
    net = Network()
    net.load_df(df_filled)
    net.cluster(run_clustering=use_dendrogram, dendro=use_dendrogram)
    viz_json = net.export_net_json("viz")
    if isinstance(viz_json, str):
        return json.loads(viz_json)
    return cast(Dict[str, Any], viz_json)


def _build_clustergrammer_mode_mat(
    matrix: List[List[Optional[float]]],
    row_labels: List[str],
    col_labels: List[str],
    ordered_row_labels: List[str],
    ordered_col_labels: List[str],
    breaks: Optional[Dict[str, float]],
) -> List[List[float]]:
    row_index = {label: idx for idx, label in enumerate(row_labels)}
    col_index = {label: idx for idx, label in enumerate(col_labels)}
    out: List[List[float]] = []
    for row_label in ordered_row_labels:
        source_row = matrix[row_index[row_label]]
        out_row: List[float] = []
        for col_label in ordered_col_labels:
            raw_value = source_row[col_index[col_label]]
            numeric = 0.0 if not _is_finite(raw_value) else cast(float, raw_value)
            selected = _value_for_breaks(numeric, breaks)
            out_row.append(float(selected) if selected is not None else 0.0)
        out.append(out_row)
    return out


def _ordered_clustergrammer_axis_labels(
    base_names: List[str],
    label_lookup: Set[str],
    trailing_labels: List[str],
) -> List[str]:
    ordered = [label for label in base_names if label in label_lookup]
    seen = set(ordered)
    for label in trailing_labels:
        if label not in seen:
            ordered.append(label)
            seen.add(label)
    return ordered


def _clone_clustergrammer_viz_for_mode(
    base_viz: Dict[str, Any],
    mode_mat: List[List[float]],
    *,
    row_nodes: Optional[List[Dict[str, Any]]] = None,
    col_nodes: Optional[List[Dict[str, Any]]] = None,
    views: Optional[List[Any]] = None,
) -> Dict[str, Any]:
    cloned = {
        "row_nodes": copy.deepcopy(
            row_nodes if row_nodes is not None else base_viz.get("row_nodes", [])
        ),
        "col_nodes": copy.deepcopy(
            col_nodes if col_nodes is not None else base_viz.get("col_nodes", [])
        ),
        "links": [],
        "mat": mode_mat,
        "cat_colors": copy.deepcopy(base_viz.get("cat_colors", {})),
        "views": copy.deepcopy(views if views is not None else base_viz.get("views", [])),
    }
    return cloned


def _render_clustergrammer_html(
    repo_root: Path,
    run_id: str,
    row_labels: List[str],
    col_labels: List[str],
    matrix: List[List[Optional[float]]],
    agg: Dict[Tuple[str, str], Dict[str, Any]],
    row_promiscuity_meta: Dict[str, Dict[str, float]],
    target_promiscuity_meta: Dict[str, Dict[str, float]],
    target_group_meta: Dict[str, Dict[str, Any]],
    ligand_group_meta: Dict[str, Dict[str, Any]],
    analysis_payload: Dict[str, Any],
    use_dendrogram: bool,
    offline_assets: bool,
    inline_assets: bool,
    asset_mode: str,
    ligand_page_prefix: str,
    scale_payload: Dict[str, Any],
    color_scale_range: List[str],
    payload_start: Optional[float] = None,
) -> Optional[str]:
    build_start = payload_start if payload_start is not None else time.perf_counter()
    try:
        base_viz_json = heatmap_call(
            "_build_clustergrammer_viz_json",
            matrix,
            row_labels,
            col_labels,
            use_dendrogram,
        )
    except Exception as exc:
        _warn_clustergrammer_failure(exc)
        return None

    chosen_breaks = cast(Dict[str, float], scale_payload.get("chosen_breaks", {}))
    fixed_breaks = cast(Dict[str, float], scale_payload.get("fixed_breaks", {}))
    quantile_breaks = cast(Optional[Dict[str, float]], scale_payload.get("quantile_breaks"))
    scale_mode = str(scale_payload.get("scale_mode") or _HEATMAP_SCALE_FIXED)
    requested_mode = str(
        scale_payload.get("requested_scale_mode") or _HEATMAP_SCALE_FIXED
    )
    fallback_reason = str(scale_payload.get("fallback_reason") or "")

    row_set = set(row_labels)
    col_set = set(col_labels)
    cell_meta: Dict[str, Dict[str, str]] = {}
    for (ligand, target), cell in agg.items():
        if ligand not in row_set or target not in col_set:
            continue
        key = f"{_clustergrammer_key(ligand)}||{_clustergrammer_key(target)}"
        raw_value = cast(Optional[float], cell.get("z_selected"))
        fixed_value = _value_for_breaks(raw_value, fixed_breaks)
        quantile_value = _value_for_breaks(raw_value, quantile_breaks)
        selected_value = fixed_value
        if scale_mode == _HEATMAP_SCALE_QUANTILE and quantile_value is not None:
            selected_value = quantile_value
        if selected_value is None:
            selected_value = raw_value
        meta_entry: Dict[str, str] = {
            "z_selected": _format_value(selected_value),
            "z_selected_raw": _format_value(raw_value),
            "z_selected_fixed": _format_value(fixed_value),
            "ligand_raw": ligand,
            "target_raw": normalize_text(cell.get("target_raw")) or target,
        }
        if quantile_value is not None:
            meta_entry["z_selected_quantile"] = _format_value(quantile_value)
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

    target_raw_by_label: Dict[str, str] = {}
    for (_ligand, target_label), cell in agg.items():
        if target_label not in col_set:
            continue
        raw_target = normalize_text(cell.get("target_raw"))
        if raw_target and target_label not in target_raw_by_label:
            target_raw_by_label[target_label] = raw_target
    row_label_meta: Dict[str, Dict[str, Any]] = {}
    for ligand_label in row_labels:
        promiscuity_meta = row_promiscuity_meta.get(ligand_label, {})
        grouping_meta = ligand_group_meta.get(ligand_label, {})
        hit_count = int(promiscuity_meta.get("significant_hits", 0.0) or 0.0)
        partner_total = int(promiscuity_meta.get("partner_total", 0.0) or 0.0)
        row_label_meta[_clustergrammer_key(ligand_label)] = {
            "ligand": ligand_label,
            "library": normalize_text(grouping_meta.get("library")),
            "chemotype_primary": normalize_text(grouping_meta.get("chemotype_primary"))
            or "Other / Unassigned",
            "chemotype_series": normalize_text(grouping_meta.get("chemotype_series"))
            or "Other / Unassigned",
            "drug_effect_primary": normalize_text(
                grouping_meta.get("drug_effect_primary")
            )
            or "Other / Unassigned",
            "drug_effect_memberships": cast(
                List[str], grouping_meta.get("drug_effect_memberships") or []
            ),
            "motif_tags": cast(List[str], grouping_meta.get("motif_tags") or []),
            "ligand_annotation_confidence": normalize_text(
                grouping_meta.get("ligand_annotation_confidence")
            )
            or "unassigned",
            "ligand_annotation_sources": cast(
                List[str], grouping_meta.get("ligand_annotation_sources") or []
            ),
            "scaffold_key": normalize_text(grouping_meta.get("scaffold_key")),
            "pairability_key": normalize_text(grouping_meta.get("pairability_key")),
            "motif_signature": normalize_text(grouping_meta.get("motif_signature")),
            **ligand_side_effect_report_fields(grouping_meta),
            "promiscuity_z": _format_value(promiscuity_meta.get("promiscuity_z")),
            "significant_hits": str(hit_count),
            "partner_total": str(partner_total),
            "significant_hits_text": f"{hit_count} / {partner_total}",
            "hit_fraction": _format_value(promiscuity_meta.get("hit_fraction")),
            "hit_fraction_pct": _format_percent(promiscuity_meta.get("hit_fraction")),
            "entropy": _format_value(promiscuity_meta.get("entropy")),
            "threshold_pct_rank": _format_value(
                promiscuity_meta.get("threshold_pct_rank")
            ),
        }
    col_label_meta: Dict[str, Dict[str, Any]] = {}
    for target_label in col_labels:
        raw_target = target_raw_by_label.get(target_label, "")
        lines = [line.strip() for line in str(target_label).split("\n") if line.strip()]
        pdb = lines[0] if lines else pdb_id_from_target_id(raw_target)
        if not pdb:
            pdb = pdb_id_from_target_id(raw_target)
        target_name = lines[1] if len(lines) > 1 else ""
        if not target_name:
            target_name = pdb
        variant_ph = lines[2] if len(lines) > 2 else ""
        promiscuity_meta = target_promiscuity_meta.get(target_label, {})
        grouping_meta = target_group_meta.get(target_label, {})
        hit_count = int(promiscuity_meta.get("significant_hits", 0.0) or 0.0)
        partner_total = int(promiscuity_meta.get("partner_total", 0.0) or 0.0)
        col_label_meta[_clustergrammer_key(target_label)] = {
            "pdb": pdb,
            "target_name": target_name,
            "variant_ph": variant_ph,
            "uniprot_accessions": cast(
                List[str], grouping_meta.get("uniprot_accessions") or []
            ),
            "protein_family": normalize_text(grouping_meta.get("protein_family"))
            or "Unassigned",
            "adme_category": normalize_text(grouping_meta.get("adme_category"))
            or ADME_DEFAULT,
            "safety_buckets": cast(List[str], grouping_meta.get("safety_buckets") or []),
            "primary_display_safety": normalize_side_effect_label(
                grouping_meta.get("primary_display_safety")
            )
            or SAFETY_DEFAULT,
            "secondary_safety_buckets": cast(
                List[str], grouping_meta.get("secondary_safety_buckets") or []
            ),
            "direct_safety_buckets": cast(
                List[str], grouping_meta.get("direct_safety_buckets") or []
            ),
            "drug_ae_buckets": cast(List[str], grouping_meta.get("drug_ae_buckets") or []),
            "safety_confidence": normalize_text(grouping_meta.get("safety_confidence"))
            or "unassigned",
            "safety_sources": cast(List[str], grouping_meta.get("safety_sources") or []),
            "safety_bucket_scores": cast(
                Dict[str, Any], grouping_meta.get("safety_bucket_scores") or {}
            ),
            "safety_evidence_summary": normalize_side_effect_label(
                grouping_meta.get("safety_evidence_summary")
            ),
            **_target_expression_report_fields(grouping_meta),
            **target_side_effect_report_fields(grouping_meta),
            "safety_evidence_terms": cast(
                List[str], grouping_meta.get("safety_evidence_terms") or []
            ),
            "direct_liability_examples": cast(
                List[str], grouping_meta.get("direct_liability_examples") or []
            ),
            "raw_adverse_event_examples": cast(
                List[str], grouping_meta.get("raw_adverse_event_examples") or []
            ),
            "pathway_memberships": cast(
                List[str], grouping_meta.get("pathway_memberships") or []
            ),
            "primary_display_pathway": normalize_text(
                grouping_meta.get("primary_display_pathway")
            )
            or "Unassigned",
            "promiscuity_z": _format_value(promiscuity_meta.get("promiscuity_z")),
            "significant_hits": str(hit_count),
            "partner_total": str(partner_total),
            "significant_hits_text": f"{hit_count} / {partner_total}",
            "hit_fraction": _format_value(promiscuity_meta.get("hit_fraction")),
            "hit_fraction_pct": _format_percent(promiscuity_meta.get("hit_fraction")),
            "entropy": _format_value(promiscuity_meta.get("entropy")),
            "threshold_pct_rank": _format_value(
                promiscuity_meta.get("threshold_pct_rank")
            ),
        }

    safe_run_id = re.sub(r"[^A-Za-z0-9_-]+", "-", str(run_id)).strip("-") or "run"
    asset_plan = heatmap_call(
        "_clustergrammer_asset_plan",
        repo_root=repo_root,
        run_id=run_id,
        offline_assets=offline_assets,
        inline_assets=inline_assets,
        asset_mode=asset_mode,
    )
    assets_block = str(asset_plan["html"])
    asset_mode = str(asset_plan["mode"])
    asset_notices = [
        str(notice).strip() for notice in cast(List[str], asset_plan.get("notices", []))
        if str(notice).strip()
    ]
    asset_notice_html = ""
    if asset_mode == "cdn":
        reason_text = asset_notices[0] if asset_notices else "local assets unavailable"
        asset_notice_html = (
            "<div class=\"heatmap-asset-note\">"
            "Asset mode: CDN fallback. "
            f"{html.escape(reason_text)}"
            "</div>"
        )
    asset_meta_json = json.dumps(
        {
            "mode": asset_mode,
            "notices": asset_notices,
            "asset_urls": asset_plan.get("asset_urls", {}),
        },
        ensure_ascii=True,
    )
    report_summary = cast(Dict[str, Any], analysis_payload.get("report_summary") or {})
    report_summary["asset_mode"] = asset_mode
    analysis_payload["report_summary"] = report_summary
    asset_meta_text = asset_meta_json.replace("</", "<\\/")

    base_row_nodes = cast(List[Dict[str, Any]], base_viz_json.get("row_nodes", []))
    base_col_nodes = cast(List[Dict[str, Any]], base_viz_json.get("col_nodes", []))
    base_row_names = [str(node.get("name", "")) for node in base_row_nodes]
    base_col_names = [str(node.get("name", "")) for node in base_col_nodes]
    row_label_lookup = set(row_labels)
    col_label_lookup = set(col_labels)
    ordered_row_labels = _ordered_clustergrammer_axis_labels(
        base_row_names, row_label_lookup, row_labels
    )
    ordered_col_labels = _ordered_clustergrammer_axis_labels(
        base_col_names, col_label_lookup, col_labels
    )
    row_node_by_name = {str(node.get("name", "")): node for node in base_row_nodes}
    col_node_by_name = {str(node.get("name", "")): node for node in base_col_nodes}
    ordered_row_nodes = [row_node_by_name.get(label, {"name": label}) for label in ordered_row_labels]
    ordered_col_nodes = [col_node_by_name.get(label, {"name": label}) for label in ordered_col_labels]
    preserve_views = (
        len(ordered_row_nodes) == len(base_row_nodes)
        and len(ordered_col_nodes) == len(base_col_nodes)
    )
    derived_views = (
        cast(List[Any], copy.deepcopy(base_viz_json.get("views", [])))
        if preserve_views
        else []
    )
    fixed_viz_json = _clone_clustergrammer_viz_for_mode(
        base_viz_json,
        _build_clustergrammer_mode_mat(
            matrix,
            row_labels,
            col_labels,
            ordered_row_labels,
            ordered_col_labels,
            fixed_breaks,
        ),
        row_nodes=ordered_row_nodes,
        col_nodes=ordered_col_nodes,
        views=derived_views,
    )
    quantile_viz_json = None
    if quantile_breaks is not None:
        quantile_viz_json = _clone_clustergrammer_viz_for_mode(
            base_viz_json,
            _build_clustergrammer_mode_mat(
                matrix,
                row_labels,
                col_labels,
                ordered_row_labels,
                ordered_col_labels,
                quantile_breaks,
            ),
            row_nodes=ordered_row_nodes,
            col_nodes=ordered_col_nodes,
            views=derived_views,
        )
    viz_payload_json_text = json.dumps(
        {"fixed": fixed_viz_json, "quantile": quantile_viz_json},
        ensure_ascii=True,
    ).replace("</", "<\\/")
    meta_json_text = json.dumps(cell_meta, ensure_ascii=True).replace("</", "<\\/")
    rowmeta_json_text = json.dumps(row_label_meta, ensure_ascii=True).replace("</", "<\\/")
    colmeta_json_text = json.dumps(col_label_meta, ensure_ascii=True).replace("</", "<\\/")
    analysis_json_text = json.dumps(analysis_payload, ensure_ascii=True).replace("</", "<\\/")
    scale_json_text = json.dumps(scale_payload, ensure_ascii=True).replace("</", "<\\/")
    mat_rows = len(matrix)
    mat_cols = len(matrix[0]) if matrix else 0
    elapsed_s = time.perf_counter() - build_start
    viz_bytes = len(viz_payload_json_text.encode("utf-8"))
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

    quantile_available = scale_payload.get("quantile_breaks") is not None
    mode_text = _scale_mode_text(scale_mode)
    breaks_line = _format_breaks_line(chosen_breaks) if chosen_breaks else ""
    fixed_line = _format_breaks_line(fixed_breaks) if fixed_breaks else ""
    return render_clustergrammer_client_html(
        assets_block=assets_block,
        safe_run_id=safe_run_id,
        asset_notice_html=asset_notice_html,
        viz_payload_json_text=viz_payload_json_text,
        meta_json_text=meta_json_text,
        rowmeta_json_text=rowmeta_json_text,
        colmeta_json_text=colmeta_json_text,
        analysis_json_text=analysis_json_text,
        scale_json_text=scale_json_text,
        asset_meta_text=asset_meta_text,
        ligand_page_prefix=ligand_page_prefix or "",
        color_scale_range=color_scale_range,
        use_dendrogram=use_dendrogram,
        mode_text=mode_text,
        breaks_line=breaks_line,
        fixed_line=fixed_line,
        fallback_reason=fallback_reason,
        scale_mode=scale_mode,
        requested_mode=requested_mode,
        quantile_available=quantile_available,
    )
