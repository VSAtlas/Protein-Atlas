from __future__ import annotations

import html
import json
from typing import Sequence

from analysis.reporting.heatmap_html_clustergrammer_template_markup import (
    build_clustergrammer_container_markup,
)
from analysis.reporting.heatmap_html_clustergrammer_template_script import (
    build_clustergrammer_client_script,
)
from analysis.reporting.heatmap_html_clustergrammer_template_styles import (
    CLUSTERGRAMMER_HEATMAP_STYLE_BLOCK,
)
from analysis.reporting.heatmap_ligand_client_extras import (
    build_heatmap_ligand_organization_extras,
)
from analysis.reporting.heatmap_report_client_extras import (
    build_heatmap_report_extras,
)


def render_clustergrammer_client_html(
    *,
    assets_block: str,
    safe_run_id: str,
    asset_notice_html: str,
    viz_payload_json_text: str,
    meta_json_text: str,
    rowmeta_json_text: str,
    colmeta_json_text: str,
    analysis_json_text: str,
    scale_json_text: str,
    asset_meta_text: str,
    ligand_page_prefix: str,
    color_scale_range: Sequence[str],
    use_dendrogram: bool,
    mode_text: str,
    breaks_line: str,
    fixed_line: str,
    fallback_reason: str,
    scale_mode: str,
    requested_mode: str,
    quantile_available: bool,
) -> str:
    container_id = f"cg-heatmap-{safe_run_id}"
    viz_json_id = f"cg-heatmap-data-{safe_run_id}"
    meta_json_id = f"cg-heatmap-meta-{safe_run_id}"
    rowmeta_json_id = f"cg-heatmap-rowmeta-{safe_run_id}"
    colmeta_json_id = f"cg-heatmap-colmeta-{safe_run_id}"
    analysis_json_id = f"cg-heatmap-analysis-{safe_run_id}"
    scale_json_id = f"cg-heatmap-scale-{safe_run_id}"
    scale_mode_label_id = f"cg-heatmap-scale-mode-{safe_run_id}"
    scale_breaks_label_id = f"cg-heatmap-scale-breaks-{safe_run_id}"
    scale_fixed_label_id = f"cg-heatmap-scale-fixed-{safe_run_id}"
    scale_legend_bar_id = f"cg-heatmap-scale-legend-bar-{safe_run_id}"
    scale_legend_min_id = f"cg-heatmap-scale-legend-min-{safe_run_id}"
    scale_legend_mid_id = f"cg-heatmap-scale-legend-mid-{safe_run_id}"
    scale_legend_mid2_id = f"cg-heatmap-scale-legend-mid2-{safe_run_id}"
    scale_legend_max_id = f"cg-heatmap-scale-legend-max-{safe_run_id}"
    scale_toggle_fixed_id = f"cg-heatmap-scale-toggle-fixed-{safe_run_id}"
    scale_toggle_quantile_id = f"cg-heatmap-scale-toggle-quantile-{safe_run_id}"
    view_summary_id = f"cg-heatmap-view-summary-{safe_run_id}"
    explanation_id = f"cg-heatmap-explanation-{safe_run_id}"
    threshold_range_id = f"cg-heatmap-threshold-range-{safe_run_id}"
    threshold_value_id = f"cg-heatmap-threshold-value-{safe_run_id}"
    threshold_pct_id = f"cg-heatmap-threshold-pct-{safe_run_id}"
    threshold_fdr_id = f"cg-heatmap-threshold-fdr-{safe_run_id}"
    threshold_scorch_fdr_id = f"cg-heatmap-threshold-scorch-fdr-{safe_run_id}"
    threshold_scorch_relaxed_fdr_id = (
        f"cg-heatmap-threshold-scorch-relaxed-fdr-{safe_run_id}"
    )
    reset_button_id = f"cg-heatmap-reset-{safe_run_id}"
    org_none_id = f"cg-heatmap-org-none-{safe_run_id}"
    org_safety_id = f"cg-heatmap-org-safety-{safe_run_id}"
    org_adme_id = f"cg-heatmap-org-adme-{safe_run_id}"
    org_family_id = f"cg-heatmap-org-family-{safe_run_id}"
    org_pathway_id = f"cg-heatmap-org-pathway-{safe_run_id}"
    safety_panel_id = f"cg-heatmap-safety-panel-{safe_run_id}"
    adme_panel_id = f"cg-heatmap-adme-panel-{safe_run_id}"
    safety_search_input_id = f"cg-heatmap-safety-search-{safe_run_id}"
    safety_search_clear_id = f"cg-heatmap-safety-search-clear-{safe_run_id}"
    org_strip_id = f"cg-heatmap-org-strip-{safe_run_id}"
    panel_id = f"cg-heatmap-panel-{safe_run_id}"
    tab_scatter_id = f"cg-heatmap-tab-scatter-{safe_run_id}"
    tab_rank_id = f"cg-heatmap-tab-rank-{safe_run_id}"
    scatter_mode_rows_id = f"cg-heatmap-scatter-rows-{safe_run_id}"
    rank_metric_prom_id = f"cg-heatmap-rank-prom-{safe_run_id}"
    rank_metric_hits_id = f"cg-heatmap-rank-hits-{safe_run_id}"
    rank_metric_frac_id = f"cg-heatmap-rank-frac-{safe_run_id}"
    host_id = f"cg-heatmap-host-{safe_run_id}"
    report_extras = build_heatmap_report_extras(safe_run_id)
    ligand_extras = build_heatmap_ligand_organization_extras(safe_run_id)

    quantile_checked = requested_mode == "quantile" and quantile_available
    fixed_checked = not quantile_checked
    fixed_reference_line = (
        f"<div class=\"heatmap-scale-fixed\" id=\"{scale_fixed_label_id}\">"
        f"Fixed reference breaks: {html.escape(fixed_line)}"
        "</div>"
        if scale_mode == "quantile" and fixed_line
        else f"<div class=\"heatmap-scale-fixed\" id=\"{scale_fixed_label_id}\"></div>"
    )
    fallback_line = (
        f"<div class=\"heatmap-scale-fixed\">Quantile fallback: {html.escape(fallback_reason)}</div>"
        if requested_mode == "quantile"
        and scale_mode != "quantile"
        and fallback_reason
        else ""
    )
    toggle_line = (
        "<div class=\"heatmap-scale-toggle\">"
        "<span>Scale display:</span>"
        f"<label><input type=\"radio\" name=\"cg-scale-mode-{safe_run_id}\" id=\"{scale_toggle_fixed_id}\" value=\"fixed\"{' checked' if fixed_checked else ''}> Fixed (absolute)</label>"
        f"<label><input type=\"radio\" name=\"cg-scale-mode-{safe_run_id}\" id=\"{scale_toggle_quantile_id}\" value=\"quantile\"{' checked' if quantile_checked else ''}{' disabled' if not quantile_available else ''}> Quantile (within-run){' (unavailable)' if not quantile_available else ''}</label>"
        "</div>"
    )
    legend_line = (
        "<div class=\"cg-scale-legend-rail\">"
        f"<div class=\"cg-scale-legend-bar\" id=\"{scale_legend_bar_id}\"></div>"
        "<div class=\"cg-scale-legend-labels\">"
        f"<span id=\"{scale_legend_max_id}\"></span>"
        f"<span id=\"{scale_legend_mid2_id}\"></span>"
        f"<span id=\"{scale_legend_mid_id}\"></span>"
        f"<span id=\"{scale_legend_min_id}\"></span>"
        "</div>"
        "</div>"
    )

    style_block = "\n".join(
        [
            CLUSTERGRAMMER_HEATMAP_STYLE_BLOCK,
            ligand_extras.style_block,
            report_extras.style_block,
        ]
    )
    container_block = build_clustergrammer_container_markup(
        report_extras=report_extras,
        ligand_extras=ligand_extras,
        container_id=container_id,
        scale_mode_label_id=scale_mode_label_id,
        mode_text=mode_text,
        scale_breaks_label_id=scale_breaks_label_id,
        breaks_line=breaks_line,
        view_summary_id=view_summary_id,
        toggle_line=toggle_line,
        threshold_range_id=threshold_range_id,
        threshold_value_id=threshold_value_id,
        threshold_pct_id=threshold_pct_id,
        threshold_fdr_id=threshold_fdr_id,
        threshold_scorch_fdr_id=threshold_scorch_fdr_id,
        threshold_scorch_relaxed_fdr_id=threshold_scorch_relaxed_fdr_id,
        org_none_id=org_none_id,
        org_safety_id=org_safety_id,
        org_adme_id=org_adme_id,
        org_family_id=org_family_id,
        org_pathway_id=org_pathway_id,
        reset_button_id=reset_button_id,
        safety_panel_id=safety_panel_id,
        adme_panel_id=adme_panel_id,
        safety_search_input_id=safety_search_input_id,
        safety_search_clear_id=safety_search_clear_id,
        explanation_id=explanation_id,
        fixed_reference_line=fixed_reference_line,
        fallback_line=fallback_line,
        asset_notice_html=asset_notice_html,
        org_strip_id=org_strip_id,
        legend_line=legend_line,
        host_id=host_id,
        panel_id=panel_id,
        tab_scatter_id=tab_scatter_id,
        tab_rank_id=tab_rank_id,
        scatter_mode_rows_id=scatter_mode_rows_id,
        rank_metric_prom_id=rank_metric_prom_id,
        rank_metric_hits_id=rank_metric_hits_id,
        rank_metric_frac_id=rank_metric_frac_id,
        viz_json_id=viz_json_id,
        meta_json_id=meta_json_id,
        rowmeta_json_id=rowmeta_json_id,
        colmeta_json_id=colmeta_json_id,
        analysis_json_id=analysis_json_id,
        scale_json_id=scale_json_id,
        safe_run_id=safe_run_id,
        asset_meta_text=asset_meta_text,
        viz_payload_json_text=viz_payload_json_text,
        meta_json_text=meta_json_text,
        rowmeta_json_text=rowmeta_json_text,
        colmeta_json_text=colmeta_json_text,
        analysis_json_text=analysis_json_text,
        scale_json_text=scale_json_text,
    )

    ligand_prefix_json = json.dumps(ligand_page_prefix or "", ensure_ascii=True)
    color_range_json = json.dumps(list(color_scale_range), ensure_ascii=True)
    row_order_js = (
        "    args.row_order = \"rank\";\n    args.col_order = \"rank\";\n"
        if not use_dendrogram
        else ""
    )
    script_block = build_clustergrammer_client_script(
        container_id=container_id,
        host_id=host_id,
        viz_json_id=viz_json_id,
        meta_json_id=meta_json_id,
        rowmeta_json_id=rowmeta_json_id,
        colmeta_json_id=colmeta_json_id,
        analysis_json_id=analysis_json_id,
        scale_json_id=scale_json_id,
        scale_mode_label_id=scale_mode_label_id,
        scale_breaks_label_id=scale_breaks_label_id,
        scale_fixed_label_id=scale_fixed_label_id,
        view_summary_id=view_summary_id,
        threshold_range_id=threshold_range_id,
        threshold_value_id=threshold_value_id,
        threshold_pct_id=threshold_pct_id,
        threshold_fdr_id=threshold_fdr_id,
        threshold_scorch_fdr_id=threshold_scorch_fdr_id,
        threshold_scorch_relaxed_fdr_id=threshold_scorch_relaxed_fdr_id,
        reset_button_id=reset_button_id,
        org_none_id=org_none_id,
        org_safety_id=org_safety_id,
        org_adme_id=org_adme_id,
        org_family_id=org_family_id,
        org_pathway_id=org_pathway_id,
        safety_panel_id=safety_panel_id,
        adme_panel_id=adme_panel_id,
        safety_search_input_id=safety_search_input_id,
        safety_search_clear_id=safety_search_clear_id,
        org_strip_id=org_strip_id,
        scale_legend_bar_id=scale_legend_bar_id,
        scale_legend_min_id=scale_legend_min_id,
        scale_legend_mid_id=scale_legend_mid_id,
        scale_legend_mid2_id=scale_legend_mid2_id,
        scale_legend_max_id=scale_legend_max_id,
        scale_toggle_fixed_id=scale_toggle_fixed_id,
        scale_toggle_quantile_id=scale_toggle_quantile_id,
        safe_run_id=safe_run_id,
        ligand_prefix_json=ligand_prefix_json,
        color_range_json=color_range_json,
        row_order_js=row_order_js,
        report_extras=report_extras,
        ligand_extras=ligand_extras,
        use_dendrogram=use_dendrogram,
    )

    return "\n".join([assets_block, style_block, container_block, script_block])
