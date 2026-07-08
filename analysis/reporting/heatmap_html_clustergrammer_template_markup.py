from __future__ import annotations

import html


def build_clustergrammer_container_markup(
    *,
    report_extras,
    ligand_extras,
    container_id: str,
    scale_mode_label_id: str,
    mode_text: str,
    scale_breaks_label_id: str,
    breaks_line: str,
    view_summary_id: str,
    toggle_line: str,
    threshold_range_id: str,
    threshold_value_id: str,
    threshold_pct_id: str,
    threshold_fdr_id: str,
    threshold_scorch_fdr_id: str,
    threshold_scorch_relaxed_fdr_id: str,
    org_none_id: str,
    org_safety_id: str,
    org_adme_id: str,
    org_family_id: str,
    org_pathway_id: str,
    reset_button_id: str,
    safety_panel_id: str,
    adme_panel_id: str,
    safety_search_input_id: str,
    safety_search_clear_id: str,
    explanation_id: str,
    fixed_reference_line: str,
    fallback_line: str,
    asset_notice_html: str,
    org_strip_id: str,
    legend_line: str,
    host_id: str,
    panel_id: str,
    tab_scatter_id: str,
    tab_rank_id: str,
    scatter_mode_rows_id: str,
    rank_metric_prom_id: str,
    rank_metric_hits_id: str,
    rank_metric_frac_id: str,
    viz_json_id: str,
    meta_json_id: str,
    rowmeta_json_id: str,
    colmeta_json_id: str,
    analysis_json_id: str,
    scale_json_id: str,
    safe_run_id: str,
    asset_meta_text: str,
    viz_payload_json_text: str,
    meta_json_text: str,
    rowmeta_json_text: str,
    colmeta_json_text: str,
    analysis_json_text: str,
    scale_json_text: str,
) -> str:
    return f"""
{report_extras.pre_html_block}
<div class="heatmap-scale-meta heatmap-controls-card" data-heatmap-controls="card">
  <div class="heatmap-scale-mode" id="{scale_mode_label_id}">{html.escape(mode_text)}</div>
  <div class="heatmap-scale-breaks" id="{scale_breaks_label_id}">Breaks used: {html.escape(breaks_line)}</div>
  <div class="heatmap-view-summary" id="{view_summary_id}"></div>
  {toggle_line}
  <div class="heatmap-analysis-controls">
    <div class="heatmap-control">
      <label for="{threshold_range_id}">Significance Threshold</label>
      <input type="range" id="{threshold_range_id}" min="0.001" max="0.100" step="0.001">
      <div class="heatmap-segmented">
        <button type="button" id="{threshold_pct_id}" class="is-active">Percentile</button>
        <button type="button" id="{threshold_fdr_id}">FDR q</button>
        <button type="button" id="{threshold_scorch_fdr_id}">SCORCH q +1</button>
        <button type="button" id="{threshold_scorch_relaxed_fdr_id}">SCORCH q relaxed</button>
      </div>
      <div class="heatmap-control-value" id="{threshold_value_id}"></div>
    </div>
    <div class="heatmap-control">
      <label>Target Organization</label>
      <div class="heatmap-segmented">
        <button type="button" id="{org_none_id}" class="is-active">None</button>
        <button type="button" id="{org_safety_id}">Safety</button>
        <button type="button" id="{org_adme_id}">ADME</button>
        <button type="button" id="{org_family_id}">Protein family</button>
        <button type="button" id="{org_pathway_id}">Pathway</button>
      </div>
    </div>
    <div class="heatmap-control">
      <label>Display Defaults</label>
      <button type="button" id="{reset_button_id}">Reset Controls</button>
    </div>
  </div>
  <div class="heatmap-org-panels">
    <div class="heatmap-org-panel" id="{safety_panel_id}"></div>
    <div class="heatmap-org-panel" id="{adme_panel_id}"></div>
  </div>
  <div class="heatmap-safety-search">
    <label for="{safety_search_input_id}">Search adverse events
      <input type="search" id="{safety_search_input_id}" placeholder="e.g. hepatic, QT, rash">
    </label>
    <button type="button" id="{safety_search_clear_id}">Clear search</button>
  </div>
  <div class="heatmap-explanation" id="{explanation_id}">
    <div><code>promiscuity_z = 0.67449 * (count - median) / MAD</code> Higher means broader than peer median.</div>
    <div><code>entropy = -sum(p_i ln p_i) / ln(k)</code> Higher means hits are more evenly spread.</div>
    <div><code>significant hit = selected statistic &lt;= threshold</code> Percentile is broad screening; FDR q is the full target-level filter; SCORCH q +1 is the finite-sample corrected conditional filter; SCORCH q relaxed is the uncorrected conditional filter.</div>
  </div>
  {fixed_reference_line}
  {fallback_line}
  {asset_notice_html}
</div>
{report_extras.start_here_block}
<div class="heatmap-org-strip" id="{org_strip_id}"></div>
<div class="cg-heatmap-layout" id="{container_id}-layout">
  <aside class="cg-heatmap-sidebar">
    {ligand_extras.html_block}
  </aside>
  <div class="cg-heatmap-main">
    <div class="cg-heatmap-wrap" id="{container_id}-wrap">
      {legend_line}
      <div class="cg-heatmap-viz" id="{host_id}">
        <div id="{container_id}" style="width:100%;height:100%;">
          <div class="wait_message">Loading heatmap...</div>
        </div>
      </div>
    </div>
  </div>
</div>
<div class="heatmap-analysis-panel" id="{panel_id}">
  <div class="heatmap-analysis-toolbar">
    <div class="heatmap-analysis-tabs">
      <button type="button" id="{tab_scatter_id}" data-analysis-tab="scatter" class="is-active">Scatter</button>
      <button type="button" id="{tab_rank_id}" data-analysis-tab="rank">Ranking</button>
    </div>
  </div>
  <div class="heatmap-analysis-view is-active" data-analysis-view="scatter">
    <div class="heatmap-analysis-head">
      <div class="heatmap-analysis-title" data-analysis-title="scatter"></div>
      <div class="heatmap-analysis-subtitle" data-analysis-subtitle="scatter"></div>
    </div>
  <div class="heatmap-analysis-subtabs" data-analysis-controls="scatter">
    <button type="button" id="{scatter_mode_rows_id}" data-scatter-mode="rows" class="is-active">Ligands</button>
  </div>
    <svg class="heatmap-analysis-canvas" data-analysis-canvas="scatter"></svg>
    <div class="heatmap-analysis-note" data-analysis-note="scatter"></div>
  </div>
  <div class="heatmap-analysis-view" data-analysis-view="rank">
    <div class="heatmap-analysis-head">
      <div class="heatmap-analysis-title" data-analysis-title="rank"></div>
      <div class="heatmap-analysis-subtitle" data-analysis-subtitle="rank"></div>
    </div>
  <div class="heatmap-analysis-subtabs" data-analysis-controls="rank">
    <button type="button" id="{rank_metric_prom_id}" data-rank-metric="promiscuity_z" class="is-active">Promiscuity z-score</button>
    <button type="button" id="{rank_metric_hits_id}" data-rank-metric="significant_hits">Hit count</button>
    <button type="button" id="{rank_metric_frac_id}" data-rank-metric="hit_fraction">Hit fraction</button>
    </div>
    <svg class="heatmap-analysis-canvas" data-analysis-canvas="rank"></svg>
    <div class="heatmap-analysis-note" data-analysis-note="rank"></div>
  </div>
</div>
{report_extras.html_block}
<details class="heatmap-debug heatmap-debug-panel" id="heatmap-debug">
  <summary>Heatmap debug</summary>
  <pre id="heatmap-debug-log"></pre>
</details>
<script type="application/json" id="{viz_json_id}">{viz_payload_json_text}</script>
<script type="application/json" id="{meta_json_id}">{meta_json_text}</script>
<script type="application/json" id="{rowmeta_json_id}">{rowmeta_json_text}</script>
<script type="application/json" id="{colmeta_json_id}">{colmeta_json_text}</script>
<script type="application/json" id="{analysis_json_id}">{analysis_json_text}</script>
<script type="application/json" id="{scale_json_id}">{scale_json_text}</script>
<script type="application/json" id="cg-heatmap-assets-{safe_run_id}">{asset_meta_text}</script>
""".strip()
