from __future__ import annotations

from dataclasses import dataclass

from analysis.reporting.heatmap_report_client_extras_core_module import (
    render_core_js_module,
)
from analysis.reporting.heatmap_report_client_extras_discovery_module import (
    render_discovery_js_module,
)
from analysis.reporting.heatmap_report_client_extras_side_effect_module import (
    render_side_effect_js_module,
)

@dataclass(frozen=True)
class HeatmapReportExtras:
    style_block: str
    pre_html_block: str
    start_here_block: str
    html_block: str
    script_block: str


def build_heatmap_report_extras(safe_run_id: str) -> HeatmapReportExtras:
    hero_id = f"cg-heatmap-hero-{safe_run_id}"
    hero_summary_id = f"cg-heatmap-hero-summary-{safe_run_id}"
    hero_search_id = f"cg-heatmap-hero-search-{safe_run_id}"
    hero_results_id = f"cg-heatmap-hero-results-{safe_run_id}"
    hero_state_id = f"cg-heatmap-hero-state-{safe_run_id}"
    hero_reset_id = f"cg-heatmap-hero-reset-{safe_run_id}"
    hero_copy_link_id = f"cg-heatmap-hero-copy-link-{safe_run_id}"
    validation_panel_id = f"cg-heatmap-validation-{safe_run_id}"
    workflow_panel_id = f"cg-heatmap-workflow-{safe_run_id}"
    discovery_panel_id = f"cg-heatmap-discovery-{safe_run_id}"
    discovery_summary_id = f"cg-heatmap-discovery-summary-{safe_run_id}"
    discovery_search_id = f"cg-heatmap-discovery-search-{safe_run_id}"
    discovery_rows_id = f"cg-heatmap-discovery-rows-{safe_run_id}"
    discovery_cols_id = f"cg-heatmap-discovery-cols-{safe_run_id}"
    discovery_table_id = f"cg-heatmap-discovery-table-{safe_run_id}"
    discovery_org_family_id = f"cg-heatmap-discovery-org-family-{safe_run_id}"
    discovery_org_pathway_id = f"cg-heatmap-discovery-org-pathway-{safe_run_id}"
    discovery_org_adme_id = f"cg-heatmap-discovery-org-adme-{safe_run_id}"
    discovery_org_safety_id = f"cg-heatmap-discovery-org-safety-{safe_run_id}"
    discovery_quick_targets_id = f"cg-heatmap-discovery-quick-targets-{safe_run_id}"
    discovery_quick_ligands_id = f"cg-heatmap-discovery-quick-ligands-{safe_run_id}"
    discovery_quick_safety_id = f"cg-heatmap-discovery-quick-safety-{safe_run_id}"
    discovery_quick_unresolved_id = f"cg-heatmap-discovery-quick-unresolved-{safe_run_id}"
    discovery_family_chips_id = f"cg-heatmap-discovery-family-chips-{safe_run_id}"
    discovery_pathway_chips_id = f"cg-heatmap-discovery-pathway-chips-{safe_run_id}"
    discovery_adme_chips_id = f"cg-heatmap-discovery-adme-chips-{safe_run_id}"
    discovery_safety_chips_id = f"cg-heatmap-discovery-safety-chips-{safe_run_id}"
    export_view_id = f"cg-heatmap-export-view-{safe_run_id}"
    export_rank_id = f"cg-heatmap-export-rank-{safe_run_id}"
    export_evidence_id = f"cg-heatmap-export-evidence-{safe_run_id}"
    save_state_id = f"cg-heatmap-save-state-{safe_run_id}"
    load_state_id = f"cg-heatmap-load-state-{safe_run_id}"
    load_state_file_id = f"cg-heatmap-load-state-file-{safe_run_id}"
    side_effect_overlap_button_id = f"cg-heatmap-side-effect-overlap-{safe_run_id}"
    side_effect_threshold_id = f"cg-heatmap-side-effect-threshold-{safe_run_id}"
    side_effect_overlap_summary_id = f"cg-heatmap-side-effect-summary-{safe_run_id}"
    side_effect_overlap_panel_id = f"cg-heatmap-side-effect-panel-{safe_run_id}"

    style_block = """
<style>
  .heatmap-report-hero { margin: 0 0 10px 0; border: 1px solid #d5dccc; border-radius: 14px; background: linear-gradient(180deg, #fcfdf9 0%, #f3f7ef 100%); box-shadow: 0 10px 24px rgba(22,34,28,0.08); overflow: hidden; }
  .heatmap-report-hero-head { padding: 12px 14px; border-bottom: 1px solid #dde5d5; display: grid; gap: 8px; }
  .heatmap-report-kicker { font-size: 0.76em; font-weight: 800; text-transform: uppercase; letter-spacing: 0.06em; color: #5e715d; }
  .heatmap-report-title { font-size: 1.28rem; font-weight: 800; color: #213629; }
  .heatmap-report-subtitle { color: #536651; line-height: 1.45; font-size: 0.94em; max-width: 90ch; }
  .heatmap-report-summary { display: flex; gap: 8px; flex-wrap: wrap; align-items: center; font-size: 0.83em; color: #435643; }
  .heatmap-report-pill { border: 1px solid #d6dfcf; border-radius: 999px; background: rgba(255,255,255,0.92); padding: 4px 9px; }
  .heatmap-report-hero-actions { padding: 12px 14px; display: grid; gap: 10px; }
  .heatmap-report-search-row { display: grid; grid-template-columns: minmax(0, 1fr) auto auto; gap: 8px; align-items: end; }
  .heatmap-report-search-row label { display: grid; gap: 4px; font-size: 0.79em; font-weight: 700; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; }
  .heatmap-report-search-row input { height: 38px; border: 1px solid #b8c5af; border-radius: 10px; padding: 0 12px; background: rgba(255,255,255,0.96); color: #243224; }
  .heatmap-report-search-row button { height: 38px; border: 1px solid #557055; border-radius: 10px; background: #fff; color: #1b3521; font-weight: 800; cursor: pointer; padding: 0 14px; }
  .heatmap-report-search-row button:hover { background: #f1f7ee; }
  .heatmap-start-here { border: 1px dashed #c8d2c1; border-radius: 10px; background: rgba(255,255,255,0.72); padding: 8px 10px; display: grid; gap: 4px; }
  .heatmap-start-here-title { font-size: 0.79em; font-weight: 800; text-transform: uppercase; letter-spacing: 0.03em; color: #334733; }
  .heatmap-start-here-copy { font-size: 0.84em; color: #536651; line-height: 1.45; }
  .heatmap-side-effect-tools { border: 1px solid #d8dfcf; border-radius: 8px; background: rgba(255,255,255,0.84); padding: 8px 10px; display: grid; gap: 7px; }
  .heatmap-side-effect-actions { display: flex; gap: 8px; flex-wrap: wrap; align-items: center; }
  .heatmap-side-effect-actions button { height: 34px; border: 1px solid #557055; border-radius: 8px; background: #fff; color: #1b3521; font-weight: 800; cursor: pointer; padding: 0 12px; }
  .heatmap-side-effect-actions button.is-active { background: #6f2f18; border-color: #6f2f18; color: #fff; box-shadow: 0 0 0 2px rgba(111,47,24,0.18); }
  .heatmap-side-effect-actions button:hover { background: #fff5ee; }
  .heatmap-side-effect-threshold { display: inline-flex; gap: 6px; align-items: center; font-size: 0.8em; font-weight: 700; color: #3f533d; }
  .heatmap-side-effect-threshold select { height: 32px; border: 1px solid #b8c5af; border-radius: 8px; background: #fff; color: #243224; padding: 0 8px; }
  .heatmap-side-effect-summary { font-size: 0.83em; color: #536651; line-height: 1.4; }
  .heatmap-side-effect-panel { display: none; border-top: 1px solid #e1e8dc; padding-top: 8px; }
  .heatmap-side-effect-panel.is-active { display: grid; gap: 8px; }
  .heatmap-side-effect-table-wrap { overflow: auto; max-height: 280px; border: 1px solid #dde5d5; border-radius: 8px; background: rgba(255,255,255,0.92); }
  .heatmap-side-effect-table { width: 100%; border-collapse: collapse; font-size: 0.81em; }
  .heatmap-side-effect-table th, .heatmap-side-effect-table td { padding: 6px 8px; border-bottom: 1px solid #edf1e8; text-align: left; vertical-align: top; }
  .heatmap-side-effect-table th { position: sticky; top: 0; background: #f7faf3; color: #314631; z-index: 1; }
  .heatmap-side-effect-table td.num { font-family: Consolas, monospace; text-align: right; white-space: nowrap; }
  .heatmap-side-effect-note { font-size: 0.8em; color: #5b6d58; line-height: 1.45; }
  .heatmap-active-state { display: flex; gap: 8px; flex-wrap: wrap; align-items: center; }
  .heatmap-active-state.is-empty { display: none; }
  .heatmap-active-chip { border: 1px solid #3e5f47; border-radius: 999px; background: rgba(239,248,240,0.98); padding: 4px 10px; display: inline-flex; gap: 8px; align-items: center; font-size: 0.82em; color: #183121; box-shadow: inset 0 0 0 1px rgba(62,95,71,0.12); }
  .heatmap-active-chip strong { color: #0f2518; }
  .heatmap-active-chip button { border: 0; background: transparent; color: #1b4f30; font-weight: 800; cursor: pointer; padding: 0; text-decoration: underline; }
  .heatmap-hero-results { border-top: 1px solid #e1e8dc; padding-top: 10px; display: none; }
  .heatmap-hero-results.is-active { display: grid; gap: 10px; }
  .heatmap-hero-results-empty { color: #5b6d58; font-size: 0.84em; }
  .heatmap-search-groups { display: grid; gap: 10px; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); }
  .heatmap-search-group { border: 1px solid #dae3d3; border-radius: 12px; background: rgba(255,255,255,0.88); overflow: hidden; }
  .heatmap-search-group-head { padding: 8px 10px; border-bottom: 1px solid #e6ece0; font-size: 0.78em; font-weight: 800; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; }
  .heatmap-search-group-body { padding: 8px; display: grid; gap: 6px; }
  .heatmap-search-result { border: 1px solid #dce5d6; border-radius: 10px; background: #fff; padding: 8px; text-align: left; cursor: pointer; display: grid; gap: 3px; }
  .heatmap-search-result-label { font-weight: 700; color: #203629; }
  .heatmap-search-result-meta { font-size: 0.82em; color: #5b6d58; line-height: 1.35; }
  .heatmap-report-grid { margin-top: 10px; display: grid; gap: 10px; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); }
  .heatmap-report-panel { border: 1px solid #d5dccc; border-radius: 12px; background: linear-gradient(180deg, #fbfcf8 0%, #f3f7ef 100%); box-shadow: 0 8px 20px rgba(22,34,28,0.06); overflow: hidden; }
  .heatmap-report-panel-head { padding: 8px 10px; border-bottom: 1px solid #dde5d5; display: grid; gap: 2px; }
  .heatmap-report-panel-title { font-weight: 700; color: #243a2a; }
  .heatmap-report-panel-subtitle { font-size: 0.82em; color: #556555; }
  .heatmap-report-panel-body { padding: 10px; display: grid; gap: 8px; }
  .heatmap-validation-list { display: grid; gap: 8px; }
  .heatmap-validation-item { border: 1px solid #d8dfcf; border-radius: 10px; background: rgba(255,255,255,0.88); padding: 8px 10px; display: grid; gap: 3px; }
  .heatmap-validation-label { font-size: 0.77em; font-weight: 800; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; }
  .heatmap-validation-value { font-weight: 700; color: #253a2a; }
  .heatmap-validation-note { font-size: 0.8em; color: #5b6d58; }
  .heatmap-workflow-actions { display: flex; gap: 8px; flex-wrap: wrap; }
  .heatmap-workflow-actions button { height: 34px; border: 1px solid #557055; border-radius: 8px; background: #fff; color: #1b3521; font-weight: 700; cursor: pointer; padding: 0 12px; }
  .heatmap-workflow-actions button:hover { background: #f1f7ee; }
  .heatmap-workflow-note { font-size: 0.8em; color: #556555; }
  .heatmap-discovery-controls { display: grid; gap: 8px; border: 1px solid #d8dfcf; border-radius: 10px; background: rgba(255,255,255,0.86); padding: 8px; }
  .heatmap-discovery-toolbar { display: flex; gap: 8px; flex-wrap: wrap; align-items: center; }
  .heatmap-discovery-toolbar input { min-width: 220px; height: 34px; border: 1px solid #b8c5af; border-radius: 8px; padding: 0 10px; background: rgba(255,255,255,0.96); color: #243224; }
  .heatmap-discovery-group-label { font-size: 0.77em; font-weight: 800; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; }
  .heatmap-discovery-summary { font-size: 0.82em; color: #556555; }
  .heatmap-discovery-quick { display: flex; gap: 8px; flex-wrap: wrap; }
  .heatmap-discovery-quick button { height: 34px; border: 1px solid #557055; border-radius: 999px; background: #fff; color: #1b3521; font-weight: 700; cursor: pointer; padding: 0 12px; }
  .heatmap-discovery-quick button:hover { background: #f1f7ee; }
  .heatmap-discovery-chip-grid { display: grid; gap: 8px; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); }
  .heatmap-discovery-chip-panel { border: 1px solid #d8dfcf; border-radius: 10px; background: rgba(255,255,255,0.88); padding: 8px; display: grid; gap: 6px; }
  .heatmap-discovery-chip-title { font-size: 0.77em; font-weight: 800; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; }
  .heatmap-discovery-chips { display: flex; gap: 6px; flex-wrap: wrap; }
  .heatmap-discovery-chip { border: 1px solid #ccd6c3; border-radius: 999px; background: rgba(255,255,255,0.95); padding: 4px 8px; display: inline-flex; gap: 6px; align-items: center; color: #294331; font-size: 0.76em; cursor: pointer; }
  .heatmap-discovery-chip.is-active { border-color: #234d35; background: rgba(228,243,232,0.98); box-shadow: 0 0 0 2px rgba(35,77,53,0.2); }
  .heatmap-discovery-chip strong { font-weight: 700; }
  .heatmap-discovery-table-wrap { overflow: auto; max-height: 420px; border: 1px solid #dde5d5; border-radius: 10px; background: rgba(255,255,255,0.9); }
  .heatmap-discovery-table { width: 100%; border-collapse: collapse; font-size: 0.82em; }
  .heatmap-discovery-table th, .heatmap-discovery-table td { padding: 7px 8px; border-bottom: 1px solid #edf1e8; text-align: left; vertical-align: top; }
  .heatmap-discovery-table th { position: sticky; top: 0; background: #f7faf3; z-index: 1; color: #314631; }
  .heatmap-discovery-table td.num { font-family: Consolas, monospace; text-align: right; }
  .heatmap-discovery-table button { height: 28px; border: 1px solid #557055; border-radius: 8px; background: #fff; color: #1b3521; font-weight: 700; cursor: pointer; padding: 0 10px; }
  .heatmap-discovery-table button:hover { background: #f1f7ee; }
  .heatmap-report-search-row button:focus-visible,
  .heatmap-report-search-row input:focus-visible,
  .heatmap-active-chip button:focus-visible,
  .heatmap-side-effect-actions button:focus-visible,
  .heatmap-workflow-actions button:focus-visible,
  .heatmap-discovery-toolbar input:focus-visible,
  .heatmap-discovery-quick button:focus-visible,
  .heatmap-discovery-chip:focus-visible,
  .heatmap-discovery-table button:focus-visible {
    outline: 3px solid #0f4c81;
    outline-offset: 2px;
  }
  .cg-heatmap-viz .cg-discovery-highlight { fill: #a14c22 !important; font-weight: 800 !important; text-decoration: underline; }
  .cg-heatmap-viz .cg-side-effect-match { stroke: #6f2f18 !important; stroke-width: 2.2px !important; vector-effect: non-scaling-stroke; filter: drop-shadow(0 0 4px rgba(111,47,24,0.55)); }
  @media (max-width: 880px) {
    .heatmap-report-search-row { grid-template-columns: minmax(0, 1fr); }
  }
</style>
""".strip()

    pre_html_block = f"""
<section class="heatmap-report-hero" id="{hero_id}">
  <div class="heatmap-report-hero-head">
    <div class="heatmap-report-kicker">Interactive Discovery Report</div>
    <div class="heatmap-report-title">Docking heatmap explorer for scientist-first discovery</div>
    <div class="heatmap-report-subtitle">Search ligands, PDBs, proteins, UniProt accessions, pathways, ADME categories, safety buckets, and side-effect/adverse-event evidence. Use search to jump into the existing heatmap interactions instead of hunting manually through the matrix.</div>
    <div class="heatmap-report-summary" id="{hero_summary_id}"></div>
  </div>
  <div class="heatmap-report-hero-actions">
    <div class="heatmap-report-search-row">
      <label for="{hero_search_id}">Global search
        <input type="search" id="{hero_search_id}" placeholder="Search ligand, PDB, target, UniProt, family, pathway, ADME, safety, side effect">
      </label>
      <button type="button" id="{hero_reset_id}">Reset all</button>
      <button type="button" id="{hero_copy_link_id}">Copy shareable link</button>
    </div>
    <div class="heatmap-active-state" id="{hero_state_id}"></div>
    <div class="heatmap-hero-results" id="{hero_results_id}"></div>
  </div>
</section>
""".strip()
    start_here_block = f"""
<div class="heatmap-start-here">
  <div class="heatmap-start-here-title">Start Here</div>
  <div class="heatmap-start-here-copy">Search for a ligand or target first. Click a ligand to sort targets left-to-right, click a PDB to focus one target and rank ligands top-to-bottom, then use target organization to regroup by family, pathway, ADME, or safety.</div>
  <div class="heatmap-side-effect-tools">
    <div class="heatmap-side-effect-actions">
      <button type="button" id="{side_effect_overlap_button_id}">Highlight Side-Effect Matches</button>
      <label class="heatmap-side-effect-threshold" for="{side_effect_threshold_id}">Overlap threshold
        <select id="{side_effect_threshold_id}">
          <option value="0.01">1%</option>
          <option value="0.02">2%</option>
          <option value="0.05">5%</option>
          <option value="0.10">10%</option>
        </select>
      </label>
      <div class="heatmap-side-effect-summary" id="{side_effect_overlap_summary_id}"></div>
    </div>
    <div class="heatmap-side-effect-panel" id="{side_effect_overlap_panel_id}"></div>
  </div>
</div>
""".strip()

    html_block = f"""
<div class="heatmap-report-grid">
  <section class="heatmap-report-panel" id="{validation_panel_id}">
    <div class="heatmap-report-panel-head">
      <div class="heatmap-report-panel-title">Validation</div>
      <div class="heatmap-report-panel-subtitle">Coverage, cache freshness, and mapping health for this report.</div>
    </div>
    <div class="heatmap-report-panel-body"></div>
  </section>
  <section class="heatmap-report-panel" id="{workflow_panel_id}">
    <div class="heatmap-report-panel-head">
      <div class="heatmap-report-panel-title">Analyst Workflows</div>
      <div class="heatmap-report-panel-subtitle">Export the current view, rankings, evidence, or a reusable state snapshot.</div>
    </div>
    <div class="heatmap-report-panel-body">
      <div class="heatmap-workflow-actions">
        <button type="button" id="{export_view_id}">Export Current View CSV</button>
        <button type="button" id="{export_rank_id}">Export Ranking CSV</button>
        <button type="button" id="{export_evidence_id}">Export Evidence CSV</button>
        <button type="button" id="{save_state_id}">Save State JSON</button>
        <button type="button" id="{load_state_id}">Load State JSON</button>
      </div>
      <div class="heatmap-workflow-note">State is mirrored into the URL hash and can also be copied directly from the top header.</div>
      <input type="file" id="{load_state_file_id}" accept="application/json,.json" style="display:none">
    </div>
  </section>
</div>
<section class="heatmap-report-panel" id="{discovery_panel_id}">
  <div class="heatmap-report-panel-head">
    <div class="heatmap-report-panel-title">Discovery</div>
    <div class="heatmap-report-panel-subtitle">Search and rank the embedded ligand universe with discovery-specific controls.</div>
  </div>
  <div class="heatmap-report-panel-body">
    <div class="heatmap-discovery-controls">
      <div class="heatmap-discovery-group-label">Discovery table controls</div>
      <div class="heatmap-discovery-toolbar">
        <div class="heatmap-segmented">
          <button type="button" id="{discovery_rows_id}" class="is-active">Ligands</button>
        </div>
        <input type="search" id="{discovery_search_id}" placeholder="Search ligand, PDB, target, UniProt, family, ADME, safety">
      </div>
      <div class="heatmap-discovery-group-label">Discovery presets</div>
      <div class="heatmap-discovery-quick">
        <button type="button" id="{discovery_quick_ligands_id}">Most promiscuous ligands</button>
      </div>
    </div>
    <div class="heatmap-discovery-summary" id="{discovery_summary_id}"></div>
    <div class="heatmap-discovery-table-wrap">
      <table class="heatmap-discovery-table" id="{discovery_table_id}"></table>
    </div>
  </div>
</section>
""".strip()

    script_block = f"""
    var validationPanelEl = document.getElementById({validation_panel_id!r});
    var discoveryPanelEl = document.getElementById({discovery_panel_id!r});
    var heroSummaryEl = document.getElementById({hero_summary_id!r});
    var heroSearchEl = document.getElementById({hero_search_id!r});
    var heroResultsEl = document.getElementById({hero_results_id!r});
    var heroStateEl = document.getElementById({hero_state_id!r});
    var heroResetEl = document.getElementById({hero_reset_id!r});
    var heroCopyLinkEl = document.getElementById({hero_copy_link_id!r});
    var discoverySummaryEl = document.getElementById({discovery_summary_id!r});
    var discoverySearchEl = document.getElementById({discovery_search_id!r});
    var discoveryRowsEl = document.getElementById({discovery_rows_id!r});
    var discoveryColsEl = document.getElementById({discovery_cols_id!r});
    var discoveryTableEl = document.getElementById({discovery_table_id!r});
    var discoveryOrgFamilyEl = document.getElementById({discovery_org_family_id!r});
    var discoveryOrgPathwayEl = document.getElementById({discovery_org_pathway_id!r});
    var discoveryOrgAdmeEl = document.getElementById({discovery_org_adme_id!r});
    var discoveryOrgSafetyEl = document.getElementById({discovery_org_safety_id!r});
    var discoveryQuickTargetsEl = document.getElementById({discovery_quick_targets_id!r});
    var discoveryQuickLigandsEl = document.getElementById({discovery_quick_ligands_id!r});
    var discoveryQuickSafetyEl = document.getElementById({discovery_quick_safety_id!r});
    var discoveryQuickUnresolvedEl = document.getElementById({discovery_quick_unresolved_id!r});
    var discoveryFamilyChipsEl = document.getElementById({discovery_family_chips_id!r});
    var discoveryPathwayChipsEl = document.getElementById({discovery_pathway_chips_id!r});
    var discoveryAdmeChipsEl = document.getElementById({discovery_adme_chips_id!r});
    var discoverySafetyChipsEl = document.getElementById({discovery_safety_chips_id!r});
    var exportViewEl = document.getElementById({export_view_id!r});
    var exportRankEl = document.getElementById({export_rank_id!r});
    var exportEvidenceEl = document.getElementById({export_evidence_id!r});
    var saveStateEl = document.getElementById({save_state_id!r});
    var loadStateEl = document.getElementById({load_state_id!r});
    var loadStateFileEl = document.getElementById({load_state_file_id!r});
    var sideEffectOverlapButtonEl = document.getElementById({side_effect_overlap_button_id!r});
    var sideEffectThresholdEl = document.getElementById({side_effect_threshold_id!r});
    var sideEffectOverlapSummaryEl = document.getElementById({side_effect_overlap_summary_id!r});
    var sideEffectOverlapPanelEl = document.getElementById({side_effect_overlap_panel_id!r});
    var activeDiscoveryMode = "rows";
    var activeDiscoveryFilter = "all";
    var activeDiscoveryHighlightMode = "";
    var activeDiscoveryHighlightValue = "";
    var activeSideEffectOverlap = false;
    var sideEffectOverlapKeySet = Object.create(null);
    var suppressHashSync = false;
    {render_core_js_module(safe_run_id)}
    {render_side_effect_js_module()}
    {render_discovery_js_module()}
    function buildHeroSearchResults(query) {{
      if (!query) return {{ ligands: [], targets: [], annotations: [] }};
      var lowered = String(query || "").trim().toLowerCase();
      var ligands = allRowNames().map(function(name) {{
        var meta = lookupRowLabelMeta(name);
        return {{
          type: "ligand",
          name: name,
          label: name,
          meta: String(meta.chemotype_primary || meta.drug_effect_primary || meta.library || ""),
          searchable: [name, meta.chemotype_primary, meta.drug_effect_primary]
            .concat(meta.motif_tags || [], meta.side_effects || [], meta.safety_buckets || [])
            .join(" ").toLowerCase()
        }};
      }}).filter(function(entry) {{ return entry.searchable.indexOf(lowered) !== -1; }}).slice(0, 6);
      var targets = allColNames().map(function(name) {{
        var meta = lookupColLabelMeta(name);
        return {{
          type: "target",
          name: name,
          pdb: String(meta.pdb || ""),
          label: String(meta.pdb || targetShortName(name, meta)),
          meta: String(meta.target_name || ""),
          searchable: targetMetaSearchText(meta)
        }};
      }}).filter(function(entry) {{ return entry.searchable.indexOf(lowered) !== -1; }}).slice(0, 6);
      var annotations = annotationEntries().filter(function(entry) {{
        return entry.searchable.indexOf(lowered) !== -1;
      }}).sort(function(a, b) {{
        if (b.count !== a.count) return b.count - a.count;
        return a.label < b.label ? -1 : (a.label > b.label ? 1 : 0);
      }}).slice(0, 6);
      return {{ ligands: ligands, targets: targets, annotations: annotations }};
    }}
    function heroGroupHtml(title, rows, renderer) {{
      if (!rows.length) return "";
      return '<div class="heatmap-search-group"><div class="heatmap-search-group-head">' + escapeHtml(title) + '</div><div class="heatmap-search-group-body">' + rows.map(renderer).join("") + "</div></div>";
    }}
    function applyAnnotationSearchEntry(entry) {{
      if (!entry) return;
      activeDiscoveryHighlightMode = "";
      activeDiscoveryHighlightValue = "";
      if (discoverySearchEl) discoverySearchEl.value = String(entry.label || "");
      activeDiscoveryMode = "rows";
      activeDiscoveryFilter = "all";
      renderDiscoveryPanel();
      renderActiveStateBar();
      syncReportStateHash();
      scrollToNode(discoveryPanelEl);
    }}
    function renderHeroSearchResults() {{
      if (!heroResultsEl) return;
      var query = heroSearchEl ? String(heroSearchEl.value || "").trim() : "";
      var results = buildHeroSearchResults(query);
      var resultCount = results.ligands.length + results.targets.length + results.annotations.length;
      heroResultsEl.classList.toggle("is-active", !!query);
      if (!query) {{
        heroResultsEl.innerHTML = "";
        return;
      }}
      if (!resultCount) {{
        heroResultsEl.innerHTML = '<div class="heatmap-hero-results-empty">No matches found. Try a ligand, PDB, target name, UniProt, family, pathway, ADME category, safety bucket, or side-effect/adverse-event term.</div>';
        return;
      }}
      heroResultsEl.innerHTML = '<div class="heatmap-search-groups">'
        + heroGroupHtml("Ligands", results.ligands, function(entry) {{
          return '<button type="button" class="heatmap-search-result" data-search-ligand="' + escapeHtml(entry.name) + '"><span class="heatmap-search-result-label">' + escapeHtml(entry.label) + '</span><span class="heatmap-search-result-meta">' + escapeHtml(entry.meta || "Ligand result") + '</span></button>';
        }})
        + heroGroupHtml("Targets", results.targets, function(entry) {{
          return '<button type="button" class="heatmap-search-result" data-search-target="' + escapeHtml(entry.name) + '" data-search-pdb="' + escapeHtml(entry.pdb) + '"><span class="heatmap-search-result-label">' + escapeHtml(entry.label) + '</span><span class="heatmap-search-result-meta">' + escapeHtml(entry.meta || "Target result") + '</span></button>';
        }})
        + heroGroupHtml("Annotations", results.annotations, function(entry) {{
          return '<button type="button" class="heatmap-search-result" data-search-annotation="' + escapeHtml(entry.kind + "||" + entry.label) + '"><span class="heatmap-search-result-label">' + escapeHtml(entry.label) + '</span><span class="heatmap-search-result-meta">' + escapeHtml(entry.subtitle + " | " + entry.count + " targets") + '</span></button>';
        }})
        + '</div>';
      var ligandButtons = heroResultsEl.querySelectorAll("[data-search-ligand]");
      for (var i = 0; i < ligandButtons.length; i++) {{
        ligandButtons[i].addEventListener("click", function(evt) {{
          evt.preventDefault();
          evt.stopPropagation();
          activateLigandSort(String(this.getAttribute("data-search-ligand") || ""));
          renderCurrentMode();
        }});
      }}
      var targetButtons = heroResultsEl.querySelectorAll("[data-search-target]");
      for (var j = 0; j < targetButtons.length; j++) {{
        targetButtons[j].addEventListener("click", function(evt) {{
          evt.preventDefault();
          evt.stopPropagation();
          activatePdbSelection(
            String(this.getAttribute("data-search-target") || ""),
            String(this.getAttribute("data-search-pdb") || "")
          );
          renderCurrentMode();
        }});
      }}
      var annotationButtons = heroResultsEl.querySelectorAll("[data-search-annotation]");
      var annotations = annotationEntries();
      for (var k = 0; k < annotationButtons.length; k++) {{
        annotationButtons[k].addEventListener("click", function(evt) {{
          evt.preventDefault();
          evt.stopPropagation();
          var key = String(this.getAttribute("data-search-annotation") || "");
          var match = annotations.find(function(item) {{ return (item.kind + "||" + item.label) === key; }});
          applyAnnotationSearchEntry(match || null);
        }});
      }}
    }}
    function targetMatchesDiscoveryHighlight(meta) {{
      var query = String(activeDiscoveryHighlightValue || "").trim().toLowerCase();
      if (!query || !meta) return false;
      if (activeDiscoveryHighlightMode === "family") return String(meta.protein_family || "").toLowerCase() === query;
      if (activeDiscoveryHighlightMode === "adme") return String(meta.adme_category || "").toLowerCase() === query;
      if (activeDiscoveryHighlightMode === "pathway") {{
        if (String(meta.primary_display_pathway || "").toLowerCase() === query) return true;
        return Array.isArray(meta.pathway_memberships) && meta.pathway_memberships.some(function(value) {{ return String(value || "").toLowerCase() === query; }});
      }}
      if (activeDiscoveryHighlightMode === "safety") {{
        if (String(meta.primary_display_safety || "").toLowerCase() === query) return true;
        if (typeof targetMetaSearchText === "function") return targetMetaSearchText(meta).indexOf(query) !== -1;
        return JSON.stringify(meta || {{}}).toLowerCase().indexOf(query) !== -1;
      }}
      return false;
    }}
    function applyDiscoveryHighlights() {{
      var chips = document.querySelectorAll(".heatmap-discovery-chip");
      for (var j = 0; j < chips.length; j++) {{
        var chip = chips[j];
        var mode = String(chip.getAttribute("data-org-mode") || "");
        var value = String(chip.getAttribute("data-org-value") || "");
        chip.classList.toggle("is-active", mode === activeDiscoveryHighlightMode && value.toLowerCase() === String(activeDiscoveryHighlightValue || "").toLowerCase());
      }}
    }}
    function buildCurrentViewCsv() {{
      var payload = activeVizData();
      var rowNames = rowNamesFromPayload(payload);
      var colNames = colNamesFromPayload(payload);
      var lines = [];
      lines.push("# heatmap_state," + serializeReportStateToHash(buildReportState()));
      lines.push(["ligand"].concat(colNames.map(function(name) {{
        var meta = lookupColLabelMeta(name);
        return meta.pdb || targetShortName(name, meta);
      }})).map(csvEscape).join(","));
      for (var r = 0; r < rowNames.length; r++) {{
        var row = [rowNames[r]];
        for (var c = 0; c < colNames.length; c++) {{
          var meta = lookupMeta(rowNames[r], colNames[c]);
          row.push(activeValueText(meta));
        }}
        lines.push(row.map(csvEscape).join(","));
      }}
      return lines.join("\\n");
    }}
    function buildRankingCsv() {{
      var rows = rankingEntries();
      var lines = [];
      lines.push(["entity_type", "label", "pdb", "target_name", "metric", "value", "threshold", "organization_mode"].map(csvEscape).join(","));
      for (var i = 0; i < rows.length; i++) {{
        var entry = rows[i];
        lines.push([
          activeRankMode === "cols" ? "target" : "ligand",
          activeRankMode === "cols" ? (entry.meta.pdb || entry.name) : entry.name,
          entry.meta && entry.meta.pdb ? entry.meta.pdb : "",
          entry.meta && entry.meta.target_name ? entry.meta.target_name : "",
          activeRankMetric,
          formatMetricValue(activeRankMetric, entry.value),
          formatMetric(activeThreshold, 6),
          targetOrganizationMode
        ].map(csvEscape).join(","));
      }}
      return lines.join("\\n");
    }}
    function buildEvidenceCsv() {{
      var cols = allColNames();
      var lines = [];
      lines.push(["pdb", "target_name", "variant_ph", "protein_family", "adme_category", "safety_primary", "safety_secondary", "safety_confidence", "safety_sources", "promiscuity_z", "significant_hits", "hit_fraction", "entropy", "side_effects"].map(csvEscape).join(","));
      for (var i = 0; i < cols.length; i++) {{
        var name = cols[i];
        var meta = lookupColLabelMeta(name);
        lines.push([
          meta.pdb || "",
          meta.target_name || "",
          meta.variant_ph || "",
          meta.protein_family || "",
          meta.adme_category || "",
          meta.primary_display_safety || "",
          Array.isArray(meta.secondary_safety_buckets) ? meta.secondary_safety_buckets.join("; ") : "",
          meta.safety_confidence || "",
          Array.isArray(meta.safety_sources) ? meta.safety_sources.join("; ") : "",
          meta.promiscuity_z || "",
          meta.significant_hits_text || "",
          meta.hit_fraction_pct || meta.hit_fraction || "",
          meta.entropy || "",
          Array.isArray(meta.side_effects) ? meta.side_effects.join("; ") : ""
        ].map(csvEscape).join(","));
      }}
      return lines.join("\\n");
    }}
    function bindReportExtrasControls() {{
      if (container.getAttribute("data-cg-report-extras-bound") === "1") return;
      container.setAttribute("data-cg-report-extras-bound", "1");
      if (heroSearchEl) heroSearchEl.addEventListener("input", function() {{ renderHeroSearchResults(); syncReportStateHash(); }});
      if (heroResetEl) heroResetEl.addEventListener("click", function() {{
        if (heroSearchEl) heroSearchEl.value = "";
        if (discoverySearchEl) discoverySearchEl.value = "";
        activeDiscoveryFilter = "all";
        activeDiscoveryHighlightMode = "";
        activeDiscoveryHighlightValue = "";
        activeSideEffectOverlap = false;
        if (resetButtonEl) resetButtonEl.click();
        renderHeroSearchResults();
      }});
      if (heroCopyLinkEl) heroCopyLinkEl.addEventListener("click", function() {{ syncReportStateHash(); copyTextToClipboard(window.location.href); }});
      if (discoveryRowsEl) discoveryRowsEl.addEventListener("click", function() {{ activeDiscoveryMode = "rows"; renderDiscoveryPanel(); syncReportStateHash(); }});
      if (discoveryQuickLigandsEl) discoveryQuickLigandsEl.addEventListener("click", function() {{ activeDiscoveryMode = "rows"; activeDiscoveryFilter = "all"; activeDiscoveryHighlightMode = ""; activeDiscoveryHighlightValue = ""; renderDiscoveryPanel(); renderActiveStateBar(); syncReportStateHash(); }});
      if (discoverySearchEl) discoverySearchEl.addEventListener("input", function() {{ renderDiscoveryPanel(); syncReportStateHash(); }});
      if (sideEffectOverlapButtonEl) sideEffectOverlapButtonEl.addEventListener("click", function(evt) {{
        evt.preventDefault();
        evt.stopPropagation();
        activeSideEffectOverlap = !activeSideEffectOverlap;
        var payload = baseVizData();
        renderSideEffectOverlapPanel(payload, derivedRowIndices(payload), derivedColIndices(payload));
        applySideEffectHighlights();
        renderActiveStateBar();
        syncReportStateHash();
      }});
      if (sideEffectThresholdEl) sideEffectThresholdEl.addEventListener("change", function() {{
        activeThreshold = Math.max(0.001, Math.min(0.1, Number(sideEffectThresholdEl.value) || activeThreshold));
        renderCurrentMode();
        syncReportStateHash();
      }});
      if (resetButtonEl) resetButtonEl.addEventListener("click", function() {{
        activeSideEffectOverlap = false;
        syncSideEffectThresholdControl();
        var payload = baseVizData();
        renderSideEffectOverlapPanel(payload, derivedRowIndices(payload), derivedColIndices(payload));
        applySideEffectHighlights();
      }});
      if (exportViewEl) exportViewEl.addEventListener("click", function() {{ downloadTextFile("heatmap_current_view.csv", buildCurrentViewCsv(), "text/csv;charset=utf-8"); }});
      if (exportRankEl) exportRankEl.addEventListener("click", function() {{ downloadTextFile("heatmap_rankings.csv", buildRankingCsv(), "text/csv;charset=utf-8"); }});
      if (exportEvidenceEl) exportEvidenceEl.addEventListener("click", function() {{ downloadTextFile("heatmap_evidence.csv", buildEvidenceCsv(), "text/csv;charset=utf-8"); }});
      if (saveStateEl) saveStateEl.addEventListener("click", function() {{ downloadTextFile("heatmap_state.json", JSON.stringify(buildReportState(), null, 2), "application/json;charset=utf-8"); }});
      if (loadStateEl && loadStateFileEl) loadStateEl.addEventListener("click", function() {{ loadStateFileEl.click(); }});
      if (loadStateFileEl) loadStateFileEl.addEventListener("change", function() {{
        var file = loadStateFileEl.files && loadStateFileEl.files[0];
        if (!file) return;
        var reader = new FileReader();
        reader.onload = function() {{
          try {{
            var parsed = JSON.parse(String(reader.result || "{{}}"));
            applyReportState(parsed);
            renderCurrentMode();
          }} catch (_err) {{
            window.alert("Could not load heatmap state JSON.");
          }} finally {{
            loadStateFileEl.value = "";
          }}
        }};
        reader.readAsText(file);
      }});
    }}
    function refreshReportExtras(payload, visibleRowIndices, visibleColIndices) {{
      renderHeroSummary();
      renderActiveStateBar();
      renderHeroSearchResults();
      renderValidationPanel(payload, visibleRowIndices || [], visibleColIndices || []);
      renderDiscoveryPanel();
      renderSideEffectOverlapPanel(payload, visibleRowIndices || [], visibleColIndices || []);
      applyDiscoveryHighlights();
      if (activeSideEffectOverlap) {{
        setTimeout(applySideEffectHighlights, 0);
      }}
      syncReportStateHash();
    }}
""".strip()

    return HeatmapReportExtras(
        style_block=style_block,
        pre_html_block=pre_html_block,
        start_here_block=start_here_block,
        html_block=html_block,
        script_block=script_block,
    )
