from __future__ import annotations


CLUSTERGRAMMER_HEATMAP_STYLE_BLOCK = """
<style>
  .heatmap-scale-meta { margin: 0 0 8px 0; border: 1px solid #d5dccc; background: linear-gradient(180deg, #fffef6 0%, #f5f8ef 100%); padding: 6px 8px; border-radius: 10px; box-shadow: 0 4px 12px rgba(32,42,28,0.08); }
  .heatmap-controls-card { margin-bottom: 10px; }
  .heatmap-scale-mode { font-weight: 700; color: #223222; }
  .heatmap-scale-breaks { margin-top: 2px; font-family: Consolas, monospace; font-size: 0.84em; color: #435043; }
  .heatmap-view-summary { margin-top: 4px; font-size: 0.82em; color: #36513b; }
  .heatmap-scale-toggle { margin-top: 4px; font-size: 0.85em; display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }
  .heatmap-scale-toggle label { display: inline-flex; align-items: center; gap: 4px; }
  .heatmap-analysis-controls { margin-top: 8px; display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 8px 14px; align-items: end; }
  .heatmap-control { display: grid; gap: 4px; }
  .heatmap-control label { font-size: 0.78em; font-weight: 700; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; }
  .heatmap-control .heatmap-control-value { font-family: Consolas, monospace; font-size: 0.84em; color: #435043; }
  .heatmap-control input[type="range"] { width: 100%; accent-color: #2c6b4d; }
  .heatmap-control button { height: 32px; border: 1px solid #4d674b; border-radius: 8px; background: #fbfcf7; color: #163020; font-weight: 800; cursor: pointer; }
  .heatmap-control button:hover { background: #eef6ea; }
  .heatmap-segmented { display: flex; gap: 6px; flex-wrap: wrap; }
  .heatmap-segmented button { height: 32px; padding: 0 10px; border: 1px solid #6a8368; border-radius: 999px; background: #fff; color: #1c3423; font-size: 0.82em; font-weight: 700; cursor: pointer; }
  .heatmap-segmented button.is-active { background: #0f4c3a; color: #fff; border-color: #0f4c3a; box-shadow: 0 0 0 2px rgba(15,76,58,0.22); }
  .heatmap-segmented button:hover { background: #eef5eb; }
  .heatmap-org-panels { margin-top: 8px; display: grid; gap: 8px; }
  .heatmap-org-panel { border: 1px solid #d7dfcf; border-radius: 10px; background: linear-gradient(180deg, #fcfdf9 0%, #f1f6ea 100%); padding: 7px 8px; display: grid; gap: 6px; }
  .heatmap-org-panel.is-active { border-color: #6d8c70; box-shadow: 0 0 0 2px rgba(60, 108, 74, 0.12); }
  .heatmap-org-panel-head { display: flex; justify-content: space-between; gap: 8px; align-items: baseline; }
  .heatmap-org-panel-title { font-size: 0.78em; font-weight: 800; color: #2f4933; text-transform: uppercase; letter-spacing: 0.03em; }
  .heatmap-org-panel-note { font-size: 0.74em; color: #5b6d58; }
  .heatmap-org-panel-chips { display: flex; gap: 6px; flex-wrap: wrap; }
  .heatmap-org-chip { border: 1px solid #ccd6c3; border-radius: 999px; background: rgba(255,255,255,0.9); padding: 4px 8px; display: inline-flex; gap: 6px; align-items: center; font-size: 0.76em; color: #2b4130; cursor: pointer; text-align: left; }
  button.heatmap-org-chip { appearance: none; -webkit-appearance: none; font: inherit; line-height: inherit; }
  .heatmap-org-chip.is-active { border-color: #234d35 !important; box-shadow: 0 0 0 2px rgba(35,77,53,0.2), inset 0 0 0 1px rgba(35,77,53,0.3) !important; background: rgba(228,243,232,0.98); }
  .heatmap-org-chip.is-search-hit { box-shadow: 0 0 0 2px rgba(154, 68, 31, 0.16); border-color: #9a441f; background: rgba(255, 248, 240, 0.98); }
  .heatmap-org-chip-name { font-weight: 700; }
  .heatmap-org-chip-count { font-family: Consolas, monospace; color: #516151; }
  .heatmap-safety-search { margin-top: 8px; display: grid; grid-template-columns: minmax(0, 1fr) auto; gap: 8px; align-items: end; }
  .heatmap-safety-search label { font-size: 0.78em; font-weight: 700; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; display: grid; gap: 4px; }
  .heatmap-safety-search input { height: 34px; border: 1px solid #b8c5af; border-radius: 8px; padding: 0 10px; background: rgba(255,255,255,0.96); color: #243224; }
  .heatmap-safety-search button { height: 34px; border: 1px solid #557055; border-radius: 8px; background: #fff; color: #1b3521; font-weight: 800; cursor: pointer; padding: 0 12px; }
  .heatmap-safety-search button:hover { background: #f1f7ee; }
  .heatmap-explanation { margin-top: 8px; padding-top: 8px; border-top: 1px dashed #c7d0bd; display: grid; gap: 3px; font-size: 0.82em; color: #334532; }
  .heatmap-explanation code { font-family: Consolas, monospace; color: #223222; background: rgba(255,255,255,0.65); padding: 1px 4px; border-radius: 4px; }
  .heatmap-asset-note { margin-top: 4px; font-size: 0.82em; color: #7c3e00; }
  .heatmap-org-strip { margin: 6px 0 8px 0; display: none; gap: 6px; align-items: stretch; width: 100%; }
  .heatmap-org-strip.is-active { display: flex; }
  .heatmap-org-block { min-width: 0; border: 1px solid #ccd7c4; border-radius: 9px; background: linear-gradient(180deg, #f9fbf5 0%, #eef4e8 100%); padding: 5px 8px; display: grid; gap: 1px; box-shadow: inset 0 1px 0 rgba(255,255,255,0.65); transition: background 120ms ease, border-color 120ms ease; }
  .heatmap-org-block-title { font-size: 0.76em; font-weight: 700; color: #28402c; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
  .heatmap-org-block-count { font-size: 0.72em; color: #5b6a58; }
  .cg-heatmap-layout { display: grid; grid-template-columns: minmax(270px, 320px) minmax(0, 1fr); gap: 10px; align-items: start; }
  .cg-heatmap-sidebar { border: 1px solid #d5dccc; border-radius: 12px; background: linear-gradient(180deg, #fbfcf8 0%, #f3f7ef 100%); box-shadow: 0 8px 20px rgba(22,34,28,0.06); padding: 8px 10px 10px; }
  .cg-heatmap-main { min-width: 0; }
  .cg-heatmap-wrap { display: flex; gap: 8px; align-items: stretch; flex-wrap: nowrap; width: 100%; min-width: 0; overflow-x: visible; }
  .cg-scale-legend-rail { flex: 0 0 68px; min-width: 68px; display: flex; align-items: flex-start; gap: 7px; padding-top: 2px; overflow: visible; }
  .cg-scale-legend-bar { flex: 0 0 12px; width: 12px; min-width: 12px; max-width: 12px; height: 180px; border: 1px solid #63745c; border-radius: 999px; box-sizing: border-box; }
  .cg-scale-legend-labels { flex: 1 1 auto; min-width: 0; height: 180px; display: flex; flex-direction: column; justify-content: space-between; font-family: Consolas, monospace; font-size: 0.73em; color: #435043; }
  .heatmap-scale-fixed { margin-top: 3px; font-family: Consolas, monospace; font-size: 0.8em; color: #435043; }
  .cg-heatmap-viz { flex: 1 1 0; width: 0; min-width: 0; min-height: clamp(560px, 72vh, 920px); border: 1px solid #d8dde3; border-radius: 12px; overflow: hidden; background: #fff; box-shadow: 0 10px 28px rgba(22,34,28,0.08); }
  .cg-heatmap-viz > div { width: 100% !important; height: 100% !important; }
  .cg-heatmap-viz .sidebar_wrapper { float: left !important; box-shadow: none !important; }
  .cg-heatmap-viz .viz_wrapper { float: left !important; }
  .cg-heatmap-viz .sidebar_text { margin-bottom: 2px !important; }
  .cg-heatmap-viz .slider_description { margin-top: 4px !important; }
  .cg-heatmap-viz .opacity_slider_container { display: none !important; }
  .cg-heatmap-viz .col_label_text,
  .cg-heatmap-viz .row_label_group,
  .cg-heatmap-viz .col_label_text text,
  .cg-heatmap-viz .row_label_group text { cursor: pointer !important; }
  .cg-heatmap-viz .col_label_text:focus-visible,
  .cg-heatmap-viz .row_label_group:focus-visible,
  .cg-heatmap-viz .col_label_text text:focus-visible,
  .cg-heatmap-viz .row_label_group text:focus-visible {
    outline: 3px solid #0f4c81 !important;
    outline-offset: 2px !important;
    border-radius: 4px;
  }
  .cg-heatmap-viz .cg-active-label text,
  .cg-heatmap-viz text.cg-active-label {
    fill: #053b2d !important;
    font-weight: 800 !important;
    text-decoration: underline;
  }
  .cg-heatmap-viz .cg-safety-search-hit {
    stroke: #9a441f !important;
    stroke-width: 0.35px !important;
  }
  .d3-tip { display: none !important; }
  .cg-cell-value-layer text { pointer-events: none; user-select: none; font-family: Consolas, monospace; font-weight: 600; letter-spacing: -0.01em; }
  .cg-hover-tip {
    position: fixed;
    max-width: 380px;
    pointer-events: none;
    z-index: 9999;
    border: 1px solid #627160;
    border-radius: 10px;
    background: rgba(249,253,247,0.99);
    padding: 9px 11px;
    font-size: 13px;
    line-height: 1.48;
    color: #17251b;
    box-shadow: 0 14px 28px rgba(0,0,0,0.2);
    white-space: normal;
  }
  .cg-hover-tip .tip-row { margin: 0; }
  .cg-hover-tip .tip-key { font-weight: 800; color: #183123; }
  .cg-hover-tip .tip-val { color: #1f3022; }
  .heatmap-debug { margin-top: 8px; }
  .heatmap-debug summary { cursor: pointer; font-weight: 600; }
  .heatmap-debug pre { margin: 6px 0 0; padding: 8px; border: 1px solid #ddd; background: #f5f5f5; max-height: 240px; overflow: auto; white-space: pre-wrap; font-family: Consolas, monospace; font-size: 0.85em; }
  .heatmap-analysis-panel { margin-top: 10px; border: 1px solid #d5dccc; border-radius: 12px; background: linear-gradient(180deg, #fbfcf8 0%, #f3f7ef 100%); box-shadow: 0 8px 20px rgba(22,34,28,0.06); overflow: hidden; }
  .heatmap-analysis-toolbar { display: flex; gap: 8px; align-items: center; justify-content: flex-start; flex-wrap: wrap; padding: 8px 10px; border-bottom: 1px solid #dde5d5; }
  .heatmap-analysis-tabs, .heatmap-analysis-subtabs { display: flex; gap: 6px; flex-wrap: wrap; }
  .heatmap-analysis-toolbar button { border: 1px solid #667e64; border-radius: 999px; background: #fff; color: #19331e; padding: 4px 10px; cursor: pointer; font-size: 0.82em; font-weight: 700; }
  .heatmap-analysis-toolbar button.is-active { background: #0f4c3a; color: #fff; border-color: #0f4c3a; box-shadow: 0 0 0 2px rgba(15,76,58,0.2); }
  .heatmap-analysis-toolbar button:hover { background: #f1f7ee; }
  .heatmap-analysis-toolbar button:focus-visible,
  .heatmap-segmented button:focus-visible,
  .heatmap-control button:focus-visible,
  .heatmap-safety-search button:focus-visible,
  .heatmap-safety-search input:focus-visible,
  .heatmap-motif-search button:focus-visible,
  .heatmap-motif-search input:focus-visible {
    outline: 3px solid #0f4c81;
    outline-offset: 2px;
  }
  .heatmap-analysis-view { display: none; padding: 8px 10px 12px; }
  .heatmap-analysis-view.is-active { display: block; }
  .heatmap-analysis-head { display: grid; gap: 2px; margin-bottom: 8px; }
  .heatmap-analysis-title { font-weight: 700; color: #243a2a; }
  .heatmap-analysis-subtitle { font-size: 0.83em; color: #556555; }
  .heatmap-analysis-subtabs { margin-bottom: 8px; }
  .heatmap-analysis-canvas { width: 100%; min-height: 330px; display: block; }
  .heatmap-analysis-empty { min-height: 330px; display: flex; align-items: center; justify-content: center; color: #566456; font-size: 0.88em; }
  .heatmap-analysis-note { font-size: 0.8em; color: #556555; padding: 0 10px 10px; }
  @media (max-width: 1100px) {
    .cg-heatmap-layout { grid-template-columns: minmax(0, 1fr); }
    .cg-heatmap-wrap { align-items: flex-start; overflow-x: auto; }
    .cg-heatmap-viz { min-height: clamp(520px, 68vh, 820px); }
  }
</style>
""".strip()
