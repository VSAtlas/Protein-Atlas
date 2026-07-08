# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.export_io import _format_num, _resolve_heatmap_source, _stage_artifacts
from analysis.reporting.run_report.html_sections import (
    _count_embedded_blocks,
    _extract_heatmap_asset_mode,
    _render_statistical_analysis_section,
)
from analysis.reporting.run_report.paths import _run_data_dir, _run_manifests_dir
from analysis.reporting.run_report.target_stats import _format_pct_display
def _write_html_report(
    report: Dict[str, Any],
    out_path: Path,
    highlight_queries: List[str],
    repo_root: Path,
    heatmap_source_override: Optional[Path] = None,
    report_asset_mode: str = "auto",
    report_html_warn_bytes: int = REPORT_HTML_WARN_BYTES_DEFAULT,
) -> None:
    logger = logging.getLogger("run-report")
    run_id = report.get("run_id", "")
    generated_at = report.get("generated_at", "")

    compact_rows = []
    compact_highlights = report.get("ligand_highlights", []) or []
    compact_cfg = report.get("ligand_highlight_config", {}) or {}
    top_targets_n = int(compact_cfg.get("top_targets_n") or 3)
    breadth_threshold_raw = as_float(compact_cfg.get("breadth_pct_threshold"))
    if breadth_threshold_raw is None or not math.isfinite(breadth_threshold_raw):
        breadth_threshold_raw = 0.01
    breadth_threshold_label = f"{(breadth_threshold_raw * 100.0):g}%"
    if compact_highlights:
        for row in compact_highlights:
            ligand_name = clean_report_text(row.get("ligand")) or clean_report_text(
                row.get("ligand_display")
            )
            if not ligand_name:
                ligand_name = clean_report_text(row.get("ligand_base")) or "—"
            best_target = clean_report_text(row.get("best_target")) or "—"
            best_score = _format_num(row.get("best_score"), ".6g") or "—"
            best_percentile = clean_report_text(row.get("best_percentile"))
            if not best_percentile:
                best_percentile = (
                    _format_pct_display(
                        row.get("best_pct_rank"),
                        2,
                        min_nonzero_pct=0.01,
                    )
                    or "—"
                )
            top_targets = clean_report_text(row.get("top_targets")) or "—"
            coverage = clean_report_text(row.get("coverage")) or "—"
            breadth = clean_report_text(row.get("breadth")) or "0"
            tail_summary = clean_report_text(row.get("tail_summary")) or "—"
            compact_rows.append(
                "<tr>"
                f"<td>{html.escape(ligand_name)}</td>"
                f"<td>{html.escape(best_target)}</td>"
                f'<td class="num">{html.escape(best_score)}</td>'
                f'<td class="num best-pct">{html.escape(best_percentile)}</td>'
                f'<td class="top-targets-col top-targets">{html.escape(top_targets)}</td>'
                f'<td class="num">{html.escape(coverage)}</td>'
                f'<td class="num">{html.escape(breadth)}</td>'
                f'<td class="num">{html.escape(tail_summary)}</td>'
                "</tr>"
            )

    highlight_rows = []
    highlight_headers = "".join(
        f'<th class="num">{html.escape(q)}</th>' for q in highlight_queries
    )
    if not compact_rows:
        targets = report.get("targets", {}) or {}
        for target_id in sorted(targets.keys()):
            target = targets[target_id]
            target_name = clean_report_text(target.get("target_name")) or "—"
            highlights = target.get("highlights", []) or []
            highlight_map = {h.get("query"): h for h in highlights if h}
            cells = []
            for query in highlight_queries:
                entry = highlight_map.get(query)
                if entry and entry.get("found"):
                    rank = entry.get("rank")
                    pct = _format_num(entry.get("pct_rank"), ".6f")
                    t_sel = _format_num(entry.get("z_selected"), ".6g")
                    if rank:
                        cell = f"{rank} ({pct}) {t_sel}".strip()
                    else:
                        cell = "not found"
                else:
                    cell = "not found"
                cells.append(f'<td class="num">{html.escape(cell)}</td>')
            row_html = f"<tr><td>{html.escape(target_id)}</td><td>{html.escape(target_name)}</td>{''.join(cells)}</tr>"
            highlight_rows.append(row_html)

    heatmap_html = ""
    heatmap_source = heatmap_source_override or _resolve_heatmap_source(
        repo_root, run_id
    )
    if heatmap_source is not None:
        try:
            heatmap_html = render_interactive_heatmap_html(
                repo_root,
                run_id,
                heatmap_source,
                report_asset_mode=report_asset_mode,
            )
        except _REPORT_HTML_BUILD_ERRORS as exc:
            heatmap_html = (
                f'<div class="meta">Heatmap unavailable: {html.escape(str(exc))}</div>'
            )
    else:
        heatmap_html = '<div class="meta">Heatmap data not available.</div>'

    heatmap_section = (
        '<section class="section section-heatmap" id="heatmap" aria-labelledby="heatmap-heading" aria-live="polite">'
        '<h2 id="heatmap-heading">Interactive Docking Heatmap</h2>'
        f"{heatmap_html}"
        "</section>"
    )

    report_suffix = ""
    if out_path.stem.startswith("report_"):
        report_suffix = out_path.stem[len("report_") :]
    data_dir = _run_data_dir(repo_root, run_id)
    if report_suffix:
        report_yaml = data_dir / f"report_{report_suffix}.yaml"
        heatmap_csv_artifact = data_dir / f"heatmap_input_{report_suffix}.csv"
        heatmap_png = data_dir / f"heatmap_{report_suffix}.png"
    else:
        report_yaml = data_dir / "report.yaml"
        heatmap_csv_artifact = data_dir / "heatmap_input.csv"
        heatmap_png = data_dir / "heatmap.png"
    statistical_analysis_html = _render_statistical_analysis_section(repo_root, data_dir)
    config_snapshot = data_dir / "config_snapshot.txt"
    run_manifest = _run_manifests_dir(repo_root, str(run_id)) / "run_manifest.yaml"
    artifacts = [
        ("report.yaml", report_yaml),
        ("heatmap_input.csv", heatmap_csv_artifact),
        ("heatmap.png", heatmap_png),
        ("config snapshot", config_snapshot),
        ("run_manifest", run_manifest),
    ]
    staged = _stage_artifacts(out_path.parent, artifacts)
    staged_links = {label: href for label, href in staged}
    link_items = [
        f'<li><a href="{html.escape(href)}">{html.escape(label)}</a></li>'
        for label, href in staged
    ]
    if link_items:
        artifacts_html = f"<ul class=\"artifact-list\">{''.join(link_items)}</ul>"
    else:
        artifacts_html = '<div class="meta">No artifacts available.</div>'
    effective_asset_mode = _extract_heatmap_asset_mode(heatmap_html, report_asset_mode)
    technical_metadata_html = (
        "<ul class=\"artifact-list\">"
        f"<li><span>asset_mode={html.escape(effective_asset_mode)}</span></li>"
        f"<li><span>asset_mode_requested={html.escape(report_asset_mode)}</span></li>"
        "</ul>"
    )
    config_snapshot_html = (
        f'<a href="{html.escape(staged_links["config snapshot"])}">config snapshot</a>'
        if "config snapshot" in staged_links
        else "not available"
    )
    run_manifest_html = (
        f'<a href="{html.escape(staged_links["run_manifest"])}">run_manifest</a>'
        if "run_manifest" in staged_links
        else "not available"
    )
    input_source_html = (
        f'<a href="{html.escape(staged_links["heatmap_input.csv"])}">heatmap_input.csv</a>'
        if "heatmap_input.csv" in staged_links
        else "heatmap input not available"
    )
    git_commit = clean_report_text(report.get("git_commit")) or clean_report_text(
        os.environ.get("ATLAS_GIT_COMMIT") or os.environ.get("GITHUB_SHA")
    )
    if git_commit:
        git_commit = git_commit[:12]
    provenance_html = (
        "<dl class=\"provenance-list\">"
        f"<dt>run_id</dt><dd>{html.escape(str(run_id))}</dd>"
        f"<dt>generated_at</dt><dd>{html.escape(str(generated_at))}</dd>"
        f"<dt>Atlas version</dt><dd>{html.escape(clean_report_text(report.get('atlas_version')) or 'unknown')}</dd>"
        f"<dt>git commit</dt><dd>{html.escape(git_commit or 'unknown')}</dd>"
        f"<dt>config snapshot</dt><dd>{config_snapshot_html}</dd>"
        f"<dt>run_manifest</dt><dd>{run_manifest_html}</dd>"
        f"<dt>input source</dt><dd>{input_source_html}</dd>"
        "</dl>"
    )

    summary = report.get("summary", {}) or {}
    visible_ligands = len(compact_rows) if compact_rows else int(
        summary.get("n_unique_ligands") or 0
    )
    target_count = int(summary.get("n_target_combos") or 0)
    significant_hits_text = ""
    breadth_values: List[int] = []
    for row in compact_highlights:
        raw_breadth = as_float(row.get("breadth"))
        if raw_breadth is not None and math.isfinite(raw_breadth):
            breadth_values.append(int(raw_breadth))
    if breadth_values:
        significant_hits_text = str(sum(breadth_values))
    selected_score_text = clean_report_text(report.get("selected_score_mode")) or clean_report_text(
        report.get("selected_score_column")
    )

    kpi_items = [
        ("Run ID", str(run_id)),
        ("Generated", str(generated_at)),
        ("Visible ligands", str(visible_ligands)),
        ("Targets", str(target_count)),
    ]
    if significant_hits_text:
        kpi_items.append(("Significant hits", significant_hits_text))
    if selected_score_text:
        kpi_items.append(("Selected score", selected_score_text))
    kpi_html = "".join(
        f'<div class="kpi-card"><div class="kpi-label">{html.escape(label)}</div><div class="kpi-value">{html.escape(value)}</div></div>'
        for label, value in kpi_items
    )

    html_body = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Atlas2 Report: {html.escape(str(run_id))}</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f6f7fb;
      --card-bg: #ffffff;
      --text: #1f2937;
      --muted: #6b7280;
      --border: #e5e7eb;
      --accent: #0f172a;
      --highlight: #eef2ff;
    }}
    * {{ box-sizing: border-box; }}
    .skip-link {{ position: absolute; left: 12px; top: 8px; transform: translateY(-140%); background: #0f172a; color: #fff; padding: 8px 12px; border-radius: 8px; z-index: 9999; }} .skip-link:focus {{ transform: translateY(0); }} a:focus-visible, button:focus-visible, input:focus-visible, select:focus-visible {{ outline: 3px solid #2563eb; outline-offset: 2px; }} @media (prefers-reduced-motion: reduce) {{ *, *::before, *::after {{ animation-duration: 0.01ms !important; animation-iteration-count: 1 !important; scroll-behavior: auto !important; transition-duration: 0.01ms !important; }} }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      margin: 0;
      background: var(--bg);
      color: var(--text);
      line-height: 1.5;
    }}
    h1, h2 {{ margin: 0; color: var(--accent); }}
    .container {{ width: 100% !important; max-width: none !important; margin: 0 auto !important; padding: 14px clamp(2px, 0.5vw, 8px) 28px !important; box-sizing: border-box; }}
    .page-header {{
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: baseline;
      gap: 12px;
      margin-bottom: 16px;
    }}
    .page-title {{ font-size: 2.1rem; }}
    .report-kpis {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 8px;
      margin-bottom: 14px;
    }}
    .kpi-card {{
      border: 1px solid var(--border);
      border-radius: 10px;
      background: #f8fafc;
      padding: 8px 10px;
    }}
    .kpi-label {{ color: var(--muted); font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.02em; }}
    .kpi-value {{ font-weight: 700; font-size: 1rem; }}
    .limitations-note {{
      margin-bottom: 14px;
      border: 1px solid #fde68a;
      border-radius: 10px;
      background: #fffbeb;
      padding: 10px 12px;
    }}
    .limitations-note h2 {{ font-size: 1rem; margin-bottom: 6px; }}
    .limitations-note ul {{ margin: 0; padding-left: 18px; }}
    .meta {{ color: var(--muted); font-size: 0.95rem; }}
    .toc {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin: 16px 0 20px;
      padding: 10px 12px;
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 12px;
    }}
    .toc a {{
      text-decoration: none;
      color: var(--accent);
      background: #f0f3f8;
      border: 1px solid #e2e8f0;
      padding: 4px 10px;
      border-radius: 999px;
      font-size: 0.9rem;
    }}
    .toc a:hover {{ background: #e8ecf4; }}
    .section {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 14px 14px;
      margin-bottom: 18px;
      box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
    }}
    #heatmap.section {{ padding: 8px 4px !important; width: 100%; max-width: none; }}
    #heatmap.section h2 {{ padding: 0 2px; }}
    .section h2 {{ margin-bottom: 12px; font-size: 1.3rem; }}
    .section h3 {{ margin: 14px 0 8px; font-size: 1rem; color: #334155; }}
    .stats-explainer {{ display: grid; gap: 8px; margin-bottom: 12px; }}
    .stats-explainer p {{ margin: 0; color: #334155; }}
    .formula-panel {{ border: 1px solid var(--border); border-radius: 10px; background: #f8fafc; padding: 10px 12px; margin: 12px 0; }}
    .formula-panel summary {{ cursor: pointer; font-weight: 700; color: #0f172a; }}
    .formula-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 10px; margin-top: 10px; }}
    .formula-grid div {{ background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 10px; }}
    .formula-grid h4 {{ margin: 0 0 6px; font-size: 0.92rem; color: #1e293b; }}
    .formula-grid p {{ margin: 4px 0; color: #475569; line-height: 1.35; }}
    .formula-note {{ margin: 8px 0 0; color: #475569; line-height: 1.35; }}
    .formula-source {{ color: #64748b !important; font-size: 0.78rem; }}
    .table-wrap {{
      overflow-x: auto;
      border: 1px solid var(--border);
      border-radius: 10px;
      background: #ffffff;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      font-size: 0.95rem;
      font-variant-numeric: tabular-nums;
    }}
    th, td {{ padding: 8px 10px; text-align: left; }}
    thead th {{
      position: sticky;
      top: 0;
      background: #f4f6fa;
      border-bottom: 1px solid var(--border);
      z-index: 1;
    }}
    tbody td {{ border-bottom: 1px solid var(--border); }}
    tbody tr:nth-child(even) {{ background: #f9fafb; }}
    tbody tr:hover {{ background: var(--highlight); }}
    th.num, td.num {{ text-align: right; }}
    .highlights-table th, .highlights-table td {{ vertical-align: top; }}
    .highlights-table td:first-child {{
      max-width: 24rem;
      white-space: normal;
      word-break: break-word;
      overflow-wrap: anywhere;
    }}
    .highlights-table .best-pct {{
      white-space: nowrap;
      min-width: 7.5rem;
      padding-right: 14px;
    }}
    .highlights-table .top-targets-col {{
      min-width: 22rem;
      padding-left: 14px;
    }}
    .top-targets {{
      white-space: pre-wrap;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
      font-size: 0.86rem;
      line-height: 1.55;
      text-align: left;
    }}
    .section-heatmap .heatmap-scale-meta,
    .section-heatmap .heatmap-report-hero,
    .section-heatmap .cg-heatmap-sidebar,
    .section-heatmap .heatmap-analysis-panel,
    .section-heatmap .heatmap-report-panel {{
      border-radius: 12px;
      box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
    }}
    .artifact-list {{
      list-style: none;
      padding-left: 0;
      margin: 0;
      display: grid;
      gap: 6px;
    }}
    .artifact-list a {{
      text-decoration: none;
      color: var(--accent);
      background: #f8fafc;
      border: 1px solid var(--border);
      padding: 6px 10px;
      border-radius: 8px;
      display: inline-block;
    }}
    .artifact-list a:hover {{ background: #eef2f7; }}
    .provenance-list {{ display: grid; grid-template-columns: max-content minmax(0, 1fr); gap: 6px 12px; margin: 0; }}
    .provenance-list dt {{ color: var(--muted); font-weight: 700; }}
    .provenance-list dd {{ margin: 0; overflow-wrap: anywhere; }}
    @media print {{
      body {{ background: #fff; }}
      .toc,
      .heatmap-report-hero,
      .heatmap-analysis-toolbar,
      .heatmap-discovery-controls,
      .heatmap-workflow-actions,
      .heatmap-side-effect-actions,
      .heatmap-debug,
      .heatmap-scale-toggle,
      button,
      input,
      select {{
        display: none !important;
      }}
      .section,
      .limitations-note,
      .report-kpis {{
        break-inside: avoid;
        page-break-inside: avoid;
      }}
      a[href]::after {{ content: " (" attr(href) ")"; font-size: 0.82em; color: #4b5563; }}
    }}
  </style>
</head>
<body>
  <a class="skip-link" href="#heatmap">Skip to heatmap</a>
  <main class="container" id="report-content">
    <header class="page-header">
      <h1 class="page-title">Atlas2 Report: {html.escape(str(run_id))}</h1>
      <div class="meta">Generated at {html.escape(str(generated_at))}</div>
    </header>
    <section class="report-kpis" aria-label="Run summary KPIs">{kpi_html}</section>
    <section class="limitations-note" aria-label="Interpretation and limitations">
      <h2>Interpretation and Limitations</h2>
      <ul>
        <li>Docking outputs are computational hypotheses and should guide prioritization, not establish biological truth.</li>
        <li>Atlas score is relative within the analyzed run and score mode; check Atlas score source for the underlying score provenance.</li>
        <li>Percent rank is computed within the screened/background distribution represented in this run.</li>
        <li>Experimental validation is required before drawing mechanistic or translational conclusions.</li>
      </ul>
    </section>
    <details class="section interpretation-guide">
      <summary>How to interpret this report</summary>
      <dl>
        <dt>Atlas score</dt><dd>Relative run-level score used to prioritize ligand-target hypotheses; Atlas score source records which upstream score produced it.</dd>
        <dt>Percentile rank</dt><dd>Rank within the screened/background distribution; lower percentiles indicate stronger relative signals.</dd>
        <dt>Coverage</dt><dd>Number of targets with usable score evidence for a ligand.</dd>
        <dt>Breadth 1%</dt><dd>Count of target pairs meeting the significant hit threshold, usually the top 1% unless configured otherwise.</dd>
        <dt>Target organization</dt><dd>Optional grouping by family, pathway, ADME category, or safety bucket to aid scanning.</dd>
        <dt>Pose validity</dt><dd>Indicates whether pose-quality checks were available and passed when exported by the run.</dd>
        <dt>MM/GBSA</dt><dd>Optional rescoring evidence when present; absence does not change the docking score definitions.</dd>
      </dl>
      <p>Rankings prioritize hypotheses, not validated binding or experimental affinity.</p>
    </details>

    <nav class="toc" aria-label="Report sections">
      <a href="#highlights">Top Ligand-Target Signals</a>
      <a href="#heatmap">Interactive Docking Heatmap</a>
      <a href="#statistics">Statistical Analysis</a>
      <a href="#artifacts">Run Artifacts</a>
      <a href="#provenance">Provenance</a>
      <a href="#technical">Technical metadata</a>
    </nav>

    <section class="section" id="highlights">
      <h2>Top Ligand-Target Signals</h2>
      <div class="table-wrap">
        <table class="highlights-table">
          <thead>
            <tr>
              {'<th>Ligand</th><th>Best target</th><th class="num">Best score</th><th class="num best-pct">Best percentile</th>'
                + f'<th class="top-targets-col">Top targets (N={top_targets_n})</th><th class="num">Coverage</th><th class="num">Breadth {breadth_threshold_label}</th><th class="num">Tail summary</th>'
                if compact_rows else '<th>target_id</th><th>target_name</th>' + highlight_headers}
            </tr>
          </thead>
          <tbody>
            {''.join(compact_rows) if compact_rows else ''.join(highlight_rows)}
          </tbody>
        </table>
      </div>
    </section>

    {heatmap_section}

    {statistical_analysis_html}

    <section class="section" id="artifacts">
      <h2>Run Artifacts</h2>
      {artifacts_html}
    </section>
    <section class="section" id="provenance">
      <h2>Provenance</h2>
      {provenance_html}
    </section>
    <section class="section" id="technical">
      <h2>Technical metadata</h2>
      {technical_metadata_html}
    </section>
  </main>
</body>
</html>
"""

    report_bytes = len(html_body.encode("utf-8"))
    embedded_blocks = _count_embedded_blocks(html_body)
    logger.info(
        "%s action=report_html_metrics path=%s bytes=%d embedded_blocks=%d asset_mode=%s",
        COMPONENT,
        out_path,
        report_bytes,
        embedded_blocks,
        effective_asset_mode,
    )
    if report_html_warn_bytes > 0 and report_bytes > report_html_warn_bytes:
        logger.warning(
            "%s action=report_html_size_warn path=%s bytes=%d threshold=%d asset_mode=%s",
            COMPONENT,
            out_path,
            report_bytes,
            report_html_warn_bytes,
            effective_asset_mode,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(html_body)

