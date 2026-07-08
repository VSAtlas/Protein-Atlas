# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.export_io import _format_num
def _extract_heatmap_asset_mode(heatmap_html: str, default_mode: str) -> str:
    match = re.search(
        r'<script type="application/json" id="cg-heatmap-assets-[^"]+">(.*?)</script>',
        heatmap_html,
        re.S,
    )
    if match is None:
        return default_mode
    try:
        payload = json.loads(match.group(1))
    except json.JSONDecodeError:
        return default_mode
    mode = str(payload.get("mode") or "").strip().lower()
    return mode or default_mode

def _read_existing_csv(path: Path, limit: int = 5000) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size <= 0:
        return []
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for idx, row in enumerate(reader):
            if idx >= limit:
                break
            rows.append(dict(row))
    return rows

def _candidate_stat_paths(repo_root: Path, data_dir: Path, filename: str) -> List[Path]:
    return [
        data_dir / filename,
        data_dir / "enrichment" / filename,
        repo_root / "outputs" / "analysis" / filename,
        repo_root / "outputs" / "analysis" / "enrichment" / filename,
        repo_root / "outputs" / "analysis" / "external_labels" / filename,
    ]

def _first_existing_rows(
    repo_root: Path, data_dir: Path, filename: str, limit: int = 5000
) -> Tuple[Path, List[Dict[str, str]]]:
    for path in _candidate_stat_paths(repo_root, data_dir, filename):
        rows = _read_existing_csv(path, limit=limit)
        if rows:
            return path, rows
    return Path(), []

def _first_existing_rows_bundle(
    repo_root: Path,
    data_dir: Path,
    filenames: Tuple[str, ...],
    limit: int = 5000,
) -> Dict[str, Tuple[Path, List[Dict[str, str]]]]:
    return {
        name: _first_existing_rows(repo_root, data_dir, name, limit=limit)
        for name in filenames
    }

def _render_compact_table(
    rows: List[Dict[str, Any]], columns: List[str], *, max_rows: int = 12
) -> str:
    if not rows:
        return '<div class="meta">No rows available.</div>'
    header = "".join(f"<th>{html.escape(_report_display_label(col))}</th>" for col in columns)
    body_rows = []
    for row in rows[:max_rows]:
        cells = "".join(
            f"<td>{html.escape(_report_display_label(clean_report_text(row.get(col)) or '—'))}</td>"
            for col in columns
        )
        body_rows.append(f"<tr>{cells}</tr>")
    if len(rows) > max_rows:
        body_rows.append(
            f'<tr><td colspan="{len(columns)}" class="meta">Showing {max_rows} of {len(rows)} rows.</td></tr>'
        )
    return f'<div class="table-wrap"><table class="stats-table"><thead><tr>{header}</tr></thead><tbody>{"".join(body_rows)}</tbody></table></div>'

def _report_display_label(value: str) -> str:
    labels = {
        "z_selected": "Atlas score",
        "z_selected_source": "Atlas score source",
        "z_selected_raw": "Raw Atlas score",
        "z_selected_fixed": "Fixed-scale Atlas score",
        "z_selected_quantile": "Quantile-scale Atlas score",
    }
    return labels.get(value, value)

def _kpi_cards_html(cards: List[Tuple[str, str, str]]) -> str:
    return "".join(
        f'<div class="kpi-card"><div class="kpi-label">{html.escape(label)}</div>'
        f'<div class="kpi-value">{html.escape(value)}</div>'
        f'<div class="meta">{html.escape(note)}</div></div>'
        for label, value, note in cards
    )

def _fdr_summary_cards(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return '<div class="meta">No target FDR summary found.</div>'
    reliable = sum(1 for row in rows if as_report_row_bool(row.get("fdr_reliable")))
    q05 = sum(int(as_float(row.get("n_hits_q05")) or 0) for row in rows)
    q10 = sum(int(as_float(row.get("n_hits_q10")) or 0) for row in rows)
    competition_q05 = sum(
        int(as_float(row.get("n_decoy_competition_hits_q05")) or 0) for row in rows
    )
    competition_q10 = sum(
        int(as_float(row.get("n_decoy_competition_hits_q10")) or 0) for row in rows
    )
    resolvable_q10 = sum(1 for row in rows if as_report_row_bool(row.get("bh_resolvable_q10")))
    raw_min_q_values = [as_float(row.get("min_possible_bh_q")) for row in rows]
    min_q_values = [
        float(v) for v in raw_min_q_values if v is not None and math.isfinite(v)
    ]
    median_min_q = median_floats(min_q_values)
    decoys = [as_float(row.get("fdr_n_decoys")) for row in rows]
    decoy_values = [float(v) for v in decoys if v is not None and math.isfinite(v)]
    median_decoys = median_floats(decoy_values)
    source_counts: Dict[str, int] = {}
    for row in rows:
        source = clean_report_text(row.get("fdr_null_source")) or "unknown"
        source_counts[source] = source_counts.get(source, 0) + 1
    source_text = "; ".join(f"{key}: {value}" for key, value in sorted(source_counts.items()))
    cards = [
        ("FDR reliable targets", f"{reliable} / {len(rows)}", "Requires sufficient decoy count and unique decoy scores."),
        ("q <= 0.05 hits", str(q05), "Target-level BH q-value hits across non-decoy rows."),
        ("q <= 0.10 hits", str(q10), "Nominal discovery set used for report triage."),
        (
            "BH resolvable targets",
            f"{resolvable_q10} / {len(rows)}",
            f"Targets where the decoy empirical p-value floor can reach q <= 0.10; median minimum q={_format_num(median_min_q, '.6g') or '—'}.",
        ),
        (
            "Decoy-competition hits",
            f"{competition_q05} / {competition_q10}",
            "Exploratory target-decoy competition hits at 5% / 10%; use alongside BH q-values.",
        ),
        ("Median decoys", _format_num(median_decoys, ".6g") or "—", source_text),
    ]
    return _kpi_cards_html(cards)

def _scorch_fdr_summary_cards(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return '<div class="meta">No SCORCH conditional FDR summary found.</div>'
    q05 = sum(
        int(as_float(row.get("n_decoy_competition_hits_q05")) or 0) for row in rows
    )
    q10 = sum(
        int(as_float(row.get("n_decoy_competition_hits_q10")) or 0) for row in rows
    )
    q05_plus1 = sum(
        int(as_float(row.get("n_decoy_competition_hits_q05_plus1")) or 0)
        for row in rows
    )
    q10_plus1 = sum(
        int(as_float(row.get("n_decoy_competition_hits_q10_plus1")) or 0)
        for row in rows
    )
    decoys = [as_float(row.get("n_decoys")) for row in rows]
    decoy_values = [float(v) for v in decoys if v is not None and math.isfinite(v)]
    tested = [as_float(row.get("n_tested")) for row in rows]
    tested_values = [float(v) for v in tested if v is not None and math.isfinite(v)]
    cards = [
        (
            "SCORCH relaxed q <= 0.10",
            str(q10),
            "Conditional target-decoy competition hits without the finite-sample +1 correction.",
        ),
        (
            "SCORCH +1 q <= 0.10",
            str(q10_plus1),
            "More conservative conditional SCORCH hits using +1 decoy-count correction.",
        ),
        (
            "SCORCH relaxed q <= 0.05",
            str(q05),
            "Exploratory 5% conditional target-decoy competition hits.",
        ),
        (
            "SCORCH +1 q <= 0.05",
            str(q05_plus1),
            "Finite-sample corrected 5% conditional SCORCH hits.",
        ),
        (
            "Median SCORCH decoys/tested",
            f"{_format_num(median_floats(decoy_values), '.6g') or '—'} / {_format_num(median_floats(tested_values), '.6g') or '—'}",
            "Rows are limited to the post-SCORCH tranche, not the full FDA/DUD library.",
        ),
    ]
    return _kpi_cards_html(cards)

def _binary_label(value: Any) -> Optional[int]:
    text = clean_report_text(value).lower()
    if text in {"1", "true", "yes", "y", "positive", "active"}:
        return 1
    if text in {"0", "false", "no", "n", "negative", "inactive"}:
        return 0
    return None

def _rank_metric_rows(
    rows: List[Dict[str, str]], score_cols: List[str], label_col: str
) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    top_fractions = (0.01, 0.05, 0.10)
    for score_col in score_cols:
        scores: List[float] = []
        labels: List[int] = []
        for row in rows:
            score = as_float(row.get(score_col))
            label = _binary_label(row.get(label_col))
            if score is None or label is None:
                continue
            scores.append(float(score))
            labels.append(int(label))
        if not scores:
            continue
        positives = sum(labels)
        if positives <= 0:
            continue
        overall = positives / len(scores)
        metrics = ranked_binary_metrics(scores, labels, top_fractions)
        metric_row: Dict[str, str] = {
            "score": score_col,
            "label": label_col,
            "metric": "derived_from_spd_benchmark",
            "n_pairs": str(len(scores)),
            "n_labeled_positive": str(positives),
            "positive_rate_overall": _format_num(overall, ".6g"),
            "AUPRC": _format_num(metrics["AUPRC"], ".6g"),
        }
        for frac in top_fractions:
            key = f"EF@{int(frac * 100)}%"
            metric_row[key] = _format_num(metrics[key], ".6g")
        out.append(metric_row)
    return out

def _derived_spd_rows(repo_root: Path, data_dir: Path) -> List[Dict[str, str]]:
    for path in _candidate_stat_paths(repo_root, data_dir, "spd_atlas_benchmark.csv"):
        rows = _read_existing_csv(path, limit=200000)
        if rows:
            return rows
    return []

def _render_statistical_analysis_section(repo_root: Path, data_dir: Path) -> str:
    stat_bundle = _first_existing_rows_bundle(
        repo_root,
        data_dir,
        (
            "target_fdr_summary.csv",
            "target_scorch_fdr_summary.csv",
            "enrichment_summary.csv",
            "ablation_summary.csv",
            "matched_control_significance.csv",
            "source_coverage_summary.csv",
            "heatmap_ml_readiness_summary.csv",
        ),
    )
    fdr_path, fdr_rows = stat_bundle["target_fdr_summary.csv"]
    scorch_fdr_path, scorch_fdr_rows = stat_bundle["target_scorch_fdr_summary.csv"]
    enrichment_path, enrichment_rows = stat_bundle["enrichment_summary.csv"]
    ablation_path, ablation_rows = stat_bundle["ablation_summary.csv"]
    matched_path, matched_rows = stat_bundle["matched_control_significance.csv"]
    coverage_path, coverage_rows = stat_bundle["source_coverage_summary.csv"]
    readiness_path, readiness_rows = stat_bundle["heatmap_ml_readiness_summary.csv"]
    spd_rows = _derived_spd_rows(repo_root, data_dir)
    if not enrichment_rows and spd_rows:
        enrichment_rows = _rank_metric_rows(
            spd_rows,
            ["z_selected", "consensus_score", "final_score", "SCORCH_score_used"],
            "spd_exposure_relevant",
        )
    if not ablation_rows and spd_rows:
        ablation_rows = [
            {
                "ablation": f"score_comparison:{row.get('score', '')}",
                "score_col": row.get("score", ""),
                "label": row.get("label", ""),
                "n_non_missing": row.get("n_pairs", ""),
                "n_labeled": row.get("n_pairs", ""),
                "n_positive": row.get("n_labeled_positive", ""),
                "EF@1%": row.get("EF@1%", ""),
                "EF@5%": row.get("EF@5%", ""),
                "AUPRC": row.get("AUPRC", ""),
                "AUROC": "",
            }
            for row in enrichment_rows
        ]
        ablation_note = (
            "No true ablation_summary.csv was found; these rows are score "
            "comparisons derived from SPD enrichment, not component-removal "
            "ablation rerankings."
        )
    else:
        ablation_note = (
            "Rows are loaded from ablation_summary.csv when available. "
            "True ablations should recompute the score after removing each "
            "component and reranking."
        )
    source_note = "; ".join(
        f"{label}: {path.name}"
        for label, path in (
            ("FDR", fdr_path),
            ("SCORCH conditional FDR", scorch_fdr_path),
            ("enrichment", enrichment_path),
            ("ablation", ablation_path),
            ("matched controls", matched_path),
            ("source coverage", coverage_path),
            ("heatmap ML readiness", readiness_path),
        )
        if path
    )
    return f"""
    <section class="section" id="statistics">
      <h2>Statistical Analysis</h2>
      <div class="stats-explainer">
        <p><strong>Atlas score source:</strong> source files may contain several score columns. The report displays the canonical master-export score as Atlas score: Stage-2/post-SCORCH z-score when available, then Stage-1 consensus z-score, then <code>consensus_score</code> only as a non-z fallback. The companion Atlas score source column records which source was used.</p>
        <p><strong>Sawada/PBAS context:</strong> PBAS-style profiles are drug-level vectors of predicted binding affinities across proteins. Atlas uses that framing as a baseline and adds pair-level significance, Vina/SCORCH consensus, MM/GBSA, exposure, tissue, target-ADR, pathway, and mechanism-graph evidence.</p>
      </div>
      {render_statistical_formula_block()}
      <div class="report-kpis">{_fdr_summary_cards(fdr_rows)}</div>
      <h3>Target FDR / Decoy Nulls</h3>
      {_render_compact_table(fdr_rows, ["pdb_id", "variant", "ph_label", "fdr_score_field", "fdr_null_source", "fdr_n_decoys", "fdr_reliable", "n_tested", "min_possible_bh_q", "bh_resolvable_q10", "n_hits_q05", "n_hits_q10", "n_decoy_competition_hits_q05", "n_decoy_competition_hits_q10", "best_decoy_competition_q"], max_rows=16)}
      <h3>SCORCH Conditional FDR</h3>
      <p class="muted">This uses only rows that were actually rescored by SCORCH and their rescored DUD rows. It is a conditional top-tranche null, not a full FDA-library null.</p>
      <div class="report-kpis">{_scorch_fdr_summary_cards(scorch_fdr_rows)}</div>
      {_render_compact_table(scorch_fdr_rows, ["pdb_id", "variant", "ph_label", "score_field", "scope", "n_decoys", "n_tested", "min_possible_bh_q", "best_bh_q", "best_decoy_competition_q", "best_decoy_competition_q_plus1", "n_decoy_competition_hits_q05", "n_decoy_competition_hits_q10", "n_decoy_competition_hits_q05_plus1", "n_decoy_competition_hits_q10_plus1"], max_rows=16)}
      <h3>Enrichment and Label-Shuffle Controls</h3>
      {_render_compact_table(enrichment_rows, ["score", "label", "metric", "observed", "EF@1%", "EF@5%", "EF@10%", "AUPRC", "ci_low", "ci_high", "drug_label_shuffle_p", "target_label_shuffle_p", "degree_preserving_label_shuffle_p"], max_rows=14)}
      <h3>Ablation Reranking</h3>
      <p class="muted">{html.escape(ablation_note)}</p>
      {_render_compact_table(ablation_rows, ["ablation", "score_col", "label", "n_non_missing", "n_labeled", "n_positive", "EF@1%", "EF@5%", "AUPRC", "AUROC"], max_rows=14)}
      <h3>Matched Controls</h3>
      {_render_compact_table(matched_rows, ["drug_id", "target_id", "score", "matched_p", "n_controls", "match_fields"], max_rows=12)}
      <h3>External Source Coverage</h3>
      {_render_compact_table(coverage_rows, ["source", "status_column", "n_rows", "n_drugs", "n_targets", "n_drug_target_pairs", "status_counts"], max_rows=8)}
      <h3>Coverage Audit / ML Usability</h3>
      {_render_compact_table(readiness_rows, ["metric", "value"], max_rows=20)}
      <div class="meta">Statistical source files: {html.escape(source_note or "none found beside this report or in outputs/analysis")}.</div>
    </section>
"""

def _count_embedded_blocks(html_body: str) -> int:
    styles = len(re.findall(r"<style(?:\s[^>]*)?>", html_body, re.IGNORECASE))
    scripts = len(
        re.findall(r"<script(?![^>]*\bsrc=)(?:\s[^>]*)?>", html_body, re.IGNORECASE)
    )
    return styles + scripts

