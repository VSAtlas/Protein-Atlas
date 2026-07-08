import csv
import json
import re
from pathlib import Path

import pytest

from analysis.reporting import heatmap_html as heatmap_html
from analysis.reporting.heatmap_html import (
    _build_clustergrammer_mode_mat,
    _build_promiscuity_meta,
    render_interactive_heatmap_html,
)


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _build_minimal_rows() -> list[dict[str, str]]:
    return [
        {
            "target_id": "T1",
            "pdb_id": "1ABC",
            "target_name": "Kinase Alpha",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "ligand_display": "triamcinolone",
            "ligand_base": "triamcinolone",
            "z_selected": "1.5",
            "rank": "1",
            "pct_rank": "0.1",
            "pose_valid_any": "true",
            "library": "lib_a",
        },
        {
            "target_id": "T2",
            "pdb_id": "2XYZ",
            "target_name": "Kinase Beta",
            "variant": "APO",
            "ph_label": "pH6_5",
            "ligand_display": "triamcinolone",
            "ligand_base": "triamcinolone",
            "z_selected": "2.0",
            "rank": "2",
            "pct_rank": "0.2",
            "pose_valid_any": "false",
            "library": "lib_a",
        },
        {
            "target_id": "T1",
            "pdb_id": "1ABC",
            "target_name": "Kinase Alpha",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "ligand_display": "itraconazole",
            "ligand_base": "itraconazole",
            "z_selected": "0.5",
            "rank": "3",
            "pct_rank": "0.3",
            "pose_valid_any": "true",
            "library": "lib_b",
        },
        {
            "target_id": "T2",
            "pdb_id": "2XYZ",
            "target_name": "Kinase Beta",
            "variant": "APO",
            "ph_label": "pH6_5",
            "ligand_display": "itraconazole",
            "ligand_base": "itraconazole",
            "z_selected": "3.5",
            "rank": "4",
            "pct_rank": "0.4",
            "pose_valid_any": "true",
            "library": "lib_b",
        },
    ]


def _mock_pathway_memberships(
    _repo_root: Path, _pdb_ids: object, **_kwargs: object
) -> dict[str, list[str]]:
    return {
        "1ABC": ["Cell Cycle", "MAPK signaling"],
        "2XYZ": ["MAPK signaling"],
        "1AAA": ["DNA repair"],
        "2BBB": [],
    }


def _mock_target_uniprots_by_pdb(
    _repo_root: Path, _pdb_ids: object
) -> dict[str, list[str]]:
    return {
        "1ABC": ["O75469"],
        "2XYZ": ["P03372"],
        "1AAA": ["P07550"],
        "2BBB": [],
    }


def test_heatmap_html_clustergrammer_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("clustergrammer")
    repo_root = tmp_path / "repo"
    run_id = "CG_TEST"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_pathway_memberships_by_pdb",
        _mock_pathway_memberships,
    )
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_target_uniprots_by_pdb",
        _mock_target_uniprots_by_pdb,
    )

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_minimal_rows())

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=5, include_decoys=False
    )

    assert "Clustergrammer(" in html_text
    safe_run_id = "CG_TEST"
    assert f'id="cg-heatmap-{safe_run_id}"' in html_text

    viz_match = re.search(
        r'<script type="application/json" id="cg-heatmap-data-CG_TEST">(.*?)</script>',
        html_text,
        re.S,
    )
    assert viz_match is not None
    viz_payloads = json.loads(viz_match.group(1))
    assert "fixed" in viz_payloads
    assert viz_payloads["fixed"]["mat"]
    assert "quantile" in viz_payloads

    meta_match = re.search(
        r'<script type="application/json" id="cg-heatmap-meta-CG_TEST">(.*?)</script>',
        html_text,
        re.S,
    )
    assert meta_match is not None
    cell_meta = json.loads(meta_match.group(1))
    alpha_key = next(
        key for key in cell_meta if key.startswith("triamcinolone||1ABC")
    )
    assert cell_meta[alpha_key]["rank"] == "1"

    rowmeta_match = re.search(
        r'<script type="application/json" id="cg-heatmap-rowmeta-CG_TEST">(.*?)</script>',
        html_text,
        re.S,
    )
    assert rowmeta_match is not None
    row_meta = json.loads(rowmeta_match.group(1))
    assert row_meta["triamcinolone"]["ligand"] == "triamcinolone"
    assert row_meta["triamcinolone"]["chemotype_primary"] == "Corticosteroid-like"
    assert row_meta["triamcinolone"]["drug_effect_primary"] == "Glucocorticoid"
    assert "steroidal core" in row_meta["triamcinolone"]["motif_tags"]
    assert row_meta["triamcinolone"]["significant_hits_text"] == "0 / 2"
    assert row_meta["itraconazole"]["chemotype_primary"] == "Azole-like"
    assert row_meta["itraconazole"]["drug_effect_primary"] == "Antifungal"

    colmeta_match = re.search(
        r'<script type="application/json" id="cg-heatmap-colmeta-CG_TEST">(.*?)</script>',
        html_text,
        re.S,
    )
    assert colmeta_match is not None
    col_meta = json.loads(colmeta_match.group(1))
    alpha_col_key = next(key for key in col_meta if key.startswith("1ABC\nKinase Alpha"))
    assert col_meta[alpha_col_key]["pdb"] == "1ABC"
    assert col_meta[alpha_col_key]["target_name"] == "Kinase Alpha"
    assert col_meta[alpha_col_key]["variant_ph"] == "HOLO | pH7_0"
    assert col_meta[alpha_col_key]["protein_family"] == "Kinase"
    assert col_meta[alpha_col_key]["adme_category"] == "Modifier/regulator"
    assert col_meta[alpha_col_key]["safety_buckets"] == ["Unassigned"]
    assert col_meta[alpha_col_key]["primary_display_safety"] == "Unassigned"
    assert col_meta[alpha_col_key]["secondary_safety_buckets"] == []
    assert col_meta[alpha_col_key]["direct_safety_buckets"] == []
    assert col_meta[alpha_col_key]["drug_ae_buckets"] == []
    assert col_meta[alpha_col_key]["safety_confidence"] == "low"
    assert col_meta[alpha_col_key]["safety_sources"] == ["Local heuristic catalog"]
    assert col_meta[alpha_col_key]["direct_liability_examples"] == []
    assert col_meta[alpha_col_key]["raw_adverse_event_examples"] == []
    assert col_meta[alpha_col_key]["pathway_memberships"] == [
        "Cell Cycle",
        "MAPK signaling",
    ]
    assert col_meta[alpha_col_key]["primary_display_pathway"] == "MAPK signaling"
    assert col_meta[alpha_col_key]["significant_hits_text"] == "0 / 2"

    analysis_match = re.search(
        r'<script type="application/json" id="cg-heatmap-analysis-CG_TEST">(.*?)</script>',
        html_text,
        re.S,
    )
    assert analysis_match is not None
    analysis_meta = json.loads(analysis_match.group(1))
    assert analysis_meta["target_organization"]["default_mode"] == "none"
    assert analysis_meta["target_organization"]["adme_blocks"][0]["label"] == "Modifier/regulator"
    assert analysis_meta["target_organization"]["safety_blocks"][-1]["label"] == "Unassigned"
    assert analysis_meta["target_organization"]["family_blocks"][0]["label"] == "Kinase"
    assert analysis_meta["target_organization"]["pathway_blocks"][0]["label"] == "MAPK signaling"
    assert analysis_meta["ligand_organization"]["default_mode"] == "none"
    assert analysis_meta["ligand_organization"]["chemotype_blocks"][0]["label"] == "Corticosteroid-like"
    assert analysis_meta["ligand_organization"]["drug_effect_blocks"][0]["label"] == "Glucocorticoid"
    assert analysis_meta["report_summary"]["run_id"] == run_id
    assert analysis_meta["coverage_summary"]["interactive_ligands"] == 2
    assert analysis_meta["coverage_summary"]["interactive_targets"] == 2
    assert "validation_summary" in analysis_meta

    assert "Cell details" not in html_text
    assert "syncShortColLabels" in html_text
    assert "renderCellValueOverlay" in html_text
    assert "variant/pH: " in html_text
    assert "promiscuity_z: " in html_text
    assert "significant_hits: " in html_text
    assert "threshold: pct_rank <= " in html_text
    assert "Ligand organization" in html_text
    assert "Search SAR motifs" in html_text
    assert 'var ligandOrganizationMode = "none";' in html_text
    assert "function organizeRowsByMode(rowNames) {" in html_text
    assert "function appendLigandAnnotationLines(lines, meta) {" in html_text
    assert "function colLabelNodes()" in html_text
    assert "function rowLabelNodes()" in html_text
    assert "function lookupRowLabelMeta(label)" in html_text
    assert "appendPromiscuityLines(lines, meta);" in html_text
    assert "container.addEventListener(\"pointermove\"" in html_text
    assert "view: aggregated zoom" not in html_text
    assert 'shortLabel = shortLabel.replace(/\\.+$/, "");' in html_text
    assert 'label.removeAttribute("textLength");' in html_text
    assert "function bindColLabelObserver()" in html_text
    assert "setTimeout(scheduleColLabelSync, 1200);" in html_text
    assert "function bindDelegatedLabelKeyHandlers() {" in html_text
    assert 'container.setAttribute("data-cg-label-key-delegated", "1");' in html_text
    assert 'if (!evt || (evt.key !== "Enter" && evt.key !== " ")) return;' in html_text
    assert 'container.addEventListener("keydown", function(evt) {' in html_text
    assert 'label.setAttribute("tabindex", "0");' in html_text
    assert 'label.setAttribute("role", "button");' in html_text
    assert 'colLabel.setAttribute("aria-pressed", colActive ? "true" : "false");' in html_text
    assert 'rowLabel.setAttribute("aria-pressed", rowActive ? "true" : "false");' in html_text
    assert "function tooltipLineHtml(line) {" in html_text
    assert '<span class="tip-key">' in html_text
    assert '<span class="tip-val">' in html_text
    assert "function activatePdbSelection(label, pdb)" in html_text
    assert "function activateLigandSort(label)" in html_text
    assert "function resetHeatmapSelection()" in html_text
    assert 'activatePdbSelection(resolvedLabel, pdb);' in html_text
    assert 'activateLigandSort(rowName);' in html_text
    assert 'selectedColLabel = sameSelection ? "" : resolvedLabel;' not in html_text
    assert 'selectedLigand = selectedLigand === rowName ? "" : rowName;' not in html_text
    assert "function updateViewSummary(payload, visibleRowIndices, visibleColIndices)" in html_text
    assert "Target organization" in html_text
    assert "Safety" in html_text
    assert "ADME" in html_text
    assert "Safety Buckets" in html_text
    assert "Biology-first grouping with AE detail on hover/search" in html_text
    assert "ADME Buckets" in html_text
    assert "Protein family" in html_text
    assert "Pathway" in html_text
    assert 'var targetOrganizationMode = "none";' in html_text
    assert "function organizeColumnsByMode(colNames) {" in html_text
    assert "function renderOrganizationStrip(payload, visibleColIndices) {" in html_text
    assert "targetOrganizationMode !== \"none\"" in html_text
    assert "protein_family: " in html_text
    assert "adme: " in html_text
    assert "safety_primary: " in html_text
    assert "safety_secondary: " in html_text
    assert "safety_confidence: " in html_text
    assert "safety_sources: " in html_text
    assert "liability_examples: " in html_text
    assert "adverse_event_examples: " in html_text
    assert "Search adverse events" in html_text
    assert "pathways: " not in html_text
    assert "display block: " not in html_text
    assert "function organizationColorMap(mode) {" in html_text
    assert '#0072B2", "#D55E00", "#009E73"' in html_text
    assert "function organizationColorForValue(value) {" in html_text
    assert 'colLabel.style.fill = orgColor;' in html_text
    assert 'var networkColor = organizationColorForValue(' in html_text
    assert 'class="heatmap-org-block"' in html_text
    assert "Promiscuity vs Significant Hit Fraction" in html_text
    assert "Significant hit fraction (%)" in html_text
    assert "Promiscuity z-score" in html_text
    assert "Ranking by " in html_text
    assert "Promiscuity z-score" in html_text
    assert 'data-analysis-tab="network"' not in html_text
    assert 'data-analysis-view="network"' not in html_text
    assert "Promiscuity z-score" in html_text
    assert "Hist Ligands" not in html_text
    assert "Hist Targets" not in html_text
    assert "Histogram" not in html_text
    assert 'data-analysis-tab="rank"' in html_text
    assert 'data-scatter-mode="cols"' not in html_text
    assert 'data-rank-mode="cols"' not in html_text
    assert 'data-rank-metric="hit_fraction"' in html_text
    assert 'data-rank-metric="entropy"' not in html_text
    assert "Ligands shown:" not in html_text
    assert "Ligands Shown" not in html_text
    assert "window.location.href = url;" not in html_text
    assert 'window.open(url, "_blank")' not in html_text
    assert "Interactive Discovery Report" in html_text
    assert "Global search" in html_text
    assert "Copy shareable link" in html_text
    assert "Start Here" in html_text
    assert "Docking heatmap explorer for scientist-first discovery" in html_text
    assert 'data-clear-state="' in html_text
    assert "Most promiscuous ligands" in html_text
    assert "Most promiscuous targets" not in html_text
    assert "High-safety targets" not in html_text
    assert 'data-discovery-filter="' not in html_text
    assert "Search and rank the embedded ligand universe with discovery-specific controls." in html_text
    assert 'data-discovery-mode="cols"' not in html_text
    assert "data-discovery-focus-col" not in html_text
    assert '<div class="heatmap-report-panel-title">Discovery</div>' in html_text
    assert "Validation" in html_text
    assert "Analyst Workflows" in html_text
    assert "Export Current View CSV" in html_text
    assert "Export Ranking CSV" in html_text
    assert "Export Evidence CSV" in html_text
    assert "Save State JSON" in html_text
    assert "Load State JSON" in html_text
    assert "function refreshReportExtras(payload, visibleRowIndices, visibleColIndices)" in html_text
    assert "function syncReportStateHash()" in html_text
    assert "function hydrateReportStateFromLocation()" in html_text
    assert "function buildCurrentViewCsv()" in html_text
    assert "function buildRankingCsv()" in html_text
    assert "function buildEvidenceCsv()" in html_text
    assert "function renderHeroSummary()" in html_text
    assert "function renderHeroSearchResults()" in html_text
    assert "function renderActiveStateBar()" in html_text
    assert "function applyDiscoveryHighlights()" in html_text
    assert "function buildHeroSearchResults(query)" in html_text
    assert "function renderDiscoveryChips()" in html_text
    assert "function navigateAnalysisTab(tabName)" in html_text
    assert "colLabelTextNodes()" not in html_text
    assert "resolveColLabelFromNode(" not in html_text
    assert 'global_search: heroSearchEl ? String(heroSearchEl.value || "") : ""' in html_text
    assert "function renderValidationPanel(payload, visibleRowIndices, visibleColIndices)" in html_text
    assert "function renderDiscoveryPanel()" in html_text
    assert 'if (typeof refreshReportExtras === "function") refreshReportExtras(payload, visibleRowIndices, visibleColIndices);' in html_text
    assert 'if (typeof bindReportExtrasControls === "function") bindReportExtrasControls();' in html_text
    assert 'if (typeof hydrateReportStateFromLocation === "function") hydrateReportStateFromLocation();' in html_text


def test_heatmap_html_clustergrammer_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo_root = tmp_path / "repo"
    run_id = "CG_FALLBACK"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_pathway_memberships_by_pdb",
        _mock_pathway_memberships,
    )
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_target_uniprots_by_pdb",
        _mock_target_uniprots_by_pdb,
    )

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_minimal_rows())

    monkeypatch.setattr(heatmap_html, "_load_clustergrammer", lambda: None)
    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=5, include_decoys=False
    )
    assert ".hm-cell" in html_text


def test_heatmap_html_clustergrammer_compact_layout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("clustergrammer")
    repo_root = tmp_path / "repo"
    run_id = "CG_COMPACT"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_pathway_memberships_by_pdb",
        _mock_pathway_memberships,
    )
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_target_uniprots_by_pdb",
        _mock_target_uniprots_by_pdb,
    )

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_minimal_rows())

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=5, include_decoys=False
    )

    assert "row_label_scale: 0.88" in html_text
    assert "col_label_scale: 0.46" in html_text
    assert "cellValueMinSide = 18" in html_text
    assert 'class="cg-heatmap-wrap"' in html_text
    assert 'class="cg-scale-legend-rail"' in html_text
    assert 'class="cg-scale-legend-bar"' in html_text
    assert 'class="cg-heatmap-viz"' in html_text
    assert 'data-heatmap-controls="card"' in html_text
    assert 'class="heatmap-debug heatmap-debug-panel"' in html_text
    assert "new ResizeObserver(" in html_text
    assert "new MutationObserver(function() { scheduleColLabelSync(); });" in html_text
    assert "function renderCurrentMode()" in html_text
    assert "network_data: activeVizData()," in html_text
    assert 'var selectedColLabel = "";' in html_text
    assert 'var targetOrganizationMode = "none";' in html_text
    assert 'data-cg-mode="safety"' not in html_text
    assert 'function applyNativeSelectionReorder() {' in html_text
    assert 'new MouseEvent("dblclick"' in html_text
    assert 'args.row_order = "rank";' in html_text
    assert 'args.col_order = "rank";' in html_text
    assert 'args.row_order = "alpha";' not in html_text
    assert 'args.col_order = "alpha";' not in html_text
    assert 'args.row_order = "ini";' not in html_text
    assert 'args.col_order = "ini";' not in html_text
    assert 'node.rank = total - i;' in html_text
    assert 'node.rankvar = total - i;' in html_text
    assert "bindDelegatedLabelClickHandlers();" in html_text
    assert "input_domain: resolveInputDomain()," in html_text
    assert ".tile, .tile_up, .tile_dn" in html_text
    assert "function syncContainerWidth()" not in html_text
    assert "function refreshTileAppearance()" not in html_text
    assert "topkRangeEl" not in html_text
    assert "topkValueEl" not in html_text
    assert 'renderOrganizationStrip(payload, colIndices);' in html_text
    assert 'targetOrganizationMode = String((analysisMeta.target_organization && analysisMeta.target_organization.default_mode) || "none");' in html_text
    assert 'renderAnnotationPanel(safetyPanelEl, "Safety Buckets"' in html_text
    assert 'renderAnnotationPanel(admePanelEl, "ADME Buckets"' in html_text
    assert 'if (mode === "safety") return "Safety";' in html_text
    assert 'if (mode === "adme") return "ADME";' in html_text
    assert "function organizationColorMap(mode) {" in html_text
    assert 'colLabel.style.fill = orgColor;' in html_text
    assert 'var networkColor = organizationColorForValue(' in html_text
    assert 'circle.setAttribute("r", "5");' in html_text
    assert "drawAxisLabel(svg, pad.left + innerW / 2, height - 16, \"Significant hit fraction (%)\", 0);" in html_text
    assert "drawAxisLabel(svg, 18, pad.top + innerH / 2, \"Promiscuity z-score\", -90);" in html_text
    assert "function renderRankingPanel()" in html_text
    assert "Showing the top 20 items for readability." in html_text
    assert "Hit Count" not in html_text


def test_heatmap_html_clustergrammer_preserves_descending_sort_direction(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pytest.importorskip("clustergrammer")
    repo_root = tmp_path / "repo"
    run_id = "CG_SORTDIR"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_pathway_memberships_by_pdb",
        _mock_pathway_memberships,
    )
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_target_uniprots_by_pdb",
        _mock_target_uniprots_by_pdb,
    )

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_minimal_rows())

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=5, include_decoys=False
    )

    assert "function sortByCellScore(names, scorer, tieOrder) {" in html_text
    assert "if (bx !== ax) return bx - ax;" in html_text
    assert "orderedNames = sortByCellScore(defaultRankedRows(), function(rowName) {" in html_text
    assert (
        "return rawScoreValue(lookupMeta(rowName, selectedColLabel));" in html_text
    )
    assert "orderedNames = sortByCellScore(colNames, function(colName) {" in html_text
    assert "return rawScoreValue(lookupMeta(selectedLigand, colName));" in html_text
    assert "node.rank = total - i;" in html_text
    assert "node.rankvar = total - i;" in html_text
    assert "node.rank = i + 1;" not in html_text
    assert "node.rankvar = i + 1;" not in html_text


def test_clustergrammer_mode_matrix_missing_values_use_neutral_fill() -> None:
    out = _build_clustergrammer_mode_mat(
        matrix=[[1.25, None], [-0.5, 2.0]],
        row_labels=["L1", "L2"],
        col_labels=["T1", "T2"],
        ordered_row_labels=["L1", "L2"],
        ordered_col_labels=["T1", "T2"],
        breaks=None,
    )

    assert out == [[1.25, 0.0], [-0.5, 2.0]]


def test_build_promiscuity_meta_uses_threshold_counts_and_entropy() -> None:
    agg = {
        ("LigA", "T1"): {"pct_rank": "0.005", "pose_valid_any": "true", "z_selected": 3.0},
        ("LigA", "T2"): {"pct_rank": "0.009", "pose_valid_any": "", "z_selected": 1.5},
        ("LigB", "T1"): {"pct_rank": "0.2", "pose_valid_any": "true", "z_selected": 2.0},
        ("LigB", "T2"): {"pct_rank": "0.004", "pose_valid_any": "false", "z_selected": 4.0},
    }

    ligand_meta, target_meta = _build_promiscuity_meta(agg, pct_threshold=0.01)

    assert ligand_meta["LigA"]["significant_hits"] == pytest.approx(2.0)
    assert ligand_meta["LigA"]["partner_total"] == pytest.approx(2.0)
    assert ligand_meta["LigA"]["hit_fraction"] == pytest.approx(1.0)
    assert ligand_meta["LigA"]["entropy"] > 0.0
    assert ligand_meta["LigB"]["significant_hits"] == pytest.approx(0.0)
    assert ligand_meta["LigB"]["promiscuity_z"] < ligand_meta["LigA"]["promiscuity_z"]

    assert target_meta["T1"]["significant_hits"] == pytest.approx(1.0)
    assert target_meta["T1"]["partner_total"] == pytest.approx(2.0)
    assert target_meta["T1"]["hit_fraction"] == pytest.approx(0.5)
    assert target_meta["T2"]["significant_hits"] == pytest.approx(1.0)


def test_heatmap_html_preserves_targets_without_top_ligand_hits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("clustergrammer")
    repo_root = tmp_path / "repo"
    run_id = "CG_TARGETS"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_pathway_memberships_by_pdb",
        _mock_pathway_memberships,
    )
    monkeypatch.setattr(
        heatmap_html,
        "_fetch_target_uniprots_by_pdb",
        _mock_target_uniprots_by_pdb,
    )

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(
        input_csv,
        [
            {
                "target_id": "T1",
                "pdb_id": "1AAA",
                "target_name": "Alpha",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "ligand_display": "L1",
                "ligand_base": "L1",
                "z_selected": "5.0",
                "rank": "1",
                "pct_rank": "0.1",
                "pose_valid_any": "true",
                "library": "lib_a",
            },
            {
                "target_id": "T2",
                "pdb_id": "2BBB",
                "target_name": "Beta",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "ligand_display": "L2",
                "ligand_base": "L2",
                "z_selected": "1.0",
                "rank": "1",
                "pct_rank": "0.1",
                "pose_valid_any": "true",
                "library": "lib_a",
            },
        ],
    )

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=1, include_decoys=False
    )

    viz_match = re.search(
        r'<script type="application/json" id="cg-heatmap-data-CG_TARGETS">(.*?)</script>',
        html_text,
        re.S,
    )
    assert viz_match is not None
    viz_payloads = json.loads(viz_match.group(1))
    col_nodes = viz_payloads["fixed"]["col_nodes"]
    assert [node["name"] for node in col_nodes] == [
        "1AAA\nAlpha\nHOLO | pH7_0",
        "2BBB\nBeta\nHOLO | pH7_0",
    ]
