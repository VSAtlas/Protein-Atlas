import csv
import json
import os
from pathlib import Path

import pytest

from analysis.reporting import heatmap_html
from analysis.reporting.heatmap_html import render_interactive_heatmap_html


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _mock_pathway_memberships(
    _repo_root: Path, _pdb_ids: object, **_kwargs: object
) -> dict[str, list[str]]:
    return {
        "1ABC": ["Cell Cycle"],
        "2XYZ": ["Cardiac conduction"],
    }


def _mock_target_uniprots_by_pdb(
    _repo_root: Path, _pdb_ids: object
) -> dict[str, list[str]]:
    return {
        "1ABC": ["P03372"],
        "2XYZ": ["Q12809"],
    }


def _build_rows() -> list[dict[str, str]]:
    return [
        {
            "target_id": "T1",
            "pdb_id": "1ABC",
            "target_name": "Progesterone receptor",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "ligand_display": "LigA",
            "ligand_base": "LigA",
            "z_selected": "1.5",
            "rank": "1",
            "pct_rank": "0.005",
            "pose_valid_any": "true",
            "library": "lib_a",
        },
        {
            "target_id": "T2",
            "pdb_id": "2XYZ",
            "target_name": "Cardiac receptor",
            "variant": "APO",
            "ph_label": "pH7_4",
            "ligand_display": "LigA",
            "ligand_base": "LigA",
            "z_selected": "2.0",
            "rank": "2",
            "pct_rank": "0.004",
            "pose_valid_any": "true",
            "library": "lib_a",
        },
        {
            "target_id": "T1",
            "pdb_id": "1ABC",
            "target_name": "Progesterone receptor",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "ligand_display": "LigB",
            "ligand_base": "LigB",
            "z_selected": "0.5",
            "rank": "3",
            "pct_rank": "0.03",
            "pose_valid_any": "true",
            "library": "lib_b",
        },
        {
            "target_id": "T2",
            "pdb_id": "2XYZ",
            "target_name": "Cardiac receptor",
            "variant": "APO",
            "ph_label": "pH7_4",
            "ligand_display": "LigB",
            "ligand_base": "LigB",
            "z_selected": "3.5",
            "rank": "4",
            "pct_rank": "0.006",
            "pose_valid_any": "true",
            "library": "lib_b",
        },
    ]


@pytest.mark.skipif(
    os.environ.get("RUN_PLAYWRIGHT_BROWSER_TESTS") != "1",
    reason="set RUN_PLAYWRIGHT_BROWSER_TESTS=1 to enable browser-level heatmap regression tests",
)
@pytest.mark.skipif(
    not Path("/stor/home/mpg2352/playwright-browsers").exists(),
    reason="local playwright browsers not installed",
)
def test_heatmap_browser_regressions(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("clustergrammer")
    playwright = pytest.importorskip("playwright.sync_api")
    repo_root = tmp_path / "repo"
    run_id = "PW_CG"
    data_dir = repo_root / "data" / run_id
    cache_dir = repo_root / "pathways" / "cache"
    data_dir.mkdir(parents=True)
    cache_dir.mkdir(parents=True)
    current_repo_root = Path(__file__).resolve().parents[2]
    (repo_root / "report_assets").write_text(
        str(current_repo_root / "report_assets"),
        encoding="utf-8",
    )
    (repo_root / "config.txt").write_text(
        "USE_DENDROGRAM=false\nREPORT_INLINE_ASSETS=true\n",
        encoding="utf-8",
    )
    (cache_dir / "target_safety_aggregated.json").write_text(
        json.dumps(
            {
                "generated_at": "2026-03-10T12:00:00+00:00",
                "entries": [
                    {
                        "pdb_id": "1ABC",
                        "target_name": "Progesterone receptor",
                        "uniprots": ["P03372"],
                        "primary_display_safety": "Endocrine",
                        "secondary_safety_buckets": ["GI"],
                        "direct_safety_buckets": ["Endocrine"],
                        "drug_ae_buckets": ["GI"],
                        "safety_buckets": ["Endocrine", "GI"],
                        "safety_confidence": "high",
                        "safety_sources": ["Open Targets safety"],
                        "direct_liability_examples": ["progesterone signaling"],
                        "raw_adverse_event_examples": ["nausea"],
                    },
                    {
                        "pdb_id": "2XYZ",
                        "target_name": "Cardiac receptor",
                        "uniprots": ["Q12809"],
                        "primary_display_safety": "Cardiotoxicity / QT",
                        "secondary_safety_buckets": ["Hepatotoxicity"],
                        "direct_safety_buckets": ["Cardiotoxicity / QT"],
                        "drug_ae_buckets": ["Hepatotoxicity"],
                        "safety_buckets": ["Cardiotoxicity / QT", "Hepatotoxicity"],
                        "safety_confidence": "high",
                        "safety_sources": ["Open Targets safety", "Open Targets drug adverse events"],
                        "direct_liability_examples": ["qt prolongation"],
                        "raw_adverse_event_examples": ["qt prolongation", "ventricular tachycardia"],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    (cache_dir / "target_safety_drift_summary.json").write_text(
        json.dumps({"generated_at": "2026-03-11T12:00:00+00:00", "changed_count": 1}),
        encoding="utf-8",
    )
    _write_csv(
        data_dir / "master_rows.csv",
        _build_rows()
        + [
            {
                "target_id": "T3",
                "pdb_id": "3HIJ",
                "target_name": "Hidden target",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "ligand_display": "LigC",
                "ligand_base": "LigC",
                "z_selected": "0.1",
                "rank": "5",
                "pct_rank": "0.2",
                "pose_valid_any": "true",
                "library": "lib_c",
            }
        ],
    )
    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_rows())
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

    html_path = data_dir / "report.html"
    html_path.write_text(
        render_interactive_heatmap_html(repo_root, run_id, input_csv, top_k=2, include_decoys=False),
        encoding="utf-8",
    )
    container_id = f"cg-heatmap-{run_id}"
    summary_id = f"cg-heatmap-view-summary-{run_id}"

    with playwright.sync_playwright() as pw:
        local_tmp = Path("/stor/home/mpg2352/playwright_tmp")
        local_tmp.mkdir(parents=True, exist_ok=True)
        browser = pw.chromium.launch(
            headless=True,
            env={
                **os.environ,
                "TMPDIR": str(local_tmp),
                "TMP": str(local_tmp),
                "TEMP": str(local_tmp),
            },
        )
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()
        page_errors: list[str] = []
        page.on("pageerror", lambda err: page_errors.append(str(err)))
        page.set_default_timeout(15000)
        page.goto(html_path.as_uri(), wait_until="domcontentloaded")
        page.wait_for_function(
            """(containerId) => {
                const labels = document.querySelectorAll('#' + containerId + ' .col_label_text text, #' + containerId + ' text.col_label_text');
                return labels.length >= 2;
            }""",
            container_id,
        )
        assert not page_errors

        page.click(f"#cg-heatmap-ligand-org-chemotype-{run_id}")
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('Ligand organization: Chemotype')""",
            summary_id,
        )
        page.fill(f"#cg-heatmap-motif-search-{run_id}", "lig")
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('Motif search: lig')""",
            summary_id,
        )
        motif_hits = page.evaluate(
            """(containerId) => document.querySelectorAll('#' + containerId + ' .cg-motif-search-hit').length""",
            container_id,
        )
        assert motif_hits >= 1
        page.click(f"#cg-heatmap-motif-search-clear-{run_id}")
        page.wait_for_function(
            """(summaryId) => !(document.getElementById(summaryId)?.innerText || '').includes('Motif search:')""",
            summary_id,
        )

        page.fill(f"#cg-heatmap-hero-search-{run_id}", "progesterone")
        page.wait_for_selector(f"#cg-heatmap-hero-results-{run_id} [data-search-target]")
        page.click(f"#cg-heatmap-hero-results-{run_id} [data-search-target]")
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('PDB selected: 1ABC')""",
            summary_id,
        )
        page.click(f"#cg-heatmap-hero-reset-{run_id}")
        page.wait_for_function(
            """(summaryId) => {
                const text = document.getElementById(summaryId)?.innerText || '';
                return !text.includes('PDB selected:') && !text.includes('Ligand sort:');
            }""",
            summary_id,
        )
        page.evaluate(
            """(containerId) => {
                const row = Array.from(document.querySelectorAll('#' + containerId + ' .row_label_group'))
                  .find((el) => (el.textContent || '').trim() === 'LigA');
                row.setAttribute('tabindex', '0');
                row.focus();
                row.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true, cancelable: true }));
            }""",
            container_id,
        )
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('Ligand sort: LigA')""",
            summary_id,
        )
        page.click(f"#cg-heatmap-hero-reset-{run_id}")
        page.wait_for_function(
            """(summaryId) => !(document.getElementById(summaryId)?.innerText || '').includes('Ligand sort:')""",
            summary_id,
        )

        initial_cols = page.evaluate(
            """(containerId) => Array.from(document.querySelectorAll('#' + containerId + ' .col_label_text text, #' + containerId + ' text.col_label_text')).map((n) => n.textContent.trim()).filter(Boolean)""",
            container_id,
        )
        assert initial_cols[:2] == ["1ABC", "2XYZ"]

        page.evaluate(
            """(containerId, rowName) => {
                const node = Array.from(document.querySelectorAll('#' + containerId + ' .row_label_group')).find((el) => (el.textContent || '').trim() === rowName);
                node.dispatchEvent(new MouseEvent('click', { bubbles: true, cancelable: true, view: window }));
            }""",
            container_id,
            "LigA",
        )
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('Ligand sort: LigA')""",
            summary_id,
        )
        ligand_sorted_cols = page.evaluate(
            """(containerId) => Array.from(document.querySelectorAll('#' + containerId + ' .col_label_text text, #' + containerId + ' text.col_label_text')).map((n) => n.textContent.trim()).filter(Boolean)""",
            container_id,
        )
        assert ligand_sorted_cols[0] == "2XYZ"

        page.evaluate(
            """(containerId, colName) => {
                const node = Array.from(document.querySelectorAll('#' + containerId + ' .col_label_text text, #' + containerId + ' text.col_label_text')).find((el) => (el.textContent || '').trim() === colName);
                node.dispatchEvent(new MouseEvent('click', { bubbles: true, cancelable: true, view: window }));
            }""",
            container_id,
            "2XYZ",
        )
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('PDB selected: 2XYZ')""",
            summary_id,
        )
        focused_cols = page.evaluate(
            """(containerId) => Array.from(document.querySelectorAll('#' + containerId + ' .col_label_text text, #' + containerId + ' text.col_label_text')).map((n) => n.textContent.trim()).filter(Boolean)""",
            container_id,
        )
        assert focused_cols == ["2XYZ"]
        focused_rows = page.evaluate(
            """(containerId) => Array.from(document.querySelectorAll('#' + containerId + ' .row_label_group')).map((n) => (n.textContent || '').trim()).filter(Boolean)""",
            container_id,
        )
        assert focused_rows[0] == "LigB"

        page.click(f"#cg-heatmap-org-safety-{run_id}")
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('Target organization: Safety')""",
            summary_id,
        )
        org_strip_text = page.locator(f"#cg-heatmap-org-strip-{run_id}").inner_text()
        assert "Cardiotoxicity / QT" in org_strip_text

        tooltip_text = page.evaluate(
            """(containerId, label) => {
                const node = Array.from(document.querySelectorAll('#' + containerId + ' .col_label_text text, #' + containerId + ' text.col_label_text')).find((el) => (el.textContent || '').trim() === label);
                const rect = node.getBoundingClientRect();
                node.dispatchEvent(new PointerEvent('pointerenter', { bubbles: true, clientX: rect.left + 3, clientY: rect.top + 3 }));
                node.dispatchEvent(new PointerEvent('pointermove', { bubbles: true, clientX: rect.left + 3, clientY: rect.top + 3 }));
                return document.getElementById('cg-hover-tip')?.innerText || '';
            }""",
            container_id,
            "2XYZ",
        )
        assert "safety_primary: Cardiotoxicity / QT" in tooltip_text
        assert "safety_secondary: Hepatotoxicity" in tooltip_text
        assert "safety_confidence: high" in tooltip_text

        page.fill(f"#cg-heatmap-safety-search-{run_id}", "qt")
        page.wait_for_function(
            """(summaryId) => (document.getElementById(summaryId)?.innerText || '').includes('AE search: qt')""",
            summary_id,
        )
        hit_count = page.evaluate(
            """(containerId) => document.querySelectorAll('#' + containerId + ' .cg-safety-search-hit').length""",
            container_id,
        )
        assert hit_count >= 1
        assert "Target organization: Safety" in page.locator(f"#{summary_id}").inner_text()

        current_hash = page.evaluate("window.location.hash")
        assert "target_organization=safety" in current_hash
        assert "safety_search=qt" in current_hash
        assert "global_search=" not in current_hash

        with page.expect_download() as download_info:
            page.click(f"#{'cg-heatmap-export-view-' + run_id}")
        download = download_info.value
        download_path = tmp_path / "current_view.csv"
        download.save_as(str(download_path))
        csv_text = download_path.read_text(encoding="utf-8")
        assert "# heatmap_state," in csv_text
        assert "LigB" in csv_text
        assert "2XYZ" in csv_text

        page.reload()
        page.wait_for_function(
            """(summaryId) => {
                const text = document.getElementById(summaryId)?.innerText || '';
                return text.includes('Target organization: Safety') && text.includes('AE search: qt');
            }""",
            summary_id,
        )
        assert not page_errors

        browser.close()
