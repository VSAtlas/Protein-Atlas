from pathlib import Path

from analysis.reporting import run_report_core as run_report


def test_report_html_has_accessibility_landmarks(tmp_path, monkeypatch) -> None:
    repo_root = tmp_path / "repo"
    run_id = "ACCESS_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (data_dir / "report.yaml").write_text("report: ok\n", encoding="utf-8")
    (data_dir / "heatmap_input.csv").write_text("placeholder\n", encoding="utf-8")

    monkeypatch.setattr(
        run_report,
        "render_interactive_heatmap_html",
        lambda *_args, **_kwargs: '<button type="button">Reset</button>',
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        {
            "run_id": run_id,
            "generated_at": "2026-05-02T00:00:00Z",
            "targets": {},
            "summary": {"n_unique_ligands": 0, "n_target_combos": 0},
        },
        out_path,
        [],
        repo_root,
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert '<main class="container" id="report-content">' in html_text
    assert '<header class="page-header">' in html_text
    assert '<nav class="toc" aria-label="Report sections">' in html_text
    assert 'class="skip-link" href="#heatmap"' in html_text
    assert 'id="heatmap" aria-labelledby="heatmap-heading" aria-live="polite"' in html_text
    assert 'id="heatmap-heading"' in html_text
    assert ":focus-visible" in html_text
    assert "@media (prefers-reduced-motion: reduce)" in html_text
