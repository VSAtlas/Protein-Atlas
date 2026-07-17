from analysis.reporting import run_report_core as run_report


def test_report_html_includes_scientist_interpretation_guide(tmp_path, monkeypatch) -> None:
    repo_root = tmp_path / "repo"
    run_id = "GUIDE_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (data_dir / "report.yaml").write_text("report: ok\n", encoding="utf-8")
    (data_dir / "heatmap_input.csv").write_text("placeholder\n", encoding="utf-8")

    monkeypatch.setattr(
        run_report,
        "render_interactive_heatmap_html",
        lambda *_args, **_kwargs: "<div>heatmap</div>",
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
    assert "How to interpret this report" in html_text
    assert "<dt>Atlas score</dt>" in html_text
    assert "Atlas score source" in html_text
    assert "Percentile rank" in html_text
    assert "Coverage" in html_text
    assert "Breadth 1%" in html_text
    assert "Target organization" in html_text
    assert "significant hit threshold" in html_text
    assert "Pose validity" in html_text
    assert "MM/GBSA" in html_text
    assert "not validated binding" in html_text
