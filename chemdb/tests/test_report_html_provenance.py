import re

from analysis.reporting import run_report_core as run_report


def test_report_html_provenance_uses_relative_links(tmp_path, monkeypatch) -> None:
    repo_root = tmp_path / "repo"
    run_id = "PROV_RUN"
    data_dir = repo_root / "data" / run_id
    manifest_dir = repo_root / "manifests" / run_id
    data_dir.mkdir(parents=True)
    manifest_dir.mkdir(parents=True)
    (data_dir / "report.yaml").write_text("report: ok\n", encoding="utf-8")
    (data_dir / "heatmap_input.csv").write_text("placeholder\n", encoding="utf-8")
    (data_dir / "config_snapshot.txt").write_text("OVERALL_DIR=/redacted\n", encoding="utf-8")
    (manifest_dir / "run_manifest.yaml").write_text("run_id: PROV_RUN\n", encoding="utf-8")

    monkeypatch.setenv("GITHUB_SHA", "abcdef1234567890")
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
            "atlas_version": "test-version",
            "targets": {},
            "summary": {"n_unique_ligands": 0, "n_target_combos": 0},
        },
        out_path,
        [],
        repo_root,
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "Provenance" in html_text
    assert "test-version" in html_text
    assert "abcdef123456" in html_text
    assert "artifacts/config_snapshot.txt" in html_text
    assert "artifacts/run_manifest.yaml" in html_text
    assert "artifacts/heatmap_input.csv" in html_text
    assert "/home/" not in html_text
    assert str(tmp_path) not in html_text
    assert re.search(r"[A-Za-z]:\\\\", html_text) is None
