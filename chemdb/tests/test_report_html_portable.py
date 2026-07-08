import csv
import json
import re
import shutil
from pathlib import Path

from analysis.reporting import heatmap_html, run_report_core as run_report


def _write_config(repo_root: Path) -> None:
    config = repo_root / "config.txt"
    config.write_text(
        "\n".join(
            [
                "HEATMAP_COLOR_MIN=#1b9e77",
                "HEATMAP_COLOR_MID=#d95f02",
                "HEATMAP_COLOR_MID2=#7570b3",
                "HEATMAP_COLOR_MAX=#e7298a",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_heatmap_input(path: Path) -> None:
    rows = [
        {
            "target_id": "T1|V1|P1",
            "ligand_display": "L1",
            "z_selected": "1.0",
            "rank": "1",
            "pct_rank": "0.1",
            "pose_valid_any": "1",
            "library": "LIB",
        },
        {
            "target_id": "T2|V1|P1",
            "ligand_display": "L1",
            "z_selected": "1.5",
            "rank": "2",
            "pct_rank": "0.2",
            "pose_valid_any": "0",
            "library": "LIB",
        },
        {
            "target_id": "T1|V1|P1",
            "ligand_display": "L2",
            "z_selected": "2.0",
            "rank": "3",
            "pct_rank": "0.3",
            "pose_valid_any": "1",
            "library": "LIB",
        },
        {
            "target_id": "T2|V1|P1",
            "ligand_display": "L2",
            "z_selected": "2.5",
            "rank": "4",
            "pct_rank": "0.4",
            "pose_valid_any": "1",
            "library": "LIB",
        },
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _make_report(run_id: str) -> dict:
    return {
        "run_id": run_id,
        "generated_at": "2025-01-01T00:00:00Z",
        "targets": {
            "T1|V1|P1": {
                "highlights": [
                    {
                        "query": "q1",
                        "found": True,
                        "rank": 1,
                        "pct_rank": 0.001,
                        "z_selected": 2.5,
                    }
                ]
            }
        },
        "multi_target_hits": [
            {
                "ligand_display": "L1",
                "ligand_base": "L1",
                "targets_qualified": 2,
                "worst_pct": 0.2,
                "mean_pct": 0.1,
                "best_rank": 1,
                "best_z_selected": 2.5,
            }
        ],
    }


def _stub_clustergrammer_viz(*_args, **_kwargs):
    return {
        "row_nodes": [{"name": "L1"}, {"name": "L2"}],
        "col_nodes": [{"name": "T1|V1|P1"}, {"name": "T2|V1|P1"}],
        "mat": [[1.0, 1.5], [2.0, 2.5]],
        "links": [],
        "views": [],
        "cat_colors": {},
    }


def _resolve_clustergrammer_assets_root(source_root: Path) -> Path | None:
    direct = source_root / "report_assets" / "clustergrammer"
    if direct.is_dir():
        return direct

    pointer = source_root / "report_assets"
    if pointer.is_file():
        try:
            resolved = Path(pointer.read_text(encoding="utf-8").strip())
        except Exception:
            return None
        candidate = resolved / "clustergrammer"
        if candidate.is_dir():
            return candidate
    return None


def test_report_html_stages_artifacts_and_relative_links(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "PORTABLE_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)
    (data_dir / "report.yaml").write_text("report: ok\n", encoding="utf-8")
    (data_dir / "heatmap_input.csv").write_text("placeholder\n", encoding="utf-8")
    (data_dir / "heatmap.png").write_bytes(b"PNG")

    monkeypatch.setattr(
        run_report,
        "render_interactive_heatmap_html",
        lambda *_args, **_kwargs: "<div>heatmap</div>",
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
    )

    html_text = out_path.read_text(encoding="utf-8")
    artifacts_dir = data_dir / "artifacts"
    assert (artifacts_dir / "report.yaml").exists()
    assert (artifacts_dir / "heatmap_input.csv").exists()
    assert (artifacts_dir / "heatmap.png").exists()
    assert "artifacts/report.yaml" in html_text
    assert "artifacts/heatmap_input.csv" in html_text
    assert "artifacts/heatmap.png" in html_text
    assert "/home/" not in html_text
    assert re.search(r"[A-Za-z]:\\\\", html_text) is None


def test_report_html_inline_assets_default(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "INLINE_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)
    _write_heatmap_input(data_dir / "heatmap_input.csv")

    monkeypatch.setattr(
        heatmap_html,
        "_build_clustergrammer_viz_json",
        _stub_clustergrammer_viz,
    )
    monkeypatch.setattr(
        heatmap_html,
        "_clustergrammer_inline_assets",
        lambda *_args, **_kwargs: "<style>/* inline */</style><script>/* inline */</script>",
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
        report_asset_mode="inline",
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "/* inline */" in html_text
    assert "<style>" in html_text
    assert "<script>" in html_text
    assert "asset_mode=inline" in html_text


def test_report_html_inline_mode_falls_back_to_cdn_without_local_assets(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "INLINE_FALLBACK_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)
    _write_heatmap_input(data_dir / "heatmap_input.csv")

    monkeypatch.setattr(
        heatmap_html,
        "_build_clustergrammer_viz_json",
        _stub_clustergrammer_viz,
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
        report_asset_mode="inline",
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "Asset mode: CDN fallback." in html_text
    assert "inline_assets_failed" in html_text
    assert "https://unpkg.com/clustergrammer@1.19.5/clustergrammer.js" in html_text
    assert "asset_mode=cdn" in html_text


def test_report_html_relative_assets_mode(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "RELATIVE_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)
    _write_heatmap_input(data_dir / "heatmap_input.csv")

    monkeypatch.setattr(
        heatmap_html,
        "_build_clustergrammer_viz_json",
        _stub_clustergrammer_viz,
    )
    monkeypatch.setattr(
        heatmap_html,
        "_clustergrammer_asset_urls",
        lambda *_args, **_kwargs: {
            "d3": "../../report_assets/clustergrammer/lib/js/d3.js",
            "jquery": "../../report_assets/clustergrammer/lib/js/jquery-1.11.2.min.js",
            "underscore": "../../report_assets/clustergrammer/lib/js/underscore-min.js",
            "bootstrap_js": "../../report_assets/clustergrammer/lib/js/bootstrap.min.js",
            "bootstrap_css": "../../report_assets/clustergrammer/lib/css/bootstrap.css",
            "font_awesome_css": "../../report_assets/clustergrammer/lib/css/font-awesome.min.css",
            "custom_css": "../../report_assets/clustergrammer/css/custom.css",
            "clustergrammer_js": "../../report_assets/clustergrammer/clustergrammer.js",
        },
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
        report_asset_mode="relative",
    )
    html_text = out_path.read_text(encoding="utf-8")
    assert "../../report_assets/clustergrammer/" in html_text
    assert "asset_mode=relative" in html_text


def test_report_html_cdn_assets_mode(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "CDN_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)
    _write_heatmap_input(data_dir / "heatmap_input.csv")

    monkeypatch.setattr(
        heatmap_html,
        "_build_clustergrammer_viz_json",
        _stub_clustergrammer_viz,
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
        report_asset_mode="cdn",
    )
    html_text = out_path.read_text(encoding="utf-8")
    assert "https://unpkg.com/clustergrammer@1.19.5/clustergrammer.js" in html_text
    assert "asset_mode=cdn" in html_text


def test_report_html_toc_and_table_styling(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "STYLE_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)

    monkeypatch.setattr(
        run_report,
        "render_interactive_heatmap_html",
        lambda *_args, **_kwargs: "<div>heatmap</div>",
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "href=\"#highlights\"" in html_text
    if 'id="multitarget"' in html_text:
        assert "href=\"#multitarget\"" in html_text
    assert "href=\"#heatmap\"" in html_text
    assert "href=\"#artifacts\"" in html_text
    assert "href=\"#technical\"" in html_text
    assert "id=\"highlights\"" in html_text
    assert "id=\"heatmap\"" in html_text
    assert "id=\"artifacts\"" in html_text
    assert "id=\"technical\"" in html_text
    assert "Atlas2 Report: STYLE_RUN" in html_text
    assert "Top Ligand-Target Signals" in html_text
    assert "Interactive Docking Heatmap" in html_text
    assert "Run Artifacts" in html_text
    assert 'class="report-kpis"' in html_text
    assert 'class="limitations-note"' in html_text
    assert "Interpretation and Limitations" in html_text
    assert "computational hypotheses" in html_text
    assert "position: sticky" in html_text
    assert ".highlights-table td:first-child" in html_text
    assert ".top-targets {" in html_text
    assert "@media print" in html_text


def test_report_html_inline_assets_from_pointer_file(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "POINTER_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)
    _write_heatmap_input(data_dir / "heatmap_input.csv")

    source_root = Path(__file__).resolve().parents[2]
    assets_src = _resolve_clustergrammer_assets_root(source_root)
    if assets_src is None:
        return

    ext_root = tmp_path / "external_assets"
    ext_clustergrammer = ext_root / "clustergrammer"
    shutil.copytree(assets_src, ext_clustergrammer)
    (repo_root / "report_assets").write_text(str(ext_root), encoding="utf-8")

    monkeypatch.setattr(
        heatmap_html,
        "_build_clustergrammer_viz_json",
        _stub_clustergrammer_viz,
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "../../report_assets/clustergrammer/" not in html_text
    assert "clustergrammer v1.19.5" in html_text
    assert "<style>" in html_text
    assert "<script>" in html_text


def test_report_html_dead_pointer_emits_cdn_asset_retry_metadata(tmp_path, monkeypatch):
    repo_root = tmp_path / "repo"
    run_id = "BROKEN_POINTER_RUN"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    _write_config(repo_root)
    _write_heatmap_input(data_dir / "heatmap_input.csv")
    (repo_root / "report_assets").write_text(
        str(tmp_path / "missing_assets_root"), encoding="utf-8"
    )

    monkeypatch.setattr(
        heatmap_html,
        "_build_clustergrammer_viz_json",
        _stub_clustergrammer_viz,
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "Asset mode: CDN fallback." in html_text
    assert "https://cdn.jsdelivr.net/npm/clustergrammer@1.19.5/clustergrammer.js" in html_text
    asset_match = re.search(
        r'<script type="application/json" id="cg-heatmap-assets-BROKEN_POINTER_RUN">(.*?)</script>',
        html_text,
        re.S,
    )
    assert asset_match is not None
    asset_meta = json.loads(asset_match.group(1))
    assert asset_meta["mode"] == "cdn"
    assert asset_meta["notices"]
    assert asset_meta["asset_urls"]["clustergrammer_js"][0].endswith("clustergrammer.js")
    assert "cdn.jsdelivr.net" in " ".join(asset_meta["asset_urls"]["clustergrammer_js"])
