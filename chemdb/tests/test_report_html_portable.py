import csv
import re
import shutil
from pathlib import Path

from analysis import heatmap_html, run_report


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
            "t_selected": "1.0",
            "rank": "1",
            "pct_rank": "0.1",
            "pose_valid_any": "1",
            "library": "LIB",
        },
        {
            "target_id": "T2|V1|P1",
            "ligand_display": "L1",
            "t_selected": "1.5",
            "rank": "2",
            "pct_rank": "0.2",
            "pose_valid_any": "0",
            "library": "LIB",
        },
        {
            "target_id": "T1|V1|P1",
            "ligand_display": "L2",
            "t_selected": "2.0",
            "rank": "3",
            "pct_rank": "0.3",
            "pose_valid_any": "1",
            "library": "LIB",
        },
        {
            "target_id": "T2|V1|P1",
            "ligand_display": "L2",
            "t_selected": "2.5",
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
                        "t_selected": 2.5,
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
                "best_t_selected": 2.5,
            }
        ],
    }


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

    source_root = Path(__file__).resolve().parents[2]
    assets_src = source_root / "report_assets" / "clustergrammer"
    assets_dst = repo_root / "report_assets" / "clustergrammer"
    shutil.copytree(assets_src, assets_dst)

    monkeypatch.setattr(
        heatmap_html,
        "_build_clustergrammer_viz_json",
        lambda *_args, **_kwargs: {
            "row_nodes": [{"name": "L1"}],
            "col_nodes": [{"name": "T1"}],
            "mat": [[1.0]],
            "links": [],
            "views": [],
            "cat_colors": {},
        },
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        _make_report(run_id),
        out_path,
        ["q1"],
        repo_root,
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "report_assets/" not in html_text
    assert "clustergrammer v1.19.5" in html_text
    assert "<style>" in html_text
    assert "<script>" in html_text


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
    assert "href=\"#multitarget\"" in html_text
    assert "href=\"#heatmap\"" in html_text
    assert "href=\"#artifacts\"" in html_text
    assert "id=\"highlights\"" in html_text
    assert "id=\"multitarget\"" in html_text
    assert "id=\"heatmap\"" in html_text
    assert "id=\"artifacts\"" in html_text
    assert "tbody tr:nth-child(even)" in html_text
    assert "tbody tr:hover" in html_text
    assert "position: sticky" in html_text
