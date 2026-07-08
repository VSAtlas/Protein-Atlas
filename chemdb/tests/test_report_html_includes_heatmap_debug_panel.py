import csv
from pathlib import Path

from analysis.reporting import heatmap_html as heatmap_html
from analysis.reporting.heatmap_html import render_interactive_heatmap_html


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def test_report_html_includes_heatmap_debug_panel(tmp_path: Path, monkeypatch) -> None:
    repo_root = tmp_path / "repo"
    run_id = "HEATMAP_DEBUG"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text(
        "\n".join(
            [
                "USE_DENDROGRAM=false",
                "REPORT_OFFLINE_ASSETS=false",
                "REPORT_INLINE_ASSETS=false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(
        input_csv,
        [
            {
                "target_id": "T1",
                "ligand_display": "L1",
                "ligand_base": "L1",
                "z_selected": "1.0",
                "rank": "1",
                "pct_rank": "0.1",
                "pose_valid_any": "1",
                "library": "LIB",
            }
        ],
    )

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

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=5, include_decoys=False
    )

    assert "id=\"heatmap-debug\"" in html_text
    assert "heatmap-debug-log" in html_text
    assert "window.onerror" in html_text
    assert "unhandledrejection" in html_text
    assert "performance.now" in html_text
    assert "clustergrammer_init_ms" in html_text
