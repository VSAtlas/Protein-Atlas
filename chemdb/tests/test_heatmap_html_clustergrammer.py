import csv
import json
import re
from pathlib import Path

import pytest

import analysis.heatmap_html as heatmap_html
from analysis.heatmap_html import render_interactive_heatmap_html


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
            "ligand_display": "L1",
            "ligand_base": "L1",
            "t_selected": "1.5",
            "rank": "1",
            "pct_rank": "0.1",
            "pose_valid_any": "true",
            "library": "lib_a",
        },
        {
            "target_id": "T2",
            "ligand_display": "L1",
            "ligand_base": "L1",
            "t_selected": "2.0",
            "rank": "2",
            "pct_rank": "0.2",
            "pose_valid_any": "false",
            "library": "lib_a",
        },
        {
            "target_id": "T1",
            "ligand_display": "L2",
            "ligand_base": "L2",
            "t_selected": "0.5",
            "rank": "3",
            "pct_rank": "0.3",
            "pose_valid_any": "true",
            "library": "lib_b",
        },
        {
            "target_id": "T2",
            "ligand_display": "L2",
            "ligand_base": "L2",
            "t_selected": "3.5",
            "rank": "4",
            "pct_rank": "0.4",
            "pose_valid_any": "true",
            "library": "lib_b",
        },
    ]


def test_heatmap_html_clustergrammer_output(tmp_path: Path) -> None:
    pytest.importorskip("clustergrammer")
    repo_root = tmp_path / "repo"
    run_id = "CG_TEST"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")

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
    json.loads(viz_match.group(1))

    meta_match = re.search(
        r'<script type="application/json" id="cg-heatmap-meta-CG_TEST">(.*?)</script>',
        html_text,
        re.S,
    )
    assert meta_match is not None
    cell_meta = json.loads(meta_match.group(1))
    assert cell_meta["L1||T1"]["rank"] == "1"


def test_heatmap_html_clustergrammer_fallback(tmp_path: Path, monkeypatch) -> None:
    repo_root = tmp_path / "repo"
    run_id = "CG_FALLBACK"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_minimal_rows())

    monkeypatch.setattr(heatmap_html, "_load_clustergrammer", lambda: None)
    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=5, include_decoys=False
    )
    assert ".hm-cell" in html_text
