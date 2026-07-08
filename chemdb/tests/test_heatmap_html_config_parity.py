import csv
import re
from pathlib import Path

from analysis.reporting import heatmap_html as heatmap_html
from analysis.reporting.heatmap_html import render_interactive_heatmap_html


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _extract_row_labels(html_text: str) -> list[str]:
    tbody_match = re.search(r"<tbody>(.*?)</tbody>", html_text, re.S)
    assert tbody_match is not None
    return re.findall(r'<th class="hm-row-label">(.*?)</th>', tbody_match.group(1))


def _extract_col_labels(html_text: str) -> list[str]:
    thead_match = re.search(r"<thead>(.*?)</thead>", html_text, re.S)
    assert thead_match is not None
    return re.findall(r'<th class="hm-col-label"><div>(.*?)</div></th>', thead_match.group(1))


def _force_legacy(monkeypatch) -> None:
    def _raise(*_args, **_kwargs) -> None:
        raise ImportError("clustergrammer disabled for legacy tests")

    monkeypatch.setattr(heatmap_html, "_build_clustergrammer_viz_json", _raise)


def test_heatmap_html_order_and_na_white(tmp_path: Path, monkeypatch) -> None:
    _force_legacy(monkeypatch)
    repo_root = tmp_path / "repo"
    run_id = "HM_ORDER"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")

    rows = [
        {"target_id": "T1", "ligand_display": "L1", "ligand_base": "L1", "z_selected": "6"},
        {"target_id": "T2", "ligand_display": "L1", "ligand_base": "L1", "z_selected": "1"},
        {"target_id": "T1", "ligand_display": "L2", "ligand_base": "L2", "z_selected": "5"},
        {"target_id": "T1", "ligand_display": "L3", "ligand_base": "L3", "z_selected": "4"},
        {"target_id": "T2", "ligand_display": "L3", "ligand_base": "L3", "z_selected": "2"},
    ]
    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, rows)

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=2, include_decoys=False
    )

    assert _extract_row_labels(html_text) == ["L1", "L2"]
    assert _extract_col_labels(html_text) == ["T1", "T2"]

    match = re.search(
        r'<button[^>]*data-ligand="L2"[^>]*data-target="T2"[^>]*>',
        html_text,
    )
    assert match is not None
    assert "background-color: white" in match.group(0)


def test_heatmap_html_config_palette_bins(tmp_path: Path, monkeypatch) -> None:
    _force_legacy(monkeypatch)
    repo_root = tmp_path / "repo"
    run_id = "HM_COLORS"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)

    config_path = repo_root / "config.txt"
    config_path.write_text(
        "\n".join(
            [
                "USE_DENDROGRAM=false",
                "HEATMAP_SCALE_MIN=-10",
                "HEATMAP_SCALE_MID=-2",
                "HEATMAP_SCALE_MID2=1",
                "HEATMAP_SCALE_MAX=5",
                "HEATMAP_COLOR_MIN=#00ff00",
                "HEATMAP_COLOR_MID=#000000",
                "HEATMAP_COLOR_MID2=#ff8800",
                "HEATMAP_COLOR_MAX=#ff0000",
                "",
            ]
        ),
        encoding="utf-8",
    )

    rows = [
        {"target_id": "T1", "ligand_display": "L1", "ligand_base": "L1", "z_selected": "-10"},
        {"target_id": "T1", "ligand_display": "L2", "ligand_base": "L2", "z_selected": "5"},
    ]
    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, rows)

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=10, include_decoys=False
    )

    min_match = re.search(
        r'<button[^>]*data-ligand="L1"[^>]*data-target="T1"[^>]*>',
        html_text,
    )
    max_match = re.search(
        r'<button[^>]*data-ligand="L2"[^>]*data-target="T1"[^>]*>',
        html_text,
    )
    assert min_match is not None
    assert max_match is not None
    assert "background-color: #00ff00" in min_match.group(0)
    assert "background-color: #ff0000" in max_match.group(0)

    config_path.write_text(
        "\n".join(
            [
                "USE_DENDROGRAM=false",
                "HEATMAP_SCALE_MIN=-10",
                "HEATMAP_SCALE_MID=-2",
                "HEATMAP_SCALE_MID2=1",
                "HEATMAP_SCALE_MAX=5",
                "HEATMAP_COLOR_MIN=#112233",
                "HEATMAP_COLOR_MID=#000000",
                "HEATMAP_COLOR_MID2=#ff8800",
                "HEATMAP_COLOR_MAX=#ff0000",
                "",
            ]
        ),
        encoding="utf-8",
    )
    html_text_alt = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=10, include_decoys=False
    )
    assert html_text_alt != html_text
    alt_match = re.search(
        r'<button[^>]*data-ligand="L1"[^>]*data-target="T1"[^>]*>',
        html_text_alt,
    )
    assert alt_match is not None
    assert "background-color: #112233" in alt_match.group(0)
