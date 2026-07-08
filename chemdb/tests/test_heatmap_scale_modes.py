import csv
import json
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


def _force_legacy(monkeypatch) -> None:
    def _raise(*_args, **_kwargs) -> None:
        raise ImportError("clustergrammer disabled for scale-mode tests")

    monkeypatch.setattr(heatmap_html, "_build_clustergrammer_viz_json", _raise)


def _extract_scale_meta(html_text: str, run_id: str) -> dict:
    safe_run_id = re.sub(r"[^A-Za-z0-9_-]+", "-", run_id).strip("-") or "run"
    script_id = f"heatmap-scale-meta-{safe_run_id}"
    match = re.search(
        rf'<script type="application/json" id="{re.escape(script_id)}">(.*?)</script>',
        html_text,
        re.S,
    )
    assert match is not None
    return json.loads(match.group(1))


def _build_rows(values: list[float], target_id: str = "T1") -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for idx, value in enumerate(values, start=1):
        ligand = f"L{idx:02d}"
        rows.append(
            {
                "target_id": target_id,
                "ligand_display": ligand,
                "ligand_base": ligand,
                "z_selected": f"{value:.6f}",
            }
        )
    return rows


def test_heatmap_scale_default_fixed_mode(tmp_path: Path, monkeypatch) -> None:
    _force_legacy(monkeypatch)
    repo_root = tmp_path / "repo"
    run_id = "HM_SCALE_FIXED"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_rows([0.2, 0.5, 0.9, 1.2, 1.8, 2.2, 2.8, 3.1, -0.2, -1.0]))

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=50, include_decoys=False
    )
    meta = _extract_scale_meta(html_text, run_id)

    assert "Scale mode: FIXED (absolute)" in html_text
    assert "Scale display:" in html_text
    assert f'id="heatmap-scale-toggle-fixed-{run_id}" value="fixed" checked' in html_text
    assert f'id="heatmap-scale-toggle-quantile-{run_id}" value="quantile"' in html_text
    assert (
        f'id="heatmap-scale-toggle-quantile-{run_id}" value="quantile" disabled'
        not in html_text
    )
    assert meta["scale_mode"] == "fixed"
    assert meta["chosen_breaks"] == meta["fixed_breaks"]
    assert meta["quantile_breaks"] is not None


def test_heatmap_scale_quantile_changes_breaks(tmp_path: Path, monkeypatch) -> None:
    _force_legacy(monkeypatch)
    repo_root = tmp_path / "repo"
    run_id = "HM_SCALE_QUANTILE"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text(
        "\n".join(
            [
                "USE_DENDROGRAM=false",
                "HEATMAP_SCALE_MODE=quantile",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    input_csv = data_dir / "heatmap_input.csv"
    values = [0.82, 0.9, 0.95, 1.01, 1.05, 1.1, 1.2, 1.25, 1.35, 1.5, 1.72, 1.95]
    _write_csv(input_csv, _build_rows(values))

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=50, include_decoys=False
    )
    meta = _extract_scale_meta(html_text, run_id)

    assert meta["scale_mode"] == "quantile"
    assert meta["quantile_breaks"] is not None
    assert meta["chosen_breaks"] != meta["fixed_breaks"]
    assert (
        f'id="heatmap-scale-toggle-quantile-{run_id}" value="quantile" checked'
        in html_text
    )


def test_heatmap_scale_html_labels_and_metadata(tmp_path: Path, monkeypatch) -> None:
    _force_legacy(monkeypatch)
    repo_root = tmp_path / "repo"
    run_id = "HM_SCALE_LABELS"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text(
        "\n".join(
            [
                "USE_DENDROGRAM=false",
                "HEATMAP_SCALE_MODE=quantile",
                "HEATMAP_QUANTILE_LOW=0.05",
                "HEATMAP_QUANTILE_MID2=0.9",
                "HEATMAP_QUANTILE_HIGH=0.97",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    input_csv = data_dir / "heatmap_input.csv"
    values = [0.8, 0.9, 1.0, 1.05, 1.12, 1.17, 1.2, 1.3, 1.4, 1.55, 1.7, 1.88]
    _write_csv(input_csv, _build_rows(values))

    html_text = render_interactive_heatmap_html(
        repo_root, run_id, input_csv, top_k=50, include_decoys=False
    )
    meta = _extract_scale_meta(html_text, run_id)

    assert "Scale mode: QUANTILE (within-run; colors are relative)" in html_text
    assert "Breaks used:" in html_text
    assert "Fixed reference breaks:" in html_text
    assert "chosen_breaks" in html_text
    assert "fixed_breaks" in html_text
    assert meta["scale_mode"] == "quantile"
    assert "chosen_breaks" in meta
    assert "fixed_breaks" in meta
    assert "quantile_breaks" in meta
    assert "n_finite" in meta
    assert "q_low" in meta and "q_mid2" in meta and "q_high" in meta


def test_heatmap_scale_quantile_guardrail_fallback(tmp_path: Path, monkeypatch, caplog) -> None:
    _force_legacy(monkeypatch)
    repo_root = tmp_path / "repo"
    run_id = "HM_SCALE_GUARDRAIL"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    (repo_root / "config.txt").write_text(
        "\n".join(
            [
                "USE_DENDROGRAM=false",
                "HEATMAP_SCALE_MODE=quantile",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    input_csv = data_dir / "heatmap_input.csv"
    _write_csv(input_csv, _build_rows([0.8, 0.9, 1.0, 1.1, 1.15, 1.22, 1.3, 1.4]))

    with caplog.at_level("WARNING", logger="heatmap-html"):
        html_text = render_interactive_heatmap_html(
            repo_root, run_id, input_csv, top_k=50, include_decoys=False
        )
    meta = _extract_scale_meta(html_text, run_id)

    assert "Scale mode: FIXED (absolute)" in html_text
    assert "Quantile fallback:" in html_text
    assert (
        f'id="heatmap-scale-toggle-fixed-{run_id}" value="fixed" checked'
        in html_text
    )
    assert (
        f'id="heatmap-scale-toggle-quantile-{run_id}" value="quantile" disabled'
        in html_text
    )
    assert meta["scale_mode"] == "fixed"
    assert meta["fallback_reason"] == "n_finite_lt_10"
    assert "action=scale_breaks_fallback" in caplog.text


def test_heatmap_r_script_references_scale_mode_keys() -> None:
    r_path = Path("analysis/HeatMap.R")
    text = r_path.read_text(encoding="utf-8")

    assert "HEATMAP_SCALE_MODE" in text
    assert "HEATMAP_QUANTILE_LOW" in text
    assert "HEATMAP_QUANTILE_MID2" in text
    assert "HEATMAP_QUANTILE_HIGH" in text
    assert "scale_sidecar" in text
