import csv
import shutil
from pathlib import Path

from analysis.reporting import run_report_core as run_report


FIXTURE = Path(__file__).parent / "fixtures" / "report_demo" / "heatmap_input.csv"


def test_report_demo_fixture_generates_html_without_external_tools(
    tmp_path, monkeypatch
) -> None:
    repo_root = tmp_path / "repo"
    run_id = "DEMO"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    shutil.copyfile(FIXTURE, data_dir / "heatmap_input.csv")
    (data_dir / "report.yaml").write_text("run_id: DEMO\n", encoding="utf-8")
    (repo_root / "config.txt").write_text("USE_DENDROGRAM=false\n", encoding="utf-8")

    rows = list(csv.DictReader(FIXTURE.open(newline="", encoding="utf-8")))
    assert len({row["ligand_base"] for row in rows}) >= 3
    assert len({row["target_id"] for row in rows}) >= 3
    assert any(float(row["pct_rank"]) <= 0.01 for row in rows)
    assert any(float(row["pct_rank"]) >= 0.90 for row in rows)
    assert any(str(row["pose_valid_any"]).lower() == "false" for row in rows)

    monkeypatch.setattr(
        run_report,
        "render_interactive_heatmap_html",
        lambda *_args, **_kwargs: "<div>Ligand Strong | Hydrolase Gamma</div>",
    )

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        {
            "run_id": run_id,
            "generated_at": "2026-05-02T00:00:00Z",
            "targets": {},
            "summary": {"n_unique_ligands": 3, "n_target_combos": 3},
        },
        out_path,
        [],
        repo_root,
        heatmap_source_override=data_dir / "heatmap_input.csv",
        report_asset_mode="cdn",
    )

    html_text = out_path.read_text(encoding="utf-8")
    assert "Atlas2 Report: DEMO" in html_text
    assert "Ligand Strong" in html_text
    assert "Interactive Docking Heatmap" in html_text
