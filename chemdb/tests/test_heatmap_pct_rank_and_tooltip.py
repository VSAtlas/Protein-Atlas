import csv
from pathlib import Path

import pytest

from analysis.reporting import heatmap_html, run_report_core as run_report


def _write_master_rows(path: Path) -> None:
    ligands = [f"Lig{chr(ord('A') + i)}" for i in range(10)]
    rows = []
    for pdb_id in ("PDB1", "PDB2"):
        for ligand in ligands:
            z_selected = ""
            if pdb_id == "PDB1" and ligand == "LigA":
                z_selected = "10"
            elif pdb_id == "PDB1" and ligand == "LigB":
                z_selected = "9"
            elif pdb_id == "PDB2" and ligand == "LigA":
                z_selected = "8"
            rows.append(
                {
                    "run_id": "testrun",
                    "pdb_id": pdb_id,
                    "variant": "V1",
                    "ph_label": "7.4",
                    "ligand_base": ligand,
                    "ligand": f"{ligand}.pdbqt",
                    "ligand_file": f"{ligand}.pdbqt",
                    "library": "FDA",
                    "is_decoy": "false",
                    "is_control": "false",
                    "pose_valid_any": "true",
                    "z_selected": z_selected,
                    "z_selected_source": "consensus",
                }
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _fake_clustergrammer_viz(
    matrix, row_labels, col_labels, _use_dendrogram
):  # type: ignore[no-untyped-def]
    return {
        "row_nodes": [{"name": name} for name in row_labels],
        "col_nodes": [{"name": name} for name in col_labels],
        "mat": [[0.0 if val is None else float(val) for val in row] for row in matrix],
        "links": [],
        "views": [],
        "cat_colors": {},
    }


def _target_entry(report: dict, target_id: str, ligand_base: str) -> dict:
    entries = report["targets"][target_id]["top_ligands_extended"]
    for entry in entries:
        if entry.get("ligand_base") == ligand_base:
            return entry
    raise AssertionError(f"Missing ligand {ligand_base} in target {target_id}")


def test_heatmap_pct_rank_and_tooltip(tmp_path: Path, monkeypatch) -> None:
    repo_root = tmp_path / "repo"
    run_id = "testrun"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True, exist_ok=True)
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

    master_rows = data_dir / "master_rows.csv"
    _write_master_rows(master_rows)

    report = run_report.build_report(
        run_id=run_id,
        repo_root=repo_root,
        top_n=10,
        extended_top_n=20,
        highlight_ligands=None,
        multi_target_min_targets=2,
        multi_target_max_pct=0.2,
        multi_target_max_hits=20,
        filter_invalid=False,
    )

    target1 = "PDB1|V1|7.4"
    target2 = "PDB2|V1|7.4"
    lig_a_t1 = _target_entry(report, target1, "LigA")
    lig_b_t1 = _target_entry(report, target1, "LigB")
    lig_a_t2 = _target_entry(report, target2, "LigA")

    assert lig_a_t1["rank"] == 1
    assert lig_a_t1["pct_rank"] == pytest.approx(0.1)
    assert lig_b_t1["rank"] == 2
    assert lig_b_t1["pct_rank"] == pytest.approx(0.2)
    assert lig_a_t2["rank"] == 1
    assert lig_a_t2["pct_rank"] == pytest.approx(0.1)

    multi_hits = report.get("multi_target_hits") or []
    assert multi_hits, "Expected at least one multi-target hit"
    lig_a_hit = next(
        (hit for hit in multi_hits if hit.get("ligand_base") == "LigA"),
        None,
    )
    assert lig_a_hit is not None, "Expected LigA in multi_target_hits"
    assert lig_a_hit["worst_pct"] == pytest.approx(0.1)
    assert lig_a_hit["mean_pct"] == pytest.approx(0.1)

    heatmap_input = data_dir / "heatmap_input.csv"
    run_report._write_heatmap_input_csv(
        repo_root=repo_root,
        run_id=run_id,
        out_path=heatmap_input,
        top_k=None,
        fda_mapping_csv=None,
        filter_invalid=False,
    )
    with heatmap_input.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    t1_lig_a = next(
        row
        for row in rows
        if row["target_id"] == target1 and row["ligand_base"] == "LigA"
    )
    t1_lig_b = next(
        row
        for row in rows
        if row["target_id"] == target1 and row["ligand_base"] == "LigB"
    )
    assert t1_lig_a["rank"] == "1"
    assert float(t1_lig_a["pct_rank"]) == pytest.approx(0.1)
    assert t1_lig_b["rank"] == "2"
    assert float(t1_lig_b["pct_rank"]) == pytest.approx(0.2)

    monkeypatch.setattr(
        heatmap_html, "_build_clustergrammer_viz_json", _fake_clustergrammer_viz
    )
    html_out = heatmap_html.render_interactive_heatmap_html(
        repo_root=repo_root,
        run_id=run_id,
        input_csv=heatmap_input,
        allowed_pdb_ids=None,
    )
    assert "pose_valid_any:" not in html_out
    assert "library:" not in html_out
    assert "pct_rank:" in html_out
    assert (
        "pct_rank: 0.1" in html_out
        or "pct_rank: 0.100000" in html_out
        or '"pct_rank":"0.1"' in html_out
        or '"pct_rank":"0.100000"' in html_out
        or '"pct_rank": "0.1"' in html_out
        or '"pct_rank": "0.100000"' in html_out
    )
