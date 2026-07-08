import csv
from pathlib import Path

from analysis.reporting import heatmap_html, run_report_core as run_report


def _write_input_pdb(path: Path) -> None:
    lines = [
        "HEADER    TEST TARGET\n",
        "COMPND    MOL_ID: 1;\n",
        "COMPND   2 MOLECULE: EGFR KINASE;\n",
        "DBREF  TEST A    1   10  UNP    P00533   EGFR_HUMAN      1    10\n",
        "ATOM      1  N   GLY A   1      11.104  13.207   8.668  1.00 20.00           N\n",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="utf-8")


def _write_master_rows(path: Path) -> None:
    rows = [
        {
            "run_id": "testrun",
            "pdb_id": "TEST",
            "variant": "V1",
            "ph_label": "7.4",
            "ligand_base": "LigA",
            "ligand": "LigA.pdbqt",
            "ligand_file": "LigA.pdbqt",
            "library": "FDA",
            "is_decoy": "false",
            "is_control": "false",
            "pose_valid_any": "true",
            "z_selected": "10",
            "z_selected_source": "consensus",
        },
        {
            "run_id": "testrun",
            "pdb_id": "TEST",
            "variant": "V1",
            "ph_label": "7.4",
            "ligand_base": "LigB",
            "ligand": "LigB.pdbqt",
            "ligand_file": "LigB.pdbqt",
            "library": "FDA",
            "is_decoy": "false",
            "is_control": "false",
            "pose_valid_any": "true",
            "z_selected": "9",
            "z_selected_source": "consensus",
        },
    ]
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


def test_heatmap_target_name_labels(tmp_path: Path, monkeypatch) -> None:
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
                "TARGET_NAME_PREFER=auto",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    _write_input_pdb(repo_root / "input_pdbs" / "TEST.pdb")
    _write_master_rows(data_dir / "master_rows.csv")

    report = run_report.build_report(
        run_id=run_id,
        repo_root=repo_root,
        top_n=5,
        extended_top_n=5,
        highlight_ligands="LigA",
        highlight_max_per_target=1,
        filter_invalid=False,
    )
    target_id = "TEST|V1|7.4"
    target_block = report["targets"][target_id]
    assert "target_name" in target_block
    assert "EGFR KINASE" in str(target_block["target_name"])

    heatmap_csv = data_dir / "heatmap_input.csv"
    run_report._write_heatmap_input_csv(
        repo_root=repo_root,
        run_id=run_id,
        out_path=heatmap_csv,
        top_k=None,
        fda_mapping_csv=None,
        filter_invalid=False,
    )

    with heatmap_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        assert "target_name" in (reader.fieldnames or [])
    assert any("EGFR KINASE" in (row.get("target_name") or "") for row in rows)

    monkeypatch.setattr(
        heatmap_html, "_build_clustergrammer_viz_json", _fake_clustergrammer_viz
    )
    html_with_name = heatmap_html.render_interactive_heatmap_html(
        repo_root=repo_root,
        run_id=run_id,
        input_csv=heatmap_csv,
    )
    assert "EGFR KINASE" in html_with_name
    assert "TEST" in html_with_name
    assert "EGFR KINASE (TEST)" not in html_with_name

    no_name_csv = data_dir / "heatmap_input_no_target_name.csv"
    no_name_fields = [f for f in (reader.fieldnames or []) if f != "target_name"]
    with no_name_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=no_name_fields)
        writer.writeheader()
        for row in rows:
            filtered = {k: row.get(k, "") for k in no_name_fields}
            writer.writerow(filtered)

    html_without_name = heatmap_html.render_interactive_heatmap_html(
        repo_root=repo_root,
        run_id=run_id,
        input_csv=no_name_csv,
    )
    assert "TEST" in html_without_name
    assert "EGFR KINASE" not in html_without_name
