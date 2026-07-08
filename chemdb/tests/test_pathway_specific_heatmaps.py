import csv
from pathlib import Path

from analysis import pathway_resolver
from analysis.reporting import heatmap_html, run_report_core as run_report


def _write_heatmap_input(path: Path) -> None:
    rows = [
        {
            "pdb_id": "PDB1",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "target_id": "PDB1|HOLO|pH7_0",
            "ligand_display": "LIG_A",
            "ligand_base": "LIG_A",
            "z_selected": "1.5",
            "rank": "1",
            "pct_rank": "0.01",
        },
        {
            "pdb_id": "PDB1",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "target_id": "PDB1|HOLO|pH7_0",
            "ligand_display": "LIG_B",
            "ligand_base": "LIG_B",
            "z_selected": "2.0",
            "rank": "2",
            "pct_rank": "0.02",
        },
        {
            "pdb_id": "PDB2",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "target_id": "PDB2|HOLO|pH7_0",
            "ligand_display": "LIG_A",
            "ligand_base": "LIG_A",
            "z_selected": "1.0",
            "rank": "1",
            "pct_rank": "0.01",
        },
        {
            "pdb_id": "PDB2",
            "variant": "HOLO",
            "ph_label": "pH7_0",
            "target_id": "PDB2|HOLO|pH7_0",
            "ligand_display": "LIG_C",
            "ligand_base": "LIG_C",
            "z_selected": "2.3",
            "rank": "2",
            "pct_rank": "0.02",
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


def test_pathway_specific_heatmaps_for_dualinhibitortest8(
    tmp_path: Path, monkeypatch
) -> None:
    repo_root = tmp_path / "repo"
    run_id = "dualinhibitortest8"
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
    heatmap_input = data_dir / "heatmap_input.csv"
    _write_heatmap_input(heatmap_input)

    monkeypatch.setattr(
        heatmap_html, "_build_clustergrammer_viz_json", _fake_clustergrammer_viz
    )

    def _fake_map_pdb_to_uniprots(
        pdb_id: str, _cache: pathway_resolver.Cache, _http: pathway_resolver.HttpClient
    ) -> list[str]:
        return {
            "PDB1": ["UNIPROTA"],
            "PDB2": ["UNIPROTB"],
        }.get(pdb_id.upper(), [])

    def _fake_map_uniprot_to_reactome_pathways(
        uniprot: str,
        _cache: pathway_resolver.Cache,
        _http: pathway_resolver.HttpClient,
        organism: str = pathway_resolver.DEFAULT_ORGANISM,
    ) -> list[dict[str, str]]:
        _ = organism
        if uniprot.upper() == "UNIPROTA":
            return [{"name": "Glycolysis", "stId": "R-HSA-70171"}]
        if uniprot.upper() == "UNIPROTB":
            return [{"name": "Pentose phosphate pathway", "stId": "R-HSA-71336"}]
        return []

    monkeypatch.setattr(pathway_resolver, "map_pdb_to_uniprots", _fake_map_pdb_to_uniprots)
    monkeypatch.setattr(
        pathway_resolver,
        "map_uniprot_to_reactome_pathways",
        _fake_map_uniprot_to_reactome_pathways,
    )

    report = {
        "run_id": run_id,
        "generated_at": "2026-01-01T00:00:00Z",
        "targets": {},
        "multi_target_hits": [],
    }
    run_report._write_html_reports(report, data_dir / "report.html", [], repo_root)

    overall_report = data_dir / "report.html"
    assert overall_report.exists()

    pathway_reports = sorted(
        p.name for p in data_dir.glob("report_*.html") if p.is_file()
    )
    assert len(pathway_reports) == 2
    assert pathway_reports == [
        "report_r_hsa_70171_human_glycolysis.html",
        "report_r_hsa_71336_human_pentose_phosphate_pathway.html",
    ]

    glycolysis_text = (
        data_dir / "report_r_hsa_70171_human_glycolysis.html"
    ).read_text(encoding="utf-8")
    pentose_text = (
        data_dir / "report_r_hsa_71336_human_pentose_phosphate_pathway.html"
    ).read_text(encoding="utf-8")
    assert "PDB1|HOLO|pH7_0" in glycolysis_text
    assert "PDB2|HOLO|pH7_0" not in glycolysis_text
    assert "PDB2|HOLO|pH7_0" in pentose_text
    assert "PDB1|HOLO|pH7_0" not in pentose_text

    csv_or_parquet = {
        path.relative_to(data_dir).as_posix()
        for path in data_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in {".csv", ".parquet"}
    }
    assert csv_or_parquet.issubset(
        {
            "heatmap_input.csv",
            "artifacts/heatmap_input.csv",
        }
    )
