from __future__ import annotations

from pathlib import Path

from analysis.generate_heatmap import build_heatmap_matrix, write_png, write_svg


def test_heatmap_outputs_preserve_z_score_hover_text(tmp_path: Path) -> None:
    rows = [
        {
            "drug_id": "drug_a",
            "target_id": "target_a",
            "atlas_score": "2.5",
            "protein_class": "Kinase",
            "ligand_chemotype": "Aromatic",
        }
    ]

    row_labels, col_labels, matrix, _row_meta, _col_meta = build_heatmap_matrix(
        rows,
        value="atlas_score",
        row_order="protein_class",
        col_order="ligand_chemotype",
    )
    svg_path = tmp_path / "atlas_heatmap.svg"
    png_path = tmp_path / "atlas_heatmap.png"
    write_svg(svg_path, row_labels, col_labels, matrix, value="atlas_score", markers={("target_a", "drug_a"): "*"})
    write_png(png_path, row_labels, col_labels, matrix, value="atlas_score", markers={})

    svg_text = svg_path.read_text(encoding="utf-8")
    assert "atlas_score z score=2.5" in svg_text
    assert png_path.exists()
    assert png_path.stat().st_size > 0

