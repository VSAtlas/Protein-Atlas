from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis._common import format_float, parse_float, write_csv_rows
from analysis.build_pair_table import (
    DERIVED_COLUMNS,
    REQUIRED_COLUMNS,
    build_pair_feature_table,
)
from analysis.generate_heatmap import build_heatmap_matrix, write_matrix_csv, write_png, write_svg
from analysis.run_enrichment import (
    PERMUTATION_COLUMNS,
    SUMMARY_COLUMNS,
    run_enrichment,
    write_enrichment_curve,
)
from analysis.run_pair_significance import OUTPUT_COLUMNS, compute_pair_significance
from analysis.schemas import write_pair_table_sidecars


def _write_significance_heatmap(
    pair_rows: Sequence[Mapping[str, Any]],
    significance_rows: Sequence[Mapping[str, Any]],
    out_path: Path,
) -> None:
    q_by_pair = {
        (str(row.get("target_id") or ""), str(row.get("drug_id") or "")): parse_float(row.get("fdr_q_value"))
        for row in significance_rows
    }
    sig_rows: list[dict[str, object]] = []
    for row in pair_rows:
        key = (str(row.get("target_id") or ""), str(row.get("drug_id") or ""))
        q_value = q_by_pair.get(key)
        score = -math.log10(max(q_value, 1e-300)) if q_value is not None and q_value > 0 else None
        sig_rows.append(
            {
                "target_id": key[0],
                "drug_id": key[1],
                "significance_score": format_float(score),
                "protein_class": row.get("protein_class", ""),
                "ligand_chemotype": row.get("ligand_chemotype", ""),
            }
        )
    rows, cols, matrix, _row_meta, _col_meta = build_heatmap_matrix(
        sig_rows,
        value="significance_score",
        row_order="protein_class",
        col_order="ligand_chemotype",
    )
    write_png(out_path, rows, cols, matrix, value="-log10(q_value)", markers={})
    write_svg(out_path.with_suffix(".svg"), rows, cols, matrix, value="-log10(q_value)", markers={})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the Atlas pair-table, heatmap, significance, and enrichment pipeline.")
    parser.add_argument("--screening-results", required=True, type=Path)
    parser.add_argument("--annotations", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--heatmap-value", default="atlas_score")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1337)
    args = parser.parse_args(argv)

    analysis_dir = args.out_dir / "analysis"
    figures_dir = args.out_dir / "figures"
    pair_table_path = analysis_dir / "pair_feature_table.csv"
    matrix_path = analysis_dir / "heatmap_matrix.csv"
    significance_path = analysis_dir / "pair_significance.csv"
    enrichment_path = analysis_dir / "enrichment_summary.csv"
    permutation_path = analysis_dir / "permutation_results.csv"

    pair_rows = build_pair_feature_table(args.screening_results, args.annotations)
    write_csv_rows(pair_table_path, pair_rows, [*REQUIRED_COLUMNS, *DERIVED_COLUMNS])
    write_pair_table_sidecars(pair_table_path, pair_rows)
    pair_table_path.with_suffix(".json").write_text(json.dumps(pair_rows, indent=2), encoding="utf-8")

    significance_rows = compute_pair_significance(pair_rows, score_field="atlas_score", direction="higher")
    write_csv_rows(significance_path, significance_rows, OUTPUT_COLUMNS)

    rows, cols, matrix, _row_meta, _col_meta = build_heatmap_matrix(
        pair_rows,
        value=args.heatmap_value,
        row_order="protein_class",
        col_order="ligand_chemotype",
    )
    write_matrix_csv(matrix_path, rows, cols, matrix)
    marker_rows = {(str(row["target_id"]), str(row["drug_id"])): str(row["fdr_q_value"]) for row in significance_rows}
    markers = {}
    for key, q_value in marker_rows.items():
        parsed_q = parse_float(q_value)
        if parsed_q is not None and parsed_q < 0.05:
            markers[key] = "***"
        elif parsed_q is not None and parsed_q < 0.10:
            markers[key] = "**"
    write_png(figures_dir / "atlas_heatmap.png", rows, cols, matrix, value=args.heatmap_value, markers=markers)
    write_svg(figures_dir / "atlas_heatmap.svg", rows, cols, matrix, value=args.heatmap_value, markers=markers)
    _write_significance_heatmap(pair_rows, significance_rows, figures_dir / "significance_heatmap.png")

    enrichment_rows, permutation_rows = run_enrichment(
        pair_rows,
        score="atlas_score",
        label="literature_supported_label",
        n_permutations=args.permutations,
        n_bootstraps=args.bootstraps,
        top_ks=[10, 25, 50],
        seed=args.seed,
    )
    write_csv_rows(enrichment_path, enrichment_rows, SUMMARY_COLUMNS)
    write_csv_rows(permutation_path, permutation_rows, PERMUTATION_COLUMNS)
    write_enrichment_curve(
        figures_dir / "enrichment_curve.png",
        pair_rows,
        score="atlas_score",
        label="literature_supported_label",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
