from __future__ import annotations

import argparse
from pathlib import Path

from analysis.build_pair_table import DERIVED_COLUMNS, REQUIRED_COLUMNS, build_pair_feature_table
from analysis.enrichment import DEFAULT_LABEL_COLS, DEFAULT_SCORE_COLS, run_ablation_analysis, run_enrichment_suite
from analysis.external.labels import build_external_label_tables
from analysis.external.spd import build_spd_benchmark, write_spd_plots
from analysis.external.toxcast import build_toxcast_benchmark
from analysis.generate_heatmap import build_heatmap_matrix, write_matrix_csv, write_png, write_svg
from analysis.graph.graph_scores import compute_pair_mechanism_scores
from analysis.graph.mechanism_graph import build_mechanism_graph
from analysis.io import load_config, output_dirs
from analysis.run_enrichment import write_enrichment_curve
from analysis.run_pair_significance import OUTPUT_COLUMNS, compute_pair_significance
from analysis.schemas import write_pair_table_sidecars
from analysis._common import read_csv_rows, write_csv_rows


def _path(config: dict, *keys: str, default: str = "") -> Path:
    value = config
    for key in keys:
        value = value.get(key, {}) if isinstance(value, dict) else {}
    return Path(str(value or default))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run expanded Atlas analysis from config.")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    analysis_dir, figures_dir, _models_dir = output_dirs(config)
    analysis_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    inputs = config.get("inputs", {})
    seed = int(config.get("project", {}).get("random_seed", 42))
    pair_path = Path(inputs.get("pair_table", analysis_dir / "pair_feature_table.csv"))
    screening = Path(inputs.get("docking_results", "outputs/docking/ranked_pairs.csv"))
    annotations = Path(inputs.get("annotations_dir", "data/annotations"))
    if not pair_path.exists() and screening.exists():
        pair_rows = build_pair_feature_table(screening, annotations, config=config)
        write_csv_rows(pair_path, pair_rows, [*REQUIRED_COLUMNS, *DERIVED_COLUMNS])
        write_pair_table_sidecars(pair_path, pair_rows)
    pair_rows = read_csv_rows(pair_path)
    sig = compute_pair_significance(pair_rows, score_field="atlas_score", direction="higher")
    sig_path = analysis_dir / "pair_significance.csv"
    write_csv_rows(sig_path, sig, OUTPUT_COLUMNS)
    rows, cols, matrix, _row_meta, _col_meta = build_heatmap_matrix(pair_rows, value="atlas_score", row_order="protein_class", col_order="ligand_chemotype")
    write_matrix_csv(analysis_dir / "heatmap_matrix.csv", rows, cols, matrix)
    write_png(figures_dir / "atlas_heatmap.png", rows, cols, matrix, value="atlas_score", markers={})
    write_svg(figures_dir / "atlas_heatmap.svg", rows, cols, matrix, value="atlas_score", markers={})
    external_label_path = analysis_dir / "external_labels" / "atlas_external_labels.csv"
    build_external_label_tables(pair_path, config, analysis_dir / "external_labels")
    enrichment_pair_path = external_label_path if external_label_path.exists() else pair_path
    write_enrichment_curve(figures_dir / "enrichment_curve.png", pair_rows, score="atlas_score", label="literature_supported_label")
    run_enrichment_suite(enrichment_pair_path, DEFAULT_SCORE_COLS, DEFAULT_LABEL_COLS, config.get("significance", {}).get("top_fractions", [0.01, 0.05, 0.10]), int(config.get("permutation", {}).get("n_permutations", 1000)), analysis_dir / "enrichment", seed, config.get("score_directions", {}))
    run_ablation_analysis(pair_path, analysis_dir / "ablation_summary.csv")
    spd_file = inputs.get("spd_file")
    if spd_file and Path(spd_file).exists():
        spd = build_spd_benchmark(pair_path, spd_file, inputs.get("mapping_file"), analysis_dir / "spd_atlas_benchmark.csv")
        write_spd_plots(spd, figures_dir)
    toxcast_file = inputs.get("toxcast_file")
    if toxcast_file and Path(toxcast_file).exists():
        build_toxcast_benchmark(pair_path, toxcast_file, inputs.get("mapping_file"), analysis_dir / "toxcast_atlas_benchmark.csv")
    graph_dir = analysis_dir / "graph"
    graph_dir.mkdir(parents=True, exist_ok=True)
    nodes = graph_dir / "mechanism_graph_nodes.csv"
    edges = graph_dir / "mechanism_graph_edges.csv"
    build_mechanism_graph(pair_path, inputs.get("sider_file"), inputs.get("ctd_file"), inputs.get("opentargets_safety_file"), inputs.get("reactome_file"), analysis_dir / "spd_atlas_benchmark.csv", edges, nodes)
    compute_pair_mechanism_scores(pair_path, edges, graph_dir / "pair_mechanism_scores.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
