from __future__ import annotations

import argparse
from pathlib import Path

from analysis.graph.graph_scores import compute_pair_mechanism_scores
from analysis.graph.mechanism_graph import build_mechanism_graph


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Atlas mechanism evidence graph.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--sider", type=Path, default=None)
    parser.add_argument("--ctd", type=Path, default=None)
    parser.add_argument("--opentargets-safety", type=Path, default=None)
    parser.add_argument("--adrecs-target", type=Path, default=None)
    parser.add_argument("--reactome", type=Path, default=None)
    parser.add_argument("--spd", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    nodes_path = args.out_dir / "mechanism_graph_nodes.csv"
    edges_path = args.out_dir / "mechanism_graph_edges.csv"
    scores_path = args.out_dir / "pair_mechanism_scores.csv"
    build_mechanism_graph(
        args.pair_table,
        args.sider,
        args.ctd,
        args.opentargets_safety,
        args.reactome,
        args.spd,
        edges_path,
        nodes_path,
        adrecs_target_path=args.adrecs_target,
    )
    compute_pair_mechanism_scores(args.pair_table, edges_path, scores_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
