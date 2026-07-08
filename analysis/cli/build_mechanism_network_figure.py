from __future__ import annotations

import argparse
from pathlib import Path

from analysis.network.figure3_network import build_mechanism_network_figure


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build a Sawada Figure 3-style Atlas mechanism subgraph.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--mechanism-nodes", required=True, type=Path)
    parser.add_argument("--mechanism-edges", required=True, type=Path)
    parser.add_argument("--pbas-pairs", type=Path, default=None)
    parser.add_argument("--drug", default=None)
    parser.add_argument("--adr", default=None)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    build_mechanism_network_figure(
        args.pair_table,
        args.mechanism_nodes,
        args.mechanism_edges,
        args.out,
        pbas_pairs_path=args.pbas_pairs,
        drug=args.drug,
        adr=args.adr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
