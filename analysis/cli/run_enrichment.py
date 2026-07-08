from __future__ import annotations

import argparse
from pathlib import Path

from analysis.enrichment import DEFAULT_LABEL_COLS, DEFAULT_SCORE_COLS, run_enrichment_suite
from analysis.io import load_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Atlas enrichment suite.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--scores", nargs="*", default=DEFAULT_SCORE_COLS)
    parser.add_argument("--labels", nargs="*", default=DEFAULT_LABEL_COLS)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    top_fractions = config.get("significance", {}).get("top_fractions", [0.01, 0.05, 0.10])
    seed = int(config.get("project", {}).get("random_seed", 42))
    run_enrichment_suite(args.pair_table, args.scores, args.labels, top_fractions, args.permutations, args.out_dir, seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

