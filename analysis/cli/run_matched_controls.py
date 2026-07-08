from __future__ import annotations

import argparse
from pathlib import Path

from analysis.controls.matched_controls import run_matched_controls, write_matched_control_plot
from analysis.io import load_config, output_dirs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run matched-control empirical significance.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--score", default="atlas_score")
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--n-controls", type=int, default=100)
    parser.add_argument("--match-on", nargs="*", default=["protein_class", "ligand_chemotype"])
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    seed = int(config.get("project", {}).get("random_seed", 42))
    run_matched_controls(args.pair_table, args.score, args.top_n, args.n_controls, args.match_on, args.out, seed)
    _analysis_dir, figures_dir, _models_dir = output_dirs(config)
    write_matched_control_plot(args.out, args.out.with_name("matched_controls.csv"), figures_dir / "matched_control_score_distributions.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

