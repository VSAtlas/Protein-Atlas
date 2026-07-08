from __future__ import annotations

import argparse
from pathlib import Path

from analysis.calibration.calibrate_atlas_score import run_atlas_calibration
from analysis.io import load_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Calibrate Atlas scores against external labels.")
    parser.add_argument("--benchmark-table", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--score", default="atlas_score")
    parser.add_argument("--methods", nargs="*", default=["logistic", "isotonic"])
    parser.add_argument("--split", default="random")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    seed = int(config.get("project", {}).get("random_seed", 42))
    run_atlas_calibration(args.benchmark_table, args.label, args.score, args.out_dir, args.methods, args.split, seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

