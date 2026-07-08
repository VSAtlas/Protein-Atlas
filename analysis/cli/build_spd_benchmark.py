from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.spd import build_spd_benchmark, write_spd_plots
from analysis.io import load_config, output_dirs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build SPD exposure-relevance benchmark.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--spd", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    thresholds = config.get("exposure", {}).get("margin_thresholds", {})
    strong = float(thresholds.get("strong", 10.0))
    weak = float(thresholds.get("weak", 100.0))
    df = build_spd_benchmark(args.pair_table, args.spd, args.mapping, args.out, strong, weak)
    _analysis_dir, figures_dir, _models_dir = output_dirs(config)
    write_spd_plots(df, figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

