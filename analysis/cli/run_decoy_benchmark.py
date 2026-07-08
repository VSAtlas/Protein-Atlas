from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.controls.decoy_benchmark import run_decoy_benchmark


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run a separated decoy/control benchmark panel.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--scores", nargs="*", default=None)
    parser.add_argument("--score-directions-json", default=None)
    parser.add_argument("--positive-col", default="is_control")
    parser.add_argument("--decoy-col", default="is_decoy")
    parser.add_argument("--group-col", default="pdb_id")
    parser.add_argument(
        "--include-fda-background",
        action="store_true",
        help="Sensitivity mode only: include non-control FDA rows as background negatives.",
    )
    args = parser.parse_args(argv)
    score_directions = json.loads(args.score_directions_json) if args.score_directions_json else None
    run_decoy_benchmark(
        args.input,
        args.out_dir,
        score_cols=args.scores,
        score_directions=score_directions,
        positive_col=args.positive_col,
        decoy_col=args.decoy_col,
        group_col=args.group_col,
        include_fda_background=args.include_fda_background,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
