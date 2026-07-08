from __future__ import annotations

import argparse
from pathlib import Path

from analysis.baselines.sawada_mode import run_sawada_side_effect_baseline


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Sawada/PBAS-style drug-level side-effect baseline.")
    parser.add_argument("--profile", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--model", default="logistic_regression", choices=("logistic_regression", "random_forest"))
    parser.add_argument("--split", default="drug_holdout")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)
    run_sawada_side_effect_baseline(args.profile, args.labels, args.out_dir, model_type=args.model, split_mode=args.split, seed=args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
