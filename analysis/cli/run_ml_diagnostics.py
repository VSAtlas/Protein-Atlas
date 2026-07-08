from __future__ import annotations

import argparse
from pathlib import Path

from analysis.io import load_config
from analysis.ml.model_diagnostics import run_ml_diagnostics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Atlas ML diagnostics: learning curves, ablations, baselines.")
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--feature-set", default="pilot_binding_only")
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--split", default="drug_holdout")
    parser.add_argument("--exclude-features", nargs="*", default=None)
    parser.add_argument("--train-fractions", nargs="*", type=float, default=None)
    parser.add_argument("--class-weight", choices=["balanced", "none"], default="balanced")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    seed = int(config.get("project", {}).get("random_seed", 42))
    run_ml_diagnostics(
        args.dataset,
        args.label,
        args.feature_set,
        args.model,
        args.split,
        args.out_dir,
        seed=seed,
        exclude_features=args.exclude_features,
        train_fractions=args.train_fractions,
        class_weight=None if args.class_weight == "none" else args.class_weight,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
