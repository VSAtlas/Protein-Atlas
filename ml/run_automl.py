#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from ml.automl import run_automl


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run AutoML search for BigBind ML workflow.")
    parser.add_argument(
        "--config",
        default="ml/config.json",
        help="Path to ml config file (JSON or YAML).",
    )
    parser.add_argument(
        "--mode",
        choices=("grid", "optuna"),
        default="grid",
        help="Search mode.",
    )
    parser.add_argument(
        "--max_trials",
        type=int,
        default=50,
        help="Max Optuna trials (used only in --mode optuna).",
    )
    parser.add_argument(
        "--top_k_report",
        type=int,
        default=5,
        help="Top-K inner-ranked trials to evaluate on holdout for final report (K<=5).",
    )
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    run_dir = run_automl(
        config_path=args.config,
        mode=args.mode,
        max_trials=args.max_trials,
        top_k_report=args.top_k_report,
    )
    best_trial_path = run_dir / "best_trial.json"
    best_trial = {}
    if best_trial_path.exists():
        best_trial = json.loads(best_trial_path.read_text(encoding="utf-8"))
    print(
        "[ml.automl] "
        f"saved={run_dir} "
        f"best_trial={best_trial.get('trial_id', 'n/a')} "
        f"variant={best_trial.get('feature_variant', 'n/a')} "
        f"model={best_trial.get('model_family', 'n/a')} "
        f"inner_EF@1%={best_trial.get('inner_metrics', {}).get('inner_mean_EF@1%', 'n/a')} "
        f"inner_PR_AUC={best_trial.get('inner_metrics', {}).get('inner_mean_PR_AUC', 'n/a')}"
    )


if __name__ == "__main__":
    main()
