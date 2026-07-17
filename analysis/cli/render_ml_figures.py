from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.figure_diagnostics import render_model_figure_suite


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render Atlas ML figures, SHAP summaries, and failure-analysis tables from model predictions."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--model-dir", type=Path, help="Trained model directory containing model_predictions.csv.")
    source.add_argument("--predictions", type=Path, help="Explicit model_predictions.csv path.")
    parser.add_argument("--dataset", type=Path, help="Optional full training/model-ready dataset for label distribution plots.")
    parser.add_argument("--label", dest="label_col", help="Binary label column. Auto-detected when omitted.")
    parser.add_argument("--score-col", default="ml_prediction_score")
    parser.add_argument("--out-dir", type=Path, help="Output figure directory. Defaults to <model-dir>/figures.")
    parser.add_argument("--group-cols", nargs="*", help="Group columns for label distribution plots.")
    parser.add_argument("--target-cols", nargs="*", help="Columns for target/PDB positive-rate plots.")
    parser.add_argument(
        "--top-ks",
        nargs="*",
        type=int,
        help="K values for Top-K recovery; clipped to row count and deduplicated.",
    )
    parser.add_argument("--max-groups", type=int, default=20)
    parser.add_argument("--min-group-size", type=int, default=5)
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Decision threshold used to filter FP/FN candidate outputs.",
    )
    parser.add_argument("--top-n-errors", type=int, default=50, help="Rows to write in top FP/FN candidate tables.")
    parser.add_argument("--shap-sample-rows", type=int, default=250)
    parser.add_argument("--shap-background-rows", type=int, default=100)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--skip-shap", action="store_true")
    parser.add_argument("--force-generic-shap", action="store_true")
    parser.add_argument("--max-generic-shap-features", type=int, default=200)
    parser.add_argument("--title-prefix")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    manifest = render_model_figure_suite(
        model_dir=args.model_dir,
        predictions=args.predictions,
        dataset=args.dataset,
        label_col=args.label_col,
        score_col=args.score_col,
        out_dir=args.out_dir,
        group_cols=args.group_cols,
        target_cols=args.target_cols,
        top_ks=args.top_ks,
        max_groups=args.max_groups,
        min_group_size=args.min_group_size,
        shap_sample_rows=args.shap_sample_rows,
        shap_background_rows=args.shap_background_rows,
        random_state=args.random_state,
        skip_shap=args.skip_shap,
        force_generic_shap=args.force_generic_shap,
        max_generic_shap_features=args.max_generic_shap_features,
        threshold=args.threshold,
        top_n_errors=args.top_n_errors,
        title_prefix=args.title_prefix,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
