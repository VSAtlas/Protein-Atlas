from __future__ import annotations

import argparse
import json

from analysis.ml.pair_residual_binding import (
    SUPPORTED_VARIANTS,
    run_pair_residual_binding,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run an exploratory grouped OOF binding architecture evaluation with "
            "cross-fitted drug propensity offsets and pair-level residual logits."
        )
    )
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--ligand-feature",
        action="append",
        required=True,
        help="Drug-level ligand descriptor column; repeat for each descriptor.",
    )
    parser.add_argument(
        "--pair-feature",
        action="append",
        required=True,
        help="Pair-level residual feature column; repeat for each feature.",
    )
    parser.add_argument(
        "--outer-group",
        action="append",
        default=None,
        help=(
            "One-axis cold group, or a two-axis true double-cold specification such as "
            "drug_id+target_id; repeat to evaluate multiple grouped OOF schemes. "
            "Default: drug_id."
        ),
    )
    parser.add_argument(
        "--variant",
        action="append",
        choices=[
            *SUPPORTED_VARIANTS,
            *(value.replace("_", "-") for value in SUPPORTED_VARIANTS),
        ],
        default=None,
        help=(
            "Architecture variant; repeat as needed. Default: pair-only, ligand-only, "
            "and combined. Shortcut-reduced requires --shortcut-feature."
        ),
    )
    parser.add_argument(
        "--shortcut-feature",
        action="append",
        default=[],
        help=(
            "Selected ligand feature removed from the shortcut-reduced propensity model; "
            "repeat as needed."
        ),
    )
    parser.add_argument(
        "--residual-model",
        action="append",
        choices=["logistic", "lightgbm"],
        default=None,
        help="Stage-two residual model; repeat to run both. Default: logistic.",
    )
    parser.add_argument("--label", default="spd_binding_label")
    parser.add_argument("--drug-col", default="drug_id")
    parser.add_argument("--target-col", default="target_id")
    parser.add_argument("--family-col", default="target_family")
    parser.add_argument("--source-col", default="label_source")
    parser.add_argument("--outer-splits", type=int, default=5)
    parser.add_argument("--inner-splits", type=int, default=5)
    parser.add_argument("--top-k", action="append", type=int, default=None)
    parser.add_argument("--propensity-l2", type=float, default=1.0)
    parser.add_argument("--residual-l2", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--model-n-jobs", type=int, default=4)
    parser.add_argument("--run-id")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = run_pair_residual_binding(
        args.dataset,
        args.out_dir,
        ligand_features=args.ligand_feature,
        pair_features=args.pair_feature,
        outer_group_specs=args.outer_group or ["drug_id"],
        variants=args.variant or ["pair_only", "ligand_only", "combined"],
        shortcut_features=args.shortcut_feature,
        residual_models=args.residual_model or ["logistic"],
        label_col=args.label,
        drug_col=args.drug_col,
        target_col=args.target_col,
        family_col=args.family_col,
        source_col=args.source_col,
        outer_splits=args.outer_splits,
        inner_splits=args.inner_splits,
        top_ks=args.top_k or [5, 10, 20],
        propensity_l2=args.propensity_l2,
        residual_l2=args.residual_l2,
        seed=args.seed,
        model_n_jobs=args.model_n_jobs,
        run_id=args.run_id,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
