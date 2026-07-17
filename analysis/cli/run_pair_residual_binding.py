from __future__ import annotations

import argparse
import json

from analysis.ml.feature_sets import get_feature_set
from analysis.ml.pair_residual_binding import (
    SUPPORTED_VARIANTS,
    run_pair_residual_binding,
)


def _resolve_features(
    feature_set: str | None,
    explicit_features: list[str] | None,
) -> list[str]:
    features = get_feature_set(feature_set) if feature_set else []
    features.extend(explicit_features or [])
    return list(dict.fromkeys(features))


def _parse_permutation_block(value: str) -> tuple[str, tuple[str, ...]]:
    name, separator, raw_features = value.partition("=")
    name = name.strip()
    features = tuple(
        dict.fromkeys(
            feature.strip() for feature in raw_features.split(",") if feature.strip()
        )
    )
    if not separator or not name or not features:
        raise argparse.ArgumentTypeError(
            "permutation blocks must use NAME=feature1,..."
        )
    return name, features


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
        help="Additional drug-level descriptor; repeat for each variable.",
    )
    parser.add_argument(
        "--ligand-feature-set",
        help="Registered analysis.ml.feature_sets name for ligand inputs.",
    )
    parser.add_argument(
        "--pair-feature",
        action="append",
        help="Additional pair-level residual feature; repeat for each variable.",
    )
    parser.add_argument(
        "--pair-feature-set",
        help="Registered analysis.ml.feature_sets name for pair inputs.",
    )
    parser.add_argument(
        "--protein-feature",
        action="append",
        default=[],
        help="Target/pocket feature column; repeat for each feature.",
    )
    parser.add_argument(
        "--permutation-repeats",
        type=int,
        default=0,
        help="Grouped pocket permutation repeats. Default: 0 (disabled).",
    )
    parser.add_argument(
        "--permutation-group",
        help="Whole-profile permutation group. Default: the selected target column.",
    )
    parser.add_argument(
        "--permutation-block",
        action="append",
        type=_parse_permutation_block,
        default=[],
        metavar="NAME=FEATURE1,...",
        help=(
            "Named pocket feature block; repeat for additional blocks. Features must "
            "also be supplied with --protein-feature."
        ),
    )
    parser.add_argument(
        "--permutation-all-pairs",
        action="store_true",
        help="Evaluate every pair of supplied protein features in addition to singles.",
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
    parser = build_parser()
    args = parser.parse_args(argv)
    permutation_blocks: dict[str, tuple[str, ...]] = {}
    for name, features in args.permutation_block:
        if name in permutation_blocks:
            parser.error(f"duplicate --permutation-block name: {name}")
        permutation_blocks[name] = features
    try:
        ligand_features = _resolve_features(
            args.ligand_feature_set,
            args.ligand_feature,
        )
        pair_features = _resolve_features(args.pair_feature_set, args.pair_feature)
    except ValueError as exc:
        parser.error(str(exc))
    if not ligand_features:
        parser.error("provide --ligand-feature-set and/or --ligand-feature")
    if not pair_features:
        parser.error("provide --pair-feature-set and/or --pair-feature")
    manifest = run_pair_residual_binding(
        args.dataset,
        args.out_dir,
        ligand_features=ligand_features,
        pair_features=pair_features,
        protein_features=args.protein_feature,
        permutation_repeats=args.permutation_repeats,
        permutation_group=args.permutation_group,
        permutation_blocks=permutation_blocks,
        permutation_all_pairs=args.permutation_all_pairs,
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
