from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.feature_sets import get_feature_set
from analysis.ml.spd_exposure_grouped_oof import run_spd_exposure_grouped_oof


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def _thread_count(value: str) -> int:
    parsed = _positive_int(value)
    if parsed > 4:
        raise argparse.ArgumentTypeError("value must be at most 4")
    return parsed


def _json_object(value: str) -> dict[str, object]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise argparse.ArgumentTypeError(f"invalid JSON: {exc.msg}") from exc
    if not isinstance(parsed, dict):
        raise argparse.ArgumentTypeError("value must be a JSON object")
    return parsed


def _resolve_features(
    feature_set: str | None,
    explicit_features: list[str] | None,
) -> list[str]:
    features = get_feature_set(feature_set) if feature_set else []
    features.extend(explicit_features or [])
    return list(dict.fromkeys(features))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run pooled grouped out-of-fold evaluation for the two-stage SPD "
            "potency/free-Cmax exposure architecture."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label", default="spd_exposure_label")
    parser.add_argument(
        "--potency-feature",
        action="append",
        help=(
            "Additional potency-stage input column; repeat for each variable. "
            "May be combined with --potency-feature-set."
        ),
    )
    parser.add_argument(
        "--potency-feature-set",
        help="Registered analysis.ml.feature_sets name for potency-stage inputs.",
    )
    parser.add_argument(
        "--pk-feature",
        action="append",
        help=(
            "Additional free-Cmax-stage input column; repeat for each variable. "
            "May be combined with --pk-feature-set."
        ),
    )
    parser.add_argument(
        "--pk-feature-set",
        help="Registered analysis.ml.feature_sets name for free-Cmax-stage inputs.",
    )
    parser.add_argument("--group-col", default="drug_id")
    parser.add_argument("--n-splits", type=_positive_int, default=5)
    parser.add_argument(
        "--model",
        choices=["ridge", "random_forest", "lightgbm"],
        default="ridge",
    )
    parser.add_argument("--seed", type=int, default=42, help="Base model seed.")
    parser.add_argument(
        "--model-params-json",
        type=_json_object,
        default=None,
        help=(
            "Estimator overrides as one JSON object. Seed and thread parameters are "
            "controlled separately."
        ),
    )
    parser.add_argument(
        "--max-threads",
        type=_thread_count,
        default=4,
        help="Hard estimator/BLAS thread cap (1-4).",
    )
    parser.add_argument(
        "--fold-assignments",
        type=Path,
        default=None,
        help="Replay this saved fold_assignments.csv instead of generating folds.",
    )
    parser.add_argument(
        "--fold-seed",
        type=int,
        default=None,
        help="Fold-generation seed; defaults to --seed when folds are generated.",
    )
    parser.add_argument(
        "--censored-policy",
        choices=["exclude", "bound"],
        default="exclude",
        help=(
            "Exclude right-censored AC50 rows from potency fitting, or treat "
            "their bounds as exact in a sensitivity analysis."
        ),
    )
    parser.add_argument("--permutation-repeats", type=int, default=10)
    parser.add_argument(
        "--top-k",
        action="append",
        type=_positive_int,
        default=None,
        help="Top-K cutoff; repeat for multiple cutoffs (default: 10, 20, 100).",
    )
    args = parser.parse_args(argv)
    if args.permutation_repeats < 0:
        parser.error("--permutation-repeats must be non-negative")
    try:
        potency_features = _resolve_features(
            args.potency_feature_set,
            args.potency_feature,
        )
        pk_features = _resolve_features(args.pk_feature_set, args.pk_feature)
    except ValueError as exc:
        parser.error(str(exc))
    if not potency_features:
        parser.error("provide --potency-feature-set and/or --potency-feature")
    if not pk_features:
        parser.error("provide --pk-feature-set and/or --pk-feature")

    manifest = run_spd_exposure_grouped_oof(
        args.dataset,
        args.out_dir,
        label_col=args.label,
        potency_features=potency_features,
        pk_features=pk_features,
        group_col=args.group_col,
        n_splits=args.n_splits,
        model_type=args.model,
        model_params=args.model_params_json,
        seed=args.seed,
        max_threads=args.max_threads,
        fold_assignments_path=args.fold_assignments,
        fold_seed=args.fold_seed,
        censored_policy=args.censored_policy,
        permutation_repeats=args.permutation_repeats,
        top_k=args.top_k or (10, 20, 100),
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
