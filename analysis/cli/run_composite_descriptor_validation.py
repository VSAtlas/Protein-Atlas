from __future__ import annotations

import argparse
import json

from analysis.ml.composite_group_validation import run_composite_descriptor_validation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run exploratory RDKit-only cold-drug and composite double-cold validation "
            "with drug-block descriptor permutation diagnostics."
        )
    )
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--label", default="spd_binding_label")
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--min-test-positives", type=int, default=10)
    parser.add_argument("--min-test-negatives", type=int, default=10)
    parser.add_argument("--permutation-repeats", type=int, default=5)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--run-id")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = run_composite_descriptor_validation(
        args.dataset,
        args.out_dir,
        label_col=args.label,
        n_splits=args.n_splits,
        min_test_positives=args.min_test_positives,
        min_test_negatives=args.min_test_negatives,
        permutation_repeats=args.permutation_repeats,
        top_k=args.top_k,
        seed=args.seed,
        run_id=args.run_id,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
