from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.leakage_overlap_audit import DEFAULT_SPLITS
from analysis.ml.audit_suite import run_ml_audit_suite


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the reusable Atlas ML data-card, leakage, source-transfer, and benchmark-independence audit suite."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--feature-set")
    parser.add_argument("--exclude-features", nargs="*", default=[])
    parser.add_argument("--splits", nargs="*", default=DEFAULT_SPLITS)
    parser.add_argument("--source-col", default="label_source")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--allow-partial-rescoring-features", action="store_true")
    parser.add_argument(
        "--allow-label-definition-features",
        action="store_true",
        help="Allow label-defining evidence fields for explicitly named sensitivity audits only.",
    )
    parser.add_argument(
        "--candidate",
        nargs="*",
        default=[],
        help="Optional benchmark/calibration tables as NAME=PATH or PATH for train-vs-benchmark independence checks.",
    )
    parser.add_argument("--strict-source-overlap", action="store_true")
    parser.add_argument("--source-transfer-train", nargs="*", default=None)
    parser.add_argument("--source-transfer-test", nargs="*", default=None)
    parser.add_argument("--predictions", type=Path, default=None)
    parser.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    args = parser.parse_args(argv)
    manifest = run_ml_audit_suite(
        args.dataset,
        args.label,
        args.out_dir,
        feature_set=args.feature_set,
        exclude_features=args.exclude_features,
        splits=args.splits,
        source_col=args.source_col,
        seed=args.seed,
        test_fraction=args.test_fraction,
        allow_partial_rescoring_features=args.allow_partial_rescoring_features,
        allow_label_definition_features=args.allow_label_definition_features,
        candidates=args.candidate,
        strict_source_overlap=args.strict_source_overlap,
        source_transfer_train=args.source_transfer_train,
        source_transfer_test=args.source_transfer_test,
        predictions_path=args.predictions,
        claim_mode=args.claim_mode,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
