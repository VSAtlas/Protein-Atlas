from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.source_transfer_audit import audit_source_transfer


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit Atlas ML source-transfer failure modes.")
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--train-sources", nargs="+", required=True)
    parser.add_argument("--test-sources", nargs="+", required=True)
    parser.add_argument("--predictions", type=Path, default=None)
    parser.add_argument("--features", nargs="*", default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    audit_source_transfer(
        args.dataset,
        args.label,
        args.out_dir,
        train_sources=args.train_sources,
        test_sources=args.test_sources,
        predictions_path=args.predictions,
        feature_cols=args.features,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
