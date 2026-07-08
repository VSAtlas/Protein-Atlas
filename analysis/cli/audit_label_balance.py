from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.label_balance_audit import DEFAULT_GROUP_COLS, audit_label_balance


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit label balance by target, source, family, and context for Atlas ML tables."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--group-cols", nargs="*", default=DEFAULT_GROUP_COLS)
    parser.add_argument("--source-cols", nargs="*", default=["label_source", "source_family", "upstream_source"])
    parser.add_argument(
        "--context-cols",
        nargs="*",
        default=["label_source", "source_family", "upstream_source", "assay_type", "endpoint_type"],
    )
    parser.add_argument("--min-positive", type=int, default=10)
    parser.add_argument("--min-negative", type=int, default=10)
    parser.add_argument("--extreme-rate-low", type=float, default=0.05)
    parser.add_argument("--extreme-rate-high", type=float, default=0.95)
    args = parser.parse_args(argv)
    manifest = audit_label_balance(
        args.dataset,
        label_col=args.label,
        out_dir=args.out_dir,
        group_cols=args.group_cols,
        source_cols=args.source_cols,
        context_cols=args.context_cols,
        min_positive=args.min_positive,
        min_negative=args.min_negative,
        extreme_rate_low=args.extreme_rate_low,
        extreme_rate_high=args.extreme_rate_high,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
