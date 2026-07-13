from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.merge_matched_inactive_scores import merge_matched_inactive_scores


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Strictly append audited, fully scored matched-inactive rows to an "
            "ML table without modifying any base row."
        )
    )
    parser.add_argument(
        "--base-table",
        required=True,
        type=Path,
        help="Authoritative base ML table; every row and nonmissing value is preserved.",
    )
    parser.add_argument(
        "--selected-pairs",
        required=True,
        type=Path,
        help="Audited selected_matched_inactive_pairs.csv artifact.",
    )
    parser.add_argument(
        "--candidate-scores",
        required=True,
        type=Path,
        help=(
            "Candidate score table containing consensus and SCORCH comparison-run "
            "z-scores with provenance."
        ),
    )
    parser.add_argument("--out", required=True, type=Path, help="Merged output CSV.")
    parser.add_argument(
        "--audit-dir",
        type=Path,
        default=None,
        help="Audit directory (default: <out stem>_audit beside --out).",
    )
    parser.add_argument(
        "--allow-base-pair-collisions",
        action="store_true",
        help=(
            "Explicitly append candidate rows whose normalized pair exists in the "
            "base table. Base rows are still never overwritten."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = merge_matched_inactive_scores(
        args.base_table,
        args.selected_pairs,
        args.candidate_scores,
        args.out,
        audit_dir=args.audit_dir,
        allow_base_pair_collisions=args.allow_base_pair_collisions,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
