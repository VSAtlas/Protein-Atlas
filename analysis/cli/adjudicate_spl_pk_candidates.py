from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.external.spl_pk_adjudication import adjudicate_spl_pk_candidates


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Adjudicate cached SPL clearance, absolute-bioavailability, and "
            "maximum-recommended-dose candidates with same-clause semantic gates."
        )
    )
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--review-decisions",
        type=Path,
        default=None,
        help=(
            "Optional canonical-candidate review CSV. Generate and edit "
            "review_decisions_template.csv before supplying it here."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.candidates.is_file():
        raise FileNotFoundError(f"candidate CSV not found: {args.candidates}")
    if args.review_decisions is not None and not args.review_decisions.is_file():
        raise FileNotFoundError(
            f"review decision CSV not found: {args.review_decisions}"
        )
    manifest = adjudicate_spl_pk_candidates(
        candidates=args.candidates,
        out_dir=args.out_dir,
        review_decisions=args.review_decisions,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
