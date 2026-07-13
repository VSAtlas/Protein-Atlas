from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.materialize_final_score import (
    DEFAULT_EXISTING_SCORE_COLUMN,
    DEFAULT_FALLBACK_COLUMN,
    materialize_final_score,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Materialize a provenance-controlled final score from an existing "
            "post-docked score or a comparison-run standardized z-score."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument(
        "--existing-score-column",
        default=DEFAULT_EXISTING_SCORE_COLUMN,
        help="Existing score to preserve only with accepted provenance.",
    )
    parser.add_argument(
        "--fallback-column",
        default=DEFAULT_FALLBACK_COLUMN,
        help="Finite standardized z-score used when no accepted existing score remains.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = materialize_final_score(
        args.dataset,
        args.out,
        args.audit_dir,
        existing_score_column=args.existing_score_column,
        fallback_column=args.fallback_column,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=True, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
