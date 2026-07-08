from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.tissue_site_scoring import build_tissue_site_scores


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build deterministic target tissue/site relevance scores from expression evidence."
    )
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--target-expression", type=Path, default=None)
    args = parser.parse_args(argv)
    build_tissue_site_scores(
        args.pair_table,
        args.out,
        target_expression_path=args.target_expression,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
