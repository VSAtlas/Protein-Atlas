from __future__ import annotations

import argparse
import os
from pathlib import Path

from analysis.external.faers import build_faers_disproportionality_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a cached openFDA FAERS disproportionality table for drug-ADR pairs."
    )
    parser.add_argument("--pairs", nargs="+", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--negative-only", action="store_true")
    parser.add_argument("--api-key", default=None)
    parser.add_argument("--api-key-env", default="OPENFDA_API_KEY")
    parser.add_argument("--sleep-sec", type=float, default=0.2)
    parser.add_argument("--max-pairs", type=int, default=None)
    args = parser.parse_args(argv)
    api_key = args.api_key or os.environ.get(args.api_key_env)
    build_faers_disproportionality_table(
        args.pairs,
        args.out,
        cache_dir=args.cache_dir,
        negative_only=args.negative_only,
        api_key=api_key,
        sleep_sec=args.sleep_sec,
        max_pairs=args.max_pairs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
