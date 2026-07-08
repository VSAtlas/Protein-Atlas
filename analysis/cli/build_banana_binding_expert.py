from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.banana import build_banana_binding_expert


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a BANANA/BigBind binding expert table for Atlas MoE scoring."
    )
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--banana-scores", type=Path, default=None)
    args = parser.parse_args(argv)
    build_banana_binding_expert(
        args.pair_table,
        args.out,
        banana_scores_path=args.banana_scores,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
