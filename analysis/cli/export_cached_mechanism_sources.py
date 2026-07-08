from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.mechanism_cache_sources import export_cached_mechanism_sources


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Export Atlas report caches into graph-ready mechanism source TSVs.")
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    export_cached_mechanism_sources(args.cache_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
