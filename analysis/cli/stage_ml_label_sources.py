from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.ml_source_catalog import stage_ml_source_catalog


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage or audit external ML label source files.")
    parser.add_argument("--out-dir", default=Path("data/external/ml_label_sources"), type=Path)
    parser.add_argument("--download-small", action="store_true")
    parser.add_argument("--download-large", action="store_true")
    parser.add_argument("--sources", nargs="*", default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    stage_ml_source_catalog(
        args.out_dir,
        download_small=args.download_small,
        download_large=args.download_large,
        sources=args.sources,
        overwrite=args.overwrite,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
