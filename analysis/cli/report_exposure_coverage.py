from __future__ import annotations

import argparse
from pathlib import Path

from analysis.exposure_coverage import write_exposure_coverage_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Report free-Cmax/Cmax/fraction-unbound coverage without imputation.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--group-cols", nargs="*", default=None)
    args = parser.parse_args(argv)
    write_exposure_coverage_report(args.input, args.out_dir, group_cols=args.group_cols)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
