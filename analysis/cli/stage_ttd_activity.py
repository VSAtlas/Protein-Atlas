from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.ttd import stage_ttd_activity


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage TTD target-compound activity for four-state labels.")
    parser.add_argument("--activity", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--crossmatch", required=True, type=Path)
    parser.add_argument("--structures", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--feature-table", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    stage_ttd_activity(
        args.activity,
        args.targets,
        args.crossmatch,
        args.structures,
        args.mapping,
        args.out,
        feature_table_path=args.feature_table,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
