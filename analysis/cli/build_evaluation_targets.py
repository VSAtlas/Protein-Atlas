from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.evaluation_targets import build_evaluation_targets


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build separated Atlas ML/evaluation target tables."
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--pair-table", type=Path, default=None)
    parser.add_argument("--bioactivity-table", type=Path, default=None)
    parser.add_argument("--spd-table", type=Path, default=None)
    parser.add_argument("--spd-full-panel", type=Path, default=None)
    parser.add_argument("--spd-mapping", type=Path, default=None)
    parser.add_argument("--mechanism-scores", type=Path, default=None)
    parser.add_argument("--mechanism-edges", type=Path, default=None)
    parser.add_argument("--bigbind-table", type=Path, default=None)
    parser.add_argument("--bigbind-mapping", type=Path, default=None)
    args = parser.parse_args(argv)
    build_evaluation_targets(
        args.out_dir,
        pair_table=args.pair_table,
        bioactivity_table=args.bioactivity_table,
        spd_table=args.spd_table,
        spd_full_panel=args.spd_full_panel,
        spd_mapping=args.spd_mapping,
        mechanism_scores=args.mechanism_scores,
        mechanism_edges=args.mechanism_edges,
        bigbind_table=args.bigbind_table,
        bigbind_mapping=args.bigbind_mapping,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
