from __future__ import annotations

import argparse
from pathlib import Path

from analysis.evaluation.sawada_vs_atlas import build_sawada_vs_atlas_summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build a manuscript-style Sawada/PBAS vs Atlas comparison table.")
    parser.add_argument("--pbas-pairs", required=True, type=Path)
    parser.add_argument("--atlas-pairs", required=True, type=Path)
    parser.add_argument("--labels", type=Path, default=None)
    parser.add_argument("--spd", type=Path, default=None)
    parser.add_argument("--sawada-metrics", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    build_sawada_vs_atlas_summary(args.pbas_pairs, args.atlas_pairs, args.out, sawada_metrics_path=args.sawada_metrics, spd_path=args.spd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
