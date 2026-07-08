from __future__ import annotations

import argparse
from pathlib import Path

from analysis.evaluation.profile_comparison import compare_pbas_atlas_profiles


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare PBAS/Sawada and Atlas pair-score profiles.")
    parser.add_argument("--pbas-pairs", required=True, type=Path)
    parser.add_argument("--atlas-pairs", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--atlas-score", default="atlas_score")
    parser.add_argument("--top-k", type=int, default=20)
    args = parser.parse_args(argv)
    compare_pbas_atlas_profiles(args.pbas_pairs, args.atlas_pairs, args.out_dir, atlas_score_col=args.atlas_score, top_k=args.top_k)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
