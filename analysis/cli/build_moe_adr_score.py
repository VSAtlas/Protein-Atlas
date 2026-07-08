from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.moe_evidence import build_moe_adr_score


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a conservative mixture-of-evidence ADR implication score."
    )
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--bioactivity-predictions", type=Path, default=None)
    parser.add_argument("--binding-expert", type=Path, default=None)
    parser.add_argument("--exposure-predictions", type=Path, default=None)
    parser.add_argument("--tissue-predictions", type=Path, default=None)
    parser.add_argument("--mechanism-scores", type=Path, default=None)
    args = parser.parse_args(argv)
    build_moe_adr_score(
        args.pair_table,
        args.out,
        bioactivity_predictions=args.bioactivity_predictions,
        binding_expert=args.binding_expert,
        exposure_predictions=args.exposure_predictions,
        tissue_predictions=args.tissue_predictions,
        mechanism_scores=args.mechanism_scores,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
