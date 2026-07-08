from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.adr_mechanism_panel import build_adr_mechanism_panel_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build drug-ADR keyed mechanism-control ML table.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--four-state-labels", required=True, type=Path)
    parser.add_argument("--mechanism-scores", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    build_adr_mechanism_panel_table(
        args.pair_table,
        args.four_state_labels,
        args.out,
        mechanism_scores_path=args.mechanism_scores,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
