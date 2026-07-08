from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.mechanism_targets import build_four_state_mechanism_labels


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build four-state Atlas mechanism labels.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--mechanism-scores", type=Path, default=None)
    parser.add_argument("--negative-evidence", type=Path, default=None)
    parser.add_argument("--positive-evidence", type=Path, default=None)
    args = parser.parse_args(argv)
    build_four_state_mechanism_labels(
        args.pair_table,
        args.out,
        mechanism_scores_path=args.mechanism_scores,
        negative_evidence_path=args.negative_evidence,
        positive_evidence_path=args.positive_evidence,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
