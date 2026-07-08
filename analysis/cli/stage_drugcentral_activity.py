from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.drugcentral import stage_drugcentral_activity


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage DrugCentral target interactions for four-state labels.")
    parser.add_argument("--interactions", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    stage_drugcentral_activity(args.interactions, args.mapping, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
