from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.external.mitotox import stage_mitotox


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage MitoTox API tables into repo-local external data.")
    parser.add_argument("--out-dir", type=Path, default=Path("data/external/mitotox"))
    parser.add_argument("--sleep-sec", type=float, default=0.1)
    args = parser.parse_args(argv)
    print(json.dumps(stage_mitotox(args.out_dir, sleep_sec=args.sleep_sec), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
