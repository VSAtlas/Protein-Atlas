from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.leaderboard import write_ml_leaderboard


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build a run-scoped Atlas ML leaderboard from manifests.")
    parser.add_argument("--ml-root", required=True, type=Path)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    manifest = write_ml_leaderboard(args.ml_root, args.out_dir)
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
