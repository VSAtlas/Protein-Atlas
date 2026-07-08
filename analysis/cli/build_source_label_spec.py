from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.source_label_spec import build_source_label_spec


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Write Atlas source verification, label-policy, and leakage-control specification tables."
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    build_source_label_spec(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
