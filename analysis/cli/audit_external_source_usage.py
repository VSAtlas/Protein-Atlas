from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.source_usage_audit import audit_external_source_usage


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Classify local data/external delimited files by Atlas training use.")
    parser.add_argument("--external-dir", default=Path("data/external"), type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    audit_external_source_usage(args.external_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
