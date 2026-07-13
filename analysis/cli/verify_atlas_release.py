"""Offline verification command for a generated Docking Atlas site."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.atlas_database.release_bundle import verify_release_bundle


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify a generated Atlas publication bundle offline."
    )
    parser.add_argument("--site-dir", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    report = verify_release_bundle(
        args.site_dir.resolve(),
        report_path=args.report.resolve() if args.report else None,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
