from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.external.spd_cmax_context import (
    SMIT_COMMIT,
    SPD_WORKBOOK_SHA256,
    recover_spd_cmax_context,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Recover auditable Cmax-source provenance from the SPD supplement and "
            "the pinned public Smit source release. Aggregate sources remain "
            "aggregate-only and are never assigned synthetic study context."
        )
    )
    parser.add_argument(
        "--spd-workbook",
        type=Path,
        default=Path(
            "data/external/spd/sutherland_2023_spd_supplementary_data_1_15.xlsx"
        ),
    )
    parser.add_argument(
        "--spd-workbook-sha256",
        default=SPD_WORKBOOK_SHA256,
        help="Expected SHA-256 for the exact SPD workbook release.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/external/spd/cmax_source_recovery"),
    )
    parser.add_argument(
        "--smit-source-dir",
        type=Path,
        default=Path("data/external/spd") / f"smit_cmax_source_{SMIT_COMMIT[:12]}",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Require already-cached pinned Smit files instead of downloading them.",
    )
    parser.add_argument("--value-tolerance", type=float, default=0.01)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.value_tolerance < 0:
        raise ValueError("--value-tolerance must be non-negative")
    manifest = recover_spd_cmax_context(
        workbook=args.spd_workbook,
        out_dir=args.out_dir,
        source_dir=args.smit_source_dir,
        download=not args.offline,
        value_tolerance=args.value_tolerance,
        expected_workbook_sha256=args.spd_workbook_sha256,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
