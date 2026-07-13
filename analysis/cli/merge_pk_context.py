from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.merge_pk_context import DEFAULT_KEY_COLUMNS, merge_pk_context


def _parse_key(value: str) -> tuple[str, ...]:
    columns = tuple(column.strip() for column in value.split(",") if column.strip())
    if not columns:
        raise argparse.ArgumentTypeError(
            "key must be a comma-separated list of column names"
        )
    return columns


def build_parser() -> argparse.ArgumentParser:
    default_keys = "; ".join("+".join(columns) for columns in DEFAULT_KEY_COLUMNS)
    parser = argparse.ArgumentParser(
        description=(
            "Losslessly backfill PK/context fields from an enriched SPD CSV into "
            "a corrected primary model-ready CSV."
        )
    )
    parser.add_argument(
        "--primary",
        required=True,
        type=Path,
        help="Corrected primary SPD model-ready CSV whose rows and order are authoritative.",
    )
    parser.add_argument(
        "--pk-enriched",
        required=True,
        type=Path,
        help="SPD PK-enriched CSV used only as a backfill source.",
    )
    parser.add_argument("--out", required=True, type=Path, help="Merged output CSV.")
    parser.add_argument(
        "--audit-dir",
        type=Path,
        default=None,
        help="Audit output directory (default: <out stem>_audit beside --out).",
    )
    parser.add_argument(
        "--column",
        dest="columns",
        action="append",
        default=None,
        help=(
            "Exact source column to backfill; repeat as needed. By default, "
            "pk_context_* and direct PK measurement columns are selected."
        ),
    )
    parser.add_argument(
        "--key",
        dest="keys",
        action="append",
        type=_parse_key,
        default=None,
        metavar="COL1,COL2",
        help=(
            "Ordered composite key tier; repeat from strongest to weakest. "
            f"Defaults: {default_keys}"
        ),
    )
    parser.add_argument(
        "--overwrite-conflicts",
        action="store_true",
        help=(
            "Explicitly replace conflicting nonmissing primary cells with source "
            "values; all replacements remain listed in conflicts.csv."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = merge_pk_context(
        args.primary,
        args.pk_enriched,
        args.out,
        audit_dir=args.audit_dir,
        backfill_columns=args.columns,
        key_columns=args.keys,
        overwrite_conflicts=args.overwrite_conflicts,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
