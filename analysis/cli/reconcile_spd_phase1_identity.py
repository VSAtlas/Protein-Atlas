"""CLI for fail-closed Phase 1 SPD identity reconciliation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from analysis.ml.spd_identity_reconciliation import (
    DEFAULT_IDENTITY_KEY_COLUMNS,
    DEFAULT_SPD_INCHIKEY_COLUMNS,
    IdentityReconciliationError,
    default_fda_mapping_path,
    reconcile_spd_phase1_identity,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reconcile Phase 1 SPD rows to canonical FDA identities by exact "
            "normalized RDK key, with structure/hash auditing and fail-closed output."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="Phase 1 SPD CSV.")
    parser.add_argument("--out", required=True, type=Path, help="Reconciled CSV.")
    parser.add_argument(
        "--mapping",
        type=Path,
        default=default_fda_mapping_path(),
        help="Canonical FDA mapping CSV (default: chemdb/data/fda_mapping_from_pdbqt.csv).",
    )
    parser.add_argument(
        "--row-audit",
        type=Path,
        default=None,
        help="Row audit CSV (default: <out stem>_identity_row_audit.csv).",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Summary JSON (default: <out stem>_identity_summary.json).",
    )
    parser.add_argument(
        "--identity-key-column",
        action="append",
        default=None,
        help=(
            "Exact RDK key column; repeat to cross-check aliases. Defaults to "
            + ", ".join(DEFAULT_IDENTITY_KEY_COLUMNS)
            + "."
        ),
    )
    parser.add_argument(
        "--spd-inchikey-column",
        action="append",
        default=None,
        help=(
            "SPD structure InChIKey column; repeat to cross-check aliases. Defaults to "
            + ", ".join(DEFAULT_SPD_INCHIKEY_COLUMNS)
            + "."
        ),
    )
    parser.add_argument(
        "--additional-spd-null-column",
        action="append",
        default=[],
        help=(
            "Additional SPD-owned value to null only on demonstrated structure/SPD "
            "mismatch; repeat as needed."
        ),
    )
    parser.add_argument(
        "--allow-blocked-output",
        action="store_true",
        help=(
            "Return success even when the summary remains release_status=blocked. "
            "The default returns exit code 2 for any blocking row."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    fail_closed = not args.allow_blocked_output
    try:
        summary = reconcile_spd_phase1_identity(
            args.input,
            args.out,
            mapping_csv=args.mapping,
            row_audit_csv=args.row_audit,
            summary_json=args.summary,
            identity_key_columns=args.identity_key_column,
            spd_inchikey_columns=args.spd_inchikey_column,
            additional_spd_null_columns=args.additional_spd_null_column,
            fail_closed=fail_closed,
        )
    except IdentityReconciliationError as exc:
        print(f"identity reconciliation failed: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, indent=2, sort_keys=True))
    return int(summary["recommended_exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
