"""Build evidence-only receptor annotations for a docking-atlas release."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.atlas_database.receptor_evidence import write_receptor_evidence


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect exact prepared receptors, retain structured metal audits, and "
            "write APO/HOLO evidence without qualifying any receptor."
        )
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--expected-context-count", type=int)
    parser.add_argument(
        "--expected-context-key",
        action="append",
        help="Repeat exact run_id|PDB|VARIANT|pH_or_base inventory keys.",
    )
    parser.add_argument(
        "--expected-inventory-sha256",
        help=(
            "Require the SHA-256 of the sorted exact receptor-context inventory."
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.expected_context_count is not None and args.expected_context_count <= 0:
        raise SystemExit("--expected-context-count must be positive")
    summary = write_receptor_evidence(
        args.manifest,
        args.repo_root,
        args.out,
        expected_context_count=args.expected_context_count,
        expected_context_keys=args.expected_context_key,
        expected_inventory_sha256=args.expected_inventory_sha256,
        overwrite=args.overwrite,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
