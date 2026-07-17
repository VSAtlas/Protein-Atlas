"""CLI for the post-run frozen-SPD holo/preparation integrity audit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from analysis.reporting.spd_holo_integrity import (
    DEFAULT_SELECTION,
    DEFAULT_STRICT_MANIFEST,
    HoloIntegrityInputError,
    audit_spd_holo_integrity,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit the frozen 93-target SPD cohort after protein preparation using "
            "run-scoped structured artifacts from one manifest-bound ligand, "
            "receptor, Meeko, and retention layout."
        )
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument(
        "--strict-manifest", type=Path, default=DEFAULT_STRICT_MANIFEST
    )
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--strict", dest="strict", action="store_true", default=True)
    parser.add_argument("--no-strict", dest="strict", action="store_false")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = audit_spd_holo_integrity(
            repo_root=args.repo_root,
            run_id=args.run_id,
            selection_csv=args.selection,
            strict_manifest_csv=args.strict_manifest,
            output_dir=args.out_dir,
            strict=bool(args.strict),
        )
    except HoloIntegrityInputError as exc:
        print(f"spd holo integrity: {exc}", file=sys.stderr)
        return 2
    print("spd_holo_integrity: complete")
    print(f"counts: {json.dumps(result.counts, sort_keys=True)}")
    print(f"csv: {result.csv_path}")
    print(f"json: {result.json_path}")
    if result.exit_code:
        print("certification: incomplete", file=sys.stderr)
    else:
        print("certification: complete")
    return result.exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["main"]
