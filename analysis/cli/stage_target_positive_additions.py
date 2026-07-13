from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.target_positive_addition_prep import (
    DEFAULT_FAMILIES,
    stage_target_positive_additions,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Select and prepare a scaffold-diverse measured-positive ligand add-on "
            "for sparse Atlas target families."
        )
    )
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--library-name", required=True)
    parser.add_argument("--family", action="append", default=None)
    parser.add_argument(
        "--allow-non-fda",
        action="store_true",
        help="Allow externally measured probe compounds; they remain sensitivity-only.",
    )
    parser.add_argument(
        "--pubchem",
        choices=("auto", "never"),
        default="auto",
        help="Resolve missing structures from cached/online PubChem CID properties.",
    )
    parser.add_argument("--retry-buffer", type=int, default=2)
    parser.add_argument("--target-cap", type=int, default=6)
    parser.add_argument("--active-nm", type=float, default=1000.0)
    parser.add_argument("--config", type=Path, default=Path("config.txt"))
    parser.add_argument("--force-prep", action="store_true")
    parser.add_argument("--no-prepare", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = stage_target_positive_additions(
        candidates_path=args.candidates,
        out_dir=args.out_dir,
        library_name=args.library_name,
        families=args.family or DEFAULT_FAMILIES,
        allow_non_fda=args.allow_non_fda,
        resolve_pubchem=args.pubchem == "auto",
        retry_buffer=args.retry_buffer,
        target_cap=args.target_cap,
        active_nm=args.active_nm,
        prepare=not args.no_prepare,
        force=args.force_prep,
        config_path=args.config,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
