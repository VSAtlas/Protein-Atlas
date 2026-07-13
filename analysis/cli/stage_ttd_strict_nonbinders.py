from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.ttd_nonbinder_additions import stage_ttd_strict_nonbinders


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Stage strict TTD >200 uM same-source nonbinders as external "
            "sensitivity controls; this command never changes SPD labels."
        )
    )
    parser.add_argument("--ttd", type=Path, required=True)
    parser.add_argument("--positive-pairs", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--library-name", default="ttd_strict_nonbinders")
    parser.add_argument("--nonbinder-nm", type=float, default=200_000.0)
    parser.add_argument("--active-nm", type=float, default=1_000.0)
    parser.add_argument("--negatives-per-positive", type=int, default=2)
    parser.add_argument("--max-per-target", type=int, default=10)
    parser.add_argument(
        "--resolve-pubchem",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--config", type=Path, default=Path("config.txt"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = stage_ttd_strict_nonbinders(
        ttd_path=args.ttd,
        positive_pairs_path=args.positive_pairs,
        out_dir=args.out_dir,
        library_name=args.library_name,
        nonbinder_nm=args.nonbinder_nm,
        active_nm=args.active_nm,
        negatives_per_positive=args.negatives_per_positive,
        max_per_target=args.max_per_target,
        resolve_pubchem=args.resolve_pubchem,
        prepare=args.prepare,
        force=args.force,
        config_path=args.config,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
