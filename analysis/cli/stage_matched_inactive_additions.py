from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.matched_inactive_additions import (
    parse_source_spec,
    stage_matched_inactive_additions,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Stage censor-aware measured inactive FDA/RDK pairs matched to "
            "same-target SPD positives. Uses local sources only and does not launch docking."
        )
    )
    parser.add_argument(
        "--dataset",
        "--spd-dataset",
        dest="spd_dataset",
        required=True,
        type=Path,
        help="Corrected SPD model-ready CSV/TSV/parquet table.",
    )
    parser.add_argument(
        "--source",
        "--source-path",
        dest="sources",
        action="append",
        required=True,
        metavar="[NAME=]PATH",
        help=(
            "Local normalized ChEMBL, Papyrus, ToxCast, or BindingDB table. "
            "Repeat for each source."
        ),
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label-col", default=None)
    parser.add_argument("--max-per-target", type=int, default=10)
    parser.add_argument("--active-nm", type=float, default=1_000.0)
    parser.add_argument("--inactive-nm", type=float, default=10_000.0)
    parser.add_argument("--min-tanimoto", type=float, default=0.0)
    parser.add_argument("--max-descriptor-distance", type=float, default=None)
    parser.add_argument("--descriptor-weight", type=float, default=0.25)
    parser.add_argument(
        "--include-existing",
        action="store_true",
        help="Retain existing non-positive SPD drug-target pairs as candidates.",
    )
    parser.add_argument("--library-name", default="matched_inactive_additions")
    parser.add_argument(
        "--prepare-pdbqt",
        action="store_true",
        help="Run local Meeko preparation after SDF generation; never launches docking.",
    )
    parser.add_argument("--force-prep", action="store_true")
    parser.add_argument("--config", type=Path, default=Path("config.txt"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = stage_matched_inactive_additions(
        spd_dataset_path=args.spd_dataset,
        source_paths=[parse_source_spec(value) for value in args.sources],
        out_dir=args.out_dir,
        max_per_target=args.max_per_target,
        label_col=args.label_col,
        active_nm=args.active_nm,
        inactive_nm=args.inactive_nm,
        min_tanimoto=args.min_tanimoto,
        max_descriptor_distance=args.max_descriptor_distance,
        descriptor_weight=args.descriptor_weight,
        include_existing=args.include_existing,
        library_name=args.library_name,
        prepare=args.prepare_pdbqt,
        force=args.force_prep,
        config_path=args.config,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
