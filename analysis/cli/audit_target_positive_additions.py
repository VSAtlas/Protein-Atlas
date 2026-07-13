from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.target_positive_addition_audit import (
    DEFAULT_ACTIVE_NM,
    DEFAULT_FAMILIES,
    DEFAULT_FAMILY_POSITIVE_FLOOR,
    DEFAULT_LABEL,
    DEFAULT_SOURCE_PATHS,
    DEFAULT_TARGET_POSITIVE_FLOOR,
    audit_target_positive_additions,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit sparse Atlas target families for measured-positive additions "
            "from local BindingDB/ChEMBL/Papyrus/TTD/DrugCentral/IUPHAR/cardiac-panel sources."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label", default=DEFAULT_LABEL)
    parser.add_argument("--family", action="append", default=None, help="Target family to audit; repeatable.")
    parser.add_argument("--active-nm", type=float, default=DEFAULT_ACTIVE_NM)
    parser.add_argument("--family-positive-floor", type=int, default=DEFAULT_FAMILY_POSITIVE_FLOOR)
    parser.add_argument("--target-positive-floor", type=int, default=DEFAULT_TARGET_POSITIVE_FLOOR)
    parser.add_argument("--top-per-target", type=int, default=25)
    parser.add_argument(
        "--include-existing-pairs",
        action="store_true",
        help="Keep candidates already represented as target-ligand rows in the dataset.",
    )
    for source_name, default_path in DEFAULT_SOURCE_PATHS.items():
        flag = "--" + source_name.lower().replace("_", "-") + "-file"
        parser.add_argument(flag, type=Path, default=Path(default_path))
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # argparse converts dashes to underscores, while source names preserve case.
    source_paths: dict[str, str | Path] = {
        "BindingDB": args.bindingdb_file,
        "ChEMBL": args.chembl_file,
        "Papyrus": args.papyrus_file,
        "TTD": args.ttd_file,
        "DrugCentral": args.drugcentral_file,
        "IUPHAR_GtoPdb": args.iuphar_gtopdb_file,
        "CardiacSafety": args.cardiacsafety_file,
    }
    manifest = audit_target_positive_additions(
        dataset_path=args.dataset,
        out_dir=args.out_dir,
        label_col=args.label,
        families=args.family or DEFAULT_FAMILIES,
        source_paths=source_paths,
        active_nm=args.active_nm,
        family_positive_floor=args.family_positive_floor,
        target_positive_floor=args.target_positive_floor,
        top_per_target=args.top_per_target,
        include_existing_pairs=args.include_existing_pairs,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
