from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.external.hpa_gtex import stage_tissue_expression_sources


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Download and normalize HPA/GTEx/Bgee/OpenTargets expression context for Atlas targets."
    )
    parser.add_argument(
        "--target-table",
        type=Path,
        default=Path("data/spd_full_panel_latest/spd_full_panel_relaxed/spd_target_selected.csv"),
        help="Target table with a gene/target_gene column. Defaults to the relaxed full SPD target selection.",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("data/external/hpa_gtex"))
    parser.add_argument("--no-download", action="store_true", help="Use already staged raw files only.")
    parser.add_argument(
        "--download-opentargets-baseline",
        action="store_true",
        help="Download/filter Open Targets baseline-expression parquet shards. This is several GB.",
    )
    parser.add_argument(
        "--opentargets-part-limit",
        type=int,
        default=0,
        help="Debug limit for Open Targets parquet shards. 0 means all discovered shards.",
    )
    parser.add_argument(
        "--opentargets-safety",
        type=Path,
        default=Path("data/external/opentargets/safety.tsv"),
        help="Optional local Open Targets Safety table used as site-safety context, not expression truth.",
    )
    args = parser.parse_args(argv)
    manifest = stage_tissue_expression_sources(
        target_table_path=args.target_table,
        out_dir=args.out_dir,
        download=not args.no_download,
        download_opentargets=args.download_opentargets_baseline,
        opentargets_part_limit=max(0, args.opentargets_part_limit),
        opentargets_safety_path=args.opentargets_safety,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
