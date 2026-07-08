from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.pilot_bioactivity_sources import stage_pilot_bioactivity_sources


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build focused Papyrus/ChEMBL/ToxCast bioactivity extracts for an Atlas pilot run."
    )
    parser.add_argument("--pilot-dir", required=True, type=Path)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument(
        "--external-dir",
        type=Path,
        default=None,
        help="Directory for downloaded source archives/extracts. Defaults to <repo-root>/data/external.",
    )
    parser.add_argument(
        "--pk-table",
        type=Path,
        default=None,
        help="Optional dedicated PK table with drug_id/free_cmax_um columns. Missing exposure is not inferred.",
    )
    parser.add_argument(
        "--fda-mapping",
        type=Path,
        default=None,
        help="Optional FDA/RDK mapping CSV override. Defaults to repo-root mapping discovery.",
    )
    parser.add_argument("--overwrite-downloads", action="store_true")
    args = parser.parse_args(argv)
    outputs = stage_pilot_bioactivity_sources(
        args.pilot_dir,
        args.repo_root,
        external_dir=args.external_dir,
        fda_mapping_path=args.fda_mapping,
        pk_table_path=args.pk_table,
        overwrite_downloads=args.overwrite_downloads,
    )
    for key, value in outputs.items():
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
