from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.paper_dti_sources import stage_paper_dti_sources


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stage downloaded ADRtarget/T-ARDIS/DTIAM/MoseDTI/EnsDTI-style sources into Atlas-compatible audit tables."
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=None,
        help="Atlas run data directory used to resolve default pair/mapping inputs.",
    )
    parser.add_argument("--ensdti-root", type=Path, default=Path("data/external/paper_sources/ensdti/EnsDTI-master"))
    parser.add_argument("--mosedti-root", type=Path, default=Path("data/external/paper_sources/mosedti/MoseDTI-main"))
    parser.add_argument("--adrtarget-root", type=Path, default=Path("data/external/paper_sources/adrtarget/ADRtarget-master"))
    parser.add_argument("--tardis-root", type=Path, default=Path("data/external/paper_sources/tardis/T-ARDIS-master"))
    parser.add_argument("--doctor-root", type=Path, default=Path("data/external/paper_sources/doctor/DocTOR-master"))
    parser.add_argument("--dtiam-root", type=Path, default=Path("data/external/paper_sources/dtiam/github/DTIAM-main"))
    parser.add_argument("--bindingdb-bioactivity", type=Path, default=Path("data/external/bindingdb/bioactivity.tsv"))
    parser.add_argument(
        "--pair-table",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--atlas-compound-map",
        type=Path,
        default=None,
    )
    args = parser.parse_args(argv)
    run_dir = args.run_dir or args.out_dir.parent
    pair_table = (
        args.pair_table
        or run_dir / "moe_experts" / "banana_four_state_with_scorch_backfill_model_ready.csv"
    )
    atlas_compound_map = (
        args.atlas_compound_map
        or run_dir / "source_overlap_deep_sources_v3" / "spd_drug_atlas_mapping_audit.csv"
    )
    stage_paper_dti_sources(
        args.out_dir,
        ensdti_root=args.ensdti_root,
        mosedti_root=args.mosedti_root,
        adrtarget_root=args.adrtarget_root,
        tardis_root=args.tardis_root,
        doctor_root=args.doctor_root,
        dtiam_root=args.dtiam_root,
        bindingdb_bioactivity_path=args.bindingdb_bioactivity,
        pair_table_path=pair_table,
        atlas_compound_map_path=atlas_compound_map,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
