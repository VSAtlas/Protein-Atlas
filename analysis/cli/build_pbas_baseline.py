from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.pbas import build_pbas_pair_scores


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Atlas-compatible Sawada/PBAS pair-score baseline tables.")
    parser.add_argument("--pbas-matrix", required=True, type=Path)
    parser.add_argument("--pbas-drug-metadata", type=Path, default=None)
    parser.add_argument("--pbas-target-metadata", type=Path, default=None)
    parser.add_argument("--pbas-pocket-info", type=Path, default=None)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    build_pbas_pair_scores(
        args.pbas_matrix,
        args.out,
        drug_metadata_path=args.pbas_drug_metadata,
        target_metadata_path=args.pbas_target_metadata,
        pocket_info_path=args.pbas_pocket_info,
        mapping_path=args.mapping,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
