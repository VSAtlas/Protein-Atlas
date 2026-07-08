from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.pubchem_bioassay import stage_pubchem_aid_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage a downloaded PubChem BioAssay AID CSV for Atlas evidence.")
    parser.add_argument("--aid", required=True, type=int)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--summary-json", type=Path, default=None)
    parser.add_argument("--mapping", type=Path, default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"))
    parser.add_argument("--cid-properties-cache", type=Path, default=None)
    parser.add_argument("--drop-unmapped", action="store_true")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    stage_pubchem_aid_table(
        args.csv,
        args.mapping,
        args.out,
        aid=args.aid,
        summary_json_path=args.summary_json,
        cid_properties_cache=args.cid_properties_cache,
        keep_unmapped=not args.drop_unmapped,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
