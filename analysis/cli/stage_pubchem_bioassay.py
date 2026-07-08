from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.pubchem_bioassay import stage_pubchem_bioassay_for_targets


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage targeted PubChem BioAssay active/inactive outcomes for Atlas targets.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"))
    parser.add_argument("--out", type=Path, default=Path("data/external/pubchem_bioassay/pubchem_bioassay.tsv"))
    parser.add_argument("--max-aids-per-target", type=int, default=25)
    parser.add_argument("--sleep-sec", type=float, default=0.2)
    parser.add_argument("--resolve-cid-inchikey", action="store_true")
    parser.add_argument("--cid-properties-cache", type=Path, default=None)
    args = parser.parse_args(argv)
    stage_pubchem_bioassay_for_targets(
        args.pair_table,
        args.mapping,
        args.out,
        max_aids_per_target=args.max_aids_per_target,
        sleep_sec=args.sleep_sec,
        resolve_cid_inchikey=args.resolve_cid_inchikey,
        cid_properties_cache=args.cid_properties_cache,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
