from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.pubchem_spd_overlap import audit_pubchem_aid_spd_overlap


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit normalized PubChem AID rows against the full local SPD workbook.")
    parser.add_argument("--pubchem-aid-table", required=True, type=Path)
    parser.add_argument("--spd-workbook", required=True, type=Path)
    parser.add_argument("--target-map", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    audit_pubchem_aid_spd_overlap(
        args.pubchem_aid_table,
        args.spd_workbook,
        args.out_dir,
        target_map_path=args.target_map,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
