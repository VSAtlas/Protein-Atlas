from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.family_panels import stage_kiba_family_panel


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Normalize family-specific bioactivity panels into Atlas evidence tables.")
    parser.add_argument("--kiba", type=Path, default=Path("data/external/family_panels/kiba.zip"))
    parser.add_argument("--mapping", type=Path, default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"))
    parser.add_argument("--chembl-mapping", type=Path, default=None)
    parser.add_argument("--resolve-chembl", action="store_true")
    parser.add_argument("--sleep-sec", type=float, default=0.05)
    parser.add_argument("--out", type=Path, default=Path("data/external/family_panels/kiba_bioactivity.tsv"))
    args = parser.parse_args(argv)
    stage_kiba_family_panel(
        args.kiba,
        args.mapping,
        args.out,
        chembl_mapping_path=args.chembl_mapping,
        resolve_chembl=args.resolve_chembl,
        sleep_sec=args.sleep_sec,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
