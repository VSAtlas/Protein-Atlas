from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.family_panels import stage_cardiac_ion_channel_panel


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Normalize extracted cardiac ion-channel panel CSVs into Atlas evidence tables.")
    parser.add_argument("--extracted-root", type=Path, default=Path("data/external/cardiac_safety/raw_extracted/raw"))
    parser.add_argument("--mapping", type=Path, default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"))
    parser.add_argument("--out", type=Path, default=Path("data/external/cardiac_safety/cardiac_ion_channel_bioactivity.tsv"))
    args = parser.parse_args(argv)
    stage_cardiac_ion_channel_panel(args.extracted_root, args.mapping, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
