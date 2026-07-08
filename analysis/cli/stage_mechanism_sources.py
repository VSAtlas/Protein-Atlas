from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.external.mechanism_sources import stage_mechanism_sources


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Download and normalize SIDER/CTD/OpenTargets/Reactome mechanism evidence into repo-local data/external paths."
    )
    parser.add_argument("--base-dir", type=Path, default=Path("data/external"))
    parser.add_argument("--pair-table", type=Path, default=None)
    parser.add_argument("--fda-mapping", type=Path, default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"))
    parser.add_argument("--sources", nargs="+", default=["sider", "ctd", "opentargets_safety", "reactome"])
    parser.add_argument("--no-download", action="store_true", help="Normalize already downloaded raw files without network access.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    manifest = stage_mechanism_sources(
        base_dir=args.base_dir,
        pair_table=args.pair_table,
        fda_mapping=args.fda_mapping,
        download=not args.no_download,
        overwrite=args.overwrite,
        sources=args.sources,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
