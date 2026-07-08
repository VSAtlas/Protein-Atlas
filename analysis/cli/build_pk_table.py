from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.pk_sources import build_combined_pk_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a local structured PK/exposure table for Atlas ML joins."
    )
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--existing-pk", type=Path, default=None)
    parser.add_argument("--spd-pk", type=Path, default=None)
    parser.add_argument("--openfda-pk", type=Path, default=None)
    parser.add_argument("--pkdb", type=Path, default=None)
    parser.add_argument("--ncats-inxight-pk", type=Path, default=None)
    parser.add_argument("--drugbank-cmax", type=Path, default=None)
    parser.add_argument("--drugbank-protein-binding", type=Path, default=None)
    build_args = parser.parse_args(argv)
    build_combined_pk_table(
        build_args.out,
        existing_pk=build_args.existing_pk,
        spd_pk=build_args.spd_pk,
        openfda_pk=build_args.openfda_pk,
        pkdb=build_args.pkdb,
        ncats_inxight_pk=build_args.ncats_inxight_pk,
        drugbank_cmax=build_args.drugbank_cmax,
        drugbank_protein_binding=build_args.drugbank_protein_binding,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
