from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.adrecs_target import stage_adrecs_drug_target_adr, stage_adrecs_target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Normalize ADReCS-Target target-ADR mechanism positives.")
    parser.add_argument("--associations", required=True, type=Path)
    parser.add_argument("--protein-info", type=Path, default=None)
    parser.add_argument("--adr-info", type=Path, default=None)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--drug-target-adr-out", type=Path, default=None)
    args = parser.parse_args(argv)
    stage_adrecs_target(
        args.associations,
        args.out,
        protein_info_path=args.protein_info,
        adr_info_path=args.adr_info,
    )
    if args.mapping is not None and args.drug_target_adr_out is not None:
        stage_adrecs_drug_target_adr(args.out, args.mapping, args.drug_target_adr_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
