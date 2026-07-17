"""Command-line wrapper for the non-mutating ligand stereo audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from prep_ligands.ligand_stereo_audit import (
    LigandStereoAuditConfig,
    LigandStereoSource,
    audit_ligand_stereochemistry,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="atlas ligands audit stereo",
        description=(
            "Report specified and unspecified potential stereochemistry from "
            "source SDF/CSV records; never assign missing stereo."
        ),
    )
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sdf", action="append", default=[], type=Path)
    parser.add_argument("--sdf-id-property", default="_Name")
    parser.add_argument("--sdf-parent-id-property")
    parser.add_argument("--csv", action="append", default=[], type=Path)
    parser.add_argument("--csv-id-column")
    parser.add_argument("--csv-structure-column")
    parser.add_argument(
        "--csv-structure-format",
        choices=("smiles", "inchi", "molblock"),
        default="smiles",
    )
    parser.add_argument("--csv-parent-id-column")
    parser.add_argument("--csv-prepared-path-column")
    parser.add_argument("--csv-prepared-sha256-column")
    parser.add_argument("--csv-prepared-path-root", type=Path)
    parser.add_argument("--fail-on-parse-error", action="store_true")
    parser.add_argument("--fail-on-unspecified", action="store_true")
    parser.add_argument("--fail-on-unresolved", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if not args.sdf and not args.csv:
        parser.error("at least one --sdf or --csv input is required")
    if args.csv and (not args.csv_id_column or not args.csv_structure_column):
        parser.error("CSV inputs require --csv-id-column and --csv-structure-column")

    sources = [
        LigandStereoSource(
            path=path,
            kind="sdf",
            id_field=args.sdf_id_property,
            parent_id_field=args.sdf_parent_id_property,
        )
        for path in args.sdf
    ]
    sources.extend(
        LigandStereoSource(
            path=path,
            kind="csv",
            id_field=args.csv_id_column,
            parent_id_field=args.csv_parent_id_column,
            structure_field=args.csv_structure_column,
            structure_format=args.csv_structure_format,
            prepared_path_field=args.csv_prepared_path_column,
            prepared_sha256_field=args.csv_prepared_sha256_column,
            prepared_path_root=args.csv_prepared_path_root,
        )
        for path in args.csv
    )
    result = audit_ligand_stereochemistry(
        LigandStereoAuditConfig(
            library_id=args.library_id,
            output_dir=args.out_dir,
            sources=tuple(sources),
        )
    )
    print(
        json.dumps(
            {
                "evidence_jsonl": str(result.evidence_jsonl),
                "parse_failed_records": result.parse_failed_records,
                "parsed_records": result.parsed_records,
                "records_csv": str(result.records_csv),
                "records_with_unspecified_stereo": result.records_with_unspecified_stereo,
                "records_with_unresolved_stereo": result.records_with_unresolved_stereo,
                "summary_json": str(result.summary_json),
                "total_records": result.total_records,
            },
            sort_keys=True,
        )
    )
    if args.fail_on_parse_error and result.parse_failed_records:
        return 2
    if args.fail_on_unspecified and result.records_with_unspecified_stereo:
        return 3
    if args.fail_on_unresolved and result.records_with_unresolved_stereo:
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
