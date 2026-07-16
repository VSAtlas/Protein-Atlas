from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.external.pkdb_recovery import (
    build_canonical_pk_identity_table,
    recover_pkdb_context,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Recover Phase 1 PK-DB context through public study source TSVs when "
            "the documented output export is unavailable. Public PK-DB access does "
            "not grant original-source ML-training rights."
        )
    )
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--fda-mapping",
        type=Path,
        default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"),
        help=(
            "Validated canonical ligand-base identity mapping. Set to an empty "
            "path only for a non-FDA table without ligand_base identifiers."
        ),
    )
    parser.add_argument(
        "--source-rights-manifest",
        "--source-rights",
        dest="source_rights_manifest",
        type=Path,
        default=None,
        help=(
            "Audited CSV with study_sid, training_allowed, and rights_reference; "
            "study_reference/source_file_url may narrow a grant."
        ),
    )
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--max-studies", type=int, default=0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not args.dataset.is_file():
        raise FileNotFoundError(f"dataset does not exist: {args.dataset}")
    if args.workers < 1 or args.workers > 8:
        raise ValueError("--workers must be between 1 and 8")
    model_table = pd.read_csv(args.dataset, low_memory=False)
    if args.fda_mapping is not None:
        if not args.fda_mapping.is_file():
            raise FileNotFoundError(
                f"canonical FDA mapping does not exist: {args.fda_mapping}"
            )
        model_table = build_canonical_pk_identity_table(
            model_table,
            pd.read_csv(args.fda_mapping, low_memory=False),
        )
    source_rights = None
    if args.source_rights_manifest is not None:
        if not args.source_rights_manifest.is_file():
            raise FileNotFoundError(
                f"source-rights manifest does not exist: {args.source_rights_manifest}"
            )
        source_rights = pd.read_csv(args.source_rights_manifest, dtype="object")
    manifest = recover_pkdb_context(
        model_table=model_table,
        out_dir=args.out_dir,
        timeout=args.timeout,
        workers=args.workers,
        refresh=args.refresh,
        max_studies=args.max_studies,
        source_rights=source_rights,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
