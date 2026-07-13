from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.external.openfda_pk_review import audit_openfda_pk_cache


DEFAULT_MODEL_TABLE = Path(
    "data/AtlasSPD_phase1/combined_activity_source_matched_20260709/"
    "model_ready/spd_binding_deduplicated.csv"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit cached SPL Cmax/protein-binding source text conservatively."
    )
    parser.add_argument("--model-table", type=Path, default=DEFAULT_MODEL_TABLE)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("data/external/dailymed_spl/phase1_openfda/records"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(
            "data/AtlasSPD_phase1/pk_context_v0_0_03/openfda_source_text_review"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.model_table.exists():
        raise FileNotFoundError(f"model table not found: {args.model_table}")
    model_table = pd.read_csv(args.model_table, low_memory=False)
    manifest = audit_openfda_pk_cache(
        cache_dir=args.cache_dir,
        model_table=model_table,
        out_dir=args.out_dir,
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
