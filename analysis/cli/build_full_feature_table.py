from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.ml.feature_metadata import enrich_ml_feature_metadata


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Backfill Atlas ML feature columns into an existing model/PU/panel table."
    )
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument("--chemical-cluster", choices=["auto", "scaffold", "smiles", "chemotype", "none"], default="auto")
    parser.add_argument("--target-family", choices=["auto", "protein_class", "gene_heuristic", "none"], default="auto")
    parser.add_argument("--source-lineage", choices=["auto", "none"], default="auto")
    parser.add_argument("--drop-column", nargs="*", default=None)
    parser.add_argument("--feature-table", action="append", type=Path, default=[], help="Additional local feature/score table to join by drug/target/PDB keys; can be repeated.")
    args = parser.parse_args(argv)

    source = pd.read_csv(args.table, low_memory=False)
    enriched, summary = enrich_ml_feature_metadata(
        source,
        chemical_cluster=args.chemical_cluster,
        target_family=args.target_family,
        source_lineage=args.source_lineage,
        drop_columns=args.drop_column,
        run_dir=args.run_dir,
        repo_root=args.repo_root.resolve(),
        extra_feature_tables=args.feature_table,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(args.out, index=False)
    summary.update({"input": str(args.table), "output": str(args.out), "columns": int(len(enriched.columns))})
    args.out.with_suffix(".manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
