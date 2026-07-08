from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.feature_metadata import refresh_tables


DEFAULT_PATTERNS = [
    "*model_ready.csv",
    "moe_experts/*model_ready.csv",
    "four_state_evidence*/four_state_joined_source_balanced*.csv",
    "negative_evidence/mechanism_four_state_labels*.csv",
]


def _discover(run_dir: Path, patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in run_dir.glob(pattern):
            if not path.is_file() or path in seen:
                continue
            seen.add(path)
            paths.append(path)
    return sorted(paths)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Refresh Atlas ML feature metadata such as chemical clusters, target families, and source lineage."
    )
    parser.add_argument("--run-dir", type=Path, default=Path("data/pilotstudy"))
    parser.add_argument("--dataset", nargs="*", type=Path, default=None)
    parser.add_argument("--include-glob", nargs="*", default=DEFAULT_PATTERNS)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--in-place", action="store_true")
    parser.add_argument("--chemical-cluster", choices=["auto", "scaffold", "smiles", "chemotype", "none"], default="auto")
    parser.add_argument("--target-family", choices=["auto", "protein_class", "gene_heuristic", "none"], default="auto")
    parser.add_argument("--source-lineage", choices=["auto", "none"], default="auto")
    parser.add_argument("--drop-column", nargs="*", default=None)
    args = parser.parse_args(argv)

    run_dir = args.run_dir.resolve()
    paths = [path.resolve() for path in args.dataset] if args.dataset else _discover(run_dir, list(args.include_glob))
    out_dir = args.out_dir.resolve() if args.out_dir else (run_dir / "ml_feature_metadata")
    manifest = refresh_tables(
        paths,
        out_dir=out_dir,
        in_place=bool(args.in_place),
        run_dir=run_dir,
        chemical_cluster=args.chemical_cluster,
        target_family=args.target_family,
        source_lineage=args.source_lineage,
        drop_columns=args.drop_column,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
