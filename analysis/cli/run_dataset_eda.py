from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.dataset_eda import DEFAULT_TOOLS, find_dataset_candidates, run_dataset_eda


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run read-only exploratory analysis on an Atlas ML/tabular dataset."
    )
    parser.add_argument("--dataset", type=Path, help="Dataset path: CSV/TSV/Parquet/JSONL/JSON.")
    parser.add_argument(
        "--dataset-name",
        help="Case-insensitive path substring to resolve under data/, outputs/data/, analysis/, or outputs/analysis/.",
    )
    parser.add_argument(
        "--list-matches",
        action="store_true",
        help="List dataset-name matches and exit without reading any dataset.",
    )
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("outputs/data/dataset_eda"),
        help="Output directory for EDA artifacts. The input dataset is never modified.",
    )
    parser.add_argument("--label-col", help="Discrete target column for mutual information ranking.")
    parser.add_argument(
        "--tool",
        action="append",
        default=[],
        help=f"Tool to run; repeatable or comma-separated. Defaults: {','.join(DEFAULT_TOOLS)}.",
    )
    parser.add_argument(
        "--tools",
        help=f"Comma-separated tool list. Defaults: {','.join(DEFAULT_TOOLS)}.",
    )
    parser.add_argument("--read-rows", type=int, help="Optional cap while loading rows from the dataset.")
    parser.add_argument(
        "--sample-rows",
        type=int,
        default=10_000,
        help="Rows sampled for expensive profile/association tools; <=0 disables sampling.",
    )
    parser.add_argument("--random-state", type=int, default=13)
    parser.add_argument("--max-association-columns", type=int, default=80)
    parser.add_argument("--max-category-levels", type=int, default=200)
    parser.add_argument("--max-mi-features", type=int, default=300)
    parser.add_argument("--network-threshold", type=float, default=0.35)
    parser.add_argument("--network-max-edges", type=int, default=300)
    parser.add_argument(
        "--full-ydata",
        action="store_true",
        help="Run the fuller ydata-profiling report instead of minimal mode.",
    )
    parser.add_argument("--dython-nominal-assoc", default="cramer")
    parser.add_argument(
        "--include-id-like",
        action="store_true",
        help="Allow ID/name-like columns in mutual-information feature ranking.",
    )
    parser.add_argument(
        "--include-high-cardinality",
        action="store_true",
        help="Allow high-cardinality categorical columns in association/MI inputs.",
    )
    parser.add_argument(
        "--fail-on-missing",
        action="store_true",
        help="Fail instead of recording skipped optional tools when packages are missing.",
    )
    return parser


def _selected_tools(args: argparse.Namespace) -> list[str] | None:
    values: list[str] = []
    if args.tools:
        values.append(str(args.tools))
    values.extend(str(value) for value in args.tool or [])
    return values or None


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list_matches:
        query = args.dataset_name or (str(args.dataset) if args.dataset else "")
        if not query:
            parser.error("--list-matches requires --dataset-name or a bare --dataset query")
        matches = find_dataset_candidates(query, repo_root=args.repo_root)
        for path in matches:
            print(path)
        return 0

    manifest = run_dataset_eda(
        args.dataset,
        args.out_dir,
        repo_root=args.repo_root,
        dataset_name=args.dataset_name,
        tools=_selected_tools(args),
        label_col=args.label_col,
        read_rows=args.read_rows,
        sample_rows=args.sample_rows,
        random_state=args.random_state,
        max_association_columns=args.max_association_columns,
        max_category_levels=args.max_category_levels,
        max_mi_features=args.max_mi_features,
        network_threshold=args.network_threshold,
        network_max_edges=args.network_max_edges,
        ydata_minimal=not args.full_ydata,
        dython_nominal_assoc=args.dython_nominal_assoc,
        include_id_like=args.include_id_like,
        include_high_cardinality=args.include_high_cardinality,
        fail_on_missing=args.fail_on_missing,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
