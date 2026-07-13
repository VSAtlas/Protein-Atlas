from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.post_docked_final_score_backfill import (
    backfill_post_docked_final_scores,
)
from src.config.output_paths import run_output_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Backfill a selected post-docked numeric score and its provenance "
            "into an ML CSV using a strict (pdb_id, ligand_base) join."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    reference = parser.add_mutually_exclusive_group(required=True)
    reference.add_argument(
        "--reference-root",
        "--post-docked-reference-root",
        dest="reference_root",
        type=Path,
        help="Exact post-docked run root containing per-PDB score artifacts.",
    )
    reference.add_argument(
        "--reference-run-id",
        help="Atlas run ID whose canonical post_docked output should be used.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root used only to resolve --reference-run-id.",
    )
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument(
        "--reference-score-column",
        default="final_score",
        help="Numeric score column read from each reference artifact.",
    )
    parser.add_argument(
        "--output-score-column",
        help="Destination dataset column (default: same as --reference-score-column).",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help=(
            "Replace conflicting existing score/provenance cells. The default "
            "fills missing cells only and preserves all nonmissing input cells."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.reference_run_id:
        reference_root = run_output_dir(
            args.repo_root.expanduser().resolve(),
            "post_docked",
            args.reference_run_id,
        )
    else:
        reference_root = args.reference_root
    manifest = backfill_post_docked_final_scores(
        args.dataset,
        reference_root,
        args.out,
        args.audit_dir,
        overwrite_existing=args.overwrite_existing,
        reference_score_column=args.reference_score_column,
        output_score_column=args.output_score_column,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=True, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
