from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.compare_run_consensus import compare_consensus_to_run_decoys
from analysis.ml.consensus_z_repair import repair_consensus_z_scores
from src.config.output_paths import run_output_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Write decoy-standardized score comparisons while preserving "
            "the canonical Atlas score columns by default."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    reference = parser.add_mutually_exclusive_group(required=True)
    reference.add_argument(
        "--reference-run-root",
        type=Path,
        help=(
            "Reference docked run root. Flat and nested variant/pH layouts are "
            "accepted; multiple strata for one PDB fail as ambiguous."
        ),
    )
    reference.add_argument(
        "--compare-runid",
        help=(
            "Compare a selected score to same-PDB decoys from this Atlas run ID "
            "without replacing canonical score fields."
        ),
    )
    parser.add_argument(
        "--compare-stream",
        choices=("docked_consensus", "post_docked_scorch"),
        default="docked_consensus",
        help=(
            "Comparison artifact stream for --compare-runid "
            "(default: docked_consensus)."
        ),
    )
    parser.add_argument(
        "--compare-score-column",
        default="consensus_score",
        help=(
            "Dataset and reference score column to standardize for --compare-runid "
            "(default: consensus_score; use scorch_composite for the post stream)."
        ),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--allow-reference-fallback",
        action="store_true",
        help=(
            "Allow rows absent from the reference run's FDA table to use the "
            "same-PDB frozen DUD null; provenance marks these rows explicitly."
        ),
    )
    parser.add_argument(
        "--replace-model-score-columns",
        action="store_true",
        help=(
            "Explicit legacy opt-in: replace consensus_score/atlas aliases with "
            "the repaired z-score. Do not use for --compare-runid."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.compare_runid:
        if args.replace_model_score_columns:
            raise ValueError(
                "--replace-model-score-columns cannot be used with --compare-runid"
            )
        runtime_name = (
            "post_docked" if args.compare_stream == "post_docked_scorch" else "docked"
        )
        comparison_root = run_output_dir(
            args.repo_root.resolve(), runtime_name, args.compare_runid
        )
        manifest = compare_consensus_to_run_decoys(
            args.dataset,
            comparison_root,
            args.out,
            comparison_run_id=args.compare_runid,
            compare_stream=args.compare_stream,
            compare_score_column=args.compare_score_column,
        )
    else:
        manifest = repair_consensus_z_scores(
            args.dataset,
            args.reference_run_root,
            args.out,
            allow_reference_fallback=args.allow_reference_fallback,
            replace_model_score_columns=args.replace_model_score_columns,
        )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
