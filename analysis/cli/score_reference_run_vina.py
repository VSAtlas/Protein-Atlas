from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from analysis.ml.addon_run_safety import (
    active_atlas_runs,
    active_runs_as_dicts,
    exclusive_addon_lock,
    require_idle_atlas,
)
from analysis.ml.reference_vina_compare import score_selected_against_reference_vina


def _worker_count(value: str) -> int:
    workers = int(value)
    if not 1 <= workers <= 32:
        raise argparse.ArgumentTypeError("workers must be between 1 and 32")
    return workers


def _minimum_exhaustiveness(value: str) -> int:
    exhaustiveness = int(value)
    if exhaustiveness < 2:
        raise argparse.ArgumentTypeError(
            "minimum exhaustiveness must be at least 2 for ML add-ons"
        )
    return exhaustiveness


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Dock selected pairs with a frozen reference run's stage-1 Vina "
            "context and standardize raw energies against its explicit DUD null."
        )
    )
    parser.add_argument(
        "--selected-pairs",
        "--pair-manifest",
        "--library-manifest",
        dest="selected_pairs",
        required=True,
        type=Path,
        help="Exact PDB-ligand pair table emitted by add-on staging.",
    )
    parser.add_argument(
        "--compare-runid",
        "--compare-run",
        dest="compare_runid",
        required=True,
        help="Completed full-accuracy run supplying receptor, grid, and DUD null context.",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--vina-exe", default="vina")
    parser.add_argument("--workers", type=_worker_count, default=8)
    parser.add_argument(
        "--minimum-exhaustiveness",
        type=_minimum_exhaustiveness,
        default=2,
        help="Reject comparison contexts below this Vina exhaustiveness (default: 2).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate concurrency policy and print the plan without docking.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        repo_root = args.repo_root.resolve()
        selected_pairs = args.selected_pairs.resolve()
        if not selected_pairs.is_file():
            raise FileNotFoundError(selected_pairs)
        require_idle_atlas(repo_root)
        lock_metadata = {
            "selected_pairs": str(selected_pairs),
            "comparison_run_id": str(args.compare_runid),
            "out_dir": str(args.out_dir.resolve()),
        }
        with exclusive_addon_lock(repo_root, metadata=lock_metadata) as lock_path:
            if args.dry_run:
                manifest = {
                    "status": "planned",
                    "selected_pairs": str(selected_pairs),
                    "comparison_run_id": str(args.compare_runid),
                    "out_dir": str(args.out_dir.resolve()),
                    "workers": int(args.workers),
                    "minimum_exhaustiveness": int(args.minimum_exhaustiveness),
                    "active_runs": active_runs_as_dicts(active_atlas_runs(repo_root)),
                    "exclusive_lock": str(lock_path),
                    "fast_mode": False,
                    "dud_only_mode": False,
                }
            else:
                manifest = score_selected_against_reference_vina(
                    selected_pairs_path=selected_pairs,
                    comparison_run_id=args.compare_runid,
                    out_dir=args.out_dir,
                    repo_root=repo_root,
                    vina_exe=args.vina_exe,
                    workers=args.workers,
                    minimum_reference_exhaustiveness=args.minimum_exhaustiveness,
                )
                manifest["exclusive_lock"] = str(lock_path)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
