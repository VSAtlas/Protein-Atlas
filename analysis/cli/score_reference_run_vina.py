from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping

from analysis.ml.addon_run_safety import (
    DEFAULT_ACTIVITY_GRACE_SECONDS,
    exclusive_addon_lock,
    require_idle_atlas,
)
from analysis.ml.reference_vina_compare import (
    plan_reference_addon_scoring,
    score_reference_vina_candidates_with_scorch,
    score_selected_against_reference_vina,
)


def _worker_count(value: str) -> int:
    workers = int(value)
    if not 1 <= workers <= 32:
        raise argparse.ArgumentTypeError("workers must be between 1 and 32")
    return workers


def _persist_final_manifest(
    manifest: dict[str, Any],
    *,
    out_dir: Path,
    lock_path: Path,
    dry_run: bool,
) -> dict[str, Any]:
    """Add execution provenance, then persist the exact returned manifest."""

    final_path = out_dir.resolve() / "reference_addon_manifest.json"
    raw_outputs = manifest.get("outputs")
    outputs = dict(raw_outputs) if isinstance(raw_outputs, Mapping) else {}
    outputs["final_manifest"] = str(final_path)
    if manifest.get("schema") == "atlas.reference-addon-score.v1":
        outputs["manifest"] = str(final_path)
    manifest["outputs"] = outputs
    manifest["exclusive_lock"] = str(lock_path)
    manifest["bypassed_stale_manifests"] = []
    manifest["fast_mode"] = False
    manifest["dud_only_mode"] = False
    manifest["dry_run"] = bool(dry_run)
    manifest["active_run_policy"] = {
        "check_result": "passed",
        "manifest_bypass_supported": False,
        "recent_manifest_grace_seconds": DEFAULT_ACTIVITY_GRACE_SECONDS,
    }
    final_path.parent.mkdir(parents=True, exist_ok=True)
    final_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Dock selected pairs with a frozen reference run's stage-1 Vina "
            "context and standardize raw energies against its explicit DUD null."
        )
    )
    parser.add_argument(
        "--run-id",
        required=True,
        help="Unique run ID used for provenance and the isolated SCORCH workspace.",
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
        "--dry-run",
        action="store_true",
        help="Validate concurrency policy and print the plan without docking.",
    )
    parser.add_argument(
        "--scorch",
        action="store_true",
        help=(
            "SCORCH-rescore every selected candidate pose and standardize the "
            "SCORCH composite against the same-PDB DUD null from --compare-run."
        ),
    )
    parser.add_argument(
        "--require-full-scorch-null",
        action="store_true",
        help=(
            "Require near-complete same-PDB DUD SCORCH coverage before any "
            "candidate docking. By default, a selection-truncated SCORCH null "
            "continues with a prominent sensitivity-only warning."
        ),
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
            "run_id": str(args.run_id),
            "selected_pairs": str(selected_pairs),
            "comparison_run_id": str(args.compare_runid),
            "out_dir": str(args.out_dir.resolve()),
            "scorch": bool(args.scorch),
        }
        with exclusive_addon_lock(repo_root, metadata=lock_metadata) as lock_path:
            plan = plan_reference_addon_scoring(
                selected_pairs_path=selected_pairs,
                comparison_run_id=args.compare_runid,
                out_dir=args.out_dir,
                run_id=args.run_id,
                repo_root=repo_root,
                vina_exe=args.vina_exe,
                workers=args.workers,
                scorch=bool(args.scorch),
            )
            for warning in plan.get("warnings", []):
                print(f"WARNING: {warning}", file=sys.stderr)
            if (
                args.scorch
                and args.require_full_scorch_null
                and not bool(plan.get("scorch_reference_full_coverage"))
            ):
                raise ValueError(
                    "SCORCH execution blocked by --require-full-scorch-null: "
                    "comparison-run DUD SCORCH coverage is selection-truncated"
                )
            if args.dry_run:
                manifest = plan
            else:
                stage1_manifest = score_selected_against_reference_vina(
                    selected_pairs_path=selected_pairs,
                    comparison_run_id=args.compare_runid,
                    out_dir=args.out_dir,
                    repo_root=repo_root,
                    vina_exe=args.vina_exe,
                    workers=args.workers,
                    run_id=args.run_id,
                )
                if args.scorch:
                    scorch_manifest = score_reference_vina_candidates_with_scorch(
                        stage1_manifest=stage1_manifest,
                        selected_pairs_path=selected_pairs,
                        comparison_run_id=args.compare_runid,
                        out_dir=args.out_dir,
                        run_id=args.run_id,
                        repo_root=repo_root,
                        workers=args.workers,
                        require_full_scorch_null=bool(
                            args.require_full_scorch_null
                        ),
                    )
                    manifest = {
                        "schema": "atlas.reference-addon-score.v1",
                        "status": "complete",
                        "run_id": str(args.run_id),
                        "comparison_run_id": str(args.compare_runid),
                        "stage1": stage1_manifest,
                        "scorch": scorch_manifest,
                    }
                else:
                    manifest = stage1_manifest
            manifest = _persist_final_manifest(
                manifest,
                out_dir=args.out_dir,
                lock_path=lock_path,
                dry_run=bool(args.dry_run),
            )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
