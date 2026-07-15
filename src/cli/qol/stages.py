"""Canonical Atlas wrappers for modular single-stage repeat planning."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from cli.qol import _bindings
from cli.stage_repeat_plan import (
    MAX_STAGE_CPUS,
    SUPPORTED_STAGES,
    PlanRequest,
    StageRepeatPlanError,
    build_stage_repeat_plan,
)


def _cmd_stages(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas stages",
        description="Plan an exact, auditable repeat of one Atlas pipeline stage.",
    )
    sub = parser.add_subparsers(dest="stages_cmd", required=True)
    repeat = sub.add_parser(
        "repeat",
        help="Persist a deterministic dry-run pair plan; never launch workloads.",
    )
    repeat.add_argument("--run-id", required=True)
    repeat.add_argument("--stage", required=True, choices=SUPPORTED_STAGES)
    receptor_scope = repeat.add_mutually_exclusive_group(required=True)
    receptor_scope.add_argument("--receptor", action="append", default=[])
    receptor_scope.add_argument("--all-receptors", action="store_true")
    ligand_scope = repeat.add_mutually_exclusive_group(required=True)
    ligand_scope.add_argument("--ligand", action="append", default=[])
    ligand_scope.add_argument("--all-ligands", action="store_true")
    repeat.add_argument(
        "--variant",
        action="append",
        default=[],
        help="Exact receptor variant; repeat for more than one.",
    )
    repeat.add_argument(
        "--ph",
        action="append",
        default=[],
        help="Exact pH label; use 'base' for the empty/base context.",
    )
    repeat.add_argument(
        "--all-contexts",
        action="store_true",
        help="Explicitly accept every variant/pH context selected by the other filters.",
    )
    selection = repeat.add_mutually_exclusive_group()
    selection.add_argument(
        "--fraction",
        type=float,
        help="Deterministically select this fraction of exact pair contexts (default: 1).",
    )
    selection.add_argument(
        "--count",
        type=int,
        help="Deterministically select exactly this many pair contexts.",
    )
    repeat.add_argument("--seed", type=int, default=0)
    repeat.add_argument(
        "--selection-strategy",
        choices=("hash", "top-score"),
        default="",
        help=(
            "Required with --count or --fraction below 1. 'hash' is an unbiased, "
            "seeded subset; 'top-score' fails closed until score direction/provenance/ties "
            "have an approved contract."
        ),
    )
    repeat.add_argument(
        "--source-stage",
        default="",
        help=(
            "Exact docking stage. Required for Vina/GNINA; otherwise defaults to "
            "the causally selected result stage."
        ),
    )
    repeat.add_argument(
        "--cpus",
        type=int,
        default=1,
        help=f"Resource request recorded in the plan (hard maximum: {MAX_STAGE_CPUS}).",
    )
    source = repeat.add_mutually_exclusive_group()
    source.add_argument(
        "--database",
        type=Path,
        help="Atlas schema-v5+ SQLite snapshot with selected-attempt and typed-artifact lineage.",
    )
    source.add_argument(
        "--master-rows",
        type=Path,
        help="Explicit master_rows.csv; defaults to the run-scoped report export.",
    )
    repeat.add_argument("--output", type=Path)
    repeat.add_argument(
        "--dry-run",
        action="store_true",
        default=True,
        help="Plan only (always enabled in this release; no workloads are launched).",
    )
    repeat.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args(list(argv))

    if args.stages_cmd != "repeat":
        return 1
    fraction = args.fraction
    if fraction is None and args.count is None:
        fraction = 1.0
    request = PlanRequest(
        run_id=str(args.run_id),
        stage=str(args.stage),
        receptors=tuple(str(value).upper() for value in args.receptor),
        all_receptors=bool(args.all_receptors),
        ligands=tuple(str(value) for value in args.ligand),
        all_ligands=bool(args.all_ligands),
        variants=tuple(str(value) for value in args.variant),
        ph_labels=tuple(str(value) for value in args.ph),
        all_contexts=bool(args.all_contexts),
        fraction=fraction,
        count=args.count,
        seed=int(args.seed),
        cpus=int(args.cpus),
        source_stage=str(args.source_stage),
        selection_strategy=str(args.selection_strategy),
    )
    try:
        result = build_stage_repeat_plan(
            _bindings.repo_root(),
            request,
            database=args.database,
            master_rows=args.master_rows,
            output=args.output,
        )
    except StageRepeatPlanError as exc:
        print(f"atlas stages repeat: {exc}", file=sys.stderr)
        return 2

    plan = result["plan"]
    summary = plan["summary"]
    output_payload = {
        "plan_id": plan["plan_id"],
        "plan_path": str(result["plan_path"]),
        "plan_sha256": result["plan_sha256"],
        "pairs_path": str(result["pairs_path"]),
        "ledger_paths": [str(path) for path in result["ledger_paths"]],
        "selected": summary["selected"],
        "ready": summary["ready"],
        "blocked": summary["blocked"],
        "dry_run": True,
        "execution_status": "not_executed",
        "executor_status": plan["execution"]["executor_status"],
        "execution_gap": plan["execution"]["gap"],
        "full_matrix_supported": plan["scalability"]["full_matrix_supported"],
        "reused_existing": bool(result["reused_existing"]),
    }
    if args.as_json:
        print(json.dumps(output_payload, indent=2, sort_keys=True))
    else:
        print(
            "Atlas stage-repeat dry-run: "
            f"selected={summary['selected']} ready={summary['ready']} "
            f"blocked={summary['blocked']}"
        )
        print(f"Plan: {result['plan_path']}")
        print(f"Pair index: {result['pairs_path']}")
        print(f"Detailed ledger shards: {len(result['ledger_paths'])}")
        print(f"Executor: planner_only (gap={plan['execution']['gap']})")
        print("No docking, rescoring, validation, or MM/GBSA workload was launched.")
    return 0


__all__ = ["_cmd_stages"]
