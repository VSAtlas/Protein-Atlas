from __future__ import annotations

import argparse
import os
import shlex
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence


from cli.qol.constants import (
    SLURM_BENCH2_CANARY_ARRAY,
    SLURM_BENCH2_CANARY_CPUS_PER_TASK,
    SLURM_BENCH2_CANARY_JOB_NAME,
    SLURM_BENCH2_CANARY_MAIN_ARGS,
    SLURM_BENCH2_CANARY_SCRIPT,
    SLURM_BENCH2_CANARY_TIME,
    SLURM_SUBMIT_DEFAULT_ARRAY,
    SLURM_SUBMIT_DEFAULT_CPUS_PER_TASK,
    SLURM_SUBMIT_DEFAULT_JOB_NAME,
    SLURM_SUBMIT_DEFAULT_SCRIPT,
)
from cli.qol.status import _cmd_status
from cli.qol import _bindings

def _cmd_slurm(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas slurm",
        description="Slurm submission, progress, and finalization wrappers.",
    )
    sub = parser.add_subparsers(dest="slurm_cmd", required=True)

    submit = sub.add_parser("submit", help="Submit the maintained Slurm worker/finalizer scripts.")
    submit.add_argument("--run-id", default="")
    submit.add_argument("--array", default=SLURM_SUBMIT_DEFAULT_ARRAY)
    submit.add_argument(
        "--cpus-per-task",
        default=str(SLURM_SUBMIT_DEFAULT_CPUS_PER_TASK),
        help="CPUs requested per Slurm task, or 'auto' to omit --cpus-per-task and use node CPU detection.",
    )
    submit.add_argument("--script", default=SLURM_SUBMIT_DEFAULT_SCRIPT)
    submit.add_argument("--finalize-script", default="tools/slurm/run_spr_finalize.sh")
    submit.add_argument("--partition")
    submit.add_argument("--account")
    submit.add_argument("--time")
    submit.add_argument("--job-name", default=SLURM_SUBMIT_DEFAULT_JOB_NAME)
    submit.add_argument("--mail-user", default="")
    submit.add_argument("--mail-type", default="BEGIN,END,FAIL")
    submit.add_argument(
        "--bench2-canary",
        action="store_true",
        help=(
            "Preset a low-SU 3-node/15-minute bench2-fast utilization canary "
            "(array 0-2%%3, exclusive whole-node workers, -bench2 -fast). Explicit --array, "
            "--cpus-per-task, --time, --job-name, and --main-args overrides are honored."
        ),
    )
    submit.add_argument(
        "--exclusive-node",
        action="store_true",
        help="Request exclusive nodes. Use --cpus-per-task auto to make Atlas detect the full node CPU count.",
    )
    submit.add_argument(
        "--execution-mode",
        default="",
        help="Optional ATLAS_EXECUTION_MODE export for worker scheduling policy.",
    )
    submit.add_argument(
        "--main-args",
        default="",
        help="Worker main args. Use {run_id} as a placeholder for the resolved run id.",
    )
    submit.add_argument(
        "--all-dirs",
        default="",
        help="Value for tools/run_relocated_mode.py --all-dirs; defaults to repo root.",
    )
    submit.add_argument(
        "--prepped-root",
        default="",
        help=(
            "Prepared ligand root exported as ATLAS_PREPPED_ROOT. When omitted for "
            "relocated runs, Atlas tries sensible run/repo/scratch defaults."
        ),
    )
    submit.add_argument("--dry-run", action="store_true")
    submit.add_argument("--no-finalize", action="store_true")

    finalize = sub.add_parser("finalize", help="Run distributed manifest reconciliation/final hooks.")
    finalize.add_argument("run_id")
    finalize.add_argument("--reconcile-only", action="store_true")
    finalize.add_argument("--no-hooks", action="store_true")

    repair = sub.add_parser(
        "repair-manifest",
        help="Repair/rebuild run_manifest.yaml from distributed chunk state.",
    )
    repair.add_argument("run_id")
    repair.add_argument(
        "--all-dirs",
        "--root",
        dest="all_dirs",
        default="",
        help="Relocated scratch root containing outputs/manifests, outputs/docked, and outputs/post_docked.",
    )
    repair.add_argument("--dry-run", action="store_true")
    repair.add_argument("--repair-scorch-selection", action="store_true")
    repair.add_argument("--preserve-incomplete-scorch-outputs", action="store_true")
    repair_scorch = repair.add_mutually_exclusive_group()
    repair_scorch.add_argument("--require-scorch", action="store_true")
    repair_scorch.add_argument("--no-require-scorch", action="store_true")

    progress = sub.add_parser("progress", help="Show run progress with live Slurm probes.")
    progress.add_argument("run_id")
    progress.add_argument("--watch", type=float, default=0.0)
    progress.add_argument("--errors", action="store_true")
    progress.add_argument("--deep", action="store_true")
    progress.add_argument("--html", nargs="?", const="")
    progress.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args(list(argv))

    if args.slurm_cmd == "submit":
        return _cmd_slurm_submit(args)
    if args.slurm_cmd == "finalize":
        return _cmd_slurm_finalize(args)
    if args.slurm_cmd == "repair-manifest":
        from cli.qol.manifest_repair import _cmd_repair_manifest

        repair_args = [args.run_id]
        if args.all_dirs:
            repair_args.extend(["--all-dirs", str(args.all_dirs)])
        if args.dry_run:
            repair_args.append("--dry-run")
        if getattr(args, "repair_scorch_selection", False):
            repair_args.append("--repair-scorch-selection")
        if getattr(args, "preserve_incomplete_scorch_outputs", False):
            repair_args.append("--preserve-incomplete-scorch-outputs")
        if getattr(args, "require_scorch", False):
            repair_args.append("--require-scorch")
        if getattr(args, "no_require_scorch", False):
            repair_args.append("--no-require-scorch")
        return _cmd_repair_manifest(repair_args)
    if args.slurm_cmd == "progress":
        status_args = [args.run_id, "--live"]
        if args.watch:
            status_args.extend(["--watch", str(args.watch)])
        if args.errors:
            status_args.append("--errors")
        if args.deep:
            status_args.append("--deep")
        if args.as_json:
            status_args.append("--json")
        if args.html is not None:
            status_args.append("--html")
            if args.html:
                status_args.append(args.html)
        return _cmd_status(status_args)
    return 1


def _main_args_are_benchmark(main_args: str) -> bool:
    try:
        tokens = set(shlex.split(str(main_args or "")))
    except ValueError:
        tokens = set(str(main_args or "").split())
    return bool(
        tokens
        & {
            "-bench",
            "--bench",
            "-bench2",
            "--bench2",
            "-bench-small",
            "--bench-small",
            "-bench-micro",
            "--bench-micro",
            "-water-bench",
            "--water-bench",
            "--waterbench",
        }
    )


def _extract_sbatch_job_id(stdout: str) -> str:
    for token in reversed(str(stdout or "").split()):
        candidate = token.split(";", 1)[0].strip()
        if candidate.isdigit():
            return candidate
    return str(stdout or "").strip().split(";", 1)[0]


def _path_has_prepped_payload(path: Path) -> bool:
    if not path.is_dir():
        return False
    try:
        for current, _dirnames, filenames in os.walk(path):
            del current
            for name in filenames:
                if name == "_manifest.json" or name.endswith(".pdbqt"):
                    return True
    except OSError:
        return False
    return False


def _scratch_mirror_prepped_root(root: Path) -> Path | None:
    parts = root.resolve().parts
    if len(parts) < 5:
        return None
    first = parts[1] if len(parts) > 1 else ""
    if not first.startswith("work"):
        return None
    return Path("/scratch", *parts[2:], "prepped_ligands")


def _resolve_submit_prepped_root(
    *,
    explicit: str,
    root: Path,
    all_dirs: str,
) -> str:
    token = str(explicit or "").strip()
    if token:
        return str(Path(token).expanduser().resolve())
    env_token = str(os.environ.get("ATLAS_PREPPED_ROOT") or "").strip()
    if env_token:
        return str(Path(env_token).expanduser().resolve())
    candidates: list[Path] = []
    if str(all_dirs or "").strip():
        candidates.append(Path(all_dirs).expanduser().resolve() / "prepped_ligands")
    candidates.append(root / "prepped_ligands")
    scratch_env = str(os.environ.get("SCRATCH") or "").strip()
    if scratch_env:
        candidates.append(Path(scratch_env).expanduser().resolve() / "atlas" / "code" / "protein_automation" / "prepped_ligands")
    mirror = _scratch_mirror_prepped_root(root)
    if mirror is not None:
        candidates.append(mirror)
    for candidate in candidates:
        if _path_has_prepped_payload(candidate):
            return str(candidate)
    return ""


def _cmd_slurm_submit(args: argparse.Namespace) -> int:
    _apply_slurm_submit_presets(args)
    try:
        cpus_per_task = _slurm_cpus_per_task_value(getattr(args, "cpus_per_task", ""))
    except ValueError as exc:
        print(f"atlas slurm submit: {exc}", file=sys.stderr)
        return 2
    if cpus_per_task is not None and cpus_per_task <= 0:
        print("atlas slurm submit: --cpus-per-task must be > 0", file=sys.stderr)
        return 2

    root = _bindings.repo_root()
    run_id = str(args.run_id or f"atlas_slurm_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}")
    worker_script = Path(args.script)
    if not worker_script.is_absolute():
        worker_script = root / worker_script
    finalizer_arg = str(args.finalize_script or "tools/slurm/run_spr_finalize.sh")
    if (
        finalizer_arg == "tools/slurm/run_spr_finalize.sh"
        and (
            worker_script.name.startswith("run_relocated_")
            or bool(str(getattr(args, "all_dirs", "") or "").strip())
        )
    ):
        finalizer_arg = "tools/slurm/run_relocated_finalize.sh"
    finalizer_script = Path(finalizer_arg)
    if not finalizer_script.is_absolute():
        finalizer_script = root / finalizer_script
    if not worker_script.exists():
        print(f"atlas slurm submit: worker script not found: {worker_script}", file=sys.stderr)
        return 1
    if not args.no_finalize and not finalizer_script.exists():
        print(f"atlas slurm submit: finalizer script not found: {finalizer_script}", file=sys.stderr)
        return 1

    main_args = str(args.main_args or "--resume --run-id {run_id}").replace("{run_id}", run_id)
    all_dirs = str(args.all_dirs or root)
    prepped_root = _resolve_submit_prepped_root(
        explicit=str(getattr(args, "prepped_root", "") or ""),
        root=root,
        all_dirs=all_dirs,
    )
    export_items = [
        "ALL",
        f"RUN_ID={run_id}",
        f"ATLAS_MAIN_ARGS={main_args}",
        f"ATLAS_ALL_DIRS={all_dirs}",
        f"ATLAS_REPO_ROOT={root}",
    ]
    if prepped_root:
        export_items.append(f"ATLAS_PREPPED_ROOT={prepped_root}")
    use_node_cpus = cpus_per_task is None
    if use_node_cpus:
        export_items.append("ATLAS_USE_NODE_CPUS=1")
    bench_seed = str(os.environ.get("ATLAS_THROUGHPUT_BENCH_SEED", "")).strip()
    if bench_seed:
        if os.environ.get("ATLAS_THROUGHPUT_BENCH") == "1" and _main_args_are_benchmark(main_args):
            export_items.append("ATLAS_THROUGHPUT_BENCH=1")
            export_items.append(f"ATLAS_THROUGHPUT_BENCH_SEED={bench_seed}")
            prep_required = str(os.environ.get("ATLAS_BENCH_CACHED_PREP_REQUIRED", "")).strip()
            if prep_required:
                export_items.append(f"ATLAS_BENCH_CACHED_PREP_REQUIRED={prep_required}")
        else:
            bench_seed = ""
    if str(args.execution_mode or "").strip():
        export_items.append(f"ATLAS_EXECUTION_MODE={str(args.execution_mode).strip()}")
    worker_cmd = [
        "sbatch",
        "--parsable",
        f"--chdir={root}",
        f"--array={args.array}",
        f"--job-name={args.job_name}",
        "--export=" + ",".join(export_items),
    ]
    if cpus_per_task is not None:
        worker_cmd.append(f"--cpus-per-task={int(cpus_per_task)}")
    if bool(getattr(args, "exclusive_node", False)):
        worker_cmd.append("--exclusive")
    if args.partition:
        worker_cmd.append(f"--partition={args.partition}")
    if args.account:
        worker_cmd.append(f"--account={args.account}")
    if args.time:
        worker_cmd.append(f"--time={args.time}")
    if str(getattr(args, "mail_user", "") or "").strip():
        worker_cmd.append(f"--mail-user={str(args.mail_user).strip()}")
        worker_cmd.append(f"--mail-type={str(args.mail_type or 'BEGIN,END,FAIL').strip()}")
    worker_cmd.append(str(worker_script))

    if args.dry_run:
        print("RUN_ID=" + run_id)
        print("ATLAS_MAIN_ARGS=" + main_args)
        print("ATLAS_ALL_DIRS=" + all_dirs)
        print("ATLAS_PREPPED_ROOT=" + (prepped_root or "auto"))
        print("SLURM_CPUS_PER_TASK=" + (str(cpus_per_task) if cpus_per_task is not None else "auto"))
        if use_node_cpus:
            print("ATLAS_USE_NODE_CPUS=1")
        if bool(getattr(args, "exclusive_node", False)):
            print("SLURM_EXCLUSIVE_NODE=1")
        if bench_seed:
            print("ATLAS_THROUGHPUT_BENCH=1")
            print("ATLAS_THROUGHPUT_BENCH_SEED=" + bench_seed)
            prep_required = str(os.environ.get("ATLAS_BENCH_CACHED_PREP_REQUIRED", "")).strip()
            if prep_required:
                print("ATLAS_BENCH_CACHED_PREP_REQUIRED=" + prep_required)
        if str(args.execution_mode or "").strip():
            print("ATLAS_EXECUTION_MODE=" + str(args.execution_mode).strip())
        if str(getattr(args, "mail_user", "") or "").strip():
            print("SLURM_MAIL_USER=" + str(args.mail_user).strip())
            print("SLURM_MAIL_TYPE=" + str(args.mail_type or "BEGIN,END,FAIL").strip())
        print("worker: " + " ".join(shlex.quote(part) for part in worker_cmd))
        if not args.no_finalize:
            print("finalizer_dependency=afterany:<ARRAY_JOB_ID>")
            finalizer_preview = [
                "sbatch",
                "--parsable",
                f"--chdir={root}",
                "--dependency=afterany:<ARRAY_JOB_ID>",
                f"--job-name={args.job_name}_finalize",
                "--export="
                + ",".join(
                    [
                        "ALL",
                        f"RUN_ID={run_id}",
                        f"ATLAS_ALL_DIRS={all_dirs}",
                        f"ATLAS_REPO_ROOT={root}",
                        *([f"ATLAS_PREPPED_ROOT={prepped_root}"] if prepped_root else []),
                    ]
                ),
            ]
            if args.partition:
                finalizer_preview.append(f"--partition={args.partition}")
            if args.account:
                finalizer_preview.append(f"--account={args.account}")
            if args.time:
                finalizer_preview.append(f"--time={args.time}")
            if str(getattr(args, "mail_user", "") or "").strip():
                finalizer_preview.append(f"--mail-user={str(args.mail_user).strip()}")
                finalizer_preview.append(f"--mail-type={str(args.mail_type or 'BEGIN,END,FAIL').strip()}")
            finalizer_preview.append(str(finalizer_script))
            print("finalizer: " + " ".join(shlex.quote(part) for part in finalizer_preview))
        return 0

    (root / "logs" / "slurm").mkdir(parents=True, exist_ok=True)
    worker_result = _bindings.subprocess().run(worker_cmd, cwd=root, check=False, text=True, capture_output=True)
    if worker_result.returncode != 0:
        sys.stderr.write(worker_result.stderr)
        return int(worker_result.returncode)
    array_job_id = _extract_sbatch_job_id(worker_result.stdout)
    print(f"Submitted array job: {array_job_id}")
    print(f"RUN_ID={run_id}")
    if args.no_finalize:
        return 0

    finalizer_cmd = [
        "sbatch",
        "--parsable",
        f"--chdir={root}",
        f"--dependency=afterany:{array_job_id}",
        f"--job-name={args.job_name}_finalize",
        "--export="
        + ",".join(
            [
                "ALL",
                f"RUN_ID={run_id}",
                f"ATLAS_ALL_DIRS={all_dirs}",
                f"ATLAS_REPO_ROOT={root}",
                *([f"ATLAS_PREPPED_ROOT={prepped_root}"] if prepped_root else []),
            ]
        ),
    ]
    if args.partition:
        finalizer_cmd.append(f"--partition={args.partition}")
    if args.account:
        finalizer_cmd.append(f"--account={args.account}")
    if args.time:
        finalizer_cmd.append(f"--time={args.time}")
    if str(getattr(args, "mail_user", "") or "").strip():
        finalizer_cmd.append(f"--mail-user={str(args.mail_user).strip()}")
        finalizer_cmd.append(f"--mail-type={str(args.mail_type or 'BEGIN,END,FAIL').strip()}")
    finalizer_cmd.append(str(finalizer_script))
    finalizer_result = _bindings.subprocess().run(finalizer_cmd, cwd=root, check=False, text=True, capture_output=True)
    if finalizer_result.returncode != 0:
        sys.stderr.write(finalizer_result.stderr)
        return int(finalizer_result.returncode)
    final_job_id = _extract_sbatch_job_id(finalizer_result.stdout)
    print(f"Submitted finalize job: {final_job_id} (afterany:{array_job_id})")
    return 0


def _apply_slurm_submit_presets(args: argparse.Namespace) -> None:
    if not bool(getattr(args, "bench2_canary", False)):
        return

    if str(getattr(args, "array", "") or "") == SLURM_SUBMIT_DEFAULT_ARRAY:
        args.array = SLURM_BENCH2_CANARY_ARRAY
    args.exclusive_node = True
    if str(getattr(args, "cpus_per_task", "") or "") == str(SLURM_SUBMIT_DEFAULT_CPUS_PER_TASK):
        args.cpus_per_task = SLURM_BENCH2_CANARY_CPUS_PER_TASK
    if str(getattr(args, "script", "") or "") == SLURM_SUBMIT_DEFAULT_SCRIPT:
        args.script = SLURM_BENCH2_CANARY_SCRIPT
    if not str(getattr(args, "time", "") or "").strip():
        args.time = SLURM_BENCH2_CANARY_TIME
    if str(getattr(args, "job_name", "") or "") == SLURM_SUBMIT_DEFAULT_JOB_NAME:
        args.job_name = SLURM_BENCH2_CANARY_JOB_NAME
    if not str(getattr(args, "main_args", "") or "").strip():
        args.main_args = SLURM_BENCH2_CANARY_MAIN_ARGS


def _slurm_cpus_per_task_value(raw: Any) -> int | None:
    token = str(raw or "").strip().lower()
    if token in {"auto", "node", "whole-node", "whole_node", "all"}:
        return None
    try:
        return int(token)
    except Exception:
        raise ValueError(f"invalid --cpus-per-task value: {raw!r}") from None


def _cmd_slurm_finalize(args: argparse.Namespace) -> int:
    cmd = [
        sys.executable,
        str(_bindings.repo_root() / "tools" / "finalize_distributed_run.py"),
        "--run-id",
        str(args.run_id),
    ]
    if args.reconcile_only:
        cmd.append("--reconcile-only")
    if args.no_hooks:
        cmd.append("--no-hooks")
    return _bindings.subprocess().run(cmd, cwd=_bindings.repo_root(), check=False).returncode
