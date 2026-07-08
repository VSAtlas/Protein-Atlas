from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence

from config.output_paths import run_output_dir

from cli.qol.reporting import _cmd_report, _run_master_export
from cli.qol.status import _cmd_status
from cli.qol import _bindings

def _cmd_debug(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas debug",
        description="Summarize structured run failures without reading raw logs.",
    )
    parser.add_argument("run_id")
    parser.add_argument("--deep", action="store_true")
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument("--html", nargs="?", const="")
    args = parser.parse_args(list(argv))

    status_args = [args.run_id, "--errors", "--explain"]
    if args.deep:
        status_args.append("--deep")
    if args.as_json:
        status_args.append("--json")
    if args.html is not None:
        status_args.append("--html")
        if args.html:
            status_args.append(args.html)
    return _cmd_status(status_args)


def _cmd_artifacts(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas artifacts",
        description="Artifact storage and measurement workflows.",
    )
    sub = parser.add_subparsers(dest="artifacts_cmd", required=True)
    measure = sub.add_parser("measure", help="Measure run-scoped output artifacts.")
    measure.add_argument("--run-id", required=True)
    measure.add_argument("--root", action="append", default=[])
    measure.add_argument("--capture", choices=["baseline", "final"])
    measure.add_argument("--diff", action="store_true")
    measure.add_argument("--out")
    measure.add_argument("--top", type=int, default=50)
    measure.add_argument("--dir-depth", type=int, default=3)
    measure.add_argument("--max-records", type=int)
    measure.add_argument("--verbose", action="store_true")
    args = parser.parse_args(list(argv))

    if args.artifacts_cmd != "measure":
        return 1

    root = _bindings.repo_root()
    data_dir = run_output_dir(root, "data", args.run_id)
    forwarded = ["--run-id", args.run_id, "--top", str(args.top), "--dir-depth", str(args.dir_depth)]
    if args.capture:
        forwarded.extend(["--capture", args.capture])
    if args.diff:
        forwarded.append("--diff")
    if args.max_records:
        forwarded.extend(["--max-records", str(args.max_records)])
    if args.verbose:
        forwarded.append("--verbose")
    out_dir = Path(args.out).expanduser() if args.out else data_dir / "artifact_reports"
    if not out_dir.is_absolute():
        out_dir = root / out_dir
    forwarded.extend(["--out", str(out_dir)])

    roots = [Path(item).expanduser() for item in args.root]
    if not roots:
        roots = [
            run_output_dir(root, "docked", args.run_id),
            run_output_dir(root, "post_docked", args.run_id),
            data_dir,
        ]
    existing_roots = [path if path.is_absolute() else root / path for path in roots]
    existing_roots = [path for path in existing_roots if path.exists()]
    if not existing_roots:
        print(f"atlas artifacts measure: no run output roots found for {args.run_id}", file=sys.stderr)
        return 1
    for path in existing_roots:
        forwarded.extend(["--root", str(path)])

    from tools import measure_artifacts

    return int(measure_artifacts.main(forwarded))


def _cmd_screenshot(argv: Sequence[str]) -> int:
    from cli.screenshot_cli import run_screenshot

    return run_screenshot(argv)


def _cmd_dev(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas dev",
        description="Developer maintenance workflows.",
    )
    sub = parser.add_subparsers(dest="dev_cmd", required=True)
    verify = sub.add_parser("verify", help="Run the repo quality gate.")
    verify.add_argument("--fix", action="store_true")
    verify.add_argument("--smoke", action="store_true")
    verify.add_argument("--full", action="store_true")
    args = parser.parse_args(list(argv))

    if args.dev_cmd == "verify":
        cmd = ["bash", str(_bindings.repo_root() / "tools" / "quality_gate.sh")]
        if args.fix:
            cmd.append("--fix")
        if args.smoke:
            cmd.append("--smoke")
        if args.full:
            cmd.append("--full")
        return _bindings.subprocess().run(cmd, cwd=_bindings.repo_root(), check=False).returncode
    return 1


def _cmd_reproduce(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas reproduce",
        description="Reproducibility helpers that reuse run-scoped output directories.",
    )
    sub = parser.add_subparsers(dest="reproduce_cmd", required=True)
    bundle = sub.add_parser("bundle", help="Write a run artifact index beside existing outputs.")
    bundle.add_argument("run_id")
    bundle.add_argument(
        "--out",
        help="Output directory for the bundle index; defaults to outputs/data/<RUN_ID>/.",
    )
    bundle.add_argument("--skip-report", action="store_true")
    bundle.add_argument("--skip-master-export", action="store_true")
    bundle.add_argument("--status-html", action="store_true", default=True)
    bundle.add_argument("--no-status-html", action="store_false", dest="status_html")
    args = parser.parse_args(list(argv))

    if args.reproduce_cmd != "bundle":
        return 1

    root = _bindings.repo_root()
    run_id = str(args.run_id)
    data_dir = run_output_dir(root, "data", run_id)
    out_dir = Path(args.out).expanduser() if args.out else data_dir
    if not out_dir.is_absolute():
        out_dir = root / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rc = 0
    if not args.skip_report:
        report_args = [run_id]
        if args.skip_master_export:
            report_args.append("--skip-master-export")
        rc = _cmd_report(report_args)
        if rc != 0:
            return rc
    elif not args.skip_master_export:
        rc = _run_master_export(run_id)
        if rc != 0:
            return rc

    status_html_path = data_dir / "status.html"
    if args.status_html:
        rc = _cmd_status([run_id, "--html", str(status_html_path)])
        if rc != 0:
            return rc

    from cli.status_dashboard import build_status_dashboard

    status = build_status_dashboard(root, run_id, include_slurm=True, include_errors=True)
    status_json_path = out_dir / "status_dashboard.json"
    status_json_path.write_text(json.dumps(status, indent=2, sort_keys=True, default=str), encoding="utf-8")

    candidates = {
        "manifest": run_output_dir(root, "manifests", run_id) / "run_manifest.yaml",
        "config_snapshot": run_output_dir(root, "configs", run_id),
        "master_rows_csv": data_dir / "master_rows.csv",
        "report_html": data_dir / "report.html",
        "report_yaml": data_dir / "report.yaml",
        "status_html": status_html_path,
        "status_json": status_json_path,
        "throughput_integrity_json": data_dir / "throughput_integrity.json",
        "throughput_integrity_csv": data_dir / "throughput_integrity.csv",
    }
    index = {
        "run_id": run_id,
        "created_utc": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "repo_root": str(root),
        "bundle_dir": str(out_dir),
        "artifacts": {
            key: str(path)
            for key, path in candidates.items()
            if path.exists()
        },
        "missing": [
            key
            for key, path in candidates.items()
            if not path.exists()
        ],
    }
    index_path = out_dir / "reproducibility_bundle.json"
    index_path.write_text(json.dumps(index, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Reproducibility bundle index: {index_path}")
    return 0


def _cmd_throughput(argv: Sequence[str]) -> int:
    from cli import throughput_bench

    try:
        return throughput_bench.run(argv, _bindings.repo_root())
    except ValueError as exc:
        print(f"atlas throughput: {exc}", file=sys.stderr)
        return 2
