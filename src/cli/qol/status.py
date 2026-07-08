from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

from config.output_paths import output_root, run_output_dir
from cli.qol import _bindings


def _cmd_status(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas status",
        description="Summarize an Atlas run and point to the useful outputs.",
    )
    parser.add_argument("run_id", nargs="?")
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument("--explain", action="store_true", help="Show plain-language operator actions.")
    parser.add_argument(
        "--html",
        nargs="?",
        const="",
        metavar="PATH",
        help="Write a standalone browser dashboard HTML file.",
    )
    parser.add_argument("--live", action="store_true", help="Include live scheduler/progress probes.")
    parser.add_argument("--errors", action="store_true", help="Include structured root-cause error summary.")
    parser.add_argument("--deep", action="store_true", help="Deep-scan run logs when summarizing errors.")
    parser.add_argument("--no-slurm", action="store_true", help="Skip Slurm squeue/sacct probes.")
    parser.add_argument(
        "--watch",
        type=float,
        default=0.0,
        metavar="SECONDS",
        help="Refresh status repeatedly at this interval.",
    )
    args = parser.parse_args(list(argv))

    root = _bindings.repo_root()
    run_id = args.run_id or _latest_run_id(root)
    if not run_id:
        print("No run id found. Start with `atlas demo` or `atlas --pdb <ID> --fast`.")
        return 1

    from cli.status_dashboard import (
        build_status_dashboard,
        default_status_html_path,
        format_status_dashboard,
        write_status_html,
    )

    include_errors = args.errors or args.deep
    include_slurm = not args.no_slurm

    def render_once() -> dict[str, Any]:
        status = build_status_dashboard(
            root,
            run_id,
            include_slurm=include_slurm,
            include_errors=include_errors,
            deep_errors=args.deep,
        )
        html_path: Path | None = None
        if args.html is not None:
            html_path = Path(args.html).expanduser() if args.html else default_status_html_path(root, run_id)
            if not html_path.is_absolute():
                html_path = root / html_path
            write_status_html(status, html_path)
            status["html_path"] = str(html_path)
        if args.as_json:
            print(json.dumps(status, indent=2, sort_keys=True, default=str))
        else:
            print(format_status_dashboard(status, explain=args.explain))
        if html_path is not None and not args.as_json:
            print(f"html: {html_path}")
        return status

    if args.watch and args.watch > 0:
        rc = 0
        try:
            while True:
                status = render_once()
                rc = 0 if status["found_any"] else 1
                time.sleep(max(1.0, float(args.watch)))
                if not args.as_json:
                    print("")
        except KeyboardInterrupt:
            return rc

    status = render_once()
    return 0 if status["found_any"] else 1


def _cmd_runs(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas runs",
        description="List recent Atlas runs and their dashboard status.",
    )
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args(list(argv))

    from cli.status_dashboard import build_runs_index, format_runs_table

    rows = build_runs_index(_bindings.repo_root(), limit=max(1, int(args.limit)))
    if args.as_json:
        print(json.dumps(rows, indent=2, sort_keys=True, default=str))
    else:
        print(format_runs_table(rows))
    return 0 if rows else 1


def _latest_run_id(root: Path) -> str | None:
    candidates: list[tuple[float, str]] = []
    for base in (root / "manifests", root / "data", root / "post_docked", root / "docked"):
        if not base.exists():
            continue
        for child in base.iterdir():
            if not child.is_dir():
                continue
            try:
                mtime = child.stat().st_mtime
            except OSError:
                continue
            candidates.append((mtime, child.name))
    if not candidates:
        return None
    candidates.sort()
    return candidates[-1][1]


def _read_manifest(root: Path, run_id: str) -> tuple[dict[str, Any] | None, Path | None]:
    try:
        from analysis.reporting.manifest_utils import load_run_manifest
    except Exception:
        return None, None
    return load_run_manifest(root, run_id)


def _build_status(root: Path, run_id: str) -> dict[str, Any]:
    manifest, manifest_path = _read_manifest(root, run_id)
    data_dir = run_output_dir(root, "data", run_id)
    paths = {
        "manifest": manifest_path,
        "report_html": data_dir / "report.html",
        "report_yaml": data_dir / "report.yaml",
        "heatmap_csv": data_dir / "heatmap_input.csv",
        "docked": run_output_dir(root, "docked", run_id),
        "post_docked": run_output_dir(root, "post_docked", run_id),
        "logs": output_root(root, "logs"),
        "configs": run_output_dir(root, "configs", run_id),
    }
    summary: dict[str, Any] = {}
    protein_status: Counter[str] = Counter()
    failed_examples: list[str] = []
    if isinstance(manifest, dict):
        raw_summary = manifest.get("summary")
        if isinstance(raw_summary, dict):
            summary.update(raw_summary)
        proteins = manifest.get("proteins")
        if isinstance(proteins, dict):
            for key, entry in proteins.items():
                if not isinstance(entry, dict):
                    continue
                status = str(entry.get("status") or "unknown")
                protein_status[status] += 1
                if status.lower() == "failed" and len(failed_examples) < 5:
                    failed_examples.append(str(key))

    existing_paths = {
        key: str(path)
        for key, path in paths.items()
        if path is not None and Path(path).exists()
    }
    return {
        "run_id": run_id,
        "found_any": bool(existing_paths or manifest),
        "manifest": str(manifest_path) if manifest_path else "",
        "summary": summary,
        "protein_status": dict(protein_status),
        "failed_examples": failed_examples,
        "paths": {key: str(value) for key, value in paths.items() if value is not None},
        "existing_paths": existing_paths,
        "next": _status_next_steps(paths),
    }


def _status_next_steps(paths: dict[str, Path | None]) -> list[str]:
    steps: list[str] = []
    html_path = paths.get("report_html")
    heatmap_path = paths.get("heatmap_csv")
    manifest_path = paths.get("manifest")
    if html_path is not None and html_path.exists():
        steps.append(f"Open report: {html_path}")
    elif heatmap_path is not None and heatmap_path.exists():
        steps.append("Generate report: atlas report <RUN_ID>")
    elif manifest_path is not None and manifest_path.exists():
        steps.append("Resume or inspect failures: atlas --resume --run-id <RUN_ID>")
    else:
        steps.append("Create a demo: atlas demo")
    return steps


def _print_status(status: dict[str, Any]) -> None:
    print(f"run_id: {status['run_id']}")
    if not status["found_any"]:
        print("status: not found")
        return
    summary = status.get("summary") or {}
    for key in (
        "total_proteins_scheduled",
        "total_proteins_completed",
        "total_proteins_failed",
    ):
        if key in summary:
            print(f"{key}: {summary[key]}")
    protein_status = status.get("protein_status") or {}
    if protein_status:
        print("protein_status: " + ", ".join(f"{k}={v}" for k, v in sorted(protein_status.items())))
    if status.get("failed_examples"):
        print("failed_examples: " + ", ".join(status["failed_examples"]))
    existing = status.get("existing_paths") or {}
    if existing:
        print("outputs:")
        for key in sorted(existing):
            print(f"  {key}: {existing[key]}")
    for step in status.get("next") or []:
        print(f"next: {step}")
