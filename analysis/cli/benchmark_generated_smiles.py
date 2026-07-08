from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from analysis.ml.generative_adapters import (
    collect_result_files,
    default_guacamol_root,
    default_moses_root,
    environment_snapshot,
    extract_metrics,
    file_record,
    get_git_sha,
    normalize_command_remainder,
    repo_root,
    run_external_command,
    utc_now,
    write_benchmark_smiles,
    write_json,
    write_metrics_csv,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare and capture GuacaMol/MOSES benchmark artifacts for generated SMILES."
    )
    parser.add_argument("--suite", required=True, choices=["guacamol", "moses"])
    parser.add_argument("--generated-smiles", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--benchmark-input", type=Path, default=None)
    parser.add_argument("--benchmark-results", nargs="*", type=Path, default=[])
    parser.add_argument("--copy-results", action="store_true")
    parser.add_argument("--metrics-csv", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--smiles-column", default=None)
    parser.add_argument("--cwd", type=Path, default=None, help="Working directory for an explicitly executed command.")
    parser.add_argument("--timeout-sec", type=int, default=None)
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually launch the command after --command. Without this flag the CLI only captures existing files.",
    )
    parser.add_argument(
        "--command",
        nargs=argparse.REMAINDER,
        help="External benchmark command to run. Must be last; use '--command -- executable ...' if needed.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    command = normalize_command_remainder(args.command)
    if command and not args.execute:
        parser.error("--command was provided, but external execution requires --execute")
    if args.execute and not command:
        parser.error("--execute requires an explicit --command")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    benchmark_input = args.benchmark_input or (out_dir / f"{args.suite}_benchmark_input.smi")
    manifest_path = args.manifest or (out_dir / f"{args.suite}_benchmark_manifest.json")
    metrics_csv = args.metrics_csv or (out_dir / f"{args.suite}_metrics_summary.csv")
    default_root = default_guacamol_root() if args.suite == "guacamol" else default_moses_root()

    manifest: dict[str, Any] = {
        "adapter": "generated_smiles_benchmark_capture",
        "suite": args.suite,
        "created_at": utc_now(),
        "repo_root": str(repo_root()),
        "git_sha": get_git_sha(),
        "environment": environment_snapshot(),
        "external_import_policy": "stdlib_only_no_guacamol_or_moses_import",
        "tool_hints": {"default_root": file_record(default_root)},
        "input_preparation": None,
        "execution": None,
        "results": [],
        "metrics_summary": None,
    }

    capture_returncode = 0
    try:
        manifest["input_preparation"] = write_benchmark_smiles(
            args.generated_smiles,
            benchmark_input,
            smiles_column=args.smiles_column,
        )
    except Exception as exc:
        capture_returncode = 1
        manifest["input_preparation"] = {
            "status": "failed",
            "error": str(exc),
            "source": file_record(args.generated_smiles),
            "benchmark_smiles": file_record(benchmark_input),
        }

    execution_returncode = 0
    if args.execute:
        execution = run_external_command(
            command,
            out_dir=out_dir,
            stem=f"{args.suite}_benchmark",
            cwd=args.cwd or (default_root if default_root.exists() else None),
            timeout_sec=args.timeout_sec,
        )
        execution_returncode = int(execution["returncode"])
        manifest["execution"] = execution

    copy_to = out_dir / "captured_results" if args.copy_results else None
    result_records = collect_result_files(args.benchmark_results, copy_to=copy_to)
    manifest["results"] = result_records

    metric_sources: list[Path] = []
    for record in result_records:
        copied = record.get("copied_to")
        candidate = copied if isinstance(copied, dict) else record
        candidate_path = candidate.get("path")
        if candidate.get("exists") and candidate_path:
            metric_sources.append(Path(str(candidate_path)))
    metrics = extract_metrics(metric_sources, suite=args.suite)
    metrics_path = write_metrics_csv(metrics_csv, metrics)
    manifest["metrics_summary"] = {
        "metrics_csv": file_record(metrics_path),
        "metric_rows": len(metrics),
    }
    manifest["status"] = "failed" if execution_returncode != 0 or capture_returncode != 0 else "captured"

    write_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return execution_returncode or capture_returncode


if __name__ == "__main__":
    raise SystemExit(main())
