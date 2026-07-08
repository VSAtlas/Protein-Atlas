from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from analysis.ml.generative_adapters import (
    default_generative_python,
    default_reinvent_root,
    environment_snapshot,
    file_record,
    get_git_sha,
    normalize_command_remainder,
    repo_root,
    run_external_command,
    standardize_smiles_file,
    utc_now,
    write_json,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Capture REINVENT4 generated molecules into Atlas ML manifests without importing REINVENT4."
    )
    parser.add_argument("--generated-smiles", type=Path, default=None, help="REINVENT4 SMILES/CSV/JSON output to capture.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--output-smiles", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--source-name", default="reinvent4")
    parser.add_argument("--smiles-column", default=None)
    parser.add_argument("--score-column", default=None)
    parser.add_argument("--config", type=Path, default=None, help="Optional REINVENT4 config recorded in the manifest.")
    parser.add_argument("--reinvent-root", type=Path, default=default_reinvent_root())
    parser.add_argument("--python-executable", type=Path, default=default_generative_python())
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
        help="External REINVENT4 command to run. Must be last; use '--command -- executable ...' if needed.",
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
    if not args.generated_smiles and not args.execute:
        parser.error("capture mode requires --generated-smiles")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    output_smiles = args.output_smiles or (out_dir / "generated_smiles.csv")
    manifest_path = args.manifest or (out_dir / "reinvent_generation_manifest.json")

    manifest: dict[str, Any] = {
        "adapter": "reinvent4_generation_capture",
        "created_at": utc_now(),
        "repo_root": str(repo_root()),
        "git_sha": get_git_sha(),
        "environment": environment_snapshot(),
        "external_import_policy": "stdlib_only_no_reinvent_import",
        "inputs": {
            "generated_smiles": file_record(args.generated_smiles),
            "config": file_record(args.config),
        },
        "tool_hints": {
            "reinvent_root": file_record(args.reinvent_root),
            "python_executable": file_record(args.python_executable),
        },
        "execution": None,
        "capture": None,
    }

    execution_returncode = 0
    if args.execute:
        cwd = args.cwd or (args.reinvent_root if args.reinvent_root.exists() else None)
        execution = run_external_command(
            command,
            out_dir=out_dir,
            stem="reinvent4_generation",
            cwd=cwd,
            timeout_sec=args.timeout_sec,
        )
        execution_returncode = int(execution["returncode"])
        manifest["execution"] = execution

    capture_returncode = 0
    if args.generated_smiles:
        try:
            manifest["capture"] = standardize_smiles_file(
                args.generated_smiles,
                output_smiles,
                source_name=args.source_name,
                smiles_column=args.smiles_column,
                score_column=args.score_column,
            )
        except Exception as exc:
            capture_returncode = 1
            manifest["capture"] = {
                "status": "failed",
                "error": str(exc),
                "standardized_smiles": file_record(output_smiles),
            }
    else:
        manifest["capture"] = {"status": "skipped", "reason": "no --generated-smiles provided"}

    manifest["status"] = (
        "failed"
        if execution_returncode != 0 or capture_returncode != 0
        else str((manifest.get("capture") or {}).get("status", "completed"))
    )
    write_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return execution_returncode or capture_returncode


if __name__ == "__main__":
    raise SystemExit(main())
