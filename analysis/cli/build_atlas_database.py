"""CLI for auditing or building an Atlas release database."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.atlas_database import (
    audit_release_inputs,
    build_release_database,
    load_release_manifest,
)
from config.output_paths import output_root


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _default_out_dir(repo_root: Path, manifest_path: Path) -> Path:
    manifest = load_release_manifest(manifest_path)
    return output_root(repo_root, "data") / str(manifest["release_id"])


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit or build a failure-complete Atlas release database."
    )
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("audit", "build"):
        child = sub.add_parser(command)
        child.add_argument("--manifest", required=True, type=Path)
        child.add_argument("--repo-root", type=Path, default=Path("."))
        child.add_argument("--out-dir", type=Path)
    sub.choices["audit"].add_argument("--strict", action="store_true")
    sub.choices["build"].add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    repo_root = args.repo_root.resolve()
    manifest_path = args.manifest.resolve()
    out_dir = (
        args.out_dir.resolve()
        if args.out_dir is not None
        else _default_out_dir(repo_root, manifest_path)
    )

    if args.command == "audit":
        summary = audit_release_inputs(manifest_path, repo_root)
        _write_json(out_dir / "audit.json", summary)
        strict_failure = bool(args.strict) and (
            bool(summary.get("errors"))
            or summary.get("failure_complete_run_count") != summary.get("run_count")
        )
        return_code = 2 if strict_failure or summary.get("errors") else 0
    else:
        summary = build_release_database(
            manifest_path,
            out_dir / "docking_atlas.sqlite",
            repo_root,
            overwrite=bool(args.overwrite),
            summary_path=out_dir / "build_summary.json",
        )
        return_code = 0

    print(json.dumps(summary, indent=2, sort_keys=True))
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
