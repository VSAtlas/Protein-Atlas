from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any, Sequence


from cli.qol.pipeline import _cmd_run
from cli.qol.utils import _display_path
from cli.qol import _bindings

def _cmd_ligands(argv: Sequence[str]) -> int:
    tokens = list(argv)
    if tokens[:2] == ["audit", "stereo"]:
        from analysis.cli.audit_ligand_stereo import main as audit_stereo

        return audit_stereo(tokens[2:])

    from prep_ligands.ligand_library_manager import known_sources, resolve_source

    source_choices = sorted(
        {
            token
            for source in known_sources()
            for token in (source.name, *source.aliases)
        }
    )
    parser = argparse.ArgumentParser(
        prog="atlas ligands",
        description="Download, prepare, and run Atlas ligand libraries.",
    )
    sub = parser.add_subparsers(dest="ligand_cmd", required=True)

    sub.add_parser("sources", help="List built-in ligand library sources.")
    audit = sub.add_parser(
        "audit",
        help="Audit a local ligand library without installing or downloading data.",
    )
    audit_sub = audit.add_subparsers(dest="audit_source", required=True)
    audit_sub.add_parser(
        "stereo",
        add_help=False,
        help="Report source-defined and unresolved ligand stereochemistry.",
    )

    for name in ("fetch", "prep", "install", "run"):
        cmd = sub.add_parser(name, help=f"{name} a ligand library source.")
        cmd.add_argument("source", choices=source_choices)
        cmd.add_argument("--source-sdf", help="Use a local SDF/SDF.GZ/ZIP instead of downloading.")
        cmd.add_argument("--source-url", help="Override the built-in source URL.")
        cmd.add_argument("--force-fetch", action="store_true")
        cmd.add_argument("--force-prep", action="store_true")
        cmd.add_argument("--limit", type=int, default=0, help="Prepare only the first N SDF records.")
        cmd.add_argument("--dry-run", action="store_true", help="Show planned paths and commands without fetching, prepping, or running.")
        if name == "run":
            cmd.add_argument("--no-install", action="store_true")

    args, pipeline_args = parser.parse_known_args(tokens)
    if args.ligand_cmd != "run" and pipeline_args:
        parser.error("unrecognized arguments: " + " ".join(pipeline_args))
    if args.ligand_cmd == "sources":
        return _print_ligand_sources()

    from prep_ligands.ligand_library_manager import (
        fetch_source,
        install_source,
        paths_for_source,
        prepare_source,
    )

    cfg = _bindings.load_effective_config()
    source = resolve_source(args.source)
    limit = int(args.limit or 0) or None
    source_sdf = Path(args.source_sdf).expanduser() if args.source_sdf else None
    if args.dry_run:
        paths = paths_for_source(cfg, source)
        _print_ligand_paths("planned", source.name, paths)
        print(f"action: {args.ligand_cmd}")
        if source_sdf:
            print(f"source_sdf: {_display_path(source_sdf)}")
        if args.source_url:
            print(f"source_url: {args.source_url}")
        if limit is not None:
            print(f"limit: {limit}")
        if args.ligand_cmd == "run":
            forwarded = _forwarded_pipeline_args(pipeline_args)
            print("planned_pipeline_args: " + (" ".join(forwarded) if forwarded else "<none>"))
        return 0

    try:
        if args.ligand_cmd == "fetch":
            paths = fetch_source(
                source.name,
                cfg,
                source_sdf=source_sdf,
                source_url=args.source_url,
                force=args.force_fetch,
                limit=limit,
            )
            _print_ligand_paths("fetched", source.name, paths)
            return 0
        if args.ligand_cmd == "prep":
            if source_sdf or args.source_url:
                fetch_source(
                    source.name,
                    cfg,
                    source_sdf=source_sdf,
                    source_url=args.source_url,
                    force=args.force_fetch,
                    limit=limit,
                )
            paths = prepare_source(source.name, cfg, force=args.force_prep, limit=limit)
            _print_ligand_paths("prepared", source.name, paths)
            return 0
        if args.ligand_cmd == "install":
            paths = install_source(
                source.name,
                cfg,
                source_sdf=source_sdf,
                source_url=args.source_url,
                force_fetch=args.force_fetch,
                force_prep=args.force_prep,
                limit=limit,
            )
            _print_ligand_paths("installed", source.name, paths)
            return 0
        if args.ligand_cmd == "run":
            paths = paths_for_source(cfg, source)
            if not args.no_install:
                paths = install_source(
                    source.name,
                    cfg,
                    source_sdf=source_sdf,
                    source_url=args.source_url,
                    force_fetch=args.force_fetch,
                    force_prep=args.force_prep,
                    limit=limit,
                )
                _print_ligand_paths("installed", source.name, paths)
            return _run_with_library(source.name, list(pipeline_args or []))
    except RuntimeError as exc:
        print(f"atlas ligands: {exc}", file=sys.stderr)
        return 2
    return 1


def _print_ligand_sources() -> int:
    from prep_ligands.ligand_library_manager import known_sources

    for source in known_sources():
        aliases = f" aliases={','.join(source.aliases)}" if source.aliases else ""
        print(f"{source.name}: {source.description}{aliases}")
        print(f"  library_subdir: {source.library_subdir}")
        print(f"  url: {source.url}")
        print(f"  license: {source.license_note}")
        print(f"  citation: {source.citation}")
    print("Use `atlas ligands install chembl`, `chebi`, `coconut`, `fda`, or `hmdb` to download and prepare.")
    return 0


def _print_ligand_paths(action: str, source_name: str, paths: Any) -> None:
    root = _bindings.repo_root()
    print(f"ligand_source: {source_name}")
    print(f"status: {action}")
    print(f"raw_sdf: {_display_path(paths.raw_sdf, root=root)}")
    print(f"library_dir: {_display_path(paths.library_dir, root=root)}")
    print(f"source_manifest: {_display_path(paths.source_manifest, root=root)}")


def _cmd_library_alias(source_name: str, argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog=f"atlas {source_name}",
        description=(
            f"Install the {source_name.upper()} ligand library if needed, then run Atlas "
            f"with TEST_MODE_ENABLE={source_name}."
        ),
    )
    parser.add_argument("--no-install", action="store_true")
    parser.add_argument("--source-sdf")
    parser.add_argument("--source-url")
    parser.add_argument("--force-fetch", action="store_true")
    parser.add_argument("--force-prep", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    args, pipeline_args = parser.parse_known_args(list(argv))

    from prep_ligands.ligand_library_manager import install_source, paths_for_source, resolve_source

    cfg = _bindings.load_effective_config()
    source = resolve_source(source_name)
    if args.dry_run:
        _print_ligand_paths("planned", source.name, paths_for_source(cfg, source))
        forwarded = _forwarded_pipeline_args(pipeline_args)
        print("planned_pipeline_args: " + (" ".join(forwarded) if forwarded else "<none>"))
        return 0
    if not args.no_install:
        try:
            paths = install_source(
                source.name,
                cfg,
                source_sdf=Path(args.source_sdf).expanduser() if args.source_sdf else None,
                source_url=args.source_url,
                force_fetch=args.force_fetch,
                force_prep=args.force_prep,
                limit=int(args.limit or 0) or None,
            )
        except RuntimeError as exc:
            print(f"atlas {source.name}: {exc}", file=sys.stderr)
            return 2
        _print_ligand_paths("installed", source.name, paths)
    else:
        _print_ligand_paths("using_existing", source.name, paths_for_source(cfg, source))
    return _run_with_library(source.name, list(pipeline_args or []))


def _forwarded_pipeline_args(pipeline_args: Sequence[str] | None) -> list[str]:
    """Drop an optional argparse-style ``--`` separator between known and forwarded args."""
    forwarded = list(pipeline_args or [])
    if forwarded and forwarded[0] == "--":
        forwarded.pop(0)
    return forwarded


def _run_with_library(source_name: str, pipeline_args: Sequence[str]) -> int:
    env_updates = {"TEST_MODE_ENABLE": source_name}
    old_env = {key: os.environ.get(key) for key in env_updates}
    try:
        for key, value in env_updates.items():
            os.environ[key] = value
        forwarded = _forwarded_pipeline_args(pipeline_args)
        print(f"running_with_library: {source_name}")
        print("pipeline_args: " + (" ".join(forwarded) if forwarded else "<none>"))
        return _cmd_run(forwarded)
    finally:
        for key, old_value in old_env.items():
            if old_value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = old_value
