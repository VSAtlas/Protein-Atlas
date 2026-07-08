from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Sequence


from cli.qol.ligands import _print_ligand_paths, _run_with_library
from cli.qol.targets import _add_target_query_filter_args, _print_target_install_result
from cli.qol import _bindings

def _cmd_run_panel(argv: Sequence[str]) -> int:
    parser = _build_run_panel_parser()
    args, pipeline_args = parser.parse_known_args(list(argv))

    cfg = _bindings.load_effective_config()
    try:
        result = _run_panel_install_targets(cfg, args)
        ligand_source = _run_panel_install_ligands(cfg, args)
    except (RuntimeError, ValueError) as exc:
        print(f"atlas run-panel: {exc}", file=sys.stderr)
        return 2

    forwarded = _run_panel_pipeline_args(pipeline_args, result)
    if args.dry_run:
        print(f"running_with_library: {ligand_source.name}")
        print("planned_pipeline_args: " + (" ".join(forwarded) if forwarded else "<none>"))
        return 0
    return _run_with_library(ligand_source.name, forwarded)


def _build_run_panel_parser() -> argparse.ArgumentParser:
    from cli.target_install import known_target_panels

    parser = argparse.ArgumentParser(
        prog="atlas run-panel",
        description="Install a target panel, install/prep a ligand library, then run Atlas.",
    )
    parser.add_argument("panel", choices=[panel.name for panel in known_target_panels()])
    parser.add_argument("--ligands", default="fda", help="Built-in ligand source token, for example chembl.")
    parser.add_argument("--dry-run", action="store_true", help="Plan target selection but do not download, prep, or run.")
    parser.add_argument("--no-target-install", action="store_true")
    parser.add_argument("--no-ligand-install", action="store_true")
    parser.add_argument("--target-out-dir")
    parser.add_argument("--selected-out")
    parser.add_argument("--target-manifest-out")
    parser.add_argument("--max-install", type=int)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--max-return", type=int)
    parser.add_argument("--resolution-max", type=float)
    parser.add_argument("--species", "--organism")
    parser.add_argument("--taxonomy-id", "--taxon", dest="taxonomy_id", type=int)
    parser.add_argument("--entity-type")
    parser.add_argument("--methods")
    parser.add_argument("--allow-apo", action="store_true")
    _add_target_query_filter_args(parser, ligand_dest="target_ligands")
    parser.add_argument("--category-limits-json")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--limit-ligands", type=int, default=0)
    parser.add_argument("--source-sdf")
    parser.add_argument("--source-url")
    parser.add_argument("--force-fetch", action="store_true")
    parser.add_argument("--force-prep", action="store_true")
    return parser


def _run_panel_install_targets(cfg: dict[str, Any], args: Any) -> Any:
    from cli.target_install import install_targets

    if args.no_target_install:
        return None
    from tools import pdb_query

    dedupe_sequence_identity = pdb_query.parse_dedupe_option(
        args.dedupe,
        args.dedupe_sequence_identity,
        no_dedupe=args.no_dedupe,
    )
    result = install_targets(
        cfg,
        panel_name=args.panel,
        selected_out=Path(args.selected_out).expanduser() if args.selected_out else None,
        manifest_out=Path(args.target_manifest_out).expanduser() if args.target_manifest_out else None,
        out_dir=Path(args.target_out_dir).expanduser() if args.target_out_dir else None,
        dry_run=args.dry_run,
        max_install=args.max_install,
        max_candidates=args.max_candidates,
        max_return=args.max_return,
        resolution_max=args.resolution_max,
        species=args.species,
        taxonomy_id=args.taxonomy_id,
        entity_type=args.entity_type,
        methods=args.methods,
        allow_apo=args.allow_apo,
        ligand_comp_ids=tuple(args.target_ligands or ()),
        dedupe_sequence_identity=dedupe_sequence_identity,
        no_dedupe=args.no_dedupe,
        ligand_filter_mode=args.ligand_filter_mode,
        quality=args.quality,
        category_limits_json=args.category_limits_json,
        no_cache=args.no_cache,
    )
    _print_target_install_result(result)
    return result


def _run_panel_install_ligands(cfg: dict[str, Any], args: Any) -> Any:
    from prep_ligands.ligand_library_manager import install_source, paths_for_source, resolve_source

    ligand_source = resolve_source(args.ligands)
    if args.dry_run:
        return ligand_source
    if args.no_ligand_install:
        _print_ligand_paths("using_existing", ligand_source.name, paths_for_source(cfg, ligand_source))
        return ligand_source
    ligand_paths = install_source(
        ligand_source.name,
        cfg,
        source_sdf=Path(args.source_sdf).expanduser() if args.source_sdf else None,
        source_url=args.source_url,
        force_fetch=args.force_fetch,
        force_prep=args.force_prep,
        limit=int(args.limit_ligands or 0) or None,
    )
    _print_ligand_paths("installed", ligand_source.name, ligand_paths)
    return ligand_source


def _run_panel_pipeline_args(pipeline_args: Sequence[str], result: Any) -> list[str]:
    forwarded = list(pipeline_args or [])
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    panel_ids = list((result.installed_ids or result.selected_ids) if result is not None else [])
    if panel_ids and not _pipeline_has_pdb_selection(forwarded):
        for pdb_id in panel_ids:
            forwarded.extend(["--pdb", pdb_id])
    return forwarded


def _pipeline_has_pdb_selection(args: Sequence[str]) -> bool:
    return any(token in {"--pdb", "--pdbs"} or token.startswith("--pdb=") or token.startswith("--pdbs=") for token in args)
