from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Sequence

from config.output_paths import run_output_dir

from cli.qol.ligands import _print_ligand_paths, _run_with_library
from cli.qol.pipeline import _cmd_run
from cli.qol.targets import _print_target_install_result
from cli.qol import _bindings

def _cmd_new_run(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas new-run",
        description="Validate first-run inputs and build a runnable Atlas command.",
    )
    parser.add_argument("--pdb", required=True, help="PDB ID matching input_pdbs/<ID>.pdb")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--fast", action="store_true", default=True)
    parser.add_argument("--no-fast", action="store_false", dest="fast")
    parser.add_argument("--small-library", action="store_true")
    parser.add_argument("--single")
    parser.add_argument("--dry-run", action="store_true", default=True)
    parser.add_argument("--no-dry-run", action="store_false", dest="dry_run")
    parser.add_argument("--yes", action="store_true", help="Run the generated command.")
    args = parser.parse_args(list(argv))

    pdb_id = str(args.pdb).strip().upper()
    root = _bindings.repo_root()
    cfg = _bindings.load_effective_config()
    input_dir = Path(str(cfg.get("INPUT_DIR") or root / "input_pdbs"))
    if not input_dir.is_absolute():
        input_dir = root / input_dir
    pdb_path = input_dir / f"{pdb_id}.pdb"
    command = ["atlas", "--pdb", pdb_id]
    if args.run_id:
        command.extend(["--run-id", args.run_id])
    if args.fast:
        command.append("--fast")
    if args.small_library:
        command.append("-test-fda")
    if args.single:
        command.extend(["--single", args.single])

    print(f"target: {pdb_id}")
    print(f"expected_input: {pdb_path}")
    if not pdb_path.exists():
        print("status: input_missing")
        print(f"next: place the structure at {pdb_path}")
        print("next: rerun this command before launching docking")
        return 1

    setup = _bindings.build_setup_report()
    if not setup["ready"]["vina_docking"]:
        print("status: tool_check_failed")
        print("next: run `atlas setup-report` and configure Vina/Open Babel")
        return 1
    if not setup["ready"]["full_receptor_prep"]:
        print("status: receptor_prep_tools_missing")
        print("next: install Meeko in the docking environment")
        return 1

    print("status: ready")
    print("command: " + " ".join(command))
    if args.yes and not args.dry_run:
        return _cmd_run(command[1:])
    print("next: rerun with `--yes --no-dry-run` when you want to launch it")
    return 0


def _cmd_first_run(argv: Sequence[str]) -> int:
    from cli.target_install import install_targets, known_target_panels
    from prep_ligands.ligand_library_manager import install_source, paths_for_source, resolve_source
    from tools import pdb_query

    parser = argparse.ArgumentParser(
        prog="atlas first-run",
        description="Select/install a first target, plan or prepare a ligand library, and build the launch command.",
    )
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--pdb", help="Use an already available input PDB ID.")
    target.add_argument("--panel", choices=[panel.name for panel in known_target_panels()])
    target.add_argument("--genes", nargs="+")
    target.add_argument("--genes-file")
    target.add_argument("--uniprot", action="append", default=[])
    parser.add_argument("--ligands", default="fda", help="Built-in ligand library source.")
    parser.add_argument("--no-ligand-install", action="store_true")
    parser.add_argument("--small-library", action="store_true")
    parser.add_argument("--single")
    parser.add_argument("--fast", action="store_true", default=True)
    parser.add_argument("--no-fast", action="store_false", dest="fast")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--dry-run", action="store_true", default=True)
    parser.add_argument("--no-dry-run", action="store_false", dest="dry_run")
    parser.add_argument("--yes", action="store_true", help="Launch after required install/prep steps.")
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--max-return", type=int, default=1)
    parser.add_argument("--max-install", type=int, default=1)
    parser.add_argument("--resolution-max", type=float)
    parser.add_argument("--species", "--organism")
    parser.add_argument("--taxonomy-id", "--taxon", dest="taxonomy_id", type=int)
    parser.add_argument("--entity-type")
    parser.add_argument("--methods")
    parser.add_argument("--allow-apo", action="store_true")
    parser.add_argument(
        "--contains-ligand",
        action="append",
        default=[],
        dest="target_ligands",
        help="Require a co-crystallized PDB chemical component ID, e.g. ATP or HEM.",
    )
    parser.add_argument("--dedupe", nargs=2, metavar=("MODE", "CUTOFF"))
    parser.add_argument(
        "--dedupe-sequence-identity",
        type=int,
        choices=pdb_query.SEQUENCE_IDENTITY_LEVELS,
    )
    parser.add_argument("--no-dedupe", action="store_true")
    parser.add_argument("--quality", choices=pdb_query.QUALITY_PRESETS, default="any")
    parser.add_argument("--category", "--gene-type", dest="category")
    parser.add_argument("--category-limits-json")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument(
        "--ligand-filter-mode",
        default="strict",
        choices=pdb_query.LIGAND_FILTER_MODES,
        help="How to treat borderline co-crystallized non-polymer ligands.",
    )
    parser.add_argument("--source-sdf")
    parser.add_argument("--source-url")
    parser.add_argument("--force-fetch", action="store_true")
    parser.add_argument("--force-prep", action="store_true")
    parser.add_argument("--limit-ligands", type=int, default=0)
    args = parser.parse_args(list(argv))

    cfg = _bindings.load_effective_config()
    root = _bindings.repo_root()
    target_result: Any | None = None
    if args.pdb:
        selected_ids = [str(args.pdb).strip().upper()]
    else:
        try:
            dedupe_sequence_identity = pdb_query.parse_dedupe_option(
                args.dedupe,
                args.dedupe_sequence_identity,
                no_dedupe=args.no_dedupe,
            )
            target_result = install_targets(
                cfg,
                panel_name=args.panel,
                genes=tuple(args.genes or ()),
                genes_file=Path(args.genes_file).expanduser() if args.genes_file else None,
                selected_out=_first_run_plan_selected_out(root, args.run_id)
                if args.dry_run
                else None,
                manifest_out=_first_run_plan_manifest_out(root, args.run_id)
                if args.dry_run
                else None,
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
                uniprots=tuple(args.uniprot or ()),
                ligand_comp_ids=tuple(args.target_ligands or ()),
                ligand_filter_mode=args.ligand_filter_mode,
                dedupe_sequence_identity=dedupe_sequence_identity,
                no_dedupe=args.no_dedupe,
                quality=args.quality,
                category=args.category,
                category_limits_json=args.category_limits_json,
                no_cache=args.no_cache,
            )
        except (RuntimeError, ValueError) as exc:
            print(f"atlas first-run: {exc}", file=sys.stderr)
            return 2
        _print_target_install_result(target_result)
        selected_ids = list(target_result.installed_ids or target_result.selected_ids)
    if not selected_ids:
        print("atlas first-run: no target PDB IDs were selected", file=sys.stderr)
        return 1

    pdb_id = selected_ids[0]
    source = resolve_source(args.ligands)
    ligand_paths = paths_for_source(cfg, source)
    if args.dry_run:
        _print_ligand_paths("planned", source.name, ligand_paths)
    elif args.no_ligand_install:
        _print_ligand_paths("using_existing", source.name, ligand_paths)
    else:
        try:
            ligand_paths = install_source(
                source.name,
                cfg,
                source_sdf=Path(args.source_sdf).expanduser() if args.source_sdf else None,
                source_url=args.source_url,
                force_fetch=args.force_fetch,
                force_prep=args.force_prep,
                limit=int(args.limit_ligands or 0) or None,
            )
        except RuntimeError as exc:
            print(f"atlas first-run: {exc}", file=sys.stderr)
            return 2
        _print_ligand_paths("installed", source.name, ligand_paths)

    pipeline_args = ["--pdb", pdb_id]
    if args.run_id:
        pipeline_args.extend(["--run-id", args.run_id])
    if args.fast:
        pipeline_args.append("--fast")
    if args.small_library:
        pipeline_args.append("-test-fda")
    if args.single:
        pipeline_args.extend(["--single", args.single])

    command = ["atlas", source.name, *pipeline_args]
    print("first_target: " + pdb_id)
    print("command: " + " ".join(command))
    if args.yes and not args.dry_run:
        return _run_with_library(source.name, pipeline_args)
    if args.dry_run:
        print("next: rerun with `--yes --no-dry-run` to install inputs and launch")
    else:
        print("next: rerun with `--yes` to launch the prepared run")
    return 0


def _first_run_plan_dir(root: Path, run_id: str | None) -> Path:
    token = str(run_id or "").strip() or "first_run_plan"
    return run_output_dir(root, "data", token) / "first_run_plan"


def _first_run_plan_selected_out(root: Path, run_id: str | None) -> Path:
    return _first_run_plan_dir(root, run_id) / "atlas_targets_custom.csv"


def _first_run_plan_manifest_out(root: Path, run_id: str | None) -> Path:
    return _first_run_plan_dir(root, run_id) / "atlas_target_install_manifest.json"
