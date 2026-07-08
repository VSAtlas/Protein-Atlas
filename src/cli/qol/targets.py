from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Sequence


from cli.qol.utils import _display_path
from cli.qol import _bindings

def _cmd_targets(argv: Sequence[str]) -> int:
    parser = _build_targets_parser()
    args = parser.parse_args(list(argv))
    if args.target_cmd == "panels":
        return _print_target_panels()
    if args.target_cmd == "species":
        return _print_target_species()
    if args.target_cmd == "from-genes":
        return _cmd_targets_from_genes(args)
    if args.target_cmd == "search":
        return _cmd_targets_search(args)
    if args.target_cmd == "genes":
        return _cmd_targets_genes(args)
    if args.target_cmd == "guide":
        return _cmd_targets_guide(args)
    if args.target_cmd == "install":
        return _cmd_targets_install(args)
    return _explain_target_csv(Path(args.csv_path))


def _build_targets_parser() -> argparse.ArgumentParser:
    from cli.target_install import known_target_panels
    from cli.target_guide import add_targets_genes_arguments, add_targets_guide_arguments
    from tools import pdb_query

    parser = argparse.ArgumentParser(
        prog="atlas targets",
        description="Wetlab-facing helpers for choosing structures before docking.",
    )
    sub = parser.add_subparsers(dest="target_cmd", required=True)
    sub.add_parser("panels", help="List built-in target panels.")
    sub.add_parser("species", help="List organism aliases for target queries.")
    from_genes = sub.add_parser("from-genes", help="Query RCSB by gene symbols.")
    from_genes.add_argument("genes", nargs="*")
    from_genes.add_argument("--genes-file")
    from_genes.add_argument("--out", default="analysis/gene_list/atlas_targets.csv")
    from_genes.add_argument("--max-candidates", type=int, default=200)
    from_genes.add_argument("--max-return", type=int, default=10)
    from_genes.add_argument("--resolution-max", type=float, default=3.2)
    from_genes.add_argument("--species", "--organism", default=None)
    from_genes.add_argument("--taxonomy-id", "--taxon", dest="taxonomy_id", type=int)
    from_genes.add_argument("--entity-type", default="Protein")
    from_genes.add_argument("--methods", default=pdb_query.DEFAULT_METHOD_FILTER)
    from_genes.add_argument("--allow-apo", action="store_true")
    from_genes.add_argument("--uniprot", action="append", default=[])
    _add_target_query_filter_args(from_genes)
    from_genes.add_argument("--no-cache", action="store_true")
    search = sub.add_parser("search", help="Full-text RCSB target search.")
    search.add_argument("query", nargs="+", help="Search text, e.g. BRCA DNA repair")
    search.add_argument("--out", default="analysis/gene_list/atlas_target_search.csv")
    search.add_argument("--max-candidates", type=int, default=200)
    search.add_argument("--max-return", type=int, default=10)
    search.add_argument("--resolution-max", type=float, default=3.2)
    search.add_argument("--species", "--organism", default=None)
    search.add_argument("--taxonomy-id", "--taxon", dest="taxonomy_id", type=int)
    search.add_argument("--entity-type", default="Protein")
    search.add_argument("--methods", default=pdb_query.DEFAULT_METHOD_FILTER)
    search.add_argument("--allow-apo", action="store_true")
    _add_target_query_filter_args(search)
    search.add_argument("--no-cache", action="store_true")
    genes = sub.add_parser("genes", help="Write a reusable gene/category CSV.")
    add_targets_genes_arguments(genes)
    guide = sub.add_parser("guide", aliases=["guided"], help="Prompt through target query/install setup.")
    add_targets_guide_arguments(guide)
    install = sub.add_parser("install", help="Select structures and download PDB files.")
    install.add_argument("gene_symbols", nargs="*", metavar="GENE")
    install_source = install.add_mutually_exclusive_group()
    install_source.add_argument("--panel", choices=[panel.name for panel in known_target_panels()])
    install_source.add_argument("--genes", nargs="+", dest="genes_option")
    install_source.add_argument("--genes-file")
    install.add_argument("--dry-run", action="store_true", help="Select targets but do not download PDB files.")
    install.add_argument("--install", action="store_true", help="Accepted for compatibility; install is the default.")
    install.add_argument("--out-dir", help="Destination for downloaded PDB files; defaults to configured INPUT_DIR.")
    install.add_argument("--selected-out", help="Selection CSV output path.")
    install.add_argument("--manifest-out", help="Target install manifest JSON path.")
    install.add_argument("--missing-out", help="Failed-download CSV path.")
    install.add_argument("--max-install", type=int)
    install.add_argument("--max-candidates", type=int)
    install.add_argument("--max-return", type=int)
    install.add_argument("--resolution-max", type=float)
    install.add_argument("--species", "--organism")
    install.add_argument("--taxonomy-id", "--taxon", dest="taxonomy_id", type=int)
    install.add_argument("--entity-type")
    install.add_argument("--methods")
    install.add_argument("--allow-apo", action="store_true")
    install.add_argument("--uniprot", action="append", default=[])
    _add_target_query_filter_args(install)
    install.add_argument("--category", "--gene-type", dest="category")
    install.add_argument("--category-limits-json")
    install.add_argument("--no-cache", action="store_true")
    explain = sub.add_parser("explain-selection", help="Summarize a target CSV.")
    explain.add_argument("csv_path")
    return parser


def _add_target_query_filter_args(
    parser: argparse.ArgumentParser,
    *,
    ligand_dest: str = "ligands",
) -> None:
    from tools import pdb_query

    parser.add_argument(
        "--ligand",
        "--contains-ligand",
        action="append",
        default=[],
        dest=ligand_dest,
        help="Require a co-crystallized PDB chemical component ID, e.g. ATP or HEM.",
    )
    parser.add_argument(
        "--dedupe",
        nargs=2,
        metavar=("MODE", "CUTOFF"),
        help="De-duplicate search hits, e.g. `--dedupe sequence-identity 90`.",
    )
    parser.add_argument(
        "--ligand-filter-mode",
        default="strict",
        choices=pdb_query.LIGAND_FILTER_MODES,
        help="How to treat borderline co-crystallized non-polymer ligands.",
    )
    parser.add_argument(
        "--dedupe-sequence-identity",
        type=int,
        choices=pdb_query.SEQUENCE_IDENTITY_LEVELS,
    )
    parser.add_argument(
        "--no-dedupe",
        action="store_true",
        help="Disable default sequence-identity de-duplication.",
    )
    parser.add_argument(
        "--quality",
        choices=pdb_query.QUALITY_PRESETS,
        default="any",
    )


def _cmd_targets_genes(args: Any) -> int:
    from cli.target_guide import cmd_targets_genes

    return cmd_targets_genes(args)


def _cmd_targets_guide(args: Any) -> int:
    from cli.target_guide import cmd_targets_guide

    return cmd_targets_guide(args)


def _cmd_targets_from_genes(args: Any) -> int:
    from tools import pdb_query

    try:
        species_name, taxonomy_id = pdb_query.resolve_species_filter(
            args.species,
            args.taxonomy_id,
        )
        ligand_comp_ids = pdb_query.normalize_comp_ids(args.ligands)
        dedupe_sequence_identity = pdb_query.parse_dedupe_option(
            args.dedupe,
            args.dedupe_sequence_identity,
            no_dedupe=args.no_dedupe,
        )
        quality = pdb_query.normalize_quality_preset(args.quality)
    except ValueError as exc:
        print(f"atlas targets from-genes: {exc}", file=sys.stderr)
        return 2

    genes = list(args.genes or [])
    if args.genes_file:
        genes.extend(pdb_query._read_gene_list(args.genes_file))
    genes = list(dict.fromkeys(gene.strip() for gene in genes if gene.strip()))
    uniprots = list(dict.fromkeys(pdb_query._upper_list(args.uniprot or [])))
    if not genes and not uniprots:
        print(
            "atlas targets from-genes: provide gene names, --genes-file, or --uniprot",
            file=sys.stderr,
        )
        return 2

    client = pdb_query.PDBQueryClient(use_cache=not args.no_cache)
    rows: list[dict[str, str]] = []
    for gene in genes:
        result = pdb_query.get_ranked_entries_for_gene(
            gene,
            max_candidates=args.max_candidates,
            max_return=args.max_return,
            resolution_max=args.resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=args.entity_type,
            experimental_methods=args.methods,
            require_ligand=not args.allow_apo,
            ligand_comp_ids=ligand_comp_ids,
            dedupe_sequence_identity=dedupe_sequence_identity,
            quality=quality,
            ligand_filter_mode=args.ligand_filter_mode,
            client=client,
        )
        rows.extend(pdb_query._entries_to_rows(result.entries))
    for uniprot_id in uniprots:
        result = pdb_query.get_ranked_entries_for_uniprot(
            uniprot_id,
            max_candidates=args.max_candidates,
            max_return=args.max_return,
            resolution_max=args.resolution_max,
            species_name=species_name,
            taxonomy_id=taxonomy_id,
            entity_type=args.entity_type,
            experimental_methods=args.methods,
            require_ligand=not args.allow_apo,
            ligand_comp_ids=ligand_comp_ids,
            dedupe_sequence_identity=dedupe_sequence_identity,
            quality=quality,
            ligand_filter_mode=args.ligand_filter_mode,
            client=client,
        )
        rows.extend(pdb_query._entries_to_rows(result.entries))
    pdb_query._write_output_csv(args.out, rows)
    print(f"Target CSV: {Path(args.out).resolve()}")
    print("Next: review with `atlas targets explain-selection <csv>`.")
    return 0


def _cmd_targets_search(args: Any) -> int:
    from tools import pdb_query

    try:
        species_name, taxonomy_id = pdb_query.resolve_species_filter(
            args.species,
            args.taxonomy_id,
        )
        ligand_comp_ids = pdb_query.normalize_comp_ids(args.ligands)
        dedupe_sequence_identity = pdb_query.parse_dedupe_option(
            args.dedupe,
            args.dedupe_sequence_identity,
            no_dedupe=args.no_dedupe,
        )
        quality = pdb_query.normalize_quality_preset(args.quality)
    except ValueError as exc:
        print(f"atlas targets search: {exc}", file=sys.stderr)
        return 2

    client = pdb_query.PDBQueryClient(use_cache=not args.no_cache)
    query_text = " ".join(str(token) for token in args.query).strip()
    result = pdb_query.get_ranked_entries_for_text(
        query_text,
        max_candidates=args.max_candidates,
        max_return=args.max_return,
        resolution_max=args.resolution_max,
        species_name=species_name,
        taxonomy_id=taxonomy_id,
        entity_type=args.entity_type,
        experimental_methods=args.methods,
        require_ligand=not args.allow_apo,
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        quality=quality,
        ligand_filter_mode=args.ligand_filter_mode,
        client=client,
    )
    rows = pdb_query._entries_to_rows(result.entries)
    pdb_query._write_output_csv(args.out, rows)
    print(f"Target CSV: {Path(args.out).resolve()}")
    print(
        "Search summary: "
        f"candidates={result.stats['candidates_before_filters']} "
        f"selected={result.stats['selected_final']}"
    )
    print("Next: review with `atlas targets explain-selection <csv>`.")
    return 0


def _cmd_targets_install(args: Any) -> int:
    from cli.target_install import install_targets

    positional_genes = tuple(args.gene_symbols or ())
    if positional_genes and (args.panel or args.genes_option or args.genes_file):
        print(
            "atlas targets install: positional genes cannot be combined with --panel, --genes, or --genes-file",
            file=sys.stderr,
        )
        return 2
    genes = positional_genes or tuple(args.genes_option or ())
    try:
        from tools import pdb_query

        dedupe_sequence_identity = pdb_query.parse_dedupe_option(
            args.dedupe,
            args.dedupe_sequence_identity,
            no_dedupe=args.no_dedupe,
        )
        result = install_targets(
            _bindings.load_effective_config(),
            panel_name=args.panel,
            genes=genes,
            genes_file=Path(args.genes_file).expanduser() if args.genes_file else None,
            selected_out=Path(args.selected_out).expanduser() if args.selected_out else None,
            manifest_out=Path(args.manifest_out).expanduser() if args.manifest_out else None,
            missing_out=Path(args.missing_out).expanduser() if args.missing_out else None,
            out_dir=Path(args.out_dir).expanduser() if args.out_dir else None,
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
            ligand_comp_ids=tuple(args.ligands or ()),
            dedupe_sequence_identity=dedupe_sequence_identity,
            no_dedupe=args.no_dedupe,
            quality=args.quality,
            category=args.category,
            category_limits_json=args.category_limits_json,
            no_cache=args.no_cache,
        )
    except (RuntimeError, ValueError) as exc:
        print(f"atlas targets install: {exc}", file=sys.stderr)
        return 2
    _print_target_install_result(result)
    return 0


def _print_target_panels() -> int:
    from cli.target_install import known_target_panels

    for panel in known_target_panels():
        print(f"{panel.name}: {panel.description}")
        print(f"  category: {panel.category}")
        print(f"  genes: {', '.join(panel.genes)}")
        print(f"  citation: {panel.citation}")
    print("Use `atlas targets install --panel kinases` to select and download panel targets.")
    return 0


def _print_target_species() -> int:
    from tools.pdb_query import known_species_aliases

    for item in known_species_aliases():
        aliases = ", ".join(item.aliases[:4])
        print(f"{item.canonical}: taxonomy_id={item.taxonomy_id} aliases={aliases}")
    print("Use `--species mouse` or `--organism zebrafish`; pass `--species any` to disable organism filtering.")
    return 0


def _print_target_install_result(result: Any) -> None:
    root = _bindings.repo_root()
    print(f"target_panel: {result.panel}")
    print(f"status: {'selected' if result.dry_run else 'installed'}")
    print(f"selected_count: {len(result.selected_ids)}")
    print(f"installed_count: {len(result.installed_ids)}")
    print(f"failure_count: {len(result.failures)}")
    print(f"genes_file: {_display_path(result.genes_file, root=root)}")
    print(f"selected_out: {_display_path(result.selected_out, root=root)}")
    print(f"output_dir: {_display_path(result.output_dir, root=root)}")
    print(f"manifest: {_display_path(result.manifest_path, root=root)}")
    if result.selected_ids:
        print("selected_ids: " + ", ".join(result.selected_ids[:20]))


def _explain_target_csv(path: Path) -> int:
    if not path.exists():
        print(f"Target CSV not found: {path}", file=sys.stderr)
        return 1
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        print(f"No rows in {path}")
        return 1
    genes = Counter(row.get("gene", "") for row in rows)
    methods = Counter(row.get("method", "") for row in rows)
    print(f"rows: {len(rows)}")
    print("genes: " + ", ".join(f"{gene}={count}" for gene, count in sorted(genes.items())))
    if methods:
        print("methods: " + ", ".join(f"{method or 'unknown'}={count}" for method, count in sorted(methods.items())))
    print("top structures:")
    for row in rows[:10]:
        print(
            "  {gene} {pdb_id} resolution={resolution} ligands={ligands}".format(
                gene=row.get("gene", ""),
                pdb_id=row.get("pdb_id", ""),
                resolution=row.get("resolution", ""),
                ligands=row.get("nontrivial_comp_ids", ""),
            )
        )
    return 0
