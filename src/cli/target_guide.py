"""Guided target gene-list, query, and install helpers for wetlab users."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from cli.target_install import (
    TargetPanel,
    clean_gene_tokens,
    install_targets,
    known_target_panels,
    resolve_target_panel,
    write_target_genes_csv,
)
from config.runtime_config import load_config
from tools import massinstall, pdb_query


DEFAULT_TARGET_DIR = Path("analysis") / "gene_list"
GUIDE_MODES = ("genes", "query", "dry-run", "install")


@dataclass(frozen=True)
class TargetGuideSpec:
    mode: str
    panel: TargetPanel | None
    genes: tuple[str, ...]
    category: str
    genes_out: Path
    selected_out: Path
    manifest_out: Path | None
    out_dir: Path | None
    species: str | None
    taxonomy_id: int | None
    entity_type: str
    methods: str
    allow_apo: bool
    ligand_comp_ids: tuple[str, ...]
    ligand_filter_mode: str
    dedupe_sequence_identity: int | None
    no_dedupe: bool
    quality: str
    resolution_max: float
    max_candidates: int
    max_return: int
    max_install: int | None
    no_cache: bool


def add_targets_genes_arguments(parser: argparse.ArgumentParser) -> None:
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--panel", choices=[panel.name for panel in known_target_panels()])
    source.add_argument("--genes", nargs="+")
    parser.add_argument("--category", "--gene-type", dest="category")
    parser.add_argument("--out")


def add_targets_guide_arguments(parser: argparse.ArgumentParser) -> None:
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--panel", choices=[panel.name for panel in known_target_panels()])
    source.add_argument("--genes", nargs="+")
    parser.add_argument("--mode", choices=GUIDE_MODES)
    parser.add_argument("--category", "--gene-type", dest="category")
    parser.add_argument("--genes-out")
    parser.add_argument("--selected-out")
    parser.add_argument("--manifest-out")
    parser.add_argument("--out-dir")
    parser.add_argument("--species", "--organism")
    parser.add_argument("--taxonomy-id", "--taxon", dest="taxonomy_id", type=int)
    parser.add_argument("--entity-type")
    parser.add_argument("--methods")
    parser.add_argument("--allow-apo", action="store_true")
    parser.add_argument("--ligand", "--contains-ligand", action="append", default=[], dest="ligands")
    parser.add_argument(
        "--ligand-filter-mode",
        default="strict",
        choices=pdb_query.LIGAND_FILTER_MODES,
    )
    parser.add_argument("--dedupe", nargs=2, metavar=("MODE", "CUTOFF"))
    parser.add_argument("--dedupe-sequence-identity", type=int)
    parser.add_argument("--no-dedupe", action="store_true")
    parser.add_argument("--quality", choices=pdb_query.QUALITY_PRESETS, default="any")
    parser.add_argument("--resolution-max", type=float)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--max-return", type=int)
    parser.add_argument("--max-install", type=int)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Use defaults for omitted guide prompts; useful for notebooks/scripts.",
    )


def cmd_targets_genes(args: Any) -> int:
    panel = resolve_target_panel(args.panel) if args.panel else None
    genes = panel.genes if panel else tuple(args.genes or ())
    category = panel.category if panel else args.category
    label = panel.name if panel else _category_label(category)
    out = Path(args.out).expanduser() if args.out else _default_genes_out(label)
    try:
        path = write_target_genes_csv(out, genes, category=category)
    except ValueError as exc:
        print(f"atlas targets genes: {exc}", file=sys.stderr)
        return 2
    print(f"Gene CSV: {path.resolve()}")
    print(f"Next: atlas targets install --genes-file {path}")
    return 0


def cmd_targets_guide(args: Any) -> int:
    try:
        spec = _resolve_guide_spec(args)
    except ValueError as exc:
        print(f"atlas targets guide: {exc}", file=sys.stderr)
        return 2

    genes_path = _write_guide_genes_csv(spec)
    print(f"Gene CSV: {genes_path.resolve()}")
    if spec.mode == "genes":
        print(f"Next: atlas targets from-genes --genes-file {genes_path}")
        return 0
    if spec.mode == "query":
        return _run_guide_query(spec, genes_path)
    return _run_guide_install(spec, genes_path)


def _resolve_guide_spec(args: Any) -> TargetGuideSpec:
    interactive = sys.stdin.isatty() and not bool(args.yes)
    panel = _resolve_guide_panel(args, interactive)
    genes = _resolve_guide_genes(args, panel, interactive)
    category = _resolve_category(args, panel, interactive)
    label = panel.name if panel else _category_label(category)
    mode = _resolve_mode(args, interactive)
    species, taxonomy_id = _resolve_species(args, interactive)
    return TargetGuideSpec(
        mode=mode,
        panel=panel,
        genes=tuple(genes),
        category=category,
        genes_out=_resolve_path(args.genes_out, _default_genes_out(label)),
        selected_out=_resolve_path(args.selected_out, _default_selected_out(label)),
        manifest_out=_optional_path(args.manifest_out),
        out_dir=_optional_path(args.out_dir),
        species=species,
        taxonomy_id=taxonomy_id,
        entity_type=_resolve_text(
            args.entity_type,
            "RCSB entity type",
            pdb_query.DEFAULT_ENTITY_TYPE,
            interactive,
        ),
        methods=_resolve_text(
            args.methods,
            "Experimental methods",
            pdb_query.DEFAULT_METHOD_FILTER,
            interactive,
        ),
        allow_apo=bool(args.allow_apo),
        ligand_comp_ids=pdb_query.normalize_comp_ids(args.ligands),
        ligand_filter_mode=str(args.ligand_filter_mode),
        dedupe_sequence_identity=pdb_query.parse_dedupe_option(
            args.dedupe,
            args.dedupe_sequence_identity,
            no_dedupe=bool(args.no_dedupe),
        ),
        no_dedupe=bool(args.no_dedupe),
        quality=pdb_query.normalize_quality_preset(args.quality),
        resolution_max=_resolve_float(
            args.resolution_max,
            "Maximum structure resolution (Angstrom)",
            pdb_query.DEFAULT_RESOLUTION_MAX,
            interactive,
        ),
        max_candidates=_resolve_int(args.max_candidates, "Candidate cap", 200, interactive),
        max_return=_resolve_int(args.max_return, "PDBs per gene", 10, interactive),
        max_install=args.max_install,
        no_cache=bool(args.no_cache),
    )


def _resolve_guide_panel(args: Any, interactive: bool) -> TargetPanel | None:
    if args.panel:
        return resolve_target_panel(args.panel)
    if args.genes or not interactive:
        return None
    panel_names = ", ".join(panel.name for panel in known_target_panels())
    answer = _prompt("Built-in panel or blank for custom genes", "", f"Options: {panel_names}")
    if not answer:
        return None
    return resolve_target_panel(answer)


def _resolve_guide_genes(
    args: Any,
    panel: TargetPanel | None,
    interactive: bool,
) -> list[str]:
    if panel is not None:
        return list(panel.genes)
    genes = clean_gene_tokens(args.genes or ())
    if genes:
        return genes
    if not interactive:
        raise ValueError("provide --genes or --panel, or run guide in an interactive shell")
    answer = _prompt("Gene symbols", "", "Comma or space separated, e.g. EGFR ABL1")
    genes = clean_gene_tokens([answer])
    if not genes:
        raise ValueError("provide at least one gene symbol")
    return genes


def _resolve_category(args: Any, panel: TargetPanel | None, interactive: bool) -> str:
    if panel is not None:
        return panel.category
    return _resolve_text(
        args.category,
        "Gene type/category label",
        massinstall.DEFAULT_CATEGORY,
        interactive,
    ).upper()


def _resolve_mode(args: Any, interactive: bool) -> str:
    if args.mode:
        return str(args.mode)
    if not interactive:
        return "dry-run"
    prompt = "Action: genes, query, dry-run, or install"
    return _prompt_choice(prompt, GUIDE_MODES, "dry-run")


def _write_guide_genes_csv(spec: TargetGuideSpec) -> Path:
    return write_target_genes_csv(spec.genes_out, spec.genes, category=spec.category)


def _run_guide_query(spec: TargetGuideSpec, genes_path: Path) -> int:
    forwarded = [
        "--genes-file",
        str(genes_path),
        "--out",
        str(spec.selected_out),
        "--max-candidates",
        str(spec.max_candidates),
        "--max-return",
        str(spec.max_return),
        "--resolution-max",
        str(spec.resolution_max),
        "--entity-type",
        spec.entity_type,
        "--methods",
        spec.methods,
    ]
    if spec.species is not None:
        forwarded.extend(["--species", spec.species])
    if spec.taxonomy_id is not None:
        forwarded.extend(["--taxonomy-id", str(spec.taxonomy_id)])
    if spec.allow_apo:
        forwarded.append("--allow-apo")
    forwarded.extend(["--ligand-filter-mode", spec.ligand_filter_mode])
    for ligand in spec.ligand_comp_ids:
        forwarded.extend(["--ligand", ligand])
    if spec.dedupe_sequence_identity is not None:
        forwarded.extend(["--dedupe-sequence-identity", str(spec.dedupe_sequence_identity)])
    if spec.no_dedupe:
        forwarded.append("--no-dedupe")
    forwarded.extend(["--quality", spec.quality])
    if spec.no_cache:
        forwarded.append("--no-cache")
    rc = int(pdb_query.main(forwarded))
    if rc == 0:
        print(f"Target CSV: {spec.selected_out.resolve()}")
        print(f"Next: atlas targets install --genes-file {genes_path}")
    return rc


def _run_guide_install(spec: TargetGuideSpec, genes_path: Path) -> int:
    cfg = load_config("config.txt", base_dir=_repo_root())
    try:
        result = install_targets(
            cfg,
            panel_name=spec.panel.name if spec.panel else None,
            genes_file=genes_path,
            selected_out=spec.selected_out,
            manifest_out=spec.manifest_out,
            out_dir=spec.out_dir,
            dry_run=spec.mode == "dry-run",
            max_install=spec.max_install,
            max_candidates=spec.max_candidates,
            max_return=spec.max_return,
            resolution_max=spec.resolution_max,
            species=spec.species,
            taxonomy_id=spec.taxonomy_id,
            entity_type=spec.entity_type,
            methods=spec.methods,
            allow_apo=spec.allow_apo,
            ligand_filter_mode=spec.ligand_filter_mode,
            ligand_comp_ids=spec.ligand_comp_ids,
            dedupe_sequence_identity=spec.dedupe_sequence_identity,
            no_dedupe=spec.no_dedupe,
            quality=spec.quality,
            category=spec.category,
            no_cache=spec.no_cache,
        )
    except (RuntimeError, ValueError) as exc:
        print(f"atlas targets guide: {exc}", file=sys.stderr)
        return 2
    _print_guide_result(result)
    return 0


def _print_guide_result(result: Any) -> None:
    print(f"status: {'selected' if result.dry_run else 'installed'}")
    print(f"selected_count: {len(result.selected_ids)}")
    print(f"installed_count: {len(result.installed_ids)}")
    print(f"selected_out: {result.selected_out}")
    print(f"output_dir: {result.output_dir}")
    print(f"manifest: {result.manifest_path}")
    if result.selected_ids:
        print("selected_ids: " + ", ".join(result.selected_ids[:20]))


def _resolve_species(args: Any, interactive: bool) -> tuple[str | None, int | None]:
    species = _resolve_text(
        args.species,
        "Species/organism",
        pdb_query.DEFAULT_SPECIES,
        interactive,
    )
    inferred_species, inferred_taxonomy_id = pdb_query.resolve_species_filter(species, None)
    if args.taxonomy_id is not None:
        return inferred_species, int(args.taxonomy_id)
    if not interactive:
        return inferred_species, inferred_taxonomy_id
    default = "" if inferred_taxonomy_id is None else str(inferred_taxonomy_id)
    while True:
        raw = _prompt("NCBI taxonomy id", default, "blank uses organism alias/default")
        if not raw:
            return inferred_species, inferred_taxonomy_id
        try:
            return inferred_species, int(raw)
        except ValueError:
            print("Enter an integer.", file=sys.stderr)


def _resolve_text(value: str | None, label: str, default: str, interactive: bool) -> str:
    if value is not None:
        clean = str(value).strip()
        if clean:
            return clean
    if interactive:
        return _prompt(label, default)
    return default


def _resolve_int(value: int | None, label: str, default: int, interactive: bool) -> int:
    if value is not None:
        return int(value)
    if not interactive:
        return default
    while True:
        raw = _prompt(label, str(default))
        try:
            return int(raw)
        except ValueError:
            print("Enter an integer.", file=sys.stderr)


def _resolve_float(value: float | None, label: str, default: float, interactive: bool) -> float:
    if value is not None:
        return float(value)
    if not interactive:
        return default
    while True:
        raw = _prompt(label, str(default))
        try:
            return float(raw)
        except ValueError:
            print("Enter a number.", file=sys.stderr)


def _prompt(label: str, default: str, hint: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    extra = f" ({hint})" if hint else ""
    answer = input(f"{label}{suffix}{extra}: ").strip()
    return answer or default


def _prompt_choice(label: str, choices: Sequence[str], default: str) -> str:
    choice_set = {choice.lower(): choice for choice in choices}
    while True:
        answer = _prompt(label, default, "/".join(choices)).lower()
        if answer in choice_set:
            return choice_set[answer]
        print(f"Choose one of: {', '.join(choices)}", file=sys.stderr)


def _resolve_path(value: str | None, default: Path) -> Path:
    return Path(value).expanduser() if value else default


def _optional_path(value: str | None) -> Path | None:
    return Path(value).expanduser() if value else None


def _default_genes_out(label: str) -> Path:
    return DEFAULT_TARGET_DIR / f"atlas_genes_{label}.csv"


def _default_selected_out(label: str) -> Path:
    return DEFAULT_TARGET_DIR / f"atlas_targets_{label}.csv"


def _category_label(category: str | None) -> str:
    clean = (category or massinstall.DEFAULT_CATEGORY).strip().lower()
    return clean.replace(" ", "_").replace("/", "_")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


__all__ = [
    "add_targets_genes_arguments",
    "add_targets_guide_arguments",
    "cmd_targets_genes",
    "cmd_targets_guide",
]
