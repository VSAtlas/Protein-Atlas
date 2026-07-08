"""User-facing target panel selection and PDB installation helpers."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from tools import massinstall


@dataclass(frozen=True)
class TargetPanel:
    name: str
    description: str
    category: str
    genes: tuple[str, ...]
    citation: str


@dataclass(frozen=True)
class TargetInstallResult:
    panel: str
    genes_file: Path
    selected_out: Path
    manifest_path: Path
    output_dir: Path
    selected_ids: list[str]
    installed_ids: list[str]
    failures: list[tuple[str, int]]
    dry_run: bool


TARGET_PANELS: dict[str, TargetPanel] = {
    "kinases": TargetPanel(
        name="kinases",
        description="Common druggable kinase and kinase-safety targets.",
        category="KINASE",
        genes=("EGFR", "ERBB2", "ABL1", "SRC", "JAK2", "BRAF", "MAPK1", "CDK2", "KIT", "MET"),
        citation="Atlas built-in kinase panel; RCSB structures selected at install time.",
    ),
    "gpcrs": TargetPanel(
        name="gpcrs",
        description="Representative GPCR pharmacology and safety targets.",
        category="GPCR",
        genes=("ADRB2", "DRD2", "HRH1", "HTR2A", "HTR2B", "OPRM1", "ADORA2A", "CXCR4", "CCR5", "AGTR1"),
        citation="Atlas built-in GPCR panel; RCSB structures selected at install time.",
    ),
    "ion-channels": TargetPanel(
        name="ion-channels",
        description="Ion-channel and excitable-tissue safety targets.",
        category="ION_CHANNEL",
        genes=("KCNH2", "SCN5A", "CACNA1C", "GABRA1", "GRIN1", "HTR3A", "KCNA1", "TRPV1"),
        citation="Atlas built-in ion-channel panel; RCSB structures selected at install time.",
    ),
    "nuclear-receptors": TargetPanel(
        name="nuclear-receptors",
        description="Nuclear receptor pharmacology and endocrine-safety targets.",
        category="NUCLEAR_RECEPTOR",
        genes=("ESR1", "AR", "PPARG", "NR3C1", "VDR", "RXRA", "RARA", "PGR"),
        citation="Atlas built-in nuclear receptor panel; RCSB structures selected at install time.",
    ),
    "adme": TargetPanel(
        name="adme",
        description="ADME-relevant enzymes, transporters, and plasma binding targets.",
        category="ADME",
        genes=("CYP3A4", "CYP2D6", "CYP2C9", "ABCB1", "SLCO1B1", "SLC22A1", "UGT1A1", "CES1", "ALB"),
        citation="Atlas built-in ADME panel; RCSB structures selected at install time.",
    ),
}


def known_target_panels() -> list[TargetPanel]:
    return [TARGET_PANELS[key] for key in sorted(TARGET_PANELS)]


def resolve_target_panel(name: str) -> TargetPanel:
    token = str(name).strip().lower().replace("_", "-")
    aliases = {
        "ionchannels": "ion-channels",
        "ion-channels": "ion-channels",
        "nuclear": "nuclear-receptors",
        "nuclear-receptor": "nuclear-receptors",
        "nuclear-receptors": "nuclear-receptors",
    }
    token = aliases.get(token, token)
    try:
        return TARGET_PANELS[token]
    except KeyError as exc:
        choices = ", ".join(sorted(TARGET_PANELS))
        raise ValueError(f"unknown target panel '{name}'. Known panels: {choices}") from exc


def clean_gene_tokens(genes: Sequence[str]) -> list[str]:
    clean: list[str] = []
    for raw_gene in genes:
        for token in str(raw_gene).replace(",", " ").replace(";", " ").split():
            gene = token.strip().upper()
            if gene and gene not in clean:
                clean.append(gene)
    return clean


def write_target_genes_csv(
    path: Path,
    genes: Sequence[str],
    *,
    category: str | None = None,
) -> Path:
    clean_genes = clean_gene_tokens(genes)
    if not clean_genes:
        raise ValueError("provide at least one gene symbol")
    category_label = (category or massinstall.DEFAULT_CATEGORY).strip().upper()
    out = path.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "category"])
        writer.writeheader()
        for gene in clean_genes:
            writer.writerow({"gene": gene, "category": category_label})
    return out


def install_targets(
    cfg: Mapping[str, Any],
    *,
    panel_name: str | None = None,
    genes: Sequence[str] = (),
    genes_file: Path | None = None,
    selected_out: Path | None = None,
    manifest_out: Path | None = None,
    missing_out: Path | None = None,
    out_dir: Path | None = None,
    dry_run: bool = False,
    max_install: int | None = None,
    max_candidates: int | None = None,
    max_return: int | None = None,
    resolution_max: float | None = None,
    species: str | None = None,
    taxonomy_id: int | None = None,
    entity_type: str | None = None,
    methods: str | None = None,
    allow_apo: bool = False,
    ligand_filter_mode: str = "strict",
    uniprots: Sequence[str] = (),
    ligand_comp_ids: Sequence[str] = (),
    dedupe_sequence_identity: int | None = None,
    no_dedupe: bool = False,
    quality: str | None = None,
    category: str | None = None,
    category_limits_json: str | None = None,
    no_cache: bool = False,
) -> TargetInstallResult:
    panel = resolve_target_panel(panel_name) if panel_name else None
    label = panel.name if panel else "custom"
    selected_path = _resolve_selected_out(selected_out, label)
    genes_path = _resolve_genes_file(
        cfg=cfg,
        panel=panel,
        genes=genes,
        uniprots=uniprots,
        genes_file=genes_file,
        selected_out=selected_path,
        category=category,
    )
    output_dir = _resolve_output_dir(cfg, out_dir)
    manifest_path = _resolve_manifest_out(manifest_out, output_dir)
    missing_path = missing_out or manifest_path.with_name("atlas_target_install_missing.csv")

    args = _massinstall_args(
        genes_file=genes_path,
        selected_out=selected_path,
        missing_out=missing_path,
        out_dir=output_dir,
        dry_run=dry_run,
        max_install=max_install,
        max_candidates=max_candidates,
        max_return=max_return,
        resolution_max=resolution_max,
        species=species,
        taxonomy_id=taxonomy_id,
        entity_type=entity_type,
        methods=methods,
        allow_apo=allow_apo,
        ligand_filter_mode=ligand_filter_mode,
        uniprots=uniprots,
        ligand_comp_ids=ligand_comp_ids,
        dedupe_sequence_identity=dedupe_sequence_identity,
        no_dedupe=no_dedupe,
        quality=quality,
        category_limits_json=category_limits_json,
        no_cache=no_cache,
    )
    selected_rows, selected_ids = massinstall._select_from_pdb_query(cfg=cfg, args=args)
    massinstall._write_selection_csv(str(selected_path), selected_rows)

    installed_ids: list[str] = []
    failures: list[tuple[str, int]] = []
    if not dry_run:
        installed_ids, failures = massinstall._download_pdbs(
            selected_ids,
            output_dir=output_dir,
            max_install=max_install,
        )
        if failures:
            massinstall._write_missing_csv(str(missing_path), failures)

    _write_manifest(
        manifest_path,
        cfg=cfg,
        panel=panel,
        label=label,
        genes_path=genes_path,
        selected_path=selected_path,
        output_dir=output_dir,
        selected_rows=selected_rows,
        selected_ids=selected_ids,
        installed_ids=installed_ids,
        failures=failures,
        dry_run=dry_run,
        args=args,
    )
    return TargetInstallResult(
        panel=label,
        genes_file=genes_path,
        selected_out=selected_path,
        manifest_path=manifest_path,
        output_dir=output_dir,
        selected_ids=list(selected_ids),
        installed_ids=list(installed_ids),
        failures=list(failures),
        dry_run=dry_run,
    )


def _resolve_selected_out(selected_out: Path | None, label: str) -> Path:
    if selected_out is not None:
        return selected_out.expanduser()
    return Path("analysis") / "gene_list" / f"atlas_targets_{label}.csv"


def _resolve_genes_file(
    *,
    cfg: Mapping[str, Any],
    panel: TargetPanel | None,
    genes: Sequence[str],
    uniprots: Sequence[str],
    genes_file: Path | None,
    selected_out: Path,
    category: str | None = None,
) -> Path:
    if genes_file is not None:
        return genes_file.expanduser()
    rows: list[tuple[str, str]]
    if panel is not None:
        rows = [(gene, panel.category) for gene in panel.genes]
    else:
        clean_genes = clean_gene_tokens(genes)
        clean_uniprots = clean_gene_tokens(uniprots)
        if not clean_genes and clean_uniprots:
            rows = [
                (accession, massinstall.DEFAULT_CATEGORY)
                for accession in clean_uniprots
            ]
        elif not clean_genes:
            args = _massinstall_args(
                genes_file="",
                selected_out=selected_out,
                missing_out=selected_out.with_name("atlas_target_install_missing.csv"),
                out_dir=_resolve_output_dir(cfg, None),
                dry_run=True,
                max_install=None,
                max_candidates=None,
                max_return=None,
                resolution_max=None,
                species=None,
                taxonomy_id=None,
                entity_type=None,
                methods=massinstall.DEFAULT_METHOD_FILTER,
                allow_apo=False,
                uniprots=(),
                ligand_comp_ids=(),
                dedupe_sequence_identity=None,
                no_dedupe=False,
                quality=None,
                ligand_filter_mode="strict",
                category_limits_json=None,
                no_cache=False,
            )
            return massinstall._resolve_genes_file(cfg, args)
        else:
            category_label = (category or massinstall.DEFAULT_CATEGORY).strip().upper()
            rows = [(gene, category_label) for gene in clean_genes]
    out = selected_out.with_name(f"{selected_out.stem}_genes.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "category"])
        writer.writeheader()
        for gene, category in rows:
            writer.writerow({"gene": gene, "category": category})
    return out


def _resolve_output_dir(cfg: Mapping[str, Any], out_dir: Path | None) -> Path:
    if out_dir is not None:
        return out_dir.expanduser()
    raw = cfg.get("PDB_QUERY_DEFAULT_OUT_DIR") or cfg.get("INPUT_DIR") or "input_pdbs"
    return Path(str(raw)).expanduser()


def _resolve_manifest_out(manifest_out: Path | None, output_dir: Path) -> Path:
    if manifest_out is not None:
        return manifest_out.expanduser()
    return output_dir / "atlas_target_install_manifest.json"


def _massinstall_args(
    *,
    genes_file: Path | str,
    selected_out: Path,
    missing_out: Path,
    out_dir: Path,
    dry_run: bool,
    max_install: int | None,
    max_candidates: int | None,
    max_return: int | None,
    resolution_max: float | None,
    species: str | None,
    taxonomy_id: int | None,
    entity_type: str | None,
    methods: str | None,
    allow_apo: bool,
    ligand_filter_mode: str,
    uniprots: Sequence[str],
    ligand_comp_ids: Sequence[str],
    dedupe_sequence_identity: int | None,
    no_dedupe: bool,
    quality: str | None,
    category_limits_json: str | None,
    no_cache: bool,
) -> massinstall.MassinstallArgs:
    return massinstall.MassinstallArgs(
        pdb_query=True,
        genes_file=str(genes_file),
        gene_uniprot_file=None,
        forced_csv=None,
        trivial_ligands_file=None,
        keep_ligands_file=None,
        max_candidates=max_candidates,
        max_return=max_return,
        resolution_max=resolution_max,
        species=species,
        taxonomy_id=taxonomy_id,
        entity_type=entity_type,
        methods=methods,
        allow_apo=bool(allow_apo),
        ligand_filter_mode=ligand_filter_mode,
        uniprot=list(uniprots),
        ligands=list(ligand_comp_ids),
        dedupe=None,
        dedupe_sequence_identity=dedupe_sequence_identity,
        no_dedupe=bool(no_dedupe),
        quality=quality or "any",
        include_forced=None,
        category_limits_json=category_limits_json or "",
        cache_dir=None,
        cache_ttl_days=None,
        no_cache=bool(no_cache),
        rate_limit=None,
        timeout=None,
        max_retries=None,
        out_dir=str(out_dir),
        selected_out=str(selected_out),
        missing_out=str(missing_out),
        dry_run=bool(dry_run),
        install=not bool(dry_run),
        max_install=max_install,
    )


def _write_manifest(
    path: Path,
    *,
    cfg: Mapping[str, Any],
    panel: TargetPanel | None,
    label: str,
    genes_path: Path,
    selected_path: Path,
    output_dir: Path,
    selected_rows: Sequence[Mapping[str, str]],
    selected_ids: Sequence[str],
    installed_ids: Sequence[str],
    failures: Sequence[tuple[str, int]],
    dry_run: bool,
    args: Any,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "panel": asdict(panel) if panel else {"name": label, "description": "Custom gene list"},
        "genes_file": str(genes_path),
        "selected_out": str(selected_path),
        "output_dir": str(output_dir),
        "dry_run": dry_run,
        "selected_count": len(selected_ids),
        "installed_count": len(installed_ids),
        "failure_count": len(failures),
        "selected_ids": list(selected_ids),
        "installed_ids": list(installed_ids),
        "failures": [{"pdb_id": pdb_id, "status_code": status} for pdb_id, status in failures],
        "selection_rows": list(selected_rows),
        "query": {
            "max_candidates": args.max_candidates,
            "max_return": args.max_return,
            "resolution_max": args.resolution_max,
            "species": args.species,
            "taxonomy_id": args.taxonomy_id,
            "entity_type": args.entity_type or massinstall.DEFAULT_ENTITY_TYPE,
            "methods": args.methods or massinstall.DEFAULT_METHOD_FILTER,
            "require_ligand": not bool(args.allow_apo),
            "ligand_filter_mode": getattr(args, "ligand_filter_mode", "strict"),
            "uniprots": list(getattr(args, "uniprot", []) or []),
            "ligand_comp_ids": list(getattr(args, "ligands", []) or []),
            "dedupe_sequence_identity": args.dedupe_sequence_identity,
            "no_dedupe": bool(getattr(args, "no_dedupe", False)),
            "quality": args.quality,
            "category_limits_json": args.category_limits_json,
            "no_cache": args.no_cache,
        },
        "sources": {
            "selector": "tools.massinstall._select_from_pdb_query",
            "pdb_download_url_template": "https://files.rcsb.org/download/{pdb_id}.pdb",
            "input_dir": str(cfg.get("INPUT_DIR", "")),
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


__all__ = [
    "TARGET_PANELS",
    "TargetInstallResult",
    "TargetPanel",
    "clean_gene_tokens",
    "install_targets",
    "known_target_panels",
    "resolve_target_panel",
    "write_target_genes_csv",
]
