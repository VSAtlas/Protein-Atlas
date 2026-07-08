from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import requests  # type: ignore[import-untyped]
import gemmi

from config.runtime_config import load_inputs, validate_config
from tools.pdb_query import (
    DEFAULT_CACHE_DIR,
    DEFAULT_CACHE_TTL_SECONDS,
    DEFAULT_ENTITY_TYPE,
    DEFAULT_DEDUPE_SEQUENCE_IDENTITY,
    DEFAULT_METHOD_FILTER,
    DEFAULT_MAX_RETRIES,
    DEFAULT_RATE_LIMIT_RPS,
    DEFAULT_RESOLUTION_MAX,
    DEFAULT_TIMEOUT_SECONDS,
    FORCED_PDBS_BY_GENE,
    PDBQueryClient,
    RankedEntry,
    get_ranked_entries_for_gene,
    normalize_comp_ids,
    normalize_quality_preset,
    parse_dedupe_option,
    resolve_species_filter,
)

DEFAULT_CATEGORY = "DEFAULT"


@dataclass(slots=True)
class MassinstallArgs:
    pdb_query: bool = False
    genes_file: str = ""
    gene_uniprot_file: str | None = None
    forced_csv: str | None = None
    trivial_ligands_file: str | None = None
    keep_ligands_file: str | None = None
    max_candidates: int | None = None
    max_return: int | None = None
    resolution_max: float | None = None
    species: str | None = None
    taxonomy_id: int | None = None
    entity_type: str | None = None
    methods: str | None = None
    allow_apo: bool = False
    ligand_filter_mode: str = "strict"
    uniprot: list[str] = field(default_factory=list)
    ligands: list[str] = field(default_factory=list)
    dedupe: Sequence[str] | None = None
    dedupe_sequence_identity: int | None = None
    no_dedupe: bool = False
    quality: str = "any"
    include_forced: bool | None = None
    category_limits_json: str = ""
    cache_dir: str | None = None
    cache_ttl_days: float | None = None
    no_cache: bool = False
    rate_limit: float | None = None
    timeout: float | None = None
    max_retries: int | None = None
    out_dir: str = ""
    selected_out: str = ""
    missing_out: str = ""
    dry_run: bool = False
    install: bool = False
    max_install: int | None = None

    @classmethod
    def from_namespace(cls, args: argparse.Namespace) -> "MassinstallArgs":
        return cls(
            pdb_query=bool(getattr(args, "pdb_query", False)),
            genes_file=str(getattr(args, "genes_file", "") or ""),
            gene_uniprot_file=getattr(args, "gene_uniprot_file", None),
            forced_csv=getattr(args, "forced_csv", None),
            trivial_ligands_file=getattr(args, "trivial_ligands_file", None),
            keep_ligands_file=getattr(args, "keep_ligands_file", None),
            max_candidates=getattr(args, "max_candidates", None),
            max_return=getattr(args, "max_return", None),
            resolution_max=getattr(args, "resolution_max", None),
            species=getattr(args, "species", None),
            taxonomy_id=getattr(args, "taxonomy_id", None),
            entity_type=getattr(args, "entity_type", None),
            methods=getattr(args, "methods", None),
            allow_apo=bool(getattr(args, "allow_apo", False)),
            ligand_filter_mode=getattr(args, "ligand_filter_mode", "strict"),
            uniprot=list(getattr(args, "uniprot", []) or []),
            ligands=list(getattr(args, "ligands", []) or []),
            dedupe=getattr(args, "dedupe", None),
            dedupe_sequence_identity=getattr(args, "dedupe_sequence_identity", None),
            no_dedupe=bool(getattr(args, "no_dedupe", False)),
            quality=str(getattr(args, "quality", "any") or "any"),
            include_forced=getattr(args, "include_forced", None),
            category_limits_json=str(getattr(args, "category_limits_json", "") or ""),
            cache_dir=getattr(args, "cache_dir", None),
            cache_ttl_days=getattr(args, "cache_ttl_days", None),
            no_cache=bool(getattr(args, "no_cache", False)),
            rate_limit=getattr(args, "rate_limit", None),
            timeout=getattr(args, "timeout", None),
            max_retries=getattr(args, "max_retries", None),
            out_dir=str(getattr(args, "out_dir", "") or ""),
            selected_out=str(getattr(args, "selected_out", "") or ""),
            missing_out=str(getattr(args, "missing_out", "") or ""),
            dry_run=bool(getattr(args, "dry_run", False)),
            install=bool(getattr(args, "install", False)),
            max_install=getattr(args, "max_install", None),
        )

# Fallback list preserved for backwards compatibility.
FALLBACK_PDB_IDS = [
    "3EML",
    "2HZI",
    "3BKL",
    "1E66",
    "2E1W",
    "2OI0",
    "2VT4",
    "3NY8",
    "3CQW",
    "3D0E",
    "2HV5",
    "1L2S",
    "2AM9",
    "1S3B",
    "3L5D",
    "3D4Q",
    "1BCD",
    "2CNK",
    "1H00",
    "3BWM",
    "1R9O",
    "3NXU",
    "3KRJ",
    "3ODU",
    "1LRU",
    "3FRJ",
    "2I78",
    "3PBL",
    "3NXO",
    "2RGP",
    "1SJ0",
    "2FSZ",
    "3KL6",
    "1W7X",
    "2NNQ",
    "3BZ3",
    "3C4F",
    "1J4H",
    "3E37",
    "1ZW5",
    "3BQD",
    "2V3F",
    "3KGC",
    "1VSO",
    "3MAX",
    "3F07",
    "3NF7",
    "1XL2",
    "3LAN",
    "3CCW",
    "1UYG",
    "3F9M",
    "2OJ9",
    "2H7L",
    "2ICA",
    "3LPB",
    "3CJO",
    "3G0E",
    "2B8T",
    "2I0E",
    "2OF2",
    "3CHP",
    "3M2W",
    "2AA2",
    "3LQ8",
    "2OJG",
    "2ZDT",
    "2QD9",
    "830C",
    "3EQH",
    "1QW6",
    "1B9V",
    "1KVO",
    "3L3M",
    "1UDT",
    "2OYU",
    "3LN1",
    "2OWB",
    "3BGS",
    "2P54",
    "2ZNP",
    "2GTK",
    "3KBA",
    "2AZR",
    "1NJS",
    "1C8K",
    "1D3G",
    "3G6Z",
    "2ETR",
    "1MV9",
    "1LI4",
    "3EL8",
    "3HMM",
    "1Q4X",
    "1YPE",
    "2AYW",
    "2ZEC",
    "1SYN",
    "1SQT",
    "2P2I",
    "3BIZ",
    "3HL5",
]


def _upper_list(values: Sequence[str]) -> list[str]:
    return [value.strip().upper() for value in values if value.strip()]


def _merge_forced_maps(*maps: Mapping[str, Sequence[str]]) -> dict[str, list[str]]:
    merged: dict[str, list[str]] = {}
    for mapping in maps:
        for gene, pdb_ids in mapping.items():
            gene_u = gene.strip().upper()
            if not gene_u:
                continue
            merged.setdefault(gene_u, [])
            for pdb_id in _upper_list(list(pdb_ids)):
                if pdb_id not in merged[gene_u]:
                    merged[gene_u].append(pdb_id)
    return merged


def _read_list_file(path: Optional[str]) -> list[str]:
    if not path:
        return []
    file_path = Path(path)
    if not file_path.exists():
        return []
    rows = file_path.read_text(encoding="utf-8").splitlines()
    return _upper_list(rows)


def _read_comp_id_file(path: Optional[str]) -> set[str]:
    if not path:
        return set()
    file_path = Path(path)
    if not file_path.exists():
        return set()
    comp_ids: set[str] = set()
    for raw_line in file_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [token.strip() for token in re.split(r"[,\s;|]+", line) if token.strip()]
        comp_ids.update(_upper_list(parts))
    return comp_ids


def _read_gene_category_rows(path: str) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("genes-file must contain a CSV header")
        gene_col = next((col for col in reader.fieldnames if col and col.lower() == "gene"), None)
        category_col = next(
            (col for col in reader.fieldnames if col and col.lower() in {"category", "pdb_type", "type"}),
            None,
        )
        if gene_col is None:
            raise ValueError("genes-file must include a `gene` column")

        for row in reader:
            gene = (row.get(gene_col) or "").strip()
            if not gene:
                continue
            category = (row.get(category_col) or DEFAULT_CATEGORY).strip() if category_col else DEFAULT_CATEGORY
            normalized = (gene.upper(), category.upper() or DEFAULT_CATEGORY)
            if normalized not in seen:
                rows.append(normalized)
                seen.add(normalized)
    return rows


def _read_gene_to_uniprots(path: Optional[str]) -> dict[str, list[str]]:
    if not path:
        return {}
    file_path = Path(path)
    if not file_path.exists():
        return {}
    mapping: dict[str, list[str]] = {}
    with file_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        gene_col = next((c for c in reader.fieldnames if c and c.lower() == "gene"), None)
        uniprot_col = next(
            (
                c
                for c in reader.fieldnames
                if c and c.lower() in {"uniprot", "uniprot_id", "uniprots"}
            ),
            None,
        )
        if gene_col is None or uniprot_col is None:
            return {}
        for row in reader:
            gene = (row.get(gene_col) or "").strip().upper()
            raw = (row.get(uniprot_col) or "").strip()
            if not gene or not raw:
                continue
            ids = [tok for tok in re.split(r"[,\s;|]+", raw) if tok]
            mapping.setdefault(gene, [])
            for uniprot_id in _upper_list(ids):
                if uniprot_id not in mapping[gene]:
                    mapping[gene].append(uniprot_id)
    return mapping


def _read_forced_map(path: Optional[str]) -> dict[str, list[str]]:
    if not path:
        return {}
    file_path = Path(path)
    if not file_path.exists():
        return {}
    forced: dict[str, list[str]] = {}
    with file_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        gene_col = next((c for c in reader.fieldnames if c and c.lower() == "gene"), None)
        pdb_col = next(
            (c for c in reader.fieldnames if c and c.lower() in {"pdb_id", "pdb", "entry"}),
            None,
        )
        if gene_col is None or pdb_col is None:
            return {}
        for row in reader:
            gene = (row.get(gene_col) or "").strip().upper()
            pdb = (row.get(pdb_col) or "").strip().upper()
            if not gene or not pdb:
                continue
            forced.setdefault(gene, [])
            if pdb not in forced[gene]:
                forced[gene].append(pdb)
    return forced


def _parse_category_limits(raw_value: Any) -> dict[str, int]:
    if raw_value is None:
        return {}
    if isinstance(raw_value, dict):
        parsed = raw_value
    else:
        text = str(raw_value).strip()
        if not text:
            return {}
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            parsed = {}
            for token in re.split(r"[,\n]+", text):
                chunk = token.strip()
                if not chunk:
                    continue
                if "=" in chunk:
                    key, value = chunk.split("=", 1)
                elif ":" in chunk:
                    key, value = chunk.split(":", 1)
                else:
                    raise ValueError(
                        "Category limits must be JSON or KEY=VALUE pairs"
                    )
                parsed[key.strip()] = value.strip()
    if not isinstance(parsed, dict):
        raise ValueError("PDB_QUERY_CATEGORY_LIMITS must be a JSON object")
    limits: dict[str, int] = {}
    for key, value in parsed.items():
        category = str(key).strip().upper()
        if not category:
            continue
        limit = int(value)
        if limit < 0:
            raise ValueError(f"Category limit cannot be negative: {key}={value}")
        limits[category] = limit
    return limits


def _cfg_str(cfg: Mapping[str, Any], key: str) -> str:
    return str(cfg.get(key, "")).strip()


def _pick_str(arg_value: Optional[str], cfg: Mapping[str, Any], key: str) -> Optional[str]:
    if arg_value is not None:
        value = arg_value.strip()
    else:
        value = _cfg_str(cfg, key)
    return value or None


def _pick_int(
    arg_value: Optional[int],
    cfg: Mapping[str, Any],
    key: str,
    *,
    default: int,
) -> int:
    if arg_value is not None:
        return arg_value
    raw = _cfg_str(cfg, key)
    return int(raw) if raw else default


def _pick_float(
    arg_value: Optional[float],
    cfg: Mapping[str, Any],
    key: str,
    *,
    default: float,
) -> float:
    if arg_value is not None:
        return arg_value
    raw = _cfg_str(cfg, key)
    return float(raw) if raw else default


def _load_legacy_pdb_ids(cfg: Mapping[str, Any]) -> list[str]:
    pdb_ids = _read_list_file(str(cfg.get("PDB_ID_LIST", "")).strip())
    return pdb_ids if pdb_ids else list(FALLBACK_PDB_IDS)


def _resolve_genes_file(cfg: Mapping[str, Any], args: MassinstallArgs) -> Path:
    if args.genes_file:
        candidate = Path(args.genes_file.strip()).expanduser()
        if candidate.exists():
            return candidate
        raise ValueError(f"`--genes-file` does not exist: {candidate}")

    direct_uniprots = _upper_list(getattr(args, "uniprot", []) or [])
    if direct_uniprots:
        selected_out = Path(
            str(
                getattr(args, "selected_out", "")
                or "analysis/gene_list/atlas_targets_uniprot.csv"
            )
        )
        out = selected_out.with_name(f"{selected_out.stem}_genes.csv")
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["gene", "category"])
            writer.writeheader()
            for accession in direct_uniprots:
                writer.writerow({"gene": accession, "category": DEFAULT_CATEGORY})
        return out

    cfg_genes_file = _cfg_str(cfg, "PDB_QUERY_GENES_FILE")
    if cfg_genes_file:
        candidate = Path(cfg_genes_file).expanduser()
        if candidate.exists():
            return candidate
        raise ValueError(f"PDB_QUERY_GENES_FILE does not exist: {candidate}")

    default_candidate = Path("analysis/gene_list/pdb_query_genes_30.csv")
    if default_candidate.exists():
        return default_candidate

    raise ValueError(
        "`--pdb_query` needs a genes source. Provide --genes-file, set PDB_QUERY_GENES_FILE, "
        "use --uniprot, or add analysis/gene_list/pdb_query_genes_30.csv."
    )


def _download_pdbs(
    pdb_ids: Sequence[str],
    *,
    output_dir: Path,
    max_install: Optional[int] = None,
) -> tuple[list[str], list[tuple[str, int]]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    successes: list[str] = []
    failures: list[tuple[str, int]] = []

    install_ids = list(pdb_ids)
    if max_install is not None:
        install_ids = install_ids[:max_install]

    with requests.Session() as session:
        for pdb_id in install_ids:
            pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            try:
                response = session.get(pdb_url, timeout=60)
            except requests.RequestException:
                print(f"Failed to download {pdb_id}: request_error")
                failures.append((pdb_id, 0))
                continue

            if response.status_code == 200:
                (output_dir / f"{pdb_id}.pdb").write_bytes(response.content)
                print(f"Downloaded {pdb_id}")
                successes.append(pdb_id)
            elif response.status_code == 404:
                cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
                try:
                    cif_response = session.get(cif_url, timeout=60)
                except requests.RequestException:
                    print(f"Failed to download {pdb_id}: request_error")
                    failures.append((pdb_id, 0))
                    continue

                if cif_response.status_code != 200:
                    print(f"Failed to download {pdb_id}: {response.status_code} / {cif_response.status_code}")
                    failures.append((pdb_id, response.status_code))
                    continue

                try:
                    doc = gemmi.cif.read_string(cif_response.text)
                    structure = gemmi.make_structure_from_block(doc.sole_block())
                except Exception as exc:
                    print(f"Failed to convert CIF to PDB for {pdb_id}: {exc}")
                    failures.append((pdb_id, response.status_code))
                    continue

                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_handle:
                    tmp_path = tmp_handle.name

                try:
                    structure.write_pdb(tmp_path)
                    (output_dir / f"{pdb_id}.pdb").write_bytes(Path(tmp_path).read_bytes())
                    print(f"Downloaded {pdb_id} from CIF")
                    successes.append(pdb_id)
                finally:
                    if os.path.exists(tmp_path):
                        Path(tmp_path).unlink()
            else:
                print(f"Failed to download {pdb_id}: {response.status_code}")
                failures.append((pdb_id, response.status_code))
    return successes, failures


def _write_selection_csv(path: str, rows: Sequence[dict[str, str]]) -> None:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "category",
        "gene",
        "pdb_id",
        "rank_score",
        "resolution",
        "method",
        "forced_included",
        "valid_entry",
        "nontrivial_comp_ids",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Install PDB files from static list or PDB query")
    parser.add_argument("--pdb_query", "--pdb-query", action="store_true", dest="pdb_query")
    parser.add_argument(
        "--genes-file",
        help="CSV with `gene` and optional `category` columns",
    )
    parser.add_argument("--uniprot", action="append", default=[], help="Direct UniProt accession to query")
    parser.add_argument("--gene-uniprot-file", help="CSV with columns gene,uniprot")
    parser.add_argument("--forced-csv", help="CSV with columns gene,pdb_id")
    parser.add_argument("--trivial-ligands-file", help="Override/additional trivial ligand comp_ids")
    parser.add_argument("--keep-ligands-file", help="Allowlist comp_ids always treated non-trivial")
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--max-return", type=int)
    parser.add_argument("--resolution-max", type=float)
    parser.add_argument("--species", "--organism", default=None)
    parser.add_argument("--taxonomy-id", "--taxon", dest="taxonomy_id", type=int)
    parser.add_argument("--entity-type", default=None)
    parser.add_argument("--methods", default=None)
    parser.add_argument("--allow-apo", action="store_true")
    parser.add_argument(
        "--ligand-filter-mode",
        default="strict",
        choices=("strict", "relaxed"),
        help="How to treat borderline non-polymer ligands.",
    )
    parser.add_argument(
        "--ligand",
        "--contains-ligand",
        action="append",
        default=[],
        dest="ligands",
        help="Require a co-crystallized PDB chemical component ID, e.g. ATP or HEM.",
    )
    parser.add_argument(
        "--dedupe",
        nargs=2,
        metavar=("MODE", "CUTOFF"),
        help="De-duplicate search hits, e.g. `--dedupe sequence-identity 90`.",
    )
    parser.add_argument(
        "--dedupe-sequence-identity",
        type=int,
        help="Expert alias for `--dedupe sequence-identity <cutoff>`.",
    )
    parser.add_argument(
        "--no-dedupe",
        action="store_true",
        help="Disable default sequence-identity de-duplication.",
    )
    parser.add_argument("--quality", default="any", help="Target quality preset.")
    forced_group = parser.add_mutually_exclusive_group()
    forced_group.add_argument("--include-forced", dest="include_forced", action="store_true")
    forced_group.add_argument("--exclude-forced", dest="include_forced", action="store_false")
    parser.set_defaults(include_forced=None)
    install_group = parser.add_mutually_exclusive_group()
    install_group.add_argument("--dry-run", action="store_true", help="Select IDs but do not download")
    install_group.add_argument("--install", action="store_true", help="Install selected IDs without prompt")
    parser.add_argument("--max-install", type=int, help="Optional cap on number of downloads")
    parser.add_argument("--out-dir", help="Install destination folder")
    parser.add_argument("--selected-out", help="Write selected PDB list/metadata CSV")
    parser.add_argument("--missing-out", help="Write failed download statuses (CSV)")
    parser.add_argument("--category-limits-json", help="Override category limits JSON")
    parser.add_argument("--cache-dir")
    parser.add_argument("--cache-ttl-days", type=float)
    parser.add_argument("--no-cache", action="store_true", default=False)
    parser.add_argument("--rate-limit", type=float)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--max-retries", type=int)
    return parser


def _select_from_pdb_query(
    *,
    cfg: Mapping[str, Any],
    args: MassinstallArgs,
) -> tuple[list[dict[str, str]], list[str]]:
    genes_file = _resolve_genes_file(cfg, args)

    max_candidates = args.max_candidates if args.max_candidates is not None else 200
    max_return = args.max_return if args.max_return is not None else 10
    resolution_max = (
        args.resolution_max
        if args.resolution_max is not None
        else DEFAULT_RESOLUTION_MAX
    )
    args_species = getattr(args, "species", None)
    args_taxonomy_id = getattr(args, "taxonomy_id", None)
    args_entity_type = getattr(args, "entity_type", None)
    args_methods = getattr(args, "methods", None)
    args_allow_apo = getattr(args, "allow_apo", False)
    ligand_filter_mode = getattr(args, "ligand_filter_mode", "strict")
    species, taxonomy_id = resolve_species_filter(args_species, args_taxonomy_id)
    args.species = species
    args.taxonomy_id = taxonomy_id
    entity_type = args_entity_type if args_entity_type is not None else DEFAULT_ENTITY_TYPE
    methods = args_methods if args_methods is not None else DEFAULT_METHOD_FILTER
    require_ligand = not bool(args_allow_apo)
    ligand_comp_ids = normalize_comp_ids(getattr(args, "ligands", []) or [])
    dedupe_sequence_identity = parse_dedupe_option(
        getattr(args, "dedupe", None),
        getattr(args, "dedupe_sequence_identity", DEFAULT_DEDUPE_SEQUENCE_IDENTITY),
        no_dedupe=bool(getattr(args, "no_dedupe", False)),
    )
    quality = normalize_quality_preset(getattr(args, "quality", "any"))
    include_forced = args.include_forced if args.include_forced is not None else True
    cache_dir = Path(args.cache_dir).expanduser() if args.cache_dir else DEFAULT_CACHE_DIR
    cache_ttl_days = (
        args.cache_ttl_days
        if args.cache_ttl_days is not None
        else DEFAULT_CACHE_TTL_SECONDS / (24 * 60 * 60)
    )
    no_cache = args.no_cache
    rate_limit = args.rate_limit if args.rate_limit is not None else DEFAULT_RATE_LIMIT_RPS
    timeout = args.timeout if args.timeout is not None else DEFAULT_TIMEOUT_SECONDS
    max_retries = args.max_retries if args.max_retries is not None else DEFAULT_MAX_RETRIES

    category_limits_raw = args.category_limits_json if args.category_limits_json is not None else ""
    category_limits = _parse_category_limits(category_limits_raw)
    gene_rows = _read_gene_category_rows(str(genes_file))

    gene_to_uniprots = _read_gene_to_uniprots(args.gene_uniprot_file)
    for accession in _upper_list(getattr(args, "uniprot", []) or []):
        gene_to_uniprots.setdefault(accession, [])
        if accession not in gene_to_uniprots[accession]:
            gene_to_uniprots[accession].append(accession)
    forced_map = _merge_forced_maps(
        FORCED_PDBS_BY_GENE,
        _read_forced_map(args.forced_csv),
    )
    trivial_ligands = _read_comp_id_file(args.trivial_ligands_file)
    keep_ligands = _read_comp_id_file(args.keep_ligands_file)

    client = PDBQueryClient(
        cache_dir=cache_dir,
        use_cache=not no_cache,
        cache_ttl_seconds=max(1, int(cache_ttl_days * 24 * 60 * 60)),
        rate_limit_rps=rate_limit,
        timeout_seconds=timeout,
        max_retries=max_retries,
    )

    category_order: list[str] = []
    bucket: dict[str, list[tuple[str, RankedEntry]]] = defaultdict(list)
    for gene, category in gene_rows:
        category_u = category.upper() or DEFAULT_CATEGORY
        if category_u not in category_order:
            category_order.append(category_u)
        selection = get_ranked_entries_for_gene(
            gene,
            max_candidates=max_candidates,
            max_return=max_return,
            resolution_max=resolution_max,
            species_name=species,
            taxonomy_id=taxonomy_id,
            entity_type=entity_type,
            experimental_methods=methods,
            require_ligand=require_ligand,
            ligand_comp_ids=ligand_comp_ids,
            dedupe_sequence_identity=dedupe_sequence_identity,
            quality=quality,
            ligand_filter_mode=ligand_filter_mode,
            gene_to_uniprots=gene_to_uniprots,
            forced_pdbs_by_gene=forced_map,
            include_forced=include_forced,
            client=client,
            trivial_ligands=trivial_ligands,
            keep_ligands=keep_ligands,
        )
        for entry in selection.entries:
            bucket[category_u].append((gene, entry))

    selected_rows: list[dict[str, str]] = []
    selected_ids: list[str] = []
    seen: set[str] = set()
    for category in category_order:
        candidates = sorted(
            bucket.get(category, []),
            key=lambda item: (-item[1].rank_score, item[1].pdb_id, item[0]),
        )
        limit = category_limits.get(category)
        used = 0
        for gene, entry in candidates:
            if entry.pdb_id in seen:
                continue
            if limit is not None and used >= limit:
                break
            seen.add(entry.pdb_id)
            used += 1
            selected_ids.append(entry.pdb_id)
            selected_rows.append(
                {
                    "category": category,
                    "gene": gene,
                    "pdb_id": entry.pdb_id,
                    "rank_score": f"{entry.rank_score:.3f}",
                    "resolution": ""
                    if entry.resolution is None
                    else f"{entry.resolution:.3f}",
                    "method": entry.method,
                    "forced_included": str(entry.forced_included).lower(),
                    "valid_entry": str(entry.valid_entry).lower(),
                    "nontrivial_comp_ids": ";".join(entry.nontrivial_comp_ids),
                }
            )
    return selected_rows, selected_ids


def _resolve_output_dir(cfg: Mapping[str, Any], args: MassinstallArgs) -> Path:
    if args.out_dir:
        return Path(args.out_dir)
    cfg_out_dir = str(cfg.get("PDB_QUERY_DEFAULT_OUT_DIR", "")).strip()
    if cfg_out_dir:
        return Path(cfg_out_dir)
    return Path(cfg["INPUT_DIR"])


def _maybe_prompt_install(selected_count: int) -> bool:
    if not sys.stdin.isatty() or not sys.stdout.isatty():
        return False
    answer = input(
        f"Dry-run selected {selected_count} PDB IDs. Install now? [y/N]: "
    ).strip().lower()
    return answer in {"y", "yes"}


def _write_missing_csv(path: str, rows: Sequence[tuple[str, int]]) -> None:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_id", "status_code"])
        writer.writeheader()
        for pdb_id, status in rows:
            writer.writerow({"pdb_id": pdb_id, "status_code": str(status)})


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = MassinstallArgs.from_namespace(parser.parse_args(argv))

    cfg = load_inputs()
    validate_config(cfg)
    output_dir = _resolve_output_dir(cfg, args)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.pdb_query:
        selected_rows, pdb_ids = _select_from_pdb_query(cfg=cfg, args=args)
        print(f"Selected {len(pdb_ids)} unique PDB IDs via --pdb_query")
        counts_by_category: dict[str, int] = defaultdict(int)
        for row in selected_rows:
            counts_by_category[row["category"]] += 1
        for category in sorted(counts_by_category):
            print(f"  {category}: {counts_by_category[category]}")

        selected_out = _pick_str(args.selected_out, cfg, "PDB_QUERY_SELECTED_OUT")
        if selected_out:
            _write_selection_csv(selected_out, selected_rows)

        should_install = False
        if args.install:
            should_install = True
        elif not args.dry_run:
            should_install = _maybe_prompt_install(len(pdb_ids))
        if not should_install:
            print("Dry-run complete. No downloads performed.")
            return 0

        max_install = args.max_install
        if max_install is None:
            max_install_raw = _cfg_str(cfg, "PDB_QUERY_MAX_INSTALL")
            if max_install_raw:
                max_install = int(max_install_raw)
        success, failures = _download_pdbs(
            pdb_ids, output_dir=output_dir, max_install=max_install
        )
        missing_out = _pick_str(args.missing_out, cfg, "PDB_QUERY_MISSING_OUT")
        if missing_out:
            _write_missing_csv(missing_out, failures)
        print(
            f"Download summary: success={len(success)} failures={len(failures)} output_dir={output_dir}"
        )
        return 0

    # Legacy mode: keep behavior of using PDB_ID_LIST/fallback and downloading.
    pdb_ids = _load_legacy_pdb_ids(cfg)
    max_install = args.max_install
    if max_install is None:
        legacy_max_raw = _cfg_str(cfg, "MASSINSTALL_MAX_INSTALL")
        if legacy_max_raw:
            max_install = int(legacy_max_raw)
    success, failures = _download_pdbs(
        pdb_ids, output_dir=output_dir, max_install=max_install
    )
    missing_out = _pick_str(args.missing_out, cfg, "MASSINSTALL_MISSING_OUT")
    if missing_out:
        _write_missing_csv(missing_out, failures)
    print(f"Download summary: success={len(success)} failures={len(failures)} output_dir={output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
