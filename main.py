# -*- coding: utf-8 -*-
# High-level pipeline for multi-stage docking.
#
# Phases per protein:
# 1) Setup & logging
# 2) Extract ligands and generate a ligand-free PDB
# 3) Protein preparation (clean PDB + receptor PDBQT), reuse if cached
# 4) Active-site detection (center, box size)
# 5) Ligand preparation & filtering
# 6) Multi-stage docking with early/fallback recenter heuristics + CenterSelector
# 7) Final pose validation & optional screenshots
# 8) Write per-protein score CSV

from __future__ import annotations

import sys, hashlib, re, logging, json, time, os, shutil, re
from dataclasses import dataclass
import datetime
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, Any, Mapping
import numpy as np
from tqdm import tqdm
import traceback
from logging_topics import (
    _tee_stdio_to,
    make_protein_logger,
    bootstrap_root_logging,
)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
from input_and_export_functions import (
    load_inputs, validate_config, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, emit_vina_config as _emit_vina_config_impl, record_score, score_key, _to_bool, init_config_run_dir
)
from protein_functions import detect_active_site
from activesite import extract_and_remove_ligands, get_atom_rules
from prep_ligands import prep_ligands_from_pdb, is_valid_ligand, enumerate_ligands_for_docking
from pose_validation import (
    validate_pose_pdbqt, extract_surface_atoms, attempt_fallback_recenter,
    filter_and_rewrite_poses_by_rmsd, compute_self_rmsd
)
from run_vina import run_docking_task, validate_all_poses
from path_router import (
    expand_variants,
    make_paths,
    Paths as RouterPaths,
    receptor_file,
    docked_dir,
    config_dir as router_config_dir,
)
from multi_stage_docking import run_one_stage, RetryManager
from ph_ensemble_docking import (
    enumerate_ligands_for_ph_context,
    init_ph_tags_and_manifest,
    prewarm_ph_ligand_microstates,
)
from fallback_recenter import (
    BudgetGuard,
    RecenterParams,
    GlobalCenterGuard,
    validate_first_valid_pose,
    fallback_recentering_if_empty,
)
from apo_holo_mode import (
    resolve_apo_holo_mode,
    _debug_normalize_mode_token,
    _variant_receptor_path,
    file_sha1,
    delete_variant_trees,
    _record_apo_holo_usage,
    _record_apo_holo_decision,
)
from library_index import LibraryIndex
import automate_protein_prep as protein_prep

# [ions] audit classification tokens
_ION_AUDIT_METALS = {
    "ZN",
    "MG",
    "MN",
    "FE",
    "CU",
    "CO",
    "NI",
    "CA",
    "HG",
}
_ION_AUDIT_SALTS = {
    "NA",
    "K",
    "CL",
    "BR",
    "I",
    "LI",
    "RB",
    "CS",
}


def _resolve_run_id(argv: list[str]) -> str:
    cli_run_id = _cli_val(argv, "--run-id")
    env_run_id = (os.environ.get("ATLAS_RUN_ID") or "").strip()
    if cli_run_id:
        return cli_run_id
    if env_run_id:
        return env_run_id
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def _prepare_run_logfile(run_id: str) -> str:
    logs_dir = Path.cwd() / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    return str(logs_dir / f"main_{run_id}.log")




class ConfigDict(dict):
    __slots__ = ()
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc
    def __setattr__(self, key, value):
        self[key] = value
    def copy(self):
        return ConfigDict(super().copy())


def _log_cfg_emit_path_check(pdb_id, receptor_path, variant, legacy):
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
    variant_label = variant_token or "legacy"
    contains_variant = bool(variant_token and variant_token in str(receptor_path))
    logging.info(
        "[cfg.emit.check] pdb=%s variant=%s receptor_path=%s path_contains_variant=%s",
        pdb_id,
        variant_label,
        receptor_path,
        contains_variant,
    )
    if variant_token and not contains_variant and not legacy:
        logging.warning(
            "[variant.mismatch] expected_variant=%s wrote_legacy_path=%s action=fail_ci",
            variant_token,
            receptor_path,
        )


def emit_vina_config(
    cfg,
    pdb_id,
    receptor_pdbqt,
    center,
    box_size,
    ligand_path,
    stage_name,
    stage_info,
    cpu_per_job,
    logger=None,
    *,
    variant=None,
    ph_token=None,
    legacy=False,
    skip_manifest_if_exists=False,
):
    if not skip_manifest_if_exists:
        result = _emit_vina_config_impl(
            cfg,
            pdb_id,
            receptor_pdbqt,
            center,
            box_size,
            ligand_path,
            stage_name,
            stage_info,
            cpu_per_job,
            logger,
            variant=variant,
            ph_token=ph_token,
            legacy=legacy,
        )
    else:
        result = _emit_vina_config_skip_manifest(
            cfg,
            pdb_id,
            receptor_pdbqt,
            center,
            box_size,
            ligand_path,
            stage_name,
            stage_info,
            cpu_per_job,
            logger,
            variant=variant,
            ph_token=ph_token,
            legacy=legacy,
        )
    _log_cfg_emit_path_check(pdb_id, receptor_pdbqt, variant, legacy)
    return result


def _emit_vina_config_skip_manifest(
    cfg,
    pdb_id,
    receptor_pdbqt,
    center,
    box_size,
    ligand_path,
    stage_name,
    stage_info,
    cpu_per_job,
    logger,
    *,
    variant=None,
    ph_token=None,
    legacy=False,
):
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
    ph_label = (str(ph_token).strip() or None) if ph_token is not None else None
    legacy_mode = bool(legacy)

    lig_base = Path(ligand_path).stem
    run_id = cfg["RUN_ID"]

    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    cfg_dir.mkdir(parents=True, exist_ok=True)

    cfg_path = cfg_dir / f"{lig_base}_{stage_name}.txt"

    stage_root = docked_dir(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    stage_root.mkdir(parents=True, exist_ok=True)
    out_dir = stage_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{lig_base}_{stage_name}.pdbqt"

    expected_receptor = receptor_file(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    receptor_exists = expected_receptor.exists()
    receptor_for_config = str(expected_receptor)

    variant_display = variant_token or "None"
    ph_display = ph_label or "None"
    breadcrumb = (
        "[cfg.emit] run=%s pdb=%s stage=%s variant=%s ph=%s\n"
        "           cfg_dir=%s receptor=%s out_root=%s"
    )
    breadcrumb_args = (
        run_id,
        pdb_id,
        stage_name,
        variant_display,
        ph_display,
        str(cfg_dir),
        receptor_for_config,
        str(stage_root),
    )
    if logger:
        logger.info(breadcrumb, *breadcrumb_args)
    else:
        print(breadcrumb % breadcrumb_args)

    if not receptor_exists:
        msg = (
            f"[router.error] missing receptor for pdb={pdb_id} variant={variant_display} "
            f"ph={ph_display} -> {expected_receptor}"
        )
        if logger:
            logger.error(msg)
        else:
            print(msg)

    lines = [
        f"receptor = {receptor_for_config}",
        f"ligand   = {ligand_path}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {box_size[0]:.3f}",
        f"size_y   = {box_size[1]:.3f}",
        f"size_z   = {box_size[2]:.3f}",
        f"cpu      = {int(cpu_per_job)}",
        f"exhaustiveness = {int(stage_info.get('exhaustiveness', 8))}",
        f"energy_range   = {int(stage_info.get('energy_range', 4))}",
        f"num_modes      = {int(stage_info.get('num_modes', 4))}",
        f"verbosity      = {int(stage_info.get('verbosity', 0))}",
        f"out = {out_path}",
    ]

    if logger:
        cx, cy, cz = center
        sx, sy, sz = box_size
        lig_name = os.path.basename(str(ligand_path))
        logger.info(
            "[vina.cfg] lig=%s center=(%.3f,%.3f,%.3f) size=(%.1f,%.1f,%.1f)",
            lig_name,
            cx,
            cy,
            cz,
            sx,
            sy,
            sz,
        )

    payload = ("\n".join(lines)).encode("utf-8")
    overwrite = cfg_path.exists()

    tmp = cfg_path.with_suffix(".part")
    with open(tmp, "wb") as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, cfg_path)

    emit_msg = (
        "[cfg.emit] run=%s pdb=%s variant=%s ph=%s stage=%s ligand=%s "
        "cfg_dir=%s docked_root=%s path=%s overwrite=%s bytes=%d"
    )
    emit_args = (
        run_id,
        pdb_id,
        variant_display,
        ph_display,
        stage_name,
        lig_base,
        str(cfg_dir),
        str(stage_root),
        str(cfg_path),
        str(overwrite).lower(),
        len(payload),
    )
    if logger:
        logger.info(emit_msg, *emit_args)
    else:
        print(emit_msg % emit_args)

    return str(cfg_path), str(out_path)

# --- Debug wrappers to locate legacy/incorrect folder creation ---
import re, traceback

_orig_mkdir = Path.mkdir
def _dbg_mkdir(self, *a, **k):
    path_str = str(self)
    # Log legacy cleaned_ligands creations
    if path_str.lower().endswith("_cleaned_ligands"):
        logging.error("[DBG] Path.mkdir for legacy path: %s\n%s",
                      path_str, "".join(traceback.format_stack(limit=6)))
    # Log unwanted processed_pdbs/<PDB>_CLEANED creations (case-insensitive)
    if "/processed_pdbs/" in path_str and re.search(r"(?i)_cleaned/?$", path_str):
        logging.error("[DBG] Path.mkdir for processed_pdbs CLEANED dir: %s\n%s",
                      path_str, "".join(traceback.format_stack(limit=6)))
    return _orig_mkdir(self, *a, **k)
Path.mkdir = _dbg_mkdir

_orig_makedirs = os.makedirs
def _dbg_makedirs(name, *a, **k):
    p = str(name)
    if p.lower().endswith("_cleaned_ligands"):
        logging.error("[DBG] os.makedirs for legacy path: %s\n%s",
                      p, "".join(traceback.format_stack(limit=6)))
    if "/processed_pdbs/" in p and re.search(r"(?i)_cleaned/?$", p):
        logging.error("[DBG] os.makedirs for processed_pdbs CLEANED dir: %s\n%s",
                      p, "".join(traceback.format_stack(limit=6)))
    return _orig_makedirs(name, *a, **k)
os.makedirs = _dbg_makedirs





import os, re, hashlib
from pathlib import Path
from typing import Iterable, Set

_SANITIZED_RUN = re.compile(r'(?:\.sanitized){2,}')

def _collapse_sanitized_token(fn: str) -> str:
    """Collapse any repeated '.sanitized' tokens anywhere in the stem."""
    stem, ext = os.path.splitext(fn)
    new_stem = _SANITIZED_RUN.sub('.sanitized', stem)
    return new_stem + ext

def _sha1(path: str, bufsize: int = 1 << 20) -> str:
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        while True:
            b = f.read(bufsize)
            if not b: break
            h.update(b)
    return h.hexdigest()

def _files_identical(a: str, b: str) -> bool:
    try:
        if os.path.getsize(a) != os.path.getsize(b):
            return False
        return _sha1(a) == _sha1(b)
    except Exception:
        return False

def collapse_sanitized_names(root_dirs: Iterable[str],
                             exts: Set[str] = {'.pdb', '.sdf', '.mol2', '.pdbqt'},
                             logger=None) -> None:
    """
    Walk given roots and collapse repeated '.sanitized' in filenames.
    If the canonical name exists:
      - if byte-identical, delete the redundant file
      - if different, keep the canonical (shortest run) and delete the longer-run file; warn.
    """
    for root in root_dirs:
        if not root or not os.path.isdir(root):
            continue
        for dirpath, _, files in os.walk(root):
            for fn in files:
                ext = os.path.splitext(fn)[1].lower()
                if ext not in exts:
                    continue
                new_fn = _collapse_sanitized_token(fn)
                if new_fn == fn:
                    continue
                src = os.path.join(dirpath, fn)
                dst = os.path.join(dirpath, new_fn)
                rel_src = os.path.relpath(src, root)
                rel_dst = os.path.relpath(dst, root)
                try:
                    if os.path.exists(dst):
                        if _files_identical(src, dst):
                            os.remove(src)
                            if logger:
                                logger.info(f"[sanitize-collapse] dedup: removed duplicate '{rel_src}' (kept '{rel_dst}')")
                        else:
                            # Prefer the shorter '.sanitized' run (i.e., dst). Remove the longer one.
                            os.remove(src)
                            if logger:
                                logger.warning(f"[sanitize-collapse] conflict: kept '{rel_dst}', removed longer-run '{rel_src}'")
                    else:
                        os.rename(src, dst)
                        if logger:
                            logger.info(f"[sanitize-collapse] rename: '{rel_src}' -> '{rel_dst}'")
                except Exception as e:
                    if logger:
                        logger.error(f"[sanitize-collapse] failed on '{rel_src}' -> '{rel_dst}': {e}")






# ======================
# ======================
# Data models & utilities
# ======================
# >>> PATHS CLASS START
Paths = RouterPaths
# >>> PATHS CLASS END


def norm(p: str | Path) -> str:
    """Normalize path to a clean, forward-slash string for logs & keys."""
    return os.path.abspath(str(p)).replace("\\", "/")

# --- Single-ligand ---
def _parse_single_from_cli(argv) -> str:
    """
    Minimal CLI parser for: --single <pattern>
    Returns the pattern string or "" if not provided.
    """
    try:
        if "--single" in argv:
            i = argv.index("--single")
            if i + 1 < len(argv) and not argv[i+1].startswith("-"):
                return argv[i+1]
    except Exception:
        pass
    return ""
def _parse_fast_flag(argv) -> bool:
    """Return True if argv includes fast/-fast/--fast (case-insensitive)."""
    try:
        return any(tok.lower().lstrip("-") == "fast" for tok in argv)
    except Exception:
        return False
def _iter_pdbqt_dirfirst(root: Path, allowed_subdirs: Optional[set[str]] = None):
    """
    Yield .pdbqt files with a directory-first strategy:
      - list files directly under `root`
      - then list files under first-level subdirs, optionally restricted by `allowed_subdirs`
    Falls back to rglob if listing fails (robustness over speed).
    """
    try:
        if not root or not root.exists():
            return
        # files at root
        for p in root.glob("*.pdbqt"):
            yield p
        # first-level subdirs (dir-name filter first, then files)
        for d in root.iterdir():
            if not d.is_dir():
                continue
            if allowed_subdirs is not None and d.name not in allowed_subdirs:
                continue
            for p in d.glob("*.pdbqt"):
                yield p
    except Exception:
        # robust fallback
        for p in root.rglob("*.pdbqt"):
            yield p


# [single-index] Deduplicate manifest roots and retain stable ordering.
def _dedupe_manifest_roots(seq) -> list[Path]:
    deduped: list[Path] = []
    seen: set[str] = set()
    for entry in seq:
        if not entry:
            continue
        path_obj = Path(entry)
        try:
            key = str(path_obj.resolve())
        except Exception:
            key = str(path_obj)
        if key not in seen:
            seen.add(key)
            deduped.append(path_obj)
    return deduped


# [single-index] Prime the manifest-backed index before single-ligand lookups.
def _ensure_single_ligand_index(cfg: Dict, paths: Paths, logger: logging.Logger) -> None:
    cfg.setdefault(
        "OUTPUT_LIGANDS_DIR",
        cfg.get("PREPPED_LIGANDS_DIR") or cfg.get("PREPPED_LIGANDS_ROOT"),
    )

    per_roots = _dedupe_manifest_roots([
        getattr(paths, "prepped_ligands_dir", None),
    ])
    library_roots = _dedupe_manifest_roots([
        cfg.get("OUTPUT_LIGANDS_DIR"),
        cfg.get("PREPPED_LIGANDS_ROOT"),
        cfg.get("PREPPED_LIGANDS_DIR"),
    ])

    cfg["_LIB_INDEX_PER_ROOTS"] = [str(p) for p in per_roots]
    cfg["_LIB_INDEX_LIBRARY_ROOTS"] = [str(p) for p in library_roots]

    if not per_roots and not library_roots:
        cfg["_LIB_INDEX"] = None
        return

    manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
    index = cfg.get("_LIB_INDEX")
    if not isinstance(index, LibraryIndex):
        index = LibraryIndex(manifest_filename=manifest_filename, logger=logger)
        cfg["_LIB_INDEX"] = index

    index_roots = _dedupe_manifest_roots([*per_roots, *library_roots])
    if index_roots:
        status_parts = []
        for root in index_roots:
            manifest_path = Path(root) / manifest_filename
            manifest_state = "present" if manifest_path.exists() else "missing"
            status_parts.append(f"{manifest_path}={manifest_state}")
        index.load(index_roots)
        logger.info(
            "[single.index] per_roots=%s lib_roots=%s manifests=%s",
            [str(p) for p in per_roots],
            [str(p) for p in library_roots],
            ";".join(status_parts) or "none",
        )
    else:
        logger.info(
            "[single.index] per_roots=%s lib_roots=%s manifests=none",
            [str(p) for p in per_roots],
            [str(p) for p in library_roots],
        )


def _resolve_single_ligand(selector: str, pdb_id: str, cfg: Dict, logger: logging.Logger) -> Optional[Path]:
    """Resolve SINGLE_LIGAND selector using manifest-backed lookups (if enabled)."""
    if not selector:
        return None

    allow_prefix = _to_bool(str(cfg.get("SINGLE_LIGAND_ALLOW_PREFIX", "false")))
    raw_order = str(cfg.get("SINGLE_LIGAND_SEARCH_ORDER", "fda_library,per_protein,global"))
    order = ["fda_library", "per_protein", "global"]
    skip_global = bool(cfg.get("SINGLE_LIGAND_SKIP_GLOBAL", True))
    if skip_global:
        order = [scope for scope in order if scope != "global"]
    manifest_only = _to_bool(str(cfg.get("SINGLE_LIGAND_MANIFEST_ONLY", True)))
    suggestions_cap = int(cfg.get("SINGLE_LIGAND_SUGGESTIONS", 5) or 5)

    per_roots = _dedupe_manifest_roots(cfg.get("_LIB_INDEX_PER_ROOTS", []))
    library_roots = _dedupe_manifest_roots(cfg.get("_LIB_INDEX_LIBRARY_ROOTS", []))

    lib_index = cfg.get("_LIB_INDEX")
    has_index = isinstance(lib_index, LibraryIndex)

    output_ligands_dir = cfg.get("OUTPUT_LIGANDS_DIR")
    fda_root = Path(output_ligands_dir).joinpath("fda_library") if output_ligands_dir else None
    name_map: Optional[dict[str, set[str]]] = None
    fda_logged = False
    key = _norm_name_key(selector)

    effective_order = order
    use_manifest = has_index and (bool(per_roots) or (not skip_global and bool(library_roots)))

    logger.info(
        "[single.debug] selector=%s raw_order=%s order=%s skip_global=%s manifest_only=%s has_index=%s use_manifest=%s per_roots=%s lib_roots=%s allow_prefix=%s",
        selector,
        raw_order,
        effective_order,
        skip_global,
        manifest_only,
        str(has_index).lower(),
        str(use_manifest).lower(),
        [str(p) for p in per_roots],
        [str(p) for p in library_roots],
        str(allow_prefix).lower(),
    )

    def _ensure_name_map() -> dict[str, set[str]]:
        nonlocal name_map
        if name_map is None:
            name_map = _load_fda_name_map(cfg, logger)
        return name_map

    def _resolve_fda_scope() -> Optional[Path]:
        nonlocal fda_logged
        if fda_root is None or not fda_root.exists():
            return None
        mapping = _ensure_name_map()
        basenames: list[str] = []
        if key:
            if key in mapping:
                basenames.extend(sorted(mapping.get(key, set())))
            if not basenames and allow_prefix:
                for map_key, values in mapping.items():
                    if map_key.startswith(key):
                        basenames.extend(sorted(values))
        basenames = list(dict.fromkeys(basenames))
        probe_target: Path = fda_root if not basenames else fda_root / basenames[0]
        if not fda_logged:
            logger.info("[single.name] key=%s basenames=%s probe=%s", key, basenames, norm(probe_target))
            fda_logged = True
        for basename in basenames:
            candidate = fda_root / basename
            if candidate.exists():
                logger.info(
                    "[single.fda.hit] key=%s basename=%s path=%s",
                    key,
                    basename,
                    norm(candidate),
                )
                logger.info("[single.lookup.hit] source=fda selector=%s path=%s", selector, norm(candidate))
                return candidate
        return None

    def _suggest_from_index() -> list[str]:
        if isinstance(lib_index, LibraryIndex):
            ordered: list[Path] = []
            seen: set[str] = set()
            for candidate in [*per_roots, *library_roots]:
                try:
                    key_str = str(candidate.resolve())
                except Exception:
                    key_str = str(candidate)
                if key_str not in seen:
                    seen.add(key_str)
                    ordered.append(candidate)
            return lib_index.suggest(selector, ordered, suggestions_cap)
        return []

    def _log_miss(suggestions: list[str]) -> None:
        logger.error(
            "[single.lookup.miss] selector=%s suggestions=[%s]",
            selector,
            ",".join(suggestions),
        )

    hit = _resolve_fda_scope()
    if hit:
        return hit

    if manifest_only and not use_manifest:
        _log_miss(_suggest_from_index())
        return None

    if not has_index:
        _log_miss(_suggest_from_index())
        return None

    if not use_manifest:
        _log_miss(_suggest_from_index())
        return None

    if per_roots:
        hit = lib_index.lookup(selector, per_roots, allow_prefix=allow_prefix)
        if hit:
            logger.info("[single.lookup.hit] source=manifest scope=per_protein selector=%s path=%s", selector, norm(hit))
            return hit

    if not skip_global and library_roots:
        hit = lib_index.lookup(selector, library_roots, allow_prefix=allow_prefix)
        if hit:
            logger.info("[single.lookup.hit] source=manifest scope=global selector=%s path=%s", selector, norm(hit))
            return hit

    _log_miss(_suggest_from_index())
    return None
# --- FDA name mapping (CSV) ---------------------------------------------------
# Lets SINGLE_LIGAND resolve by generic/brand/synonym (e.g., "imatinib", "Gleevec").
_FDA_NAME_MAP_CACHE = None

def _norm_name_key(s: str) -> str:
    import re
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())

def _split_multi_names(v: str) -> list[str]:
    import re
    parts = re.split(r"[|;,/]", v or "")
    return [p.strip() for p in parts if p and p.strip()]

def _load_fda_name_map(cfg: Dict, logger: logging.Logger) -> dict[str, set[str]]:
    """
    Build dict: normalized_name -> {pdbqt_basename, ...}
    CSV must have at least: column 'path' pointing to a *.pdbqt, plus name columns.
    """
    global _FDA_NAME_MAP_CACHE
    if isinstance(_FDA_NAME_MAP_CACHE, dict):
        return _FDA_NAME_MAP_CACHE

    import csv
    from pathlib import Path

    csv_path = os.environ.get("FDA_MAPPING_CSV", "").strip() or str(cfg.get("FDA_MAPPING_CSV", "")).strip()
    mapping: dict[str, set[str]] = {}
    if not csv_path:
        _FDA_NAME_MAP_CACHE = {}
        return _FDA_NAME_MAP_CACHE

    p = Path(csv_path)
    if not p.exists():
        logger.info(f"[single:name] FDA_MAPPING_CSV not found at {csv_path} (name lookup disabled).")
        _FDA_NAME_MAP_CACHE = {}
        return _FDA_NAME_MAP_CACHE

    try:
        with open(p, newline="", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f)
            for row in reader:
                path = (row.get("path") or "").strip()
                if not path.endswith(".pdbqt"):
                    continue
                base = os.path.basename(path)

                # "single" name fields
                singles = [
                    row.get("display_name", ""),
                    row.get("generic_name", ""),
                    row.get("rxnorm_generic_name", ""),
                    row.get("drugcentral_generic_name", ""),
                    row.get("pubchem_name", ""),
                    row.get("pubchem_record_title", ""),
                ]
                # multi-value fields (split)
                multis = []
                for col in ("brand_names", "rxnorm_brand_names", "drugcentral_brand_names", "pubchem_synonyms"):
                    v = row.get(col, "")
                    if v:
                        multis.extend(_split_multi_names(v))

                for nm in [*singles, *multis]:
                    key = _norm_name_key(nm)
                    if key:
                        mapping.setdefault(key, set()).add(base)

        logger.info(f"[single:name] Loaded FDA name map ({len(mapping)} keys) from {p}")
    except Exception as e:
        logger.warning(f"[single:name] Failed to load name map: {e}")
        mapping = {}

    _FDA_NAME_MAP_CACHE = mapping
    return mapping

def _cli_val(argv, flag):
    try:
        if flag in argv:
            i = argv.index(flag)
            if i + 1 < len(argv) and not argv[i+1].startswith("-"):
                return argv[i+1]
    except Exception:
        pass
    return None

def _cli_has(argv, flag):
    try:
        return flag in argv
    except Exception:
        return False

# --- Specified Proteins Mode helpers (NEW) -----------------------------------
from typing import Iterable

def _norm_pdb_id(token: str) -> Optional[str]:
    """
    Normalize a user token to a 4-char PDB ID (uppercase).
    Accepts bare IDs (2HYY), QoL flags (--2HYY or -2HYY), and filenames (2HYY.pdb).
    Returns None if it cannot produce a 4-char alnum ID.
    """
    if not token:
        return None
    t = str(token).strip()
    # Strip any leading dashes (one or two)
    while t.startswith("-"):
        t = t[1:]
    t = os.path.basename(t)
    if t.lower().endswith(".pdb"):
        t = t[:-4]
    t = t.replace("_cleaned", "")
    t = t.upper()
    if len(t) >= 4:
        cand = t[:4]
        return cand if cand.isalnum() else None
    return None


def _split_ids(s: str) -> list[str]:
    """Split a comma/whitespace separated string into normalized 4-char IDs."""
    if not s:
        return []
    parts = s.replace(",", " ").split()
    out = []
    for p in parts:
        nid = _norm_pdb_id(p)
        if nid:
            out.append(nid)
    return out

def _dedupe_order(seq: Iterable[str]) -> list[str]:
    """De-duplicate while preserving first-seen order."""
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

def _parse_specified_proteins(argv, cfg) -> tuple[list[str], str]:
    """
    Resolve requested PDB IDs with precedence CLI > ENV > CFG.
    CLI:
      --pdb 2HYY        (repeatable)
      --pdbs 2HYY,3ERT  (comma/space separated)
      --2HYY            (QoL: any --<4char> alnum)
    ENV: ONLY_PDBS="2HYY 3ERT"
    CFG: SPECIFIED_PROTEINS: JSON list or string "2HYY, 3ERT"
    Returns: (normalized_ids, source or "")
    """
    # --- CLI ---
    cli_ids: list[str] = []

    # --pdb (repeatable)
    i = 0
    while i < len(argv):
        if argv[i] == "--pdb" and i + 1 < len(argv) and not argv[i + 1].startswith("-"):
            nid = _norm_pdb_id(argv[i + 1])
            if nid:
                cli_ids.append(nid)
            i += 2
            continue
        i += 1

    # --pdbs "A B,C"
    try:
        if "--pdbs" in argv:
            j = argv.index("--pdbs")
            if j + 1 < len(argv) and not argv[j + 1].startswith("-"):
                cli_ids.extend(_split_ids(argv[j + 1]))
    except Exception:
        pass

    # QoL: --2HYY / -2HYY style (exact length, starts with '-' or '--', next 4 alnum)
    for tok in argv:
        low = tok.lower()
        # don't treat fast/-fast/--fast as a PDB short-form token
        if low in ("fast", "-fast", "--fast"):
            continue
        if (tok.startswith("--") and len(tok) == 6) or (tok.startswith("-") and len(tok) == 5):
            nid = _norm_pdb_id(tok)
            if nid:
                cli_ids.append(nid)



    if cli_ids:
        return _dedupe_order(cli_ids), "CLI"

    env_val = os.environ.get("ONLY_PDBS", "").strip()
    if env_val:
        return _dedupe_order(_split_ids(env_val)), "ENV"

    # --- CFG ---
    val = cfg.get("SPECIFIED_PROTEINS", "")
    if isinstance(val, list):
        cfg_ids = [_norm_pdb_id(x) for x in val]
        cfg_ids = [x for x in cfg_ids if x]
        if cfg_ids:
            return _dedupe_order(cfg_ids), "CFG"
        return [], ""
    s = str(val or "").strip()
    if s:
        return _dedupe_order(_split_ids(s)), "CFG"

    return [], ""



def get_recenter_params(cfg: Dict) -> RecenterParams:
    """Load recenter parameters from config with safe defaults."""
    return RecenterParams(
        EARLY_RECENTER_RATIO=float(cfg.get("EARLY_RECENTER_RATIO", 0.70)),
        EARLY_RECENTER_MIN_EVAL=int(cfg.get("EARLY_RECENTER_MIN_EVAL", 10)),
        EARLY_RECENTER_FAR_A=float(cfg.get("EARLY_RECENTER_FAR_A", 15.0)),
        EARLY_RECENTER_MEDIAN_A=float(cfg.get("EARLY_RECENTER_MEDIAN_A", 10.0)),
        ALLOW_BOX_EXPAND=bool(cfg.get("ALLOW_BOX_EXPAND", True)),
        MAX_RECENTER_ATTEMPTS=int(cfg.get("MAX_RECENTER_ATTEMPTS", 1)),
    )





# >>> MAKE_PATHS SHIM START
# make_paths is imported from path_router above (legacy helper removed).
# >>> MAKE_PATHS SHIM END


# --------- control ligand lookup (prefer SDF > MOL2 > PDB; search ligands_raw + reference) ---------
def build_control_lookup(paths: Paths) -> dict:
    """
    Map base extracted-ligand stem -> crystal file path.
    Prefer a readable .sdf > .mol2 > .pdb, searching ligands_raw and optional reference folder.
    """
    prefs = [".sdf", ".mol2", ".pdb"]
    by_base: dict[str, dict[str, Path]] = {}

    search_dirs = [paths.ligand_output_dir, paths.ligand_output_dir.parent / "reference"]
    for root in search_dirs:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            ext = p.suffix.lower()
            if ext not in prefs:
                continue
            base = p.stem.split("_stage")[0]
            by_base.setdefault(base, {})
            by_base[base][ext] = p

    chosen: dict[str, Path] = {}
    for base, candidates in by_base.items():
        # try in preference order, but only accept if RDKit can read it
        picked = None
        for ext in prefs:
            p = candidates.get(ext)
            if p and _is_readable_ref(p):
                picked = p
                break
        if picked:
            chosen[base] = picked

    return chosen


# --------- extra helpers (easy-win features) ---------
def _map_reason_to_category(reason: str) -> str:
    if not reason:
        return "no_valid_pose"
    r = str(reason).lower()
    if "timeout" in r:
        return "timeout"
    if "too far" in r or "distance" in r or "centroid" in r:
        return "too_far_from_pocket"
    if "malformed" in r or "parse" in r or "format" in r:
        return "malformed"
    if "no pose" in r or "no_valid" in r or "all_poses_invalid" in r:
        return "no_valid_pose"
    return "no_valid_pose"


# --------- improved checkpointing (fingerprinted) ---------
def _file_md5(path: str, blocksize: int = 1 << 20) -> Optional[str]:
    try:
        h = hashlib.md5()
        with open(path, "rb") as f:
            while True:
                b = f.read(blocksize)
                if not b:
                    break
                h.update(b)
        return h.hexdigest()
    except Exception:
        return None


def _round_tuple(t: Tuple[float, float, float], ndp: int = 1) -> Tuple[float, float, float]:
    return tuple(None if (x is None) else round(float(x), ndp) for x in t)


def _fingerprint_stage(cfg: Dict,
                       receptor_pdbqt: str,
                       center: Tuple[float, float, float],
                       box_size: Tuple[float, float, float],
                       stage: Dict) -> Dict[str, Any]:
    rec_hash = _file_md5(receptor_pdbqt) if receptor_pdbqt else None
    stage_keys = ["name", "size", "exhaustiveness", "energy_range", "num_modes", "seed"]
    stage_core = {k: stage.get(k) for k in stage_keys if k in stage}
    return {
        "receptor_md5": rec_hash,
        "center": _round_tuple(center, 1),
        "box_size": _round_tuple(box_size, 1),
        "stage": stage_core,
        "vina_exe": str(cfg.get("VINA_EXE", "")),
        "threads_per_vina": int(cfg.get("THREADS_PER_VINA", 1)),
        "version_tag": "ckpt_v2",
    }


def _checkpoint_path(cfg: Dict, pdb_id: str, stage_name: str) -> Path:
    ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None
    variant = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    root = docked_dir(pdb_id, variant=variant, ph_tag=ph_label, legacy=legacy_mode)
    return root / f".ckpt_{stage_name}.json"


def checkpoint_should_skip(cfg: Dict,
                           pdb_id: str,
                           stage_name: str,
                           fingerprint: Dict[str, Any]) -> bool:
    p = _checkpoint_path(cfg, pdb_id, stage_name)
    if not p.exists():
        return False
    try:
        prev = json.loads(p.read_text())
    except Exception:
        return False
    return prev == fingerprint


def checkpoint_mark_done(cfg: Dict,
                         pdb_id: str,
                         stage_name: str,
                         fingerprint: Dict[str, Any]) -> None:
    p = _checkpoint_path(cfg, pdb_id, stage_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        p.write_text(json.dumps(fingerprint, indent=2))
    except Exception:
        pass


def checkpoint_invalidate_from(cfg: Dict, pdb_id: str, stages: List[Dict], start_index: int) -> None:
    for j in range(start_index, len(stages)):
        try:
            _checkpoint_path(cfg, pdb_id, stages[j]["name"]).unlink(missing_ok=True)
        except Exception:
            pass


def _write_audit_json(cfg: Dict, pdb_id: str, summary: Dict):
    try:
        if not cfg.get("AUDIT_JSON", True):
            return
        ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None
        variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
        legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
        out = docked_dir(pdb_id, variant=variant_env, ph_tag=ph_label, legacy=legacy_mode) / "audit.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, indent=2))
    except Exception:
        pass


def receptor_sanity_check(receptor_pdbqt: str, min_atoms: int = 10) -> bool:
    try:
        atoms = 0
        any_nonzero = False
        with open(receptor_pdbqt, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if (ln.startswith("ATOM") or ln.startswith("HETATM")) and len(ln) >= 54:
                    atoms += 1
                    try:
                        x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                        if (abs(x) + abs(y) + abs(z)) > 0.0:
                            any_nonzero = True
                    except Exception:
                        continue
        return (atoms >= min_atoms) and any_nonzero
    except Exception:
        return False


# ======================
# Phase 1 5: Prep steps
# ======================
def extract_ligands_to_nolig(paths: Paths, logger: logging.Logger) -> Tuple[int, set]:
    """
    Run crystallographic ligand extraction into processed_pdbs/<PDB>/ligands_raw/,
    create/update nolig/<PDB>_nolig.pdb, and return (n_controls, control_stems).
    Compatible with multiple legacy signatures of activesite.extract_and_remove_ligands().
    """
    from activesite import extract_and_remove_ligands  # authoritative extractor

    # Ensure intermediate dirs exist and clear the malformed log for a fresh run
    paths.ligand_output_dir.mkdir(parents=True, exist_ok=True)
    paths.ligands_mol2_dir.mkdir(parents=True, exist_ok=True)
    malformed_log = paths.ligands_mol2_dir / "malformed_ligands.txt"
    if malformed_log.exists():
        malformed_log.unlink()

    src_pdb = paths.input_pdb_path
    nolig_dst = paths.nolig_pdb_path
    ligands_dir = paths.ligand_output_dir

    logger.debug("[extract.debug] in=%s nolig=%s ldir=%s", src_pdb, nolig_dst, ligands_dir)

    ligands_dict, _ = extract_and_remove_ligands(
        str(src_pdb), str(nolig_dst), str(ligands_dir)
    )
    logger.info(f"Extracted {len(ligands_dict)} ligands -> {ligands_dir}")
    try:
        counts = {".pdb": 0, ".mol2": 0, ".sdf": 0}
        samples = []
        for ext in (".pdb", ".mol2", ".sdf"):
            for p in paths.ligand_output_dir.glob(f"*{ext}"):
                counts[ext] += 1
                if len(samples) < 6:
                    samples.append(p.name)
        logger.info("[extract.audit] counts=%s samples=%s", counts, samples)
    except Exception as e:
        logger.debug("[extract.audit] listing failed: %s", e)

    # Build control stems from what was actually written
    control_stems: set[str] = set()
    for ext in (".mol2", ".pdb", ".sdf"):
        for p in paths.ligand_output_dir.rglob(f"*{ext}"):
            control_stems.add(Path(p).stem)

    logger.info(f"Extracted {len(control_stems)} ligands -> {paths.ligand_output_dir}")
    return len(control_stems), control_stems




def _ph_control_centroid(paths: Paths) -> Optional[Tuple[float, float, float]]:
    try:
        primary = paths.ligand_output_dir
    except Exception:
        primary = None
    roots = []
    if primary is not None:
        roots.append(primary)
        legacy = primary.parent.parent / f"{paths.pdb_id}_NOLIG" / 'ligands_raw'
        roots.append(legacy)
    else:
        roots.append(paths.root_pdb_dir / 'ligands_raw')
        roots.append(paths.root_pdb_dir.parent / f"{paths.pdb_id}_NOLIG" / 'ligands_raw')
    centroids = []
    seen = set()
    for root in roots:
        if not root:
            continue
        root = Path(root)
        key = str(root)
        if key in seen or not root.exists():
            continue
        seen.add(key)
        for pdb_path in sorted(root.glob('*.pdb')):
            xs = ys = zs = count = 0.0
            try:
                with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as fh:
                    for line in fh:
                        if not line.startswith(('ATOM  ', 'HETATM')):
                            continue
                        try:
                            xs += float(line[30:38])
                            ys += float(line[38:46])
                            zs += float(line[46:54])
                            count += 1.0
                        except Exception:
                            continue
            except Exception:
                continue
            if count:
                centroids.append((xs / count, ys / count, zs / count))
    if centroids:
        n = float(len(centroids))
        return (sum(x for x, _, _ in centroids) / n,
                sum(y for _, y, _ in centroids) / n,
                sum(z for _, _, z in centroids) / n)
    return None


# pH helpers  ----------------------
def _resolve_ph_scope(scope_cfg: str, radius_nominal: float, paths: Paths, cleaned_pdb: str, log: logging.Logger) -> Tuple[str, Tuple[float, float, float], float]:
    scope = (scope_cfg or '').strip().lower()
    try:
        radius = float(radius_nominal)
    except Exception:
        radius = 10.0
    if radius <= 0:
        radius = 10.0
    if scope == 'pocket':
        center = _ph_control_centroid(paths)
        if center is None:
            try:
                detect_res = detect_active_site(cleaned_pdb)
            except Exception as exc:
                log.debug('[ph_ensemble.scope] detect_active_site failed: %s', exc)
                detect_res = None
            if detect_res and detect_res[0]:
                center = tuple(float(x) for x in detect_res[0])
        if center is None:
            log.warning('[ph_ensemble.scope] pocket requested but no center found; fallback=global')
            return 'global', (0.0, 0.0, 0.0), 1_000_000.0
        cx, cy, cz = (float(center[0]), float(center[1]), float(center[2]))
        return 'pocket', (cx, cy, cz), radius
    return 'global', (0.0, 0.0, 0.0), 1_000_000.0


def _ph_values_from_context(pdb_path: str) -> list[float]:
    """
    Ask context_ph for the full pH list (ensemble if present, else [target]),
    then round to 0.1 and clamp to [3.0, 10.5].
    """
    vals = []
    try:
        from context_ph import select_ph_values_for_protonation
        raw = select_ph_values_for_protonation(pdb_path)  # returns ensemble or [target]
        logging.info(f"[ph.ctx.list] taken_from_context={raw}")

        for x in (raw or []):
            # round & clamp
            v = max(3.0, min(10.5, round(float(x), 1)))
            vals.append(v)
        # dedupe + sort for stability
        vals = sorted({v for v in vals})
    except Exception as e:
        logging.warning(f"[ph.context] failed to resolve; falling back to [7.0]: {e}")
        vals = [7.0]
    logging.info(f"[ph.list] n={len(vals)} values={vals}")
    return vals


def _ph_ligand_mode(cfg: Mapping[str, Any]) -> str:
    """
    Interpret PH_LIGAND_MODE from config.

    Recognized values (case-insensitive):
    - off / none / false / 0 / "" -> "off"
    - anything else -> "context_window"

    This keeps the behavior opt-in while allowing future modes later.
    """
    raw = str(cfg.get("PH_LIGAND_MODE", "off")).strip().lower()
    if raw in ("", "off", "none", "false", "0"):
        return "off"
    return "context_window"


# ----------------------


# ----------------------



def prepare_receptor(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    center: Optional[Tuple[float, float, float]] = None,
    box_size: Optional[Tuple[float, float, float]] = None,
) -> Tuple[Optional[str], Optional[str]]:
    import automate_protein_prep
    from distutils.util import strtobool

    prepare_receptor.last_provenance = "unknown"
    force_reprocess = bool(strtobool(str(cfg.get("FORCE_REPROCESS", False))))
    log = logging.getLogger("ph_ensemble")
    # >>> RECEPTOR PATHS PATCH START
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant = var if var else None
    ph_token = None
    cleaned_pdb_path = paths.receptor_cleaned_pdb(variant)
    receptor_pdbqt_path = paths.receptor_pdbqt(variant, ph_token)
    # >>> RECEPTOR PATHS PATCH END
    logger.info(
        f"FORCE_REPROCESS={force_reprocess} | "
        f"cleaned_exists={cleaned_pdb_path.exists()} "
        f"receptor_exists={receptor_pdbqt_path.exists()}"
    )

    def _build_ph_ensemble(cleaned_path: str) -> Optional[str]:
        if not bool(cfg.get("PH_ENSEMBLE")):
            return None
        cleaned_path = str(cleaned_path)
        log.info("[ph_ensemble.anchor] cleaned_receptor_pdb=%s", cleaned_path)
        log.info("[ph_ensemble.begin] pdb_id=%s path=%s", paths.pdb_id, cleaned_path)

        def _collect_dock_targets(manifest_path: str) -> Optional[list[tuple[str, str]]]:
            manifest_file = Path(manifest_path)
            try:
                payload = json.loads(manifest_file.read_text())
            except Exception as exc:
                log.error("[ph_ensemble.manifest.read.error] path=%s err=%s", manifest_path, exc)
                return None

            members = payload.get("members") or []
            # Diagnostics: enumerate keys and canonical counts.
            try:
                key_universe = sorted({k for m in members for k in (m.keys() if isinstance(m, dict) else [])})
            except Exception:
                key_universe = []
            log.info(
                "[ph_ensemble.manifest.stats] members=%d canonical=%d keys=%s",
                len(members),
                sum(1 for m in members if isinstance(m, dict) and bool(m.get("canonical", False))),
                ",".join(key_universe),
            )

            canonical = [m for m in members if isinstance(m, dict) and bool(m.get("canonical", False))]
            if canonical:
                members = canonical

            prefix = f"{paths.pdb_id}_"
            targets: list[tuple[str, str]] = []
            for entry in members:
                if not isinstance(entry, dict):
                    continue
                receptor_path = (
                    entry.get("pdbqt")
                    or entry.get("receptor_pdbqt")
                    or entry.get("output_pdbqt")
                    or entry.get("path")
                    or entry.get("receptor")
                )
                if not receptor_path:
                    log.warning("[ph_ensemble.manifest.entry.missing_pdbqt] keys=%s", list(entry.keys()))
                    continue

                ph_label = entry.get("label") or entry.get("ph_label")
                if not ph_label:
                    stem = Path(receptor_path).stem
                    ph_label = stem[len(prefix):] if stem.startswith(prefix) else stem

                targets.append((str(ph_label), str(receptor_path)))

            if not targets:
                log.error("[ph_ensemble.manifest.no_targets] path=%s members=%d", manifest_path, len(members))
            return targets

        def _bridge_manifest_targets(manifest_path: str) -> None:
            targets = _collect_dock_targets(manifest_path)
            if targets is None:
                return
            cfg.setdefault("_PH_ENSEMBLE_CANONICAL", {})[paths.pdb_id] = targets
            # --- mapping debug (anchor: [ph_ensemble.map]) ---
            log.info(
                "[ph_ensemble.map] pdb_id=%s canonical=%s",
                paths.pdb_id,
                ";".join(f"{lbl}:{Path(p).name}" for lbl, p in targets) if targets else ""
            )
            variant_label = variant or "legacy"
            log.info(
                "[ph_ensemble.dock.begin] pdb_id=%s variant=%s n=%d labels=%s",
                paths.pdb_id,
                variant_label,
                len(targets),
                ",".join(lbl for lbl, _ in targets) or "",
            )

            if not targets:
                log.info(
                    "[ph_ensemble.dock.done] pdb_id=%s variant=%s n=0 targets=",
                    paths.pdb_id,
                    variant_label,
                )
                return
            legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
            for ph_label, _ in targets:
                ph_root = docked_dir(paths.pdb_id, variant=variant, ph_tag=ph_label, legacy=legacy_mode)
                log.info(
                    "[ph_ensemble.dock.root] pdb_id=%s variant=%s ph=%s dock_root=%s",
                    paths.pdb_id,
                    variant_label,
                    ph_label,
                    str(ph_root),
                )
            summary = ";".join(f"{ph}:{rec}" for ph, rec in targets)
            log.info(
                "[ph_ensemble.dock.done] pdb_id=%s variant=%s n=%d targets=%s",
                paths.pdb_id,
                variant_label,
                len(targets),
                summary,
            )

        try:
            from context_ph import select_ph_from_pdb
            ctx_result = select_ph_from_pdb(cleaned_path)
        except Exception as exc:
            log.warning("[ph_ensemble.ctx.error] %s", exc)
            ctx_result = {"target_pH": 7.0, "ensemble": None}
        target_pH = float(ctx_result.get("target_pH", 7.0) or 7.0)
        ensemble_from_context = ctx_result.get("ensemble")
        log.info("[ph_ensemble.ctx] target_pH=%.2f raw_ensemble=%s", target_pH, repr(ensemble_from_context))
        raw_values = list(ensemble_from_context or [target_pH])
        ph_values: list[float] = []
        for value in raw_values:
            try:
                ph = float(value)
            except Exception:
                continue
            ph = round(ph, 1)
            if ph < 3.0:
                ph = 3.0
            if ph > 10.5:
                ph = 10.5
            ph_values.append(ph)
        if not ph_values:
            fallback = round(target_pH, 1)
            if fallback < 3.0:
                fallback = 3.0
            if fallback > 10.5:
                fallback = 10.5
            ph_values = [fallback]
        ph_values = sorted({round(p, 1) for p in ph_values})
        log.info("[ph_ensemble.list] canonical=%s", ",".join(f"{p:.1f}" for p in ph_values))
        radius_nominal = getattr(cfg, "PH_RADIUS", 10.0)
        try:
            radius_nominal = float(radius_nominal)
        except Exception:
            radius_nominal = 10.0
        scope_cfg = getattr(cfg, "PH_SCOPE", "")
        scope, center, eff_radius = _resolve_ph_scope(scope_cfg, radius_nominal, paths, cleaned_path, log)
        log.info("[ph_ensemble.pick.scope] scope=%s", scope)
        if scope == "pocket":
            log.info("[ph_ensemble.pick.center] center=(%.3f,%.3f,%.3f) radius=%.1f", center[0], center[1], center[2], radius_nominal)
        else:
            log.info("[ph_ensemble.pick.center] center=GLOBAL radius=ALL")
        log.info("[ph_ensemble.call] building pH ensemble for %s", paths.pdb_id)
        prev_cfg = getattr(automate_protein_prep, "config", None)
        try:
            automate_protein_prep.config = cfg
        except Exception:
            prev_cfg = None
        try:
            import ph_ensemble
            manifest_path = ph_ensemble.build_ph_ensemble(
                pdb_id=paths.pdb_id,
                cleaned_receptor_pdb=cleaned_path,
                out_dir=str(Path(cfg["OUTPUT_DIR"])),
                center=center,
                radius=eff_radius,
                ph_values=ph_values,
                variant=variant,
                legacy=bool(cfg.get("_ROUTER_LEGACY", False)),
            )

            log.info("[ph_ensemble.manifest] path=%s", manifest_path)
            log.info("[ph.manifest.json] written=%s", manifest_path)

            # --- DEBUG: measure map size before/after bridge ---
            try:
                _pre = len((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).get(paths.pdb_id, []))
            except Exception:
                _pre = -1
            log.info("[ph_ensemble.debug] before-bridge map_len[%s]=%d", paths.pdb_id, _pre)

            if manifest_path:
                _bridge_manifest_targets(str(manifest_path))

            try:
                _post = len((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).get(paths.pdb_id, []))
            except Exception:
                _post = -1
            log.info("[ph_ensemble.debug] after-bridge map_len[%s]=%d", paths.pdb_id, _post)
            # ---------------------------------------------------

            log.info("[ph_ensemble.done] ok=True")
            return manifest_path




        except Exception as exc:
            log.error("[ph_ensemble.error] %s", exc)
            log.info("[ph_ensemble.done] ok=False")
            return None
        finally:
            if prev_cfg is not None:
                automate_protein_prep.config = prev_cfg

    if cleaned_pdb_path.exists() and receptor_pdbqt_path.exists() and not force_reprocess:
        logger.info("Reusing existing cleaned PDB and receptor PDBQT.")
        try:
            if bool(cfg.get("RECEPTOR_SANITY_CHECK", True)) and not receptor_sanity_check(str(receptor_pdbqt_path)):
                logger.warning("Receptor sanity check failed (cached receptor).")
                prepare_receptor.last_provenance = "cache_reuse_failed"
                return None, None
        except Exception as _e:
            logger.warning(f"Receptor sanity check skipped due to error: {_e}")
        cleaned_norm = norm(cleaned_pdb_path)
        receptor_norm = norm(receptor_pdbqt_path)
        if bool(cfg.get("PH_ENSEMBLE_IN_PREP", False)):
            _build_ph_ensemble(cleaned_norm)
            mode = _ph_ligand_mode(cfg)
            if mode != "off":
                try:
                    context_ph_values = _ph_values_from_context(cleaned_pdb_path)
                except Exception as e:
                    logger.warning("[ph_ligand] failed to load context pH values; skipping ligand enumeration: %s", e)
                    context_ph_values = []

                if context_ph_values:
                    window: set[float] = set()
                    for p in context_ph_values:
                        if p is None:
                            continue
                        try:
                            p_val = float(p)
                        except Exception:
                            continue
                        for delta in (-1.0, 0.0, +1.0):
                            v = p_val + delta
                            if v < 3.0 or v > 10.5:
                                continue
                            window.add(round(v, 1))

                    ligand_ph_values = sorted(window)
                    if ligand_ph_values:
                        logger.info(
                            "[ph_ligand] mode=%s context_pH=%s ligand_pH_window=%s",
                            mode,
                            ",".join(f"{p:.1f}" for p in sorted(context_ph_values)),
                            ",".join(f"{p:.1f}" for p in ligand_ph_values),
                        )
                        _ctrl_roots, noncontrol_roots = _lib_roots_for_pdb(
                            cfg, paths.pdb_id.upper(), paths, logger
                        )
                        ph_ligand_root = noncontrol_roots[0] if noncontrol_roots else None
                        logger.info(
                            "[ph_ligand.context.bridge] pdb=%s root_dir=%s requested_ph=%s",
                            paths.pdb_id,
                            ph_ligand_root,
                            ligand_ph_values,
                        )
                        if ph_ligand_root is None or not ph_ligand_root.exists():
                            logger.info(
                                "[ph_ligand.context.bridge.skip] no valid ligand root; skipping microstate priming",
                            )
                        else:
                            try:
                                enumerate_ligands_for_docking(
                                    requested_ph_values=ligand_ph_values,
                                    root_dir=ph_ligand_root,
                                    microstate_dedup=True,
                                    force=False,
                                    cfg=cfg,
                                    pdb_id=paths.pdb_id,
                                )
                            except Exception as e:
                                logger.warning("[ph_ligand] ligand enumeration failed (non-fatal): %s", e)
                    else:
                        logger.info("[ph_ligand] context pH values present but window is empty after clamping; skipping ligand enumeration")
                else:
                    logger.info("[ph_ligand] no context pH values available; ligand enumeration skipped")
        prepare_receptor.last_provenance = "cache_reuse"
        return cleaned_norm, receptor_norm
    
    # --- PH_ENSEMBLE gating of legacy protonation ---
    if bool(cfg.get("PH_ENSEMBLE", False)):
        os.environ["A2_SKIP_PDB2PQR"] = "1"
        logger.info("[ph_ensemble] enabling ensemble mode: A2_SKIP_PDB2PQR=1 for cleaning stage")
    else:
        os.environ.pop("A2_SKIP_PDB2PQR", None)
        logger.info("[ph_ensemble] disabled; legacy cleaning path unchanged")
    # ------------------------------------------------
    # Fresh prep path: clean PDB then create PDBQT into the variant-aware target
    try:
        cleaned_pdb = automate_protein_prep.clean_pdb(
            pdb_file=str(paths.input_pdb_path),
            output_root=str(Path(cfg["OUTPUT_DIR"])),  # processed_pdbs root; module lays out subdirs
            logger=logger,
        )
    except Exception as e:
        logger.warning(f"Protein cleaning failed: {e}")
        prepare_receptor.last_provenance = "clean_failed"
        return None, None

    if not cleaned_pdb or not Path(cleaned_pdb).exists():
        logger.warning("Protein cleaning did not produce a cleaned PDB.")
        prepare_receptor.last_provenance = "clean_failed"
        return None, None

    # Relocate cleaned PDB into the variant receptor dir if needed
    try:
        if Path(cleaned_pdb).resolve() != cleaned_pdb_path.resolve():
            cleaned_pdb_path.parent.mkdir(parents=True, exist_ok=True)
            from shutil import copy2
            copy2(str(cleaned_pdb), str(cleaned_pdb_path))
            cleaned_pdb = str(cleaned_pdb_path)
        else:
            cleaned_pdb = str(cleaned_pdb_path)
    except Exception as e:
        logger.warning(f"Could not relocate cleaned PDB: {e}")

    # Build protein pH ensemble if requested
    ph_ensemble_in_prep = bool(cfg.get("PH_ENSEMBLE_IN_PREP", False))
    if ph_ensemble_in_prep:
        _build_ph_ensemble(cleaned_pdb)

    # Optional ligand microstate priming based on context pH
    mode = _ph_ligand_mode(cfg)
    if ph_ensemble_in_prep and mode != "off":
        try:
            # Use the existing context-pH helper to recover the canonical pH values
            # used to build the ensemble, instead of parsing any filenames.
            context_ph_values = _ph_values_from_context(cleaned_pdb_path)
        except Exception as e:
            logger.warning("[ph_ligand] failed to load context pH values; skipping ligand enumeration: %s", e)
            context_ph_values = []

        if context_ph_values:
            # Build a ±1 pH window union across all context pHs, clamped to [3.0, 10.5]
            window: set[float] = set()
            for p in context_ph_values:
                if p is None:
                    continue
                try:
                    p_val = float(p)
                except Exception:
                    continue
                for delta in (-1.0, 0.0, +1.0):
                    v = p_val + delta
                    if v < 3.0 or v > 10.5:
                        continue
                    # round to 1 decimal place to match microstate registry convention
                    window.add(round(v, 1))

            ligand_ph_values = sorted(window)
            if ligand_ph_values:
                logger.info(
                    "[ph_ligand] mode=%s context_pH=%s ligand_pH_window=%s",
                    mode,
                    ",".join(f"{p:.1f}" for p in sorted(context_ph_values)),
                    ",".join(f"{p:.1f}" for p in ligand_ph_values),
                )
                _ctrl_roots, noncontrol_roots = _lib_roots_for_pdb(cfg, paths.pdb_id.upper(), paths, logger)
                ph_ligand_root = noncontrol_roots[0] if noncontrol_roots else None

                if ligand_ph_values and ph_ligand_root is not None and ph_ligand_root.exists():
                    logger.info(
                        "[ph_ligand.context.bridge] pdb=%s root_dir=%s requested_ph=%s",
                        paths.pdb_id,
                        ph_ligand_root,
                        ligand_ph_values,
                    )
                    try:
                        enumerate_ligands_for_docking(
                            requested_ph_values=ligand_ph_values,
                            root_dir=ph_ligand_root,
                            microstate_dedup=True,
                            force=False,
                        )
                    except Exception as e:
                        logger.warning("[ph_ligand] ligand enumeration failed (non-fatal): %s", e)
                else:
                    logger.info(
                        "[ph_ligand.context.bridge.skip] no valid ligand root or empty window; skipping microstate priming"
                    )
            else:
                logger.info("[ph_ligand] context pH values present but window is empty after clamping; skipping ligand enumeration")
        else:
            logger.info("[ph_ligand] no context pH values available; ligand enumeration skipped")

    try:
        provenance = getattr(automate_protein_prep, "get_clean_provenance", lambda: "clean_pdb")()
    except Exception:
        provenance = "clean_pdb"
    prepare_receptor.last_provenance = provenance

    # Generate receptor PDBQT directly at the variant-aware path
    try:
        ok = automate_protein_prep.run_prepare_receptor(
            input_pdb=cleaned_pdb,
            output_pdbqt=str(receptor_pdbqt_path),
            cfg=cfg
        )
        receptor_pdbqt = str(receptor_pdbqt_path) if ok else None
    except Exception as e:
        logger.warning(f"Receptor PDBQT prep failed: {e}")
        receptor_pdbqt = None

    if not receptor_pdbqt or not Path(receptor_pdbqt).exists():
        logger.warning("Receptor PDBQT was not created.")
        prepare_receptor.last_provenance = "clean_failed"
        return None, None

    try:
        if Path(receptor_pdbqt).resolve() != receptor_pdbqt_path.resolve():
            from shutil import copy2
            receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
            copy2(receptor_pdbqt, receptor_pdbqt_path)
            receptor_pdbqt = str(receptor_pdbqt_path)
    except Exception as e:
        logger.warning(f"Could not relocate receptor PDBQT: {e}")

    try:
        if bool(cfg.get("RECEPTOR_SANITY_CHECK", True)):
            ok = receptor_sanity_check(receptor_pdbqt)
            if not ok:
                logger.warning("Receptor sanity check failed (too few atoms or zero coords).")
                prepare_receptor.last_provenance = "clean_failed"
                return None, None
    except Exception as _e:
        logger.warning(f"Receptor sanity check skipped due to error: {_e}")


    return norm(cleaned_pdb), norm(receptor_pdbqt)


# ---- Multi-control center selection via crystallographic controls ----


def _ensure_ctrl_vina_manifest(cfg, paths, stage_name, variant_token, ph_label, legacy_mode, logger):
    run_id = cfg["RUN_ID"]
    stage_dir = paths.configs_stage_dir(run_id, variant_token, stage_name, ph_label)
    manifest_path = stage_dir / "vina.json"
    stage_dir.mkdir(parents=True, exist_ok=True)
    precreated = 0
    if not manifest_path.exists():
        payload = {
            "run_id": run_id,
            "pdb_id": paths.pdb_id,
            "stage": stage_name,
            "variant": variant_token,
            "ph": ph_label,
            "legacy": bool(legacy_mode),
            "entries": [],
        }
        tmp = manifest_path.with_suffix(".part")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, manifest_path)
        precreated = 1
    logger.info("[ctrl.vina.json] precreated=%d path=%s", precreated, manifest_path)
    return manifest_path


def _finalize_ctrl_vina_manifest(
    cfg,
    paths,
    stage_name,
    variant_token,
    ph_label,
    legacy_mode,
    manifest_path,
    jobs,
):
    if manifest_path is None:
        return
    run_id = cfg["RUN_ID"]
    stage_dir = manifest_path.parent
    receptor_path = receptor_file(
        paths.pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_root = docked_dir(
        paths.pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_root.mkdir(parents=True, exist_ok=True)
    out_dir = out_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    for job in jobs:
        lig_base = Path(job["ligand_path"]).stem
        entry = {
            "ligand": lig_base,
            "config": str(stage_dir / f"{lig_base}_{stage_name}.txt"),
            "out": str(out_dir / f"{lig_base}_{stage_name}.pdbqt"),
            "receptor": str(receptor_path),
        }
        entries.append(entry)

    entries.sort(key=lambda e: e["ligand"])

    payload = {
        "run_id": run_id,
        "pdb_id": paths.pdb_id,
        "stage": stage_name,
        "variant": variant_token,
        "ph": ph_label,
        "legacy": bool(legacy_mode),
        "entries": entries,
    }

    tmp = manifest_path.with_suffix(".part")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, manifest_path)


def _ctrl_quick_file_sig(pth: str) -> str:
    try:
        p = Path(pth)
        if not p.exists():
            return "exists=False size=-1 sha=0000000000"
        sz = p.stat().st_size
        h = hashlib.sha1()
        with open(p, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith(("ATOM", "HETATM")):
                    h.update(ln[12:54].encode("utf-8", "ignore"))
        return f"exists=True size={sz} sha={h.hexdigest()[:10]}"
    except Exception:
        return "sig=unavailable"


def _ctrl_best_model_to_pdb(pdbqt_file: str, obabel: str):
    from pathlib import Path as _Path
    import tempfile
    import subprocess
    import shutil as _sh

    best_e = None
    best_chunk = None
    pdbqt_path = _Path(pdbqt_file)
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
        chunk = []
        in_model = False
        for ln in fh:
            u = ln.strip().upper()
            if u.startswith("MODEL"):
                chunk = [ln]
                in_model = True
            elif u.startswith("ENDMDL"):
                chunk.append(ln)
                in_model = False
                for cl in chunk:
                    if "REMARK VINA RESULT" in cl.upper():
                        try:
                            e = float(cl.strip().split()[3])
                            if (best_e is None) or (e < best_e):
                                best_e = e
                                best_chunk = chunk[:]
                        except Exception:
                            pass
            else:
                if in_model:
                    chunk.append(ln)
    if best_chunk is None:
        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if "REMARK VINA RESULT" in ln.upper():
                        best_e = float(ln.strip().split()[3])
                        break
            best_chunk = None
        except Exception:
            return None, None

    td = _Path(tempfile.mkdtemp(prefix="ctrl_redock_"))
    best_pdbqt = td / "best.pdbqt"
    if best_chunk:
        with open(best_pdbqt, "w", encoding="utf-8") as out:
            out.writelines(best_chunk)
    else:
        _sh.copy2(pdbqt_path, best_pdbqt)
    out_pdb = td / "best.pdb"
    try:
        subprocess.run(
            [obabel, "-ipdbqt", str(best_pdbqt), "-opdb", "-O", str(out_pdb)],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None, best_e
    return (out_pdb if out_pdb.exists() else None), best_e


def _ctrl_redock_job(payload: dict) -> dict:
    cfg = payload["cfg"]
    stage_info = dict(payload["stage_info"])
    conf_path, out_path = emit_vina_config(
        cfg,
        payload["pdb_id"],
        payload["receptor_pdbqt"],
        tuple(payload["center"]),
        tuple(payload["box_size"]),
        payload["ligand_path"],
        payload["stage_name"],
        stage_info,
        int(payload["threads_per_job"]),
        logger=None,
        variant=payload["variant"],
        ph_token=payload["ph"],
        legacy=payload["legacy"],
        skip_manifest_if_exists=bool(payload.get("skip_manifest", False)),
    )
    try:
        _, score = run_docking_task(
            payload["vina_exe"],
            str(conf_path),
            payload["ligand_name"],
            str(out_path),
        )
    except Exception:
        score = None
    best_pdb, best_e = _ctrl_best_model_to_pdb(str(out_path), payload["obabel"])
    return {
        "order": payload["order"],
        "base": payload["base"],
        "center": tuple(payload["center"]),
        "ligand_name": payload["ligand_name"],
        "best_pdb": str(best_pdb) if best_pdb else None,
        "best_e": best_e,
        "score": score,
        "out_path": str(out_path),
    }


def select_center_via_control_redock(
    cfg,
    paths,
    receptor_pdbqt,
    logger,
    *,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    legacy: bool = False,
):
    """
    Returns (center_tuple, (24.0,24.0,24.0)) or (None, None).
    - If multiple controls and max pairwise centroid distance <= CONTROL_CENTER_CLOSE_MAX_A -> average (consensus).
    - Else (far apart), redock each prepped control and pick lowest RMSD vs its crystal.
    - If no controls or all redocks fail, returns (None, None) to signal P2Rank fallback.
    """
    import numpy as _np
    import math as _math
    from pathlib import Path as _Path
    from run_vina import run_docking_task as _run_dock

    def _find_control_pdbs(d: _Path) -> list[_Path]:
        return sorted([p for p in d.glob("*.pdb") if p.is_file()])

    def _centroid_from_pdb(p: _Path):
        xs, ys, zs = [], [], []
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    if ln.startswith(("ATOM","HETATM")):
                        try:
                            xs.append(float(ln[30:38])); ys.append(float(ln[38:46])); zs.append(float(ln[46:54]))
                        except Exception:
                            continue
        except Exception:
            return None
        if not xs: return None
        return (float(_np.mean(xs)), float(_np.mean(ys)), float(_np.mean(zs)))

    # discover crystal controls in preferred locations (current + legacy sibling)
    ctrl_pdbs = _find_control_pdbs(paths.ligand_output_dir)
    logger.info(f"[control-redock] controls_found={len(ctrl_pdbs)} dir={paths.ligand_output_dir}")
    if not ctrl_pdbs:
        legacy = paths.ligand_output_dir.parent.parent / f"{paths.pdb_id}_NOLIG" / "ligands_raw"
        if legacy.exists():
            ctrl_pdbs = _find_control_pdbs(legacy)

    if not ctrl_pdbs:
        logger.warning("[control-redock] No extracted control PDBs present; skipping redock.")
        return None, None  # let caller go to P2Rank directly


    policy = str(cfg.get("CONTROL_CENTER_POLICY", "best_redock")).lower().strip()
    variant_env = variant if variant is not None else (os.environ.get("APO_HOLO_VARIANT", "") or "")
    variant_token = (str(variant_env).strip().upper() or None)
    ph_label = str(ph_token).strip() if ph_token is not None else None
    legacy_mode = bool(legacy)
    thr = float(cfg.get("CONTROL_CENTER_CLOSE_MAX_A", 8.0))

    # compute centroids + pairwise spread
    centroids = {}
    for p in ctrl_pdbs:
        c = _centroid_from_pdb(p)
        if c: centroids[p.stem.split("_stage")[0]] = c

    bases = list(centroids.keys())
    coords = [centroids[b] for b in bases]
    def _dist(a,b):
        return float(((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5)
    max_delta = 0.0
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            d = _dist(coords[i], coords[j])
            if d > max_delta: max_delta = d

    if not coords:
        raise RuntimeError("[control-centers] No control centroids available; cannot select center.")
    logger.info(f"[control-centers] n={len(coords)} max?={max_delta:.2f}A policy={policy} thr={thr:.2f}A")
    for b, c in zip(bases, coords):
        logger.debug(f"[control-centers] {b}: ({c[0]:.3f},{c[1]:.3f},{c[2]:.3f})")

    # single-control or simple policies
    if len(coords) == 1 and policy != "best_redock":
        center = coords[0]
        logger.info(f"[Control-center] chosen={bases[0]} center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f}) box=(24,24,24)")
        return center, (24.0,24.0,24.0)

    if policy == "first":
        center = coords[0]
        logger.info(f"[Control-center] chosen={bases[0]} center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f}) box=(24,24,24)")
        return center, (24.0,24.0,24.0)
    # if best_redock, skip consensus short-circuit:
    if policy == "best_redock":
        pass  # fall through to redock block below

    elif max_delta <= thr:
        # consensus average when controls are close
        c = (float(_np.mean([x for x,_,_ in coords])),
             float(_np.mean([y for _,y,_ in coords])),
             float(_np.mean([z for _,_,z in coords])))
        logger.info(f"[Control-center] chosen=consensus center=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}) box=(24,24,24)")
        return c, (24.0,24.0,24.0)

    elif policy == "average_when_close":
        # far apart ? fall back to first per spec
        center = coords[0]
        logger.info(f"[Control-center] chosen={bases[0]} center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f}) box=(24,24,24)")
        return center, (24.0,24.0,24.0)

    # best_redock path (controls far apart)
    control_lookup = build_control_lookup(paths)  # base -> reference path
    # collect prepped control pdbqts (per-protein dir and global output dir)
    prepped_dirs = [paths.prepped_ligands_dir, _Path(str(cfg.get("OUTPUT_LIGANDS_DIR", ""))) / paths.pdb_id]
    cand_pdbqts = []
    logger.info(f"[control-redock] search_prepped_dirs={[str(d) for d in prepped_dirs]}")

    seen = set()
    for root in prepped_dirs:
        if not root or not _Path(root).exists(): continue
        for p in _Path(root).glob("*.pdbqt"):
            base = p.stem.split("_stage")[0].split(".sanitized")[0]
            if base in centroids and base in control_lookup and base not in seen:
                cand_pdbqts.append(p); seen.add(base)
    logger.info(f"[control-redock] candidates={len(cand_pdbqts)}")

    if not cand_pdbqts:
        logger.warning("[control-redock] No prepped control PDBQTs found; redock impossible (will fall back).")
        return None, None

    ex = int(cfg.get("CTRL_REDOCK_EXHAUSTIVENESS", 24)) #AAA CHANGE FOR test RUNS
    nm = int(cfg.get("CTRL_REDOCK_NMODES", 9)) # AAA CHANGE FOR test RUNS
    threads_per_vina = int(cfg.get("THREADS_PER_VINA_CTRL", 8))  # control redock uses its own threads default=8
    vina_exe = str(cfg.get("VINA_EXE") or cfg.get("VINA_PATH") or "vina")
    obabel = str(cfg.get("OPENBABEL_PATH") or "obabel")


    best = None  # (rmsd, score, base, center_tuple)

    cpu_total = int(cfg.get("CPU", os.cpu_count() or 1) or 1)
    cpu_total = max(1, cpu_total)
    max_parallel = int(cfg.get("MAX_PARALLEL_JOBS", 1) or 1)
    max_parallel = max(1, max_parallel)
    workers = max(1, min(cpu_total, max_parallel))
    job_threads = max(1, min(threads_per_vina, max(1, cpu_total // workers)))
    total_threads = workers * job_threads
    logger.info(
        "[ctrl.parallel] CPU=%d MAX_PARALLEL_JOBS=%d workers=%d threads_per_job=%d total_threads=%d",
        cpu_total,
        max_parallel,
        workers,
        job_threads,
        total_threads,
    )
    if job_threads < threads_per_vina:
        logger.info(
            "[ctrl.parallel] adjust threads_per_vina old=%d new=%d",
            threads_per_vina,
            job_threads,
        )

    stage_name = "ctrl_redock"
    parallel_enabled = workers > 1 and len(cand_pdbqts) > 1

    manifest_path = None

    if not parallel_enabled:
        for lig_pdbqt in cand_pdbqts:
            base = lig_pdbqt.stem.split("_stage")[0].split(".sanitized")[0]
            center = centroids.get(base)
            if not center:
                ref = control_lookup.get(base)
                if ref and ref.suffix.lower() == ".pdb":
                    center = _centroid_from_pdb(ref)
            if not center:
                continue

            stage_info = {"exhaustiveness": ex, "num_modes": nm}
            if cfg.get("FAST_MODE"):
                stage_info["exhaustiveness"] = 1
            logger.info(
                "[emit.debug] pdb=%s stage=%s variant=%s ph=%s",
                paths.pdb_id,
                stage_name,
                variant_token or "None",
                ph_label or "None",
            )
            conf_path, out_path = emit_vina_config(
                cfg,
                paths.pdb_id,
                receptor_pdbqt,
                center,
                (24.0, 24.0, 24.0),
                str(lig_pdbqt),
                stage_name,
                stage_info,
                job_threads,
                logger=None,
                variant=variant_token,
                ph_token=ph_label,
                legacy=legacy_mode,
            )
            try:
                _, score = _run_dock(vina_exe, conf_path, lig_pdbqt.name, out_path)
            except Exception:
                score = None

            best_pdb, best_e = _ctrl_best_model_to_pdb(str(out_path), obabel)
            ref_path = control_lookup.get(base)
            rmsd = float("inf")

            if best_pdb and ref_path:
                logger.info(f"[rmsd.debug] ref={ref_path} | {_ctrl_quick_file_sig(str(ref_path))}")
                logger.info(f"[rmsd.debug] dock={best_pdb} | {_ctrl_quick_file_sig(str(best_pdb))}")
                same_file = (Path(ref_path).resolve() == Path(best_pdb).resolve())
                if same_file:
                    logger.warning("[rmsd.debug] ref and dock paths resolve to the same file! RMSD=0.0 is expected.")
                try:
                    rmsd = compute_rmsd(str(ref_path), str(best_pdb))
                except Exception as e:
                    rmsd = float("inf")
                    logger.exception(f"[rmsd.debug] compute_rmsd failed: {e}")

            e_print = best_e if (best_e is not None) else (score if score is not None else float("nan"))
            logger.info(f"[control-redock] lig={lig_pdbqt.name} rmsd={rmsd:.2f}A score={e_print if e_print is not None else float('nan')} kcal/mol")

            if _math.isfinite(rmsd):
                if (best is None) or (rmsd < best[0]) or (rmsd == best[0] and (e_print is not None) and (best[1] is None or e_print < best[1])):
                    best = (rmsd, e_print if e_print is not None else None, base, center)
    else:
        manifest_path = _ensure_ctrl_vina_manifest(
            cfg,
            paths,
            stage_name,
            variant_token,
            ph_label,
            legacy_mode,
            logger,
        )
        logger.info("[ctrl.parallel] write_once=on path=%s", manifest_path)
        cfg_payload = dict(cfg)
        jobs = []
        for order, lig_pdbqt in enumerate(cand_pdbqts):
            base = lig_pdbqt.stem.split("_stage")[0].split(".sanitized")[0]
            center = centroids.get(base)
            if not center:
                ref = control_lookup.get(base)
                if ref and ref.suffix.lower() == ".pdb":
                    center = _centroid_from_pdb(ref)
            if not center:
                continue
            stage_info = {"exhaustiveness": ex, "num_modes": nm}
            if cfg.get("FAST_MODE"):
                stage_info["exhaustiveness"] = 1
            logger.info(
                "[emit.debug] pdb=%s stage=%s variant=%s ph=%s",
                paths.pdb_id,
                stage_name,
                variant_token or "None",
                ph_label or "None",
            )
            jobs.append(
                {
                    "order": order,
                    "cfg": cfg_payload,
                    "pdb_id": paths.pdb_id,
                    "receptor_pdbqt": str(receptor_pdbqt),
                    "center": tuple(center),
                    "box_size": (24.0, 24.0, 24.0),
                    "ligand_path": str(lig_pdbqt),
                    "ligand_name": lig_pdbqt.name,
                    "stage_name": stage_name,
                    "stage_info": stage_info,
                    "threads_per_job": job_threads,
                    "variant": variant_token,
                    "ph": ph_label,
                    "legacy": legacy_mode,
                    "vina_exe": vina_exe,
                    "obabel": obabel,
                    "base": base,
                    "skip_manifest": True,
                }
            )

        if not jobs:
            return None, None

        results = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            future_map = {pool.submit(_ctrl_redock_job, job): job for job in jobs}
            for fut in as_completed(future_map):
                res = fut.result()
                results.append(res)

        results.sort(key=lambda r: r["order"])
        _finalize_ctrl_vina_manifest(
            cfg,
            paths,
            stage_name,
            variant_token,
            ph_label,
            legacy_mode,
            manifest_path,
            jobs,
        )
        for res in results:
            base = res["base"]
            center = tuple(res["center"])
            ref_path = control_lookup.get(base)
            best_pdb = res["best_pdb"]
            rmsd = float("inf")
            if best_pdb and ref_path:
                logger.info(f"[rmsd.debug] ref={ref_path} | {_ctrl_quick_file_sig(str(ref_path))}")
                logger.info(f"[rmsd.debug] dock={best_pdb} | {_ctrl_quick_file_sig(str(best_pdb))}")
                same_file = (Path(ref_path).resolve() == Path(best_pdb).resolve())
                if same_file:
                    logger.warning("[rmsd.debug] ref and dock paths resolve to the same file! RMSD=0.0 is expected.")
                try:
                    rmsd = compute_rmsd(str(ref_path), str(best_pdb))
                except Exception as e:
                    rmsd = float("inf")
                    logger.exception(f"[rmsd.debug] compute_rmsd failed: {e}")

            best_e = res["best_e"]
            score = res["score"]
            e_print = best_e if (best_e is not None) else (score if score is not None else float("nan"))
            logger.info(f"[control-redock] lig={res['ligand_name']} rmsd={rmsd:.2f}A score={e_print if e_print is not None else float('nan')} kcal/mol")

            if _math.isfinite(rmsd):
                if (best is None) or (rmsd < best[0]) or (rmsd == best[0] and (e_print is not None) and (best[1] is None or e_print < best[1])):
                    best = (rmsd, e_print if e_print is not None else None, base, center)

    if best is None:
        return None, None

    chosen_center = best[3]
    logger.info(f"[Control-center] chosen={best[2]} center=({chosen_center[0]:.3f},{chosen_center[1]:.3f},{chosen_center[2]:.3f})")
    #BOX SIZE SPECIFIED HERE, NEED TO EDIT THIS TO CALCULATE BOX SIZE, LARGE BOX  SIZES DECREASE VINA  ACCURACY 
    return chosen_center, (24.0,24.0,24.0)

def detect_pocket(cleaned_pdb: str,
                  ligand_dir: Path,
                  logger: logging.Logger) -> Tuple[
    Optional[Tuple[float,float,float]],
    Optional[Tuple[float,float,float]],
    str  # source ("control" | "p2rank" | "none")
]:
    """
    Prefer control ligands for docking center/box. If none, fall back to P2Rank.
    """
    def _his_counts_within(pdb_path: str, center_xyz: tuple[float,float,float], r: float = 6.0) -> tuple[int,int,int]:
        HID = HIE = HIP = 0
        try:
            with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
                seen = set()
                cx, cy, cz = center_xyz
                for ln in fh:
                    if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                        continue
                    res = ln[17:20].strip().upper()  # residue name
                    if res not in {"HID","HIE","HIP"}:
                        continue
                    try:
                        x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                    except Exception:
                        continue
                    if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 <= r*r:
                        # key by (chain, resseq, resname) so we count each residue once
                        key = (ln[21].strip(), ln[22:26].strip(), res)
                        if key in seen:
                            continue
                        seen.add(key)
                        if res == "HID": HID += 1
                        elif res == "HIE": HIE += 1
                        elif res == "HIP": HIP += 1
        except Exception:
            pass
        return HID, HIE, HIP

    def find_control_pdbs(d: Path) -> list[Path]:
        return sorted([p for p in d.glob("*.pdb") if p.is_file()])

    # 1) Controls check in canonical ligands_raw
    ctrl_files = find_control_pdbs(ligand_dir)
    # --- AUDIT: summarize control centroids & policy ---
    def _centroid_of_pdb(p: Path) -> tuple[float,float,float] | None:
        xs, ys, zs = [], [], []
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if ln.startswith(("ATOM","HETATM")) and len(ln) >= 54:
                        xs.append(float(ln[30:38])); ys.append(float(ln[38:46])); zs.append(float(ln[46:54]))
        except Exception:
            return None
        if xs:
            return (float(np.mean(xs)), float(np.mean(ys)), float(np.mean(zs)))
        return None

    centers = [c for c in (_centroid_of_pdb(p) for p in ctrl_files) if c]
    dmax = 0.0
    if len(centers) >= 2:
        for i in range(len(centers)):
            for j in range(i+1, len(centers)):
                dx = centers[i][0] - centers[j][0]
                dy = centers[i][1] - centers[j][1]
                dz = centers[i][2] - centers[j][2]
                d = float((dx*dx + dy*dy + dz*dz) ** 0.5)
                if d > dmax: dmax = d
    policy = "first" if ctrl_files else "p2rank"
    logger.info("[control-centers] n=%d max?=%.2f A policy=%s", len(ctrl_files), dmax, policy)

    # Back-compat (read-only): if none found, check legacy sibling <PDB>_NOLIG/ligands_raw
    if not ctrl_files:
        # ligand_dir = .../<PDB>/ligands_raw
        pdb_root = ligand_dir.parent  # .../<PDB>
        legacy = pdb_root.parent / f"{pdb_root.name}_NOLIG" / "ligands_raw"
        if legacy.exists():
            ctrl_files = find_control_pdbs(legacy)
            if ctrl_files:
                logger.info(f"[Control-center] Found controls in legacy sibling: {legacy}")

    if ctrl_files:
        p = ctrl_files[0]
        xs, ys, zs = [], [], []
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith(("ATOM","HETATM")):
                    try:
                        x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                        xs.append(x); ys.append(y); zs.append(z)
                    except ValueError:
                        continue
        if xs:
            ctrl_center = (float(np.mean(xs)), float(np.mean(ys)), float(np.mean(zs)))
            box_size = (24.0, 24.0, 24.0)
            hid, hie, hip = _his_counts_within(cleaned_pdb, ctrl_center, r=6.0)
            logger.info("[reduce] his={'HID':%d,'HIE':%d,'HIP':%d} flips_near_box=%d", hid, hie, hip, 0)
            return ctrl_center, box_size, "control"

    # 2) Fallback to P2Rank
    center, box_size = detect_active_site(cleaned_pdb)
    if center:
        box_size = tuple(min(28.0, float(s)) for s in box_size)
        logger.info(f"[P2Rank] Using P2Rank center {center} with box {box_size}")
        hid, hie, hip = _his_counts_within(cleaned_pdb, center, r=6.0)
        logger.info("[reduce] his={'HID':%d,'HIE':%d,'HIP':%d} flips_near_box=%d", hid, hie, hip, 0)
        return center, box_size, "p2rank"
    else:
        logger.error("Active-site detection failed (no controls, P2Rank returned None).")
        return None, None, "none"


def _count_heavy_atoms_from_pdbqt(pdbqt_path: Path) -> int:
    heavy = 0
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            if element:
                if element.upper() != "H":
                    heavy += 1
            else:
                atom_name = line[12:16].strip()
                if not atom_name.upper().startswith("H"):
                    heavy += 1
    return heavy


from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize
import subprocess, tempfile


def _load_mol_any(pdbqt_path: Path, obabel_exe: str | None) -> Chem.Mol | None:
    base = pdbqt_path.with_suffix("")
    # prefer SDF, then MOL2, then PDB
    sdf = base.with_suffix(".sdf"); mol2 = base.with_suffix(".mol2"); pdb = base.with_suffix(".pdb")
    if sdf.exists():
        supp = Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True)
        for m in supp:
            if m: return m
    for fp, reader in [(mol2, Chem.MolFromMol2File), (pdb, Chem.MolFromPDBFile)]:
        if fp.exists():
            m = reader(str(fp), sanitize=True, removeHs=False)
            if m: return m
    # fallback: PDBQT -> SDF via obabel (Linux-friendly)
    if obabel_exe:
        with tempfile.TemporaryDirectory() as td:
            out = Path(td, "tmp.sdf")
            try:
                subprocess.check_call([obabel_exe, "-ipdbqt", str(pdbqt_path), "-osdf", "-O", str(out), "--retype", "--addh"])
                supp = Chem.SDMolSupplier(str(out), removeHs=False, sanitize=True)
                for m in supp:
                    if m: return m
            except Exception:
                return None
    return None


def _standardize(m: Chem.Mol) -> Chem.Mol:
    parent = rdMolStandardize.ChargeParent(m)   # neutralize/parent
    rdMolStandardize.Normalize(parent)          # FG normalization
    Chem.SanitizeMol(parent)
    return parent


# Build catalog with PAINS A/B/C
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
pains_catalog = FilterCatalog.FilterCatalog(params)


def _coerce_test_map(m) -> Dict[str, str]:
    import json as _json, ast as _ast

    if isinstance(m, dict):
        return {str(k).upper(): str(v) for k, v in m.items()}
    s = str(m).strip()
    if not s:
        return {}
    parsed = None
    try:
        parsed = _json.loads(s)
    except Exception:
        try:
            parsed = _ast.literal_eval(s)
        except Exception:
            parsed = {}
    return {str(k).upper(): str(v) for k, v in (parsed if isinstance(parsed, dict) else {}).items()}


def _resolve_test_mode(cfg) -> str:
    """
    Normalize TEST_MODE_ENABLE to one of: "off", "dud", "fda+dud".

    Accepts:
      - Booleans / bool-like strings for backwards compatibility:
          True  / "true" / "yes" / "on" / "1"  -> "fda+dud"
          False / "false" / "no"  / "off" / "0" / "" / None -> "off"
      - Explicit string modes:
          "off"        -> "off"
          "dud"        -> "dud"
          "fda+dud"    -> "fda+dud"
          "both"       -> "fda+dud"
          "fda_dud"    -> "fda+dud"
    Any unrecognized string should log a warning and fall back to "off".
    """
    #env has priority over config
    raw = os.environ.get("TEST_MODE_ENABLE", cfg.get("TEST_MODE_ENABLE", "off"))

    if isinstance(raw, bool):
        return "fda+dud" if raw else "off"

    s = str(raw).strip().lower()
    if s in ("", "0", "false", "no", "off", "none", "null"):
        return "off"
    if s in ("true", "yes", "on", "1"):
        return "fda+dud"
    if s in ("dud", "dud-only", "dud_only"):
        return "dud"
    if s in ("fda+dud", "dud+fda", "both", "fda_and_dud", "fda_dud"):
        return "fda+dud"

    print(f"[test-mode] WARNING: Unknown TEST_MODE_ENABLE={raw!r}; treating as 'off'.")
    return "off"

def _dedup_index_roots(seq: list[Path]) -> list[Path]:
    deduped: list[Path] = []
    seen: set[str] = set()
    for candidate in seq:
        if not candidate:
            continue
        path_obj = Path(candidate)
        try:
            key = str(path_obj.resolve())
        except Exception:
            key = str(path_obj)
        if key not in seen:
            seen.add(key)
            deduped.append(path_obj)
    return deduped


def _lib_roots_for_pdb(cfg: Dict, pdb_id: str, paths: Paths, logger: logging.Logger) -> tuple[list[Path], list[Path]]:
    extra_dirs = str(cfg.get("LIBRARY_EXTRA_DIRS", "")).strip()
    extra_paths: list[Path] = []
    if extra_dirs:
        for d in extra_dirs.split(";"):
            d = d.strip()
            if not d:
                continue
            p = Path(d)
            if p.exists():
                extra_paths.append(p)

    subdir_default = str(cfg.get("LIBRARY_SUBDIR_DEFAULT", "fda_library"))
    test_mode = _resolve_test_mode(cfg)

    maybe_map = cfg.get("TEST_LIBRARY_MAP", {})
    test_map: Dict[str, str] = {}
    test_map = _coerce_test_map(maybe_map)
    logger.info(f"[lib-roots.map] raw_type={type(maybe_map).__name__} keys={len(test_map)}")
    mapped_value = (test_map or {}).get(pdb_id)

    base_root = Path(
        cfg.get("OUTPUT_LIGANDS_DIR")
        or cfg.get("PREPPED_LIGANDS_ROOT")
        or "prepped_ligands"
    )

    roots_for_mode: list[Path]
    if test_mode == "off" or not mapped_value:
        roots_for_mode = [base_root / subdir_default]
        if test_mode in ("dud", "fda+dud") and not mapped_value:
            logger.warning(
                "[test-mode] PDB %s missing from TEST_LIBRARY_MAP; using default library=%s",
                pdb_id,
                subdir_default,
            )
    elif test_mode == "dud":
        roots_for_mode = [base_root / mapped_value]
    elif test_mode == "fda+dud":
        roots_for_mode = [base_root / mapped_value, base_root / subdir_default]
    else:
        roots_for_mode = [base_root / subdir_default]

    deduped_roots = _dedup_index_roots(roots_for_mode)
    allowed_noncontrol_roots: list[Path] = []
    for r in deduped_roots:
        if r.exists():
            allowed_noncontrol_roots.append(r)
        else:
            logger.warning("[ligands.test-roots] missing=%s", r)

    extra_paths_cfg: list[str] = cfg.get("EXTRA_LIGAND_ROOTS", []) or []
    for d in extra_paths_cfg:
        p = Path(d)
        if p.exists():
            allowed_noncontrol_roots.append(p)
        else:
            logger.warning("[ligands.extra-roots] missing=%s", p)

    allowed_noncontrol_roots.extend(extra_paths)
    allowed_noncontrol_roots = _dedup_index_roots(allowed_noncontrol_roots)

    logger.info(
        "[ligands.allowed-roots] test_mode=%s pdb=%s roots=%d",
        test_mode,
        pdb_id,
        len(allowed_noncontrol_roots),
    )
    if allowed_noncontrol_roots:
        cfg["_PH_LIGAND_ROOT"] = str(allowed_noncontrol_roots[0])
    else:
        cfg.pop("_PH_LIGAND_ROOT", None)
    cfg["_ALLOWED_NONCONTROL_ROOTS"] = [str(p) for p in allowed_noncontrol_roots]
    logger.info(
        "[ph_ligand.roots] test_mode=%s pdb=%s ph_root=%s noncontrol_roots=%s",
        test_mode,
        pdb_id,
        cfg.get("_PH_LIGAND_ROOT"),
        cfg.get("_ALLOWED_NONCONTROL_ROOTS"),
    )
    cfg["_TEST_MODE_EFFECTIVE"] = test_mode
    return [], allowed_noncontrol_roots


def prepare_and_filter_ligands(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[List[str], Dict[str, int], Dict[str, bool]]:
    """
    Gathers candidate ligands, keeps existing validation/PAINS logic, and
    filters the *non-control* pool to allowed library roots:

      - TEST_MODE_ENABLE="off"      -> OUTPUT_LIGANDS_DIR/<LIBRARY_SUBDIR_DEFAULT>
      - TEST_MODE_ENABLE="dud"      -> OUTPUT_LIGANDS_DIR/<mapped_subdir>
      - TEST_MODE_ENABLE="fda+dud"  -> OUTPUT_LIGANDS_DIR/<mapped_subdir> + OUTPUT_LIGANDS_DIR/<LIBRARY_SUBDIR_DEFAULT>

    Controls are *never* filtered out here.
    LIBRARY_EXTRA_DIRS remain included (unchanged).
    """
    # Keep existing prep step for extracted controls (harmless if nothing to do)
    prep_ligands_from_pdb(
        ligand_output_dir=paths.ligand_output_dir,
        ligands_mol2_dir=paths.ligands_mol2_dir,
        prepped_ligands_dir=paths.prepped_ligands_dir,
    )

    _control_roots, allowed_noncontrol_roots = _lib_roots_for_pdb(cfg, paths.pdb_id.upper(), paths, logger)
    per_index_roots: list[Path] = []
    if paths.prepped_ligands_dir:
        per_index_roots.append(paths.prepped_ligands_dir)

    index_library_roots: list[Path] = _dedup_index_roots(allowed_noncontrol_roots)
    per_index_roots = _dedup_index_roots(per_index_roots)
    index_roots = _dedup_index_roots(per_index_roots + index_library_roots)

    manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
    lib_index = cfg.get("_LIB_INDEX")
    if not isinstance(lib_index, LibraryIndex):
        lib_index = LibraryIndex(manifest_filename=manifest_filename, logger=logger)
        cfg["_LIB_INDEX"] = lib_index
    if index_roots:
        lib_index.load(index_roots)
    cfg["_LIB_INDEX_PER_ROOTS"] = [str(p) for p in per_index_roots]
    cfg["_LIB_INDEX_LIBRARY_ROOTS"] = [str(p) for p in index_library_roots]

    def _under(p: Path, root: Path) -> bool:
        try:
            p.resolve().relative_to(root.resolve())
            return True
        except Exception:
            return False

    def _enumerate_noncontrol_candidates_via_index(
        cfg: Dict,
        allowed_roots: list[Path],
        logger: logging.Logger,
    ) -> list[Path]:
        resolved_roots = _dedup_index_roots([Path(r) for r in allowed_roots if r])
        resolved_existing = [r for r in resolved_roots if r.exists()]
        logger.info(
            "[lib-roots] non-control roots = %s",
            [str(p) for p in resolved_existing],
        )

        candidates: list[Path] = []
        microstate_roots = [r for r in resolved_existing if (r / "microstates.json").exists()]
        if microstate_roots:
            try:
                all_ms: list[Path] = []
                for root in microstate_roots:
                    ms_paths = enumerate_ligands_for_docking(
                        requested_ph_values=None,
                        root_dir=root,
                        microstate_dedup=True,
                        force=False,
                    )
                    all_ms.extend(ms_paths)
                ms_filtered = [p for p in all_ms if any(_under(p, root) for root in resolved_existing)]
                candidates.extend(ms_filtered)
                logger.info(
                    "[lib-index.microstate] roots=%d ligands=%d",
                    len(microstate_roots),
                    len(ms_filtered),
                )
            except Exception as exc:
                logger.warning("[lib-index.microstate] error=%s", exc)

        manifest_candidates: list[Path] = []
        if isinstance(cfg.get("_LIB_INDEX"), LibraryIndex):
            index_obj: LibraryIndex = cfg["_LIB_INDEX"]
            for root in resolved_existing:
                manifest = getattr(index_obj, "_cache", {}).get(Path(root))
                if not manifest:
                    continue
                for rel in manifest.entries.values():
                    manifest_candidates.append(Path(root) / rel)
        if manifest_candidates:
            candidates.extend(manifest_candidates)
            logger.info(
                "[lib-index.manifest] roots=%d ligands=%d",
                len(resolved_existing),
                len(manifest_candidates),
            )

        if not candidates:
            logger.warning(
                "[lib-index] no usable index detected; falling back to filesystem scan under %d roots (this is slow)",
                len(resolved_existing),
            )
            seen: set[str] = set()
            for root in resolved_existing:
                for p in _iter_pdbqt_dirfirst(root):
                    pn = norm(p)
                    if pn not in seen:
                        seen.add(pn)
                        candidates.append(p)

        deduped: list[Path] = []
        seen_keys: set[str] = set()
        for p in candidates:
            key = norm(p)
            if key not in seen_keys:
                seen_keys.add(key)
                deduped.append(Path(p))
        return deduped

    scan_roots: list[Path] = []
    prepped_lig_root = paths.prepped_ligands_dir
    if prepped_lig_root and prepped_lig_root.exists():
        scan_roots.append(prepped_lig_root)

    scan_roots.extend(allowed_noncontrol_roots)
    scan_roots = _dedup_index_roots(scan_roots)

    single_selector = cfg.get("_EFFECTIVE_SINGLE_LIGAND")
    crawl_allowed = not bool(single_selector)
    logger.info(
        "[ligands.scan.guard] roots=%d crawl_allowed=%s selector=%s",
        len(scan_roots),
        crawl_allowed,
        single_selector,
    )

    seen: set[str] = set()
    controls: list[Path] = []
    per_protein_allow = {"controls", "reference"}
    if prepped_lig_root and prepped_lig_root.exists():
        logger.info(
            "[ligands.scan.root] root=%s allowed_subdirs=%s crawl=%s",
            prepped_lig_root,
            per_protein_allow,
            crawl_allowed,
        )
        for p in _iter_pdbqt_dirfirst(prepped_lig_root, allowed_subdirs=per_protein_allow):
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                controls.append(p)

    noncontrol_candidates = _enumerate_noncontrol_candidates_via_index(cfg, allowed_noncontrol_roots, logger)

    all_pdbqt_paths: list[Path] = []
    all_pdbqt_paths.extend(controls)
    all_pdbqt_paths.extend(noncontrol_candidates)

    logger.info("[ligands.scan] roots=%d found=%d", len(scan_roots), len(all_pdbqt_paths))
    cfg["ALL_LIGAND_PATHS"] = [str(p) for p in all_pdbqt_paths]

    if not all_pdbqt_paths:
        logger.warning("No .pdbqt ligands were found under the configured roots.")
        return [], {}, {}

    global_root = Path(cfg["OUTPUT_LIGANDS_DIR"]) if cfg.get("OUTPUT_LIGANDS_DIR") else None

    valid_pdbqt: Dict[str, Path] = {}
    for p in all_pdbqt_paths:
        try:
            if p.exists() and p.stat().st_size > 100:
                valid_pdbqt[norm(p)] = p
            else:
                logger.debug(f"Excluded malformed ligand (missing or too small): {p}")
        except Exception:
            logger.debug(f"Excluded malformed ligand (exception): {p}")

    cfg["ALL_LIGAND_PATHS_VALID"] = [str(p) for p in valid_pdbqt.values()]
    logger.info("[ligands.valid] count=%d", len(valid_pdbqt))

    controls_valid: list[Path] = []
    noncontrols: list[Path] = []
    for p in valid_pdbqt.values():
        if prepped_lig_root and _under(p, prepped_lig_root):
            controls_valid.append(p)
        else:
            noncontrols.append(p)

    filtered_noncontrols: list[Path] = []
    for p in noncontrols:
        keep = False
        for root in allowed_noncontrol_roots:
            if root.exists() and _under(p, root):
                keep = True
                break
        if keep:
            filtered_noncontrols.append(p)

    # --- Optional: build library manifests from scan results ---
    if cfg.get("LIBRARY_MANIFEST_BUILD_ON_SCAN"):
        try:
            manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
            by_root: Dict[Path, list[Path]] = {}

            for root in allowed_noncontrol_roots:
                root = Path(root)
                if not root.exists():
                    continue
                for lig in filtered_noncontrols:
                    lig_path = Path(lig)
                    if not lig_path.exists():
                        continue
                    if not _under(lig_path, root):
                        continue
                    by_root.setdefault(root, []).append(lig_path)

            for root, ligs in by_root.items():
                manifest_path = root / manifest_filename
                if manifest_path.exists():
                    logger.info(
                        "[lib-manifest.scan.skip] root=%s reason=exists path=%s",
                        str(root),
                        str(manifest_path),
                    )
                    continue

                if not ligs:
                    continue

                logger.info(
                    "[lib-manifest.scan.build] root=%s ligands=%d manifest=%s",
                    str(root),
                    len(ligs),
                    str(manifest_path),
                )

                entries: list[str] = []
                for lig in ligs:
                    try:
                        rel = Path(lig).resolve().relative_to(root.resolve())
                        entries.append(rel.as_posix())
                    except Exception:
                        entries.append(os.path.relpath(str(lig), str(root)))

                tmp_path = manifest_path.with_suffix(manifest_path.suffix + ".tmp")

                try:
                    lib_index = LibraryIndex(
                        manifest_filename=manifest_filename,
                        logger=logger,
                    )
                    if hasattr(lib_index, "write_manifest_for_root"):
                        lib_index.write_manifest_for_root(root, entries, tmp_path)
                    else:
                        import json

                        tmp_path.parent.mkdir(parents=True, exist_ok=True)
                        data = {"root": str(root), "entries": entries}
                        tmp_path.write_text(json.dumps(data, indent=2, sort_keys=True))

                    tmp_path.replace(manifest_path)
                except Exception:
                    logger.exception(
                        "[lib-manifest.scan.error] root=%s manifest=%s",
                        str(root),
                        str(manifest_path),
                    )
        except Exception:
            logger.exception("[lib-manifest.scan.error] unexpected failure during build_on_scan")
    # ------------------------------------------------------------

    # Merge back: controls (unaltered) + filtered non-controls
    if cfg.get("_EFFECTIVE_SINGLE_LIGAND") and cfg.get("_SINGLE_RESOLVED_PATH"):
        resolved_path = Path(cfg["_SINGLE_RESOLVED_PATH"])
        final_paths = controls_valid + [resolved_path]
        filtered_noncontrols = [resolved_path]
        logger.info("[single.fuel] resolved=%s controls=%d (blocking non-control pool)", cfg["_SINGLE_RESOLVED_PATH"], len(controls_valid))
    else:
        final_paths = controls_valid + filtered_noncontrols

    # --- PAINS flags (keep as before; default to {}) ---
    pains_flags: Dict[str, bool] = {}
    try:
        from rdkit import Chem
        from rdkit.Chem import FilterCatalog, rdMolStandardize
        # Build catalog once at module-level if you prefer; safe inline here too
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        pains_catalog = FilterCatalog.FilterCatalog(params)

        def _has_pains(pdbqt_path: Path) -> bool:
            try:
                # Try to locate a mol2 (or similar) neighbor if your original logic requires it.
                # Fallback: False (non-blocking).
                return False
            except Exception:
                return False

        for p in final_paths:
            pains_flags[p.stem] = _has_pains(p)
    except Exception:
        pains_flags = {}

    # Heavy atom counts (reuse your existing helper)
    heavy_atom_counts: Dict[str, int] = {}
    for p in final_paths:
        try:
            heavy_atom_counts[p.stem] = _count_heavy_atoms_from_pdbqt(p)
        except Exception:
            heavy_atom_counts[p.stem] = 0

    # Final return (stringify paths)
    ligands = [str(p) for p in final_paths]
    non_control_count = max(0, len(final_paths) - len(controls))
    logger.info(f"Selected ligands -> controls={len(controls)} + non-controls={non_control_count} = total={len(ligands)}")
    # GUARD: enforce .pdbqt-only pool
    bad = [p for p in ligands if not str(p).lower().endswith(".pdbqt")]
    if bad:
        raise ValueError(f"Ligand is not a .pdbqt file: {bad[0]}")
    return ligands, heavy_atom_counts, pains_flags



# ======================
# Phase 6: Docking loop
# ======================
def early_recenter_decision(
        i: int,
        scores: Dict[str, float],
        all_distances: List[float],
        box_size: Tuple[float, float, float],
        center: Tuple[float, float, float],
        stage1_original: List[str],
        attempts_used: int,
        params: RecenterParams,
        cfg: Dict,
        pdb_id: str,
        receptor_pdbqt: str,
        logger: logging.Logger,
        raw_docked: Dict[str, str],
        guard: GlobalCenterGuard,
        control_anchor_hit: bool
) -> Tuple[bool, Tuple[float, float, float], Tuple[float, float, float], List[str], int]:
    """
    Stage-1 heuristic for expanding box or recentering when everything docks far from the pocket.
    De-duped and control-anchored: will not fire if (a) a control validated this stage, (b) a global switch already
    occurred this stage, or (c) global switch cap reached.
    """
    if i != 0:
        return False, center, box_size, [], attempts_used
    if control_anchor_hit:
        logger.info("Early recenter skipped: control-anchored validation present.")
        return False, center, box_size, [], attempts_used
    if not guard.can_switch():
        logger.info("Early recenter skipped: global switch guard disallows further switches this stage/cap reached.")
        return False, center, box_size, [], attempts_used

    evaluated = len(all_distances)
    valid_count = len(scores)
    if evaluated < max(params.EARLY_RECENTER_MIN_EVAL, 15):
        logger.info(f"Early recenter skipped: evaluated={evaluated} < threshold.")
        return False, center, box_size, [], attempts_used

    far = sum(1 for d in all_distances if isinstance(d, (int, float)) and d > params.EARLY_RECENTER_FAR_A)
    far_ratio = far / evaluated if evaluated else 0.0
    med_dist = float(np.median(all_distances)) if all_distances else 0.0

    # Prefer a single mild box expand over recenter
    if params.ALLOW_BOX_EXPAND and (0.55 <= far_ratio < params.EARLY_RECENTER_RATIO) and (9.0 <= med_dist < params.EARLY_RECENTER_MEDIAN_A) and (valid_count == 0):
        box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
        new_box = tuple(min(box_cap, s + 4.0) for s in box_size)
        if new_box != box_size:
            logger.info(
                f"Borderline far_ratio={far_ratio:.2f}, median={med_dist:.1f} A -> "
                f"expand box to {new_box} and redo stage1."
            )
            # Note: not counted as a global switch
            return True, center, new_box, stage1_original[:], attempts_used

    if (far_ratio >= params.EARLY_RECENTER_RATIO) and (med_dist >= params.EARLY_RECENTER_MEDIAN_A) and (valid_count == 0):
        if attempts_used >= params.MAX_RECENTER_ATTEMPTS:
            logger.warning("Early recenter max attempts reached; proceeding without recenter.")
            return False, center, box_size, [], attempts_used

        logger.warning(f"Early recenter trigger: far_ratio={far_ratio:.2f}, median={med_dist:.1f}  , valid=0 -> recentering.")
        # >>> DOCKED PATHS PATCH START
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        # >>> DOCKED PATHS PATCH END
        fb_pose, new_center, _best_score, _chosen = attempt_fallback_recenter(
            fallback_ligands=raw_docked,
            receptor_pdbqt=receptor_pdbqt,
            docking_dir=str(paths.docked_pdb_root()),
            stage_name="stage1",
            pocket_center=center,
            logger=logger,
            exclude_basenames=set(),
        )
        if new_center is not None:
            attempts_used += 1
            new_box = tuple(min(float(cfg.get("BOX_SIZE_MAX_A", 28.0)), s) for s in box_size)
            guard.mark_switch()  # counts as a global switch
            logger.info("Re-running stage1 with new center and tightened box. [global switch]")
            return True, new_center, new_box, stage1_original[:], attempts_used
        logger.warning("Fallback could not produce a new center; proceeding without recenter.")

    return False, center, box_size, [], attempts_used


def select_ligands_for_next(
        docking_mode: str,
        i: int,
        stages: List[Dict],
        scores: Dict[str, float],
        logger: logging.Logger,
        base_pool_n: Optional[int] = None,          #  if provided, select % of this
        force_include: Optional[set] = None         #  always add these
) -> List[str]:
    if not scores:
        # Still allow force-carry if provided and next stage exists
        return sorted(force_include) if force_include else []

    schedule = {
        "discovery": [1.0, 0.1, 0.01, 0.001, 0.001],
        "polypharmacology": [1.0, 0.05, 0.005],   # i=1 -> next is stage3 uses 0.5%
    }.get(docking_mode, [1.0] * len(stages))

    pct = schedule[i + 1] if i + 1 < len(schedule) else 0.01

    # Use provided base if given (e.g., Stage1 pool size) -- otherwise fall back to valid-count
    pool_n = base_pool_n if (base_pool_n is not None) else len(scores)

    # Select K by the base pool, but cap at the number of valid scores available
    k_target = max(1, int(pool_n * pct))
    k = max(1, min(k_target, len(scores)))

    # take best k from valid scores
    next_list = [l for l, _ in sorted(scores.items(), key=lambda kv: kv[1])[:k]]

    # Force-carry: add any requested ligands (e.g., extracted controls) to the next stage
    if force_include:
        # maintain stable order: extend with any forced ligands not already selected
        in_set = set(next_list)
        forced_add = [l for l in sorted(force_include) if l not in in_set]
        next_list.extend(forced_add)
        if forced_add:
            logger.info(f"[Force-carry] Added {len(forced_add)} extracted ligands to next stage.")

    logger.info(
        f"Selected {k} by score (+{len(force_include or [])} forced) "
        f"= {len(next_list)} total ({pct * 100:.5f}% of base={pool_n})."
    )
    return next_list

# ======================
# Center selection (controls + discovery)
# ======================
@dataclass
class CenterDecision:
    new_center: Optional[Tuple[float, float, float]]
    reason: str = ""
    switchscore: float = 0.0
    promoted: bool = False


class CenterSelector:
    """
    Decides when to keep the current center (control-anchored) vs. switch to a newly discovered pocket.
    Uses clustering of valid pose centroids + SwitchScore, now with stronger control anchoring and guard.
    """

    def __init__(self, cfg: Dict, logger: logging.Logger,
                 control_stems: set,
                 heavy_atom_counts: Dict[str, int],
                 initial_center: Tuple[float, float, float]):
        self.cfg = cfg
        self.log = logger
        self.control_stems = {s.lower() for s in control_stems}
        self.heavy = heavy_atom_counts
        self.eps = float(cfg.get("CLUSTER_EPS_ANG", 3.5))
        self.mode = cfg.get("CENTER_MODE", "control-first").lower()
        self.ctrl_blacklist = {s.strip().upper() for s in cfg.get("CONTROL_BLACKLIST", "").split(",") if s.strip()}
        self.ctrl_min_heavy = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))
        self.allow_switch_from_control = bool(cfg.get("ALLOW_SWITCH_FROM_CONTROL", True))
        self.require_control_failure = bool(cfg.get("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False))
        self.threshold = float(cfg.get("SWITCH_SCORE_THRESHOLD", 0.7))
        self.away_from_control_boost = float(cfg.get("SWITCH_AWAY_FROM_CONTROL_MIN_BOOST", 2.5))
        self.hysteresis = float(cfg.get("SWITCH_SCORE_HYSTERESIS", 0.5))
        self.switch_history = []  # keep last few SwitchScores
        self.current_center = np.array(initial_center, float)

        self.curr_stats = {"median_score": None, "median_le": None, "valid_rate": 0.0}
        self.last_decision_had_control_anchor = False

    @staticmethod
    def _strip_stage(name: str) -> str:
        stem = Path(name).stem
        return stem.split("_stage")[0].lower()

    @staticmethod
    def _pdbqt_centroid(pdbqt_path: str) -> Optional[np.ndarray]:
        if not pdbqt_path or not os.path.exists(pdbqt_path):
            return None
        xs, ys, zs = [], [], []
        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        try:
                            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                        except Exception:
                            parts = line.split()
                            if len(parts) >= 8:
                                x = float(parts[5]); y = float(parts[6]); z = float(parts[7])
                            else:
                                continue
                        xs.append(x); ys.append(y); zs.append(z)
            if xs:
                return np.array([np.mean(xs), np.mean(ys), np.mean(zs)], dtype=float)
        except Exception:
            return None
        return None

    def _is_blacklisted_control_name(self, stem_upper: str) -> bool:
        return any(stem_upper.startswith(bad) for bad in self.ctrl_blacklist)

    def _is_control(self, lig_path: str) -> bool:
        stem = self._strip_stage(os.path.basename(lig_path))
        if stem.upper() and self._is_blacklisted_control_name(stem.upper()):
            return False
        ha = self.heavy.get(lig_path)
        if ha is not None and ha < self.ctrl_min_heavy:
            return False
        return (stem in self.control_stems)

    def _cluster(self, points: List[np.ndarray]) -> List[List[int]]:
        clusters: List[List[int]] = []
        for i, p in enumerate(points):
            placed = False
            for cl in clusters:
                c = np.mean([points[j] for j in cl], axis=0)
                if np.linalg.norm(p - c) <= self.eps:
                    cl.append(i); placed = True; break
            if not placed:
                clusters.append([i])
        return clusters

    def _median(self, arr):
        return float(np.median(arr)) if arr else None

    def _cluster_stats(self,
                       member_idxs: List[int],
                       ligs: List[str],
                       scores_map: Dict[str, float],
                       centroids: List[np.ndarray],
                       total_docked: int) -> dict:
        members = [ligs[i] for i in member_idxs]
        scores = [scores_map[m] for m in members if m in scores_map]
        les = []
        for m in members:
            s = scores_map.get(m)
            if s is None:
                continue
            ha = self.heavy.get(m)
            if ha and ha > 0:
                les.append((-s) / ha)  # s is negative kcal/mol
        center = np.mean([centroids[i] for i in member_idxs], axis=0)
        valid_rate = len(members) / max(1, total_docked)
        ctrl_hits = sum(1 for m in members if self._is_control(m))
        return {
            "center": tuple(center.tolist()),
            "n": len(members),
            "valid_rate": float(valid_rate),
            "median_score": self._median(scores),
            "median_le": self._median(les),
            "control_hits": ctrl_hits,
        }

    def _switch_score(self, valid_rate, score_boost_kcal, le_gain, pocketability=0.5) -> float:
        consensus_boost = max(0.0, min(1.0, score_boost_kcal / 3.0))
        le_norm = max(0.0, min(1.0, le_gain / 0.05))
        return 0.40*valid_rate + 0.20*consensus_boost + 0.20*pocketability + 0.20*le_norm

    def consider_switch(self,
                        stage_name: str,
                        scores: Dict[str, float],
                        validated_ligands: List[str],
                        raw_docked: Dict[str, str],
                        receptor_pdbqt: str,
                        current_center: Tuple[float, float, float],
                        guard: GlobalCenterGuard) -> CenterDecision:
        # If the guard forbids a switch this stage or we're capped out, exit early
        if not guard.can_switch():
            return CenterDecision(None, "guard_disallowed")

        if not validated_ligands:
            self.last_decision_had_control_anchor = False
            return CenterDecision(None, "no_validated_poses")

        ligs = sorted(set(validated_ligands))
        centroids = []
        ligs_kept = []
        for lig in ligs:
            pose = raw_docked.get(lig)
            c = self._pdbqt_centroid(pose) if pose else None
            if c is not None:
                centroids.append(c); ligs_kept.append(lig)

        if len(centroids) < 3:
            self.last_decision_had_control_anchor = False
            return CenterDecision(None, "too_few_centroids")

        clusters = self._cluster(centroids)
        total_docked = max(1, len(raw_docked))
        stats = [self._cluster_stats(cl, ligs_kept, scores, centroids, total_docked) for cl in clusters]

        curr_idx = None
        curr_center = np.array(current_center, float)
        for idx, st in enumerate(stats):
            if np.linalg.norm(np.array(st["center"]) - curr_center) <= self.eps:
                curr_idx = idx; break

        curr_median = stats[curr_idx]["median_score"] if curr_idx is not None else None
        curr_le = stats[curr_idx]["median_le"] if curr_idx is not None else None
        curr_valid = stats[curr_idx]["valid_rate"] if curr_idx is not None else 0.0

        curr_ctrl_hits = stats[curr_idx]["control_hits"] if curr_idx is not None else 0
        self.curr_stats = {"median_score": curr_median, "median_le": curr_le, "valid_rate": curr_valid}

        # Strong control anchoring (compute for THIS stage before using it):
        self.last_decision_had_control_anchor = (
                curr_ctrl_hits > 0 and curr_valid >= float(self.cfg.get("CONTROL_ANCHOR_MIN_VALID_RATE", 0.10))
        )

        # If policy requires control failure to switch, stop here while control anchors
        if self.require_control_failure and self.last_decision_had_control_anchor:
            return CenterDecision(None, "control_anchor_lock")

        # Choose best alternative cluster
        best_cand = None
        best_score = -1.0
        for idx, st in enumerate(stats):
            if idx == curr_idx:
                continue
            score_boost = 0.0
            le_gain = 0.0
            if curr_median is not None and st["median_score"] is not None:
                score_boost = (curr_median - st["median_score"])
            if curr_le is not None and st["median_le"] is not None:
                le_gain = st["median_le"] - curr_le

            if st["valid_rate"] < float(self.cfg.get("SWITCH_VALID_RATE_MIN", 0.40)):
                continue
            # Require bigger improvement if leaving a control-anchored cluster
            required_boost = float(self.cfg.get("SWITCH_SCORE_IMPROVE_MIN", 1.5))
            if self.last_decision_had_control_anchor:
                required_boost = max(required_boost, self.away_from_control_boost)

            if score_boost < required_boost:
                continue
            if le_gain < float(self.cfg.get("SWITCH_LE_GAIN_MIN", 0.02)):
                continue

            sw = self._switch_score(st["valid_rate"], score_boost, le_gain, pocketability=0.5)
            if sw > best_score:
                best_score = sw
                best_cand = st

        if not best_cand:
            return CenterDecision(None, "no_candidate_passed_gates")

        self.switch_history.append(best_score)
        promoted = best_score >= self.threshold
        if not promoted:
            return CenterDecision(None, f"below_threshold:{best_score:.2f}", best_score, False)

        # mark guard switch outside (caller), but flag that we want to promote
        return CenterDecision(new_center=tuple(best_cand["center"]),
                              reason=f"promote_new_center score={best_score:.2f}",
                              switchscore=best_score,
                              promoted=True)


# ======================
# Phase 7-8: Finalization
# ======================
def final_pose_validation_and_screenshots(
        cfg: Dict,
        pdb_id: str,
        stages: List[Dict],
        receptor_pdbqt: str,
        center: Tuple[float, float, float],
        validated_ligands_last: List[str],
        score_history: Dict[str, Dict[str, Dict]],
        cleaned_pdb: str,
        docking_mode: str,
        logger: logging.Logger,
        ph_label: Optional[str] = None,
) -> None:
    if not validated_ligands_last:
        return

    # >>> DOCKED PATHS PATCH START
    variant = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # >>> DOCKED PATHS PATCH END
    last_stage = stages[-1]["name"]
    stage_dir = paths.docked_stage_dir(variant, last_stage, ph_label)
    final_surface = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)

    if docking_mode == "polypharmacology":
        final_scores = score_history.get(last_stage, {})
        top_ligs = sorted(final_scores.items(), key=score_key)[:20]
        validated_ligands_last = [lig for lig, _ in top_ligs]
        logger.info(f"[Polypharmacology] Selected top {len(validated_ligands_last)} ligands for images/validation.")

    for lig in validated_ligands_last:
        out_path = stage_dir / f"{Path(lig).stem}_{last_stage}.pdbqt"
        if not out_path.exists():
            logger.warning(f"Pose file not found for {os.path.basename(lig)} -- likely filtered earlier.")
            continue
        try:
            filter_and_rewrite_poses_by_rmsd(
                str(out_path),
                rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
            )
        except Exception as e:
            logger.warning(f"Final RMSD filtering failed: {e}")

        best_model, best_valid_score = validate_all_poses(
            pdbqt_path=str(out_path),
            receptor_pdbqt=receptor_pdbqt,
            center=center,
            surface_coords=final_surface,
            validate_fn=validate_pose_pdbqt
        )
        if best_model:
            record_score(score_history, last_stage, lig, best_valid_score, True, reason="rescued_best_pose")
            print(f"{Path(lig).name} | {last_stage} rescued: {best_valid_score:.2f} kcal/mol (valid)")
        else:
            try:
                fallback_score = extract_best_score(str(out_path))
                record_score(score_history, last_stage, lig, fallback_score, False, reason="all_poses_invalid")
            except Exception:
                record_score(score_history, last_stage, lig, None, False, reason="all_poses_invalid_no_score")
            print(f"{Path(lig).name} | all poses invalid (kept for logs)")

        try:
            import subprocess
            top = validated_ligands_last[0]
            pose = stage_dir / f"{Path(top).stem}_{last_stage}.pdbqt"
            if pose.exists():
                out_prefix_root = paths.docked_variant_root(variant, ph_label)
                out_prefix = out_prefix_root / "top_pose"
                out_prefix.parent.mkdir(parents=True, exist_ok=True)

                # -- PyMOL screenshot block (Option A: -r + -d python) --
                cap_py = Path(__file__).with_name("capture_pose.py")

                py_cfg = str(cfg.get("PYMOL_PATH", "")).strip()
                pymol_exe = py_cfg if (py_cfg and Path(py_cfg).is_file()) else (shutil.which("pymol") or "pymol")

                d_arg = f"""python
                from __main__ import capture_pose
                capture_pose({repr(cleaned_pdb)}, {repr(str(pose))}, {repr(str(out_prefix))})
                python end
                quit
                """

                # -cq keeps PyMOL headless/quiet; keep -r to load helper script
                cmd = [pymol_exe, "-cq", "-r", str(cap_py), "-d", d_arg]
                print("Running PyMOL:", cmd)
                res = subprocess.run(cmd, capture_output=True, text=True)
                print("PyMOL stdout:", res.stdout)
                print("PyMOL stderr:", res.stderr)

        except Exception as e:
            logger.warning(f"Screenshot generation failed: {e}")


def _pose_path_for(csv_cfg: Dict, pdb_id: str, stage_name: str, lig_path: str) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    from pathlib import Path
    # >>> DOCKED PATHS PATCH START
    paths = make_paths(csv_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_label = (csv_cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None
    variant = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    stage_dir = paths.docked_stage_dir(variant, stage_name, ph_label)
    return str(stage_dir / f"{Path(lig_path).stem}_{stage_name}.pdbqt")
    # >>> DOCKED PATHS PATCH END


def write_scores_csv(cfg: Dict, pdb_id: str, score_history: Dict[str, Dict[str, Dict]]) -> str:
    import csv, math

    # >>> DOCKED PATHS PATCH START
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    dock_dir = paths.docked_variant_root(variant_token, ph_label)
    dock_dir.mkdir(parents=True, exist_ok=True)
    # >>> DOCKED PATHS PATCH END

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

    # --- Wide summary (unchanged shape) ---
    csv_out_wide = str(dock_dir / "docking_score_summary.csv")
    flat = {}
    for stage_name, stage_map in score_history.items():
        flat[stage_name] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(lig)
            s = rec.get("score", None)
            if rec.get("valid", False):
                flat[stage_name][lig_key] = s if s is not None else ""
            else:
                flat[stage_name][lig_key] = f"{s:.2f} (invalid)" if isinstance(s, (int, float)) else "(invalid)"
    write_score_summary_to_csv(
        flat,
        output_path=csv_out_wide,
        run_id=run_id_value,
        variant=variant_value if include_variant else None,
    )

    # --- Long format with self_rmsd added ---
    csv_out_long = str(dock_dir / "docking_score_long.csv")
    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        header = ["run_id"]
        if include_variant:
            header.append("variant")
        header.extend(["stage", "ligand", "score", "valid", "reason", "heavy_atoms", "le", "self_rmsd", "pains_flag"])
        writer.writerow(header)

        for stage_name, stage_map in score_history.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(lig)
                score = rec.get("score", None)
                valid = bool(rec.get("valid", False))
                reason = rec.get("reason", "")

                ha = rec.get("heavy_atoms", None)
                le = rec.get("le", None)
                pains_hit = rec.get("pains_flag", False)

                # Compute self-RMSD from the saved pose for this stage (if present)
                pose_path = _pose_path_for(cfg, pdb_id, stage_name, lig)
                if os.path.exists(pose_path):
                    try:
                        sr = compute_self_rmsd(pose_path)
                        sr_str = f"{sr:.2f}" if isinstance(sr, (int, float)) and math.isfinite(sr) else ""
                    except Exception:
                        sr_str = ""
                else:
                    sr_str = ""

                score_str = f"{score:.2f}" if isinstance(score, (int, float)) else ""
                ha_str = str(int(ha)) if isinstance(ha, (int, float)) else ""
                le_str = f"{le:.4f}" if isinstance(le, (int, float)) else ""
                reason_str = str(reason) if reason is not None else ""

                row = [run_id_value]
                if include_variant:
                    row.append(variant_value)
                row.extend([stage_name, lig_key, score_str, int(valid), reason_str, ha_str, le_str, sr_str, int(pains_hit)])
                writer.writerow(row)

    return csv_out_wide


def compute_ligand_efficiency(score: Optional[float], heavy_atoms: Optional[int]) -> Optional[float]:
    try:
        if score is None or heavy_atoms is None or int(heavy_atoms) <= 0:
            return None
        return float(-float(score) / int(heavy_atoms))
    except Exception:
        return None


def record_le(score_history: Dict[str, Dict[str, Dict]],
              stage_name: str,
              lig_path: str,
              score: Optional[float],
              heavy_atom_counts: Dict[str, int]) -> Optional[float]:
    ha = heavy_atom_counts.get(lig_path)
    le = compute_ligand_efficiency(score, ha)
    stage_map = score_history.setdefault(stage_name, {})
    rec = stage_map.setdefault(lig_path, {})
    rec["heavy_atoms"] = int(ha) if isinstance(ha, (int, float)) else None
    rec["le"] = le
    return le
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFMCS,  AllChem



def _read_any_lig(path: str):
    """
    Load ligand from SDF/MOL2/PDB with consistent settings.
    Returns an RDKit Mol or None.
    """
    mol = None
    loader = "unknown"
    ext = os.path.splitext(path)[1].lower()

    try:
        if ext in (".sdf", ".sd"):
            loader = "SDMolSupplier"
            suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
            mol = next((m for m in suppl if m is not None), None)
        elif ext in (".mol2",):
            loader = "MolFromMol2File"
            mol = Chem.MolFromMol2File(path, sanitize=True, removeHs=False)
        elif ext in (".pdb",):
            loader = "MolFromPDBFile"
            # If you use proximityBonding or flavor flags elsewhere, keep them consistent here.
            mol = Chem.MolFromPDBFile(path, sanitize=True, removeHs=False, proximityBonding=True)
        else:
            loader = "auto"
            mol = Chem.MolFromMolFile(path, sanitize=True, removeHs=False)  # last-ditch; or return None
    except Exception as e:
        logging.getLogger("rmsd").info(f"[read_any] loader={loader} path='{path}' load_failed={e}")
        mol = None

    # ──  single debug line about what we actually loaded ───────────────────
    try:
        _rlog = logging.getLogger("rmsd")
        if _rlog and mol is not None:
            from rdkit.Chem import rdMolDescriptors
            # formula = e.g., "C20H25N3O"
            formula = rdMolDescriptors.CalcMolFormula(mol)
            # InChIKey may be unavailable if RDKit was built without InChI; guard it.
            try:
                from rdkit.Chem import inchi
                inchikey = inchi.MolToInchiKey(mol)
            except Exception:
                inchikey = "NA"
            _rlog.info(f"[read_any] loader={loader} path='{path}' atoms={mol.GetNumAtoms()} "
                       f"heavy={mol.GetNumHeavyAtoms()} formula={formula} inchikey={inchikey}")
        elif _rlog:
            _rlog.info(f"[read_any] loader={loader} path='{path}' mol=None")
    except Exception:
        pass
    # ───────────────────────────────────────────────────────────────────────────

    return mol


def _is_readable_ref(pth: Path) -> bool:
    try:
        m = _read_any_lig(str(pth))
        return (m is not None) and (m.GetNumHeavyAtoms() > 0)
    except Exception:
        return False

def compute_rmsd(ref_path: str, docked_path: str) -> float:
    """Heavy-atom RMSD using best mapping; supports PDB/SDF/MOL2 refs and adds a minimal MCS fallback."""
    ref = _read_any_lig(ref_path)
    dock = _read_any_lig(docked_path)
    if not ref or not dock:
        return float("inf")

    # 1) Fast path: RDKit best alignment
    try:
        return float(rdMolAlign.GetBestRMS(ref, dock))
    except Exception:
        pass

    # 2) Tiny, robust fallback via MCS
    try:
        mcs = rdFMCS.FindMCS([ref, dock],
                             ringMatchesRingOnly=True,
                             completeRingsOnly=True,
                             matchValences=True)
        patt = Chem.MolFromSmarts(mcs.smartsString)
        if patt is None:
            return float("inf")
        ref_match = ref.GetSubstructMatch(patt)
        dock_match = dock.GetSubstructMatch(patt)
        if not ref_match or not dock_match or (len(ref_match) != len(dock_match)):
            return float("inf")
        amap = list(zip(dock_match, ref_match))  # (probe->ref)
        return float(rdMolAlign.AlignMol(dock, ref, atomMap=amap))
    except Exception:
        return float("inf")


def validate_ligand(
        ligand_name: str,
        docked_path: str,
        crystal_path: str = None,
        rmsd_thresh: float = 2.0,
        self_rmsd: float = None,
        logger=None
) -> bool:
    """
    Validate ligand docking.
      * If crystal structure available ? use redocking RMSD (hard gate).
      * Otherwise (non-controls) ? self-RMSD is *log-only* (never reject).
    """
    if crystal_path and Path(crystal_path).exists():
        redock_rmsd = compute_rmsd(crystal_path, docked_path)
        if logger:
            sr = f"{self_rmsd:.2f}" if isinstance(self_rmsd, (int, float)) else "n/a"
            logger.info(f"[validate] {ligand_name}: redock_RMSD={redock_rmsd:.2f} A, self_RMSD={sr}")
        if redock_rmsd <= rmsd_thresh:
            return True
        else:
            if logger:
                logger.warning(
                    f"[validate] {ligand_name}: redocking failed (RMSD {redock_rmsd:.2f} A > {rmsd_thresh:.2f})"
                )
            return False

    # Non-controls: log self-RMSD but do not gate on it
    try:
        sr_val = float(self_rmsd) if self_rmsd is not None else None
    except Exception:
        sr_val = None
    if logger:
        sr_txt = f"{sr_val:.2f}" if isinstance(sr_val, (int, float)) else "n/a"
        logger.info(f"[validate] {ligand_name}: self_RMSD={sr_txt} A (LOG-ONLY)")
    return True



# ======================
# Per-protein driver
# ======================
def _summarize_ions_file(file_path: Path | str) -> dict[str, object]:
    path = Path(file_path)
    if not path.exists():
        return {
            "hist": "missing",
            "counts": {},
            "metals_present": False,
            "salts_present": False,
            "error": "missing",
        }
    counts = Counter()
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith("HETATM"):
                    continue
                res = ln[17:20].strip().upper()
                elem = (ln[76:78].strip() or res).upper()
                token = elem if elem.isalpha() and 1 <= len(elem) <= 2 else res
                if token and token.isalpha() and len(token) <= 3:
                    counts[token] += 1
    except Exception as exc:  # pragma: no cover - diagnostics
        return {
            "hist": "error",
            "counts": {},
            "metals_present": False,
            "salts_present": False,
            "error": str(exc),
        }

    hist = ",".join(f"{tok}:{counts[tok]}" for tok in sorted(counts)) if counts else "none"
    metals_present = any(token in _ION_AUDIT_METALS and counts[token] > 0 for token in counts)
    salts_present = any(token in _ION_AUDIT_SALTS and counts[token] > 0 for token in counts)
    return {
        "hist": hist,
        "counts": dict(counts),
        "metals_present": metals_present,
        "salts_present": salts_present,
        "error": None,
    }


def process_one_protein(cfg: Dict, pdb_file: str, stages: List[Dict], params: RecenterParams) -> None:
    # ======================
    # Phase 0 – ID & path setup
    # ======================
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = re.sub(r'(?i)_cleaned$', '', base_id)
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=os.path.basename(pdb_file))

    logger = make_protein_logger(str(paths.docked_pdb_root()), pdb_id, cfg)
    logger.info(f"[paths] base_id={base_id} -> pdb_id={pdb_id}")
    logger.info(f"Processing protein: {pdb_file} (id={pdb_id})")
    # --- Canonicalize any runaway '.sanitized' filenames before we touch them ---
    lig_raw_dir  = os.path.join(cfg['OUTPUT_DIR'], pdb_id, 'ligands_raw')
    prepped_dir  = os.path.join(cfg['PREPPED_LIGANDS_DIR'], pdb_id)
    collapse_sanitized_names([lig_raw_dir, prepped_dir], logger=logger)

    # ======================
    # Phase 1 – Variant & ion context
    # ======================
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    variant_label = variant_env or "legacy"
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    active_ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None

    ion_audit_root: dict = cfg.setdefault("_ION_AUDIT", {})
    pdb_audit: dict = ion_audit_root.setdefault(paths.pdb_id, {})
    clean_audit: dict = pdb_audit.setdefault("clean_counts", {})

    input_summary = _summarize_ions_file(paths.input_pdb_path)
    input_hist = str(input_summary.get("hist", "none"))
    input_error = input_summary.get("error")
    if input_error not in (None, "missing"):
        logger.warning(
            "[ions.input.counts] pdb=%s file=%s action=skip err=%s",
            paths.pdb_id,
            paths.input_pdb_path,
            input_error,
        )
    else:
        logger.info(
            "[ions.input.counts] pdb=%s file=%s present_pdb=%s",
            paths.pdb_id,
            paths.input_pdb_path,
            input_hist,
        )
    pdb_audit["input_counts"] = {
        "hist": input_hist,
        "counts": dict(input_summary.get("counts", {})),
        "metals_present": bool(input_summary.get("metals_present", False)),
        "salts_present": bool(input_summary.get("salts_present", False)),
        "file": str(paths.input_pdb_path),
        "error": input_error,
    }

    cfg["_CURRENT_VARIANT"] = variant_env
    cleaned_target = paths.receptor_cleaned_pdb(variant_token)
    receptor_target = paths.receptor_pdbqt(variant_token, None)
    logger.info(
        "[receptor.path] pdb=%s variant=%s cleaned_pdb=%s exists=%s",
        paths.pdb_id,
        variant_label,
        cleaned_target,
        cleaned_target.exists(),
    )
    logger.info(
        "[receptor.path] pdb=%s variant=%s receptor_pdbqt=%s exists=%s",
        paths.pdb_id,
        variant_label,
        receptor_target,
        receptor_target.exists(),
    )

    # ======================
    # Phase 2 – Receptor prep (with ion summary)
    # ======================

    # 2) Protein prep (re-use if cached)
    logger.info("[ph.debug] calling prepare_receptor; PH_ENSEMBLE=%s", cfg.get("PH_ENSEMBLE", False))
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    provenance = getattr(prepare_receptor, "last_provenance", None)
    if provenance is None:
        try:
            provenance = getattr(protein_prep, "get_clean_provenance", lambda: "unknown")()
        except Exception:
            provenance = "unknown"
    logger.info("[receptor.clean.provenance] created_by=%s", provenance)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("Skipping protein due to prep failure.")
        return

    cleaned_hist = "none"
    if cleaned_pdb:
        clean_summary = _summarize_ions_file(cleaned_pdb)
        cleaned_hist = str(clean_summary.get("hist", "none"))
        clean_error = clean_summary.get("error")
        if clean_error not in (None, "missing"):
            logger.warning(
                "[ions.clean.counts] pdb=%s variant=%s action=skip err=%s",
                paths.pdb_id,
                variant_label,
                clean_error,
            )
        else:
            logger.info(
                "[ions.clean.counts] pdb=%s variant=%s file=%s present_pdb=%s",
                paths.pdb_id,
                variant_label,
                cleaned_pdb,
                cleaned_hist,
            )
        clean_audit[variant_label] = {
            "hist": cleaned_hist,
            "counts": dict(clean_summary.get("counts", {})),
            "metals_present": bool(clean_summary.get("metals_present", False)),
            "salts_present": bool(clean_summary.get("salts_present", False)),
            "file": str(cleaned_pdb),
            "error": clean_error,
        }
    else:
        clean_audit[variant_label] = {
            "hist": "missing",
            "counts": {},
            "metals_present": False,
            "salts_present": False,
            "file": "",
            "error": "missing",
        }

    try:
        probe_map = protein_prep.get_ion_probe_map(paths.pdb_id)  # type: ignore[attr-defined]
    except Exception as exc:
        logger.warning("[ions.summary] pdb=%s variant=%s action=skip err=%s", paths.pdb_id, variant_label, exc)
    else:
        before_counts = probe_map.get("strip_nonstandard:before", {})
        final_counts = probe_map.get("receptor_write", {})
        metals_before = sum(before_counts.values())
        metals_kept = sum(final_counts.values())
        metals_stripped = max(0, metals_before - metals_kept)
        logger.info("[ions.summary] pdb=%s variant=%s metals_kept=%d metals_stripped=%d", paths.pdb_id, variant_label, metals_kept, metals_stripped)

    if cleaned_pdb:
        if variant_env == "HOLO":
            skip_reason = "variant"
        elif not variant_env:
            skip_reason = "legacy"
        else:
            skip_reason = "disabled"
        logger.info(
            "[ions.prep-early] pdb=%s variant=%s action=skip reason=%s file=%s",
            paths.pdb_id,
            variant_label,
            skip_reason,
            cleaned_pdb,
        )

    # ======================
    # Phase 3 – Crystallographic ligand extraction / control setup
    # ======================
    _lig_count, control_stems = extract_ligands_to_nolig(paths, logger)
    try:
        from prep_ligands import prep_ligands_from_pdb
        prep_ligands_from_pdb(
            ligand_output_dir=paths.ligand_output_dir,
            ligands_mol2_dir=paths.ligands_mol2_dir,
            prepped_ligands_dir=paths.prepped_ligands_dir,
        )
        logger.info("[Controls] Prepped extracted controls ahead of redock.")
    except Exception as e:
        logger.warning(f"[Controls] Prepping extracted controls failed: {e}")

    ctrl_pdbqts: list[Path] = []
    for root in {paths.prepped_ligands_dir, Path(cfg["OUTPUT_LIGANDS_DIR"])}:
        if root.exists():
            ctrl_pdbqts.extend(root.glob("*.pdbqt"))

    logger.info(f"[Controls] Prepped control PDBQTs found (union): {len(ctrl_pdbqts)}")
    for p in ctrl_pdbqts[:10]:
        logger.info(f"[Controls]   {p.name}")

    control_lookup = build_control_lookup(paths)

    # ======================
    # Phase 4 – Pocket detection and initial center/box
    # ======================
    # 3) Pocket detection
    center, box_size, center_source = None, None, "none"
    try:
        sel_center, sel_box = select_center_via_control_redock(
            cfg,
            paths,
            receptor_pdbqt,
            logger,
            variant=variant_token,
            ph_token=active_ph_label,
            legacy=legacy_mode,
        )
    except Exception as _e:
        sel_center, sel_box = (None, None)
        logger.debug(f"[control-centers] helper errored: {_e}")
    if sel_center is not None:
        center, box_size, center_source = sel_center, sel_box, "control"
        logger.info(f"[control-redock] Using control-derived center {center} with box {box_size}")
    else:
        # P2Rank last resort (controls absent or all redocks failed)
        c2, b2 = detect_active_site(cleaned_pdb)
        if c2:
            box_size = tuple(min(28.0, float(s)) for s in b2)
            center = c2
            center_source = "p2rank"
            logger.info(f"[P2Rank] Using P2Rank center {center} with box {box_size}")
        else:
            logger.error("Active-site detection failed (no usable controls, P2Rank returned None).")
            return
    if center is None:
        return

    try:
        import automate_protein_prep as _auto_prep_mod
    except Exception as ions_err:
        logger.warning("[ions] pocket_refine_skip err=%s", ions_err)
    else:
        if cleaned_pdb and variant_env == "HOLO":
            logger.info(
                "[ions.pocket-pass] pdb=%s variant=%s action=refine_with_center file=%s",
                paths.pdb_id,
                variant_label,
                cleaned_pdb,
            )
            try:
                _auto_prep_mod._maybe_strip_ions(
                    Path(cleaned_pdb),
                    cfg=cfg,
                    variant=variant_token,
                    pocket_center=center,
                )
            except Exception as pocket_err:
                logger.warning(
                    "[ions] pocket_refine_skip err=%s",
                    pocket_err,
                )

    # Override control-box size from config (keeps existing 24 A default)
    if center_source == "control":
        side = float(cfg.get("CONTROL_BOX_A", 24.0))
        box_size = (side, side, side)

    # clamp initial box once to keep Vina happy (detect_pocket already caps P2Rank path)
    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
    box_size = tuple(min(box_cap, float(s)) for s in box_size)
    logger.info(f"Initial box clamped to {box_size} (cap={box_cap} A)")

    # explicit console breadcrumb so you don't need to open logs
    try:
        c_print = tuple(round(float(x), 3) for x in center)
        b_print = tuple(round(float(x), 1) for x in box_size)
        print(f"[CENTER] source={center_source} center={c_print} box={b_print}")
    except Exception:
        pass

    # HOLO-only restore of metals/cofactors, now that center/box_size are known
    metals_added = 0
    cofactors_added = 0
    regen = False
    try:
        logger.info("[holo.restore.call] invoking for pdb=%s", paths.pdb_id)
        metals_added, cofactors_added, regen = protein_prep._holo_restore_from_input_if_needed(
            pdb_id=paths.pdb_id,
            cleaned_pdb=cleaned_pdb,
            output_pdbqt=receptor_pdbqt,
            config=cfg,
            center=center,
            box_size=box_size,
        )

    except Exception as _restore_err:
        logger.warning("[holo.restore] action=skip reason=%s", _restore_err)

    if variant_env == "HOLO":
        logger.info(
            "[holo.restore.summary] pdb=%s variant=%s metals_added=%d cofactors_added=%d regen=%s",
            paths.pdb_id,
            variant_env,
            metals_added,
            cofactors_added,
            regen,
        )

        if regen:
            logger.warning(
                "[holo.restore.regen] pdb=%s variant=%s regen=True; rebuilding receptor PDBQT from %s -> %s",
                paths.pdb_id,
                variant_env,
                cleaned_pdb,
                receptor_target,
            )
            try:
                ok_after = protein_prep.run_prepare_receptor(
                    input_pdb=cleaned_pdb,
                    output_pdbqt=str(receptor_target),
                    cfg=cfg,
                )
                if not ok_after or not receptor_target.exists():
                    logger.warning(
                        "[holo.restore.regen] status=failed pdb=%s; keeping previous receptor PDBQT=%s",
                        paths.pdb_id,
                        receptor_pdbqt,
                    )
                else:
                    receptor_pdbqt = str(receptor_target)
                    logger.info(
                        "[holo.restore.regen] status=ok pdb=%s receptor_pdbqt=%s",
                        paths.pdb_id,
                        receptor_pdbqt,
                    )
            except Exception as regen_err:
                logger.warning(
                    "[holo.restore.regen] status=error pdb=%s err=%s; keeping previous receptor PDBQT=%s",
                    paths.pdb_id,
                    regen_err,
                    receptor_pdbqt,
                )
        else:
            logger.info(
                "[holo.restore.regen] pdb=%s variant=%s regen=False; skipping receptor PDBQT rebuild",
                paths.pdb_id,
                variant_env,
            )

        if receptor_pdbqt:
            try:
                protein_prep.run_metal_site_audit(
                    pdb_id=paths.pdb_id,
                    router_paths=paths,
                    input_pdb_path=str(paths.input_pdb_path),
                    receptor_pdb_path=cleaned_pdb,
                    receptor_pdbqt_path=receptor_pdbqt,
                    center=center,
                    variant_label=variant_label,
                    ph_label=active_ph_label,
                )
            except Exception as audit_err:
                logger.warning(
                    "[holo.metal_audit] action=skip pdb=%s reason=%s",
                    paths.pdb_id,
                    audit_err,
                )


    # Preflight HOLO skip: avoid redundant HOLO work when receptors are byte-identical to APO
    resolved_mode = (str(cfg.get("_RESOLVED_APO_HOLO_MODE")) or "").strip().lower() or "legacy"
    if variant_env == "HOLO" and resolved_mode == "apo_vs_holo":
        apo_clean = _variant_receptor_path(pdb_id, "APO", cfg)
        holo_clean = cleaned_pdb or _variant_receptor_path(pdb_id, "HOLO", cfg)
        apo_path = Path(apo_clean) if apo_clean else None
        holo_path = Path(holo_clean) if holo_clean else None
        apo_exists = apo_path.exists() if apo_path else False
        holo_exists = holo_path.exists() if holo_path else False

        if not apo_exists or not holo_exists:
            logger.warning(
                "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=continue reason=missing_paths apo=%s holo=%s",
                pdb_id,
                apo_clean,
                holo_clean,
            )
            _record_apo_holo_decision(cfg, pdb_id, "HOLO", "missing_paths")
        else:
            # >>> path+exists breadcrumb just before SHA calculation <<<
            logger.info(
                "[apo-vs-holo] compare.preflight apo=%s exists=%s holo=%s exists=%s",
                norm(apo_path), ("T" if apo_exists else "F"),
                norm(holo_path), ("T" if holo_exists else "F"),
            )
            try:
                apo_sha = file_sha1(str(apo_path))
                holo_sha = file_sha1(str(holo_path))
            except Exception as hash_err:
                logger.warning(
                    "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=skip reason=sha_error err=%s",
                    pdb_id,
                    hash_err,
                )
                _record_apo_holo_decision(cfg, pdb_id, "HOLO", "sha_error")
            else:
                if apo_sha == holo_sha:
                    logger.info(
                        "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=skip reason=identical apo_sha=%s holo_sha=%s",
                        pdb_id,
                        apo_sha,
                        holo_sha,
                    )
                    # [ions] dedup audit guard
                    audit_root = cfg.get("_ION_AUDIT", {})
                    pdb_entry = audit_root.get(paths.pdb_id) or audit_root.get(pdb_id)
                    warn_needed = False
                    if isinstance(pdb_entry, dict):
                        input_info = pdb_entry.get("input_counts", {})
                        clean_map = pdb_entry.get("clean_counts", {}) or {}
                        holo_info = clean_map.get("HOLO") or clean_map.get(variant_label) or {}
                        if input_info.get("metals_present") or input_info.get("salts_present"):
                            warn_needed = True
                        if holo_info.get("metals_present") or holo_info.get("salts_present"):
                            warn_needed = True
                    if warn_needed:
                        logger.warning(
                            "[apo-vs-holo] unexpected_identical_after_ion_policy pdb=%s apo_sha=%s holo_sha=%s",
                            pdb_id,
                            apo_sha,
                            holo_sha,
                        )
                    if receptor_pdbqt:
                        _record_apo_holo_usage(cfg, pdb_id, variant_token, None, receptor_pdbqt)
                    _record_apo_holo_decision(cfg, pdb_id, "HOLO", "skipped_preflight")
                    try:
                        delete_variant_trees(pdb_id, "HOLO", cfg)
                    except Exception as cleanup_err:
                        logger.warning(
                            "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=cleanup_warn err=%s",
                            pdb_id,
                            cleanup_err,
                        )
                    return
                else:
                    logger.info(
                        "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=continue reason=not_identical apo_sha=%s holo_sha=%s",
                        pdb_id,
                        apo_sha,
                        holo_sha,
                    )
                    _record_apo_holo_decision(cfg, pdb_id, "HOLO", "not_identical")
                    # (regen-aware rebuild happens earlier when HOLO restore requests it)




    # ======================
    # Phase 5 – pH ensemble manifest (global protein-level)
    # ======================
    # >>> PH ENSEMBLE (GLOBAL) START
    if bool(cfg.get("PH_ENSEMBLE", False)):
        try:
            from context_ph import select_ph_values_for_protonation
            from ph_ensemble import build_ph_ensemble

            # Pull raw list from context_ph on the **raw input PDB** (header intact)
            raw_vals = select_ph_values_for_protonation(str(paths.input_pdb_path))
            logger.info("[ph.ctx.raw] path=%s values=%s", str(paths.input_pdb_path),
                        ",".join(f"{v:.2f}" for v in (raw_vals or [])))

            # Round to 0.1 and clamp to [3.0, 10.5]; dedupe + sort
            ph_values = sorted({max(3.0, min(10.5, round(float(x), 1))) for x in (raw_vals or [])})
            if not ph_values:
                logger.warning("[ph.ctx.fallback] context list empty -> using [7.0]")
                ph_values = [7.0]

            logger.info("[ph.list] n=%d values=%s", len(ph_values),
                        ",".join(f"{v:.1f}" for v in ph_values))

            # GLOBAL scope: use the propka_wire sentinel (radius >= 1e6)
            manifest_path = build_ph_ensemble(
                pdb_id=paths.pdb_id,
                cleaned_receptor_pdb=str(Path(cleaned_pdb)),
                out_dir=str(Path(cfg["OUTPUT_DIR"])),
                center=(0.0, 0.0, 0.0),
                radius=1_000_000.0,
                ph_values=ph_values,
                member_index_start=0,
                variant=variant_token,
                legacy=legacy_mode,
            )
            logger.info("[ph_ensemble.manifest] path=%s", manifest_path)

        except Exception as e:
            logger.warning("[ph_ensemble.skip] error=%s", e)
    # >>> PH ENSEMBLE (GLOBAL) END

    # ======================
    # Phase 6 – Ligand prep & filtering
    # ======================
    # 4) Ligand prep & filtering
    cfg.setdefault("_EFFECTIVE_SINGLE_LIGAND", "")
    single_ligand_hit: Optional[Path] = None
    cfg.pop("_SINGLE_RESOLVED_PATH", None)
    # --- Single-ligand mode (if active) --------------------------------------
    if cfg["_EFFECTIVE_SINGLE_LIGAND"]:
        _ensure_single_ligand_index(cfg, paths, logger)
        # Provide per-protein paths to resolver
        cfg.setdefault("paths", {})
        cfg["paths"]["prepped_ligands_dir"] = str(paths.prepped_ligands_dir)

        hit = _resolve_single_ligand(cfg["_EFFECTIVE_SINGLE_LIGAND"], pdb_id, cfg, logger)
        if hit:
            cfg["_SINGLE_RESOLVED_PATH"] = str(hit)
            single_ligand_hit = hit
        else:
            selector_token = cfg["_EFFECTIVE_SINGLE_LIGAND"]
            suggestions: list[str] = []
            try:
                import difflib

                fda_map = _load_fda_name_map(cfg, logger)
                suggestions = difflib.get_close_matches(
                    selector_token,
                    list(fda_map.keys()),
                    n=5,
                    cutoff=0.7,
                )
            except Exception:
                suggestions = []
            if suggestions:
                logger.error("[single.miss.suggest] did_you_mean=%s", ", ".join(suggestions))

            allow_flag = os.environ.get("ALLOW_FDA_FALLBACK")
            if allow_flag is None:
                allow_flag = cfg.get("ALLOW_FDA_FALLBACK", False)
            if not _to_bool(allow_flag):
                logger.error(
                    "[single.block] selector '%s' not found in fda_library via FDA_MAPPING_CSV; aborting instead of fallback.",
                    selector_token,
                )
                raise SystemExit(2)
            logger.warning(
                "[single.block] selector '%s' not found; ALLOW_FDA_FALLBACK enabled, continuing with fallback flow.",
                selector_token,
            )
            cfg["_EFFECTIVE_SINGLE_LIGAND"] = ""

    if cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
        if not single_ligand_hit:
            return
        ligands = [str(single_ligand_hit)]
        ha = _count_heavy_atoms_from_pdbqt(single_ligand_hit)
        heavy_atom_counts = {str(single_ligand_hit): ha}
        pains_flags = {}
        logger.info(f"[single] Active ? docking only: {single_ligand_hit.name} (heavy={ha})")

        # --- PH-ligand support for single-ligand mode ---
        if cfg.get("PH_LIGAND_MODE", "").lower() == "context_window" and cfg.get("PH_ENSEMBLE_IN_PREP"):
            try:
                from prep_ligands import enumerate_ligands_for_docking

                ph_values = [6.0, 8.0]
                if "_PH_CONTEXT_VALUES" in cfg:
                    ph_values = cfg["_PH_CONTEXT_VALUES"]
                ligand_window = sorted(
                    {round(p, 1) for ph in ph_values for p in (float(ph) - 1.0, float(ph), float(ph) + 1.0)}
                )
                logger.info(f"[single.ph_ligand] Using ligand window {ligand_window}")

                ph_root_cfg = cfg.get("_PH_LIGAND_ROOT", "")
                try:
                    ph_root_path = Path(ph_root_cfg) if ph_root_cfg else None
                except Exception:
                    ph_root_path = None

                logger.info(
                    "[single.ph_ligand.bridge] ph_root_cfg=%s ph_root_path=%s exists=%s",
                    ph_root_cfg,
                    str(ph_root_path) if ph_root_path is not None else "",
                    ph_root_path.exists() if ph_root_path is not None else False,
                )
                if ph_root_path is not None and ph_root_path.exists():
                    enumerate_ligands_for_docking(
                        requested_ph_values=ligand_window,
                        root_dir=ph_root_path,
                        microstate_dedup=True,
                        force=False,
                    )
                else:
                    logger.info(
                        "[single.ph_ligand.bridge.skip] no valid _PH_LIGAND_ROOT; "
                        "skipping microstate priming for single-ligand mode"
                    )
            except Exception as e:
                logger.warning(f"[single.ph_ligand.skip] Could not run PH-ligand window for single mode: {e}")
        # ------------------------------------------------
    else:
        ligands, heavy_atom_counts, pains_flags = prepare_and_filter_ligands(cfg, paths, logger)

    # (skipped in single-ligand mode)
    if not cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
        # Force-inject control PDBQTs if they exist on disk but weren't selected
        ctrl_stems_lower = {s.lower() for s in control_stems}

        prepped_control_pdbqts = []
        # Re-scan now that prep_ligands_from_pdb has run
        scan_roots = [paths.prepped_ligands_dir]
        if cfg.get("OUTPUT_LIGANDS_DIR"):
            try:
                out_root = Path(cfg["OUTPUT_LIGANDS_DIR"])
                if out_root.exists():
                    scan_roots.append(out_root)
            except Exception:
                pass

        for root in scan_roots:
            if root and root.exists():
                for p in root.glob("*.pdbqt"):
                    stem0 = p.stem.split("_stage")[0].lower()
                    if stem0 in ctrl_stems_lower:
                        prepped_control_pdbqts.append(p)

        lig_set = {norm(x) for x in ligands}
        missing_controls = [p for p in prepped_control_pdbqts if norm(p) not in lig_set]

        if missing_controls:
            logger.info(f"[Controls] Adding {len(missing_controls)} prepared control(s) to Stage1.")
            # Front-load controls
            ligands = [str(p) for p in missing_controls] + ligands
            for p in missing_controls:
                try:
                    heavy_atom_counts.setdefault(str(p), _count_heavy_atoms_from_pdbqt(p))
                except Exception:
                    heavy_atom_counts.setdefault(str(p), 0)



    # --- Normalize & de-dupe Stage1 ligand list (keep order) ---
    def _norm_dedupe(seq):
        seen = set()
        out = []
        for p in seq:
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                out.append(pn)
        return out

    ligands = _norm_dedupe(ligands)

    base_ligands = ligands[:]
    base_heavy_atoms = dict(heavy_atom_counts)
    base_pains_flags = dict(pains_flags)
    base_center = tuple(center)
    base_box = tuple(box_size)

    # ======================
    # Phase 7 – pH/variant loop (core docking)
    # ======================
    ph_log = logging.getLogger("ph_ensemble")
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))
    plan_only = os.environ.get("A2_PLAN_ONLY") == "1"

    ph_tags = init_ph_tags_and_manifest(cfg, paths.pdb_id, variant_token, legacy_mode)
    if ph_enabled and not ph_tags:
        # keep the early return behavior for empty ensembles
        return

    ph_ligand_root = None
    ph_ligand_root_str = cfg.get("_PH_LIGAND_ROOT")
    if ph_ligand_root_str:
        try:
            ph_ligand_root = Path(ph_ligand_root_str)
        except Exception:
            ph_ligand_root = None

    prewarm_ph_ligand_microstates(cfg, ph_tags, ph_ligand_root)

    for ph_label in ph_tags:
        rec_path = receptor_file(paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode)
        out_root = docked_dir(paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode)
        ph_print = ph_label or "(none)"
        rec_exists = rec_path.exists()
        logger.info(
            "[router] pdb=%s variant=%s ph=%s receptor_file=%s docked_dir=%s exists=%s",
            paths.pdb_id,
            variant_label,
            ph_print,
            str(rec_path),
            str(out_root),
            rec_exists,
        )

        if plan_only:
            print(
                f"pdb={paths.pdb_id} variant={variant_label} ph={ph_print} "
                f"receptor_file={rec_path} docked_dir={out_root} exists={rec_exists}"
            )
            continue

        _record_apo_holo_usage(cfg, paths.pdb_id, variant_token, ph_label, rec_path)

        if not rec_exists:
            ph_log.warning(
                "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s receptor_missing=%s",
                paths.pdb_id,
                variant_label,
                ph_print,
                str(rec_path),
            )
            continue

        if ph_label:
            cfg["_ACTIVE_PH_LABEL"] = ph_label
        else:
            cfg.pop("_ACTIVE_PH_LABEL", None)

        ligands = base_ligands[:]
        heavy_atom_counts = dict(base_heavy_atoms)
        pains_flags = dict(base_pains_flags)
        center = tuple(base_center)
        box_size = tuple(base_box)
        receptor_pdbqt = str(rec_path)

        enumerated = enumerate_ligands_for_ph_context(
            cfg=cfg,
            pdb_id=paths.pdb_id,
            ph_label=ph_label,
            ph_ligand_root=ph_ligand_root,
        )

        if enumerated:
            ligands = [str(p) for p in enumerated]
            heavy_atom_counts = {
                str(p): _count_heavy_atoms_from_pdbqt(p) for p in enumerated
            }
            pains_flags = {
                k: base_pains_flags.get(
                    k,
                    base_pains_flags.get(Path(k).stem, False),
                )
                for k in ligands
            }
        else:
            # keep ligands, heavy_atom_counts, pains_flags at their base values
            ligands = base_ligands[:]
            heavy_atom_counts = dict(base_heavy_atoms)
            pains_flags = dict(base_pains_flags)
        # --------------------------------------------------------


        ctrl_stems_lower = {s.lower() for s in control_stems}
        ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
        min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

        def _is_control_path(p: str) -> bool:
            stem = Path(p).stem.split("_stage")[0]
            if stem.upper() in ctrl_blacklist:
                return False
            if stem.lower() not in ctrl_stems_lower:
                return False
            ha = heavy_atom_counts.get(p)
            return (ha is None) or (ha >= min_ha)

        ctrls = [p for p in ligands if _is_control_path(p)]
        non_ctrls = [p for p in ligands if not _is_control_path(p)]
        if ctrls:
            ligands = ctrls + non_ctrls
            logger.info(
                f"[Controls] Front-loading {len(ctrls)} controls. "
                f"First wave: {[Path(x).name for x in ligands[:int(cfg.get('MAX_PARALLEL_JOBS', 1))]]}"
            )

        present_ctrls = [
            Path(l).stem.split("_stage")[0].lower()
            for l in ligands
            if Path(l).stem.split("_stage")[0].lower() in ctrl_stems_lower
        ]

        if not present_ctrls:
            logger.warning(
                "[Controls] No control ligands present in Stage1 ligand list -- "
                "self-RMSD/locking will not be possible. (Check prep errors above.)"
            )
        if not ligands:
            logger.warning("No valid ligands after filtering; skipping protein.")
            continue

        selector = CenterSelector(cfg, logger, control_stems, heavy_atom_counts, center)
        guard = GlobalCenterGuard(
            max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2))
        )

        stage1_original = ligands[:]

        control_stems_lower = {s.lower() for s in control_stems}
        forced_extracted_for_stage3 = {
            lig for lig in stage1_original
            if Path(lig).stem.split("_stage")[0].lower() in control_stems_lower
        }
        logger.info(
            f"[Force-carry] Extracted ligands earmarked for Stage3: {len(forced_extracted_for_stage3)}"
        )

        score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
        validated_ligands_last: List[str] = []
        recenter_attempts = 0
        docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()

        retry_mgr = RetryManager()

        i = 0
        while i < len(stages):
            guard.reset_stage()
            stage = stages[i]

            if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                if checkpoint_should_skip(cfg, paths.pdb_id, stage["name"], fp):
                    logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                    i += 1
                    continue

            if not ligands:
                logger.warning(f"No ligands to dock at {stage['name']}; stopping for this protein.")
                break

            logger.info(f"Starting {stage['name']} with {len(ligands)} ligands...")
            if ph_label:
                stage_dir = paths.docked_stage_dir(variant_env or None, stage["name"], ph_label)
                ph_log.info(
                    "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s stage=%s receptor=%s out=%s",
                    paths.pdb_id,
                    variant_label,
                    ph_label,
                    stage["name"],
                    receptor_pdbqt,
                    str(stage_dir),
                )

            if i == 0 and ctrls and non_ctrls:
                logger.info(
                    f"Stage1 two-wave: {len(ctrls)} controls first, then {len(non_ctrls)} others."
                )

                s1, v1, d1, rd1, inv1 = run_one_stage(
                    cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                    ctrls, logger, retry_mgr, control_lookup
                )

                try:
                    dec = selector.consider_switch(stage['name'], s1, v1, rd1, receptor_pdbqt, center, guard)
                    if dec.promoted and dec.new_center is not None:
                        old = center
                        center = dec.new_center
                        guard.mark_switch()
                        logger.info(
                            f"[CENTER] Switched before library run: {old} -> {center} ({dec.reason}) [global switch]"
                        )
                except Exception as e:
                    logger.warning(f"CenterSelector (controls-only) failed gracefully: {e}")

                lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", -6.0))
                lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))
                qualified_controls = []

                for lig in v1:
                    stem = Path(lig).stem.split("_stage")[0].lower()
                    ha = heavy_atom_counts.get(lig)
                    if stem in control_stems_lower and (ha is None or ha >= min_ha):
                        sc = s1.get(lig)
                        if sc is not None and np.isfinite(sc) and sc <= lock_score_max:
                            pose_path = rd1.get(lig)
                            c = CenterSelector._pdbqt_centroid(pose_path) if pose_path else None
                            if c is not None and np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                                qualified_controls.append(lig)

                if len(qualified_controls) >= lock_min_hits and not guard.locked:
                    guard.lock()
                    logger.info(
                        "[CONTROL-LOCK] Early lock from controls-only wave "
                        f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} A); "
                        "future center switches disabled."
                    )

                s2, v2, d2, rd2, inv2 = run_one_stage(
                    cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                    non_ctrls, logger, retry_mgr, control_lookup
                )

                scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
                raw_docked = {**rd1, **rd2}
                invalids = {**inv1, **inv2}
            else:
                scores, validated, distances, raw_docked, invalids = run_one_stage(
                    cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                    ligands, logger, retry_mgr, control_lookup
                )

            validated_ligands_last = validated

            def _is_control(lig: str) -> bool:
                stem = Path(lig).stem.split("_stage")[0].lower()
                if stem.upper() in {s.strip().upper() for s in cfg.get("CONTROL_BLACKLIST", "").split(",") if s.strip()}:
                    return False
                ha = heavy_atom_counts.get(lig)
                if ha is not None and ha < int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10)):
                    return False
                return stem in {s.lower() for s in control_stems}

            control_anchor_hit = any(_is_control(lig) for lig in validated)

            lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", float("inf")))
            lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
            lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))

            qualified_controls = []
            for lig in validated:
                if not _is_control(lig):
                    continue
                sc = scores.get(lig)
                if sc is None or not np.isfinite(sc):
                    continue
                if sc > lock_score_max:
                    continue
                pose_path = raw_docked.get(lig)
                if not pose_path:
                    continue
                c = CenterSelector._pdbqt_centroid(pose_path)
                if c is None:
                    continue
                if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                    qualified_controls.append(lig)

            if len(qualified_controls) >= lock_min_hits and not guard.locked:
                guard.lock()
                logger.info(
                    "[CONTROL-LOCK] Control(s) validated with strong confidence "
                    f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Ang); "
                    "center is now anchored; future center switches are disabled."
                )

            try:
                processed = {norm(x) for x in ligands}
                valid_set = {norm(x) for x in scores.keys()}
                invalid_set = {norm(x) for x in invalids.keys()}
                both = valid_set & invalid_set
                missing = processed - (valid_set | invalid_set)
                if both or missing:
                    logger.error(
                        f"Invariant violation at {stage['name']}: both={len(both)}, missing={len(missing)}"
                    )
                    if both:
                        logger.error(
                            "Ligands marked both valid & invalid: "
                            + ", ".join(os.path.basename(x) for x in list(both)[:10])
                        )
                    if missing:
                        logger.error(
                            "Ligands missing from results: "
                            + ", ".join(os.path.basename(x) for x in list(missing)[:10])
                        )
            except Exception as _e:
                logger.warning(f"Invariant check failed: {_e}")

            for lig, sc in scores.items():
                record_score(score_history, stage['name'], lig, sc, True)
                record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)
            for lig, (sc, reason) in invalids.items():
                record_score(score_history, stage['name'], lig, sc, False, reason=reason)
                record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)

            promoted_this_stage = False
            try:
                decision = selector.consider_switch(
                    stage['name'], scores, validated, raw_docked, receptor_pdbqt, center, guard
                )
                if decision.promoted and decision.new_center is not None:
                    old = center
                    center = decision.new_center
                    promoted_this_stage = True
                    guard.mark_switch()
                    if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                        checkpoint_invalidate_from(cfg, paths.pdb_id, stages, start_index=i)
                    logger.info(
                        f"[CENTER] Switched from {old} -> {center} ({decision.reason}, "
                        f"SwitchScore={decision.switchscore:.2f}) [global switch]"
                    )
            except Exception as e:
                logger.warning(f"CenterSelector failed gracefully: {e}")

            if ph_label:
                ph_log.info(
                    "[ph_ensemble.dock.scores] pdb_id=%s variant=%s ph=%s stage=%s valid=%d invalid=%d",
                    paths.pdb_id,
                    variant_label,
                    ph_label,
                    stage["name"],
                    len(scores),
                    len(invalids),
                )

            if not promoted_this_stage:
                restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                    i, scores, distances, box_size, center, stage1_original, recenter_attempts, params,
                    cfg, paths.pdb_id, receptor_pdbqt, logger, raw_docked, guard, control_anchor_hit
                )
                if restart:
                    ligands = redo_ligands
                    if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                        checkpoint_invalidate_from(cfg, paths.pdb_id, stages, start_index=0)
                    i = 0
                    continue

            try:
                if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                    med = (
                        float(np.median([d for d in distances if isinstance(d, (int, float))]))
                        if distances
                        else None
                    )
                    if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                        dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                        min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                        new_box = tuple(max(min_box, s - dec) for s in box_size)
                        if new_box != box_size:
                            logger.info(
                                f"Adaptive shrink: median dist {med:.2f} A -> box {box_size} -> {new_box}"
                            )
                            box_size = new_box
            except Exception as _e:
                logger.warning(f"Adaptive shrink skipped: {_e}")

            if i < len(stages) - 1:
                if not scores:
                    restart, center, box_size, redo_ligands = fallback_recentering_if_empty(
                        cfg, paths.pdb_id, stage['name'], scores, raw_docked,
                        receptor_pdbqt, center, box_size, stage1_original, logger, guard, control_anchor_hit
                    )
                    if restart:
                        ligands = redo_ligands
                        if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                            checkpoint_invalidate_from(cfg, paths.pdb_id, stages, start_index=0)
                        i = 0
                        continue
                    else:
                        break

                use_stage1_base = (docking_mode == "polypharmacology" and i == 1)

                rescue = []
                if i < len(stages) - 1:
                    for lig, (sc, reason) in invalids.items():
                        if sc is not None and "self_rmsd_" in str(reason).lower() and sc <= float(
                                cfg.get("RESCUE_SELF_RMSD_SCORE_MAX", -8.0)):
                            rescue.append((sc, lig))
                    rescue = [lig for _, lig in sorted(rescue)[:int(cfg.get("RESCUE_SELF_RMSD_TOP_N", 10))]]

                selected = select_ligands_for_next(
                    docking_mode,
                    i,
                    stages,
                    scores,
                    logger,
                    base_pool_n=(len(stage1_original) if use_stage1_base else None),
                    force_include=(forced_extracted_for_stage3 if use_stage1_base else None)
                )

                if rescue:
                    sel_set = set(selected)
                    rescue_unique = [r for r in rescue if r not in sel_set]
                    ligands = rescue_unique + selected
                else:
                    ligands = selected
                if not ligands:
                    logger.warning(f"No ligands selected for {stages[i + 1]['name']}; stopping.")
                    break

            if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                try:
                    fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                    checkpoint_mark_done(cfg, paths.pdb_id, stage["name"], fp)
                except Exception:
                    pass

            i += 1

        # ======================
        # Phase 8 – Summary outputs & cleanup
        # ======================
        final_pose_validation_and_screenshots(
            cfg, paths.pdb_id, stages, receptor_pdbqt, center, validated_ligands_last,
            score_history, cleaned_pdb, docking_mode, logger, ph_label
        )

        csv_path = write_scores_csv(cfg, paths.pdb_id, score_history)
        logger.info(
            "[Scores] ph_label=%s summary=%s",
            ph_label if ph_label else "base",
            csv_path,
        )

        try:
            summary = {
                "pdb_id": paths.pdb_id,
                "center": tuple(map(float, center)) if center else None,
                "box_size": tuple(map(float, box_size)) if box_size else None,
                "n_ligands_stage1": len(stage1_original),
                "n_valid_last_stage": len(validated_ligands_last),
                "switch_history": getattr(selector, "switch_history", []),
                "global_switches": guard.global_switches,
                "stages": [s["name"] for s in stages],
                "ph_label": ph_label,
            }
            _write_audit_json(cfg, paths.pdb_id, summary)
        except Exception as _e:
            logger.warning(f"Audit JSON write failed: {_e}")

    cfg.pop("_ACTIVE_PH_LABEL", None)


def _smoke_emit_config_demo() -> None:
    """Emit a small config to exercise router paths in isolation."""
    smoke_log = logging.getLogger("smoke")
    old_variant = os.environ.get("APO_HOLO_VARIANT")
    try:
        base_cfg = ConfigDict(load_inputs())
    except Exception as exc:
        smoke_log.warning("[smoke.emit.skip] reason=%s", exc)
        return

    try:
        cfg = ConfigDict(base_cfg.copy())
        repo_root = Path(__file__).resolve().parent
        smoke_root = repo_root / "analysis" / "_smoke"
        overrides = {
            "OVERALL_DIR": smoke_root,
            "INPUT_DIR": smoke_root / "input_pdbs",
            "OUTPUT_DIR": smoke_root / "processed_pdbs",
            "DOCKED_DIR": smoke_root / "docked",
            "PREPPED_LIGANDS_DIR": smoke_root / "prepped_ligands",
            "OUTPUT_LIGANDS_DIR": smoke_root / "prepped_ligands",
            "PREPPED_LIGANDS_ROOT": smoke_root / "prepped_ligands",
            "LIGANDS_MOL2_DIR": smoke_root / "ligands_mol2",
            "CONFIGS_DIR": smoke_root / "configs",
        }
        for key, path_value in overrides.items():
            cfg[key] = str(path_value)
            Path(path_value).mkdir(parents=True, exist_ok=True)

        cfg["RUN_ID"] = "smoke_demo"
        cfg["RESET_CONFIGS"] = False
        init_config_run_dir(cfg, run_id=cfg["RUN_ID"], reset=False, logger=smoke_log)

        mode, variants = resolve_apo_holo_mode(cfg)
        smoke_log.info("[smoke.emit] mode=%s variants=%s", mode, variants)

        pdb_id = "3CS9"
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        lig_dir = paths.prepped_ligands_dir
        lig_dir.mkdir(parents=True, exist_ok=True)
        lig_path = lig_dir / "smoke_ligand.pdbqt"
        if not lig_path.exists():
            lig_path.write_text("SMOKE", encoding="utf-8")

        ph_label = "pH6_7"
        cfg["_ACTIVE_PH_LABEL"] = ph_label
        os.environ["APO_HOLO_VARIANT"] = "APO"

        receptor_path = receptor_file(paths.pdb_id, variant="APO", ph_tag=ph_label, legacy=False)
        receptor_path.parent.mkdir(parents=True, exist_ok=True)
        if not receptor_path.exists():
            receptor_path.write_text("RECEPTOR", encoding="utf-8")

        stage_info = {"name": "smoke_stage", "exhaustiveness": 8, "num_modes": 9, "verbosity": 0}
        conf_path, out_path = emit_vina_config(
            cfg,
            paths.pdb_id,
            str(receptor_path),
            (0.0, 0.0, 0.0),
            (20.0, 20.0, 20.0),
            str(lig_path),
            stage_info["name"],
            stage_info,
            1,
            smoke_log,
            variant="APO",
            ph_token=ph_label,
            legacy=False,
        )
        smoke_log.info("[smoke.emit.done] config=%s out=%s", conf_path, out_path)
    except Exception as exc:
        smoke_log.warning("[smoke.emit.skip] reason=%s", exc)
    finally:
        if old_variant is None:
            os.environ.pop("APO_HOLO_VARIANT", None)
        else:
            os.environ["APO_HOLO_VARIANT"] = old_variant


# ======================
# Program entry point
# ======================
def main() -> None:
    print("MODELLER is working with license.")
    run_id = _resolve_run_id(sys.argv)
    os.environ["ATLAS_RUN_ID"] = run_id
    log_path = _prepare_run_logfile(run_id)
    os.environ["ATLAS_LOG_FILE"] = log_path
    _tee_stdio_to(log_path)
    cfg = ConfigDict(load_inputs())
    bootstrap_root_logging(cfg, log_path)
    logging.info("[probe.root] root-logger INFO now visible")
    print(f"[run] log_file={log_path} run_id={run_id}")
    validate_config(cfg)

    rules = get_atom_rules()
    alias_sets = getattr(rules, "alias_sets", None)
    waters_set = set(getattr(rules, "waters", set()))
    if not waters_set and alias_sets is not None:
        waters_set = {str(tok).strip().upper() for tok in getattr(alias_sets, "waters", set()) if str(tok).strip()}
    cofactors_set = {str(tok).strip().upper() for tok in getattr(rules, "cofactors", set()) if str(tok).strip()}
    elements_set = {str(tok).strip().upper() for tok in getattr(rules, "elem_tokens_canonical", set()) if str(tok).strip()}
    logging.info(
        "[aliases.summary] mode=%s keep={waters:%d, cofactors:%d, elements:%d}",
        getattr(rules, "policy_mode", "LEGACY"),
        len(waters_set),
        len(cofactors_set),
        len(elements_set),
    )

    cfg.setdefault("PH_ENSEMBLE", False)
    cfg.setdefault("PH_RADIUS", 10.0)

    env_ph_flag = os.environ.get("PH_ENSEMBLE")
    if env_ph_flag is not None:
        cfg.PH_ENSEMBLE = _to_bool(env_ph_flag)
    else:
        cfg.PH_ENSEMBLE = _to_bool(cfg.PH_ENSEMBLE)

    scope_env = os.environ.get("PH_SCOPE")
    if scope_env is not None:
        cfg.PH_SCOPE = scope_env.strip()
    elif "PH_SCOPE" in cfg:
        cfg.PH_SCOPE = cfg["PH_SCOPE"]

    radius_env = os.environ.get("PH_RADIUS")
    if radius_env is not None:
        try:
            cfg.PH_RADIUS = float(radius_env)
        except Exception:
            cfg.PH_RADIUS = 10.0
    else:
        try:
            cfg.PH_RADIUS = float(cfg.PH_RADIUS)
        except Exception:
            cfg.PH_RADIUS = 10.0
    if cfg.PH_RADIUS <= 0:
        cfg.PH_RADIUS = 10.0

    log = logging.getLogger("ph_ensemble")
    log.info("[ph_ensemble.mode] enabled=%s scope=%s radius=%s", cfg.PH_ENSEMBLE, getattr(cfg, "PH_SCOPE", "auto"), getattr(cfg, "PH_RADIUS", 10.0))
    # --- PH-ligand mode (CLI > ENV > CFG) ---
    cfg.setdefault("PH_LIGAND_MODE", "off")
    cli_ph_mode = _cli_val(sys.argv, "--ph-ligand-mode")
    env_ph_mode = os.environ.get("PH_LIGAND_MODE")
    cfg_ph_mode = str(cfg.get("PH_LIGAND_MODE", "off"))

    effective_ph_mode = next(
        (
            m
            for m in (cli_ph_mode, env_ph_mode, cfg_ph_mode)
            if m is not None and str(m).strip() != ""
        ),
        "off",
    )
    cfg["PH_LIGAND_MODE"] = effective_ph_mode

    log.info(
        "[ph_ligand.mode] effective=%r (source=%s)",
        effective_ph_mode,
        "CLI" if cli_ph_mode else "ENV" if env_ph_mode else "CFG",
    )

    # --- Per-run configs (RUN_DIR) ---
    cfg.setdefault("CONFIGS_DIR", str(Path(cfg["OVERALL_DIR"]) / "configs"))
    cfg.setdefault("RESET_CONFIGS", True)


    # CLI > ENV > CFG
    cli_cfg_dir = _cli_val(sys.argv, "--configs-dir")
    cli_no_reset = _cli_has(sys.argv, "--no-reset-configs")

    if cli_cfg_dir:  cfg["CONFIGS_DIR"] = cli_cfg_dir
    cfg["RUN_ID"] = run_id
    if cli_no_reset: cfg["RESET_CONFIGS"] = False

    init_config_run_dir(cfg, run_id=cfg.get("RUN_ID"), reset=cfg.get("RESET_CONFIGS"),
                        logger=logging.getLogger("run"))
    print(f"[cfg.run] run_id={cfg['RUN_ID']} run_dir={cfg['CONFIG_RUN_DIR']}")
    print(f"[ph.mode] PH_ENSEMBLE={cfg.get('PH_ENSEMBLE', False)}")

    # --- Single-ligand config (ported) ---------------------------------------
    cfg.setdefault("SINGLE_LIGAND", "")
    cfg.setdefault("SINGLE_LIGAND_SEARCH_ORDER", "fda_library,per_protein,global")
    cfg.setdefault("SINGLE_LIGAND_ALLOW_PREFIX", False)
    cfg.setdefault("SINGLE_LIGAND_MANIFEST_ONLY", True)
    cfg.setdefault("SINGLE_LIGAND_SUGGESTIONS", 5)
    cfg.setdefault("ALLOW_FDA_FALLBACK", False)
    cfg.setdefault("LIBRARY_MANIFEST_FILENAME", "_manifest.json")
    # Default + env override for building manifests during fallback scans
    cfg.setdefault("LIBRARY_MANIFEST_BUILD_ON_SCAN", True)
    env_build_flag = os.environ.get("LIBRARY_MANIFEST_BUILD_ON_SCAN")
    if env_build_flag is not None:
        try:
            cfg["LIBRARY_MANIFEST_BUILD_ON_SCAN"] = _to_bool(env_build_flag)
        except Exception:
            # If parsing fails, keep the config/default
            pass

    logging.getLogger("lib-manifest").info(
        "[lib-manifest.scan] build_on_scan=%s",
        str(bool(cfg.get("LIBRARY_MANIFEST_BUILD_ON_SCAN", True))).lower(),
    )
    cfg.setdefault("FDA_MAPPING_CSV", str(Path(__file__).with_name("fda_mapping_from_pdbqt.csv")))

    # CLI > ENV > CFG precedence
    cli_single = _parse_single_from_cli(sys.argv)
    env_single = os.environ.get("SINGLE_LIGAND", "").strip()
    cfg_single = str(cfg.get("SINGLE_LIGAND", "")).strip()

    effective_single = next((x for x in (cli_single, env_single, cfg_single) if x), "")
    cfg["_EFFECTIVE_SINGLE_LIGAND"] = effective_single
    if effective_single:
        print(f"[config] SINGLE_LIGAND effective='{effective_single}' "
              f"(order=CLI>{'ENV' if env_single else ''}>{'CFG' if cfg_single else ''})")

    # --- Library subfolder selection -----------------------------------
    cfg.setdefault("LIBRARY_SUBDIR_DEFAULT", "fda_library")
    cfg.setdefault("TEST_MODE_ENABLE", "off")
    # Accept dict or JSON-ish string
    if "TEST_LIBRARY_MAP" not in cfg:
        cfg["TEST_LIBRARY_MAP"] = {}
    # --- Specified Proteins Mode ---------------------------------------
    cfg.setdefault("SPECIFIED_PROTEINS", "")
    requested_ids, _sel_src = _parse_specified_proteins(sys.argv, cfg)
    cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] = requested_ids
    print(f"[config] SPECIFIED_PROTEINS effective={requested_ids} (precedence: CLI>ENV>CFG)")
    # --- Fast mode: force exhaustiveness=1 everywhere ---
    cfg["FAST_MODE"] = _parse_fast_flag(sys.argv) or bool(cfg.get("FAST_MODE", False))
    if cfg["FAST_MODE"]:
        print("[config] FAST_MODE effective=True (exhaustiveness=1)")

    # --- Center selection knobs (safe defaults) ---
    cfg.setdefault("CENTER_MODE", "control-first")  # ["control-first","hybrid","library-first"]
    cfg.setdefault("CONTROL_BLACKLIST", "GOL,EDO,PG4,MPD,ACT,SO4,PO4,CL,NA,CA")
    cfg.setdefault("CONTROL_MIN_HEAVY_ATOMS", 10)
    cfg.setdefault("CONTROL_ANCHOR_MIN_VALID_RATE", 0.10)  # if current cluster has control hits + =10% valid, anchor
    cfg.setdefault("ALLOW_SWITCH_FROM_CONTROL", True)
    cfg.setdefault("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False)
    cfg.setdefault("SWITCH_AWAY_FROM_CONTROL_MIN_BOOST", 2.5)  # kcal/mol median boost needed to leave control
    # Optional lock score gate (kcal/mol). Use a large positive number (or remove) to lock on RMSD alone.
    cfg.setdefault("CONTROL_LOCK_SCORE_MAX", -6.0)
    cfg.setdefault("CONTROL_LOCK_MIN_HITS", 1)  # require = this many validated controls
    cfg.setdefault("CONTROL_LOCK_CENTER_MAX_DIST", 4.0)  # A; control centroid must be within this of center

    # clustering + switching thresholds
    cfg.setdefault("CLUSTER_EPS_ANG", 3.5)
    cfg.setdefault("SWITCH_VALID_RATE_MIN", 0.40)
    cfg.setdefault("SWITCH_SCORE_IMPROVE_MIN", 1.5)
    cfg.setdefault("SWITCH_LE_GAIN_MIN", 0.02)
    cfg.setdefault("SWITCH_SCORE_THRESHOLD", 0.70)
    cfg.setdefault("SWITCH_SCORE_HYSTERESIS", 0.50)

    # Hard cap on global switching (early recenter, empty-stage fallback, selector promotions)
    cfg.setdefault("MAX_GLOBAL_CENTER_SWITCHES", 2)

    # Critical defaults
    cfg.setdefault("THREADS_PER_VINA", 1)
    cfg.setdefault("RMSD_FILTER_ANG", 2.0)
    cfg.setdefault("RMSD_MAX_MODELS", 3)

    # --- safe defaults ---
    cfg.setdefault("RETRY_NEAR_MISS", True)
    cfg.setdefault("RETRY_EXHAUST_MULT", 2)
    cfg.setdefault("RETRY_DIST_THRESH", 7.0)
    cfg.setdefault("RETRY_SCORE_THRESH", -7.5)

    cfg.setdefault("ADAPTIVE_SHRINK_ENABLE", True)
    cfg.setdefault("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0)
    cfg.setdefault("ADAPTIVE_SHRINK_DEC", 4.0)
    cfg.setdefault("ADAPTIVE_SHRINK_MIN_BOX", 14.0)

    cfg.setdefault("CHECKPOINT_ENABLE", True)
    cfg.setdefault("PAINS_NAMES", "")
    cfg.setdefault("RECEPTOR_SANITY_CHECK", True)
    cfg.setdefault("AUDIT_JSON", True)

    # --- logging/noise controls ---
    cfg.setdefault("QUIET_CONSOLE", False)   # console shows WARN+ only; file keeps DEBUG
    cfg.setdefault("VINA_VERBOSITY", 2)     # 0=minimal, 1=normal, 2=verbose
    cfg.setdefault("FILTER_VINA_STDOUT", False)  # reserved if we need extra filtering later

    #RMSD PARAMETERS
    cfg.setdefault("SELF_RMSD_MAX_ANG", 2.0)
    cfg.setdefault("SELF_RMSD_REQUIRE_FOR_CONTROLS", True)  # reserved for future stricter gating
    cfg.setdefault("EARLY_EXIT_MAX_MODELS", 3)  # validate at most N poses, stop on first PASS
    cfg.setdefault("MAX_RETRY_SECONDS_PER_LIGAND", 300)  # wall-clock for retries/validation per ligand

    params = get_recenter_params(cfg)
    stages = define_docking_stages(cfg.get("DOCKING_MODE", "discovery").lower())
    print("current docking mode is ", cfg.get("DOCKING_MODE"))

    # Discover all candidate PDB files (unchanged default behavior)
    pdb_files = [
        f for f in os.listdir(cfg["INPUT_DIR"])
        if f.lower().endswith(".pdb") and "_nolig" not in f.lower()
    ]

    # Build an index: PDBID (4-char, upper) -> filename
    id_index: dict[str, str] = {}
    for f in pdb_files:
        base = os.path.splitext(f)[0].replace("_cleaned", "")
        nid = _norm_pdb_id(base)
        if nid:
            # preserve first occurrence to retain directory order
            id_index.setdefault(nid, f)

    req = list(cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS", []) or [])
    if req:
        # Compute present/missing and apply filter in user-specified order
        hits = [nid for nid in req if nid in id_index]
        miss = [nid for nid in req if nid not in id_index]

        print(f"[filter.proteins] mode=on requested={len(req)} present={len(hits)} missing={len(miss)} ? {hits}")
        for m in miss:
            print(f"WARNING: requested PDB '{m}' not found under INPUT_DIR={cfg['INPUT_DIR']} or was excluded (_nolig).")

        if not hits:
            print("ERROR: No requested proteins found. Exiting with status 2 to avoid a no-op run.")
            sys.exit(2)

        # Restrict queue to the selected files, preserving user order
        pdb_files = [id_index[nid] for nid in hits]
        print("Selected proteins (Specified Proteins Mode): " + ", ".join(hits))
    else:
        print(f"[filter.proteins] mode=off requested=0 present={len(pdb_files)} missing=0 ? []")



    # --- Test-mode protein filter: keep only PDBs listed in TEST_LIBRARY_MAP ---
    test_mode = _resolve_test_mode(cfg)
    if test_mode != "off":
        raw_map = cfg.get("TEST_LIBRARY_MAP", {})
        test_map = _coerce_test_map(raw_map)
        test_keys = {k[:4] for k in test_map.keys()}
        if test_keys:
            kept, skipped = [], []
            for f in pdb_files:
                nid = _norm_pdb_id(f)
                if nid and nid.upper()[:4] in test_keys:
                    kept.append(f)
                else:
                    skipped.append(f)

            if skipped:
                print(f"[test-mode] Enabled mode={test_mode}; restricting to {len(kept)} PDBs from TEST_LIBRARY_MAP keys.")
                for s in skipped:
                    print(f"[test-mode] Skipping {s} (not in TEST_LIBRARY_MAP).")

            pdb_files = kept
        else:
            print("[test-mode] TEST_LIBRARY_MAP empty/invalid; no extra filtering applied.")

    print("Working directory:", os.getcwd())
    print("Loaded config keys:", list(cfg.keys()))
    print(f"Proteins queued: {len(pdb_files)}")


    start = time.time()
    plan_only = os.environ.get("A2_PLAN_ONLY") == "1"
    mode, variants = resolve_apo_holo_mode(cfg)
    router_legacy = (mode == "legacy")
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))
    cfg_raw_mode = cfg.get("APO_HOLO_MODE")
    logging.info(
        "[apo-holo] cfg_token_raw=%r resolved_mode=%s variants=%s",
        cfg_raw_mode,
        mode,
        variants,
    )
    logging.info(
        "[apo-holo] router_legacy=%s ph_enabled=%s",
        router_legacy,
        ph_enabled,
    )
    logging.info(
        "[apo-holo] normalized_mode_token=%s",
        _debug_normalize_mode_token(cfg_raw_mode),
    )
    cfg["_ROUTER_LEGACY"] = router_legacy
    cfg["_RESOLVED_APO_HOLO_MODE"] = mode
    logging.info(
        "[apo-holo] resolved mode=%s variants=%s cfg_token=%r",
        mode,
        variants,
        cfg_raw_mode,
    )



    # Where to write per-PDB failure logs
    overall_dir = cfg.get("OVERALL_DIR", ".")
    failed_root = os.path.join(overall_dir, "failed")
    os.makedirs(failed_root, exist_ok=True)
    logging.info("[apo-holo] failed log directory: %s", failed_root)

    failed_entries = []  # (pdb_id, variant_label, log_path, exc_type, exc_msg)

    for variant in variants:
        # Make variant visible to any module still reading env (legacy compatibility)
        os.environ["APO_HOLO_VARIANT"] = "" if variant is None else str(variant).upper()
        label = "legacy" if variant is None else str(variant).lower()
        logging.info(
            "[apo-holo] start_variant mode=%s variant_label=%s env_token=%s",
            mode,
            label,
            os.environ.get("APO_HOLO_VARIANT", ""),
        )

        with tqdm(
                total=len(pdb_files),
                desc=f"Processing Proteins ({label})",
                unit="protein",
                position=0,
                dynamic_ncols=True,
                mininterval=0.2,
                leave=True,
                file=sys.stdout,
        ) as bar:
            cfg_v = cfg  # no per-variant mutation; variant is propagated via APO_HOLO_VARIANT env

            for pdb_file in pdb_files:
                pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()

                try:
                    # Main per-PDB work
                    process_one_protein(cfg_v, pdb_file, stages, params)


                except Exception as exc:
                    # Per-PDB failure handling
                    exc_type = type(exc).__name__
                    exc_msg = str(exc)
                    traceback_str = traceback.format_exc()

                    fail_log_path = Path(failed_root) / f"{pdb_id}.{label}.log"

                    # Write a dedicated failure log for this PDB+variant
                    with fail_log_path.open("w", encoding="utf-8") as fh:
                        fh.write(
                            f"[FAILED PDB]\n"
                            f"  pdb_id        = {pdb_id}\n"
                            f"  variant_label = {label}\n"
                            f"  mode          = {mode}\n"
                            f"  pdb_file      = {pdb_file}\n"
                            f"  exception     = {exc_type}: {exc_msg}\n\n"
                            f"[TRACEBACK]\n"
                            f"{traceback_str}\n"
                        )
                        traceback.print_exc(file=fh)

                    logging.error(
                        "[apo-holo] FAILED pdb_id=%s variant_label=%s; "
                        "see failure log at %s",
                        pdb_id,
                        label,
                        fail_log_path,
                    )

                    failed_entries.append(
                        (pdb_id, label, str(fail_log_path), exc_type, exc_msg)
                    )

                finally:
                    # Always advance the progress bar, even if this PDB failed
                    bar.update(1)

    elapsed_min = (time.time() - start) / 60.0
    print(f"\nAll proteins processed in {elapsed_min:.2f} minutes.")

    if failed_entries:
        print("\nThe following proteins failed. See per-PDB logs under:", failed_root)
        for pdb_id, label, log_path, exc_type, exc_msg in failed_entries:
            print(
                f"  - {pdb_id} ({label}): {exc_type} — {exc_msg}\n"
                f"      log: {log_path}"
            )
    else:
        print("\nNo proteins recorded as failed.")

    if plan_only:
        sys.exit(0)


def _send_run_email(status: int, start_time: str, end_time: str) -> None:
    """
    Best-effort email notification using the system `mail` command.
    Must never raise, so it is safe to call from finally blocks.
    """
    try:
        # Hostname: use uname if available (Linux/Unix), else fall back.
        try:
            host = os.uname().nodename
        except AttributeError:
            host = "unknown-host"

        subject = f"Atlas run exited with status {status} on {host}"

        # Mirror your shell script body as closely as possible
        cmd_line = "python main.py " + " ".join(sys.argv[1:])
        body_lines = [
            f"Atlas run finished on host: {host}",
            f"Start time: {start_time}",
            f"End time:   {end_time}",
            f"Exit status: {status}",
            "",
            "Command:",
            cmd_line,
            "",
        ]
        body = "\n".join(body_lines)

        # Use the same `mail` CLI you already tested in your bash wrapper
        try:
            pipe = os.popen(f'mail -s "{subject}" mpg2352@utexas.edu', "w")
            try:
                pipe.write(body)
            finally:
                pipe.close()
        except Exception:
            # If mail fails, log it but never break the run
            try:
                logging.exception("[notify] failed to send mail notification")
            except Exception:
                # Logging itself should not be able to kill the run
                pass

    except Exception:
        # Absolute last-resort guard: never let notification kill the process
        try:
            logging.exception("[notify] unexpected error while building notification email")
        except Exception:
            pass


def _run_with_email_notification() -> None:
    """
    Wrapper used ONLY when main.py is invoked as a script.

    - Calls _smoke_emit_config_demo() and main() in the same order as before.
    - Preserves all exit codes (SystemExit, unhandled exceptions).
    - Always attempts to send an email in a finally block.
    """
    # Import-time already brought in `time`, `os`, `sys`, `logging`, etc.
    start_time = time.strftime("%Y-%m-%d %H:%M:%S")
    status: int = 0

    try:
        _smoke_emit_config_demo()
        main()
        # If main returns normally, status 0
        status = 0
    except SystemExit as exc:
        # Preserve the original exit code from sys.exit()
        code = exc.code
        status = code if isinstance(code, int) else 1
        raise
    except BaseException:
        # KeyboardInterrupt and other errors → non-zero status
        status = 1
        raise
    finally:
        end_time = time.strftime("%Y-%m-%d %H:%M:%S")
        try:
            _send_run_email(status=status, start_time=start_time, end_time=end_time)
        except Exception:
            # Never let notification interfere with the original exit behavior
            try:
                logging.exception("[notify] email wrapper raised unexpectedly")
            except Exception:
                pass


if __name__ == "__main__":
    _run_with_email_notification()

