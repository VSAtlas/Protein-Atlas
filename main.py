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
from dataclasses import dataclass, field
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
    extract_best_score, emit_vina_config as _emit_vina_config_impl, _to_bool, init_config_run_dir
)
from docking import (
    norm,
    _fingerprint_stage,
    select_ligands_for_next,
    final_pose_validation_and_screenshots,
    RetryManager,
    run_one_stage,
    _coerce_test_map,
    _resolve_test_mode,
    process_one_protein,
)
from path_router import (
    expand_variants,
    make_paths,
    Paths as RouterPaths,
    receptor_file,
    docked_dir,
    config_dir as router_config_dir,
)
from fallback_recenter import RecenterParams
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






# >>> PATHS CLASS START
Paths = RouterPaths
# >>> PATHS CLASS END

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

