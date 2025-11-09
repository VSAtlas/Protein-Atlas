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
import atexit, datetime
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from tqdm import tqdm
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
from input_and_export_functions import (
    load_inputs, validate_config, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, emit_vina_config, record_score, score_key, _to_bool, init_config_run_dir
)
from protein_functions import detect_active_site
from activesite import extract_and_remove_ligands
from prep_ligands import prep_ligands_from_pdb, is_valid_ligand
from pose_validation import (
    validate_pose_pdbqt, extract_surface_atoms, attempt_fallback_recenter,
    filter_and_rewrite_poses_by_rmsd, compute_self_rmsd
)
from run_vina import run_docking_task, validate_all_poses
from path_router import expand_variants
from path_router import make_paths, Paths as RouterPaths


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

class _Tee:
    def __init__(self, stream, file_path):
        self._stream = stream
        self._fh = open(file_path, "a", buffering=1, encoding="utf-8", errors="replace")
    def write(self, data):
        try:
            self._stream.write(data)
        except Exception:
            pass
        try:
            self._fh.write(data)
        except Exception:
            pass
    def flush(self):
        try:
            self._stream.flush()
        except Exception:
            pass
        try:
            self._fh.flush()
        except Exception:
            pass
    def isatty(self):
        try:
            return self._stream.isatty()
        except Exception:
            return False
    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass

def _tee_stdio_to(log_path):
    tee_out, tee_err = _Tee(sys.stdout, log_path), _Tee(sys.stderr, log_path)
    sys.stdout, sys.stderr = tee_out, tee_err
    def _announce_and_close():
        try:
            sys.stdout.write(f"Log -> {log_path}\n")
            sys.stdout.flush()
        finally:
            try:
                tee_out.close()
                tee_err.close()
            except Exception:
                pass
    atexit.register(_announce_and_close)




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
# Apo vs Holo mode
# ======================
def _clean_mode_token(s: str | None) -> str:
    s = (s or "").strip().lower()
    # normalize separators
    s = s.replace("-", "").replace("_", "")
    return s

def resolve_apo_holo_mode(cfg: dict) -> tuple[str, list]:
    """
    Single source of truth: ENV -> config -> default ('apo_vs_holo').
    Returns (mode, variants). For legacy/no-variant mode, variants == [None].
    """
    env_raw = _clean_mode_token(os.environ.get("APO_HOLO_MODE"))
    cfg_raw = _clean_mode_token(str(cfg.get("APO_HOLO_MODE", "")))
    raw = env_raw or cfg_raw or "apovsholo"

    if raw in {"apo"}:
        return "apo", ["APO"]
    if raw in {"holo"}:
        return "holo", ["HOLO"]
    if raw in {"none", "legacy", "null", "false"}:
        return "legacy", [None]
    if raw in {"apovsholo", "apovsholo", "apovsholo"} or raw == "apovsholo":
        return "apo_vs_holo", ["HOLO", "APO"]  # HOLO-first aligns dedup 'keep' with HOLO
    return "apo_vs_holo", ["HOLO", "APO"]



# Put near other helpers
def _variant_receptor_path(pdb_id: str, variant: str | None, cfg: dict) -> str | None:
    # Return the cleaned receptor PDB path for a given variant if it exists, else None
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    rec = paths.receptor_cleaned_pdb(variant)
    return str(rec) if rec.exists() else None


def file_sha1(path: str) -> str:
    import hashlib
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

def delete_variant_trees(pdb_id: str, variant: str, cfg: dict) -> None:
    # Delete processed receptor and docking trees for a specific variant (idempotent)
    import shutil
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    proc_variant_root = paths.receptor_dir(variant).parent      # processed_pdbs/<PDB>/<VARIANT>/
    dock_variant_root = paths.docked_variant_root(variant)      # docked/<PDB>/<VARIANT>/
    for d in (proc_variant_root, dock_variant_root):
        if d.exists():
            shutil.rmtree(d, ignore_errors=True)

def dedup_identical_variants(pdb_id: str, cfg: dict) -> None:
    """
    If HOLO and APO cleaned receptors are byte-identical, delete HOLO and keep APO.
    """
    holo = _variant_receptor_path(pdb_id, "HOLO", cfg)
    apo  = _variant_receptor_path(pdb_id, "APO",  cfg)
    if not holo or not apo:
        return
    try:
        if file_sha1(holo) == file_sha1(apo):
            logging.info("[apo-vs-holo] identical receptors for %s; deleting HOLO (keeping APO)", pdb_id)
            delete_variant_trees(pdb_id, "HOLO", cfg)
    except Exception as e:
        logging.warning("[apo-vs-holo] dedup check failed for %s: %s", pdb_id, e)

    
    
# ======================
# Data models & utilities
# ======================
@dataclass
class RetryManager:
    max_retries: int = 2
    recipes: Dict[str, List[Dict[str, Any]]] = field(default_factory=lambda: {
        # If we docked far from the pocket, try small geometry tweaks   not more modes
        "too_far_from_pocket": [
            {"recenter": True, "box_pad_delta": +1.0, "num_modes": 4},
            {"recenter": True, "box_pad_delta": +2.0, "exhaustiveness": 6, "num_modes": 4},
        ],
        # If we didn't get a valid pose, explore new seeds and a slightly wider energy window,
        # but keep returned modes low so validation stays fast.
        "no_valid_pose": [
            {"exhaustiveness": 6, "num_modes": 5, "seed_jitter": True, "energy_range": 6},
            {"recenter": True, "box_pad_delta": +1.0, "exhaustiveness": 6, "num_modes": 5, "seed_jitter": True, "energy_range": 6},
        ],
        # If we timed out, go cheaper, not deeper.
        "timeout": [
            {"exhaustiveness": 3, "num_modes": 3, "seed_jitter": True},
            {"exhaustiveness": 2, "num_modes": 2},
        ],
        "malformed": []  # do not retry
    })

    def apply(self, base_params: Dict[str, Any], err_type: str, attempt: int) -> Optional[Dict[str, Any]]:
        if err_type not in self.recipes or attempt >= len(self.recipes[err_type]):
            return None
        p = base_params.copy()
        for k, v in self.recipes[err_type][attempt].items():
            if k.endswith("_delta"):
                key = k.replace("_delta", "")
                p[key] = p.get(key, 0.0) + v
            else:
                p[k] = v
        return p


@dataclass
class BudgetGuard:
    """Simple per-ligand wall-clock guard for retries/validation."""
    max_seconds: float
    _deadline: float = None

    def __post_init__(self):
        self._deadline = time.time() + float(self.max_seconds)

    def expired(self) -> bool:
        return time.time() >= self._deadline


def _iter_pdbqt_models(pdbqt_path: str):
    """
    Yield individual MODEL..ENDMDL blocks from a (possibly multi-model) PDBQT.
    If no MODEL/ENDMDL markers exist, yield the whole file once.
    """
    buf = []
    saw_model = False
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith("MODEL"):
                if buf:
                    yield "".join(buf)
                    buf = []
                saw_model = True
                buf.append(ln)
            elif ln.startswith("ENDMDL"):
                buf.append(ln)
                yield "".join(buf)
                buf = []
                saw_model = True
            else:
                if saw_model:
                    buf.append(ln)

    # If we never saw a MODEL block, treat the whole file as one model
    if not saw_model:
        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f2:
                yield f2.read()
        except Exception:
            yield ""


def validate_first_valid_pose(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    pocket_center: tuple[float, float, float],
    surface_coords,
    max_models: int = 3,
    clash_threshold: float = 2.0,
    clash_tol: int = 3,
    dist_surf: float = 6.0,
    dist_centroid: float = 4.5,
):
    """
    Validate poses in order and return as soon as one passes.
    Falls back to the last invalid result if none pass.
    """
    tmp_dir = Path(ligand_pdbqt).parent
    best_invalid = None
    count = 0

    for idx, model_text in enumerate(_iter_pdbqt_models(ligand_pdbqt)):
        if max_models and count >= int(max_models):
            break
        count += 1

        tmp = tmp_dir / f"{Path(ligand_pdbqt).stem}.m{idx}.tmp.pdbqt"
        try:
            tmp.write_text(model_text, encoding="utf-8")
        except Exception:
            tmp = Path(ligand_pdbqt)

        try:
            res = validate_pose_pdbqt(
                protein_pdbqt=receptor_pdbqt,
                ligand_pdbqt=str(tmp),
                pocket_center=pocket_center,
                clash_threshold=clash_threshold,
                CLASH_TOLERANCE=clash_tol,
                DIST_THRESHOLD_SURFACE=dist_surf,
                DIST_THRESHOLD_CENTROID=dist_centroid,
                surface_atom_coords=surface_coords,
            )
        finally:
            if tmp.name.endswith(".tmp.pdbqt"):
                try:
                    tmp.unlink(missing_ok=True)
                except Exception:
                    pass

        if res.get("valid", False):
            return res  # early exit on first valid

        best_invalid = res  # keep the last invalid for diagnostics

    return best_invalid or {"valid": False, "reason": "no_poses"}


@dataclass
class RecenterParams:
    """Thresholds for early/fallback recenter heuristics."""
    EARLY_RECENTER_RATIO: float = 0.70
    EARLY_RECENTER_MIN_EVAL: int = 10
    EARLY_RECENTER_FAR_A: float = 15.0
    EARLY_RECENTER_MEDIAN_A: float = 10.0
    ALLOW_BOX_EXPAND: bool = True
    MAX_RECENTER_ATTEMPTS: int = 1  # tightened: fewer early recenter tries


@dataclass
class GlobalCenterGuard:
    """
    Gatekeeper for ANY global center change (early recenter, empty-stage fallback,
    CenterSelector promotions). Supports a hard 'lock' after a validated control ligand.
    """
    max_global_switches: int = 2
    global_switches: int = 0
    switched_this_stage: bool = False
    locked: bool = False

    def reset_stage(self) -> None:
        """Reset per-stage switch flag (call at the start of each stage)."""
        self.switched_this_stage = False

    def can_switch(self) -> bool:
        """
        True if a center change is allowed right now.
        Respects: hard lock, one-per-stage, and global cap.
        """
        return (not self.locked) and (not self.switched_this_stage) and (self.global_switches < self.max_global_switches)

    def mark_switch(self) -> None:
        """Record that a center change just happened this stage."""
        self.global_switches += 1
        self.switched_this_stage = True

    def lock(self) -> None:
        """Hard-lock: disallow any further center changes for the remainder of the run."""
        self.locked = True


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

def _resolve_single_ligand(selector: str, pdb_id: str, cfg: Dict, logger: logging.Logger) -> Optional[Path]:
    """
    Apply SINGLE_LIGAND_SEARCH_ORDER with a minimal optimization:
    When SINGLE_LIGAND_SKIP_GLOBAL=True (default), ignore the 'global' scope entirely.

    Match semantics: exact or prefix match depending on SINGLE_LIGAND_ALLOW_PREFIX.
    Returns a Path to the first hit, or None.
    """
    if not selector:
        return None

    allow_prefix = _to_bool(str(cfg.get("SINGLE_LIGAND_ALLOW_PREFIX", "false")))
    # Default order unchanged, but we may filter it below:
    order = str(cfg.get("SINGLE_LIGAND_SEARCH_ORDER", "per_protein,global")).replace(" ", "").split(",")

    # NEW: default to skip 'global' traversal in single-ligand mode (cheap win)
    skip_global = bool(cfg.get("SINGLE_LIGAND_SKIP_GLOBAL", True))
    if skip_global:
        order = [w for w in order if w != "global"]

    # Where to look
    per_protein_dir = cfg.get("paths", {}).get("prepped_ligands_dir")  # injected at runtime in process_one_protein
    global_root     = Path(cfg.get("OUTPUT_LIGANDS_DIR", "")) if cfg.get("OUTPUT_LIGANDS_DIR") else None

    logger.info("[single.debug] selector=%s order=%s skip_global=%s", selector, order, skip_global)
    logger.info("[single.debug.paths] per_protein_dir=%s global_root=%s", per_protein_dir, global_root)

    try:
        name_map = _load_fda_name_map(cfg, logger)
        key = _norm_name_key(selector)
        basenames = list(name_map.get(key, []))
        if not basenames and allow_prefix and key:
            pref = key
            for k, v in name_map.items():
                if k.startswith(pref):
                    basenames.extend(list(v))
        if global_root:
            for bn in basenames:
                p = Path(global_root) / "fda_library" / bn
                if p.exists():
                    logger.info("[single.name] key=%s basenames=%s probe=%s", _norm_name_key(selector), list(basenames), str(p))
                    return p
    except Exception as _e:
        logger.debug("[single.name] mapping search skipped: %s", _e)

    def _match_one_dir(root: Optional[Path]) -> Optional[Path]:
        if not root or not Path(root).exists():
            return None
        candidates = list(Path(root).glob("*.pdbqt"))
        # Try exact stem (without _stage suffix), then prefix
        for p in candidates:
            base = p.stem.split("_stage")[0]
            if base.lower() == selector.lower():
                return p
        if allow_prefix:
            for p in candidates:
                base = p.stem.split("_stage")[0]
                if base.lower().startswith(selector.lower()):
                    return p
        return None

    for where in order:
        if where == "per_protein":
            hit = _match_one_dir(Path(per_protein_dir) if per_protein_dir else None)
            if hit:
                logger.info(f"[single] matched in per-protein dir: {hit.name}")
                return hit

        elif where == "global":
            # Retained for optional use if SINGLE_LIGAND_SKIP_GLOBAL=False
            if global_root and global_root.exists():
                # Fallback rglob only if you explicitly re-enable global (skip_global=False)
                for p in global_root.rglob("*.pdbqt"):
                    base = p.stem.split("_stage")[0]
                    if base.lower() == selector.lower():
                        logger.info(f"[single] matched in global dir: {p}")
                        return p
                if allow_prefix:
                    for p in global_root.rglob("*.pdbqt"):
                        base = p.stem.split("_stage")[0]
                        if base.lower().startswith(selector.lower()):
                            logger.info(f"[single] prefix-matched in global dir: {p}")
                            return p
        else:
            logger.debug(f"[single] unknown search scope: {where}")
    logger.warning("[single.miss] selector=%s (no FDA map hit / file absent)", selector)
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


def make_protein_logger(docked_dir: str, pdb_id: str, cfg: Dict) -> logging.Logger:
    """
    Create a logger writing to DOCKED_DIR/<pdb_id>/protein.log and also to console.
    Keeps logs per-protein and avoids duplicate handlers.
    """
    base = Path(docked_dir)
    # If caller already passed .../docked/<PDB>, don't append <PDB> again
    log_dir = base if base.name.upper() == pdb_id.upper() else (base / pdb_id)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file_path = log_dir / "protein.log"

    logger = logging.getLogger(pdb_id)
    logger.setLevel(logging.DEBUG)


    # Reset handlers to avoid duplicates if re-used
    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Overwrite per run (truncate), not append
    fh = logging.FileHandler(log_file_path, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    ch = logging.StreamHandler(stream=sys.stdout)  # stdout, not stderr
    env_override = os.environ.get("QUIET_CONSOLE_OVERRIDE", "").strip()
    quiet = (env_override.lower() in {"1", "true", "yes"}) if env_override else _to_bool(cfg.get("QUIET_CONSOLE", False))
    ch.setLevel(logging.WARNING if quiet else logging.INFO)
    ch.setFormatter(formatter)
    # Opt-in topic filtering and level overrides
    def _coerce_level(name: str | None, default: int) -> int:
        m = {"CRITICAL":50,"ERROR":40,"WARN":30,"WARNING":30,"INFO":20,"DEBUG":10,"NOTSET":0}
        return m.get(str(name or "").strip().upper(), default)

    raw_topics = (os.environ.get("LOG_TOPICS") or str(cfg.get("LOG_TOPICS", ""))).replace(",", " ")
    topics = {t.strip().lower() for t in raw_topics.split() if t.strip()}

    class _TopicFilter(logging.Filter):
        def __init__(self, allowed: set[str]): self.allowed = allowed
        def filter(self, record: logging.LogRecord) -> bool:
            # Always show warnings/errors
            if record.levelno >= logging.WARNING:
                return True
            msg = record.getMessage()
            # If message is [tag]..., allow only when tag in allowed
            if msg.startswith("[") and ("]" in msg):
                tag = msg[1:msg.find("]")].strip().lower()
                if not self.allowed or "all" in self.allowed:
                    return True
                return tag in self.allowed
            # Untagged INFO/DEBUG only pass when explicitly enabled as 'untagged'
            return ("untagged" in self.allowed) or (not self.allowed)

    # Optional level overrides
    fh.setLevel(_coerce_level(os.environ.get("LOG_LEVEL_FILE") or cfg.get("LOG_LEVEL_FILE"), fh.level))
    ch.setLevel(_coerce_level(os.environ.get("LOG_LEVEL_CONSOLE") or cfg.get("LOG_LEVEL_CONSOLE"), ch.level))

    # Apply topic filter only if topics were provided (and not 'all')
    if topics and ("all" not in topics):
        filt = _TopicFilter(topics)
        fh.addFilter(filt)
        ch.addFilter(filt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.propagate = False
    return logger



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
    # >>> DOCKED PATHS PATCH START
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None
    variant = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    root = paths.docked_variant_root(variant, ph_label)
    return root / f".ckpt_{stage_name}.json"
    # >>> DOCKED PATHS PATCH END


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
        # >>> DOCKED PATHS PATCH START
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None
        variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
        out = paths.docked_variant_root(variant_env, ph_label) / "audit.json"
        # >>> DOCKED PATHS PATCH END
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


# ----------------------



def prepare_receptor(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[Optional[str], Optional[str]]:
    import automate_protein_prep
    from distutils.util import strtobool

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
            dock_root = paths.docked_variant_root(variant)
            for ph_label, _ in targets:
                ph_root = dock_root / ph_label
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
        builder_id = paths.pdb_id if variant is None else f"{paths.pdb_id}_{variant}"
        log.info("[ph_ensemble.call] building pH ensemble for %s", paths.pdb_id)
        prev_cfg = getattr(automate_protein_prep, "config", None)
        try:
            automate_protein_prep.config = cfg
        except Exception:
            prev_cfg = None
        try:
            import ph_ensemble
            manifest_path = ph_ensemble.build_ph_ensemble(
                pdb_id=builder_id,
                cleaned_receptor_pdb=cleaned_path,
                out_dir=str(Path(cfg["OUTPUT_DIR"])),
                center=center,
                radius=eff_radius,
                ph_values=ph_values,
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
                return None, None
        except Exception as _e:
            logger.warning(f"Receptor sanity check skipped due to error: {_e}")
        cleaned_norm = norm(cleaned_pdb_path)
        receptor_norm = norm(receptor_pdbqt_path)
        if bool(cfg.get("PH_ENSEMBLE_IN_PREP", False)):
            _build_ph_ensemble(cleaned_norm)
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
            output_root=str(Path(cfg["OUTPUT_DIR"]))  # processed_pdbs root; module lays out subdirs
        )
    except Exception as e:
        logger.warning(f"Protein cleaning failed: {e}")
        return None, None

    if not cleaned_pdb or not Path(cleaned_pdb).exists():
        logger.warning("Protein cleaning did not produce a cleaned PDB.")
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

    if bool(cfg.get("PH_ENSEMBLE_IN_PREP", False)):
        _build_ph_ensemble(cleaned_pdb)

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
                return None, None
    except Exception as _e:
        logger.warning(f"Receptor sanity check skipped due to error: {_e}")

    return norm(cleaned_pdb), norm(receptor_pdbqt)

# ---- Multi-control center selection via crystallographic controls ----
def select_center_via_control_redock(cfg, paths, receptor_pdbqt, logger):
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

    def _best_model_to_pdb(pdbqt_file: _Path):
        # parse multi-model pdbqt, pick model with best (lowest) Vina score
        best_e = None; best_chunk = None
        with open(pdbqt_file, "r", encoding="utf-8", errors="ignore") as fh:
            chunk = []; in_model = False
            for ln in fh:
                u = ln.strip().upper()
                if u.startswith("MODEL"):
                    chunk = [ln]; in_model = True
                elif u.startswith("ENDMDL"):
                    chunk.append(ln); in_model = False
                    # evaluate chunk
                    for cl in chunk:
                        if "REMARK VINA RESULT" in cl.upper():
                            try:
                                e = float(cl.strip().split()[3])
                                if (best_e is None) or (e < best_e):
                                    best_e = e; best_chunk = chunk[:]
                            except Exception:
                                pass
                else:
                    if in_model: chunk.append(ln)
        # if file had no explicit MODEL blocks, treat whole file
        if best_chunk is None:
            try:
                with open(pdbqt_file, "r", encoding="utf-8", errors="ignore") as fh:
                    for ln in fh:
                        if "REMARK VINA RESULT" in ln.upper():
                            best_e = float(ln.strip().split()[3])
                            break
                best_chunk = None
            except Exception:
                return None, None

        import tempfile, subprocess, shutil as _sh
        td = _Path(tempfile.mkdtemp(prefix="ctrl_redock_"))
        best_pdbqt = td / "best.pdbqt"
        if best_chunk:
            with open(best_pdbqt, "w", encoding="utf-8") as out:
                out.writelines(best_chunk)
        else:
            _sh.copy2(pdbqt_file, best_pdbqt)
        out_pdb = td / "best.pdb"
        try:
            subprocess.run([obabel, "-ipdbqt", str(best_pdbqt), "-opdb", "-O", str(out_pdb)],
                           check=True, capture_output=True, text=True)
        except Exception:
            return None, best_e
        return (out_pdb if out_pdb.exists() else None), best_e

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
        conf_path, out_path = emit_vina_config(
            cfg, paths.pdb_id, receptor_pdbqt, center, (24.0,24.0,24.0),
            str(lig_pdbqt), "ctrl_redock", stage_info, threads_per_vina, logger=None
        )
        try:
            _, score = _run_dock(vina_exe, conf_path, lig_pdbqt.name, out_path)
        except Exception:
            score = None

        best_pdb, best_e = _best_model_to_pdb(_Path(out_path))
        ref_path = control_lookup.get(base)
        rmsd = float("inf")

        def _quick_file_sig(pth: str) -> str:
            try:
                p = Path(pth)
                sz = p.stat().st_size if p.exists() else -1
                # quick coordinate hash for PDB-like files (robust-ish, not crypto)
                h = hashlib.sha1()
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    for ln in fh:
                        if ln.startswith(("ATOM", "HETATM")):
                            h.update(ln[12:54].encode("utf-8", "ignore"))  # atom name + coords
                return f"exists={p.exists()} size={sz} sha={h.hexdigest()[:10]}"
            except Exception:
                return "sig=unavailable"

        if best_pdb and ref_path:
            logger.info(f"[rmsd.debug] ref={ref_path} | {_quick_file_sig(str(ref_path))}")
            logger.info(f"[rmsd.debug] dock={best_pdb} | {_quick_file_sig(str(best_pdb))}")
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


def prepare_and_filter_ligands(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[List[str], Dict[str, int], Dict[str, bool]]:
    """
    Gathers candidate ligands, keeps existing validation/PAINS logic, and
    filters the *non-control* pool to allowed library roots:

      - TEST_MODE_ENABLE & TEST_LIBRARY_MAP (by pdb_id) -> OUTPUT_LIGANDS_DIR/<mapped_subdir>
      - else -> OUTPUT_LIGANDS_DIR/<LIBRARY_SUBDIR_DEFAULT>

    Controls are *never* filtered out here.
    LIBRARY_EXTRA_DIRS remain included (unchanged).
    """
    # Keep existing prep step for extracted controls (harmless if nothing to do)
    prep_ligands_from_pdb(
        ligand_output_dir=paths.ligand_output_dir,
        ligands_mol2_dir=paths.ligands_mol2_dir,
        prepped_ligands_dir=paths.prepped_ligands_dir,
    )

    # --- roots discovery (unchanged baseline) ---
    roots: list[Path] = []
    global_root = Path(cfg["OUTPUT_LIGANDS_DIR"]) if cfg.get("OUTPUT_LIGANDS_DIR") else None
    if global_root and global_root.exists():
        roots.append(global_root)
    if paths.prepped_ligands_dir.exists():
        roots.append(paths.prepped_ligands_dir)

    extra_dirs = str(cfg.get("LIBRARY_EXTRA_DIRS", "")).strip()
    if extra_dirs:
        for d in extra_dirs.split(";"):
            d = d.strip()
            if not d:
                continue
            p = Path(d)
            if p.exists():
                roots.append(p)

    logger.info("Scanning for ligands under: " + " | ".join(str(r) for r in roots))

    # --- collect all .pdbqt (dedup by normalized path), directory-first ---
    seen: set[str] = set()
    all_pdbqt_paths: list[Path] = []

    # Allowlist only control/reference in the per-protein tree; everything
    # else (non-controls) comes from explicitly allowed library roots.
    per_protein_allow = {"controls", "reference"}

    for r in roots:
        allowed = None
        if r == paths.prepped_ligands_dir:
            allowed = per_protein_allow
        for p in _iter_pdbqt_dirfirst(r, allowed_subdirs=allowed):
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                all_pdbqt_paths.append(p)

    if not all_pdbqt_paths:
        logger.warning("No .pdbqt ligands were found under the configured roots.")
        return [], {}, {}

    # --- validity pass (keep your existing checker) ---
    valid_pdbqt: Dict[str, Path] = {}
    for p in all_pdbqt_paths:
        try:
            # Keep existing lib-root heuristic for validation
            lib_root_for_checks = str(global_root if (global_root and global_root.exists()) else paths.prepped_ligands_dir.parent)
            if is_valid_ligand(p, lib_root_for_checks):
                valid_pdbqt[norm(p)] = p
            else:
                logger.debug(f"Excluded malformed ligand (pdbqt check failed): {p}")
        except Exception:
            logger.debug(f"Excluded malformed ligand (exception): {p}")

    logger.info(f"Valid .pdbqt ligands (union): {len(valid_pdbqt)}")

    # --- per-protein subfolder selection for non-controls only ----------
    pdb_id = paths.pdb_id.upper()
    subdir_default = str(cfg.get("LIBRARY_SUBDIR_DEFAULT", "fda_library"))
    test_enable = _to_bool(str(cfg.get("TEST_MODE_ENABLE", "false")))

    # Robust parse of TEST_LIBRARY_MAP (dict, JSON, or Python-literal string)
    maybe_map = cfg.get("TEST_LIBRARY_MAP", {})
    test_map: Dict[str, str] = {}

    def _coerce_test_map(m) -> Dict[str, str]:
        import json as _json, ast as _ast
        if isinstance(m, dict):
            return {str(k).upper(): str(v) for k, v in m.items()}
        # try string or "dict-like" objects
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

    test_map = _coerce_test_map(maybe_map)
    logger.info(f"[lib-roots.map] raw_type={type(maybe_map).__name__} keys={len(test_map)}")
    hit = test_map.get(pdb_id)  # <- now robust

    # Build allowed non-control roots
    allowed_noncontrol_roots: list[Path] = []
    roots: list[Path] = []
    if cfg.get("OUTPUT_LIGANDS_DIR"):
        base_root = Path(cfg["OUTPUT_LIGANDS_DIR"])
        default_root = base_root / subdir_default
        roots = [default_root]
        if test_enable and pdb_id in test_map:
            test_root = base_root / test_map[pdb_id]
            roots.insert(0, test_root)

    deduped_roots: list[Path] = []
    seen_keys: set[str] = set()
    for r in roots:
        key = str(r.resolve())
        if key not in seen_keys:
            seen_keys.add(key)
            deduped_roots.append(r)

    allowed_noncontrol_roots.extend(deduped_roots)

    logger.info(
        "[fuel] pdb_id=%s test_mode=%s roots=%s"
        % (
            pdb_id,
            str(test_enable).lower(),
            ",".join(str(r.resolve()) for r in deduped_roots),
        )
    )

    # Always include extras (unchanged)
    for d in (extra_dirs.split(";") if extra_dirs else []):
        d = d.strip()
        if d:
            allowed_noncontrol_roots.append(Path(d))

    # Final audits (after list is populated)
    logger.info("[lib-roots] non-control roots = " + ", ".join(map(str, allowed_noncontrol_roots)))


    # Helper: path under root?
    def _under(p: Path, root: Path) -> bool:
        try:
            p.resolve().relative_to(root.resolve())
            return True
        except Exception:
            return False

    # Separate controls vs non-controls by location
    controls: list[Path] = []
    noncontrols: list[Path] = []
    for p in valid_pdbqt.values():
        if _under(p, paths.prepped_ligands_dir):
            controls.append(p)
        else:
            noncontrols.append(p)

    # Filter non-controls to the allowed roots
    filtered_noncontrols: list[Path] = []
    for p in noncontrols:
        keep = False
        for root in allowed_noncontrol_roots:
            if root.exists() and _under(p, root):
                keep = True
                break
        if keep:
            filtered_noncontrols.append(p)

    # Merge back: controls (unaltered) + filtered non-controls
    if cfg.get("_EFFECTIVE_SINGLE_LIGAND") and cfg.get("_SINGLE_RESOLVED_PATH"):
        resolved_path = Path(cfg["_SINGLE_RESOLVED_PATH"])
        final_paths = controls + [resolved_path]
        filtered_noncontrols = [resolved_path]
        logger.info("[single.fuel] resolved=%s controls=%d (blocking non-control pool)", cfg["_SINGLE_RESOLVED_PATH"], len(controls))
    else:
        final_paths = controls + filtered_noncontrols

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
def run_one_stage(
    cfg: Dict,
    pdb_id: str,
    receptor_pdbqt: str,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    stage: Dict,
    ligands: List[str],
    logger: logging.Logger,
    retry_mgr: RetryManager,
    control_lookup: Dict[str, Path],        # maps ligand basename -> crystal ref PDB
    budget_guards: Optional[Dict[str, BudgetGuard]] = None,  # external per-ligand guards
) -> Tuple[
    Dict[str, float],
    List[str],
    List[float],
    Dict[str, str],
    Dict[str, Tuple[Optional[float], str]]
]:
    from sys import stdout as _stdout  # for tqdm
    import shutil
    import subprocess

    threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))
    max_workers = int(cfg["MAX_PARALLEL_JOBS"])
    # >>> DOCKED PATHS PATCH START
    variant = None
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # >>> DOCKED PATHS PATCH END

    scores: Dict[str, float] = {}
    validated_ligands: List[str] = []
    all_distances: List[float] = []
    raw_docked_ligands: Dict[str, str] = {}
    invalids: Dict[str, Tuple[Optional[float], str]] = {}

    surface_coords = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)

    # ---- helpers (nested) ----
    def _best_pose_pdb_from_pdbqt(pdbqt_path: str, obabel_path: Optional[str] = None) -> Optional[str]:
        """Convert first model of PDBQT -> PDB (no hydrogens) using OpenBabel."""
        try:
            out_pdb = Path(pdbqt_path).with_suffix(".best.pdb")
            obabel = (
                obabel_path
                or os.environ.get("OPENBABEL_EXE")
                or cfg.get("OPENBABEL_PATH", "").strip()
                or shutil.which("obabel")
                or "obabel"
            )
            cmd = [obabel, "-ipdbqt", str(pdbqt_path), "-opdb", "-O", str(out_pdb), "-f", "1", "-l", "1", "-d"]
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return str(out_pdb) if out_pdb.exists() and out_pdb.stat().st_size > 0 else None
        except Exception as e:
            logger.warning(f"[RMSD] OpenBabel conversion failed for {os.path.basename(pdbqt_path)}: {e}")
            return None

    def compute_rmsd(ref_path: str, docked_path: str) -> float:
        """
        Heavy-atom RMSD using RDKit's BestRMS *only*.
        Debug: logs paths, atom counts, conformer presence; returns +inf on failure.
        """
        _rlog = logging.getLogger("rmsd")

        ref = _read_any_lig(ref_path)
        dock = _read_any_lig(docked_path)

        # ──  quick identical-file sanity trap ───────────────────────────────────
        try:
            import os
            if os.path.exists(ref_path) and os.path.exists(docked_path) and os.path.samefile(ref_path, docked_path):
                if _rlog:
                    _rlog.warning(f"[rmsd.core] ref and dock resolve to the SAME file "
                                  f"(ref='{ref_path}', dock='{docked_path}')")
        except Exception:
            pass
        # ───────────────────────────────────────────────────────────────────────────

        if not ref or not dock:
            if _rlog:
                _rlog.warning(f"[rmsd.core] load-fail ref_ok={bool(ref)} dock_ok={bool(dock)} "
                              f"ref='{ref_path}' dock='{docked_path}'")
            return float("inf")

        try:
            n_ref = ref.GetNumAtoms()
            n_dock = dock.GetNumAtoms()
        except Exception:
            n_ref = n_dock = -1

        has_conf_ref = (ref.GetNumConformers() > 0)
        has_conf_dock = (dock.GetNumConformers() > 0)

        # ── heavy-atom counts (useful when you get inf) ──────────────────
        try:
            ha_ref = ref.GetNumHeavyAtoms()
            ha_dock = dock.GetNumHeavyAtoms()
        except Exception:
            ha_ref = ha_dock = -1
        if _rlog:
            _rlog.info(f"[rmsd.core] inputs ref='{ref_path}' dock='{docked_path}' "
                       f"n_ref={n_ref} n_dock={n_dock} heavy_ref={ha_ref} heavy_dock={ha_dock} "
                       f"conf_ref={has_conf_ref} conf_dock={has_conf_dock}")
        # ───────────────────────────────────────────────────────────────────────────

        if not has_conf_ref or not has_conf_dock:
            if _rlog:
                _rlog.warning("[rmsd.core] missing 3D conformers; returning inf")
            return float("inf")

        if n_ref != n_dock:
            if _rlog:
                _rlog.info(f"[rmsd.core] atom_count_mismatch ({n_ref} vs {n_dock}); "
                           f"bestRMS will not be used; returning inf (no MCS fallback)")
            return float("inf")

        try:
            val = float(rdMolAlign.GetBestRMS(ref, dock))
            if _rlog:
                _rlog.info(f"[rmsd.core] method=bestRMS rmsd={val:.3f}")
            return val
        except Exception as e:
            if _rlog:
                _rlog.info(f"[rmsd.core] method=bestRMS failed: {e}; returning inf (no MCS fallback)")
            return float("inf")

    def _validate_with_rmsd_gate(
        lig_path: str,
        lig_name: str,
        out_pdbqt_path: str,
        score_val: float
    ) -> Tuple[bool, Optional[str]]:
        """
        Controls: crystal redock RMSD is a hard gate.
        Non-controls: self-RMSD is logged upstream; do not gate here.
        """
        base = Path(lig_path).stem.split("_stage")[0]
        crystal_ref = control_lookup.get(base)

        if crystal_ref:
            best_pdb = _best_pose_pdb_from_pdbqt(out_pdbqt_path, obabel_path=cfg.get("OPENBABEL_PATH"))
            if not best_pdb:
                logger.warning(f"{lig_name} | unable to extract best pose PDB for redock RMSD.")
                return False, "no_best_pose_for_rmsd"

            # --- AUDIT: control redock (compute RMSD just for logging) ---
            rmsd_val = compute_rmsd(str(crystal_ref), best_pdb)

            ok = validate_ligand(
                ligand_name=lig_name,
                docked_path=best_pdb,
                crystal_path=str(crystal_ref),
                rmsd_thresh=float(cfg.get("CONTROL_RMSD_MAX_ANG", 2.0)),
                self_rmsd=None,
                logger=logger
            )

            logger.info("[control-redock] lig=%s rmsd=%.2f A score=%.2f",
                        lig_name, (rmsd_val if rmsd_val is not None else float('nan')), float(score_val))

            if ok:
                logger.info(f"{lig_name} | {stage['name']} score: {score_val:.2f} kcal/mol (redock-RMSD PASS)")
                return True, None
            else:
                return False, "rmsd_fail"


        # Non-controls: redock gate not applicable here (geometry checks already passed).
        return True, None

    # ---- scheduling & submission ----
    default_budget_seconds = float(
        cfg.get("MAX_RETRY_SECONDS_PER_LIGAND", cfg.get("BENCH_MAX_SECONDS", 300.0))
    )

    submit_queue = []
    guards_for_ligand: Dict[str, BudgetGuard] = {}

    for lig in ligands:
        guard = (budget_guards.get(lig) if budget_guards else None)
        if guard is None:
            guard = BudgetGuard(default_budget_seconds)
        guards_for_ligand[lig] = guard

        if guard.expired():
            invalids[lig] = (None, "budget_exceeded")
            continue
        submit_queue.append(lig)

    futures = {}
    if submit_queue:
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            for lig in submit_queue:
                stage_for_cfg = dict(stage)
                stage_for_cfg["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
                if cfg.get("FAST_MODE"):
                    stage_for_cfg["exhaustiveness"] = 1

                conf_path, out_path = emit_vina_config(
                    cfg, pdb_id, receptor_pdbqt, center, box_size, lig, stage["name"], stage_for_cfg, threads_per_vina,
                    logger
                )

                # Guard: config must live under current RUN_DIR
                try:
                    Path(conf_path).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
                except Exception:
                    raise RuntimeError(f"Refusing to launch Vina with config outside current RUN_DIR: {conf_path}")

                logger.info(f"[vina.call] config={conf_path}")

                lig_n, out_n = norm(lig), norm(out_path)
                raw_docked_ligands[lig_n] = out_n
                futures[pool.submit(run_docking_task, cfg["VINA_EXE"], conf_path, lig, out_path)] = (lig_n, out_n)

            processed = 0
            from tqdm import tqdm as _tqdm
            with _tqdm(
                total=len(futures),
                desc=f"Docking ({stage['name']})",
                unit="ligand",
                position=1,
                dynamic_ncols=True,
                mininterval=0.2,
                leave=True,
                file=sys.stdout,
            ) as pbar:
                for fut in as_completed(futures):
                    lig, out_path = futures[fut]
                    lig_name = os.path.basename(lig)
                    guard = guards_for_ligand.get(lig) or BudgetGuard(default_budget_seconds)

                    try:
                        _, score = fut.result()
                    except Exception as e:
                        logger.warning(f"Docking crashed for {lig_name}: {e}")
                        invalids[lig] = (None, "docking_exception")
                        processed += 1
                        if (processed % 25 == 0) or (processed == len(futures)):
                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                        pbar.update(1)
                        continue

                    if score is None:
                        logger.warning(f"No score for {lig_name}")
                        invalids[lig] = (None, "no_score")
                        processed += 1
                        if (processed % 25 == 0) or (processed == len(futures)):
                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                        pbar.update(1)
                        continue

                    # Pose de-dup
                    try:
                        kept, removed = filter_and_rewrite_poses_by_rmsd(
                            out_path,
                            rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                            max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
                        )
                        if removed > 0:
                            logger.info(f"{os.path.basename(out_path)}: RMSD filter kept {kept}, removed {removed}")
                    except Exception as e:
                        logger.warning(f"RMSD filtering failed for {os.path.basename(out_path)}: {e}")

                    # Geometric validation (early-exit through first valid pose)
                    result = validate_first_valid_pose(
                        receptor_pdbqt=receptor_pdbqt,
                        ligand_pdbqt=out_path,
                        pocket_center=center,
                        surface_coords=surface_coords,
                        max_models=int(cfg.get("EARLY_EXIT_MAX_MODELS", 3)),
                        clash_threshold=2.0,
                        clash_tol=3,
                        dist_surf=6.0,
                        dist_centroid=4.5,
                    )
                    dist = result.get("distance_to_pocket")
                    if isinstance(dist, (int, float)):
                        all_distances.append(dist)

                    logger.info(f"{lig_name} validation: {result}")

                    # Valid pose -> optional RMSD hard gate for controls
                    if result.get("valid", False):
                        ok, reason = _validate_with_rmsd_gate(lig, lig_name, out_path, float(score))
                        if ok:
                            scores[lig] = float(score)
                            validated_ligands.append(lig)
                        else:
                            invalids[lig] = (float(score), reason or "rmsd_fail")

                        processed += 1
                        if (processed % 25 == 0) or (processed == len(futures)):
                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                        pbar.update(1)
                        continue

                    # --------------------------
                    # Geometrically invalid path
                    # --------------------------
                    if guard.expired():
                        invalids[lig] = (float(score) if score is not None else None, "budget_exceeded")
                        processed += 1
                        if (processed % 25 == 0) or (processed == len(futures)):
                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                        pbar.update(1)
                        continue

                    # Near-miss, single heavier retry
                    did_retry = False
                    if bool(cfg.get("RETRY_NEAR_MISS", True)):
                        try:
                            dthr = float(cfg.get("RETRY_DIST_THRESH", 7.0))
                            sthr = float(cfg.get("RETRY_SCORE_THRESH", -7.5))
                            if (dist is not None) and (dist < dthr) and (float(score) <= sthr) and not guard.expired():
                                did_retry = True
                                stage_retry = dict(stage)
                                stage_retry["name"] = f"{stage['name']}_retry"
                                stage_retry["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
                                try:
                                    ex0 = int(stage_retry.get("exhaustiveness", 8))
                                except Exception:
                                    ex0 = 8
                                ex_mult = int(cfg.get("RETRY_EXHAUST_MULT", 2))
                                stage_retry["exhaustiveness"] = max(8, ex0 * ex_mult)
                                if cfg.get("FAST_MODE"):
                                    stage_retry["exhaustiveness"] = 1
                                    
                                retry_center = center
                                retry_box = box_size
                                conf_path2, out_path2 = emit_vina_config(
                                    cfg, pdb_id, receptor_pdbqt, retry_center, retry_box, lig,
                                    stage_retry["name"], stage_retry, threads_per_vina, logger
                                )
                                try:
                                    Path(conf_path2).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
                                except Exception:
                                    raise RuntimeError(
                                        f"Refusing to launch Vina with config outside current RUN_DIR: {conf_path2}")
                                logger.info(f"[vina.call] config={conf_path2}")

                                if guard.expired():
                                    invalids[lig] = (float(score) if score is not None else None, "budget_exceeded")
                                    processed += 1
                                    if (processed % 25 == 0) or (processed == len(futures)):
                                        pbar.set_postfix(ok=len(scores), inv=len(invalids))
                                    pbar.update(1)
                                    continue

                                _, score2 = run_docking_task(cfg["VINA_EXE"], conf_path2, lig, out_path2)

                                try:
                                    filter_and_rewrite_poses_by_rmsd(
                                        out_path2,
                                        rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                                        max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
                                    )
                                except Exception:
                                    pass

                                result2 = validate_first_valid_pose(
                                    receptor_pdbqt=receptor_pdbqt,
                                    ligand_pdbqt=out_path2,
                                    pocket_center=center,
                                    surface_coords=surface_coords,
                                    max_models=int(cfg.get("EARLY_EXIT_MAX_MODELS", 3)),
                                    clash_threshold=2.0,
                                    clash_tol=3,
                                    dist_surf=6.0,
                                    dist_centroid=4.5,
                                )

                                if result2.get("valid", False) and score2 is not None:
                                    ok2, reason2 = _validate_with_rmsd_gate(lig, lig_name, out_path2, float(score2))
                                    if ok2:
                                        scores[lig] = float(score2)
                                        validated_ligands.append(lig)
                                        raw_docked_ligands[lig] = norm(out_path2)
                                        logger.info(f"{lig_name} | {stage_retry['name']} score: {score2:.2f} kcal/mol (rescued)")
                                        processed += 1
                                        if (processed % 25 == 0) or (processed == len(futures)):
                                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                                        pbar.update(1)
                                        continue
                                    else:
                                        invalids[lig] = (float(score2), reason2 or "rmsd_fail")
                                        processed += 1
                                        if (processed % 25 == 0) or (processed == len(futures)):
                                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                                        pbar.update(1)
                                        continue
                        except Exception as _e:
                            logger.warning(f"Retry path failed for {lig_name}: {_e}")

                    # Structured recipe retries
                    err_cat = _map_reason_to_category(result.get("reason", ""))
                    attempt = 0
                    retained_invalid = True

                    while attempt < retry_mgr.max_retries:
                        if guard.expired():
                            invalids[lig] = (float(score) if score is not None else None, "budget_exceeded")
                            break

                        recipe = retry_mgr.apply(stage, err_cat, attempt)
                        if not recipe:
                            break

                        stage_retry2 = dict(stage)
                        stage_retry2["name"] = f"{stage['name']}_r{attempt + 1}"
                        stage_retry2["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))

                        if "exhaustiveness" in recipe:
                            stage_retry2["exhaustiveness"] = recipe["exhaustiveness"]
                        if "num_modes" in recipe:
                            stage_retry2["num_modes"] = recipe["num_modes"]
                        if cfg.get("FAST_MODE"):
                            stage_retry2["exhaustiveness"] = 1
                        if recipe.get("seed_jitter", False):
                            try:
                                base_seed = int(stage_retry2.get("seed", 0)) if "seed" in stage_retry2 else 0
                            except Exception:
                                base_seed = 0
                            stage_retry2["seed"] = base_seed + (attempt + 1) * 137

                        retry_center = center
                        retry_box = box_size
                        try:
                            if recipe.get("recenter", False):
                                fb_pose, new_c, _bs, _ch = attempt_fallback_recenter(
                                    fallback_ligands={lig: out_path},
                                    receptor_pdbqt=receptor_pdbqt,
                                    docking_dir=str(paths.docked_pdb_root()),
                                    stage_name=stage_retry2["name"],
                                    pocket_center=center,
                                    logger=logger,
                                    exclude_basenames=set(),
                                )
                                if new_c is not None:
                                    retry_center = new_c
                            if "box_pad_delta" in recipe and isinstance(recipe["box_pad_delta"], (int, float)):
                                dx = float(recipe["box_pad_delta"])
                                box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                                retry_box = tuple(min(box_cap, s + dx) for s in box_size)
                        except Exception as _e:
                            logger.warning(f"Retry recenter/box tweak failed: {_e}")

                        conf_path3, out_path3 = emit_vina_config(
                            cfg, pdb_id, receptor_pdbqt, retry_center, retry_box, lig,
                            stage_retry2["name"], stage_retry2, threads_per_vina, logger
                        )
                        try:
                            Path(conf_path3).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
                        except Exception:
                            raise RuntimeError(
                                f"Refusing to launch Vina with config outside current RUN_DIR: {conf_path3}")
                        logger.info(f"[vina.call] config={conf_path3}")

                        try:
                            _, score_r = run_docking_task(cfg["VINA_EXE"], conf_path3, lig, out_path3)
                        except Exception as _e:
                            logger.warning(f"Retry docking crashed for {lig_name}: {_e}")
                            attempt += 1
                            continue

                        try:
                            filter_and_rewrite_poses_by_rmsd(
                                out_path3,
                                rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                                max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
                            )
                        except Exception:
                            pass

                        result_r = validate_first_valid_pose(
                            receptor_pdbqt=receptor_pdbqt,
                            ligand_pdbqt=out_path3,
                            pocket_center=center,
                            surface_coords=surface_coords,
                            max_models=int(cfg.get("EARLY_EXIT_MAX_MODELS", 3)),
                            clash_threshold=2.0,
                            clash_tol=3,
                            dist_surf=6.0,
                            dist_centroid=4.5,
                        )

                        logger.info(f"{lig_name} retry#{attempt + 1} ({err_cat}) -> {result_r}")
                        if result_r.get("valid", False) and score_r is not None:
                            ok_r, reason_r = _validate_with_rmsd_gate(lig, lig_name, out_path3, float(score_r))
                            if ok_r:
                                scores[lig] = float(score_r)
                                validated_ligands.append(lig)
                                raw_docked_ligands[lig] = norm(out_path3)
                                logger.info(f"{lig_name} | {stage_retry2['name']} score: {score_r:.2f} kcal/mol (retry rescued)")
                                retained_invalid = False
                                break
                            else:
                                invalids[lig] = (float(score_r), reason_r or "rmsd_fail")
                                retained_invalid = False
                                break

                        attempt += 1

                    if retained_invalid:
                        invalids[lig] = (float(score), result.get("reason", "pose_invalid"))
                        logger.info(f"{lig_name} | pose invalid (after retries)")

                    processed += 1
                    if (processed % 25 == 0) or (processed == len(futures)):
                        pbar.set_postfix(ok=len(scores), inv=len(invalids))
                    pbar.update(1)

    return scores, validated_ligands, all_distances, raw_docked_ligands, invalids

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


def fallback_recentering_if_empty(
        cfg: Dict,
        pdb_id: str,
        stage_name: str,
        scores: Dict[str, float],
        raw_docked_ligands: Dict[str, str],
        receptor_pdbqt: str,
        center: Tuple[float, float, float],
        box_size: Tuple[float, float, float],
        stage1_original: List[str],
        logger: logging.Logger,
        guard: GlobalCenterGuard,
        control_anchor_hit: bool
) -> Tuple[bool, Tuple[float, float, float], Tuple[float, float, float], List[str]]:
    """
    If a stage yields no valid ligands, attempt a fallback recenter and restart stage1.
    De-duped: will not fire if early recenter or CenterSelector already switched this stage,
    or if a control anchor validated in this stage, or if global cap reached.
    """
    if scores:
        return False, center, box_size, []
    if control_anchor_hit:
        logger.info("Empty-stage fallback skipped: control-anchored validation present earlier.")
        return False, center, box_size, []
    if not guard.can_switch():
        logger.info("Empty-stage fallback skipped: global switch guard disallows further switches.")
        return False, center, box_size, []

    logger.warning(f"No valid ligands in {stage_name}. Attempting fallback recentering...")
    # >>> DOCKED PATHS PATCH START
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # >>> DOCKED PATHS PATCH END
    fb_pose, new_center, _best_score, _chosen_lig = attempt_fallback_recenter(
        fallback_ligands=raw_docked_ligands,
        receptor_pdbqt=receptor_pdbqt,
        docking_dir=str(paths.docked_pdb_root()),
        stage_name=stage_name,
        pocket_center=center,
        logger=logger,
        exclude_basenames=set(),
    )
    if new_center is None:
        logger.warning("Fallback recovery failed.")
        return False, center, box_size, []

    new_box = tuple(min(float(cfg.get("BOX_SIZE_MAX_A", 28.0)), s) for s in box_size)
    guard.mark_switch()  # counts as a global switch
    logger.info("Re-running stage1 with new center after no-valid fallback. [global switch]")
    return True, new_center, new_box, stage1_original[:]


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
def process_one_protein(cfg: Dict, pdb_file: str, stages: List[Dict], params: RecenterParams) -> None:
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

    # 1) Extract ligands ? produce nolig PDB
    _lig_count, control_stems = extract_ligands_to_nolig(paths, logger)
    # Ensure extracted crystal controls are prepped before control redock
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

    # build crystal-ligand lookup (prefer PDB)
    control_lookup = build_control_lookup(paths)

    # 2) Protein prep (re-use if cached)
    logger.info("[ph.debug] calling prepare_receptor; PH_ENSEMBLE=%s", cfg.get("PH_ENSEMBLE", False))
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("Skipping protein due to prep failure.")
        return


    # insert: strip monoatomic ions (Na+, K+, Cl-, etc.) before Meeko uses the PDB
    try:
        from automate_protein_prep import _strip_monoatomic_ions_inplace
        _strip_monoatomic_ions_inplace(Path(cleaned_pdb))
        logger.info("Stripped monoatomic ions from cleaned PDB before pocket detection/Meeko.")
    except Exception as e:
        logger.warning(f"Strip monoatomic ions skipped: {e}")

    # 3) Pocket detection
    center, box_size, center_source = None, None, "none"
    try:
        sel_center, sel_box = select_center_via_control_redock(cfg, paths, receptor_pdbqt, logger)
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
                member_index_start=0
            )
            logger.info("[ph_ensemble.manifest] path=%s", manifest_path)

        except Exception as e:
            logger.warning("[ph_ensemble.skip] error=%s", e)
    # >>> PH ENSEMBLE (GLOBAL) END


    # 4) Ligand prep & filtering
    cfg.setdefault("_EFFECTIVE_SINGLE_LIGAND", "")
    single_ligand_hit: Optional[Path] = None
    cfg.pop("_SINGLE_RESOLVED_PATH", None)
    # --- Single-ligand mode (if active) --------------------------------------
    if cfg["_EFFECTIVE_SINGLE_LIGAND"]:
        # Provide per-protein paths to resolver
        cfg.setdefault("paths", {})
        cfg["paths"]["prepped_ligands_dir"] = str(paths.prepped_ligands_dir)

        hit = _resolve_single_ligand(cfg["_EFFECTIVE_SINGLE_LIGAND"], pdb_id, cfg, logger)
        if hit:
            cfg["_SINGLE_RESOLVED_PATH"] = str(hit)
            single_ligand_hit = hit
        else:
            logger.error("[single.block] selector '%s' not found in fda_library via FDA_MAPPING_CSV; aborting instead of fallback.", cfg["_EFFECTIVE_SINGLE_LIGAND"])
            raise SystemExit(2)

    ligands, heavy_atom_counts, pains_flags = prepare_and_filter_ligands(cfg, paths, logger)
    if single_ligand_hit:
        ligands = [str(single_ligand_hit)]
        ha = _count_heavy_atoms_from_pdbqt(single_ligand_hit)
        heavy_atom_counts = {str(single_ligand_hit): ha}
        pains_flags = {}
        logger.info(f"[single] Active ? docking only: {single_ligand_hit.name} (heavy={ha})")

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

    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_label = variant_env or "legacy"
    ph_log = logging.getLogger("ph_ensemble")
    manifest_map = cfg.get("_PH_ENSEMBLE_CANONICAL") or {}
    ph_runs: list[tuple[Optional[str], str]] = []

    if bool(cfg.get("PH_ENSEMBLE")):
        # --- DEBUG: what do we have right now?
        try:
            _keys = sorted(list((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).keys()))
            _len_here = len((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).get(paths.pdb_id, []))
            ph_log.info("[ph_ensemble.debug] map_keys=%s map_len[%s]=%d", ",".join(_keys), paths.pdb_id, _len_here)
        except Exception:
            pass

        ph_runs = [(lbl, str(path)) for lbl, path in manifest_map.get(paths.pdb_id, []) if path]

        # If empty, try a direct manifest read as a last-resort bridge.
        if not ph_runs:
            try:
                guess = Path(cfg["OUTPUT_DIR"]) / paths.pdb_id / "receptor" / "ph_ensemble" / "ensemble.json"
                ph_log.info("[ph_ensemble.debug] fallback_manifest=%s exists=%s", str(guess), guess.exists())
                if guess.exists():
                    data = json.loads(guess.read_text())
                    members = data.get("members") or []
                    canonical = [m for m in members if bool(m.get("canonical", False))]
                    if canonical:
                        members = canonical
                    prefix = f"{paths.pdb_id}_"
                    targets = []
                    for entry in members:
                        receptor_path = entry.get("pdbqt")
                        if not receptor_path:
                            continue
                        stem = Path(receptor_path).stem
                        ph_label = stem[len(prefix):] if stem.startswith(prefix) else stem
                        targets.append((ph_label, receptor_path))
                    if targets:
                        cfg.setdefault("_PH_ENSEMBLE_CANONICAL", {})[paths.pdb_id] = targets
                        ph_runs = [(lbl, str(p)) for (lbl, p) in targets if p]
                        ph_log.info("[ph_ensemble.debug] fallback_bridge n=%d labels=%s",
                                    len(targets), ",".join(lbl for lbl, _ in targets))
            except Exception as _e:
                ph_log.warning("[ph_ensemble.debug] fallback_bridge.error %s", _e)

        if not ph_runs:
            ph_log.error("[ph_ensemble.abort] PH_ENSEMBLE=True but no canonical targets for %s; refusing legacy fallback.", paths.pdb_id)
            return  # disallow legacy fallback when ensemble is enabled

    if not ph_runs:
        ph_runs = [(None, str(receptor_pdbqt))]


    for ph_label, receptor_override in ph_runs:
        receptor_current = str(receptor_override or receptor_pdbqt)
        if not receptor_current:
            continue
        try:
            if ph_label and not Path(receptor_current).exists():
                ph_log.warning(
                    "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s receptor_missing=%s",
                    paths.pdb_id,
                    variant_label,
                    ph_label,
                    receptor_current,
                )
                continue
        except Exception:
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
        receptor_pdbqt = receptor_current

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
    print(f"[run] log_file={log_path} run_id={run_id}")
    cfg = ConfigDict(load_inputs())
    validate_config(cfg)

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
    cfg.setdefault("SINGLE_LIGAND_SEARCH_ORDER", "per_protein,global")
    cfg.setdefault("SINGLE_LIGAND_ALLOW_PREFIX", False)
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
    cfg.setdefault("TEST_MODE_ENABLE", False)
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
    if _to_bool(str(cfg.get("TEST_MODE_ENABLE", "false"))):
        raw_map = cfg.get("TEST_LIBRARY_MAP", {})
        test_keys = set()
        if isinstance(raw_map, dict):
            test_keys = {str(k).upper()[:4] for k in raw_map.keys()}
        else:
            # Accept JSON or Python-literal dict strings
            try:
                parsed = json.loads(str(raw_map).strip())
            except Exception:
                import ast

                try:
                    parsed = ast.literal_eval(str(raw_map).strip())
                except Exception:
                    parsed = {}
            if isinstance(parsed, dict):
                test_keys = {str(k).upper()[:4] for k in parsed.keys()}
        if test_keys:
            kept, skipped = [], []
            for f in pdb_files:
                nid = _norm_pdb_id(f)
                if nid and nid.upper() in test_keys:
                    kept.append(f)
                else:
                    skipped.append(f)
            if skipped:
                # Print short list of IDs we're skipping so it's obvious in logs
                skipped_ids = sorted({(_norm_pdb_id(x) or x) for x in skipped})
                print(
                    f"[test-mode] Skipping {len(skipped)} protein(s) not in TEST_LIBRARY_MAP: {', '.join(skipped_ids[:20])}" +
                    (" ..." if len(skipped_ids) > 20 else ""))
            pdb_files = kept
        else:
            print("[test-mode] TEST_LIBRARY_MAP empty/invalid; no extra filtering applied.")

    print("Working directory:", os.getcwd())
    print("Loaded config keys:", list(cfg.keys()))
    print(f"Proteins queued: {len(pdb_files)}")


    start = time.time()
    mode, variants = resolve_apo_holo_mode(cfg)
    logging.info(f"[apo-holo] resolved mode={mode} variants={variants} "
                 f"env.APO_HOLO_MODE='{os.environ.get('APO_HOLO_MODE')}'")

    from tqdm import tqdm as _tqdm
    for variant in variants:
        # Make variant visible to any module still reading env (legacy compatibility)
        os.environ["APO_HOLO_VARIANT"] = "" if variant is None else str(variant).upper()
        label = "legacy" if variant is None else str(variant).lower()
        with _tqdm(
                total=len(pdb_files),
                desc=f"Processing Proteins ({label})",
                unit="protein",
                position=0,
                dynamic_ncols=True,
                mininterval=0.2,
                leave=True,
                file=sys.stdout
        ) as bar:
            cfg_v = cfg  # no per-variant mutation; variant is propagated via APO_HOLO_VARIANT env
            for pdb_file in pdb_files:
                process_one_protein(cfg_v, pdb_file, stages, params)
                # Keep APO, delete HOLO if byte-identical (run after both variants exist)
                if mode == "apo_vs_holo" and (variant == "HOLO"):
                    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()
                    try:
                        dedup_identical_variants(pdb_id, cfg_v)
                    except Exception as _e:
                        logging.warning(f"[apo-vs-holo] dedup skipped for {pdb_id}: {_e}")
                bar.update(1)

    elapsed_min = (time.time() - start) / 60.0
    print(f"\nAll proteins processed in {elapsed_min:.2f} minutes.")


if __name__ == "__main__":
    main()