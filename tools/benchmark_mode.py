# benchmark_mode.py
"""
Benchmark mode: single ultra-stage docking of likely co-crystal FDA ligands, per pocket.

This version:
- Resolves INPUT/OUTPUT/PREPPED from YAML config or CLI; no Windows paths.
- Resolves FDA mapping CSV from CLI -> YAML (FDA_MAPPING_CSV) -> alongside library
  -> chemdb/data/fda_mapping_from_pdbqt.csv under OVERALL_DIR or CWD.
- Sources CHEMCOMP_ALIAS / EXCLUDE_* / PER_PDB_HINTS / HARD_FDA_CONTROL_BY_PDB from chemdb.chem_alias_db
- Keeps optional alias augmentation from prior benchmark details CSVs
- Uses  logic from main.py for docking/validation/switching
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import threading
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import docking.capture_pose as capture_pose
from docking.capture_pose import (
    _render_native_on_original_pdb,
    _render_three_views_with_pymol,
    _safe_open_csv_for_write,
    pick_control_and_nearest_rdk,
)
# Per-run timestamp (exported later after RUN_DIR is set)
RUN_START_EPOCH: float = float(time.time())

def _bool_env(name: str, default: bool) -> bool:
    v = str(os.environ.get(name, "")).strip().lower()
    if not v:
        return default
    return v in ("1", "true", "yes", "on")

def _render_only_current_run_enabled(cfg: dict) -> bool:
    # Default true: pre-docking render enqueue is disabled
    return _bool_env("BENCH_RENDER_ONLY_CURRENT_RUN", True)

from concurrent.futures import ThreadPoolExecutor, as_completed

# Config
from config.runtime_config import load_inputs, validate_config, load_config
from config.output_paths import output_root
from cli.run_bootstrap_helpers import init_config_run_dir
import logging
import sys
cfg = load_config()
# Library resolution helpers
from prep_ligands.metabolite_resolver import (
    ensure_parent_drugs_for_controls,
    load_library_index,
    resolve_corresponding_name_for_rdk,
    resolve_corresponding_name_from_text,
)
try:
    from tools.fda_mapping_index import MappingIndex, MappingRow
except Exception:
    from fda_mapping_index import MappingIndex, MappingRow

# Fallback active-site detector
from docking.active_site_detection import detect_active_site

# Unified chem/alias config (YAML-backed)
from chemdb.chem_alias_db import (
    CHEMCOMP_ALIAS,
    EXCLUDE_HET_IDS,
    EXCLUDE_HET_NAME_KEYWORDS,
    HARD_FDA_CONTROL_BY_PDB,
    PER_PDB_HINTS,
)

# Import core pipeline pieces (updated to new package paths)
from docking.fallback_recenter import BudgetGuard, GlobalCenterGuard
from docking.docking_centering import CenterSelector
from docking.docking_controls import build_control_lookup, detect_pocket, extract_ligands_to_nolig
from docking.checkpoints import checkpoint_invalidate_from, checkpoint_mark_done, checkpoint_should_skip
from docking.docking_utils import final_pose_validation_and_screenshots, _fingerprint_stage, early_recenter_decision
from main import get_recenter_params
from cli.runtime_logging import make_protein_logger
from docking.docking_receptor import prepare_receptor
from docking.ligand_metrics import record_le
from docking.score_io import record_score
from docking.docking_stage_runner import run_one_stage, RetryManager
from docking.docking_ligands import _count_heavy_atoms_from_pdbqt
# >>> PATHS IMPORT START
from path_router import make_paths, expand_variants
# >>> PATHS IMPORT END
from prep_ligands.prep_ligands import prep_ligands_from_pdb
# ----- Deferred render queue (global) -----
from contextlib import contextmanager

# --- [AUTO-ENV/LOG GUARD for micromamba + tee logging] ---
import shutil
import datetime
import atexit

TARGET_ENV = "docking-env"

def _in_target_env():
    # Treat either micromamba or conda naming as "already in docking-env"
    return (os.environ.get("MAMBA_DEFAULT_ENV") == TARGET_ENV) or \
           (os.environ.get("CONDA_DEFAULT_ENV") == TARGET_ENV)

def _prepare_logfile_from_cwd():
    logs_dir = output_root(Path(os.getcwd()), "logs")
    logs_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    return str(logs_dir / ("bench_%s.log" % ts))

class _Tee(object):
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
        try: self._stream.flush()
        except Exception: pass
        try: self._fh.flush()
        except Exception: pass
    def isatty(self):
        try: return self._stream.isatty()
        except Exception: return False
    def close(self):
        try: self._fh.close()
        except Exception: pass

def _tee_stdio_to(log_path):
    _orig_stdout, _orig_stderr = sys.stdout, sys.stderr
    tee_out, tee_err = _Tee(_orig_stdout, log_path), _Tee(_orig_stderr, log_path)
    sys.stdout, sys.stderr = tee_out, tee_err
    def _announce_and_close():
        try:
            sys.stdout.write("Log -> %s\n" % log_path)
            sys.stdout.flush()
        finally:
            try:
                tee_out.close()
                tee_err.close()
            except Exception:
                pass
    atexit.register(_announce_and_close)

# Require Python 3+ to parse this file at all.
if sys.version_info[0] < 3:
    sys.stderr.write("[FATAL] This program requires Python 3+. Your interpreter is: %s\n" % sys.version.split()[0])
    sys.stderr.write("Use `python3 benchmark_mode.py` or `micromamba run -n %s python -u benchmark_mode.py ...`\n" % TARGET_ENV)
    sys.exit(2)

# If not in docking-env, re-exec via micromamba (no conda fallback).
if not _in_target_env():
    mm = shutil.which("micromamba")
    if not mm:
        sys.stderr.write("[FATAL] micromamba not found on PATH; cannot enter '%s'.\n" % TARGET_ENV)
        sys.stderr.write("Hint: ensure micromamba is on PATH, then try again.\n")
        sys.exit(2)
    log_path = _prepare_logfile_from_cwd()
    env = os.environ.copy()
    env["ATLAS_LOG_FILE"] = log_path
    script_path = os.path.abspath(__file__)
    cmd = [mm, "run", "-n", TARGET_ENV, "python", "-u", script_path] + sys.argv[1:]
    os.execvpe(cmd[0], cmd, env)

# Already inside docking-env: enable tee logging (use provided path if any).
_log_path = os.environ.get("ATLAS_LOG_FILE") or _prepare_logfile_from_cwd()
os.environ["ATLAS_LOG_FILE"] = _log_path
_tee_stdio_to(_log_path)
# --- [end AUTO-ENV/LOG GUARD] ---



# --- PyMOL deferral & parallelism bootstrap
_DEFER_DEFAULT = str(os.environ.get("DEFER_PYMOL", "1")).lower() in ("1","true","yes","on")
capture_pose.set_defer_mode(_DEFER_DEFAULT, queue_path=os.environ.get("PYMOL_DEFER_QUEUE"))
_RENDER_LOCK = threading.Semaphore(int(os.environ.get("PYMOL_PARALLEL", "1")))
def _parse_image_variants(cfg) -> tuple[set[str], bool]:
    """
    Returns (variants, is_subset_mode)

    variants: normalized set in {"ctrl", "rdk", "ctrl+rdk"}.
    is_subset_mode: True iff the user explicitly set PYMOL_IMAGE_VARIANTS
                    (i.e., filter strictly to those; otherwise keep back-compat extras).
    """
    raw = (os.environ.get("PYMOL_IMAGE_VARIANTS")
           or str(cfg.get("PYMOL_IMAGE_VARIANTS") or "")).strip().lower()

    def _bool_from_cfg_env(key: str, default: bool) -> bool:
        v = (os.environ.get(key) or str(cfg.get(key) or "")).strip().lower()
        if v in ("", None):
            return default
        return v not in ("0", "false", "no", "off")

    include_native = _bool_from_cfg_env("PYMOL_IMAGE_INCLUDE_NATIVE", False)

    if not raw:
        # Back-compat default: render everything we render today
        return {"ctrl", "rdk", "ctrl+rdk"}, False
    toks = {t.strip() for t in raw.split(",") if t.strip()}
    valid = {"ctrl", "rdk", "ctrl+rdk"}
    return (toks & valid) or set(), True


@contextmanager
def _acquire(sem: threading.Semaphore):
    sem.acquire()
    try:
        yield
    finally:
        sem.release()




# ---- Re-entrant patch for docking ThreadPoolExecutor ----
import threading
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor as _RealTPE

_GLOBAL_LIGAND_SEM = None  # unchanged
class _GuardedTPE(_RealTPE):
    def submit(self, fn, *args, **kwargs):
        sem = _GLOBAL_LIGAND_SEM
        if sem is None:
            return super().submit(fn, *args, **kwargs)
        def _wrapped(*a, **k):
            with _acquire(sem):
                return fn(*a, **k)
        return super().submit(_wrapped, *args, **kwargs)

# re-entrant state
_TP_PATCH_LOCK   = threading.Lock()
_TP_PATCH_COUNT  = 0
_TP_PREV_EXEC    = {}

@contextmanager
def _patch_main_threadpool(sem):
    """
    Temporarily replace the executor used by docking stage submissions with a
    guarded one.
    Re-entrant and thread-safe: multiple overlapping uses are OK.
    """
    import docking.docking_stage_runner_support as _stage_runner_support
    import main as _main_mod
    global _GLOBAL_LIGAND_SEM, _TP_PATCH_COUNT, _TP_PREV_EXEC
    with _TP_PATCH_LOCK:
        _GLOBAL_LIGAND_SEM = sem
        if _TP_PATCH_COUNT == 0:
            _TP_PREV_EXEC = {
                "benchmark_mode": globals().get("ThreadPoolExecutor"),
                "main": getattr(_main_mod, "ThreadPoolExecutor", None),
                "stage_runner_support": getattr(
                    _stage_runner_support,
                    "ThreadPoolExecutor",
                    None,
                ),
            }
            globals()["ThreadPoolExecutor"] = _GuardedTPE
            setattr(_main_mod, "ThreadPoolExecutor", _GuardedTPE)
            setattr(_stage_runner_support, "ThreadPoolExecutor", _GuardedTPE)
        _TP_PATCH_COUNT += 1
    try:
        yield
    finally:
        with _TP_PATCH_LOCK:
            _TP_PATCH_COUNT -= 1
            if _TP_PATCH_COUNT == 0:
                globals()["ThreadPoolExecutor"] = _TP_PREV_EXEC.get("benchmark_mode")
                setattr(_main_mod, "ThreadPoolExecutor", _TP_PREV_EXEC.get("main"))
                setattr(
                    _stage_runner_support,
                    "ThreadPoolExecutor",
                    _TP_PREV_EXEC.get("stage_runner_support"),
                )
                _GLOBAL_LIGAND_SEM = None
                _TP_PREV_EXEC = {}


from docking import capture_pose as _cap
from contextlib import contextmanager

_DEFERRED_PYMOL_CALLS = []
_ORIG_RENDER_THREE = getattr(_cap, "_render_three_views_with_pymol")
_ORIG_RENDER_NATIVE = getattr(_cap, "_render_native_on_original_pdb")

_PYMOL_DEFER_LOCK  = threading.Lock()
_PYMOL_DEFER_COUNT = 0



def run_deferred_captures(max_workers: int):
    print("[render] WARNING: legacy mid-run capture replay is disabled. Use capture_pose.replay_deferred_jobs_mp().")



def _enqueue_render_task(
    *,
    pdb_id: str,
    stage_name: str,
    paths,
    variant: Optional[str],
    cleaned_pdb_path: str,
    original_pdb_path: str,
    ctrl_pose_path: str | None,
    rdk_pose_path: str | None,
    exclude_resns: List[str],
) -> None:
    """
    Unified producer: enqueue directly into capture_pose's JSONL by calling its APIs with DEFER_PYMOL=1.
    """
    from docking import capture_pose as _cap
    stage_dir_target = paths.docked_stage_dir(variant, stage_name)

    # Viewport first (used by native + three-views)
    viewport_w = int(cfg.get("VIEWPORT_W", 640))
    viewport_h = int(cfg.get("VIEWPORT_H", 480))
    viewport = (viewport_w, viewport_h)

    # Resolve which three-view sets to generate
    _variants, _subset = _parse_image_variants(cfg)
    include_native = str(os.environ.get("PYMOL_IMAGE_INCLUDE_NATIVE", "") or
                         str(cfg.get("PYMOL_IMAGE_INCLUDE_NATIVE", ""))).strip().lower() in ("1", "true", "yes", "on")

    # Native view
    if include_native and original_pdb_path:
        _cap._render_native_on_original_pdb(
            original_pdb=original_pdb_path,
            outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-native__NATIVE"),
            exclude_resns=exclude_resns,
            viewport=viewport,
        )

    try:
        label_top_n_res = int(os.environ.get("PYMOL_LABEL_TOP_N_RES") or 5)
    except (TypeError, ValueError):
        label_top_n_res = 5

    try:
        label_cutoff = float(os.environ.get("PYMOL_LABEL_CUTOFF_A") or 5.0)
    except (TypeError, ValueError):
        label_cutoff = 5.0
    viewport_w = int(cfg.get("VIEWPORT_W", 640))
    viewport_h = int(cfg.get("VIEWPORT_H", 480))
    viewport = (viewport_w, viewport_h)
    # Three-view(s)
    # ctrl-only (cleaned receptor + control)
    if ctrl_pose_path and ("ctrl" in _variants):
        _cap._render_three_views_with_pymol(
            receptor_path=cleaned_pdb_path,
            ligand_paths_and_colors=[(ctrl_pose_path, "control", "green")],
            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL"),
        )

    # rdk-only (cleaned receptor + rdk)
    if rdk_pose_path and ("rdk" in _variants):
        _cap._render_three_views_with_pymol(
            receptor_path=cleaned_pdb_path,
            ligand_paths_and_colors=[(rdk_pose_path, "rdk", "magenta")],
            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__RDKclosest"),
        )

    # Historical extra: rdk on original receptor (kept only for back-compat when user didn’t set PYMOL_IMAGE_VARIANTS)
    if rdk_pose_path and not _subset:
        _cap._render_three_views_with_pymol(
            receptor_path=original_pdb_path,
            ligand_paths_and_colors=[(rdk_pose_path, "rdk", "magenta")],
            outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-rdk__RDKclosest"),
        )

    # ctrl+rdk (cleaned receptor + both ligands)
    if ctrl_pose_path and rdk_pose_path and ("ctrl+rdk" in _variants):
        _cap._render_three_views_with_pymol(
            receptor_path=cleaned_pdb_path,
            ligand_paths_and_colors=[
                (ctrl_pose_path, "control", "green"),
                (rdk_pose_path, "rdk", "magenta"),
            ],
            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL+RDKclosest"),
        )


def run_deferred_renders(max_workers: int) -> None:
    if not _DEFERRED_RENDER_TASKS:
        print("[render] nothing to render (queue empty).")
        return
    print(f"[render] starting deferred renders: {len(_DEFERRED_RENDER_TASKS)} tasks | workers={max_workers}")
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futs = [pool.submit(_execute_render_task, t) for t in _DEFERRED_RENDER_TASKS]
        for fut in as_completed(futs):
            try:
                fut.result()
            except Exception as e:
                print(f"[render] task failed: {e}")
    print("[render] all deferred renders complete.")


def _resolve_receptor_for(pdb_id: str, cfg: dict) -> Optional[Path]:
    """
    Locate the cleaned receptor PDB for <pdb_id>.
    Canonical: processed_pdbs/<PDB>/receptor/<PDB>_cleaned.pdb
    Legacy fallback (read-only): processed_pdbs/<PDB>_NOLIG/receptor/*_NOLIG_cleaned.pdb (or *cleaned.pdb)
    """
    base = Path(cfg.get("OUTPUT_DIR") or Path(cfg["OVERALL_DIR"]).parent / "processed_pdbs").resolve()
    canonical = base / pdb_id / "receptor" / f"{pdb_id}_cleaned.pdb"
    if canonical.exists():
        return canonical.resolve()

    # Back-compat (read-only): accept legacy sibling if present
    legacy_dir = base / f"{pdb_id}_NOLIG" / "receptor"
    for pat in (f"{pdb_id}_NOLIG_cleaned.pdb", "*_NOLIG_cleaned.pdb", "*cleaned.pdb"):
        hits = sorted(legacy_dir.glob(pat))
        if hits:
            return hits[0].resolve()

    return None

def _resolve_allowed_rdk_ids(mapping: MappingIndex, names: Sequence[str]) -> set[str]:
    """Resolve allow-list of drug names to concrete RDK ids via mapping CSV."""
    ids: set[str] = set()
    base: list[str] = [n for n in (names or []) if n]
    SYN = {  # minimal shim; extend as needed
        "venetoclax": ["venclexta", "abt-199", "abt199"],
    }
    for nm in base:
        q = [nm]
        extra = SYN.get(str(nm).strip().lower())
        if extra:
            q.extend(extra)
        for row, _sc, _why in mapping.search(q, inchikey=None, max_results=6):
            rid = _extract_rdk_id(Path(row.path).name)
            if rid:
                ids.add(rid.lower())
    return ids


def _find_control_files(prepped_dir: Path, control_names: list[str]) -> list[str]:
    """Return full paths to existing .pdbqt files matching the requested control base names."""
    found = []
    for name in control_names:
        # Match typical variants (prefix, suffix, stage suffixes); keep existing selection heuristics downstream.
        for p in prepped_dir.glob(f"{name}*.pdbqt"):
            try:
                if p.is_file():
                    found.append(str(p))
            except Exception:
                pass
    # Stable order helps reproducibility
    return sorted(set(found))

# =============================
# Small text/path helper utils
# =============================

def _remap_to_scopes(p: Path, scope_dirs: Sequence[Path]) -> Optional[Path]:
    """
    Try to relocate file p into one of the chosen scope directories by:
      (1) exact basename match, else
      (2) rdk_XXXXXX core match allowing suffixes (e.g., *_stage5).
    Returns the first in-scope hit or None if not found.
    """
    rid = _extract_rdk_id(p.name)  # may be None
    for root in scope_dirs:
        # (1) exact basename
        q = root / p.name
        if q.is_file():
            return q.resolve()
        # (2) rdk id core match
        if rid:
            hits = sorted(root.glob(f"*{rid}*.pdbqt"))
            if hits:
                return hits[0].resolve()
    return None

def _any_within(child: Path, parents: Sequence[Path]) -> bool:
    return any(_path_is_within(child, pr) for pr in parents)









def _collapse_name_once(p: Path) -> Path:
    # turn "...sanitized.sanitized.xyz" into "...sanitized.xyz"
    return p.with_name(re.sub(r'(?:\.sanitized)+(?=\.)', '.sanitized', p.name))
def _norm(s: Optional[str]) -> str:
    """Lowercase, trim, collapse spaces, and strip punctuation-like separators."""
    if s is None:
        return ""
    s = str(s).strip().lower()
    s = s.replace("\u00A0", " ")
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[-_/\\,;:\|\[\]\(\)\{\}\.\+\*'\"]+", " ", s)
    return s

def _tokenize(s: str) -> List[str]:
    """Tokenize to [a-z0-9+] words used by search/ranking."""
    if not s:
        return []
    toks = re.split(r"[^a-z0-9\+]+", s.lower())
    return [t for t in toks if t]

def _extract_rdk_id(text: str) -> Optional[str]:
    """Extract rdk_XXXXXX id from a filename-like string."""
    if not text:
        return None
    m = re.search(r"(rdk_\d{6,8})", Path(text).stem.lower())
    return m.group(1) if m else None

def _norm_text(s: Optional[str]) -> str:
    """Letters+digits only, for approximate comparisons."""
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())

def _looks_like_brand_or_generic(name: str) -> bool:
    n = (name or "").strip()
    if len(n) < 3 or len(n) > 64:
        return False
    if not re.search(r"[a-zA-Z]", n):
        return False
    noisy = len(re.findall(r"[0-9\[\]\(\)\/\.\-\,;:]", n))
    return (noisy / max(1, len(n))) < 0.40

def _clean_alias(name: str) -> Optional[str]:
    """Light cleanup and guardrails for text harvested from CSVs."""
    if not name:
        return None
    s = name.strip()
    if s.startswith("?"):
        s = s.lstrip("?\uFF1F").strip()
    s = s.replace("\u2019", "'").replace("\u00AE", "").replace("\u2122", "")
    s = re.sub(r"\s+", " ", s)
    return s.lower() if _looks_like_brand_or_generic(s) else None

def _expand_globs(paths: Sequence[str]) -> List[Path]:
    out: List[Path] = []
    for p in paths:
        if any(ch in p for ch in "*?[]"):
            out.extend(Path(x) for x in glob.glob(p))
        else:
            out.append(Path(p))
    return out

def _dedupe_str(seq: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out

def _path_is_within(child: Path, parent: Path) -> bool:
    """Return True if `child` is within `parent` (no strict existence requirement)."""
    try:
        child = child.resolve(strict=False)
        parent = parent.resolve(strict=False)
        child.relative_to(parent)
        return True
    except Exception:
        return False

# =============================
# PDB parsing / hint extraction
# =============================

def parse_pdb_het_hints(pdb_path: Path) -> Tuple[List[str], List[str]]:
    """
    Parse HET / HETNAM / HETATM records from a PDB to extract plausible ligand HET codes
    and human-readable names, applying EXCLUDE_* filters.
    """
    pdb_path = Path(pdb_path)
    het_ids: set[str] = set()
    hetnam_map: Dict[str, List[Tuple[int, str]]] = defaultdict(list)

    if not pdb_path.is_file():
        return [], []

    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                rec = raw[:6].strip().upper()

                if rec == "HET":
                    het = (raw[7:10].strip() if len(raw) >= 10 else "").upper()
                    if not het:
                        toks = raw.split()
                        if len(toks) >= 2:
                            het = toks[1].upper()
                    if het:
                        het_ids.add(het)

                elif rec == "HETNAM":
                    toks = raw.split()
                    cont = 1
                    het, name = "", ""
                    if len(toks) >= 3 and toks[1].isdigit():
                        cont = int(toks[1]); het = toks[2].upper(); name = " ".join(toks[3:])
                    elif len(toks) >= 2:
                        het = toks[1].upper(); name = " ".join(toks[2:])
                    if het and name:
                        hetnam_map[het].append((cont, name.strip()))

                elif rec == "HETATM":
                    resn = raw[17:20].strip().upper()
                    if resn:
                        het_ids.add(resn)
    except Exception:
        return [], []

    STANDARD_RES = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
        "PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"
    }

    het_ids = {
        h for h in het_ids
        if h not in EXCLUDE_HET_IDS and h not in STANDARD_RES and 2 <= len(h) <= 5
    }

    het_names: List[str] = []
    for het in sorted(het_ids):
        if het in hetnam_map:
            parts = [t for _, t in sorted(hetnam_map[het], key=lambda x: x[0])]
            nm = " ".join(parts).strip()
            if not nm:
                continue
            nm_up = nm.upper()
            if any(kw in nm_up for kw in EXCLUDE_HET_NAME_KEYWORDS):
                continue
            cleaned = _clean_alias(nm)
            het_names.append(cleaned if cleaned else nm)

    # dedupe preserving order
    seen = set()
    dedup_names: List[str] = []
    for n in het_names:
        nn = _norm(n)
        if nn and nn not in seen:
            seen.add(nn)
            dedup_names.append(n)

    return sorted(het_ids), dedup_names

# =============================
# Alias augmentation (details CSV)
# =============================

def load_aliases_from_details_csv(paths: Sequence[str]) -> Dict[str, List[str]]:
    """
    Harvest (HET -> alias) from benchmark_analysis_details.csv files.
    Returns aliases lowercased and deduped per HET code.
    """
    exp: Dict[str, set] = {}
    for csv_path in _expand_globs(paths):
        if not csv_path.is_file():
            continue
        try:
            with open(csv_path, "r", encoding="utf-8", errors="ignore", newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    het = (row.get("control_het") or "").strip().upper()
                    nm = _clean_alias(row.get("rdk_name") or "")
                    if het and nm:
                        exp.setdefault(het, set()).add(nm)
        except Exception as e:
            print(f"[alias] WARN: failed to parse {csv_path}: {e}")
    return {k: sorted(v) for k, v in exp.items()}

def merge_aliases_into_chemcomp(base: Dict[str, List[str]], extra: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """Merge extra aliases into an existing chemcomp map. Output values are lowercased & deduped."""
    out: Dict[str, List[str]] = {k: sorted(set(vv.lower() for vv in vs)) for k, vs in base.items()}
    for het, names in (extra or {}).items():
        merged = set(out.get(het, []))
        for nm in names:
            nm_norm = (nm or "").strip().lower()
            if nm_norm:
                merged.add(nm_norm)
        out[het] = sorted(merged)
    return out

# =============================
# Selection / orchestration
# =============================

def _resolve_name_for_path_or_text(p_or_text: str, fda_index) -> str:
    """Prefer RDK-based resolution; fallback to filename text matching."""
    rdk = _extract_rdk_id(p_or_text or "")
    if rdk:
        name = resolve_corresponding_name_for_rdk(rdk, fda_index) or ""
        if name:
            return name
    stem = Path(p_or_text).stem if os.path.exists(p_or_text) else str(p_or_text)
    name = resolve_corresponding_name_from_text(stem, fda_index) or ""
    return name

def _remap_to_prepped(p: Path, prepped_dir: Path) -> Path:
    """If p isn't within prepped_dir, try to map it to a file that *is* there (by basename or rdk_* id)."""
    if p.is_file() and p.resolve().parent == prepped_dir.resolve():
        return p
    # 1) Same basename inside prepped
    q = prepped_dir / p.name
    if q.is_file():
        return q
    # 2) rdk_* id match inside prepped
    rid = _extract_rdk_id(p.name)
    if rid:
        hits = list(prepped_dir.glob(f"*{rid}*.pdbqt"))
        if hits:
            return hits[0]
    return p  # unchanged if we couldn't remap
def select_candidates_for_protein(
    mapping: MappingIndex,
    hints: Sequence[str],
    fda_scope_dirs: Sequence[Path],
    extra_parent_ids: Sequence[str],
    max_candidates: int,
) -> List[Tuple[MappingRow, int, str]]:
    """
    Build the FDA/RDK candidate list restricted to the union of fda_scope_dirs.
    Absolute paths in the mapping CSV are treated as hints and remapped into the scope.
    """
    scope_dirs = [Path(d).resolve() for d in fda_scope_dirs or []]
    candidates = mapping.search(hints, inchikey=None, max_results=max(6, max_candidates * 4))

    cand_rows: List[Tuple[MappingRow, int, str]] = []
    dropped_not_exist = 0
    dropped_out_scope = 0

    # Expand parent (metabolite->parent) hints
    parent_rows: List[Tuple[MappingRow, int, str]] = []
    for pid in extra_parent_ids or []:
        rid = _extract_rdk_id(pid or "")
        if rid:
            parent_rows.extend((r, 99, "parent_from_metabolite") for r in mapping.rows_by_rdk_id(rid))
        else:
            for (r, sc, why) in mapping.search([pid], inchikey=None, max_results=2):
                parent_rows.append((r, max(sc, 92), f"{why}|parent_from_metabolite"))

    # Best-per-path before fencing
    best_by_path: Dict[str, Tuple[MappingRow, int, str]] = {}
    for row, sc, why in (candidates + parent_rows):
        p_abs = str(Path(row.path).resolve())
        if (p_abs not in best_by_path) or (sc > best_by_path[p_abs][1]):
            best_by_path[p_abs] = (row, sc, why)

    # Fence to scope and remap
    fenced: List[Tuple[MappingRow, int, str]] = []
    for row, sc, why in best_by_path.values():
        p = Path(row.path)
        # Use as-is if in-scope and exists
        if p.is_file() and _any_within(p.resolve(), scope_dirs):
            fenced.append((row, sc, why))
            continue
        # Otherwise, attempt remap by basename/rdk id into the scope
        remapped = _remap_to_scopes(p, scope_dirs)
        if remapped and remapped.is_file() and _any_within(remapped, scope_dirs):
            new_row = MappingRow(
                path=str(remapped),
                display_name=row.display_name,
                generic_name=row.generic_name,
                brand_names=row.brand_names,
                pubchem_name=row.pubchem_name,
                pubchem_record_title=row.pubchem_record_title,
                pubchem_iupac_name=row.pubchem_iupac_name,
                pubchem_synonyms=row.pubchem_synonyms,
                rxnorm_generic_name=row.rxnorm_generic_name,
                rxnorm_brand_names=row.rxnorm_brand_names,
                drugcentral_generic_name=row.drugcentral_generic_name,
                drugcentral_brand_names=row.drugcentral_brand_names,
                remark_name=row.remark_name,
                sdf_title=row.sdf_title,
                inchikey=row.inchikey,
            )
            fenced.append((new_row, sc, why))
            continue

        # Count drop reason
        if not p.exists():
            dropped_not_exist += 1
        else:
            dropped_out_scope += 1

    # Sort and cap
    fenced.sort(key=lambda t: t[1], reverse=True)
    kept = fenced[:max_candidates]

    # Debug summary (one line)
    try:
        scope_txt = ", ".join(str(s) for s in scope_dirs)
        print(f"[DEBUG] select_candidates scope=[{scope_txt}] | kept={len(kept)} "
              f"(dropped: not_exist={dropped_not_exist}, out_of_scope={dropped_out_scope})")
    except Exception:
        pass

    return kept


def split_controls_and_whitelist(
    whitelist_paths: Sequence[str],
    prepped_control_pdbqts: Sequence[str],
    control_stems_lower: set,
    heavy_atom_counts: Dict[str, int],
    cfg: Dict,
    logger,
) -> Tuple[List[str], List[str]]:
    """Divide ligands into control vs non-control (with caps and HA minimum)."""
    ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
    min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

    def _is_control_path(p: str) -> bool:
        stem = Path(p).stem.split("_stage")[0]
        if stem.upper() in ctrl_blacklist:
            return False
        if stem.lower() not in control_stems_lower:
            return False
        ha = heavy_atom_counts.get(p)
        return (ha is None) or (ha >= min_ha)

    ctrls = [p for p in prepped_control_pdbqts if _is_control_path(p)]
    non_ctrls = [p for p in whitelist_paths if not _is_control_path(p)]
    ctrls = _dedupe_str(ctrls)
    non_ctrls = _dedupe_str(non_ctrls)

    max_ctrls = int(cfg.get("BENCH_MAX_CONTROLS", 8))
    if len(ctrls) > max_ctrls:
        if logger:
            logger.warning(f"[SANITY] controls={len(ctrls)} looks high; capping to {max_ctrls}. Offenders will be skipped.")
            for p in ctrls[max_ctrls:]:
                logger.warning(f"[SANITY-OFFENDER] {p}")
        ctrls = ctrls[:max_ctrls]
    return ctrls, non_ctrls

def _merge_per_pdb_hints(pdb_id: str, hints: List[str]) -> List[str]:
    """Union (deduped) of existing hints with per-PDB and hard-control hints."""
    out = list(hints or [])
    add: List[str] = []
    pid = pdb_id.upper()
    if pid in PER_PDB_HINTS:
        add.extend(PER_PDB_HINTS[pid])
    if pid in HARD_FDA_CONTROL_BY_PDB:
        add.extend(HARD_FDA_CONTROL_BY_PDB[pid])
    seen = set()
    merged: List[str] = []
    for h in (out + add):
        hn = _norm(h)
        if hn and hn not in seen:
            seen.add(hn)
            merged.append(h)
    return merged

def _promote_forced_controls_to_ctrls(
    pdb_id: str,
    mapping: MappingIndex,
    prepped_dir: Path,
    names: Sequence[str],
    logger=None,
) -> Tuple[List[str], List[str]]:
    """Promote explicitly requested controls (from HARD_FDA_CONTROL_BY_PDB) to the control pool."""
    if not names:
        return [], []
    rows = mapping.search(names, inchikey=None, max_results=max(8, len(names) * 3))
    prepped_dir_resolved = prepped_dir.resolve()
    promoted_paths: List[str] = []
    promoted_stems_lower: List[str] = []
    for row, sc, why in rows:
        p = Path(row.path)
        if p.is_file() and p.resolve().parent == prepped_dir_resolved:
            promoted_paths.append(str(p.resolve()))
            stem0 = p.stem.split("_stage")[0].lower()
            promoted_stems_lower.append(stem0)
            if logger:
                logger.info(f"[forced-control] {pdb_id}: promoting '{stem0}' as control (reason={why}, score={sc})")
    promoted_paths = _dedupe_str(promoted_paths)
    promoted_stems_lower = _dedupe_str(promoted_stems_lower)
    if logger and not promoted_paths:
        logger.warning(
            f"[forced-control] {pdb_id}: requested {list(names)} but none were found in prepped library '{prepped_dir}'."
        )
    return promoted_paths, promoted_stems_lower

def collect_prepped_controls_for_protein(paths, control_stems: Sequence[str], logger=None) -> Tuple[List[str], set]:
    """Find control pdbqts colocated in the protein's prepped ligands directory."""
    control_stems_lower = {s.lower() for s in control_stems if s}
    prepped_control_pdbqts: List[str] = []
    prepped_dir_resolved = paths.prepped_ligands_dir.resolve()
    if paths.prepped_ligands_dir.exists():
        for p in paths.prepped_ligands_dir.glob("*.pdbqt"):
            stem0 = p.stem.split("_stage")[0].lower()
            if stem0 in control_stems_lower and p.parent.resolve() == prepped_dir_resolved:
                prepped_control_pdbqts.append(str(p.resolve()))
    prepped_control_pdbqts = sorted(set(prepped_control_pdbqts))
    if logger:
        leaks = [pp for pp in prepped_control_pdbqts if Path(pp).resolve().parent != prepped_dir_resolved]
        if leaks:
            logger.warning(f"[CONTROL-LEAK] Found {len(leaks)} controls outside {prepped_dir_resolved}")
            for pp in leaks[:10]:
                logger.warning(f"  leak -> {pp}")
    return prepped_control_pdbqts, control_stems_lower

# =============================
# Benchmark driver
# =============================

def run_benchmark_for_protein(
    cfg: Dict,
    mapping: MappingIndex,
    pdb_file: str,
    prepped_dir: Path,                  # controls root for this PDB
    exhaustiveness: int,
    num_modes: int,
    max_candidates: int,
    manual_hints: Optional[List[str]] = None,
    fda_index=None,
    fda_scope_dirs: Optional[Sequence[Path]] = None,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
) -> None:
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = base_id.replace("_cleaned", "").upper()

    # Logger + ASCII filter for mixed terminals
    logger = make_protein_logger(cfg["DOCKED_DIR"], pdb_id, cfg)
    import logging
    import sys

    def _sanitize_msg(s: str) -> str:
        return (
            s.replace("≤", "<=")
             .replace("≥", ">=")
             .replace("Å", " Angstrom")
             .replace("µ", "u")
             .replace("°", " deg")
        )

    class _AsciiFilter(logging.Filter):
        def filter(self, record):
            if isinstance(record.msg, str):
                record.msg = _sanitize_msg(record.msg)
            return True

    for h in logger.handlers:
        h.addFilter(_AsciiFilter())
    for stream_name in ("stdout", "stderr"):
        try:
            getattr(sys, stream_name).reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    variant_mode = variant if variant is not None else (
        cfg.get("APO_HOLO_MODE") or os.environ.get("APO_HOLO_MODE")
    )
    variant_candidates = expand_variants(variant_mode)
    variant = variant_candidates[0] if variant_candidates else None
    ph_token = ph_token or None

    # >>> PATHS INIT START
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # >>> PATHS INIT END

    prepped_lig_dir = paths.prepped_ligands_dir
    ligands_raw_dir = paths.ligand_output_dir
    ligands_mol2 = paths.ligands_mol2_dir

    # Prefer canonical cleaned receptor; if missing, allow legacy read-only fallback.
    cleaned_pdb_path = paths.receptor_cleaned_pdb(variant)
    if not cleaned_pdb_path.exists():
        base = Path(cfg.get("OUTPUT_DIR") or Path(cfg["OVERALL_DIR"]).parent / "processed_pdbs").resolve()
        legacy = base / f"{pdb_id}_NOLIG" / "receptor" / f"{pdb_id}_NOLIG_cleaned.pdb"
        if legacy.exists():
            cleaned_pdb_path = legacy

    input_pdb_path = Path(paths.input_pdb_path)

    # --- DEBUG visibility (clarify scopes) ---
    try:
        per_pdb_prepped = prepped_lig_dir.resolve()
        lib_root = Path(cfg.get("PREPPED_LIGANDS_DIR", str(per_pdb_prepped.parent))).resolve()

        logger.info("[DEBUG] per_pdb_prepped_dir=%s", str(per_pdb_prepped))
        logger.info("[DEBUG] library_root=%s", str(lib_root))

        # Mapping stats
        try:
            n_rows = len(mapping.rows)
        except Exception:
            n_rows = -1
        logger.info("[DEBUG] mapping_rows=%d", n_rows)

        if n_rows > 0:
            sample = mapping.rows[:12]  # first dozen rows as a quick sample
            bad_exist = 0
            bad_scope_lib = 0
            for r in sample:
                p = Path(r.path)
                if not p.exists():
                    bad_exist += 1
                else:
                    rp = p.resolve()
                    # check against the LIBRARY ROOT (the true selection scope)
                    if rp.parent != lib_root and lib_root not in rp.parents:
                        bad_scope_lib += 1
            logger.info("[DEBUG] mapping sample vs library_root: not_exist=%d, out_of_scope=%d (first 12)",
                        bad_exist, bad_scope_lib)
    except Exception as _e:
        logger.warning("[DEBUG] mapping-precheck failed: %s", _e)

    # 1) Extract controls & build nolig
    _, control_stems = extract_ligands_to_nolig(paths, logger)
    # Always prep the extracted crystal ligands with the same sanitizer/ADT flags (-A hydrogens)
    prep_marker = prepped_lig_dir / ".prepped.ok"
    prepped_lig_dir.mkdir(parents=True, exist_ok=True)



    # If marker exists and nothing newer is in ligands_raw, skip prep
    def _latest_mtime(globpat: str) -> float:
        import glob
        import os
        mt = 0.0
        for p in glob.glob(globpat):
            try:
                mt = max(mt, os.stat(p).st_mtime)
            except Exception:
                pass
        return mt

    # --- collapse helper (idempotent) ---
    import re
    def _collapse_sanitized_suffixes(root: Path, exts=(".pdb", ".mol2", ".pdbqt")) -> int:
        """
        Collapse runs of '.sanitized' (2 or more) just before the extension.
        Returns count of files renamed.
        """
        pat = re.compile(r'(?:\.sanitized){2,}(?=\.)')  # only collapse 2+ in a row
        fixed = 0
        for ext in exts:
            for f in root.glob(f"*{ext}"):
                new_name = pat.sub('.sanitized', f.name)
                if new_name != f.name:
                    dst = f.with_name(new_name)
                    try:
                        if dst.exists():
                            # keep the newer timestamp
                            if f.stat().st_mtime > dst.stat().st_mtime:
                                dst.unlink()
                                f.replace(dst)
                            else:
                                f.unlink()
                        else:
                            f.replace(dst)
                        fixed += 1
                    except Exception:
                        pass
        if fixed:
            logger.info("[tidy] collapsed .sanitized x%d under %s", fixed, root)
        return fixed

    # -------------------------------------

    # 0) Normalize RAW filenames first so our freshness check is stable
    _collapse_sanitized_suffixes(ligands_raw_dir, (".pdb", ".mol2"))

    # 1) Freshness check (consider both .pdb and .mol2 as raw inputs)
    raw_latest = max(
        _latest_mtime(str(ligands_raw_dir / "*.pdb")),
        _latest_mtime(str(ligands_raw_dir / "*.mol2")),
    )
    prepped_latest = _latest_mtime(str(prepped_lig_dir / "*.pdbqt"))

    if prep_marker.exists() and prepped_latest >= raw_latest:
        logger.info("[prep] prepped ligands present & up-to-date; skipping ligand prep.")
    else:
        logger.info("[prep] running ligand prep (raw newer than prepped or first run).")
        prep_ligands_from_pdb(
            ligand_output_dir=ligands_raw_dir,
            ligands_mol2_dir=ligands_mol2,
            prepped_ligands_dir=prepped_lig_dir,
        )
        try:
            prep_marker.touch()
        except Exception:
            pass

    # 2) Always normalize both trees after the decision (whether we ran prep or not)
    _collapse_sanitized_suffixes(ligands_raw_dir, (".pdb", ".mol2"))
    _collapse_sanitized_suffixes(prepped_lig_dir, (".pdbqt",))

    # 1b) metabolites -> parents
    extra_parent_ids: List[str] = []
    try:
        extra_parent_ids = ensure_parent_drugs_for_controls(
            pdb_code=pdb_id,
            ligands_raw_dir=str(ligands_raw_dir),
            fda_index=fda_index,
            max_additions=3,
        ) or []
        if extra_parent_ids:
            logger.info(f"[metabolite→parent] Will also consider parents: {', '.join(extra_parent_ids)}")
    except Exception as e:
        logger.warning(f"[metabolite→parent] resolver failed: {e}")

    # 2) Prepare receptor
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("[benchmark] receptor prep failed; skipping protein")
        return

    # 3) Detect pocket (main logic)
    center, detected_box, src = detect_pocket(cleaned_pdb, ligands_raw_dir, logger)
    # honor BOX_SIZE_MAX_A from config (default 28.0)
    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
    if center:
        pockets = [("pocket1", center)]
        box_size = tuple(min(box_cap, float(s)) for s in (detected_box or (24.0, 24.0, 24.0)))
    else:
        logger.error("[benchmark] no pocket could be detected; skipping protein")
        return
    if not pockets:
        c, _, _ = detect_active_site(cleaned_pdb)
        if c:
            pockets = [("pocket1", c)]
        else:
            logger.error("[benchmark] no pocket could be detected; skipping protein")
            return

    # 4) Build hints
    het_ids, het_names = parse_pdb_het_hints(input_pdb_path)
    hints: List[str] = []
    hints.extend(het_names)
    hints.extend(het_ids)
    for het in het_ids:
        hints.extend(CHEMCOMP_ALIAS.get(het.upper(), []))
    if manual_hints:
        hints.extend(manual_hints)
    hints.extend(extra_parent_ids)

    # dedupe normalized
    seen = set()
    final_hints: List[str] = []
    for h in hints:
        hn = _norm(h)
        if hn and hn not in seen:
            final_hints.append(h)
            seen.add(hn)
    final_hints = _merge_per_pdb_hints(pdb_id, final_hints)
    allowed_rdk_ids: set[str] = set()

    # STRICT policy (default ON): only extracted controls + hard-coded per-PDB RDK identities
    if bool(cfg.get("BENCH_ONLY_EXTRACTED_CTRL_AND_RDK", True)):
        final_hints = list(PER_PDB_HINTS.get(pdb_id, []))  # keep scope narrow
        try:
            allowed_rdk_ids = _resolve_allowed_rdk_ids(mapping, final_hints)
        except Exception:
            allowed_rdk_ids = set()
        logger.info("[policy] %s restrict_hints → PER_PDB_HINTS=%s | resolved_rdk=%s",
                    pdb_id, final_hints, sorted(list(allowed_rdk_ids))[:6])
    # --- Controls-only policy vars (used to guard candidate selection) ---
    is_hardcoded = pdb_id in HARD_FDA_CONTROL_BY_PDB
    enforce_controls_only = bool(cfg.get("BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY", True))
    allow_fallback = bool(cfg.get("BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING", False))
    hardcoded_names = HARD_FDA_CONTROL_BY_PDB.get(pdb_id, []) if is_hardcoded else []
    prepped_dir = Path(paths.prepped_ligands_dir)
    ctrl_paths_from_disk = _find_control_files(prepped_dir, hardcoded_names) if is_hardcoded else []

    logger.info("[DEBUG] hint_count=%d | examples=%s", len(final_hints), ", ".join(map(str, final_hints[:8])))
    # 5) Candidate selection (now fenced to multi-scope)
    scope_dirs = [Path(d).resolve() for d in (fda_scope_dirs or [])]
    if not scope_dirs:
        # Fallback to default fda_library under PREPPED_ROOT if not provided
        default = Path(cfg.get("PREPPED_LIGANDS_DIR", str(prepped_dir))).resolve() / "fda_library"
        scope_dirs = [default]

    cand_rows = select_candidates_for_protein(
        mapping=mapping,
        hints=final_hints,
        fda_scope_dirs=scope_dirs,
        extra_parent_ids=extra_parent_ids,
        max_candidates=max_candidates,
    )
    logger.info("[DEBUG] select_candidates scope=%s | kept=%d",
                [str(s) for s in scope_dirs], len(cand_rows))

    # --- render enqueue helpers (inside run_benchmark_for_protein) ---
    def _first(globpat: str) -> Optional[Path]:
        hits = sorted(glob.glob(globpat))
        return Path(hits[0]) if hits else None

    def _pick_best_control_and_rdk(out_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Find representative control and rdk pose .pdbqt files produced for this protein.
        We try a few common patterns; adjust if your filenames differ.
        """
        # control: crystal/control poses
        ctrl = (
                _first(str(out_dir / "*control*best*.pdbqt")) or
                _first(str(out_dir / "*CONTROL*best*.pdbqt")) or
                _first(str(out_dir / "*control*.pdbqt")) or
		_first(str(out_dir / "*best*.pdb"))                 
	)
        # rdk: library match (nearest / selected)
        rdk = (
                _first(str(out_dir / "*rdk*best*.pdbqt")) or
                _first(str(out_dir / "*RDK*best*.pdbqt")) or
                _first(str(out_dir / "*rdk*.pdbqt"))
        )
        return ctrl, rdk

    STAGE_NAME = str(cfg.get("BENCH_STAGE_NAME", "bench_pocket1_single"))

    def _pick_best_control_and_rdk_in_stage(
            paths,
            variant: Optional[str],
            stage_name: str | None = None
    ) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Controls (your examples):
          P30_B1001.sanitized_bench_pocket1_single.best.pdb
          P30_B1001_bench_pocket1_single.best.pdb
          VGH_A2346.sanitized_bench_pocket1_single.best.pdb
          1N1_B502_bench_pocket1_single.best.pdb
          1N1_A501.sanitized_bench_pocket1_single.best.pdb

        RDK (your example):
          rdk_0002967_bench_pocket1_single.pdbqt
        """
        stage_name = stage_name or STAGE_NAME
        logger.info("[render-pick] using stage_name=%s", stage_name)

        stage_dir = paths.docked_variant_root(variant) / stage_name
        if not stage_dir.is_dir():
            return None, None

        def _first(globpat: str) -> Optional[Path]:
            hits = sorted(glob.glob(globpat))
            return Path(hits[0]) if hits else None

        # control: crystal/control poses
        ctrl_patterns = [
            f"*{stage_name}.best.pdb",
            "*control*best*.pdbqt",
            "*CONTROL*best*.pdbqt",
            "*control*.pdbqt",
        ]
        ctrl = None
        for pat in ctrl_patterns:
            ctrl = _first(str(stage_dir / pat))
            if ctrl:
                break

        # rdk: library match (nearest / selected)
        rdk_patterns = [
            "*rdk*best*.pdbqt",
            "*RDK*best*.pdbqt",
            "*rdk*.pdbqt",
            f"rdk*{stage_name}*.pdbqt",
        ]
        rdk = None
        for pat in rdk_patterns:
            rdk = _first(str(stage_dir / pat))
            if rdk:
                break

        return ctrl, rdk



    # Audit
    out_dir = paths.docked_variant_root(variant)
    # Pick representative poses (best control + best rdk) from this protein's outputs
    ctrl_pose, rdk_pose = _pick_best_control_and_rdk_in_stage(paths, variant, STAGE_NAME)

    logger.info("[render-enqueue] %s %s: receptor=%s | orig=%s | ctrl=%s | rdk=%s",
                pdb_id, STAGE_NAME, cleaned_pdb_path, input_pdb_path,
                ctrl_pose, rdk_pose)
    if not _render_only_current_run_enabled(cfg):
        _enqueue_render_task(
            pdb_id=pdb_id,
            stage_name=STAGE_NAME,
            paths=paths,
            variant=variant,
            cleaned_pdb_path=str(cleaned_pdb_path.resolve()),
            original_pdb_path=str(input_pdb_path.resolve()),
            ctrl_pose_path=(str(ctrl_pose) if ctrl_pose else None),
            rdk_pose_path=(str(rdk_pose) if rdk_pose else None),
            exclude_resns=list(EXCLUDE_HET_IDS),
        )
        logger.info("[render-enqueue-final] %s %s: ctrl=%s | rdk=%s (pre-docking allowed)",
                    pdb_id, STAGE_NAME, bool(ctrl_pose), bool(rdk_pose))
    else:
        logger.info("[render-pre-skip] %s %s: BENCH_RENDER_ONLY_CURRENT_RUN=1 → "
                    "skip pre-docking enqueue; renders will come only from current-run poses.",
                    pdb_id, STAGE_NAME)

    # Enqueue consolidated renders (native + control + rdk + overlay)
    fh, cand_csv_path = _safe_open_csv_for_write(out_dir / f"benchmark_candidates_{pdb_id}.csv")
    with fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "rank",
                "score",
                "reason",
                "path",
                "display_name",
                "generic_name",
                "brand_names",
                "inchikey",
                "corresponding_name",
            ]
        )
        for i, (row, sc, why) in enumerate(cand_rows, start=1):
            corr = _resolve_name_for_path_or_text(row.path, fda_index) or row.generic_name or row.display_name
            w.writerow([i, sc, why, row.path, row.display_name, row.generic_name, row.brand_names, row.inchikey, corr])
    if str(cand_csv_path) != str(out_dir / f"benchmark_candidates_{pdb_id}.csv"):
        logger.info(f"[benchmark] candidates CSV was locked; wrote to {cand_csv_path.name}")
    # --- Controls-only policy gate for HARD_FDA_CONTROL_BY_PDB ---
    is_hardcoded = pdb_id in HARD_FDA_CONTROL_BY_PDB
    enforce_controls_only = bool(cfg.get("BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY", True))
    allow_fallback = bool(cfg.get("BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING", False))

    hardcoded_names = HARD_FDA_CONTROL_BY_PDB.get(pdb_id, []) if is_hardcoded else []
    prepped_dir = Path(paths.prepped_ligands_dir)
    ctrl_paths_from_disk = _find_control_files(prepped_dir, hardcoded_names) if is_hardcoded else []

    # Start with "controls" as discovered by today's flow (e.g., promotion/lock), then overlay if hardcoded hit:
    controls = list(controls) if "controls" in locals() else []  # keep existing discovery if already computed
    if is_hardcoded and enforce_controls_only:
        if ctrl_paths_from_disk:
            # Enforce: controls only; nuke whitelist later by feeding empty non_ctrls.
            controls = ctrl_paths_from_disk
            whitelist = []  # will be ignored; we’ll still let downstream two-wave code run with 0 whitelist
            logger.info(f"[policy] {pdb_id} policy=controls-only n_ctrl={len(controls)} n_whitelist=0 reason=hardcoded")
        else:
            looked_for = ", ".join(hardcoded_names) if hardcoded_names else "(none)"
            logger.warning(
                f"[hardcoded-control-missing-noabort] {pdb_id}: looked_for=[{looked_for}]; "
                f"continuing with discovered controls and keeping RDK whitelist.")
    else:
        # Non-hardcoded: keep current behavior. We'll print the audit once counts are known.
        pass
    # --- end controls-only policy gate ---

    whitelist = [row.path for (row, _, _) in cand_rows]

    if (not whitelist) and not (is_hardcoded and enforce_controls_only and ctrl_paths_from_disk):
        logger.warning("[benchmark] no candidate ligands matched mapping; skipping protein")
        return


    # Heavy atom counts
    heavy_atom_counts: Dict[str, int] = {}
    for p in whitelist:
        try:
            heavy_atom_counts[p] = int(_count_heavy_atoms_from_pdbqt(Path(p)))
        except Exception:
            heavy_atom_counts[p] = 0

    # Control lookup
    control_lookup = build_control_lookup(paths)

    # Controls present in THIS protein's prepped dir
    prepped_control_pdbqts, control_stems_lower = collect_prepped_controls_for_protein(
        paths=paths, control_stems=control_stems, logger=logger
    )
    if pdb_id.upper() in HARD_FDA_CONTROL_BY_PDB:
        forced_names = HARD_FDA_CONTROL_BY_PDB[pdb_id.upper()]
        forced_paths, forced_stems_lower = _promote_forced_controls_to_ctrls(
            pdb_id=pdb_id,
            mapping=mapping,
            prepped_dir=prepped_dir,
            names=forced_names,
            logger=logger,
        )
        if forced_paths:
            prepped_control_pdbqts = _dedupe_str(prepped_control_pdbqts + forced_paths)
            control_stems_lower = set(control_stems_lower) | set(forced_stems_lower)

    for p in prepped_control_pdbqts:
        try:
            heavy_atom_counts[p] = int(_count_heavy_atoms_from_pdbqt(Path(p)))
        except Exception:
            heavy_atom_counts[p] = heavy_atom_counts.get(p, 0)

    # -------- single-stage per pocket --------
    pocket_strength_rows: List[List[str]] = []

    for pocket_name, center in pockets[:2]:
        stage = {
            "name": f"bench_{pocket_name}_single",
            "exhaustiveness": int(exhaustiveness),
            "num_modes": int(num_modes),
            "energy_range": int(cfg.get("ENERGY_RANGE", 9)),
        }

        params = get_recenter_params(cfg)
        guard = GlobalCenterGuard(max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2)))
        selector = CenterSelector(
            cfg,
            logger,
            control_stems=set(control_stems),
            heavy_atom_counts=heavy_atom_counts,
            initial_center=center,
        )
        # use BOX_SIZE_MAX_A from config (default 28.0) for all clamps in this loop
        box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
        box_size: Tuple[float, float, float] = tuple(min(box_cap, float(s)) for s in (24.0, 24.0, 24.0))
        if 'detected_box' in locals():
            box_size = tuple(min(box_cap, float(s)) for s in (detected_box or box_size))

        # Ligand pool
        ctrls, non_ctrls = split_controls_and_whitelist(
            whitelist_paths=whitelist,
            prepped_control_pdbqts=prepped_control_pdbqts,
            control_stems_lower=control_stems_lower,
            heavy_atom_counts=heavy_atom_counts,
            cfg=cfg,
            logger=logger,
        )
        # When a per-PDB hardcoded list exists, filter RDK to those names (by resolved brand/generic).
        if is_hardcoded and hardcoded_names:
            allowed = {_norm(n) for n in hardcoded_names}

            def _rdk_name_ok(p: str) -> bool:
                # Prefer exact RDK id fence when available
                rid = _extract_rdk_id(Path(p).name)
                if allowed_rdk_ids and rid:
                    return rid.lower() in allowed_rdk_ids
                nm = _resolve_name_for_path_or_text(p, fda_index)
                return _norm(nm) in allowed

            filtered = [p for p in non_ctrls if _rdk_name_ok(p)]
            before = len(non_ctrls)
            if filtered:
                non_ctrls = filtered
                logger.info("[rdk-whitelist] %s: filtered RDK from %d→%d using %s",
                            pdb_id, before, len(non_ctrls), hardcoded_names)
            else:
                if bool(cfg.get("BENCH_ONLY_EXTRACTED_CTRL_AND_RDK", False)):
                    non_ctrls = []
                    logger.warning("[rdk-whitelist-empty-enforced] %s: no matches for %s; non_ctrls cleared by policy.",
                                   pdb_id, hardcoded_names)
                else:
                    logger.warning("[rdk-whitelist-empty] %s: no matches for %s; keeping original list.",
                                   pdb_id, hardcoded_names)

        # toggles to adjust
        max_ctrls = int(cfg.get("BENCH_MAX_CONTROLS", 8))

        max_total = max_ctrls+3
        if max_total > 0:
            allowed_non_ctrls = max(0, max_total - len(ctrls))
            if allowed_non_ctrls < len(non_ctrls):
                if logger:
                    logger.info(
                        f"[CAP] {pdb_id}: limiting non-controls from {len(non_ctrls)} to {allowed_non_ctrls} "
                        f"(controls={len(ctrls)}, cap={max_total})"
                    )
                non_ctrls = non_ctrls[:allowed_non_ctrls]

        ligands = _dedupe_str(ctrls + non_ctrls)
        # If strict allow-list yielded 0 RDKs, mark sentinel so analysis can skip cleanly
        if not non_ctrls:
            try:
                stage_dir_target = paths.docked_stage_dir(variant, stage["name"])
                (stage_dir_target / ".rdk_skipped").touch()
                logger.info(f"[rdk-sentinel] {pdb_id} {stage['name']}: no RDK candidates -> wrote .rdk_skipped")
            except Exception as _e:
                logger.warning(f"[rdk-sentinel] {pdb_id} {stage['name']}: could not write .rdk_skipped: {_e}")
        stage1_original = ligands[:]
        logger.info(f"[DEBUG] pools: ctrls={len(ctrls)}, whitelist={len(non_ctrls)}, total={len(ligands)}")
        if not (is_hardcoded and enforce_controls_only and ctrl_paths_from_disk):
            # Default path (or fallback enabled): print the same audit but marked as default
            try:
                logger.info(
                    f"[policy] {pdb_id} policy=controls+whitelist n_ctrl={len(controls)} "
                    f"n_whitelist={len(whitelist)} reason=default"
                )
            except Exception:
                pass

        # Per-ligand budget guards (persist across restarts/waves)
        per_ligand_seconds = float(cfg.get("BENCH_MAX_SECONDS", 9000.0))
        guards_all: Dict[str, BudgetGuard] = {lig: BudgetGuard(per_ligand_seconds) for lig in ligands}

        score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
        validated_ligands_last: List[str] = []
        recenter_attempts = 0

        # Optional checkpoint skip
        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
            fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
            if checkpoint_should_skip(cfg, pdb_id, stage["name"], fp):
                logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                continue

        restarting = True
        while restarting:
            restarting = False
            guard.reset_stage()

            if not ligands:
                logger.warning(f"No ligands to dock at {stage['name']}; skipping pocket.")
                break

            logger.info(f"[benchmark] {pdb_id} | {pocket_name} | starting {stage['name']} with {len(ligands)} ligands")

            if ctrls and non_ctrls:
                logger.info(f"Single-stage two-wave: {len(ctrls)} controls first, then {len(non_ctrls)} whitelist.")

                # Wave A - controls
                with _patch_main_threadpool(cfg.get("GLOBAL_LIGAND_SEM")):
                    s1, v1, d1, rd1, inv1 = run_one_stage(
                        cfg=cfg,
                        pdb_id=pdb_id,
                        receptor_pdbqt=receptor_pdbqt,
                        center=center,
                        box_size=box_size,
                        stage=stage,
                        ligands=ctrls,
                        logger=logger,
                        retry_mgr=RetryManager(),
                        control_lookup=control_lookup,
                        budget_guards={lig: guards_all[lig] for lig in ctrls},
                    )

                # CenterSelector after controls
                try:
                    dec = selector.consider_switch(stage["name"], s1, v1, rd1, receptor_pdbqt, center, guard)
                    if dec.promoted and dec.new_center is not None:
                        old = center
                        center = dec.new_center
                        guard.mark_switch()
                        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                            checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                        logger.info(f"[CENTER] switched (controls wave): {old} -> {center} ({dec.reason}) [global switch]")
                except Exception as e:
                    logger.warning(f"CenterSelector (controls) failed: {e}")

                # Optional control lock
                lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", -6.0))
                lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))
                qualified_controls = []
                for lig in v1:
                    stem = Path(lig).stem.split("_stage")[0].lower()
                    if stem not in control_stems_lower:
                        continue
                    sc = s1.get(lig)
                    if sc is None or not isinstance(sc, (int, float)):
                        continue
                    if sc > lock_score_max:
                        continue
                    pose_path = rd1.get(lig)
                    if not pose_path:
                        continue
                    c = CenterSelector._pdbqt_centroid(pose_path)
                    if c is None:
                        continue
                    import numpy as np
                    if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                        qualified_controls.append(lig)
                if len(qualified_controls) >= lock_min_hits and not guard.locked:
                    guard.lock()
                    logger.info(
                        f"[CONTROL-LOCK] n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Angstrom"
                    )

                # Wave B - whitelist
                with _patch_main_threadpool(cfg.get("GLOBAL_LIGAND_SEM")):
                    s2, v2, d2, rd2, inv2 = run_one_stage(
                        cfg=cfg,
                        pdb_id=pdb_id,
                        receptor_pdbqt=receptor_pdbqt,
                        center=center,
                        box_size=box_size,
                        stage=stage,
                        ligands=non_ctrls,
                        logger=logger,
                        retry_mgr=RetryManager(),
                        control_lookup=control_lookup,
                        budget_guards={lig: guards_all[lig] for lig in non_ctrls},
                    )

                # Merge
                scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
                raw_docked = {**rd1, **rd2}
                invalids = {**inv1, **inv2}
            else:
                # Single wave
                with _patch_main_threadpool(cfg.get("GLOBAL_LIGAND_SEM")):
                    scores, validated, distances, raw_docked, invalids = run_one_stage(
                        cfg=cfg,
                        pdb_id=pdb_id,
                        receptor_pdbqt=receptor_pdbqt,
                        center=center,
                        box_size=box_size,
                        stage=stage,
                        ligands=ligands,
                        logger=logger,
                        retry_mgr=RetryManager(),
                        control_lookup=control_lookup,
                        budget_guards={lig: guards_all[lig] for lig in ligands},
                    )

            validated_ligands_last = validated

            # Record
            for lig, sc in scores.items():
                record_score(score_history, stage["name"], lig, sc, True)
                record_le(score_history, stage["name"], lig, sc, heavy_atom_counts)
            for lig, (sc, reason) in invalids.items():
                record_score(score_history, stage["name"], lig, sc, False, reason=str(reason))
                record_le(score_history, stage["name"], lig, sc, heavy_atom_counts)

            # CenterSelector after full run
            try:
                decision = selector.consider_switch(
                    stage["name"], scores, validated, raw_docked, receptor_pdbqt, center, guard
                )
                if decision.promoted and decision.new_center is not None:
                    old = center
                    center = decision.new_center
                    guard.mark_switch()
                    if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                        checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                    logger.info(f"[CENTER] {pdb_id} {pocket_name}: {old} -> {center} ({decision.reason}) [global switch]")
            except Exception as e:
                logger.warning(f"CenterSelector failed gracefully: {e}")

            # Early recenter (single-stage flavor)
            restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                0,
                scores,
                distances,
                box_size,
                center,
                stage1_original,
                recenter_attempts,
                params,
                cfg,
                pdb_id,
                receptor_pdbqt,
                logger,
                raw_docked,
                guard,
                control_anchor_hit=False,
            )
            if restart:
                ligands = redo_ligands
                if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                    checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                restarting = True
                continue

            # Adaptive shrink
            try:
                if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                    import numpy as np
                    med = float(
                        np.median([d for d in distances if isinstance(d, (int, float))])
                    ) if distances else None
                    if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                        dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                        min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                        new_box = tuple(max(min_box, s - dec) for s in box_size)
                        if new_box != box_size:
                            logger.info(f"Adaptive shrink: median {med:.2f} Angstrom -> box {box_size} -> {new_box}")
                            box_size = new_box
            except Exception as _e:
                logger.warning(f"Adaptive shrink skipped: {_e}")

            # Checkpoint done
            if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                try:
                    fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                    checkpoint_mark_done(cfg, pdb_id, stage["name"], fp)
                except Exception:
                    pass

            # ---- outputs ----
            out_dir.mkdir(parents=True, exist_ok=True)
            long_csv = out_dir / f"docking_score_long_{pocket_name}.csv"
            results = score_history.get(stage["name"], {})
            ordered = sorted(
                results.items(),
                key=lambda kv: (
                    kv[1].get("score", 1e9) if isinstance(kv[1].get("score"), (int, float)) else 1e9,
                    os.path.basename(kv[0]),
                ),
            )
            with open(long_csv, "w", newline="", encoding="utf-8") as f:
                f.write(f"# pocket={pocket_name}, center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f})\n")
                w = csv.writer(f)
                w.writerow(["stage", "ligand", "library_id", "corresponding_name", "score", "valid", "reason"])
                for lig, rec in ordered:
                    sc = rec.get("score")
                    sc_str = f"{sc:.2f}" if isinstance(sc, (int, float)) else ""
                    lib_id = _extract_rdk_id(lig) or ""
                    corr = _resolve_name_for_path_or_text(lig, fda_index)
                    w.writerow(
                        [
                            stage["name"],
                            os.path.basename(lig),
                            lib_id,
                            corr,
                            sc_str,
                            int(bool(rec.get("valid", False))),
                            str(rec.get("reason", "")),
                        ]
                    )

            if validated_ligands_last:
                final_pose_validation_and_screenshots(
                    cfg,
                    pdb_id,
                    [stage],
                    receptor_pdbqt,
                    center,
                    validated_ligands_last,
                    score_history,
                    cleaned_pdb,
                    cfg.get("DOCKING_MODE", "benchmark"),
                    logger,
                )

            # Optional PyMOL renders (when DOCKING_MODE == "benchmark")
            if str(cfg.get("DOCKING_MODE", "")).lower() == "benchmark":
                stage_name = stage["name"]

                results_for_stage = score_history.get(stage_name, {})
                best_ctrl_lig, nearest_rdk_lig = pick_control_and_nearest_rdk(
                    results_for_stage, raw_docked, control_stems_lower
                )

                cleaned_pdb_str = str(cleaned_pdb)
                original_pdb_str = str(input_pdb_path)
                ctrl_pose_path = raw_docked.get(best_ctrl_lig) if best_ctrl_lig else None
                rdk_pose_path = raw_docked.get(nearest_rdk_lig) if nearest_rdk_lig else None

                if bool(cfg.get("DEFER_PYMOL", True)):
                    if not ctrl_pose_path and not rdk_pose_path:
                        logger.info(
                            f"[render-skip] {pdb_id} {stage_name}: no ctrl/rdk pose paths → no screenshots will be generated.")
                    else:
                        logger.info(f"[render-enqueue] {pdb_id} {stage_name}: "
                                    f"ctrl={'yes' if ctrl_pose_path else 'no'} | rdk={'yes' if rdk_pose_path else 'no'}")
                    _enqueue_render_task(
                        pdb_id=pdb_id,
                        stage_name=stage_name,
                        paths=paths,
                        variant=variant,
                        cleaned_pdb_path=cleaned_pdb_str,
                        original_pdb_path=original_pdb_str,
                        ctrl_pose_path=ctrl_pose_path,
                        rdk_pose_path=rdk_pose_path,
                        exclude_resns=sorted(list(EXCLUDE_HET_IDS)),
                    )
                else:
                    with _RENDER_LOCK:
                        stage_dir_target = paths.docked_stage_dir(variant, stage_name)

                        _render_native_on_original_pdb(
                            original_pdb=original_pdb_str,
                            outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-native__NATIVE"),
                            exclude_resns=sorted(list(EXCLUDE_HET_IDS)),
                        )
                        if ctrl_pose_path and Path(ctrl_pose_path).is_file():
                            _render_three_views_with_pymol(
                                receptor_path=cleaned_pdb_str,
                                ligand_paths_and_colors=[(ctrl_pose_path, "control", "blue")],
                                outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL"),
                            )
                        if rdk_pose_path and Path(rdk_pose_path).is_file():
                            _render_three_views_with_pymol(
                                receptor_path=cleaned_pdb_str,
                                ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "orange")],
                                outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__RDKclosest"),
                            )
                            _render_three_views_with_pymol(
                                receptor_path=original_pdb_str,
                                ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "orange")],
                                outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-rdk__RDKclosest"),
                            )
                            if ctrl_pose_path and Path(ctrl_pose_path).is_file():
                                _render_three_views_with_pymol(
                                    receptor_path=cleaned_pdb_str,
                                    ligand_paths_and_colors=[(ctrl_pose_path, "control", "blue"),
                                                             (rdk_pose_path, "rdk_closest", "orange")],
                                    outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL+RDKclosest"),
                                    label_top_n_res=label_top_n_res,
                                    label_cutoff=label_cutoff,
                                )

            scores_only = [rec.get("score") for rec in results.values() if isinstance(rec.get("score"), (int, float))]
            best = min(scores_only) if scores_only else None
            n_valid = sum(1 for rec in results.values() if rec.get("valid"))
            n_tested = len(results)
            pocket_strength_rows.append(
                [
                    pocket_name,
                    f"{center[0]:.3f}",
                    f"{center[1]:.3f}",
                    f"{center[2]:.3f}",
                    (f"{best:.2f}" if isinstance(best, (int, float)) else ""),
                    n_valid,
                    n_tested,
                ]
            )

    if pocket_strength_rows:
        summary_csv = out_dir / "benchmark_pocket_strength.csv"
        with open(summary_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["pocket", "center_x", "center_y", "center_z", "best_score", "n_valid", "n_tested"])
            for row in pocket_strength_rows:
                w.writerow(row)

# =============================
# Analysis wrapper
# =============================

def _run_auto_analysis(
    docked_root: Path,
    mapping_csv: Path,
    only_pdb: Optional[str],
    score_tol: float,
    center_tol: float,
    rmsd_tol: float,
    include_identity: bool,
    out_dir: Optional[Path] = None,
    since_epoch: Optional[float] = None,
    html_success_only: Optional[bool] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    try:
        from tools import benchmark_auto_analysis as ana
    except Exception as e:
        print(f"[analysis] ERROR: could not import benchmark_auto_analysis.py: {e}")
        return None, None
    if not hasattr(ana, "run_analysis"):
        print("[analysis] ERROR: benchmark_auto_analysis.py is missing run_analysis().")
        return None, None
    try:
        return ana.run_analysis(
            docked_root=docked_root,
            mapping_csv=mapping_csv,
            only_pdb=only_pdb,
            score_tol=score_tol,
            center_tol=center_tol,
            rmsd_tol=rmsd_tol,
            include_identity=include_identity,
            out_dir=out_dir,
            since_epoch=since_epoch,
            html_success_only=html_success_only,
        )
    except Exception as e:
        print(f"[analysis] ERROR running run_analysis(): {e}")
        return None, None

# =============================
# CLI
# =============================

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Benchmark mode: single ultra-stage docking of likely co-crystal FDA ligands, per pocket."
    )
    # Defaults are None; we resolve from YAML or fallbacks in main()
    p.add_argument("--input-dir", default=None)
    p.add_argument("--out-root", default=None)
    p.add_argument("--prepped", default=None)
    p.add_argument("--mapping", default=None)
    p.add_argument("--max-candidates", type=int, default=9)
    p.add_argument("--exhaustiveness", type=int, default=1)
    p.add_argument("--num-modes", type=int, default=1)
    p.add_argument("--hints", help="Optional manual comma-separated hints (e.g., 'imatinib,STI571')")
    # Controls-only policy toggles
    p.add_argument("--enforce-hardcoded-controls-only",
                   dest="enforce_controls_only", action="store_true",
                   help="If PDB has a hardcoded control, dock only that control.")
    p.add_argument("--no-enforce-hardcoded-controls-only",
                   dest="enforce_controls_only", action="store_false",
                   help="Disable controls-only enforcement (use whitelist).")
    p.set_defaults(enforce_controls_only=None)  # let YAML/ENV decide unless provided

    p.add_argument("--allow-whitelist-fallback-if-controls-missing",
                   dest="allow_whitelist_fallback", action="store_true",
                   help="If hardcoded controls are missing on disk, allow RDK whitelist fallback.")

    # Outer (proteins) parallelism
    p.add_argument(
        "--jobs",
        type=int,
        default=max(1, (os.cpu_count() or 4)),
        help="Number of proteins to process in parallel",
    )

    #  CPU/parallel shaping flags
    p.add_argument(
        "--total-cpus",
        type=int,
        default=None,
        help="Override: total logical CPUs available to this run (auto-detected if omitted).",
    )
    p.add_argument(
        "--ligand-workers",
        type=int,
        default=None,
        help="Max ligands per protein to dock concurrently. If omitted, auto = floor(total_cpus / jobs).",
    )
    p.add_argument(
        "--blas-threads",
        type=int,
        default=1,
        help="Set OMP/BLAS thread env (OMP/MKL/OPENBLAS/etc) to avoid oversubscription (default: 1).",
    )

    p.add_argument(
        "--alias-from-details",
        nargs="*",
        default=[],
        help=("Path(s) or glob(s) to benchmark_analysis_details.csv to harvest aliases"),
    )

    # Budget CLI flag
    p.add_argument(
        "--max-seconds",
        type=float,
        default=9000.0,
        help="Per-ligand wall-clock budget in seconds for retries/validation (default: 9000.0).",
    )

    # Auto-Analysis options
    p.add_argument(
        "--skip-analysis",
        dest="run_analysis",
        action="store_false",
        help="Skip the post-run Auto-Analysis.",
    )
    p.set_defaults(run_analysis=True)
    p.add_argument("--analysis-only-pdb", default=None, help="Optional: restrict analysis to one PDB (e.g., 5MO4).")
    p.add_argument("--analysis-score-tol", type=float, default=1.0, help="Score tolerance (kcal/mol) for a point.")
    p.add_argument("--analysis-center-tol", type=float, default=1.0, help="Centroid distance tolerance (Å).")
    p.add_argument("--analysis-rmsd-tol", type=float, default=3.0, help="RMSD tolerance (Å).")
    p.add_argument("--analysis-no-identity", action="store_true", help="Exclude identity-match from scoring.")
    p.add_argument("--analysis-out", default=None, help="Optional output dir for analysis CSVs.")
    p.add_argument("--only-pdb", default=None, help="Comma-separated PDB IDs to dock, e.g. '6GQO,4XUF'.")
    p.add_argument(
        "--ligands-folder",
        dest="ligands_folders",
        action="append",
        default=None,
        help=("Limit FDA/RDK search to this subfolder of PREPPED_ROOT (repeatable). "
              "Default: fda_library. Examples: "
              "--ligands-folder fda_library2  "
              "--ligands-folder fda_library --ligands-folder fda_new"),
    )
    p.add_argument("--apo-holo-mode",
                   choices=["holo", "apo", "apo_vs_holo", "apovsholo"],
                   default="apo_vs_holo",
                   help="Dock holo, apo, or both (default: apo_vs_holo).")



    return p









# =============================
# CPU detection helpers
# =============================

def _env_int(*names: str) -> Optional[int]:
    """Return the first present env var (int) among names, else None."""
    for n in names:
        v = os.environ.get(n)
        if v:
            try:
                return int(v)
            except Exception:
                pass
    return None

def _detect_available_cpus() -> int:
    """
    Conservative CPU detector for clusters/desktops.
    Priority:
      1) SLURM_CPUS_PER_TASK
      2) SLURM_CPUS_ON_NODE
      3) NSLOTS / PBS_NP
      4) OMP_NUM_THREADS
      5) os.cpu_count()
    """
    for v in (
        _env_int("SLURM_CPUS_PER_TASK"),
        _env_int("SLURM_CPUS_ON_NODE"),
        _env_int("NSLOTS", "PBS_NP"),
        _env_int("OMP_NUM_THREADS"),
        os.cpu_count() or 1,
    ):
        if v and v > 0:
            return int(v)
    return 1

def _set_thread_env(n: int) -> None:
    """Pin common threaded libs to n to avoid hidden oversubscription."""
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_MAX_THREADS",
        "BLIS_NUM_THREADS",
    ):
        os.environ.setdefault(key, str(n))


# -----------------------------
# Small CLI helper shims (match main.py behavior)
# -----------------------------
def _cli_has(argv, flag: str) -> bool:
    try:
        return any((x or "").strip() == flag for x in argv)
    except Exception:
        return False

def _cli_val(argv, flag: str):
    try:
        argv = list(argv or [])
        for i, tok in enumerate(argv):
            if (tok or "").strip() == flag and i + 1 < len(argv):
                return argv[i + 1]
    except Exception:
        pass
    return None


# =============================
# Main (Linux/BRCF friendly)
# =============================

def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_argparser().parse_args(argv)

    # Load YAML config first so we can use BRCF paths by default
    cfg = load_inputs()
    validate_config(cfg)
    # --- Per-run configs (RUN_DIR) ---  # (copied from main.py)
    cfg.setdefault("CONFIGS_DIR", str(output_root(Path(cfg["OVERALL_DIR"]), "configs")))
    cfg.setdefault("RESET_CONFIGS", True)

    # CLI > ENV > CFG
    cli_cfg_dir = _cli_val(sys.argv, "--configs-dir")
    cli_run_id = _cli_val(sys.argv, "--run-id")
    cli_no_reset = _cli_has(sys.argv, "--no-reset-configs")

    if cli_cfg_dir:  cfg["CONFIGS_DIR"] = cli_cfg_dir
    if cli_run_id:   cfg["RUN_ID"] = cli_run_id
    if cli_no_reset: cfg["RESET_CONFIGS"] = False

    init_config_run_dir(cfg,
                        run_id=cfg.get("RUN_ID"),
                        reset=cfg.get("RESET_CONFIGS"),
                        logger=logging.getLogger("run"),
                        )
    print(f"[cfg.run] run_id={cfg['RUN_ID']} run_dir={cfg['CONFIG_RUN_DIR']}")
    # Stamp the run start and export for analysis (HTML/CSV mtime filter)
    global RUN_START_EPOCH
    RUN_START_EPOCH = float(time.time())
    os.environ["BENCH_RUN_START_EPOCH"] = str(RUN_START_EPOCH)

    # CLI → cfg overrides for new knobs (only if explicitly passed)
    if getattr(args, "enforce_controls_only", None) is not None:
        cfg["BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY"] = bool(args.enforce_controls_only)
    if getattr(args, "allow_whitelist_fallback", False):
        cfg["BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING"] = True

    # Re-apply deferral and render parallelism using config + env override
    _defer = str(os.environ.get("DEFER_PYMOL", "") or cfg.get("DEFER_PYMOL", "1")).lower() in ("1", "true", "yes", "on")
    capture_pose.set_defer_mode(_defer, queue_path=os.environ.get("PYMOL_DEFER_QUEUE") or cfg.get("PYMOL_DEFER_QUEUE"))

    global _RENDER_LOCK
    _RENDER_LOCK = threading.Semaphore(
        int(os.environ.get("PYMOL_PARALLEL", str(cfg.get("PYMOL_PARALLEL", 1))))
    )

    # Helper: first existing path in a list (returns Path or None)
    def _first_existing(*cands) -> Optional[Path]:
        for c in cands:
            if not c:
                continue
            p = Path(str(c))
            if p.exists():
                return p
        return None

    # 1) Resolve INPUT / OUTPUT / PREPPED
    input_dir = Path(args.input_dir or cfg.get("INPUT_DIR", ""))
    out_root  = Path(args.out_root  or cfg.get("DOCKED_DIR") or cfg.get("OUTPUT_DIR", ""))
    prepped   = Path(args.prepped   or cfg.get("PREPPED_LIGANDS_DIR", ""))


    # FDA scope resolution from CLI (default: ['fda_library'])
    lig_folders = args.ligands_folders if args.ligands_folders else ["fda_library"]
    fda_scope_dirs: List[Path] = []

    _missing: List[str] = []
    for name in lig_folders:
        cand = (prepped / name).resolve()
        if cand.is_dir():
            fda_scope_dirs.append(cand)
        else:
            _missing.append(str(cand))

    if _missing:
        print(f"[benchmark] WARN: missing FDA scope folders ({len(_missing)}):")
        for m in _missing:
            print(f"  - {m}")

    if not fda_scope_dirs:
        print("[benchmark] WARN: no valid FDA scope dirs; will run controls only.")
    else:
        print("[DEBUG] FDA scope dirs:", [str(d) for d in fda_scope_dirs])
        # Optional quick inventory per scope dir
        for d in fda_scope_dirs:
            try:
                n = len(list(Path(d).glob("*.pdbqt")))
                print(f"[DEBUG] in-scope files: {d} -> {n} *.pdbqt")
            except Exception:
                pass

    # 2) Resolve mapping CSV in priority order
    mapping_path = (
        Path(args.mapping) if args.mapping else
        _first_existing(
            cfg.get("FDA_MAPPING_CSV"),
            prepped / "fda_mapping_from_pdbqt.csv",
            Path(cfg.get("OVERALL_DIR", Path(out_root).parent if out_root else Path.cwd()))
            / "chemdb"
            / "data"
            / "fda_mapping_from_pdbqt.csv",
            Path.cwd() / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv",
        )
    )

    # 3) Sanity checks with helpful messages
    if not input_dir.exists():
        print(f"[benchmark] INPUT_DIR does not exist: {input_dir}")
        return

    if not prepped.exists():
        print(f"[benchmark] PREPPED library does not exist: {prepped}")
        return

    if not out_root:
        out_root = Path(cfg.get("OUTPUT_DIR", "")) if cfg.get("OUTPUT_DIR") else (Path.cwd() / "benchmarks")
    out_root.mkdir(parents=True, exist_ok=True)

    if not mapping_path or not Path(mapping_path).exists():
        print(
            "[benchmark] Mapping CSV not found.\n"
            "  Tried (in order):\n"
            f"    --mapping={args.mapping}\n"
            f"    cfg['FDA_MAPPING_CSV']={cfg.get('FDA_MAPPING_CSV')}\n"
            f"    {prepped / 'fda_mapping_from_pdbqt.csv'}\n"
            f"    {Path(cfg.get('OVERALL_DIR', Path(out_root).parent)) / 'chemdb' / 'data' / 'fda_mapping_from_pdbqt.csv'}\n"
            f"    {Path.cwd() / 'chemdb' / 'data' / 'fda_mapping_from_pdbqt.csv'}\n"
            "  Fix by passing --mapping /path/to/fda_mapping_from_pdbqt.csv or setting FDA_MAPPING_CSV in your config.txt."
        )
        return

    # 4) Build MappingIndex now that we know the file exists
    mapping = MappingIndex(Path(mapping_path))

    # 5) Freeze resolved paths back into cfg for downstream calls
    cfg = dict(cfg)
    cfg["INPUT_DIR"] = str(input_dir)
    cfg["DOCKED_DIR"] = str(out_root)
    cfg["OUTPUT_DIR"] = str(Path(cfg.get("OUTPUT_DIR", Path(out_root).parent / "processed_pdbs")))
    cfg["PREPPED_LIGANDS_DIR"] = str(prepped)
    cfg["DOCKING_MODE"] = "benchmark"
    cfg.setdefault("OVERALL_DIR", str(Path(out_root).parent))
    cfg["BENCH_MAX_SECONDS"] = float(args.max_seconds)

    # Defer ALL PyMOL work until the very end (no mid-run rendering)
    cfg["DEFER_PYMOL"] = True
    # Optional: choose a queue file; default is OVERALL_DIR/deferred_pymol_jobs.jsonl
    # os.environ["PYMOL_DEFER_QUEUE"] = str(Path(cfg["OVERALL_DIR"]) / "deferred_pymol_jobs.jsonl")
    capture_pose.set_defer_mode(True, queue_path=os.environ.get("PYMOL_DEFER_QUEUE"))


    # Optionally harvest aliases from prior results
    alias_sources = args.alias_from_details or []
    if alias_sources:
        try:
            extra_aliases = load_aliases_from_details_csv(alias_sources)
            if extra_aliases:
                global CHEMCOMP_ALIAS
                CHEMCOMP_ALIAS = merge_aliases_into_chemcomp(CHEMCOMP_ALIAS, extra_aliases)
                # keep values normalized+deduped
                for _k, _vals in list(CHEMCOMP_ALIAS.items()):
                    CHEMCOMP_ALIAS[_k] = sorted(set(v.lower() for v in _vals))
                total_added = sum(len(v) for v in extra_aliases.values())
                print(f"[alias] Expanded CHEMCOMP_ALIAS with {total_added} aliases across {len(extra_aliases)} HET codes.")
        except Exception as e:
            print(f"[alias] WARN: could not augment aliases: {e}")

    # enumerate PDBs
    pdb_files = [f for f in os.listdir(cfg["INPUT_DIR"]) if f.lower().endswith(".pdb") and "_nolig" not in f.lower()]
    print(f"[benchmark] proteins queued: {len(pdb_files)} from {cfg['INPUT_DIR']}")

    # ONLY-PDB filter (optional)
    if args.only_pdb:
        wanted = {s.strip().upper() for s in str(args.only_pdb).split(",") if s.strip()}
        before = len(pdb_files)
        pdb_files = [f for f in pdb_files if os.path.splitext(f)[0].upper() in wanted]
        print(f"[benchmark] --only-pdb={sorted(wanted)} -> matched {len(pdb_files)} of {before} files:")
        for f in pdb_files:
            print(f"  - {f}")

    fda_index = load_library_index(Path(mapping_path))
    manual_hints = [h.strip() for h in (args.hints or "").split(",") if h.strip()] or None

    # ---- Dynamic CPU-saturating scheduler (outer proteins × inner ligands) ----
    from collections import deque

    total_cpus = int(args.total_cpus) if args.total_cpus else _detect_available_cpus()
    _set_thread_env(int(args.blas_threads))


    # <<< global token pool used by ALL ligand tasks (live rebalancing of cpu usage) >>>
    global_ligand_sem = threading.BoundedSemaphore(total_cpus)
    cfg["GLOBAL_LIGAND_SEM"] = global_ligand_sem
    # Quick, cheap prescan: estimate per-protein candidate counts to weight scheduling.
    # (We intentionally skip metabolite→parent expansion here for speed.)
    prescan: List[Tuple[str, int]] = []
    for pdb in pdb_files:
        pdb_path = Path(cfg["INPUT_DIR"]) / pdb
        try:
            het_ids, het_names = parse_pdb_het_hints(pdb_path)
        except Exception:
            het_ids, het_names = [], []
        hints = list(het_names) + list(het_ids)
        for het in het_ids:
            hints.extend(CHEMCOMP_ALIAS.get(het.upper(), []))
        if manual_hints:
            hints.extend(manual_hints)

        # prescan estimate uses fenced FDA scope
        try:
            est = len(select_candidates_for_protein(
                mapping=mapping,
                hints=hints,
                fda_scope_dirs=fda_scope_dirs,  # <<< fence to chosen FDA folders
                extra_parent_ids=[],  # keep fast
                max_candidates=int(args.max_candidates),
            ))
        except Exception:
            est = 8
        prescan.append((pdb, max(1, est)))

    # Big jobs first helps keep CPUs busy as small jobs finish.
    prescan.sort(key=lambda t: t[1], reverse=True)
    q = deque(prescan)

    # Outer cap: how many proteins can run at once
    max_protein_jobs = max(1, min(int(args.jobs) if args.jobs else total_cpus, len(pdb_files), total_cpus))

    print(f"[parallel] total_cpus={total_cpus} | max_proteins={max_protein_jobs} | queued={len(prescan)}")

    def _work(pdb_file: str, inner_for_this: int):
        cfg_local = dict(cfg)
        cfg_local["GLOBAL_LIGAND_SEM"] = cfg["GLOBAL_LIGAND_SEM"]
        # Keep per-protein workers aligned with the scheduler's assigned share.
        cfg_local["CPU"] = max(1, int(inner_for_this))

        # Determine which variants to run for this protein
        variant_mode = args.apo_holo_mode or os.environ.get("APO_HOLO_MODE") or cfg.get("APO_HOLO_MODE")

        for variant in expand_variants(variant_mode):
            if variant is None:
                os.environ.pop("APO_HOLO_MODE", None)
            else:
                os.environ["APO_HOLO_MODE"] = variant

            cfg_v = dict(cfg_local)
            cfg_v["APO_HOLO_MODE"] = variant

            run_benchmark_for_protein(
                cfg=cfg_v,
                mapping=mapping,
                pdb_file=pdb_file,
                prepped_dir=prepped,
                exhaustiveness=int(args.exhaustiveness),
                num_modes=int(args.num_modes),
                max_candidates=int(args.max_candidates),
                manual_hints=manual_hints,
                fda_index=fda_index,
                fda_scope_dirs=fda_scope_dirs,
                variant=variant,
            )

        return pdb_file

    import traceback

    # Submit loop: we keep track of "CPU slots" allocated to each running protein.
    running: Dict[object, Tuple[str, int]] = {}
    sum_inner = 0

    def _alloc_plan(remaining_slots: int, remaining_jobs: int) -> int:
        """Greedy: give each new protein at least floor, first few get +1 to absorb remainder."""
        base = max(1, remaining_slots // max(1, remaining_jobs))
        extra = remaining_slots - base * max(1, remaining_jobs)
        return base + (1 if extra > 0 else 0)

    with ThreadPoolExecutor(max_workers=max_protein_jobs) as pool:
        # Seed as many as we can to fill CPU budget
        while q and len(running) < max_protein_jobs and sum_inner < total_cpus:
            remaining_slots = total_cpus - sum_inner
            remaining_jobs = max_protein_jobs - len(running)
            inner = _alloc_plan(remaining_slots, remaining_jobs)
            pdb, est = q.popleft()
            fut = pool.submit(_work, pdb, inner)
            running[fut] = (pdb, inner)
            sum_inner += inner

        done = failed = 0
        total = len(pdb_files)

        while running:
            fut = next(as_completed(list(running.keys())))
            pdb, inner = running.pop(fut)
            sum_inner -= inner
            try:
                fut.result()
                done += 1
                print(f"[benchmark] ✔ finished {pdb} (used_inner={inner}) ({done}/{total})")
            except Exception as e:
                failed += 1
                tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
                print(f"[benchmark] ✖ error on {pdb} (used_inner={inner}): {e} ({done + failed}/{total})\n{tb}")

            # Refill immediately, redistributing freed CPU to the next proteins
            while q and len(running) < max_protein_jobs and sum_inner < total_cpus:
                remaining_slots = total_cpus - sum_inner
                remaining_jobs  = max_protein_jobs - len(running)
                inner = _alloc_plan(remaining_slots, remaining_jobs)
                next_pdb, _est = q.popleft()
                fut2 = pool.submit(_work, next_pdb, inner)
                running[fut2] = (next_pdb, inner)
                sum_inner += inner
    # ----- Post-run: render everything without stalling -----
    capture_pose.set_defer_mode(False)

    # ----- Post-run: render everything out-of-process (multi-process) -----
    def _int_env(name: str, default: int) -> int:
        try:
            v = int((os.environ.get(name, "") or "").strip() or str(default))
            return v if v > 0 else default
        except Exception:
            return default

    workers_cfg = _int_env("PYMOL_RENDER_WORKERS", int(cfg.get("PYMOL_RENDER_WORKERS", 30)))
    mode_cfg    = (os.environ.get("PYMOL_RENDER_MODE") or str(cfg.get("PYMOL_RENDER_MODE", "cli"))).lower()

    # Env-only diagnostics (no behavior change unless envs are set)
    soft_cap  = _int_env("PYMOL_RENDER_SOFT_CAP", 0)
    batch_sz  = _int_env("PYMOL_RENDER_BATCH", 0)
    mode_force = (os.environ.get("PYMOL_RENDER_MODE_FORCE") or "").strip().lower()
    ab_env     = (os.environ.get("PYMOL_RENDER_AB") or "").strip()
    queue_path = os.environ.get("PYMOL_DEFER_QUEUE") or cfg.get("PYMOL_DEFER_QUEUE") or "<default>"

    eff_workers = min(workers_cfg, soft_cap) if soft_cap else workers_cfg
    eff_mode    = mode_force if mode_force in ("cli","pymol2") else mode_cfg

    print("[render] replay start | "
          f"workers={workers_cfg} soft_cap={soft_cap or 'none'} -> used={eff_workers} "
          f"mode={mode_cfg} force={mode_force or 'none'} batch={batch_sz or 'none'} "
          f"AB={ab_env or 'none'} queue={queue_path}")

    # Replay (A/B/batching handled inside capture_pose.replay_deferred_jobs_mp)
    capture_pose.replay_deferred_jobs_mp(max_workers=eff_workers, mode=eff_mode)

    score_tol = float(cfg.get("SCORE_TOL_KCAL", args.analysis_score_tol))
    center_tol = float(cfg.get("CENTER_TOL_A", args.analysis_center_tol))
    rmsd_tol = float(cfg.get("RMSD_TOL_A", args.analysis_rmsd_tol))
    out_dir = Path(cfg.get("BENCH_ANALYSIS_OUTDIR") or (args.analysis_out or "")) if (
                cfg.get("BENCH_ANALYSIS_OUTDIR") or args.analysis_out) else None
    # HTML success-only default ON (env ANALYSIS_HTML_SUCCESS_ONLY=1/0)
    html_success_only = (str(os.environ.get("ANALYSIS_HTML_SUCCESS_ONLY", "1")).strip().lower()
                         in ("1", "true", "yes", "on"))

    # since_epoch default: BENCH_ANALYSIS_SINCE_EPOCH or this run's start (BENCH_RUN_START_EPOCH)
    _since_env = os.environ.get("BENCH_ANALYSIS_SINCE_EPOCH") or os.environ.get("BENCH_RUN_START_EPOCH")
    since_epoch = float(_since_env) if _since_env else None

    # ----- Post-run analysis -----
    if args.run_analysis:
        docked_root = Path(cfg["DOCKED_DIR"])
        if not docked_root.is_dir():
            print(f"[analysis] Docked root not found at {docked_root} - skipping.")
        else:
            details, summary = _run_auto_analysis(
                docked_root=Path(cfg["DOCKED_DIR"]),
                mapping_csv=Path(mapping_path),
                only_pdb=args.analysis_only_pdb,
                score_tol=score_tol,
                center_tol=center_tol,
                rmsd_tol=rmsd_tol,
                include_identity=not bool(args.analysis_no_identity),
                out_dir=out_dir,
                since_epoch=since_epoch,
                html_success_only=html_success_only,
            )

            if details and summary:
                print(f"[analysis] ✅ Details: {details}")
                print(f"[analysis] ✅ Summary: {summary}")
            else:
                print("[analysis] ❌ Analysis did not produce outputs.")


if __name__ == "__main__":
    main()
