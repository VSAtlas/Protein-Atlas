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

import sys
if hasattr(sys.stdout, "reconfigure"):  # Py3.7+
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import os
import time
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from tqdm import tqdm

from input_and_export_functions import (
    load_inputs, validate_config, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, generate_config, record_score, score_key,_to_bool
)
from protein_functions import detect_active_site
from activesite import extract_and_remove_ligands
from prep_ligands import prep_ligands_from_pdb, is_valid_ligand
from pose_validation import (
    validate_pose_pdbqt, extract_surface_atoms, attempt_fallback_recenter,
    filter_and_rewrite_poses_by_rmsd, compute_self_rmsd   # <— NEW IMPORT
)
from run_vina import run_docking_task, validate_all_poses

# ======================
# Data models & utilities
# ======================

@dataclass
class RetryManager:
    max_retries: int = 2
    recipes: Dict[str, List[Dict[str, Any]]] = field(default_factory=lambda: {
        # If we docked far from the pocket, try small geometry tweaks — not more modes
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

from dataclasses import dataclass
import time
from pathlib import Path

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
    from pose_validation import validate_pose_pdbqt  # local import to avoid cycles
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
            # If we can't write a temp file, just fall back to validating the whole file once
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

@dataclass
class Paths:
    """File/dir paths relevant to one protein."""
    pdb_id: str
    pdb_path: str
    nolig_pdb_path: str
    ligand_output_dir: Path
    ligands_mol2_dir: Path
    prepped_ligands_dir: Path
    cleaned_pdb_path: Path
    receptor_pdbqt_path: Path


def norm(p: str | Path) -> str:
    """Normalize path to a clean, forward-slash string for logs & keys."""
    return os.path.abspath(str(p)).replace("\\", "/")


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
    logger = logging.getLogger(pdb_id)
    logger.setLevel(logging.DEBUG)

    log_dir = os.path.join(docked_dir, pdb_id)
    os.makedirs(log_dir, exist_ok=True)
    log_file_path = os.path.join(log_dir, "protein.log")

    # Reset handlers to avoid duplicates if re-used
    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    fh = logging.FileHandler(log_file_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    ch = logging.StreamHandler(stream=sys.stdout)  # stdout, not stderr
    env_override = os.environ.get("QUIET_CONSOLE_OVERRIDE", "").strip()
    quiet = (env_override.lower() in {"1", "true", "yes"}) if env_override else _to_bool(cfg.get("QUIET_CONSOLE", False))
    ch.setLevel(logging.WARNING if quiet else logging.INFO)
    ch_level = logging.WARNING if quiet else logging.INFO
    ch.setLevel(ch_level)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.propagate = False
    return logger


def make_paths(cfg: Dict, base_id: str, pdb_file: str) -> Paths:
    pdb_id = base_id.upper()  # keep consistent casing
    root = Path(cfg["OUTPUT_DIR"]) / pdb_id  # e.g., processed_pdbs/1IEP

    # Canonical subfolders
    raw_dir      = root / "raw"
    lig_raw_dir  = root / "ligands_raw"
    nolig_dir    = root / "nolig"
    receptor_dir = root / "receptor"
    work_dir     = root / "work"

    # Ensure they exist
    for d in (raw_dir, lig_raw_dir, nolig_dir, receptor_dir, work_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Inputs / outputs
    pdb_path           = os.path.join(cfg["INPUT_DIR"], pdb_file)
    nolig_pdb_path     = str(nolig_dir / f"{pdb_id}_nolig.pdb")                 # intermediate nolig
    cleaned_pdb_path   = receptor_dir / f"{pdb_id}_cleaned.pdb"                 # final cleaned PDB
    receptor_pdbqt     = receptor_dir / f"{pdb_id}.pdbqt"                       # final receptor PDBQT
    ligand_output_dir  = lig_raw_dir                                            # <— IMPORTANT: was *_cleaned_ligands
    ligands_mol2_dir   = Path(cfg["LIGANDS_MOL2_DIR"]) / pdb_id                 # keep as-is
    prepped_lig_dir    = Path(cfg["OUTPUT_LIGANDS_DIR"]) / pdb_id               # prepped .pdbqt library

    prepped_lig_dir.mkdir(parents=True, exist_ok=True)

    return Paths(
        pdb_id=pdb_id,
        pdb_path=pdb_path,
        nolig_pdb_path=nolig_pdb_path,
        ligand_output_dir=ligand_output_dir,
        ligands_mol2_dir=ligands_mol2_dir,
        prepped_ligands_dir=prepped_lig_dir,
        cleaned_pdb_path=cleaned_pdb_path,
        receptor_pdbqt_path=receptor_pdbqt,
    )


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
        # if none are readable, skip this base
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
    return Path(cfg["DOCKED_DIR"]) / pdb_id / f".ckpt_{stage_name}.json"

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
        out = Path(cfg["DOCKED_DIR"]) / pdb_id / "audit.json"
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
# Phase 1–5: Prep steps
# ======================
def extract_ligands_to_nolig(paths: Paths, logger: logging.Logger) -> Tuple[int, set]:
    malformed_log = paths.ligands_mol2_dir / "malformed_ligands.txt"
    if malformed_log.exists():
        malformed_log.unlink()

    ligands_dict, _ = extract_and_remove_ligands(
        paths.pdb_path, paths.nolig_pdb_path, str(paths.ligand_output_dir)
    )
    logger.info(f"Extracted {len(ligands_dict)} ligands -> {paths.ligand_output_dir}")

    control_stems = set()
    for ext in (".mol2", ".pdb", ".sdf"):
        for p in paths.ligand_output_dir.rglob(f"*{ext}"):
            control_stems.add(Path(p).stem)

    return len(ligands_dict), control_stems
def robust_prepare_controls(paths: Paths, cfg: Dict, logger: logging.Logger) -> None:
    import subprocess
    from pathlib import Path
    from rdkit import Chem

    out_dir = paths.prepped_ligands_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    def sanitized_path(p: Path) -> Path:
        return p.parent / f"{p.stem}.sanitized.pdb"

    def sanitize_pdb(in_pdb: Path, out_pdb: Path) -> bool:
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        wrote_any = False
        with open(in_pdb, "r", encoding="utf-8", errors="ignore") as fin, \
             open(out_pdb, "w", encoding="utf-8") as fout:
            for ln in fin:
                if not ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln); continue
                altloc = ln[16].strip() if len(ln) > 16 else ""
                if altloc and altloc.upper() not in ("A", ""):
                    continue
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                    if any([(x != x), (y != y), (z != z)]) or max(abs(x),abs(y),abs(z)) > 1e6:
                        continue
                except Exception:
                    continue
                fout.write(ln); wrote_any = True
        if not wrote_any:
            logger.warning(f"[control-prep] {in_pdb.name}: no safe ATOM/HETATM lines kept.")
        return out_pdb.exists()

    mgl_py = str(Path(cfg["MGLTOOLS_PYTHON"])) if cfg.get("MGLTOOLS_PYTHON") else None
    prep_lig_script = cfg.get("PREPARE_LIGAND_SCRIPT")
    if not prep_lig_script and cfg.get("PREPARE_RECEPTOR_SCRIPT"):
        prep_lig_script = str(Path(cfg["PREPARE_RECEPTOR_SCRIPT"]).parent / "prepare_ligand4.py")

    obabel_cfg = (cfg.get("OPENBABEL_PATH") or "").strip()
    obabel_exe = obabel_cfg if obabel_cfg.lower().endswith(".exe") else str(Path(obabel_cfg) / "obabel.exe")

    for p in paths.ligand_output_dir.rglob("*.pdb"):
        base = p.stem
        san = sanitized_path(p)
        try:
            if not sanitize_pdb(p, san):
                logger.error(f"[control-prep] Failed to create sanitized: {san}")
                continue

            out_pdbqt = out_dir / f"{base}.pdbqt"
            if (not san.exists()) or (san.stat().st_size == 0):
                logger.error(f"[control-prep] Missing or empty just before MGLTools: {san}")
                continue

            if mgl_py and prep_lig_script and Path(prep_lig_script).exists():
                if not san.exists():
                    logger.error(f"[control-prep] Missing just before MGLTools: {san}")
                    continue
                lig_basename = san.name
                subprocess.check_call([mgl_py, prep_lig_script,
                                       "-l", lig_basename,
                                       "-o", str(out_pdbqt),
                                       "-U", "nphs_lps",
                                       "-A", "checkhydrogens"],
                                      cwd=str(san.parent))
            else:
                # RDKit (no sanitize explosions) + obabel (no gen3D)
                mol = Chem.MolFromPDBFile(str(san), sanitize=False, removeHs=False)
                if mol is None:
                    logger.warning(f"[control-prep] RDKit failed to read {san.name}")
                    continue
                tmp_pdb = san.with_suffix(".tmp.pdb")
                Chem.MolToPDBFile(mol, str(tmp_pdb))
                subprocess.check_call([obabel_exe, "-ipdb", str(tmp_pdb), "-opdbqt", "-O", str(out_pdbqt)])
                tmp_pdb.unlink(missing_ok=True)

            if not out_pdbqt.exists():
                logger.warning(f"[control-prep] Expected output not created: {out_pdbqt}")

        except subprocess.CalledProcessError as e:
            logger.warning(f"[control-prep] Failed for {p.name}: {e}")
        except Exception as e:
            logger.warning(f"[control-prep] Unexpected failure for {p.name}: {e}")
        finally:
            san.unlink(missing_ok=True)


def prepare_receptor(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[Optional[str], Optional[str]]:
    import automate_protein_prep
    from distutils.util import strtobool

    force_reprocess = bool(strtobool(str(cfg.get("FORCE_REPROCESS", False))))
    logger.info(
        f"FORCE_REPROCESS={force_reprocess} | "
        f"cleaned_exists={paths.cleaned_pdb_path.exists()} "
        f"receptor_exists={paths.receptor_pdbqt_path.exists()}"
    )

    if paths.cleaned_pdb_path.exists() and paths.receptor_pdbqt_path.exists() and not force_reprocess:
        logger.info("Reusing existing cleaned PDB and receptor PDBQT.")
        try:
            if bool(cfg.get("RECEPTOR_SANITY_CHECK", True)) and not receptor_sanity_check(str(paths.receptor_pdbqt_path)):
                logger.warning("Receptor sanity check failed (cached receptor).")
                return None, None
        except Exception as _e:
            logger.warning(f"Receptor sanity check skipped due to error: {_e}")
        return norm(paths.cleaned_pdb_path), norm(paths.receptor_pdbqt_path)

    result = automate_protein_prep.main(str(paths.nolig_pdb_path))
    if not result or not isinstance(result, tuple) or len(result) != 2:
        logger.warning("Protein prep failed.")
        return None, None

    cleaned_pdb, receptor_pdbqt = result

    try:
        if Path(receptor_pdbqt).resolve() != paths.receptor_pdbqt_path.resolve():
            from shutil import copy2
            paths.receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
            copy2(receptor_pdbqt, paths.receptor_pdbqt_path)
            receptor_pdbqt = str(paths.receptor_pdbqt_path)
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

from pathlib import Path  # (top-level import; you already have it)
import numpy as np        # (top-level; you already have it)

def detect_pocket(cleaned_pdb: str,
                  ligand_dir: Path,
                  logger: logging.Logger) -> Tuple[
                      Optional[Tuple[float,float,float]],
                      Optional[Tuple[float,float,float]],
                      str  # <— source ("control" or "p2rank" or "none")
                  ]:
    """
    Prefer control ligands for docking center/box. If none, fall back to P2Rank.
    """

    import numpy as np

    def find_control_pdbs(d: Path) -> list[Path]:
        return sorted([p for p in d.glob("*.pdb") if p.is_file()])

    # --- 1) Controls check ---
    ctrl_files = find_control_pdbs(ligand_dir)

    # Also check sibling if this is *_NOLIG
    if not ctrl_files and ligand_dir.parent.name.upper().endswith("_NOLIG"):
        sibling = ligand_dir.parents[1] / ligand_dir.parent.name[:-6] / "ligands_raw"
        if sibling.exists():
            ctrl_files = find_control_pdbs(sibling)
            if ctrl_files:
                logger.info(f"[Control-center] Found controls in sibling: {sibling}")

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
            logger.info(f"[Control-center] Using control centroid {ctrl_center} with box {box_size}")
            return ctrl_center, box_size, "control"

    # --- 2) Fallback to P2Rank ---
    center, box_size = detect_active_site(cleaned_pdb)
    if center:
        box_size = tuple(min(28.0, float(s)) for s in box_size)
        logger.info(f"[P2Rank] Using P2Rank center {center} with box {box_size}")
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

def prepare_and_filter_ligands(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[List[str], Dict[str, int]]:
    # Keep your existing prep step (controls/extracted); harmless if nothing to do
    prep_ligands_from_pdb(
        ligand_output_dir=paths.ligand_output_dir,
        ligands_mol2_dir=paths.ligands_mol2_dir,
        prepped_ligands_dir=paths.prepped_ligands_dir,
    )

    # --- NEW: scan multiple roots ---
    roots = []
    global_root = Path(cfg["OUTPUT_LIGANDS_DIR"])
    if global_root.exists():
        roots.append(global_root)
    # also include the per-protein folder (controls you just prepared land here)
    if paths.prepped_ligands_dir.exists():
        roots.append(paths.prepped_ligands_dir)

    # optional extras (semicolon-separated absolute paths)
    extra = str(cfg.get("LIBRARY_EXTRA_DIRS", "")).strip()
    if extra:
        for d in extra.split(";"):
            d = d.strip()
            if d:
                p = Path(d)
                if p.exists():
                    roots.append(p)

    # Gather all .pdbqt files from all roots (deduped by normalized absolute path)
    logger.info("Scanning for ligands under: " + " | ".join(str(r) for r in roots))
    seen_paths: set[str] = set()
    all_pdbqt_paths: list[Path] = []
    for r in roots:
        for p in r.rglob("*.pdbqt"):
            pn = norm(p)
            if pn not in seen_paths:
                seen_paths.add(pn)
                all_pdbqt_paths.append(p)

    if not all_pdbqt_paths:
        logger.warning("No .pdbqt ligands found in any library roots; nothing to dock.")
        return [], {}

    # Validate each ligand and collect heavy atom counts
    valid_pdbqt: Dict[str, Path] = {}
    valid_count = 0
    for p in all_pdbqt_paths:
        stem = p.stem
        try:
            # use the global root as "library root" for path-based checks
            lib_root_for_checks = str(global_root if global_root.exists() else paths.prepped_ligands_dir.parent)
            if is_valid_ligand(p, lib_root_for_checks):
                valid_pdbqt[norm(p)] = p
                valid_count += 1
            else:
                logger.debug(f"Excluded malformed ligand (pdbqt check failed): {p}")
        except Exception:
            logger.debug(f"Excluded malformed ligand (exception): {p}")

    logger.info(f"Valid .pdbqt ligands (union of roots): {valid_count}")

    # Optional blacklist by name prefix (PAINS_NAMES)
    pains_tokens = [t.strip().upper() for t in str(cfg.get("PAINS_NAMES", "")).split(",") if t.strip()]
    if pains_tokens:
        before = len(valid_pdbqt)
        valid_pdbqt = {
            k: v for k, v in valid_pdbqt.items()
            if not any(Path(k).stem.upper().startswith(tok) for tok in pains_tokens)
        }
        removed = before - len(valid_pdbqt)
        if removed > 0:
            logger.info(f"Name blacklist removed {removed} ligands (PAINS_NAMES).")

    # Heavy atoms per ligand
    heavy_atom_counts: Dict[str, int] = {}
    ligands_to_dock: List[str] = []
    for k, p in valid_pdbqt.items():
        ligands_to_dock.append(k)
        try:
            ha = _count_heavy_atoms_from_pdbqt(p)
        except Exception:
            ha = 0
        heavy_atom_counts[k] = int(ha)

    # Sampling / limits for Stage1
    import random
    sample_n = int(cfg.get("LIBRARY_SAMPLE_N", 0) or 0)
    limit_n  = int(cfg.get("LIBRARY_LIMIT", 0) or 0)

    if sample_n > 0 and sample_n < len(ligands_to_dock):
        ligands_to_dock = random.sample(ligands_to_dock, sample_n)
        heavy_atom_counts = {k: heavy_atom_counts[k] for k in ligands_to_dock}
        logger.info(f"[Debug] Sampling {sample_n} ligands from library for stage1.")

    if limit_n > 0 and limit_n < len(ligands_to_dock):
        ligands_to_dock = ligands_to_dock[:limit_n]
        heavy_atom_counts = {k: heavy_atom_counts[k] for k in ligands_to_dock}
        logger.info(f"[Debug] Limiting library to first {limit_n} ligands for stage1.")

    logger.info(f"Ligands queued for docking (union roots): {len(ligands_to_dock)}")
    return ligands_to_dock, heavy_atom_counts

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
    budget_guards: Optional[Dict[str, BudgetGuard]] = None,  # <<< NEW: external per-ligand guards (Option B)
) -> Tuple[
    Dict[str, float],
    List[str],
    List[float],
    Dict[str, str],
    Dict[str, Tuple[Optional[float], str]]
]:
    from sys import stdout as _stdout

    threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))
    max_workers = int(cfg["MAX_PARALLEL_JOBS"])

    scores: Dict[str, float] = {}
    validated_ligands: List[str] = []
    all_distances: List[float] = []
    raw_docked_ligands: Dict[str, str] = {}
    invalids: Dict[str, Tuple[Optional[float], str]] = {}

    surface_coords = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)

    # ---- helpers (unchanged in spirit) ----
    def _best_pose_pdb_from_pdbqt(pdbqt_path: str, obabel_path: Optional[str] = None) -> Optional[str]:
        try:
            from pathlib import Path
            import subprocess, shutil
            out_pdb = Path(pdbqt_path).with_suffix(".best.pdb")

            obabel = obabel_path or os.environ.get("OPENBABEL_EXE") or ""
            if os.name == "nt" and obabel and not obabel.lower().endswith(".exe"):
                obabel = str(Path(obabel) / "obabel.exe")
            if not obabel or not shutil.which(obabel):
                obabel = shutil.which("obabel")
            if not obabel:
                raise RuntimeError("OpenBabel not found; set OPENBABEL_PATH/OPENBABEL_EXE")

            # first model only; strip H so heavy-atom counts line up
            cmd = [obabel, "-ipdbqt", str(pdbqt_path), "-opdb", "-O", str(out_pdb), "-f", "1", "-l", "1", "-d"]
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return str(out_pdb) if out_pdb.exists() and out_pdb.stat().st_size > 0 else None
        except Exception as e:
            logger.warning(f"[RMSD] OpenBabel conversion failed for {os.path.basename(pdbqt_path)}: {e}")
            return None

    def _validate_with_rmsd_gate(
        lig_path: str,
        lig_name: str,
        out_pdbqt_path: str,
        score_val: float
    ) -> Tuple[bool, Optional[str]]:
        """
        Controls: crystal redock RMSD is a hard gate.
        Non-controls: self-RMSD is logged in upstream code; do not gate here.
        """
        base = Path(lig_path).stem.split("_stage")[0]
        crystal_ref = control_lookup.get(base)

        if crystal_ref:
            best_pdb = _best_pose_pdb_from_pdbqt(out_pdbqt_path, obabel_path=cfg.get("OPENBABEL_PATH"))
            if not best_pdb:
                logger.warning(f"{lig_name} | unable to extract best pose PDB for redock RMSD.")
                return False, "no_best_pose_for_rmsd"

            ok = validate_ligand(
                ligand_name=lig_name,
                docked_path=best_pdb,
                crystal_path=str(crystal_ref),
                rmsd_thresh=float(cfg.get("CONTROL_RMSD_MAX_ANG", 2.0)),
                self_rmsd=None,
                logger=logger
            )
            if ok:
                logger.info(f"{lig_name} | {stage['name']} score: {score_val:.2f} kcal/mol (redock-RMSD PASS)")
                return True, None
            else:
                return False, "rmsd_fail"

        # Non-controls: redock gate not applicable here (geometry checks already passed).
        return True, None

    # ---- scheduling & submission ----
    # default budget used only if caller didn't supply a guard for a ligand
    default_budget_seconds = float(
        cfg.get("MAX_RETRY_SECONDS_PER_LIGAND",
                cfg.get("BENCH_MAX_SECONDS", 300.0))
    )

    submit_queue = []
    guards_for_ligand: Dict[str, BudgetGuard] = {}

    for lig in ligands:
        # pick external guard if provided; else create a local one so code stays robust
        guard = (budget_guards.get(lig) if budget_guards else None)
        if guard is None:
            guard = BudgetGuard(default_budget_seconds)
        guards_for_ligand[lig] = guard

        # if the ligand's budget is already exhausted (e.g., earlier wave), don’t even submit
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

                conf_path, out_path = generate_config(
                    cfg["OVERALL_DIR"], pdb_id, receptor_pdbqt, center, box_size,
                    lig, stage["name"], stage_for_cfg, threads_per_vina
                )

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
                file=_stdout
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

                    # if pose valid → RMSD gate for controls; non-controls accepted
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
                    # stop if we blew the budget for this ligand
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

                                conf_path2, out_path2 = generate_config(
                                    cfg["OVERALL_DIR"], pdb_id, receptor_pdbqt, center, box_size,
                                    lig, stage_retry["name"], stage_retry, threads_per_vina
                                )

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

                    # Structured recipe retries (BudgetGuard enforced each loop)
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
                                    docking_dir=os.path.join(cfg['DOCKED_DIR'], pdb_id),
                                    stage_name=stage_retry2["name"],
                                    pocket_center=center,
                                    logger=logger,
                                    exclude_basenames=set(),
                                )
                                if new_c is not None:
                                    retry_center = new_c
                            if "box_pad_delta" in recipe and isinstance(recipe["box_pad_delta"], (int, float)):
                                dx = float(recipe["box_pad_delta"])
                                retry_box = tuple(min(28.0, s + dx) for s in box_size)
                        except Exception as _e:
                            logger.warning(f"Retry recenter/box tweak failed: {_e}")

                        conf_path3, out_path3 = generate_config(
                            cfg["OVERALL_DIR"], pdb_id, receptor_pdbqt, retry_center, retry_box,
                            lig, stage_retry2["name"], stage_retry2, threads_per_vina
                        )
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
    De-duped and control-anchored: will not fire if (a) a control validated this stage, (b) a global switch already occurred this stage,
    or (c) global switch cap reached.
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
    if evaluated < max(params.EARLY_RECENTER_MIN_EVAL, 15):  # slightly stricter
        logger.info(f"Early recenter skipped: evaluated={evaluated} < threshold.")
        return False, center, box_size, [], attempts_used

    far = sum(1 for d in all_distances if isinstance(d, (int, float)) and d > params.EARLY_RECENTER_FAR_A)
    far_ratio = far / evaluated if evaluated else 0.0
    med_dist = float(np.median(all_distances)) if all_distances else 0.0

    # Prefer a single mild box expand over recenter
    if params.ALLOW_BOX_EXPAND and (0.55 <= far_ratio < params.EARLY_RECENTER_RATIO) and (9.0 <= med_dist < params.EARLY_RECENTER_MEDIAN_A) and (valid_count == 0):
        new_box = tuple(min(28.0, s + 4.0) for s in box_size)
        if new_box != box_size:
            logger.info(f"Borderline far_ratio={far_ratio:.2f}, median={med_dist:.1f} Å -> expand box to {new_box} and redo stage1.")
            # Note: not counted as a global switch
            return True, center, new_box, stage1_original[:], attempts_used

    if (far_ratio >= params.EARLY_RECENTER_RATIO) and (med_dist >= params.EARLY_RECENTER_MEDIAN_A) and (valid_count == 0):
        if attempts_used >= params.MAX_RECENTER_ATTEMPTS:
            logger.warning("Early recenter max attempts reached; proceeding without recenter.")
            return False, center, box_size, [], attempts_used

        logger.warning(f"Early recenter trigger: far_ratio={far_ratio:.2f}, median={med_dist:.1f} Å, valid=0 -> recentering.")
        fb_pose, new_center, _best_score, _chosen = attempt_fallback_recenter(
            fallback_ligands=raw_docked,
            receptor_pdbqt=receptor_pdbqt,
            docking_dir=os.path.join(cfg['DOCKED_DIR'], pdb_id),
            stage_name="stage1",
            pocket_center=center,
            logger=logger,
            exclude_basenames=set(),
        )
        if new_center is not None:
            attempts_used += 1
            new_box = tuple(min(28.0, s) for s in box_size)
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

    # Use provided base if given (e.g., Stage1 pool size) — otherwise fall back to valid-count
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
    fb_pose, new_center, _best_score, _chosen_lig = attempt_fallback_recenter(
        fallback_ligands=raw_docked_ligands,
        receptor_pdbqt=receptor_pdbqt,
        docking_dir=os.path.join(cfg["DOCKED_DIR"], pdb_id),
        stage_name=stage_name,
        pocket_center=center,
        logger=logger,
        exclude_basenames=set(),
    )
    if new_center is None:
        logger.warning("Fallback recovery failed.")
        return False, center, box_size, []

    new_box = tuple(min(28.0, s) for s in box_size)
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
# Phase 7–8: Finalization
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
) -> None:
    if not validated_ligands_last:
        return

    last_stage = stages[-1]["name"]
    final_surface = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)

    if docking_mode == "polypharmacology":
        final_scores = score_history.get(last_stage, {})
        top_ligs = sorted(final_scores.items(), key=score_key)[:20]
        validated_ligands_last = [lig for lig, _ in top_ligs]
        logger.info(f"[Polypharmacology] Selected top {len(validated_ligands_last)} ligands for images/validation.")

    for lig in validated_ligands_last:
        out_path = Path(cfg["DOCKED_DIR"]) / pdb_id / last_stage / f"{Path(lig).stem}_{last_stage}.pdbqt"
        if not out_path.exists():
            logger.warning(f"Pose file not found for {os.path.basename(lig)} — likely filtered earlier.")
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
        pose = Path(cfg["DOCKED_DIR"]) / pdb_id / last_stage / f"{Path(top).stem}_{last_stage}.pdbqt"
        if pose.exists():
            out_prefix = Path(cfg["DOCKED_DIR"]) / pdb_id / "top_pose"
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
            cmd = [str(Path(cfg["PYMOL_PATH"])), "capture_pose.py", cleaned_pdb, str(pose), str(out_prefix)]
            print("Running PyMOL:", cmd)
            res = subprocess.run(cmd, capture_output=True, text=True)
            print("PyMOL stdout:", res.stdout)
            print("PyMOL stderr:", res.stderr)
    except Exception as e:
        logger.warning(f"Screenshot generation failed: {e}")


def _pose_path_for(csv_cfg: Dict, pdb_id: str, stage_name: str, lig_path: str) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    from pathlib import Path
    return str(Path(csv_cfg["DOCKED_DIR"]) / pdb_id / stage_name / f"{Path(lig_path).stem}_{stage_name}.pdbqt")


def write_scores_csv(cfg: Dict, pdb_id: str, score_history: Dict[str, Dict[str, Dict]]) -> str:
    import csv, math

    dock_dir = os.path.join(cfg["DOCKED_DIR"], pdb_id)
    os.makedirs(dock_dir, exist_ok=True)

    # --- Wide summary (unchanged shape) ---
    csv_out_wide = os.path.join(dock_dir, "docking_score_summary.csv")
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
    write_score_summary_to_csv(flat, output_path=csv_out_wide)

    # --- Long format with self_rmsd added ---
    csv_out_long = os.path.join(dock_dir, "docking_score_long.csv")
    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["stage", "ligand", "score", "valid", "reason", "heavy_atoms", "le", "self_rmsd"])

        for stage_name, stage_map in score_history.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(lig)
                score = rec.get("score", None)
                valid = bool(rec.get("valid", False))
                reason = rec.get("reason", "")

                ha = rec.get("heavy_atoms", None)
                le = rec.get("le", None)

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

                writer.writerow([stage_name, lig_key, score_str, int(valid), reason_str, ha_str, le_str, sr_str])

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
from rdkit.Chem import rdMolAlign, rdFMCS



def _read_any_lig(path: str):
    p = Path(path)
    ext = p.suffix.lower()
    m = None
    try:
        if ext in (".pdb", ".ent"):
            m = Chem.MolFromPDBFile(str(p), sanitize=False, removeHs=False)
        elif ext == ".sdf":
            # tolerant read: don't die on weird SDF headers
            sup = Chem.SDMolSupplier(str(p), sanitize=False, removeHs=False, strictParsing=False)
            m = next((x for x in sup if x is not None), None)
        elif ext == ".mol2":
            m = Chem.MolFromMol2File(str(p), sanitize=False, removeHs=False)
        else:
            m = Chem.MolFromPDBFile(str(p), sanitize=False, removeHs=False)
        if m is None:
            return None
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
        return Chem.RemoveHs(m)
    except Exception:
        return None

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
      • If crystal structure available → use redocking RMSD (hard gate).
      • Otherwise (non-controls) → self-RMSD is *log-only* (never reject).
    """
    if crystal_path and Path(crystal_path).exists():
        redock_rmsd = compute_rmsd(crystal_path, docked_path)
        if logger:
            sr = f"{self_rmsd:.2f}" if isinstance(self_rmsd, (int, float)) else "n/a"
            logger.info(f"[validate] {ligand_name}: redock_RMSD={redock_rmsd:.2f} Å, self_RMSD={sr}")
        if redock_rmsd <= rmsd_thresh:
            return True
        else:
            if logger:
                logger.warning(
                    f"[validate] {ligand_name}: redocking failed (RMSD {redock_rmsd:.2f} Å > {rmsd_thresh:.2f})"
                )
            return False

    # Non-controls: log self-RMSD but do not gate on it
    try:
        sr_val = float(self_rmsd) if self_rmsd is not None else None
    except Exception:
        sr_val = None
    if logger:
        sr_txt = f"{sr_val:.2f}" if isinstance(sr_val, (int, float)) else "n/a"
        logger.info(f"[validate] {ligand_name}: self_RMSD={sr_txt} Å (LOG-ONLY)")
    return True



# ======================
# Per-protein driver
# ======================
import subprocess
def process_one_protein(cfg: Dict, pdb_file: str, stages: List[Dict], params: RecenterParams) -> None:
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = base_id.replace("_cleaned", "")
    logger = make_protein_logger(cfg["DOCKED_DIR"], pdb_id, cfg)
    logger.info(f"Processing protein: {pdb_file} (id={pdb_id})")

    paths = make_paths(cfg, base_id, pdb_file)

    # 1) Extract ligands → produce nolig PDB
    _lig_count, control_stems = extract_ligands_to_nolig(paths, logger)
    robust_prepare_controls(paths, cfg, logger)
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
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("Skipping protein due to prep failure.")
        return

    # 3) Pocket detection
    center, box_size, center_source = detect_pocket(cleaned_pdb, paths.ligand_output_dir, logger)
    if center is None:
        return

    # clamp initial box once to keep Vina happy (detect_pocket already caps P2Rank path)
    box_size = tuple(min(28.0, float(s)) for s in box_size)
    logger.info(f"Initial box clamped to {box_size} (cap=28 Å)")

    # explicit console breadcrumb so you don't need to open logs
    try:
        c_print = tuple(round(float(x), 3) for x in center)
        b_print = tuple(round(float(x), 1) for x in box_size)
        print(f"[CENTER] source={center_source} center={c_print} box={b_print}")
    except Exception:
        pass

    # 4) Ligand prep & filtering
    ligands, heavy_atom_counts = prepare_and_filter_ligands(cfg, paths, logger)
    # Force-inject control PDBQTs if they exist on disk but weren't selected
    control_stems_lower = {s.lower() for s in control_stems}
    prepped_control_pdbqts = []
    for root in {paths.prepped_ligands_dir, Path(cfg["OUTPUT_LIGANDS_DIR"])}:
        if not root.exists():
            continue
        for p in root.glob("*.pdbqt"):
            if p.stem.split("_stage")[0].lower() in control_stems_lower:
                prepped_control_pdbqts.append(norm(p))

    # Fallback: if any extracted control stems are missing as PDBQT, try on-the-fly obabel convert
    missing_stems = {s.lower() for s in control_stems} - {Path(p).stem.split("_stage")[0].lower() for p in
                                                          prepped_control_pdbqts}
    if missing_stems:
        logger.warning(
            f"[Controls] {len(missing_stems)} extracted ligand(s) missing as PDBQT; attempting quick obabel convert.")
        obabel = cfg.get("OPENBABEL_PATH", "").strip()
        if obabel and not obabel.lower().endswith(".exe"):
            obabel = str(Path(obabel) / "obabel.exe")
        for p in paths.ligand_output_dir.glob("*.pdb"):
            stem = p.stem.split("_stage")[0].lower()
            if stem not in missing_stems:
                continue
            out_pdbqt = paths.prepped_ligands_dir / f"{p.stem}.pdbqt"
            try:
                subprocess.check_call([obabel, "-ipdb", str(p), "-opdbqt", "-O", str(out_pdbqt)])
                if out_pdbqt.exists():
                    logger.info(f"[Controls] Quick-converted {p.name} -> {out_pdbqt.name}")
            except Exception as e:
                logger.warning(f"[Controls] Quick obabel failed for {p.name}: {e}")

    lig_set = {norm(x) for x in ligands}
    missing_controls = [x for x in prepped_control_pdbqts if norm(x) not in lig_set]
    if missing_controls:
        logger.info(f"[Controls] Adding {len(missing_controls)} prepared control(s) to Stage1.")

        ligands = missing_controls + ligands  # prepend to ensure they’re seen early
        for x in missing_controls:
            if x not in heavy_atom_counts:
                try:
                    ha = _count_heavy_atoms_from_pdbqt(Path(x))
                    heavy_atom_counts[x] = int(ha)
                except Exception:
                    heavy_atom_counts[x] = 0

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
    # Also normalize keys in heavy_atom_counts to match
    heavy_atom_counts = {norm(k): v for k, v in heavy_atom_counts.items()}
    # ---- super-simple: front-load controls at the head of Stage1 ----
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
        logger.info(f"[Controls] Front-loading {len(ctrls)} controls. "
                    f"First wave: {[Path(x).name for x in ligands[:int(cfg.get('MAX_PARALLEL_JOBS', 1))]]}")

    present_ctrls = [Path(l).stem.split("_stage")[0].lower() for l in ligands
                     if Path(l).stem.split("_stage")[0].lower() in control_stems_lower]
    if not present_ctrls:
        logger.warning("[Controls] No control ligands present in Stage1 ligand list — "
                       "self-RMSD/locking will not be possible. (Check prep errors above.)")
    if not ligands:
        logger.warning("No valid ligands after filtering; skipping protein.")
        return

    # 4.5) Initialize center selector & guard
    selector = CenterSelector(cfg, logger, control_stems, heavy_atom_counts, center)
    guard = GlobalCenterGuard(
        max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2))
    )

    stage1_original = ligands[:]

    #cache extracted ligands (from the PDB) that actually exist in the prepped pool
    control_stems_lower = {s.lower() for s in control_stems}
    forced_extracted_for_stage3 = {
        lig for lig in stage1_original
        if Path(lig).stem.split("_stage")[0].lower() in control_stems_lower
    }
    logger.info(f"[Force-carry] Extracted ligands earmarked for Stage3: {len(forced_extracted_for_stage3)}")

    score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
    validated_ligands_last: List[str] = []
    recenter_attempts = 0
    docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()

    retry_mgr = RetryManager()

    i = 0
    while i < len(stages):
        guard.reset_stage()  # only one global switch allowed per stage
        stage = stages[i]

        # Optional checkpoint skip (fingerprinted)
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
        # Stage-1 two-wave run: controls first, then others (so center can update/lock early)
        if i == 0 and ctrls and non_ctrls:
            logger.info(f"Stage1 two-wave: {len(ctrls)} controls first, then {len(non_ctrls)} others.")

            # Wave A — controls only
            s1, v1, d1, rd1, inv1 = run_one_stage(
                cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                ctrls, logger, retry_mgr, control_lookup
            )

            # Let controls influence center/lock immediately
            try:
                dec = selector.consider_switch(stage['name'], s1, v1, rd1, receptor_pdbqt, center, guard)
                if dec.promoted and dec.new_center is not None:
                    old = center
                    center = dec.new_center
                    guard.mark_switch()
                    logger.info(
                        f"[CENTER] Switched before library run: {old} -> {center} ({dec.reason}) [global switch]")
            except Exception as e:
                logger.warning(f"CenterSelector (controls-only) failed gracefully: {e}")
            # Early control-lock using Wave A (controls-only) results
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
                    f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Å); "
                    "future center switches disabled."
                )

            # (Optional) reuse your existing qualified-control locking gate here using s1/v1/rd1.

            # Wave B — non-controls, using (possibly) updated/locked center
            s2, v2, d2, rd2, inv2 = run_one_stage(
                cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                non_ctrls, logger, retry_mgr, control_lookup
            )

            # Merge results to keep the rest of the pipeline unchanged
            scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
            raw_docked = {**rd1, **rd2}
            invalids = {**inv1, **inv2}
        else:
            logger.info(f"Starting {stage['name']} with {len(ligands)} ligands...")
            scores, validated, distances, raw_docked, invalids = run_one_stage(
                cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                ligands, logger, retry_mgr, control_lookup
            )

        validated_ligands_last = validated

        # --- Compute control anchor hit for this stage (used to gate recentering) ---
        def _is_control(lig: str) -> bool:
            stem = Path(lig).stem.split("_stage")[0].lower()
            if stem.upper() in {s.strip().upper() for s in cfg.get("CONTROL_BLACKLIST", "").split(",") if s.strip()}:
                return False
            ha = heavy_atom_counts.get(lig)
            if ha is not None and ha < int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10)):
                return False
            return stem in {s.lower() for s in control_stems}

        control_anchor_hit = any(_is_control(lig) for lig in validated)
        # === Confidence‑gated lock when a validated control anchors the site ===
        # Preconditions already satisfied here:
        #  - controls in `validated` passed self‑RMSD (run_one_stage demotes failures)
        # Confidence criteria we add:
        #  (a) control score ≤ CONTROL_LOCK_SCORE_MAX  (kcal/mol; negative is better)
        #  (b) centroid of the control's best pose is close to current center (≤ CONTROL_LOCK_CENTER_MAX_DIST Å)
        #  (c) at least CONTROL_LOCK_MIN_HITS such controls
        lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", float("inf")))
        lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
        lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))

        # Gather validated controls with score gate and proximity gate
        qualified_controls = []
        for lig in validated:
            if not _is_control(lig):
                continue
            sc = scores.get(lig)
            if sc is None or not np.isfinite(sc):
                continue
            if sc > lock_score_max:
                continue
            # Proximity: centroid of the actually docked pose vs current center
            pose_path = raw_docked.get(lig)
            if not pose_path:
                continue
            c = CenterSelector._pdbqt_centroid(pose_path)
            if c is None:
                continue
            if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                qualified_controls.append(lig)

        if len(qualified_controls) >= lock_min_hits:
            if not guard.locked:
                guard.lock()
                logger.info(
                    "[CONTROL-LOCK] Control(s) validated with strong confidence "
                    f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Ang); "
                    "center is now anchored; future center switches are disabled."
                )

        # ---------- Stage invariant check ----------
        try:
            processed = {norm(x) for x in ligands}
            valid_set = {norm(x) for x in scores.keys()}
            invalid_set = {norm(x) for x in invalids.keys()}
            both = valid_set & invalid_set
            missing = processed - (valid_set | invalid_set)
            if both or missing:
                logger.error(f"Invariant violation at {stage['name']}: both={len(both)}, missing={len(missing)}")
                if both:
                    logger.error("Ligands marked both valid & invalid: " + ", ".join(
                        os.path.basename(x) for x in list(both)[:10]))
                if missing:
                    logger.error(
                        "Ligands missing from results: " + ", ".join(os.path.basename(x) for x in list(missing)[:10]))
        except Exception as _e:
            logger.warning(f"Invariant check failed: {_e}")
        # ------------------------------------------

        # Record both valid and invalid into history for CSV
        for lig, sc in scores.items():
            record_score(score_history, stage['name'], lig, sc, True)
            record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)
        for lig, (sc, reason) in invalids.items():
            record_score(score_history, stage['name'], lig, sc, False, reason=reason)
            record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)

        # Consider switching the center based on clusters/controls (guarded)
        promoted_this_stage = False
        try:
            decision = selector.consider_switch(stage['name'], scores, validated, raw_docked, receptor_pdbqt, center, guard)
            if decision.promoted and decision.new_center is not None:
                old = center
                center = decision.new_center
                promoted_this_stage = True
                guard.mark_switch()  # counts as a global switch
                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    checkpoint_invalidate_from(cfg, paths.pdb_id, stages, start_index=i)
                logger.info(
                    f"[CENTER] Switched from {old} -> {center} ({decision.reason}, SwitchScore={decision.switchscore:.2f}) [global switch]"
                )
        except Exception as e:
            logger.warning(f"CenterSelector failed gracefully: {e}")

        # Stage-1 early recenter / expand box — skip if promoted or control anchored
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

        # Adaptive shrink if we have a tight cluster of valid poses
        try:
            if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                med = float(np.median([d for d in distances if isinstance(d, (int, float))])) if distances else None
                if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                    dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                    min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                    new_box = tuple(max(min_box, s - dec) for s in box_size)
                    if new_box != box_size:
                        logger.info(f"Adaptive shrink: median dist {med:.2f} Å -> box {box_size} -> {new_box}")
                        box_size = new_box
        except Exception as _e:
            logger.warning(f"Adaptive shrink skipped: {_e}")

        # Selection or generic fallback (guarded)
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

            # By default, Stage3 (i==1) in polypharmacology should use Stage1 pool size
            use_stage1_base = (docking_mode == "polypharmacology" and i == 1)

            # --- near-miss rescue (place BEFORE select_ligands_for_next) ---
            rescue = []
            if i < len(stages) - 1:
                for lig, (sc, reason) in invalids.items():
                    if sc is not None and "self_rmsd_" in str(reason).lower() and sc <= float(
                            cfg.get("RESCUE_SELF_RMSD_SCORE_MAX", -8.0)):
                        rescue.append((sc, lig))
                rescue = [lig for _, lig in sorted(rescue)[:int(cfg.get("RESCUE_SELF_RMSD_TOP_N", 10))]]

            # Select by score
            selected = select_ligands_for_next(
                docking_mode,
                i,
                stages,
                scores,
                logger,
                base_pool_n=(len(stage1_original) if use_stage1_base else None),
                force_include=(forced_extracted_for_stage3 if use_stage1_base else None)
            )

            # Union (rescue-first), preserve order, dedupe
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

    # 6) Final pose validation & screenshots
    final_pose_validation_and_screenshots(
        cfg, paths.pdb_id, stages, receptor_pdbqt, center, validated_ligands_last,
        score_history, cleaned_pdb, docking_mode, logger
    )

    # 7) Write scores CSV
    csv_path = write_scores_csv(cfg, paths.pdb_id, score_history)
    logger.info(f"Score summary written: {csv_path}")

    # 8) Audit JSON summary
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
        }
        _write_audit_json(cfg, paths.pdb_id, summary)
    except Exception as _e:
        logger.warning(f"Audit JSON write failed: {_e}")


# ======================
# Program entry point
# ======================
def main() -> None:
    print("MODELLER is working with license.")

    cfg = load_inputs()
    validate_config(cfg)

    # --- Center selection knobs (safe defaults) ---
    cfg.setdefault("CENTER_MODE", "control-first")  # ["control-first","hybrid","library-first"]
    cfg.setdefault("CONTROL_BLACKLIST", "GOL,EDO,PG4,MPD,ACT,SO4,PO4,CL,NA,CA")
    cfg.setdefault("CONTROL_MIN_HEAVY_ATOMS", 10)
    cfg.setdefault("CONTROL_ANCHOR_MIN_VALID_RATE", 0.10)  # if current cluster has control hits + ≥10% valid, anchor
    cfg.setdefault("ALLOW_SWITCH_FROM_CONTROL", True)
    cfg.setdefault("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False)
    cfg.setdefault("SWITCH_AWAY_FROM_CONTROL_MIN_BOOST", 2.5)  # kcal/mol median boost needed to leave control
    # Optional lock score gate (kcal/mol). Use a large positive number (or remove) to lock on RMSD alone.
    cfg.setdefault("CONTROL_LOCK_SCORE_MAX", -6.0)
    cfg.setdefault("CONTROL_LOCK_MIN_HITS", 1)  # require ≥ this many validated controls
    cfg.setdefault("CONTROL_LOCK_CENTER_MAX_DIST", 4.0)  # Å; control centroid must be within this of center

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

    pdb_files = [
        f for f in os.listdir(cfg["INPUT_DIR"])
        if f.endswith(".pdb") and "_nolig" not in f.lower()
    ]
    print("Working directory:", os.getcwd())
    print("Loaded config keys:", list(cfg.keys()))
    print(f"Proteins queued: {len(pdb_files)}")

    start = time.time()
    with tqdm(total=len(pdb_files), desc="Processing Proteins", unit="protein") as bar:
        for pdb_file in pdb_files:
            process_one_protein(cfg, pdb_file, stages, params)
            bar.update(1)

    elapsed_min = (time.time() - start) / 60.0
    print(f"\nAll proteins processed in {elapsed_min:.2f} minutes.")


if __name__ == "__main__":
    main()
