"""
High-level pipeline for multi-stage docking.

Phases per protein:
1) Setup & logging
2) Extract ligands and generate a ligand-free PDB
3) Protein preparation (clean PDB + receptor PDBQT), reuse if cached
4) Active-site detection (center, box size)
5) Ligand preparation & filtering
6) Multi-stage docking with early/fallback recenter heuristics + CenterSelector
7) Final pose validation & optional screenshots
8) Write per-protein score CSV
"""

from __future__ import annotations

import sys
if hasattr(sys.stdout, "reconfigure"):  # Py3.7+
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import os
import time
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from tqdm import tqdm

from input_and_export_functions import (
    load_inputs, validate_config, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, generate_config, record_score, score_key
)
from protein_functions import detect_active_site
from activesite import extract_and_remove_ligands
from prep_ligands import prep_ligands_from_pdb, is_valid_ligand
from pose_validation import (
    validate_pose_pdbqt, extract_surface_atoms, attempt_fallback_recenter, filter_and_rewrite_poses_by_rmsd
)
from run_vina import run_docking_task, validate_all_poses


# ======================
# Data models & utilities
# ======================
from dataclasses import dataclass, field
from typing import Dict, Any, List

@dataclass
class RetryManager:
    max_retries: int = 2
    recipes: Dict[str, List[Dict[str, Any]]] = field(default_factory=lambda: {
        "too_far_from_pocket": [
            {"recenter": True, "box_pad_delta": +2.0},
            {"recenter": True, "box_pad_delta": +3.0, "exhaustiveness": 6},
        ],
        "no_valid_pose": [
            {"exhaustiveness": 6, "num_modes": 12},
            {"exhaustiveness": 4, "num_modes": 16, "seed_jitter": True},
        ],
        "timeout": [
            {"num_modes": 4},
            {"exhaustiveness": 4, "num_modes": 3},
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
class RecenterParams:
    """Thresholds for early/fallback recenter heuristics."""
    EARLY_RECENTER_RATIO: float = 0.70
    EARLY_RECENTER_MIN_EVAL: int = 10
    EARLY_RECENTER_FAR_A: float = 15.0
    EARLY_RECENTER_MEDIAN_A: float = 10.0
    ALLOW_BOX_EXPAND: bool = True
    MAX_RECENTER_ATTEMPTS: int = 3


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
        MAX_RECENTER_ATTEMPTS=int(cfg.get("MAX_RECENTER_ATTEMPTS", 3)),
    )


def make_protein_logger(docked_dir: str, pdb_id: str) -> logging.Logger:
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

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.propagate = False
    return logger


def make_paths(cfg: Dict, base_id: str, pdb_file: str) -> Paths:
    """
    Build common paths for a protein and ensure necessary folders exist.
    This keeps path construction in one place for readability.
    """
    pdb_id = base_id
    ligand_output_dir = Path(cfg["OUTPUT_DIR"]) / pdb_id / f"{pdb_id}_cleaned_ligands"
    ligands_mol2_dir = Path(cfg["LIGANDS_MOL2_DIR"]) / pdb_id
    prepped_ligands_dir = Path(cfg["OUTPUT_LIGANDS_DIR"]) / pdb_id
    protein_dir = Path(cfg["OUTPUT_DIR"]) / f"{base_id}_nolig"
    cleaned_pdb_path = protein_dir / f"{base_id}_nolig_cleaned.pdb"
    receptor_pdbqt_path = Path(cfg["PDBQT_DIR"]) / f"{base_id}_receptor.pdbqt"
    pdb_path = os.path.join(cfg["INPUT_DIR"], pdb_file)
    nolig_pdb_path = os.path.join(cfg["OUTPUT_DIR"], f"{base_id}_nolig.pdb")

    ligand_output_dir.mkdir(parents=True, exist_ok=True)
    prepped_ligands_dir.mkdir(parents=True, exist_ok=True)
    protein_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)

    return Paths(
        pdb_id=pdb_id,
        pdb_path=pdb_path,
        nolig_pdb_path=nolig_pdb_path,
        ligand_output_dir=ligand_output_dir,
        ligands_mol2_dir=ligands_mol2_dir,
        prepped_ligands_dir=prepped_ligands_dir,
        cleaned_pdb_path=cleaned_pdb_path,
        receptor_pdbqt_path=receptor_pdbqt_path,
    )


# --------- extra helpers (easy-win features) ---------

# Map validator/reported reasons to retry categories understood by RetryManager
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
    # round to 0.1 Å to avoid float noise triggering spurious invalidation
    return tuple(None if (x is None) else round(float(x), ndp) for x in t)


def _fingerprint_stage(cfg: Dict,
                       receptor_pdbqt: str,
                       center: Tuple[float, float, float],
                       box_size: Tuple[float, float, float],
                       stage: Dict) -> Dict[str, Any]:
    """
    Minimal but robust fingerprint that defines the 'work' of a stage.
    If any of this changes, the stage must be re-run.
    """
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
        "version_tag": "ckpt_v2",  # bump if fingerprint schema changes
    }


def _checkpoint_path(cfg: Dict, pdb_id: str, stage_name: str) -> Path:
    return Path(cfg["DOCKED_DIR"]) / pdb_id / f".ckpt_{stage_name}.json"


def checkpoint_should_skip(cfg: Dict,
                           pdb_id: str,
                           stage_name: str,
                           fingerprint: Dict[str, Any]) -> bool:
    """
    Return True if an existing checkpoint matches the current fingerprint closely enough.
    """
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
    """
    If we change center/box mid-run (e.g., via fallback recenter or CenterSelector),
    invalidate checkpoints for this and all later stages to avoid stale skips.
    """
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
    """Quick guard against empty/degenerate receptor pdbqt files."""
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
    """
    Extract ligands and output a ligand-free PDB at `nolig_pdb_path`.
    Returns (num_extracted, control_stems) where control_stems are names/stems of extracted ligands.
    """
    malformed_log = paths.ligands_mol2_dir / "malformed_ligands.txt"
    if malformed_log.exists():
        malformed_log.unlink()

    ligands_dict, _ = extract_and_remove_ligands(
        paths.pdb_path, paths.nolig_pdb_path, str(paths.ligand_output_dir)
    )
    logger.info(f"Extracted {len(ligands_dict)} ligands -> {paths.ligand_output_dir}")

    # Build a set of control ligand stems based on what was written out
    control_stems = set()
    for ext in (".mol2", ".pdb", ".sdf"):
        for p in paths.ligand_output_dir.rglob(f"*{ext}"):
            control_stems.add(Path(p).stem)

    return len(ligands_dict), control_stems


def prepare_receptor(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[Optional[str], Optional[str]]:
    """
    Create/reuse:
      - cleaned PDB without ligands (Phenix/Reduce/etc. via automate_protein_prep)
      - receptor PDBQT (copied into canonical PDBQT_DIR if generated elsewhere)
    Returns normalized (cleaned_pdb, receptor_pdbqt) or (None, None) on failure.
    """
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
        # Sanity check receptor (optional even when reusing)
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

    # Ensure receptor is placed in canonical PDBQT_DIR
    try:
        if Path(receptor_pdbqt).resolve() != paths.receptor_pdbqt_path.resolve():
            from shutil import copy2
            paths.receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
            copy2(receptor_pdbqt, paths.receptor_pdbqt_path)
            receptor_pdbqt = str(paths.receptor_pdbqt_path)
    except Exception as e:
        logger.warning(f"Could not relocate receptor PDBQT: {e}")

    # Sanity check receptor (optional)
    try:
        if bool(cfg.get("RECEPTOR_SANITY_CHECK", True)):
            ok = receptor_sanity_check(receptor_pdbqt)
            if not ok:
                logger.warning("Receptor sanity check failed (too few atoms or zero coords).")
                return None, None
    except Exception as _e:
        logger.warning(f"Receptor sanity check skipped due to error: {_e}")

    return norm(cleaned_pdb), norm(receptor_pdbqt)


def detect_pocket(cleaned_pdb: str, logger: logging.Logger) -> Tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]]]:
    """
    Detect binding-site center and box size from the cleaned PDB.
    Returns (center, box_size) or (None, None) on failure.
    """
    center, box_size = detect_active_site(cleaned_pdb)
    if center is None:
        logger.warning("Active-site detection failed.")
    return center, box_size


def _count_heavy_atoms_from_pdbqt(pdbqt_path: Path) -> int:
    """
    Count non-hydrogen atoms from a PDBQT file without RDKit.
    Uses element column if present, otherwise falls back to atom name.
    """
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
                # fallback: parse atom name field
                atom_name = line[12:16].strip()
                if not atom_name.upper().startswith("H"):
                    heavy += 1
    return heavy


def prepare_and_filter_ligands(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[List[str], Dict[str, int]]:
    """
    Prepare ligands and build the docking set from the global prepped library:
      - Preps any ligands extracted from this protein (kept for 'controls')
      - Collects ALL valid *.pdbqt ligands under OUTPUT_LIGANDS_DIR recursively
      - Validates basic PDBQT structure with is_valid_ligand
      - Computes heavy-atom counts directly from PDBQT (no RDKit dependency)
      - Optional sampling/limiting via config to control runtime in stage1

    Returns:
      ligands_to_dock: [normalized pdbqt paths]
      heavy_atom_counts: {normalized pdbqt path -> int heavy atoms}
    """
    # Still prep per-protein extracted ligands (controls etc.)
    prep_ligands_from_pdb(
        ligand_output_dir=paths.ligand_output_dir,
        ligands_mol2_dir=paths.ligands_mol2_dir,
        prepped_ligands_dir=paths.prepped_ligands_dir,
    )

    # 1) Collect the entire prepped library (*.pdbqt) globally
    lib_root = Path(cfg["OUTPUT_LIGANDS_DIR"])
    all_pdbqt_paths = list(lib_root.rglob("*.pdbqt"))
    if not all_pdbqt_paths:
        logger.warning("No .pdbqt ligands found under OUTPUT_LIGANDS_DIR; nothing to dock.")
        return [], {}

    # 2) Validate and deduplicate by stem
    seen_stems: set = set()
    valid_pdbqt: Dict[str, Path] = {}
    valid_count = 0
    for p in all_pdbqt_paths:
        stem = p.stem
        if stem in seen_stems:
            continue
        if is_valid_ligand(p, cfg["OUTPUT_LIGANDS_DIR"]):
            valid_pdbqt[stem] = p
            valid_count += 1
        else:
            logger.debug(f"Excluded malformed ligand (pdbqt check failed): {p}")
        seen_stems.add(stem)

    logger.info(f"Valid .pdbqt ligands (global): {valid_count}")

    # 2.5) Optional light name blacklist based on stems
    pains_tokens = [t.strip().upper() for t in str(cfg.get("PAINS_NAMES","")).split(",") if t.strip()]
    if pains_tokens:
        before = len(valid_pdbqt)
        valid_pdbqt = {
            stem: p for stem, p in valid_pdbqt.items()
            if not any(stem.upper().startswith(tok) for tok in pains_tokens)
        }
        removed = before - len(valid_pdbqt)
        if removed > 0:
            logger.info(f"Name blacklist removed {removed} ligands (PAINS_NAMES).")

    # 3) Build docking list + heavy atom counts (from PDBQT directly)
    ligands_to_dock: List[str] = []
    heavy_atom_counts: Dict[str, int] = {}
    for stem, p in valid_pdbqt.items():
        lig_norm = norm(p)
        try:
            ha = _count_heavy_atoms_from_pdbqt(p)
        except Exception:
            ha = 0
        ligands_to_dock.append(lig_norm)
        heavy_atom_counts[lig_norm] = int(ha)

    # 4) Optional runtime control knobs for stage1 (purely convenience)
    import random
    sample_n = int(cfg.get("LIBRARY_SAMPLE_N", 0) or 0)
    limit_n  = int(cfg.get("LIBRARY_LIMIT", 0) or 0)

    if sample_n > 0 and sample_n < len(ligands_to_dock):
        sampled = random.sample(ligands_to_dock, sample_n)
        ligands_to_dock = sampled
        heavy_atom_counts = {k: heavy_atom_counts[k] for k in ligands_to_dock}
        logger.info(f"[Debug] Sampling {sample_n} ligands from library for stage1.")

    if limit_n > 0 and limit_n < len(ligands_to_dock):
        ligands_to_dock = ligands_to_dock[:limit_n]
        heavy_atom_counts = {k: heavy_atom_counts[k] for k in ligands_to_dock}
        logger.info(f"[Debug] Limiting library to first {limit_n} ligands for stage1.")

    logger.info(f"Ligands queued for docking: {len(ligands_to_dock)}")
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
    retry_mgr: RetryManager,  # <-- added
) -> Tuple[
    Dict[str, float],                 # scores (valid)
    List[str],                        # validated_ligands
    List[float],                      # all_distances
    Dict[str, str],                   # raw_docked_ligands {lig -> pose_path}
    Dict[str, Tuple[Optional[float], str]]  # invalids {lig -> (score_or_None, reason)}
]:
    """
    Run one docking stage in parallel:
    - Generates configs (THREADS_PER_VINA default 1 to avoid oversubscription)
    - Runs Vina tasks
    - Validates poses and returns:
        scores: {ligand_path -> float score} for valid poses
        validated_ligands: subset with valid poses
        all_distances: list of distances to pocket for heuristic checks
        raw_docked_ligands: {ligand_path -> out_pose_path}
        invalids: {ligand_path -> (score or None, reason)}
    """
    threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))
    max_workers = int(cfg["MAX_PARALLEL_JOBS"])

    scores: Dict[str, float] = {}
    validated_ligands: List[str] = []
    all_distances: List[float] = []
    raw_docked_ligands: Dict[str, str] = {}
    invalids: Dict[str, Tuple[Optional[float], str]] = {}

    surface_coords = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {}
        for lig in ligands:
            conf_path, out_path = generate_config(
                cfg["OVERALL_DIR"], pdb_id, receptor_pdbqt, center, box_size,
                lig, stage["name"], stage, threads_per_vina
            )
            lig_n, out_n = norm(lig), norm(out_path)
            raw_docked_ligands[lig_n] = out_n
            futures[pool.submit(run_docking_task, cfg["VINA_EXE"], conf_path, lig, out_path)] = (lig_n, out_n)

        # --- per-stage audit hook (first 10 ligands) ---
        try:
            audit_path = Path(cfg["DOCKED_DIR"]) / pdb_id / f"{stage['name']}_ligand_audit.txt"
            audit_path.parent.mkdir(parents=True, exist_ok=True)
            with open(audit_path, "w", encoding="utf-8") as af:
                af.write(f"Total ligands queued: {len(ligands)}\n")
                af.write("First 10 queued:\n")
                for ln in ligands[:10]:
                    af.write(f" - {os.path.basename(ln)}\n")
        except Exception as _e:
            logger.warning(f"Audit hook failed: {_e}")
        # ------------------------------------------------

        with tqdm(total=len(futures), desc=f"Docking ({stage['name']})", unit="ligand") as pbar:
            for fut in as_completed(futures):
                lig, out_path = futures[fut]
                lig_name = os.path.basename(lig)
                try:
                    _, score = fut.result()
                except Exception as e:
                    logger.warning(f"Docking crashed for {lig_name}: {e}")
                    invalids[lig] = (None, "docking_exception")
                    pbar.update(1)
                    continue

                if score is None:
                    logger.warning(f"No score for {lig_name}")
                    invalids[lig] = (None, "no_score")
                    pbar.update(1)
                    continue

                # --- prune near-identical poses in-place ---
                try:
                    kept, removed = filter_and_rewrite_poses_by_rmsd(
                        out_path,
                        rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                        max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
                    )
                    if removed > 0:
                        logger.info(
                            f"{os.path.basename(out_path)}: RMSD filter kept {kept}, removed {removed} duplicate poses.")
                except Exception as e:
                    logger.warning(f"RMSD filtering failed for {os.path.basename(out_path)}: {e}")
                # -----------------------------------------------

                result = validate_pose_pdbqt(
                    protein_pdbqt=receptor_pdbqt,
                    ligand_pdbqt=out_path,
                    pocket_center=center,
                    clash_threshold=2.0,
                    CLASH_TOLERANCE=3,
                    DIST_THRESHOLD_SURFACE=6.0,
                    DIST_THRESHOLD_CENTROID=4.5,
                    surface_atom_coords=surface_coords
                )
                dist = result.get("distance_to_pocket")
                if isinstance(dist, (int, float)):
                    all_distances.append(dist)

                logger.info(f"{lig_name} validation: {result}")
                if result.get("valid", False):
                    scores[lig] = float(score)
                    validated_ligands.append(lig)
                    print(f"{lig_name} | {stage['name']} score: {score:.2f} kcal/mol (valid)")
                    pbar.update(1)
                else:
                    # ---------- optional near-miss rescue (single retry) ----------
                    did_retry = False
                    if bool(cfg.get("RETRY_NEAR_MISS", True)):
                        try:
                            dthr = float(cfg.get("RETRY_DIST_THRESH", 7.0))
                            sthr = float(cfg.get("RETRY_SCORE_THRESH", -7.5))
                            if (dist is not None) and (dist < dthr) and (float(score) <= sthr):
                                did_retry = True
                                stage_retry = dict(stage)
                                stage_retry["name"] = f"{stage['name']}_retry"
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
                                # run docking again synchronously
                                _, score2 = run_docking_task(cfg["VINA_EXE"], conf_path2, lig, out_path2)

                                # prune & re-validate the retry pose
                                try:
                                    filter_and_rewrite_poses_by_rmsd(
                                        out_path2,
                                        rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                                        max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
                                    )
                                except Exception:
                                    pass

                                result2 = validate_pose_pdbqt(
                                    protein_pdbqt=receptor_pdbqt,
                                    ligand_pdbqt=out_path2,
                                    pocket_center=center,
                                    clash_threshold=2.0,
                                    CLASH_TOLERANCE=3,
                                    DIST_THRESHOLD_SURFACE=6.0,
                                    DIST_THRESHOLD_CENTROID=4.5,
                                    surface_atom_coords=surface_coords
                                )

                                if result2.get("valid", False):
                                    scores[lig] = float(score2)
                                    validated_ligands.append(lig)
                                    raw_docked_ligands[lig] = norm(out_path2)  # point to rescued pose
                                    print(f"{lig_name} | {stage_retry['name']} score: {score2:.2f} kcal/mol (rescued)")
                                    pbar.update(1)
                                    continue  # move on to next ligand
                        except Exception as _e:
                            logger.warning(f"Retry path failed for {lig_name}: {_e}")

                    # ---------- failure queue retries (category-driven) ----------
                    err_cat = _map_reason_to_category(result.get("reason", ""))
                    attempt = 0
                    retained_invalid = True  # assume failure stays till we flip it

                    while attempt < retry_mgr.max_retries:
                        recipe = retry_mgr.apply(stage, err_cat, attempt)
                        if not recipe:
                            break

                        # Compose a ligand-specific retry stage
                        stage_retry2 = dict(stage)
                        stage_retry2["name"] = f"{stage['name']}_r{attempt+1}"

                        # Mutate docking params from recipe
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

                        # Per-ligand recenter/box tweak (doesn't affect other ligands)
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

                        # Generate retry config & re-dock synchronously for this ligand
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

                        # prune & re-validate
                        try:
                            filter_and_rewrite_poses_by_rmsd(
                                out_path3,
                                rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                                max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
                            )
                        except Exception:
                            pass

                        result_r = validate_pose_pdbqt(
                            protein_pdbqt=receptor_pdbqt,
                            ligand_pdbqt=out_path3,
                            pocket_center=retry_center,
                            clash_threshold=2.0,
                            CLASH_TOLERANCE=3,
                            DIST_THRESHOLD_SURFACE=6.0,
                            DIST_THRESHOLD_CENTROID=4.5,
                            surface_atom_coords=surface_coords
                        )

                        logger.info(f"{lig_name} retry#{attempt+1} ({err_cat}) -> {result_r}")
                        if result_r.get("valid", False) and score_r is not None:
                            scores[lig] = float(score_r)
                            validated_ligands.append(lig)
                            raw_docked_ligands[lig] = norm(out_path3)
                            print(f"{lig_name} | {stage_retry2['name']} score: {score_r:.2f} kcal/mol (retry rescued)")
                            retained_invalid = False
                            break

                        attempt += 1

                    # Record final failure if still invalid after retries
                    if retained_invalid:
                        invalids[lig] = (float(score), result.get("reason", "pose_invalid"))
                        print(f"{lig_name} | pose invalid (after retries)")
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
) -> Tuple[bool, Tuple[float, float, float], Tuple[float, float, float], List[str], int]:
    """
    Stage-1 heuristic for expanding box or recentering when everything docks far from the pocket.
    Returns:
      (should_restart_stage1, new_center, new_box_size, redock_list, attempts_used)
    """
    if i != 0:
        return False, center, box_size, [], attempts_used

    evaluated = len(all_distances)
    valid_count = len(scores)
    if evaluated < params.EARLY_RECENTER_MIN_EVAL:
        logger.info(f"Early recenter skipped: evaluated={evaluated} < {params.EARLY_RECENTER_MIN_EVAL}.")
        return False, center, box_size, [], attempts_used

    far = sum(1 for d in all_distances if isinstance(d, (int, float)) and d > params.EARLY_RECENTER_FAR_A)
    far_ratio = far / evaluated if evaluated else 0.0
    med_dist = float(np.median(all_distances)) if all_distances else 0.0

    # Expand box slightly if borderline and no valids
    if params.ALLOW_BOX_EXPAND and (0.5 <= far_ratio < params.EARLY_RECENTER_RATIO) and (8.0 <= med_dist < params.EARLY_RECENTER_MEDIAN_A) and (valid_count == 0):
        new_box = tuple(min(28.0, s + 4.0) for s in box_size)
        if new_box != box_size:
            logger.info(f"Borderline far_ratio={far_ratio:.2f}, median={med_dist:.1f} Å -> expand box to {new_box} and redo stage1.")
            return True, center, new_box, stage1_original[:], attempts_used

    # Recenter if clearly off and no valids
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
            logger.info("Re-running stage1 with new center and tightened box.")
            return True, new_center, new_box, stage1_original[:], attempts_used
        logger.warning("Fallback could not produce a new center; proceeding without recenter.")

    return False, center, box_size, [], attempts_used


def select_ligands_for_next(docking_mode: str, i: int, stages: List[Dict], scores: Dict[str, float], logger: logging.Logger) -> List[str]:
    """
    Choose a fraction of ligands to carry forward, based on docking mode schedule.
    Returns a list of ligand paths for the next stage.
    """
    if not scores:
        return []

    schedule = {
        "discovery": [1.0, 0.1, 0.01, 0.001, 0.001],
        "polypharmacology": [1.0, 0.05, 0.005],
    }.get(docking_mode, [1.0] * len(stages))

    pct = schedule[i + 1] if i + 1 < len(schedule) else 0.01
    k = max(1, int(len(scores) * pct))
    next_list = [l for l, _ in sorted(scores.items(), key=lambda kv: kv[1])[:k]]
    logger.info(f"Selected top {len(next_list)} ligands ({pct * 100:.5f}%) for the next stage.")
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
) -> Tuple[bool, Tuple[float, float, float], Tuple[float, float, float], List[str]]:
    """
    If a stage yields no valid ligands, attempt a fallback recenter and restart stage1.
    Returns (should_restart_stage1, new_center, new_box_size, ligands_for_stage1).
    """
    if scores:
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
    logger.info("Re-running stage1 with new center after no-valid fallback.")
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
    Uses clustering of valid pose centroids + simple SwitchScore with hysteresis.
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
        self.hysteresis = float(cfg.get("SWITCH_SCORE_HYSTERESIS", 0.5))
        self.switch_history = []  # keep last few SwitchScores
        self.current_center = np.array(initial_center, float)

        # cache current-center cluster stats (median score, LE, valid_rate) after each stage
        self.curr_stats = {"median_score": None, "median_le": None, "valid_rate": 0.0}

    # ---------- utilities ----------
    @staticmethod
    def _strip_stage(name: str) -> str:
        # e.g., "GOL_A401_stage1.pdbqt" -> "GOL_A401"
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
        # skip solvents/buffers/ions
        return any(stem_upper.startswith(bad) for bad in self.ctrl_blacklist)

    def _is_control(self, lig_path: str) -> bool:
        stem = self._strip_stage(os.path.basename(lig_path))
        if stem.upper() and self._is_blacklisted_control_name(stem.upper()):
            return False
        # Enforce minimum heavy atoms for controls to avoid tiny crystallization junk
        ha = self.heavy.get(lig_path)
        if ha is not None and ha < self.ctrl_min_heavy:
            return False
        return (stem in self.control_stems)

    def _cluster(self, points: List[np.ndarray]) -> List[List[int]]:
        # simple greedy clustering by distance threshold self.eps
        clusters: List[List[int]] = []
        for i, p in enumerate(points):
            placed = False
            for cl in clusters:
                c = np.mean([points[j] for j in cl], axis=0)  # cluster centroid
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
        # ligand efficiency per member
        les = []
        for m in members:
            s = scores_map.get(m)
            if s is None:
                continue
            ha = self.heavy.get(m)
            if ha and ha > 0:
                les.append((-s) / ha)  # s is negative kcal/mol
        center = np.mean([centroids[i] for i in member_idxs], axis=0)
        valid_rate = len(members) / max(1, total_docked)  # fraction of all docked ligands that landed here
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
        # normalize to [0..1] ranges
        consensus_boost = max(0.0, min(1.0, score_boost_kcal / 3.0))  # ≥3 kcal/mol -> ~1.0
        le_norm = max(0.0, min(1.0, le_gain / 0.05))                  # ≥0.05 kcal/mol/HA -> ~1.0
        return 0.40*valid_rate + 0.20*consensus_boost + 0.20*pocketability + 0.20*le_norm

    # ---------- main entry ----------
    def consider_switch(self,
                        stage_name: str,
                        scores: Dict[str, float],
                        validated_ligands: List[str],
                        raw_docked: Dict[str, str],
                        receptor_pdbqt: str,
                        current_center: Tuple[float, float, float]) -> CenterDecision:
        if not validated_ligands:
            return CenterDecision(None, "no_validated_poses")

        # collect centroids for validated ligands
        ligs = sorted(set(validated_ligands))
        centroids = []
        ligs_kept = []
        for lig in ligs:
            pose = raw_docked.get(lig)
            c = self._pdbqt_centroid(pose) if pose else None
            if c is not None:
                centroids.append(c); ligs_kept.append(lig)

        if len(centroids) < 3:
            return CenterDecision(None, "too_few_centroids")

        clusters = self._cluster(centroids)
        total_docked = max(1, len(raw_docked))

        # compute stats per cluster
        stats = [self._cluster_stats(cl, ligs_kept, scores, centroids, total_docked) for cl in clusters]

        # identify cluster overlapping current center (if any)
        curr_idx = None
        curr_center = np.array(current_center, float)
        for idx, st in enumerate(stats):
            if np.linalg.norm(np.array(st["center"]) - curr_center) <= self.eps:
                curr_idx = idx; break

        curr_median = stats[curr_idx]["median_score"] if curr_idx is not None else None
        curr_le = stats[curr_idx]["median_le"] if curr_idx is not None else None
        self.curr_stats = {
            "median_score": curr_median,
            "median_le": curr_le,
            "valid_rate": stats[curr_idx]["valid_rate"] if curr_idx is not None else 0.0
        }

        # control anchoring: if any cluster contains validated controls, prefer that cluster
        control_clusters = [st for st in stats if st["control_hits"] > 0]
        if control_clusters and self.mode in ("control-first","hybrid"):
            ctrl_best = sorted(control_clusters, key=lambda st: (-st["control_hits"], -st["valid_rate"]))[0]
            if np.linalg.norm(np.array(ctrl_best["center"]) - curr_center) > self.eps:
                if self.allow_switch_from_control:
                    self.switch_history.append(1.0)
                    return CenterDecision(new_center=tuple(ctrl_best["center"]),
                                          reason="anchor_to_control_cluster",
                                          switchscore=1.0,
                                          promoted=True)
            if curr_idx is not None and curr_idx == stats.index(ctrl_best):
                if self.require_control_failure:
                    return CenterDecision(None, "control_validated—no_switch")

        # choose best non-current cluster by SwitchScore
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
            if score_boost < float(self.cfg.get("SWITCH_SCORE_IMPROVE_MIN", 1.5)):
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
    """
    Validate all poses for the last-stage ligands; if a best valid pose exists, record it.
    Render top pose screenshots via PyMOL.
    """
    if not validated_ligands_last:
        return

    last_stage = stages[-1]["name"]
    final_surface = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)

    # In polypharmacology mode, select top-N before images to keep artifacts focused
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
        # Deduplicate final-stage poses before deep validation
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

    # Screenshots (best ligand only)
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


def write_scores_csv(cfg: Dict, pdb_id: str, score_history: Dict[str, Dict[str, Dict]]) -> str:
    """
    Write two CSVs into DOCKED_DIR/<pdb_id>/:
      1) docking_score_summary.csv  (wide, backwards-compatible with your pipeline)
      2) docking_score_long.csv     (long-form with score, valid, reason, heavy_atoms, le)

    Returns the path of the wide CSV for convenience.
    """
    import csv

    dock_dir = os.path.join(cfg["DOCKED_DIR"], pdb_id)
    os.makedirs(dock_dir, exist_ok=True)

    # -------- Wide summary (backwards compatible) --------
    csv_out_wide = os.path.join(dock_dir, "docking_score_summary.csv")
    flat = {}
    for stage_name, stage_map in score_history.items():
        flat[stage_name] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(lig)  # normalize to basename for CSV
            s = rec.get("score", None)
            if rec.get("valid", False):
                flat[stage_name][lig_key] = s if s is not None else ""
            else:
                flat[stage_name][lig_key] = f"{s:.2f} (invalid)" if isinstance(s, (int, float)) else "(invalid)"

    write_score_summary_to_csv(flat, output_path=csv_out_wide)

    # -------- Long-form table with LE --------
    csv_out_long = os.path.join(dock_dir, "docking_score_long.csv")
    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["stage", "ligand", "score", "valid", "reason", "heavy_atoms", "le"])

        for stage_name, stage_map in score_history.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(lig)
                score = rec.get("score", None)
                valid = bool(rec.get("valid", False))
                reason = rec.get("reason", "")
                ha = rec.get("heavy_atoms", None)
                le = rec.get("le", None)

                # Normalize types for clean CSV
                score_str = f"{score:.2f}" if isinstance(score, (int, float)) else ""
                ha_str = str(int(ha)) if isinstance(ha, (int, float)) else ""
                le_str = f"{le:.4f}" if isinstance(le, (int, float)) else ""
                reason_str = str(reason) if reason is not None else ""

                writer.writerow([stage_name, lig_key, score_str, int(valid), reason_str, ha_str, le_str])

    return csv_out_wide

def compute_ligand_efficiency(score: Optional[float], heavy_atoms: Optional[int]) -> Optional[float]:
    """
    Compute ligand efficiency (LE) as -score / heavy_atoms.
    - score: docking score (kcal/mol), usually negative.
    - heavy_atoms: number of non-hydrogen atoms.
    Returns LE as a positive float (kcal/mol/HA) or None if not computable.
    """
    try:
        if score is None or heavy_atoms is None or int(heavy_atoms) <= 0:
            return None
        return float(-float(score) / int(heavy_atoms))
    except Exception:
        return None

#LE is ligand efficiency, it measures score/heavy atom count
#this basically tells you how efficient the ligand is at binding
#this prevent larger ligands from getting really high scores
#you dont want that because larger ligands are biased towards binding with more things
#hypothetically high LE means that the ligands is more specific and wont just bind
#to any old protein like a low LE might
def record_le(score_history: Dict[str, Dict[str, Dict]],
              stage_name: str,
              lig_path: str,
              score: Optional[float],
              heavy_atom_counts: Dict[str, int]) -> Optional[float]:
    """
    Store LE (and heavy atom count) into score_history for later CSV export or iteration logic.

    score_history structure (existing):
        { stage_name: { lig_path: { "score": float|None, "valid": bool, ... } } }

    After this call, you'll also have:
        { ..., "le": float|None, "heavy_atoms": int|None }

    Returns the computed LE (or None if not computable).
    """
    # look up heavy atoms by *normalized* path if that's how you keyed them
    ha = heavy_atom_counts.get(lig_path)
    le = compute_ligand_efficiency(score, ha)

    stage_map = score_history.setdefault(stage_name, {})
    rec = stage_map.setdefault(lig_path, {})
    rec["heavy_atoms"] = int(ha) if isinstance(ha, (int, float)) else None
    rec["le"] = le
    return le

# ======================
# Per-protein driver
# ======================

def process_one_protein(cfg: Dict, pdb_file: str, stages: List[Dict], params: RecenterParams) -> None:
    """
    Orchestrates all phases for a single protein.
    """
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = base_id.replace("_cleaned", "")
    logger = make_protein_logger(cfg["DOCKED_DIR"], pdb_id)
    logger.info(f"Processing protein: {pdb_file} (id={pdb_id})")

    paths = make_paths(cfg, base_id, pdb_file)

    # 1) Extract ligands → produce nolig PDB
    _lig_count, control_stems = extract_ligands_to_nolig(paths, logger)

    # 2) Protein prep (re-use if cached)
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("Skipping protein due to prep failure.")
        return

    # 3) Pocket detection
    center, box_size = detect_pocket(cleaned_pdb, logger)
    if center is None:
        return

    # 4) Ligand prep & filtering
    ligands, heavy_atom_counts = prepare_and_filter_ligands(cfg, paths, logger)
    if not ligands:
        logger.warning("No valid ligands after filtering; skipping protein.")
        return

    # 4.5) Initialize center selector
    selector = CenterSelector(cfg, logger, control_stems, heavy_atom_counts, center)

    # 5) Multi-stage docking w/ heuristics
    stage1_original = ligands[:]
    score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
    validated_ligands_last: List[str] = []
    recenter_attempts = 0
    docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()

    retry_mgr = RetryManager()  # <-- instantiate once per protein

    i = 0
    while i < len(stages):
        stage = stages[i]

        # Optional checkpoint skip (fingerprinted)
        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
            fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
            if checkpoint_should_skip(cfg, paths.pdb_id, stage["name"], fp):
                logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                i += 1
                continue

        if not ligands:
            logger.warning(f"No ligands to dock at {stage['name']}; stopping for this protein.")
            break

        logger.info(f"Starting {stage['name']} with {len(ligands)} ligands...")
        scores, validated, distances, raw_docked, invalids = run_one_stage(
            cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage, ligands, logger, retry_mgr
        )
        validated_ligands_last = validated

        # ---------- Stage invariant check (every ligand must be accounted for) ----------
        try:
            processed = set(ligands)
            valid_set = set(scores.keys())
            invalid_set = set(invalids.keys())
            both = valid_set & invalid_set
            missing = processed - (valid_set | invalid_set)
            if both or missing:
                logger.error(f"Invariant violation at {stage['name']}: both={len(both)}, missing={len(missing)}")
                if both:
                    logger.error("Ligands marked both valid & invalid: " + ", ".join(os.path.basename(x) for x in list(both)[:10]))
                if missing:
                    logger.error("Ligands missing from results: " + ", ".join(os.path.basename(x) for x in list(missing)[:10]))
        except Exception as _e:
            logger.warning(f"Invariant check failed: {_e}")
        # -------------------------------------------------------------------------------

        # Record both valid and invalid into history for CSV
        for lig, sc in scores.items():
            record_score(score_history, stage['name'], lig, sc, True)
            record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)
        for lig, (sc, reason) in invalids.items():
            record_score(score_history, stage['name'], lig, sc, False, reason=reason)
            record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)

        # Consider switching the center based on clusters/controls
        promoted_this_stage = False
        try:
            decision = selector.consider_switch(stage['name'], scores, validated, raw_docked, receptor_pdbqt, center)
            if decision.promoted and decision.new_center is not None:
                old = center
                center = decision.new_center
                promoted_this_stage = True
                if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                    checkpoint_invalidate_from(cfg, paths.pdb_id, stages, start_index=i)
                logger.info(
                    f"[CENTER] Switched from {old} -> {center} ({decision.reason}, SwitchScore={decision.switchscore:.2f})"
                )
        except Exception as e:
            logger.warning(f"CenterSelector failed gracefully: {e}")

        # Stage-1 early recenter / expand box — skip if we just promoted a center
        if not promoted_this_stage:
            restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                i, scores, distances, box_size, center, stage1_original, recenter_attempts, params,
                cfg, paths.pdb_id, receptor_pdbqt, logger, raw_docked
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

        # Selection or generic fallback (for all but last stage)
        if i < len(stages) - 1:
            if not scores:
                restart, center, box_size, redo_ligands = fallback_recentering_if_empty(
                    cfg, paths.pdb_id, stage['name'], scores, raw_docked,
                    receptor_pdbqt, center, box_size, stage1_original, logger
                )
                if restart:
                    ligands = redo_ligands
                    if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                        checkpoint_invalidate_from(cfg, paths.pdb_id, stages, start_index=0)
                    i = 0
                    continue
                else:
                    break

            ligands = select_ligands_for_next(docking_mode, i, stages, scores, logger)
            if not ligands:
                logger.warning(f"No ligands selected for {stages[i + 1]['name']}; stopping.")
                break

        # Mark checkpoint complete (fingerprinted)
        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
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
    cfg.setdefault("ALLOW_SWITCH_FROM_CONTROL", True)   # allow switching away if evidence is strong
    cfg.setdefault("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False)

    # clustering + switching thresholds
    cfg.setdefault("CLUSTER_EPS_ANG", 3.5)               # Å for pose-centroid clustering
    cfg.setdefault("SWITCH_VALID_RATE_MIN", 0.40)        # fraction of all docked ligands landing in candidate cluster
    cfg.setdefault("SWITCH_SCORE_IMPROVE_MIN", 1.5)      # kcal/mol better (median) vs current-center cluster
    cfg.setdefault("SWITCH_LE_GAIN_MIN", 0.02)           # kcal/mol per heavy atom gain (median)
    cfg.setdefault("SWITCH_SCORE_THRESHOLD", 0.70)       # final SwitchScore promote threshold
    cfg.setdefault("SWITCH_SCORE_HYSTERESIS", 0.50)      # reserved for future demotion logic

    # Critical defaults to avoid CPU oversubscription when MAX_PARALLEL_JOBS is large.
    cfg.setdefault("THREADS_PER_VINA", 1)
    cfg.setdefault("RMSD_FILTER_ANG", 2.0)  # Å; typical 1.5–2.0
    cfg.setdefault("RMSD_MAX_MODELS", 3)    # keep at most N distinct poses per ligand per stage

    # ---  safe defaults  ---
    cfg.setdefault("RETRY_NEAR_MISS", True)
    cfg.setdefault("RETRY_EXHAUST_MULT", 2)           # multiply stage exhaustiveness
    cfg.setdefault("RETRY_DIST_THRESH", 7.0)          # Å: close to pocket
    cfg.setdefault("RETRY_SCORE_THRESH", -7.5)        # kcal/mol: promising score

    cfg.setdefault("ADAPTIVE_SHRINK_ENABLE", True)
    cfg.setdefault("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0) # Å: cluster is tight
    cfg.setdefault("ADAPTIVE_SHRINK_DEC", 4.0)        # Å: shrink per dimension
    cfg.setdefault("ADAPTIVE_SHRINK_MIN_BOX", 14.0)   # Å: do not go below

    cfg.setdefault("CHECKPOINT_ENABLE", True)        # ON by default
    cfg.setdefault("PAINS_NAMES", "")                 # e.g. "EUG,HEM,CHS,PLP"
    cfg.setdefault("RECEPTOR_SANITY_CHECK", True)
    cfg.setdefault("AUDIT_JSON", True)

    params = get_recenter_params(cfg)
    stages = define_docking_stages(cfg.get("DOCKING_MODE", "discovery").lower())
    print("current docking mode is ", cfg.get("DOCKING_MODE"))

    # Discover input PDBs (skip *_nolig and non-PDB files)
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
