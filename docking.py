# docking.py
# Docking orchestration helpers extracted from main.py

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from input_and_export_functions import emit_vina_config
from rdkit.Chem import rdMolAlign
from tqdm import tqdm

import automate_protein_prep as protein_prep
from apo_holo_mode import (
    _debug_normalize_mode_token,
    _record_apo_holo_decision,
    _record_apo_holo_usage,
    _variant_receptor_path,
    delete_variant_trees,
    file_sha1,
    resolve_apo_holo_mode,
)
from checkpoints import (
    checkpoint_invalidate_from,
    checkpoint_mark_done,
    checkpoint_should_skip,
)
from fallback_recenter import (
    BudgetGuard,
    GlobalCenterGuard,
    RecenterParams,
    fallback_recentering_if_empty,
    validate_first_valid_pose,
)
from input_and_export_functions import (
    load_inputs, validate_config, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, emit_vina_config as _emit_vina_config_impl, record_score, score_key, _to_bool, init_config_run_dir, 
     _to_bool, extract_best_score, record_score, score_key,
)
from logging_topics import make_protein_logger
# collapse_sanitized_names lives in main.py and is not used here.
from path_router import (
    RouterPaths,
    Paths,
    docked_dir,
    make_paths,
    receptor_file,
)
from single_ligand_index import (
    _ensure_single_ligand_index,
    _resolve_single_ligand,
)
from ph_ensemble_docking import (
    enumerate_ligands_for_ph_context,
    init_ph_tags_and_manifest,
    prewarm_ph_ligand_microstates,
)
from pose_validation import (
    attempt_fallback_recenter,
    compute_self_rmsd,
    extract_surface_atoms,
    filter_and_rewrite_poses_by_rmsd,
    validate_pose_pdbqt,
    _pose_centroid_from_pdbqt,
)
from prep_ligands import (
    enumerate_ligands_for_docking,
    prep_ligands_from_pdb,
)
from protein_functions import detect_active_site
from record_data import record_le, write_scores_csv
from run_vina import run_docking_task, validate_all_poses

import docking_controls
from docking_controls import (
    build_control_lookup,
    receptor_sanity_check,
    extract_ligands_to_nolig,
    _ph_values_from_context,
    _ph_ligand_mode,
    _resolve_ph_scope,
    select_center_via_control_redock,
    detect_pocket,
    _summarize_ions_file,
)
from docking_ligands import (
    select_ligands_for_next,
    _count_heavy_atoms_from_pdbqt,
    _coerce_test_map,
    _resolve_test_mode,
    _dedup_index_roots,
    _lib_roots_for_pdb,
    prepare_and_filter_ligands,
    _read_any_lig,
    _is_readable_ref,
    compute_rmsd,
    validate_ligand,
)
from docking_receptor import prepare_receptor



def norm(p: str | Path) -> str:
    """Normalize path to a clean, forward-slash string for logs & keys."""
    return os.path.abspath(str(p)).replace("\\", "/")


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

def _write_audit_json(cfg: Dict, pdb_id: str, summary: Dict, ph_label: Optional[str] = None,
                     variant: Optional[str] = None):
    try:
        if not cfg.get("AUDIT_JSON", True):
            return
        variant_env = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
        legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
        out = docked_dir(pdb_id, variant=variant_env, ph_tag=ph_label, legacy=legacy_mode) / "audit.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, indent=2))
    except Exception:
        pass

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


def _pose_path_for(csv_cfg: Dict, pdb_id: str, stage_name: str, lig_path: str,
                   ph_label: Optional[str] = None, variant: Optional[str] = None) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    from pathlib import Path
    # >>> DOCKED PATHS PATCH START
    paths = make_paths(csv_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_token = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    stage_dir = paths.docked_stage_dir(variant_token, stage_name, ph_token)
    return str(stage_dir / f"{Path(lig_path).stem}_{stage_name}.pdbqt")
    # >>> DOCKED PATHS PATCH END


docking_controls._is_readable_ref = _is_readable_ref
docking_controls.compute_rmsd = compute_rmsd

# Bridge functions for docking helpers
import docking as _docking

_docking.emit_vina_config = emit_vina_config
_docking.norm = norm
_docking._read_any_lig = _read_any_lig
_docking.validate_ligand = validate_ligand


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
    # ======================
    # Phase 1 – Variant & ion context
    # ======================
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    variant_label = variant_env or "legacy"
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    active_ph_label = None

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

                _, noncontrol_roots_single = _lib_roots_for_pdb(cfg, paths.pdb_id.upper(), paths, logger)
                ph_root_path = noncontrol_roots_single[0] if noncontrol_roots_single else None
                ph_root_cfg = str(ph_root_path) if ph_root_path else ""

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
                        "[single.ph_ligand.bridge.skip] no valid ph_ligand_root; "
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

    # Resolve the primary non-control library root for pH ligands
    _ctrl_roots_ph, noncontrol_roots_ph = _lib_roots_for_pdb(cfg, paths.pdb_id.upper(), paths, logger)
    ph_ligand_root = noncontrol_roots_ph[0] if noncontrol_roots_ph else None

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
                if checkpoint_should_skip(
                    cfg,
                    paths.pdb_id,
                    stage["name"],
                    fp,
                    ph_label=ph_label,
                    variant=variant_env or None,
                ):
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
                    ctrls, logger, retry_mgr, control_lookup, ph_label=ph_label
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
                    non_ctrls, logger, retry_mgr, control_lookup, ph_label=ph_label
                )

                scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
                raw_docked = {**rd1, **rd2}
                invalids = {**inv1, **inv2}
            else:
                scores, validated, distances, raw_docked, invalids = run_one_stage(
                    cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                    ligands, logger, retry_mgr, control_lookup, ph_label=ph_label
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
                        checkpoint_invalidate_from(
                            cfg,
                            paths.pdb_id,
                            stages,
                            start_index=i,
                            ph_label=ph_label,
                            variant=variant_env or None,
                        )
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
                        checkpoint_invalidate_from(
                            cfg,
                            paths.pdb_id,
                            stages,
                            start_index=0,
                            ph_label=ph_label,
                            variant=variant_env or None,
                        )
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
                            checkpoint_invalidate_from(
                                cfg,
                                paths.pdb_id,
                                stages,
                                start_index=0,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
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
                    checkpoint_mark_done(
                        cfg,
                        paths.pdb_id,
                        stage["name"],
                        fp,
                        ph_label=ph_label,
                        variant=variant_env or None,
                    )
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

        csv_path = write_scores_csv(
            cfg,
            paths.pdb_id,
            score_history,
            ph_label=ph_label,
            variant=variant_env or None,
        )
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
            _write_audit_json(cfg, paths.pdb_id, summary, ph_label=ph_label, variant=variant_env or None)
        except Exception as _e:
            logger.warning(f"Audit JSON write failed: {_e}")

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
    ph_label: Optional[str] = None,
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
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

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

                stage_name = stage["name"]
                logger.info(
                    "[emit.debug] pdb=%s stage=%s variant=%s ph=%s",
                    pdb_id,
                    stage_name,
                    variant_token or "None",
                    ph_label or "None",
                )
                conf_path, out_path = emit_vina_config(
                    cfg,
                    pdb_id,
                    receptor_pdbqt,
                    center,
                    box_size,
                    lig,
                    stage_name,
                    stage_for_cfg,
                    threads_per_vina,
                    logger,
                    variant=variant_token,
                    ph_token=ph_label,
                    legacy=legacy_mode,
                )

                # anchor: emit_vina_config resolves variant/pH from cfg/env
                logger.info("[cfg.emit] %s -> %s", os.path.basename(lig), conf_path)

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
            with tqdm(
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
                            # --- compute self-RMSD for near-misses (log-only) ---
                            if bool(cfg.get("LOG_SELF_RMSD", True)):
                                try:
                                    self_rmsd_val = compute_self_rmsd(out_path)
                                except Exception as _e:
                                    self_rmsd_val = None
                                    logger.warning(f"self-RMSD failed for {lig_name}: {_e}")
                                try:
                                    logger.info(f"[self-rmsd] lig={lig_name} rmsd={self_rmsd_val}")
                                except Exception:
                                    pass

                            reason = result.get("reason", "") or ""
                            near_miss = (
                                ("clash" in reason) or
                                (result.get("distance_to_surface") or 0.0) < 6.5 or
                                (result.get("distance_to_centroid") or 0.0) < 4.0
                            )

                            if near_miss:
                                stage_retry = dict(stage)
                                stage_retry["name"] = f"{stage['name']}_retry"
                                stage_retry["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
                                if cfg.get("FAST_MODE"):
                                    stage_retry["exhaustiveness"] = 1

                                try:
                                    base_seed = int(stage_retry.get("seed", 0)) if "seed" in stage_retry else 0
                                except Exception:
                                    base_seed = 0
                                stage_retry["seed"] = base_seed + 137  # small jitter
                                stage_retry["num_modes"] = int(cfg.get("NEAR_MISS_NUM_MODES", stage_retry.get("num_modes", 4)))
                                stage_retry["energy_range"] = float(cfg.get("NEAR_MISS_ENERGY_RANGE", stage_retry.get("energy_range", 4.0)))

                                if bool(cfg.get("NEAR_MISS_RECENTER", True)):
                                    # If near miss, try to recenter on current pose centroid
                                    try:
                                        cent = _pose_centroid_from_pdbqt(str(out_path))
                                    except Exception as _e:
                                        logger.warning(f"[near-miss] failed to compute centroid for {lig_name}: {_e}")
                                        cent = None

                                    if cent and isinstance(cent, (list, tuple)) and len(cent) == 3:
                                        try:
                                            center_nm = tuple(float(x) for x in cent)
                                            logger.info(f"[near-miss] recentering on best pose centroid {center_nm}")
                                            try:
                                                max_box = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                                            except Exception:
                                                max_box = 28.0
                                            box_nm = tuple(min(max_box, s + 1.0) for s in box_size)
                                        except Exception:
                                            center_nm = center
                                            box_nm = box_size
                                    else:
                                        center_nm = center
                                        box_nm = box_size
                                else:
                                    center_nm = center
                                    box_nm = box_size

                                # anchor: retry config repeats variant/pH tagging
                                conf_path2, out_path2 = emit_vina_config(
                                    cfg,
                                    pdb_id,
                                    receptor_pdbqt,
                                    center_nm,
                                    box_nm,
                                    lig,
                                    stage_retry["name"],
                                    stage_retry,
                                    threads_per_vina,
                                    logger,
                                    variant=variant_token,
                                    ph_token=ph_label,
                                    legacy=legacy_mode,
                                )
                                logger.info("[cfg.emit] %s -> %s", os.path.basename(lig), conf_path2)
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

                        stage_retry2_name = stage_retry2["name"]
                        logger.info(
                            "[emit.debug] pdb=%s stage=%s variant=%s ph=%s",
                            pdb_id,
                            stage_retry2_name,
                            variant_token or "None",
                            ph_label or "None",
                        )
                        conf_path3, out_path3 = emit_vina_config(
                            cfg,
                            pdb_id,
                            receptor_pdbqt,
                            retry_center,
                            retry_box,
                            lig,
                            stage_retry2_name,
                            stage_retry2,
                            threads_per_vina,
                            logger,
                            variant=variant_token,
                            ph_token=ph_label,
                            legacy=legacy_mode,
                        )
                        # anchor: fallback config mirrors initial variant/pH discovery
                        logger.info("[cfg.emit] %s -> %s", os.path.basename(lig), conf_path3)
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


__all__ = [
    "_map_reason_to_category",
    "RetryManager",
    "run_one_stage",
    "process_one_protein",
    "_coerce_test_map",
    "_resolve_test_mode",
]
