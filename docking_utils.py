from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from path_router import docked_dir, make_paths
from pose_validation import (
    attempt_fallback_recenter,
    extract_surface_atoms,
    filter_and_rewrite_poses_by_rmsd,
    validate_pose_pdbqt,
)
from run_vina import validate_all_poses
from fallback_recenter import GlobalCenterGuard, RecenterParams
from input_and_export_functions import extract_best_score, record_score, score_key


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

