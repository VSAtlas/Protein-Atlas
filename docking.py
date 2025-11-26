# docking.py
# Docking orchestration helpers extracted from main.py

from __future__ import annotations

import hashlib
import logging
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rdkit.Chem import rdMolAlign
from tqdm import tqdm

# Copy the imports that run_one_stage relies on.
# Easiest is to copy the relevant subset from main.py:
from fallback_recenter import (
    BudgetGuard,
    RecenterParams,
    GlobalCenterGuard,
    validate_first_valid_pose,
    fallback_recentering_if_empty,
)
from checkpoints import (
    checkpoint_invalidate_from,
    checkpoint_mark_done,
    checkpoint_should_skip,
)
from input_and_export_functions import extract_best_score, record_score, score_key
from path_router import RouterPaths, make_paths
from pose_validation import (
    validate_pose_pdbqt,
    extract_surface_atoms,
    attempt_fallback_recenter,
    filter_and_rewrite_poses_by_rmsd,
    compute_self_rmsd,
)
from run_vina import run_docking_task, validate_all_poses


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
                                    # If near miss, try to recenter on best-scoring pose centroid
                                    res = validate_pose_pdbqt(
                                        receptor_pdbqt=receptor_pdbqt,
                                        ligand_pdbqt=out_path,
                                        pocket_center=center,
                                        surface_coords=surface_coords,
                                        clash_threshold=2.0,
                                        clash_tol=3,
                                        dist_surf=6.0,
                                        dist_centroid=4.5,
                                        max_models=int(cfg.get("EARLY_EXIT_MAX_MODELS", 3)),
                                        best_pose_centroid=True,
                                        early_exit=False,
                                    )
                                    cent = res.get("best_pose_centroid")
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
]
