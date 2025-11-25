from __future__ import annotations

import logging
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from tqdm import tqdm
from rdkit.Chem import rdMolAlign

from fallback_recenter import validate_first_valid_pose, BudgetGuard
from pose_validation import extract_surface_atoms, filter_and_rewrite_poses_by_rmsd, attempt_fallback_recenter, compute_self_rmsd
from path_router import make_paths
from run_vina import run_docking_task

"""
Single-stage docking engine used by the Atlas pipeline.

This module provides run_one_stage(...), which:
- runs a multi-ligand dock for a given 'stage' configuration,
- performs per-ligand pose validation and structured retries,
- applies RMSD-based filters and pocket-distance checks,
- returns scores, validated ligands, distances, raw pose paths, and invalid reasons.

High-level orchestration (center selection, checkpoints, stage-to-stage transitions)
is intentionally kept in main.py or higher-level orchestration modules.
"""


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

    from main import emit_vina_config, _read_any_lig

    threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))
    max_workers = int(cfg["MAX_PARALLEL_JOBS"])
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None
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

    def _map_reason_to_category(reason: str) -> str:
        r = (reason or "").lower()
        if "pocket" in r:
            return "too_far_from_pocket"
        if "pose" in r:
            return "no_valid_pose"
        if "timeout" in r:
            return "timeout"
        return "malformed"

    def _validate_with_rmsd_gate(lig_path: str, lig_name: str, out_path: str, score: float) -> Tuple[bool, str]:
        # Hard gate controls if experimental heavy-atom RMSD is too high.
        lig_key = Path(lig_path).stem
        control_ref = control_lookup.get(lig_key)
        if control_ref and Path(control_ref).exists():
            # Try PDBQT (preferred) or SDF/MOL2 fallback. Any valid path is ok.
            ref_path = control_ref
            if os.path.splitext(str(ref_path))[-1].lower() != ".pdbqt":
                ref_path2 = Path(ref_path).with_suffix(".pdbqt")
                if ref_path2.exists():
                    ref_path = ref_path2
            if not os.path.exists(ref_path):
                logger.warning(f"{lig_name}: reference control not found: {ref_path}")
            else:
                rmsd = compute_rmsd(str(ref_path), out_path)
                if not np.isfinite(rmsd) or rmsd > float(cfg.get("MAX_CONTROL_RMSD", 8.0)):
                    logger.warning(f"{lig_name} control RMSD {rmsd} exceeds max; marking invalid")
                    return False, "rmsd_gate_fail"
                else:
                    logger.info(f"{lig_name} control RMSD {rmsd}")

        # Ligand-only self RMSD gate (supports per-stage override of the threshold)
        self_rmsd_threshold = float(stage.get("SELF_RMSD_THRESH", cfg.get("SELF_RMSD_THRESHOLD", 0)))
        if self_rmsd_threshold > 0:
            try:
                sr = compute_self_rmsd(out_path)
            except Exception:
                sr = None
            if sr is None or (not np.isfinite(sr)) or (sr > self_rmsd_threshold):
                return False, "self_rmsd_gate_fail"

        return True, ""

    # ----- Stage-level emit + docking -----
    # Key: (future -> ligand)
    futures: Dict = {}

    # Pool within the stage
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        ligs = ligands
        ligs_with_name = [(lig, os.path.basename(lig)) for lig in ligs]
        ligs_sorted = sorted(ligs_with_name, key=lambda kv: kv[1])

        # The stage name is used to build output paths; ensure no spaces.
        stage_name = stage["name"].replace(" ", "_")
        ws_path = Path(cfg.get("WORKSPACE_DIR", cfg.get("OUTPUT_DIR", ".")))

        # For tracking and dedup of retries with jittered configs
        vinaconfs: Dict[str, str] = {}  # lig -> config path
        vinapaths: Dict[str, str] = {}

        logger.info(
            "[dock.stage] stage=%s ligands=%d threads_per_job=%d max_workers=%d",
            stage_name,
            len(ligs),
            threads_per_vina,
            max_workers,
        )

        attempt = 0
        while attempt < retry_mgr.max_retries:
            futures = {}
            processed = 0
            with tqdm(ligs_sorted, desc=f"Docking {stage_name}", file=_stdout, mininterval=5) as pbar:
                for lig, lig_name in ligs_sorted:
                    # Workdir within the workspace
                    ligand_ws = ws_path / "docking" / Path(lig_name).stem
                    ligand_ws.mkdir(parents=True, exist_ok=True)

                    guard = None
                    if budget_guards is not None:
                        guard = budget_guards.get(lig)
                    if guard is None:
                        guard = BudgetGuard(
                            dock_max_seconds=float(cfg.get("POSE_MAX_SEC", 240)),
                            idle_max_seconds=float(cfg.get("POSE_IDLE_MAX_SEC", 120)),
                            log=logger,
                        )
                    guard.reset()

                    lig_n = norm(lig)
                    out_n = norm(ligand_ws / f"{Path(lig_name).stem}_{stage_name}.pdbqt")

                    stage_params = dict(stage)
                    stage_params["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
                    if cfg.get("FAST_MODE"):
                        try:
                            ex0 = int(stage_params.get("exhaustiveness", 8))
                        except Exception:
                            ex0 = 8
                        stage_params["exhaustiveness"] = 1
                        stage_params["num_modes"] = min(3, stage_params.get("num_modes", 3))

                    try:
                        ex0 = int(stage_params.get("exhaustiveness", 8))
                    except Exception:
                        ex0 = 8
                    stage_params["exhaustiveness"] = int(cfg.get("STAGE_EXHAUST_OVERRIDE", ex0))

                    variant_specific_gates = cfg.get("_VARIANT_STAGE_GATES") or {}
                    variant_key = variant_token or "ALL"
                    if variant_key in variant_specific_gates:
                        try:
                            gates_for_variant = variant_specific_gates[variant_key]
                            stage_params.update(gates_for_variant.get(stage_name, {}))
                        except Exception:
                            pass

                    lig_wc = os.path.basename(lig)
                    if lig_wc in control_lookup:
                        # Use a relaxed energy range for control re-docks.
                        stage_params["energy_range"] = max(stage_params.get("energy_range", 3.0), 5.0)

                    lig_cfg, out_pdbqt = emit_vina_config(
                        cfg,
                        pdb_id,
                        receptor_pdbqt,
                        center,
                        box_size,
                        lig,
                        stage_name,
                        stage_params,
                        threads_per_vina,
                        logger,
                        variant=variant_token,
                        ph_token=ph_label,
                        legacy=legacy_mode,
                    )
                    vinaconfs[lig] = lig_cfg
                    vinapaths[lig] = out_pdbqt

                    guard.bind(lig_cfg, out_pdbqt)

                    futures[pool.submit(run_docking_task, cfg["VINA_EXE"], lig_cfg, lig, out_pdbqt)] = (lig_n, out_n)

                for fut in as_completed(futures):
                    lig, out_path = futures[fut]
                    lig_name = os.path.basename(lig)

                    guard = None
                    if budget_guards is not None:
                        guard = budget_guards.get(lig)
                    if guard is None:
                        guard = BudgetGuard(
                            dock_max_seconds=float(cfg.get("POSE_MAX_SEC", 240)),
                            idle_max_seconds=float(cfg.get("POSE_IDLE_MAX_SEC", 120)),
                            log=logger,
                        )

                    score = None
                    try:
                        score = fut.result()
                    except Exception as e:
                        if guard.expired():
                            invalids[lig] = (None, "budget_exceeded")
                        else:
                            invalids[lig] = (None, f"vina_failed: {e}")

                        logger.warning(f"{lig_name} | Vina failed: {e}")
                        processed += 1
                        if (processed % 25 == 0) or (processed == len(futures)):
                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                        pbar.update(1)
                        continue

                    out_path = out_path or vinapaths.get(lig)
                    lig_cfg = vinaconfs.get(lig)

                    if out_path is None or (not os.path.exists(out_path)) or (os.path.getsize(out_path) < 100):
                        invalids[lig] = (None, "vina_missing_output")
                        logger.warning(f"{lig_name} | Vina output missing or too small: {out_path}")
                        processed += 1
                        if (processed % 25 == 0) or (processed == len(futures)):
                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                        pbar.update(1)
                        continue

                    # Optional quick filter + RMSD rewrite against the best pose to reduce duplicates
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
                                stage_retry_name = stage_retry["name"]
                                logger.info(
                                    "[emit.debug] pdb=%s stage=%s variant=%s ph=%s",
                                    pdb_id,
                                    stage_retry_name,
                                    variant_token or "None",
                                    ph_label or "None",
                                )
                                conf_path2, out_path2 = emit_vina_config(
                                    cfg,
                                    pdb_id,
                                    receptor_pdbqt,
                                    retry_center,
                                    retry_box,
                                    lig,
                                    stage_retry_name,
                                    stage_retry,
                                    threads_per_vina,
                                    logger,
                                    variant=variant_token,
                                    ph_token=ph_label,
                                    legacy=legacy_mode,
                                )
                                # anchor: retry config uses same variant/pH resolution
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

                                result = validate_first_valid_pose(
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

                                logger.info(f"{lig_name} retry -> {result}")

                                if result.get("valid", False) and score2 is not None:
                                    ok2, reason2 = _validate_with_rmsd_gate(lig, lig_name, out_path2, float(score2))
                                    if ok2:
                                        scores[lig] = float(score2)
                                        validated_ligands.append(lig)
                                        raw_docked_ligands[lig] = norm(out_path2)
                                        logger.info(f"{lig_name} retry score: {score2:.2f} kcal/mol")
                                        processed += 1
                                        if (processed % 25 == 0) or (processed == len(futures)):
                                            pbar.set_postfix(ok=len(scores), inv=len(invalids))
                                        pbar.update(1)
                                        continue
                                    else:
                                        invalids[lig] = (float(score2), reason2 or "rmsd_fail")

                        except Exception:
                            pass

                    if did_retry:
                        continue

                    try:
                        fb_pose, new_c, _bs, _ch = attempt_fallback_recenter(
                            receptor_pdbqt,
                            center,
                            box_size,
                            out_path,
                            pocket_center=paths.pocket_center,
                            box_size_hint=paths.box_size,
                            logger=logger,
                            exclude_basenames=set(),
                        )
                        if new_c is not None:
                            retry_center = new_c
                        else:
                            retry_center = center
                        if fb_pose is not None:
                            retry_box = _bs
                            fallback_ligands = {lig: fb_pose}
                        else:
                            retry_box = box_size
                            fallback_ligands = {lig: out_path}
                    except Exception as e:
                        logger.warning(f"Fallback recentering failed for {lig_name}: {e}")
                        fallback_ligands = {lig: out_path}
                        retry_center = center
                        retry_box = box_size

                    # Hard-coded retry strategy (one per ligand for now)
                    retained_invalid = True
                    while attempt < retry_mgr.max_retries:
                        err_cat = _map_reason_to_category(result.get("reason", ""))
                        recipe = retry_mgr.apply(stage, err_cat, attempt)
                        if recipe is None:
                            attempt += 1
                            continue

                        stage_retry2 = dict(stage)
                        stage_retry2.update(recipe)
                        stage_retry2["name"] = f"{stage['name']}_retry{attempt + 1}"
                        stage_retry2.setdefault("verbosity", int(cfg.get("VINA_VERBOSITY", 0)))

                        try:
                            if recipe.get("recenter"):
                                _, new_c, _bs2, _ch = attempt_fallback_recenter(
                                    receptor_pdbqt,
                                    center,
                                    box_size,
                                    fallback_ligands,
                                    pocket_center=paths.pocket_center,
                                    box_size_hint=paths.box_size,
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
