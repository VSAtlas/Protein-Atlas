from __future__ import annotations

import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from tqdm import tqdm
from rdkit.Chem import rdMolAlign

from .docking_vina import emit_vina_config
from .docking_ligands import _read_any_lig, validate_ligand
from .docking_utils import norm
from .fallback_recenter import (
    BudgetGuard,
    attempt_fallback_recenter,
    validate_first_valid_pose,
)
from .pose_validation import (
    compute_self_rmsd,
    extract_surface_atoms,
    filter_and_rewrite_poses_by_rmsd,
)
from run_manifest import update_manifest_for_docking_stage
from .run_vina import run_docking_task


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
    recipes: Dict[str, List[Dict[str, Any]]] = field(
        default_factory=lambda: {
            # Near-miss rescue: small geometry tweak with a jittered seed.
            "near_miss": [
                {"recenter": True, "box_pad_delta": +1.0, "seed_jitter": True},
            ],
            # If we docked far from the pocket, try small geometry tweaks   not more modes
            "too_far_from_pocket": [
                {"recenter": True, "box_pad_delta": +1.0, "num_modes": 4},
                {
                    "recenter": True,
                    "box_pad_delta": +2.0,
                    "exhaustiveness": 6,
                    "num_modes": 4,
                },
            ],
            # If we didn't get a valid pose, explore new seeds and a slightly wider energy window,
            # but keep returned modes low so validation stays fast.
            "no_valid_pose": [
                {
                    "exhaustiveness": 6,
                    "num_modes": 5,
                    "seed_jitter": True,
                    "energy_range": 6,
                },
                {
                    "recenter": True,
                    "box_pad_delta": +1.0,
                    "exhaustiveness": 6,
                    "num_modes": 5,
                    "seed_jitter": True,
                    "energy_range": 6,
                },
            ],
            # If we timed out, go cheaper, not deeper.
            "timeout": [
                {"exhaustiveness": 3, "num_modes": 3, "seed_jitter": True},
                {"exhaustiveness": 2, "num_modes": 2},
            ],
            "malformed": [],  # do not retry
        }
    )

    def apply(
        self, base_params: Dict[str, Any], err_type: str, attempt: int
    ) -> Optional[Dict[str, Any]]:
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
    control_lookup: Dict[str, Path],  # maps ligand basename -> crystal ref PDB
    budget_guards: Optional[
        Dict[str, BudgetGuard]
    ] = None,  # external per-ligand guards
    ph_label: Optional[str] = None,
) -> Tuple[
    Dict[str, float],
    List[str],
    List[float],
    Dict[str, str],
    Dict[str, Tuple[Optional[float], str]],
]:
    import shutil
    import subprocess

    threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))
    cpu = int(cfg.get("CPU", os.cpu_count() or 1))
    max_jobs = int(cfg.get("MAX_PARALLEL_JOBS", cpu))
    max_workers = min(max_jobs, len(ligands))
    if max_workers < 1:
        max_workers = 1
    if max_jobs == 1 and len(ligands) > 0:
        logger.info(
            "[dock.parallel.single_ligand] MAX_PARALLEL_JOBS=1 forcing max_workers=1 for %d ligands",
            len(ligands),
        )
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    # paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # Note: make_paths is not available here. However, checking run_one_stage body:
    # It uses paths.docked_pdb_root() in attempt_fallback_recenter.
    # It constructs paths via make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb") inside the function?
    # Let me check the original code I read.
    # Yes, line 1888: paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    # I need to import make_paths from path_router.
    from path_router.path_router import make_paths

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    run_id = str(cfg.get("RUN_ID", "") or "")
    variant_label = variant_env or "legacy"
    raw_stage_name = (
        stage.get("name")
        or stage.get("stage_name")
        or stage.get("label")
        or stage.get("id")
        or "stage"
    )
    stage_name = str(raw_stage_name).strip() or "stage"
    stage_start = time.time()

    if cfg.get("FAST_MODE") and retry_mgr.max_retries > 1:
        retry_mgr.max_retries = 1
        logger.info("[retry] mode=fast stage=%s max_retries=1", stage_name)

    if run_id:
        try:
            update_manifest_for_docking_stage(
                cfg,
                run_id,
                pdb_id,
                variant_label,
                stage_name,
                status="running",
                ph_tag=ph_label,
            )
        except Exception:
            logger.warning(
                "[run-manifest.docking-stage] failed to record start pdb=%s variant=%s ph=%s stage=%s",
                pdb_id,
                variant_label,
                ph_label if ph_label is not None else "base",
                stage_name,
                exc_info=True,
            )

    scores: Dict[str, float] = {}
    validated_ligands: List[str] = []
    all_distances: List[float] = []
    raw_docked_ligands: Dict[str, str] = {}
    invalids: Dict[str, Tuple[Optional[float], str]] = {}

    surface_coords = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)
    processed = 0

    # ---- helpers (nested) ----
    def _best_pose_pdb_from_pdbqt(
        pdbqt_path: str, obabel_path: Optional[str] = None
    ) -> Optional[str]:
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
            cmd = [
                obabel,
                "-ipdbqt",
                str(pdbqt_path),
                "-opdb",
                "-O",
                str(out_pdb),
                "-f",
                "1",
                "-l",
                "1",
                "-d",
            ]
            subprocess.check_call(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            return (
                str(out_pdb)
                if out_pdb.exists() and out_pdb.stat().st_size > 0
                else None
            )
        except Exception as e:
            logger.warning(
                f"[RMSD] OpenBabel conversion failed for {os.path.basename(pdbqt_path)}: {e}"
            )
            return None

    # NOTE:
    # This nested compute_rmsd intentionally shadows the generic compute_rmsd imported from docking_ligands.
    # It adds extra logging and redock-specific behavior for controls in run_one_stage without changing the module-level helper.
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

            if (
                os.path.exists(ref_path)
                and os.path.exists(docked_path)
                and os.path.samefile(ref_path, docked_path)
            ):
                if _rlog:
                    _rlog.warning(
                        f"[rmsd.core] ref and dock resolve to the SAME file "
                        f"(ref='{ref_path}', dock='{docked_path}')"
                    )
        except Exception:
            pass
        # ───────────────────────────────────────────────────────────────────────────

        if not ref or not dock:
            if _rlog:
                _rlog.warning(
                    f"[rmsd.core] load-fail ref_ok={bool(ref)} dock_ok={bool(dock)} "
                    f"ref='{ref_path}' dock='{docked_path}'"
                )
            return float("inf")

        try:
            n_ref = ref.GetNumAtoms()
            n_dock = dock.GetNumAtoms()
        except Exception:
            n_ref = n_dock = -1

        has_conf_ref = ref.GetNumConformers() > 0
        has_conf_dock = dock.GetNumConformers() > 0

        # ── heavy-atom counts (useful when you get inf) ──────────────────
        try:
            ha_ref = ref.GetNumHeavyAtoms()
            ha_dock = dock.GetNumHeavyAtoms()
        except Exception:
            ha_ref = ha_dock = -1
        if _rlog:
            _rlog.info(
                f"[rmsd.core] inputs ref='{ref_path}' dock='{docked_path}' "
                f"n_ref={n_ref} n_dock={n_dock} heavy_ref={ha_ref} heavy_dock={ha_dock} "
                f"conf_ref={has_conf_ref} conf_dock={has_conf_dock}"
            )
        # ───────────────────────────────────────────────────────────────────────────

        if not has_conf_ref or not has_conf_dock:
            if _rlog:
                _rlog.warning("[rmsd.core] missing 3D conformers; returning inf")
            return float("inf")

        if n_ref != n_dock:
            if _rlog:
                _rlog.info(
                    f"[rmsd.core] atom_count_mismatch ({n_ref} vs {n_dock}); "
                    f"bestRMS will not be used; returning inf (no MCS fallback)"
                )
            return float("inf")

        try:
            val = float(rdMolAlign.GetBestRMS(ref, dock))
            if _rlog:
                _rlog.info(f"[rmsd.core] method=bestRMS rmsd={val:.3f}")
            return val
        except Exception as e:
            if _rlog:
                _rlog.info(
                    f"[rmsd.core] method=bestRMS failed: {e}; returning inf (no MCS fallback)"
                )
            return float("inf")

    def _validate_with_rmsd_gate(
        lig_path: str, lig_name: str, out_pdbqt_path: str, score_val: float
    ) -> Tuple[bool, Optional[str]]:
        """
        Controls: crystal redock RMSD is a hard gate.
        Non-controls: self-RMSD is logged upstream; do not gate here.
        """
        base = Path(lig_path).stem.split("_stage")[0]
        crystal_ref = control_lookup.get(base)

        if crystal_ref:
            best_pdb = _best_pose_pdb_from_pdbqt(
                out_pdbqt_path, obabel_path=cfg.get("OPENBABEL_PATH")
            )
            if not best_pdb:
                logger.warning(
                    f"{lig_name} | unable to extract best pose PDB for redock RMSD."
                )
                return False, "no_best_pose_for_rmsd"

            # --- AUDIT: control redock (compute RMSD just for logging) ---
            rmsd_val = compute_rmsd(str(crystal_ref), best_pdb)

            ok = validate_ligand(
                ligand_name=lig_name,
                docked_path=best_pdb,
                crystal_path=str(crystal_ref),
                rmsd_thresh=float(cfg.get("CONTROL_RMSD_MAX_ANG", 2.0)),
                self_rmsd=None,
                logger=logger,
            )

            logger.info(
                "[control-redock] lig=%s rmsd=%.2f A score=%.2f",
                lig_name,
                (rmsd_val if rmsd_val is not None else float("nan")),
                float(score_val),
            )

            if ok:
                logger.info(
                    f"{lig_name} | {stage['name']} score: {score_val:.2f} kcal/mol (redock-RMSD PASS)"
                )
                return True, None
            else:
                return False, "rmsd_fail"

        # Non-controls: redock gate not applicable here (geometry checks already passed).
        return True, None

    def _per_ligand_postprocess(
        lig: str,
        lig_name: str,
        out_path: str,
        score: float,
        guard: BudgetGuard,
        pbar: "tqdm",
    ) -> None:
        nonlocal \
            scores, \
            validated_ligands, \
            all_distances, \
            raw_docked_ligands, \
            invalids, \
            processed

        try:
            kept, removed = filter_and_rewrite_poses_by_rmsd(
                out_path,
                rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                max_models=int(cfg.get("RMSD_MAX_MODELS", 3)),
            )
            if removed > 0:
                logger.info(
                    f"{os.path.basename(out_path)}: RMSD filter kept {kept}, removed {removed}"
                )
        except Exception as e:
            logger.warning(
                f"RMSD filtering failed for {os.path.basename(out_path)}: {e}"
            )

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
            return

        if guard.expired():
            invalids[lig] = (float(score), "budget_exceeded")
            processed += 1
            if (processed % 25 == 0) or (processed == len(futures)):
                pbar.set_postfix(ok=len(scores), inv=len(invalids))
            pbar.update(1)
            return

        if _handle_structured_recipe_retries(
            lig=lig,
            lig_name=lig_name,
            out_path=out_path,
            score=score,
            result=result,
            guard=guard,
            pbar=pbar,
        ):
            return

        invalids[lig] = (float(score), result.get("reason", "pose_invalid"))
        logger.info(f"{lig_name} | pose invalid (after retries)")
        processed += 1
        if (processed % 25 == 0) or (processed == len(futures)):
            pbar.set_postfix(ok=len(scores), inv=len(invalids))
        pbar.update(1)

    def _handle_structured_recipe_retries(
        lig: str,
        lig_name: str,
        out_path: str,
        score: float,
        result: Dict[str, Any],
        guard: BudgetGuard,
        pbar: "tqdm",
    ) -> bool:
        nonlocal processed

        def _is_near_miss_failure(res: Dict[str, Any]) -> bool:
            if not bool(cfg.get("RETRY_NEAR_MISS", True)):
                return False
            reason_val = res.get("reason", "") or ""
            reason_lower = str(reason_val).lower()
            try:
                dist_surface = float(res.get("distance_to_surface", float("inf")))
            except Exception:
                dist_surface = float("inf")
            try:
                dist_centroid = float(res.get("distance_to_centroid", float("inf")))
            except Exception:
                dist_centroid = float("inf")
            return (
                ("clash" in reason_lower) or dist_surface < 6.5 or dist_centroid < 4.0
            )

        near_miss_hit = _is_near_miss_failure(result)
        if near_miss_hit and bool(cfg.get("LOG_SELF_RMSD", True)):
            try:
                self_rmsd_val = compute_self_rmsd(out_path)
                logger.info("[self-rmsd] lig=%s rmsd=%s", lig_name, self_rmsd_val)
            except Exception as _e:
                logger.warning(f"self-RMSD failed for {lig_name}: {_e}")

        err_cat = (
            "near_miss"
            if near_miss_hit
            else _map_reason_to_category(result.get("reason", ""))
        )
        attempt = 0
        retained_invalid = True

        while attempt < retry_mgr.max_retries:
            if guard.expired():
                invalids[lig] = (float(score), "budget_exceeded")
                retained_invalid = False
                break

            retry_params = retry_mgr.apply(stage, err_cat, attempt)
            if not retry_params:
                break

            stage_retry2 = dict(retry_params)
            stage_retry2["name"] = (
                f"{stage['name']}_retry"
                if attempt == 0
                else f"{stage['name']}_retry{attempt + 1}"
            )
            stage_retry2["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))

            if err_cat == "near_miss":
                if not bool(cfg.get("NEAR_MISS_RECENTER", True)):
                    stage_retry2.pop("recenter", None)
                stage_retry2["num_modes"] = int(
                    cfg.get(
                        "NEAR_MISS_NUM_MODES",
                        stage_retry2.get("num_modes", stage.get("num_modes", 4)),
                    )
                )
                stage_retry2["energy_range"] = float(
                    cfg.get(
                        "NEAR_MISS_ENERGY_RANGE",
                        stage_retry2.get(
                            "energy_range", stage.get("energy_range", 4.0)
                        ),
                    )
                )

            if cfg.get("FAST_MODE"):
                stage_retry2["exhaustiveness"] = 1
                stage_retry2["num_modes"] = 1

            seed_jitter = bool(stage_retry2.pop("seed_jitter", False))
            if seed_jitter:
                try:
                    base_seed = (
                        int(stage_retry2.get("seed", 0))
                        if "seed" in stage_retry2
                        else 0
                    )
                except Exception:
                    base_seed = 0
                stage_retry2["seed"] = base_seed + (attempt + 1) * 137

            retry_center = center
            retry_box = box_size
            try:
                if stage_retry2.get("recenter", False):
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
                box_pad = stage_retry2.get("box_pad_delta", None)
                if isinstance(box_pad, (int, float)):
                    dx = float(box_pad)
                    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                    retry_box = tuple(min(box_cap, s + dx) for s in box_size)
            except Exception as _e:
                logger.warning(f"Retry recenter/box tweak failed: {_e}")

            stage_retry2.pop("recenter", None)
            stage_retry2.pop("box_pad_delta", None)

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
            logger.info("[cfg.emit] %s -> %s", os.path.basename(lig), conf_path3)
            try:
                Path(conf_path3).resolve().relative_to(
                    Path(cfg["CONFIG_RUN_DIR"]).resolve()
                )
            except Exception:
                raise RuntimeError(
                    f"Refusing to launch Vina with config outside current RUN_DIR: {conf_path3}"
                )
            logger.info(f"[vina.call] config={conf_path3}")

            try:
                _, score_r = run_docking_task(
                    cfg["VINA_EXE"],
                    conf_path3,
                    lig,
                    out_path3,
                    write_failure_marker_flag=True,
                )
            except Exception as _e:
                logger.warning(f"Retry docking crashed for {lig_name}: {_e}")
                attempt += 1
                continue

            try:
                filter_and_rewrite_poses_by_rmsd(
                    out_path3,
                    rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                    max_models=int(cfg.get("RMSD_MAX_MODELS", 3)),
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

            logger.info(
                f"{lig_name} retry#{attempt + 1} ({stage_retry2_name}|{err_cat}) -> {result_r}"
            )
            if result_r.get("valid", False) and score_r is not None:
                ok_r, reason_r = _validate_with_rmsd_gate(
                    lig, lig_name, out_path3, float(score_r)
                )
                if ok_r:
                    scores[lig] = float(score_r)
                    validated_ligands.append(lig)
                    raw_docked_ligands[lig] = norm(out_path3)
                    logger.info(
                        f"{lig_name} | {stage_retry2['name']} score: {score_r:.2f} kcal/mol (retry rescued)"
                    )
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
        return True

    # ---- scheduling & submission ----
    default_budget_seconds = float(
        cfg.get("MAX_RETRY_SECONDS_PER_LIGAND", cfg.get("BENCH_MAX_SECONDS", 300.0))
    )

    submit_queue = []
    guards_for_ligand: Dict[str, BudgetGuard] = {}

    for lig in ligands:
        guard = budget_guards.get(lig) if budget_guards else None
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
                    stage_for_cfg["num_modes"] = 1

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
                    Path(conf_path).resolve().relative_to(
                        Path(cfg["CONFIG_RUN_DIR"]).resolve()
                    )
                except Exception:
                    raise RuntimeError(
                        f"Refusing to launch Vina with config outside current RUN_DIR: {conf_path}"
                    )

                logger.info(f"[vina.call] config={conf_path}")

                lig_n, out_n = norm(lig), norm(out_path)
                raw_docked_ligands[lig_n] = out_n
                futures[
                    pool.submit(
                        run_docking_task,
                        cfg["VINA_EXE"],
                        conf_path,
                        lig,
                        out_path,
                        write_failure_marker_flag=True,
                    )
                ] = (lig_n, out_n)

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
                    guard = guards_for_ligand.get(lig) or BudgetGuard(
                        default_budget_seconds
                    )

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

                    _per_ligand_postprocess(
                        lig=lig,
                        lig_name=lig_name,
                        out_path=out_path,
                        score=float(score),
                        guard=guard,
                        pbar=pbar,
                    )
                    continue

    stage_elapsed = time.time() - stage_start
    if run_id:
        try:
            update_manifest_for_docking_stage(
                cfg,
                run_id,
                pdb_id,
                variant_label,
                stage_name,
                status="completed",
                ph_tag=ph_label,
                elapsed_sec=stage_elapsed,
            )
        except Exception:
            logger.warning(
                "[run-manifest.docking-stage] failed to record completion pdb=%s variant=%s ph=%s stage=%s",
                pdb_id,
                variant_label,
                ph_label if ph_label is not None else "base",
                stage_name,
                exc_info=True,
            )

    return scores, validated_ligands, all_distances, raw_docked_ligands, invalids
