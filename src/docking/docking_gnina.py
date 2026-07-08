from __future__ import annotations

import logging
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from tqdm import tqdm

from docking.score_io import extract_gnina_scores
from docking.docking_gnina_support import (
    RetryManager,
    _map_reason_to_category,
    _resolve_gnina_stage_params,
    emit_gnina_config,
)
from docking.docking_ligands import validate_ligand
from docking.pose_validation import (
    compute_self_rmsd,
    extract_surface_atoms,
    filter_and_rewrite_poses_by_rmsd,
    _pose_centroid_from_pdbqt,
)
from docking.fallback_recenter import BudgetGuard, validate_first_valid_pose
from docking.docking_utils import run_completion_audit


def _run_gnina_for_ligand(
    cfg: Dict[str, Any],
    pdb_id: str,
    receptor_pdbqt: str,
    lig: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    stage: Dict[str, Any],
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
    logger: logging.Logger,
    threads_per_vina: int,
) -> Optional[str]:
    """
    Minimal GNINA follow-up docking for a single ligand.

    - Writes config under gnina_<stage_name> and runs GNINA with --no_gpu --cpu 1.
    - Returns the GNINA out_path on success, or None on failure.
    """
    gnina_exe = (cfg.get("GNINA_EXE") or "").strip()
    if not gnina_exe:
        logger.debug("[gnina.skip] GNINA_EXE not set; skipping GNINA follow-up.")
        return None

    # For fast test targets that use the small test_library_10, allow GNINA to use
    # more CPU threads to keep CI/local tests quick. This must be reflected both
    # in the emitted config and the CLI invocation.
    test_pdbs = {"TEST", "T3MP", "TEMP"}
    cpu_threads = 10 if pdb_id.strip().upper() in test_pdbs else 1

    # Build GNINA-specific stage params without mutating the shared stage dict.
    stage_for_gnina = _resolve_gnina_stage_params(cfg, stage)

    base_stage_name = stage.get("name") or "stage"
    gnina_stage_name = f"gnina_{base_stage_name}"

    logger.info(
        "[gnina.emit] pdb=%s stage=%s variant=%s ph=%s lig=%s",
        pdb_id,
        gnina_stage_name,
        variant_token or "None",
        ph_label or "None",
        os.path.basename(lig),
    )

    conf_path, out_path = emit_gnina_config(
        cfg,
        pdb_id,
        receptor_pdbqt,
        center,
        box_size,
        lig,
        gnina_stage_name,
        stage_for_gnina,
        cpu_threads,  # GNINA threads (test targets may use more)
        logger,
        variant=variant_token,
        ph_token=ph_label,
        legacy=legacy_mode,
    )

    # GNINA accepts a mostly Vina-compatible config, but rejects some Vina-only
    # options like energy_range/verbosity. Strip unsupported keys in-place to
    # keep router paths identical while ensuring GNINA runs.
    try:
        conf_p = Path(conf_path)
        raw_lines = conf_p.read_text(encoding="utf-8").splitlines()
        filtered_lines: List[str] = []
        removed_keys: List[str] = []
        for line in raw_lines:
            key = line.split("=", 1)[0].strip().lower()
            if key in {"energy_range", "verbosity"}:
                removed_keys.append(key)
                continue
            filtered_lines.append(line)
        if removed_keys:
            tmp = conf_p.with_suffix(".part")
            tmp.write_text("\n".join(filtered_lines) + "\n", encoding="utf-8")
            os.replace(tmp, conf_p)
            logger.debug(
                "[gnina.cfg.strip] removed=%s path=%s",
                ",".join(sorted(set(removed_keys))),
                conf_path,
            )
    except Exception as e:
        logger.warning("[gnina.cfg.strip] failed path=%s err=%s", conf_path, e)

    try:
        Path(conf_path).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
    except Exception:
        raise RuntimeError(
            f"Refusing to launch GNINA with config outside current RUN_DIR: {conf_path}"
        )

    # Build GNINA CLI. In FAST_MODE, explicitly disable CNN scoring to keep runs fast
    # and easier to test, and sort poses by energy instead of CNN score.
    cmd = [gnina_exe, "--no_gpu", "--cpu", str(cpu_threads), "--config", conf_path]
    if cfg.get("FAST_MODE"):
        cmd.extend(["--cnn_scoring", "none", "--pose_sort_order", "energy"])

    logger.info(
        "[gnina.call] exe=%s config=%s cmd=%s",
        gnina_exe,
        conf_path,
        " ".join(cmd),
    )

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as e:
        logger.warning(
            "[gnina.error] pdb=%s lig=%s stage=%s err=%s",
            pdb_id,
            os.path.basename(lig),
            gnina_stage_name,
            e,
        )
        return None

    return out_path


def run_gnina_for_stage(
    *,
    cfg: Dict[str, Any],
    paths: Any,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    stage_name: str,
    stage_info: Mapping[str, Any],
    ligands: List[str],
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    logger: logging.Logger,
    receptor_pdbqt: Optional[str] = None,
    control_lookup: Optional[Dict[str, Path]] = None,
) -> tuple[
    Dict[str, Optional[float]], Dict[str, Dict[str, Any]], Dict[str, Any]
]:
    """
    Run GNINA for the given stage and ligand list.

    - Writes GNINA config(s) into configs/.../gnina_<stage>/gnina.json
    - Writes docked poses into docked/.../gnina_<stage>/<ligand>.pdbqt
    - Uses GNINA_EXE from cfg
    - Uses --no_gpu and --cpu 1
    - Returns (primary_scores, metrics) where:
        primary_scores: {ligand_path: cnn_affinity or None}
        metrics: {
            ligand_path: {
                minimized_affinity_kcal, cnn_score, cnn_affinity_pK, gnina_primary_score
            }
        }
    """
    threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))
    cpu = int(cfg.get("CPU", os.cpu_count() or 1))
    max_workers = min(cpu, len(ligands)) if ligands else 1
    if max_workers < 1:
        max_workers = 1

    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )

    base_stage_name = stage_name
    gnina_stage_name = f"gnina_{base_stage_name}"
    stage_dir = paths.docked_stage_dir(variant_token, gnina_stage_name, ph_label)
    logger.info(
        "[gnina.stage] pdb=%s stage=%s variant=%s ph=%s n_lig=%d",
        pdb_id,
        gnina_stage_name,
        variant_token or "None",
        ph_label or "None",
        len(ligands),
    )

    receptor_pdbqt_path = receptor_pdbqt or str(
        paths.receptor_pdbqt(variant_token, ph_label)
    )
    try:
        surface_coords = extract_surface_atoms(
            pdbqt_path=receptor_pdbqt_path, center=center
        )
    except Exception:
        surface_coords = None

    def _best_pose_pdb_from_pdbqt(
        pdbqt_path: str, obabel_path: Optional[str] = None
    ) -> Optional[str]:
        """Convert first model of PDBQT -> PDB (no hydrogens) using OpenBabel."""
        try:
            from shutil import which

            obabel = obabel_path or which("obabel")
            if not obabel:
                return None

            out_pdb = Path(pdbqt_path).with_suffix(".best.pdb")
            cmd = [
                obabel,
                "-ipdbqt",
                pdbqt_path,
                "-opdb",
                "-O",
                str(out_pdb),
                "-d",
                "--first",
            ]
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            return str(out_pdb)
        except Exception:
            return None

    scores: Dict[str, Optional[float]] = {}
    gnina_metrics: Dict[str, Dict[str, Any]] = {}
    gnina_futures: Dict[Any, str] = {}

    default_budget_seconds = float(
        cfg.get("MAX_RETRY_SECONDS_PER_LIGAND", cfg.get("BENCH_MAX_SECONDS", 300.0))
    )
    guards_for_ligand: Dict[str, BudgetGuard] = {}
    retry_mgr = RetryManager()

    def _dock_once(
        lig: str,
        stage_params: Dict[str, Any],
        center_use: tuple[float, float, float],
        box_use: tuple[float, float, float],
        stage_label: str,
    ) -> tuple[Optional[float], Dict[str, Any], Dict[str, Any], Optional[str]]:
        """
        Run GNINA once, parse scores, and perform pose validation/self-RMSD/control gate.
        Returns (primary_score, metrics, validation_result, out_path).
        """
        out_path = _run_gnina_for_ligand(
            cfg,
            pdb_id,
            str(paths.receptor_pdbqt(variant_token, ph_label)),
            lig,
            center_use,
            box_use,
            stage_params,
            variant_token,
            ph_label,
            legacy_mode,
            logger,
            threads_per_vina,
        )

        metrics: Dict[str, Any] = {
            "minimized_affinity_kcal": None,
            "cnn_score": None,
            "cnn_affinity_pK": None,
            "gnina_primary_score": None,
            "valid": False,
            "reason": "",
            "self_rmsd": None,
        }
        validation_result: Dict[str, Any] = {"valid": False, "reason": "gnina_failed"}

        if not out_path:
            metrics["reason"] = "gnina_failed"
            return None, metrics, validation_result, None

        try:
            parsed = extract_gnina_scores(out_path)
        except Exception:
            parsed = {
                "minimized_affinity_kcal": None,
                "cnn_score": None,
                "cnn_affinity_pK": None,
            }

        minimized_affinity = parsed.get("minimized_affinity_kcal")
        cnn_score = parsed.get("cnn_score")
        cnn_affinity = parsed.get("cnn_affinity_pK")
        # Primary score is CNN affinity when available; otherwise fall back to the
        # minimized empirical affinity (e.g. when FAST_MODE disables CNN scoring).
        primary_score = cnn_affinity if cnn_affinity is not None else minimized_affinity

        metrics.update(
            {
                "minimized_affinity_kcal": minimized_affinity,
                "cnn_score": cnn_score,
                "cnn_affinity_pK": cnn_affinity,
                "gnina_primary_score": primary_score,
            }
        )
        logger.info(
            "[gnina.out] pdb=%s stage=%s lig=%s out=%s cnn_affinity=%s cnn_score=%s minimized_affinity=%s",
            pdb_id,
            stage_label,
            os.path.basename(lig),
            out_path,
            cnn_affinity if cnn_affinity is not None else "None",
            cnn_score if cnn_score is not None else "None",
            minimized_affinity if minimized_affinity is not None else "None",
        )

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

        try:
            validation_result = validate_first_valid_pose(
                receptor_pdbqt=receptor_pdbqt_path,
                ligand_pdbqt=out_path,
                pocket_center=center_use,
                surface_coords=surface_coords,
                max_models=int(cfg.get("EARLY_EXIT_MAX_MODELS", 3)),
                clash_threshold=2.0,
                clash_tol=3,
                dist_surf=6.0,
                dist_centroid=4.5,
            )
        except Exception:
            validation_result = {"valid": False, "reason": "pose_invalid"}

        valid_pose = bool(validation_result.get("valid", False))
        reason_str = validation_result.get("reason", "") or ""

        self_rmsd_val = None
        if bool(cfg.get("LOG_SELF_RMSD", True)):
            try:
                self_rmsd_val = compute_self_rmsd(out_path)
                logger.info(
                    f"[gnina.self-rmsd] lig={os.path.basename(lig)} rmsd={self_rmsd_val}"
                )
            except Exception as _e:
                logger.warning(
                    f"[gnina.self-rmsd] failed for {os.path.basename(lig)}: {_e}"
                )

        ctrl_ref = None
        if control_lookup:
            base = Path(lig).stem.split("_stage")[0]
            ctrl_ref = control_lookup.get(base)
        if ctrl_ref:
            best_pdb = _best_pose_pdb_from_pdbqt(
                out_path, obabel_path=cfg.get("OPENBABEL_PATH")
            )
            if not best_pdb:
                valid_pose = False
                reason_str = "no_best_pose_for_rmsd"
            else:
                ok = validate_ligand(
                    ligand_name=os.path.basename(lig),
                    docked_path=best_pdb,
                    crystal_path=str(ctrl_ref),
                    rmsd_thresh=float(cfg.get("CONTROL_RMSD_MAX_ANG", 2.0)),
                    self_rmsd=self_rmsd_val,
                    logger=logger,
                )
                if not ok:
                    valid_pose = False
                    reason_str = "rmsd_fail"

        metrics["valid"] = bool(valid_pose)
        metrics["reason"] = reason_str
        metrics["self_rmsd"] = self_rmsd_val
        return primary_score, metrics, validation_result, out_path

    def _near_miss_retry(
        lig: str,
        base_stage: Dict[str, Any],
        result: Dict[str, Any],
        guard: BudgetGuard,
        last_out: Optional[str],
    ) -> Optional[Dict[str, Any]]:
        if not bool(cfg.get("RETRY_NEAR_MISS", True)):
            return None

        reason = result.get("reason", "") or ""
        near_miss = (
            ("clash" in reason)
            or (result.get("distance_to_surface") or 0.0) < 6.5
            or (result.get("distance_to_centroid") or 0.0) < 4.0
        )
        if not near_miss:
            return None

        stage_retry = dict(base_stage)
        stage_retry["name"] = f"{base_stage.get('name', base_stage_name)}_retry"
        stage_retry["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
        if cfg.get("FAST_MODE"):
            stage_retry["exhaustiveness"] = 1

        try:
            base_seed = int(stage_retry.get("seed", 0)) if "seed" in stage_retry else 0
        except Exception:
            base_seed = 0
        stage_retry["seed"] = base_seed + 137
        stage_retry["num_modes"] = int(
            cfg.get("NEAR_MISS_NUM_MODES", stage_retry.get("num_modes", 4))
        )
        stage_retry["energy_range"] = float(
            cfg.get("NEAR_MISS_ENERGY_RANGE", stage_retry.get("energy_range", 4.0))
        )

        center_nm = center
        box_nm = box_size
        if bool(cfg.get("NEAR_MISS_RECENTER", True)) and last_out:
            try:
                cent = _pose_centroid_from_pdbqt(str(last_out))
            except Exception as _e:
                logger.warning(
                    f"[near-miss] failed to compute centroid for {os.path.basename(lig)}: {_e}"
                )
                cent = None
            if cent and isinstance(cent, (list, tuple)) and len(cent) == 3:
                try:
                    center_nm = (
                        float(cent[0]),
                        float(cent[1]),
                        float(cent[2]),
                    )
                    max_box = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                    box_nm = (
                        min(max_box, float(box_size[0]) + 1.0),
                        min(max_box, float(box_size[1]) + 1.0),
                        min(max_box, float(box_size[2]) + 1.0),
                    )
                except Exception:
                    center_nm = center
                    box_nm = box_size

        if guard.expired():
            return {
                "minimized_affinity_kcal": None,
                "cnn_score": None,
                "cnn_affinity_pK": None,
                "gnina_primary_score": None,
                "valid": False,
                "reason": "budget_exceeded",
                "self_rmsd": None,
            }

        primary_retry, metrics_retry, _, _ = _dock_once(
            lig,
            stage_retry,
            center_nm,
            box_nm,
            f"{stage_retry['name']}",
        )
        if metrics_retry.get("valid", False):
            logger.info(
                "%s | %s score: %s (gnina near-miss rescued)",
                os.path.basename(lig),
                stage_retry["name"],
                f"{primary_retry:.2f}"
                if isinstance(primary_retry, (int, float))
                else "None",
            )
        return metrics_retry

    def _structured_retries(
        lig: str,
        base_stage: Dict[str, Any],
        result: Dict[str, Any],
        guard: BudgetGuard,
        last_out: Optional[str],
    ) -> Optional[Dict[str, Any]]:
        err_cat = _map_reason_to_category(result.get("reason", ""))
        attempt = 0
        metrics_candidate: Optional[Dict[str, Any]] = None
        last_out_path = last_out

        while attempt < retry_mgr.max_retries:
            if guard.expired():
                metrics_candidate = {
                    "minimized_affinity_kcal": None,
                    "cnn_score": None,
                    "cnn_affinity_pK": None,
                    "gnina_primary_score": None,
                    "valid": False,
                    "reason": "budget_exceeded",
                    "self_rmsd": None,
                }
                break

            recipe = retry_mgr.apply(base_stage, err_cat, attempt)
            if not recipe:
                break

            stage_retry2 = dict(base_stage)
            stage_retry2[
                "name"
            ] = f"{base_stage.get('name', base_stage_name)}_r{attempt + 1}"
            stage_retry2["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))

            if "exhaustiveness" in recipe:
                stage_retry2["exhaustiveness"] = recipe["exhaustiveness"]
            if "num_modes" in recipe:
                stage_retry2["num_modes"] = recipe["num_modes"]
            if recipe.get("energy_range") is not None:
                stage_retry2["energy_range"] = recipe["energy_range"]
            if cfg.get("FAST_MODE"):
                stage_retry2["exhaustiveness"] = 1
            if recipe.get("seed_jitter", False):
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
                if recipe.get("recenter", False) and last_out_path:
                    cent = _pose_centroid_from_pdbqt(str(last_out_path))
                    if cent and isinstance(cent, (list, tuple)) and len(cent) == 3:
                        retry_center = (
                            float(cent[0]),
                            float(cent[1]),
                            float(cent[2]),
                        )
                if "box_pad_delta" in recipe and isinstance(
                    recipe["box_pad_delta"], (int, float)
                ):
                    dx = float(recipe["box_pad_delta"])
                    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                    retry_box = (
                        min(box_cap, float(box_size[0]) + dx),
                        min(box_cap, float(box_size[1]) + dx),
                        min(box_cap, float(box_size[2]) + dx),
                    )
            except Exception as _e:
                logger.warning(f"Retry recenter/box tweak failed: {_e}")

            primary_retry, metrics_retry, result_retry, out_retry = _dock_once(
                lig,
                stage_retry2,
                retry_center,
                retry_box,
                f"{stage_retry2['name']}",
            )
            last_out_path = out_retry or last_out_path
            metrics_candidate = metrics_retry
            if metrics_retry.get("valid", False):
                logger.info(
                    "%s | %s score: %s (gnina structured retry rescued)",
                    os.path.basename(lig),
                    stage_retry2["name"],
                    f"{primary_retry:.2f}"
                    if isinstance(primary_retry, (int, float))
                    else "None",
                )
                break
            result = result_retry
            attempt += 1

        return metrics_candidate

    submit_ligands: List[str] = []
    for lig in ligands:
        guard = guards_for_ligand.get(lig) or BudgetGuard(
            default_budget_seconds, log=logger
        )
        guards_for_ligand[lig] = guard
        if guard.expired():
            scores[lig] = None
            gnina_metrics[lig] = {
                "minimized_affinity_kcal": None,
                "cnn_score": None,
                "cnn_affinity_pK": None,
                "gnina_primary_score": None,
                "valid": False,
                "reason": "budget_exceeded",
                "self_rmsd": None,
            }
            logger.info(
                "[gnina.budget] pdb=%s lig=%s stage=%s reason=budget_exceeded",
                pdb_id,
                lig,
                gnina_stage_name,
            )
            continue
        submit_ligands.append(lig)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        for lig in submit_ligands:
            gnina_futures[
                pool.submit(
                    _dock_once,
                    lig,
                    dict(stage_info),
                    center,
                    box_size,
                    gnina_stage_name,
                )
            ] = lig

        with tqdm(
            total=len(gnina_futures),
            desc=f"GNINA ({base_stage_name})",
            unit="ligand",
            position=2,
            dynamic_ncols=True,
            mininterval=0.2,
            leave=True,
            file=sys.stdout,
        ) as pbar:
            for fut in as_completed(gnina_futures):
                lig = gnina_futures[fut]
                lig_name = os.path.basename(lig)
                try:
                    primary, metrics, result, out_gnina = fut.result()
                except Exception as e:
                    logger.warning(
                        "[gnina.error] pdb=%s lig=%s stage=%s err=%s",
                        pdb_id,
                        lig_name,
                        gnina_stage_name,
                        e,
                    )
                    _primary, metrics, result, out_gnina = (
                        None,
                        {
                            "minimized_affinity_kcal": None,
                            "cnn_score": None,
                            "cnn_affinity_pK": None,
                            "gnina_primary_score": None,
                            "valid": False,
                            "reason": "gnina_failed",
                            "self_rmsd": None,
                        },
                        {"valid": False, "reason": "gnina_failed"},
                        None,
                    )

                guard = guards_for_ligand.get(lig) or BudgetGuard(
                    default_budget_seconds, log=logger
                )
                guards_for_ligand[lig] = guard

                final_metrics = metrics
                if not metrics.get("valid", False):
                    if guard.expired():
                        final_metrics = dict(metrics)
                        final_metrics["reason"] = "budget_exceeded"
                    else:
                        near_miss_metrics = _near_miss_retry(
                            lig, dict(stage_info), result, guard, out_gnina
                        )
                        if near_miss_metrics and near_miss_metrics.get("valid", False):
                            final_metrics = near_miss_metrics
                        else:
                            structured_metrics = _structured_retries(
                                lig, dict(stage_info), result, guard, out_gnina
                            )
                            if structured_metrics and structured_metrics.get(
                                "valid", False
                            ):
                                final_metrics = structured_metrics
                            elif structured_metrics:
                                final_metrics = structured_metrics

                scores[lig] = final_metrics.get("gnina_primary_score")
                gnina_metrics[lig] = final_metrics
                pbar.update(1)

    # Completion audit: rerun missing GNINA outputs once before checkpointing.
    stage_params_for_gnina = _resolve_gnina_stage_params(cfg, stage_info)

    def _expected_gnina_path(lig: str) -> Path:
        return stage_dir / f"{Path(lig).stem}_{gnina_stage_name}.pdbqt"

    def _rerun_gnina_missing(lig: str) -> tuple[bool, Optional[str], Optional[Path]]:
        try:
            primary_r, metrics_r, _result_r, out_path = _dock_once(
                lig,
                dict(stage_params_for_gnina),
                center,
                box_size,
                gnina_stage_name,
            )
            gnina_metrics[lig] = metrics_r
            scores[lig] = metrics_r.get("gnina_primary_score")
        except Exception as exc:
            return False, f"rerun_error:{exc}", None
        if out_path and Path(out_path).exists() and Path(out_path).stat().st_size > 0:
            return True, "rerun_ok", None
        return False, "rerun_no_output", None

    completion_report_gnina = run_completion_audit(
        engine="gnina",
        pdb_id=pdb_id,
        stage_name=gnina_stage_name,
        ligands=ligands,
        expected_output_path=_expected_gnina_path,
        rerun_one=_rerun_gnina_missing,
        stage_dir=stage_dir,
        cfg=cfg,
        logger=logger,
        retries=1,
        ph_label=ph_label,
        variant=variant_token,
    )
    missing_after = completion_report_gnina.get("missing_ligands_after") or []
    if missing_after:
        for lig in missing_after:
            if lig not in gnina_metrics:
                gnina_metrics[lig] = {
                    "minimized_affinity_kcal": None,
                    "cnn_score": None,
                    "cnn_affinity_pK": None,
                    "gnina_primary_score": None,
                    "valid": False,
                    "reason": "completion_missing",
                    "self_rmsd": None,
                }
            scores.setdefault(lig, gnina_metrics[lig].get("gnina_primary_score"))

    return scores, gnina_metrics, completion_report_gnina
