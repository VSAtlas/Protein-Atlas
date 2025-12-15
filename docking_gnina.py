from __future__ import annotations

import logging
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from tqdm import tqdm

from input_and_export_functions import _to_bool, emit_gnina_config, extract_gnina_scores
from chemdb.target_difficulty import TargetDifficulty
from docking_ligands import compute_rmsd, validate_ligand
from pose_validation import compute_self_rmsd, extract_surface_atoms, filter_and_rewrite_poses_by_rmsd
from fallback_recenter import validate_first_valid_pose

HARD_BUCKETS = {"hard", "degenerate"}


def should_run_gnina_for_target(td: Optional[TargetDifficulty], cfg: Mapping[str, Any]) -> bool:
    """
    Decide whether to run GNINA based on difficulty and config.

    New behavior: honor config flag USE_GNINA / use_gnina.
    If set false-y, GNINA is globally disabled even if difficulty would request it.

    - GNINA only runs if td.difficulty is "hard" or "degenerate".
    - Also requires cfg.get("ENABLE_GNINA", True) and a valid GNINA_EXE.
    - Unknown difficulty is treated as "easy" (skip).
    """
    # --- Global GNINA on/off switch via config ---
    raw_flag = None
    if cfg is not None:
        # If the config file explicitly disables GNINA, that takes precedence
        # over any environment overrides.
        file_cfg = cfg.get("_FILE_CFG") if isinstance(cfg, dict) else None
        if isinstance(file_cfg, dict):
            raw_flag = file_cfg.get("USE_GNINA", file_cfg.get("use_gnina", None))
        if raw_flag is None:
            raw_flag = cfg.get("USE_GNINA", cfg.get("use_gnina", None))
    if raw_flag is not None and not _to_bool(raw_flag):
        logging.getLogger(__name__).info(
            "[gnina.disabled] USE_GNINA=%r -> GNINA follow-up docking is globally disabled",
            raw_flag,
        )
        return False

    if not cfg.get("ENABLE_GNINA", True):
        return False
    if not (cfg.get("GNINA_EXE") or "").strip():
        return False
    if td is None:
        return False
    diff = (getattr(td, "difficulty", "") or "").strip().lower()
    return diff in HARD_BUCKETS


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
    stage_for_gnina = dict(stage)
    stage_for_gnina["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
    if cfg.get("FAST_MODE"):
        # FAST_MODE for GNINA forces both knobs to 1; Vina remains unchanged.
        stage_for_gnina["exhaustiveness"] = 1
        stage_for_gnina["num_modes"] = 1

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

    logger.info("[gnina.call] exe=%s config=%s", gnina_exe, conf_path)

    cmd = [gnina_exe, "--no_gpu", "--cpu", str(cpu_threads), "--config", conf_path]
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
) -> tuple[Dict[str, Optional[float]], Dict[str, Dict[str, Optional[float]]]]:
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
    max_jobs = int(cfg.get("MAX_PARALLEL_JOBS", cpu))
    max_workers = min(max_jobs, len(ligands)) if ligands else 1
    if max_workers < 1:
        max_workers = 1

    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None

    base_stage_name = stage_name
    gnina_stage_name = f"gnina_{base_stage_name}"
    logger.info(
        "[gnina.stage] pdb=%s stage=%s variant=%s ph=%s n_lig=%d",
        pdb_id,
        gnina_stage_name,
        variant_token or "None",
        ph_label or "None",
        len(ligands),
    )

    receptor_pdbqt_path = receptor_pdbqt or str(paths.receptor_pdbqt(variant_token, ph_label))
    try:
        surface_coords = extract_surface_atoms(pdbqt_path=receptor_pdbqt_path, center=center)
    except Exception:
        surface_coords = None

    def _best_pose_pdb_from_pdbqt(pdbqt_path: str, obabel_path: Optional[str] = None) -> Optional[str]:
        """Convert first model of PDBQT -> PDB (no hydrogens) using OpenBabel."""
        try:
            from shutil import which

            obabel = obabel_path or which("obabel")
            if not obabel:
                return None

            out_pdb = Path(pdbqt_path).with_suffix(".best.pdb")
            cmd = [obabel, "-ipdbqt", pdbqt_path, "-opdb", "-O", str(out_pdb), "-d", "--first"]
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            return str(out_pdb)
        except Exception:
            return None

    scores: Dict[str, Optional[float]] = {}
    gnina_metrics: Dict[str, Dict[str, Optional[float]]] = {}
    gnina_futures: Dict[Any, str] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        for lig in ligands:
            gnina_futures[
                pool.submit(
                    _run_gnina_for_ligand,
                    cfg,
                    pdb_id,
                    str(paths.receptor_pdbqt(variant_token, ph_label)),
                    lig,
                    center,
                    box_size,
                    dict(stage_info),
                    variant_token,
                    ph_label,
                    legacy_mode,
                    logger,
                    threads_per_vina,
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
                    out_gnina = fut.result()
                except Exception as e:
                    logger.warning(
                        "[gnina.error] pdb=%s lig=%s stage=%s err=%s",
                        pdb_id,
                        lig_name,
                        gnina_stage_name,
                        e,
                    )
                    out_gnina = None

                best = None
                if out_gnina:
                    try:
                        parsed = extract_gnina_scores(out_gnina)
                    except Exception:
                        parsed = {
                            "minimized_affinity_kcal": None,
                            "cnn_score": None,
                            "cnn_affinity_pK": None,
                        }
                    minimized_affinity = parsed.get("minimized_affinity_kcal")
                    cnn_score = parsed.get("cnn_score")
                    cnn_affinity = parsed.get("cnn_affinity_pK")
                    best = cnn_affinity
                    logger.info(
                        "[gnina.out] pdb=%s stage=%s lig=%s out=%s cnn_affinity=%s cnn_score=%s minimized_affinity=%s",
                        pdb_id,
                        gnina_stage_name,
                        lig_name,
                        out_gnina,
                        cnn_affinity if cnn_affinity is not None else "None",
                        cnn_score if cnn_score is not None else "None",
                        minimized_affinity if minimized_affinity is not None else "None",
                    )
                    gnina_metrics[lig] = {
                        "minimized_affinity_kcal": minimized_affinity,
                        "cnn_score": cnn_score,
                        "cnn_affinity_pK": cnn_affinity,
                        "gnina_primary_score": cnn_affinity,
                        "valid": False,
                        "reason": "",
                        "self_rmsd": None,
                    }
                    # Pose filtering + validation (mirror Vina)
                    try:
                        kept, removed = filter_and_rewrite_poses_by_rmsd(
                            out_gnina,
                            rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                            max_models=int(cfg.get("RMSD_MAX_MODELS", 3)),
                        )
                        if removed > 0:
                            logger.info(f"{os.path.basename(out_gnina)}: RMSD filter kept {kept}, removed {removed}")
                    except Exception as e:
                        logger.warning(f"RMSD filtering failed for {os.path.basename(out_gnina)}: {e}")

                    try:
                        result = validate_first_valid_pose(
                            receptor_pdbqt=receptor_pdbqt_path,
                            ligand_pdbqt=out_gnina,
                            pocket_center=center,
                            surface_coords=surface_coords,
                            max_models=int(cfg.get("EARLY_EXIT_MAX_MODELS", 3)),
                            clash_threshold=2.0,
                            clash_tol=3,
                            dist_surf=6.0,
                            dist_centroid=4.5,
                        )
                    except Exception:
                        result = {"valid": False, "reason": "pose_invalid"}

                    valid_pose = bool(result.get("valid", False))
                    reason_str = result.get("reason", "") or ""

                    self_rmsd_val = None
                    if bool(cfg.get("LOG_SELF_RMSD", True)):
                        try:
                            self_rmsd_val = compute_self_rmsd(out_gnina)
                            logger.info(f"[gnina.self-rmsd] lig={lig_name} rmsd={self_rmsd_val}")
                        except Exception as _e:
                            logger.warning(f"[gnina.self-rmsd] failed for {lig_name}: {_e}")

                    # Control redock gate (hard RMSD gate)
                    ctrl_ref = None
                    if control_lookup:
                        base = Path(lig).stem.split("_stage")[0]
                        ctrl_ref = control_lookup.get(base)
                    if ctrl_ref:
                        best_pdb = _best_pose_pdb_from_pdbqt(out_gnina, obabel_path=cfg.get("OPENBABEL_PATH"))
                        if not best_pdb:
                            valid_pose = False
                            reason_str = "no_best_pose_for_rmsd"
                        else:
                            ok = validate_ligand(
                                ligand_name=lig_name,
                                docked_path=best_pdb,
                                crystal_path=str(ctrl_ref),
                                rmsd_thresh=float(cfg.get("CONTROL_RMSD_MAX_ANG", 2.0)),
                                self_rmsd=self_rmsd_val,
                                logger=logger,
                            )
                            if not ok:
                                valid_pose = False
                                reason_str = "rmsd_fail"

                    gnina_metrics[lig]["valid"] = bool(valid_pose)
                    gnina_metrics[lig]["reason"] = reason_str
                    gnina_metrics[lig]["self_rmsd"] = self_rmsd_val
                else:
                    gnina_metrics[lig] = {
                        "minimized_affinity_kcal": None,
                        "cnn_score": None,
                        "cnn_affinity_pK": None,
                        "gnina_primary_score": None,
                        "valid": False,
                        "reason": "gnina_failed",
                        "self_rmsd": None,
                    }
                scores[lig] = best
                pbar.update(1)

    return scores, gnina_metrics
