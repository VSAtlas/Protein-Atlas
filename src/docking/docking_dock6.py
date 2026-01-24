from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import math
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from input_and_export_functions import _to_bool, write_score_summary_to_csv
from path_router.path_router import make_paths, ph_ensemble_dir
from prep_docking.prep_for_ledock import ensure_mol2_for_ledock, map_pdbqt_to_mol2_path


def _get_dock6_exe(cfg: Dict[str, Any], logger: Optional[logging.Logger] = None) -> str:
    log = logger or logging.getLogger(__name__)
    candidate = None
    if isinstance(cfg, dict):
        candidate = str(cfg.get("DOCK6_EXE") or "").strip() or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = str(file_cfg.get("DOCK6_EXE") or "").strip() or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/bin/dock6"
    if candidate and not Path(candidate).exists():
        log.warning("[dock6.exe.missing] candidate=%s falling_back_to=dock6", candidate)
        candidate = "dock6"
    log.info("[dock6.exe] path=%s", candidate)
    return candidate


def _get_vdw_defn_file(cfg: Dict[str, Any], logger: logging.Logger) -> str:
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("DOCK6_VDW_DEFN_FILE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("DOCK6_VDW_DEFN_FILE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/parameters/vdw_AMBER_parm99.defn"
    logger.info("[dock6.vdw_defn] path=%s", candidate)
    return str(candidate)


def _get_flex_defn_file(cfg: Dict[str, Any], logger: logging.Logger) -> str:
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("DOCK6_FLEX_DEFN_FILE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("DOCK6_FLEX_DEFN_FILE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/parameters/flex.defn"
    logger.info("[dock6.flex_defn] path=%s", candidate)
    return str(candidate)


def _get_flex_drive_file(cfg: Dict[str, Any], logger: logging.Logger) -> str:
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("DOCK6_FLEX_DRIVE_FILE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("DOCK6_FLEX_DRIVE_FILE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/parameters/flex_drive.tbl"
    logger.info("[dock6.flex_drive] path=%s", candidate)
    return str(candidate)


def _normalize_ph_label(ph_label: Optional[str]) -> Optional[str]:
    if ph_label is None:
        return None
    token = str(ph_label).strip()
    return token or None


def _variant_for_ph(variant: Optional[str], legacy_mode: bool) -> Optional[str]:
    v = (str(variant).strip().upper() or None) if variant is not None else None
    if v:
        return v
    if legacy_mode:
        return None
    return "HOLO"


def should_run_dock6_for_target(cfg: Dict[str, Any]) -> bool:
    raw_flag = None
    if isinstance(cfg, dict):
        file_cfg = cfg.get("_FILE_CFG") if isinstance(cfg, dict) else None
        raw_flag = cfg.get("USE_DOCK6", cfg.get("use_dock6", None))
        if raw_flag is None and isinstance(file_cfg, dict):
            raw_flag = file_cfg.get("USE_DOCK6", file_cfg.get("use_dock6", None))
    if raw_flag is None:
        env_raw = os.environ.get("USE_DOCK6", os.environ.get("use_dock6"))
        raw_flag = env_raw

    if not _to_bool(raw_flag):
        logging.getLogger(__name__).info(
            "[dock6.disabled] USE_DOCK6=%r -> DOCK6 docking is disabled", raw_flag
        )
        return False

    if not bool(cfg.get("PH_ENSEMBLE")):
        logging.getLogger(__name__).info("[dock6.skip] reason=no_ph_ensemble")
        return False
    canonical = cfg.get("_PH_ENSEMBLE_CANONICAL")
    if isinstance(canonical, dict) and not canonical:
        logging.getLogger(__name__).info("[dock6.skip] reason=empty_ph_ensemble")
        return False

    exe = _get_dock6_exe(cfg)
    if not exe:
        logging.getLogger(__name__).warning("[dock6.skip] reason=missing_exe")
        return False
    return True


def get_dock6_stage_preset(cfg: Dict[str, Any], stage_key: str) -> Dict[str, Any]:
    fast_mode = bool(cfg.get("FAST_MODE"))
    fast_preset = {
        "grid_prefix": "grid_fast",
        "grid_spacing": 0.60,
        "max_orientations": 150,
        "minimize_ligand": False,
        "simplex_max_iterations": 1000,
        "simplex_max_cycles": 1,
        "num_final_scored_poses": 1,
    }
    if fast_mode:
        return fast_preset

    stage = (stage_key or "").lower()
    presets = {
        "stage1": {
            "grid_prefix": "grid_s1",
            "grid_spacing": 0.50,
            "max_orientations": 300,
            "minimize_ligand": False,
            "simplex_max_iterations": 1000,
            "simplex_max_cycles": 1,
            "num_final_scored_poses": 1,
        },
        "stage2": {
            "grid_prefix": "grid_s2",
            "grid_spacing": 0.35,
            "max_orientations": 1000,
            "minimize_ligand": True,
            "simplex_max_iterations": 500,
            "simplex_max_cycles": 1,
            "num_final_scored_poses": 1,
        },
        "stage3": {
            "grid_prefix": "grid_s3",
            "grid_spacing": 0.30,
            "max_orientations": 3000,
            "minimize_ligand": True,
            "simplex_max_iterations": 1500,
            "simplex_max_cycles": 2,
            "num_final_scored_poses": 3,
        },
    }
    return presets.get(stage, presets["stage2"])


def _build_multi_mol2_for_stage(
    ligand_pairs: List[Tuple[Path, Path]],
    out_path: Path,
    logger: logging.Logger,
) -> None:
    """
    Build a multi-ligand MOL2 where each molecule name matches the mol2 filename
    (avoids leaking .pdbqt names into DOCK6 outputs).
    """
    parts: List[str] = []
    for pdbqt_path, mol2_path in ligand_pairs:
        try:
            text = mol2_path.read_text(encoding="utf-8", errors="ignore")
        except Exception as exc:
            logger.warning("[dock6.mol2.skip] ligand=%s reason=%s", mol2_path, exc)
            continue

        blocks = text.split("@<TRIPOS>MOLECULE")
        if len(blocks) < 2:
            logger.warning(
                "[dock6.mol2.skip] ligand=%s reason=no_molecule_block", mol2_path
            )
            continue
        # Drop any prefix before first molecule
        for block in blocks[1:]:
            lines = [l for l in block.splitlines() if l.strip()]
            if not lines:
                continue
            # first line is molecule name; replace with mol2 filename to keep DOCK6 outputs .mol2-based
            lines[0] = mol2_path.name
            parts.append("@<TRIPOS>MOLECULE\n" + "\n".join(lines) + "\n")
            break

    out_path.write_text("\n".join(parts), encoding="utf-8")
    logger.info("[dock6.multimol2] path=%s ligands=%d", out_path, len(ligand_pairs))


def _parse_dock6_ranked_mol2(
    ranked_path: Path,
    ligand_pairs: List[Tuple[Path, Path]],
    logger: logging.Logger,
) -> Dict[Path, Dict[str, Any]]:
    """
    Parse a DOCK6 *_ranked.mol2 file and extract grid scores per ligand.

    Supports header lines of the form:
      '########## Name:ligand.mol2'
      '########## Name: ligand.mol2'
      '########## Name ligand.mol2'
    by extracting the ligand name after 'Name' or after the first ':'.
    """
    mol2_by_pdbqt: Dict[Path, Path] = {pdbqt: mol2 for pdbqt, mol2 in ligand_pairs}

    # Build a name map that includes both pdbqt and mol2 basenames/stems so either can match.
    name_map: Dict[str, Path] = {}
    for pdbqt_path, mol2_path in ligand_pairs:
        for token in {
            pdbqt_path.name,
            pdbqt_path.stem,
            mol2_path.name,
            mol2_path.stem,
        }:
            if token:
                name_map[token] = pdbqt_path

    metrics: Dict[Path, Dict[str, Any]] = {}
    parsed: List[Tuple[str, Optional[float]]] = []

    try:
        current_name: Optional[str] = None
        current_score: Optional[float] = None

        with ranked_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                # Extract the ligand name
                if "Name" in line:
                    tokens = line.replace("#", " ").replace("\t", " ").split()
                    for idx, tok in enumerate(tokens):
                        low = tok.lower()
                        if low.startswith("name"):
                            name_val: Optional[str] = None

                            # Case 1: token looks like "Name:ligand" or "Name:ligand.mol2"
                            if ":" in tok:
                                parts = tok.split(":", 1)
                                if len(parts) == 2 and parts[1].strip():
                                    name_val = parts[1].strip()

                            # Case 2: "Name" or "Name:" followed by the ligand name as next token
                            if not name_val and idx + 1 < len(tokens):
                                candidate = tokens[idx + 1].strip()
                                # Avoid treating a bare ":" as the name
                                if candidate and candidate != ":":
                                    name_val = candidate

                            if name_val:
                                current_name = name_val
                            break

                # Extract the grid score
                if "Grid_Score" in line:
                    tokens = line.replace("#", " ").replace("\t", " ").split()
                    for idx, tok in enumerate(tokens):
                        if tok.lower().startswith("grid_score"):
                            if idx + 1 < len(tokens):
                                try:
                                    current_score = float(tokens[idx + 1])
                                except Exception:
                                    current_score = None
                            break

                    # Only record a row when we have both name and score
                    if current_name is not None and current_score is not None:
                        parsed.append((current_name, current_score))
                        current_name = None
                        current_score = None

    except Exception as exc:
        logger.warning("[dock6.parse.error] path=%s reason=%s", ranked_path, exc)

    used: set[Path] = set()

    # First pass: map by molecule name
    for name, score in parsed:
        match = name_map.get(name)
        if match and match not in used:
            metrics[match] = {
                "grid_score": score,
                "n_poses": 1 if score is not None else 0,
                "valid": score is not None,
                "reason": None if score is not None else "parse_error",
                "mol2_name": name,
            }
            used.add(match)

    # Second pass: fall back to positional matching for anything not seen in the ranked file
    for idx, (pdbqt_path, _mol2_path) in enumerate(ligand_pairs):
        if pdbqt_path in metrics:
            continue
        score: Optional[float] = None
        if idx < len(parsed):
            score = parsed[idx][1]

        metrics[pdbqt_path] = {
            "grid_score": score,
            "n_poses": 1 if score is not None else 0,
            "valid": score is not None,
            "reason": None if score is not None else "missing_from_ranked",
            "mol2_name": mol2_by_pdbqt.get(pdbqt_path, pdbqt_path).name,
        }

    # Ensure mol2_name is populated for all metrics
    for lig, rec in metrics.items():
        if "mol2_name" not in rec or not rec.get("mol2_name"):
            rec["mol2_name"] = mol2_by_pdbqt.get(lig, lig).name

    return metrics


def run_dock6_for_stage(
    cfg: Dict[str, Any],
    paths: Any,
    pdb_id: str,
    *,
    variant: Optional[str],
    ph_label: Optional[str],
    stage_name: str,
    stage_info: Dict[str, Any],
    ligands: List[Path],
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    logger: logging.Logger,
    receptor_pdb: Optional[Path] = None,
    skip_completion: bool = False,
    output_root: Optional[Path] = None,
) -> Tuple[Dict[Path, Optional[float]], Dict[Path, Dict[str, Any]]]:
    scores: Dict[Path, Optional[float]] = {}
    metrics: Dict[Path, Dict[str, Any]] = {}

    if not ligands or center is None or box_size is None:
        logger.info("[dock6.skip] reason=missing_ligands_or_center")
        return scores, metrics

    if not should_run_dock6_for_target(cfg):
        return scores, metrics

    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    variant_for_ph = _variant_for_ph(variant_token, legacy_mode)
    ph_token = _normalize_ph_label(ph_label)

    # Ensure mol2 prep for all ligands up front so we always draw from prepped_ligands/<lib>/<lib>_mol2
    try:
        ensure_mol2_for_ledock(cfg, ligands, logger)
    except Exception as exc:
        logger.warning("[dock6.mol2.prep.warn] n=%d reason=%s", len(ligands), exc)

    ensemble_dir = ph_ensemble_dir(pdb_id, variant=variant_for_ph, legacy=legacy_mode)
    dock6_root = ensemble_dir / "dock6"
    rec_ms = dock6_root / "rec.ms"
    if not rec_ms.exists() or rec_ms.stat().st_size == 0:
        logger.warning(
            "[dock6.skip] reason=site_prep_missing pdb=%s variant=%s ph=%s dir=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token or "base",
            dock6_root,
        )
        return scores, metrics

    ligand_pairs: List[Tuple[Path, Path]] = []
    missing_mol2: List[Path] = []
    for lig in ligands:
        mol2_path = map_pdbqt_to_mol2_path(lig, logger=logger)
        if mol2_path is None or not mol2_path.exists():
            missing_mol2.append(lig)
            continue
        ligand_pairs.append((lig, mol2_path))

    for lig in missing_mol2:
        metrics[lig] = {
            "grid_score": None,
            "n_poses": 0,
            "valid": False,
            "reason": "dock6_missing_mol2",
        }

    if not ligand_pairs:
        logger.info("[dock6.skip] reason=no_mol2_ligands")
        return scores, metrics

    mol2_by_pdbqt = {pdbqt: mol2 for pdbqt, mol2 in ligand_pairs}
    stage_key = stage_info.get("key") or stage_name
    stage_preset = get_dock6_stage_preset(cfg, stage_key)
    grid_prefix = stage_preset["grid_prefix"]
    grid_nrg = dock6_root / f"{grid_prefix}.nrg"
    if not grid_nrg.exists() or grid_nrg.stat().st_size == 0:
        logger.warning(
            "[dock6.skip] reason=missing_grid pdb=%s prefix=%s dir=%s",
            pdb_id,
            grid_prefix,
            dock6_root,
        )
        return scores, metrics
    prefix = f"dock6_{stage_key}"
    ligand_items = list(ligand_pairs)
    n_lig = len(ligand_items)

    cpu = int(cfg.get("CPU", os.cpu_count() or 1))
    max_jobs = int(cfg.get("MAX_PARALLEL_JOBS", cpu))
    n_workers_cap = max(1, min(cpu, max_jobs, n_lig))

    raw_chunks = int(cfg.get("DOCK6_CHUNK_SIZE", 0) or 0)
    if raw_chunks > 0:
        target_batches = max(1, min(raw_chunks, n_lig))
    else:
        target_batches = max(1, min(n_workers_cap, n_lig))

    chunk_size = max(1, math.ceil(n_lig / target_batches))

    batches: List[List[Tuple[Path, Path]]] = []
    for i in range(0, n_lig, chunk_size):
        batches.append(ligand_items[i : i + chunk_size])

    n_batches = len(batches)
    max_workers = min(n_workers_cap, n_batches)

    logger.info(
        "[dock6.stage] pdb=%s stage=%s variant=%s ph=%s n_lig=%d n_batches=%d max_workers=%d",
        pdb_id,
        stage_name,
        variant_for_ph or "HOLO",
        ph_token or "ph_ensemble",
        n_lig,
        n_batches,
        max_workers,
    )

    dock6_exe = _get_dock6_exe(cfg)
    vdw_defn = _get_vdw_defn_file(cfg, logger)
    flex_defn = _get_flex_defn_file(cfg, logger)
    flex_drive = _get_flex_drive_file(cfg, logger)
    max_orientations = int(stage_preset.get("max_orientations", 1000))
    minimize_flag = "yes" if stage_preset.get("minimize_ligand", True) else "no"
    simplex_max_iterations = int(stage_preset.get("simplex_max_iterations", 1000))
    simplex_max_cycles = int(stage_preset.get("simplex_max_cycles", 1))
    num_final_scored_poses = int(stage_preset.get("num_final_scored_poses", 1))

    logger.info(
        "[dock6.run] pdb=%s stage=%s ligands=%d exe=%s",
        pdb_id,
        stage_name,
        len(ligand_pairs),
        dock6_exe,
    )

    start_ts = time.time()
    all_artifacts: List[Path] = []

    def _run_dock6_batch(
        batch_idx: int,
        batch_pairs: List[Tuple[Path, Path]],
    ) -> Tuple[Dict[Path, Dict[str, Any]], List[Path]]:
        local_prefix = prefix if n_batches == 1 else f"{prefix}_b{batch_idx}"
        multi_mol2 = dock6_root / f"{local_prefix}_ligands.mol2"
        input_path = dock6_root / f"{local_prefix}.in"
        output_path = dock6_root / f"{local_prefix}.out"
        ranked_path = dock6_root / f"{local_prefix}_ranked.mol2"

        _build_multi_mol2_for_stage(batch_pairs, multi_mol2, logger)

        local_config_lines = [
            "conformer_search_type rigid",
            "use_internal_energy no",
            f"ligand_atom_file {multi_mol2.name}",
            "limit_max_ligands yes",
            f"max_ligands {len(batch_pairs)}",
            "skip_molecule no",
            "read_mol_solvation no",
            "calculate_rmsd no",
            "use_database_filter no",
            "orient_ligand yes",
            "automated_matching yes",
            "receptor_site_file selected_spheres.sph",
            f"max_orientations {max_orientations}",
            "critical_points no",
            "chemical_matching no",
            "use_ligand_spheres no",
            "bump_filter yes",
            f"bump_grid_prefix {grid_prefix}",
            "max_bumps_anchor 2",
            "max_bumps_growth 2",
            "score_molecules yes",
            "contact_score_primary no",
            "grid_score_primary yes",
            "grid_score_secondary no",
            "grid_score_rep_rad_scale 1.0",
            "grid_score_vdw_scale 1.0",
            "grid_score_es_scale 1.0",
            "grid_lig_efficiency no",
            f"grid_score_grid_prefix {grid_prefix}",
            f"minimize_ligand {minimize_flag}",
            f"simplex_max_iterations {simplex_max_iterations}",
            "simplex_tors_premin_iterations 0",
            f"simplex_max_cycles {simplex_max_cycles}",
            "simplex_score_converge 0.1",
            "simplex_cycle_converge 1.0",
            "simplex_trans_step 1.0",
            "simplex_rot_step 0.1",
            "simplex_tors_step 10.0",
            "simplex_final_min no",
            "simplex_random_seed 0",
            "simplex_restraint_min no",
            "atom_model all",
            f"vdw_defn_file {vdw_defn}",
            f"flex_defn_file {flex_defn}",
            f"flex_drive_file {flex_drive}",
            f"ligand_outfile_prefix {local_prefix}",
            "write_mol_solvation no",
            "write_orientations yes",
            f"num_final_scored_poses {num_final_scored_poses}",
            "num_preclustered_conformers 1",
            "score_threshold 100.0",
            "rank_ligands yes",
            f"max_ranked_ligands {max(len(batch_pairs), 1)}",
        ]
        input_path.write_text("\n".join(local_config_lines) + "\n", encoding="utf-8")

        logger.info(
            "[dock6.cmd] pdb=%s stage=%s batch=%d n_lig=%d exe=%s cwd=%s",
            pdb_id,
            stage_name,
            batch_idx,
            len(batch_pairs),
            dock6_exe,
            str(dock6_root),
        )

        try:
            subprocess.run(
                [dock6_exe, "-i", input_path.name, "-o", output_path.name],
                check=True,
                cwd=str(dock6_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except Exception as exc:
            logger.warning(
                "[dock6.error] pdb=%s stage=%s batch=%d reason=%s",
                pdb_id,
                stage_name,
                batch_idx,
                exc,
            )
            batch_metrics: Dict[Path, Dict[str, Any]] = {}
            for lig, _mol2 in batch_pairs:
                batch_metrics[lig] = {
                    "grid_score": None,
                    "n_poses": 0,
                    "valid": False,
                    "reason": "dock6_run_failed",
                    "mol2_name": mol2_by_pdbqt.get(lig, lig).name,
                }
            return batch_metrics, [input_path, output_path, multi_mol2, ranked_path]

        if not ranked_path.exists() or ranked_path.stat().st_size == 0:
            batch_metrics = {}
            for lig, _mol2 in batch_pairs:
                batch_metrics[lig] = {
                    "grid_score": None,
                    "n_poses": 0,
                    "valid": False,
                    "reason": "no_ranked_output",
                    "mol2_name": mol2_by_pdbqt.get(lig, lig).name,
                }
        else:
            batch_metrics = _parse_dock6_ranked_mol2(ranked_path, batch_pairs, logger)

        return batch_metrics, [input_path, output_path, multi_mol2, ranked_path]

    if n_batches == 1:
        batch_metrics, artifact_paths = _run_dock6_batch(0, batches[0])
        metrics.update(batch_metrics)
        all_artifacts.extend(artifact_paths)
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures: Dict[Any, int] = {}
            for idx, batch_pairs in enumerate(batches):
                futures[pool.submit(_run_dock6_batch, idx, batch_pairs)] = idx

            for fut in as_completed(futures):
                idx = futures[fut]
                try:
                    batch_metrics, artifact_paths = fut.result()
                    metrics.update(batch_metrics)
                    all_artifacts.extend(artifact_paths)
                except Exception as exc:
                    logger.warning(
                        "[dock6.error] pdb=%s stage=%s batch=%d reason=%s",
                        pdb_id,
                        stage_name,
                        idx,
                        exc,
                    )
                    for lig, _mol2 in batches[idx]:
                        if lig in metrics:
                            continue
                        metrics[lig] = {
                            "grid_score": None,
                            "n_poses": 0,
                            "valid": False,
                            "reason": "dock6_run_failed",
                            "mol2_name": mol2_by_pdbqt.get(lig, lig).name,
                        }

    scores = {lig: rec.get("grid_score") for lig, rec in metrics.items()}

    dock_root = (
        Path(output_root)
        if output_root is not None
        else Path(paths.docked_variant_root(variant_for_ph, ph_token))
    )
    if output_root is not None:
        dock_dest = dock_root
    else:
        stage_dir_name = (
            f"dock6_{stage_name}"
            if stage_name
            else f"dock6_{stage_key or 'stage'}"
        )
        dock_dest = dock_root / stage_dir_name
    dock_dest.mkdir(parents=True, exist_ok=True)

    for fname in all_artifacts:
        if fname.exists():
            try:
                shutil.copy2(fname, dock_dest / fname.name)
            except Exception:
                logger.debug("[dock6.copy.skip] src=%s", fname)

    elapsed = time.time() - start_ts
    logger.info(
        "[dock6.done] pdb=%s stage=%s variant=%s ph=%s elapsed_sec=%.2f",
        pdb_id,
        stage_name,
        variant_for_ph or "HOLO",
        ph_token or "base",
        elapsed,
    )
    return scores, metrics


def write_dock6_scores_csv(
    cfg: Dict,
    pdb_id: str,
    dock6_metrics_by_stage: Dict[str, Dict[Path, Dict[str, Any]]],
    *,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    csv_prefix: str = "",
) -> str:
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    dock_dir = Path(paths.docked_variant_root(variant, ph_label))

    dock_dir.mkdir(parents=True, exist_ok=True)

    flat: Dict[str, Dict[str, str]] = {}
    for stage_name, stage_map in dock6_metrics_by_stage.items():
        stage_dict: Dict[str, str] = {}
        for lig, rec in stage_map.items():
            lig_key = rec.get("mol2_name") or Path(lig).name
            score = rec.get("grid_score")
            if isinstance(score, (int, float)):
                stage_dict[lig_key] = f"{score:.2f}"
            else:
                stage_dict[lig_key] = ""
        flat[stage_name] = stage_dict

    summary_path = dock_dir / f"{csv_prefix}dock6_docking_score_summary.csv"
    write_score_summary_to_csv(
        flat,
        output_path=summary_path,
        run_id=str(cfg.get("RUN_ID") or ""),
        variant=(variant or "").strip() or None,
    )

    import csv as _csv

    long_path = dock_dir / f"{csv_prefix}dock6_docking_score_long.csv"
    with open(long_path, "w", newline="", encoding="utf-8") as f:
        writer = _csv.writer(f)
        header = ["run_id"]
        if variant:
            header.append("variant")
        header += ["stage", "ligand", "dock6_grid_score", "dock6_n_poses"]
        writer.writerow(header)
        run_id_val = str(cfg.get("RUN_ID") or "")
        for stage_name, stage_map in dock6_metrics_by_stage.items():
            for lig, rec in stage_map.items():
                row = [run_id_val]
                if variant:
                    row.append(variant)
                row.extend(
                    [
                        stage_name,
                        rec.get("mol2_name") or Path(lig).name,
                        f"{rec.get('grid_score'):.2f}"
                        if isinstance(rec.get("grid_score"), (int, float))
                        else "",
                        rec.get("n_poses", ""),
                    ]
                )
                writer.writerow(row)

    return str(summary_path)
