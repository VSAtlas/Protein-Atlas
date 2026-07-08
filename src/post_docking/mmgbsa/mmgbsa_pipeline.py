# -*- coding: utf-8 -*-
from __future__ import annotations

import os
import re
import csv
import json
import hashlib
import logging
import random
import shutil
import statistics
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
from pathlib import Path

from config.output_paths import output_root
from config.normalize import _to_bool
from post_docking.mmgbsa.prep_for_mmgbsa import prep_mmgbsa_from_sdfs
import post_docking.mmgbsa.protein_prep_mmgbsa as ppm
from post_docking.mmgbsa.protein_prep_mmgbsa import (
    prep_mmgbsa_receptor,
    run_tleap,
    write_leap_for_ligands,
)
from post_docking.mmgbsa.mmgbsa_trajectory import (
    make_mmgbsa_trajectory,
    run_implicit_md,
)
from post_docking.mmgbsa.mmgbsa_pose_selection import (
    _mmgbsa_canonical_ligand_id,
    _mmgbsa_group_pose_sdfs,
    _mmgbsa_pose_group_regexes,
    _mmgbsa_sdf_base_candidates,
    _mmgbsa_select_sdfs,
)
from post_docking.mmgbsa.mmgbsa_methods import write_methods_json
from post_docking.mmgbsa.mmgbsa_explicit import build_and_run_explicit_workflow
from post_docking.mmgbsa.mmgbsa_qc import (
    detect_outliers,
    summarize_values,
    write_frame_qc_artifacts,
)
from post_docking.mmgbsa.mmgbsa_review import (
    default_review_record,
    review_record_path,
    validate_review_record,
    write_review_record,
)
from post_docking.mmgbsa.run_mmgbsa import (
    run_mmgbsa,
    parse_mmpbsa_delta_total,
    write_aggregated_mmpbsa_results,
)
from path_router.path_router import (
    make_paths,
    config_dir as router_config_dir,
    ph_ensemble_dir as router_ph_ensemble_dir,
)
from post_docking.mmgbsa.operation_exceptions import (
    FLOAT_COERCE_ERRORS,
    INT_COERCE_ERRORS,
    OPTIONAL_QC_STEP_ERRORS,
    STAT_AGG_ERRORS,
    TEXT_FILE_READ_ERRORS,
)


def _mmgbsa_find_vina_config(
    run_id: str,
    pdb_id: str,
    stage_dir: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
) -> Optional[Path]:
    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage_dir,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    if not cfg_dir.exists():
        cfg_dir = router_config_dir(
            run_id,
            pdb_id,
            stage_dir,
            variant=variant_token,
            ph_tag=None,
            legacy=legacy_mode,
        )
        if not cfg_dir.exists():
            return None
    txt_files = sorted(cfg_dir.glob("*.txt"))
    return txt_files[0] if txt_files else None


def _mmgbsa_parse_vina_config(
    cfg_path: Path,
) -> tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]]]:
    center_vals: dict[str, float] = {}
    size_vals: dict[str, float] = {}
    try:
        lines = cfg_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except TEXT_FILE_READ_ERRORS as exc:
        logging.getLogger("mmgbsa.pipeline").debug(
            "[mmgbsa.vina_config] action=read_failed path=%s err=%s",
            cfg_path,
            exc,
            exc_info=True,
        )
        return None, None

    for line in lines:
        if "=" not in line:
            continue
        key, raw_val = line.split("=", 1)
        key = key.strip().lower()
        val_text = raw_val.strip()
        try:
            val = float(val_text)
        except FLOAT_COERCE_ERRORS:
            continue
        if key in {"center_x", "center_y", "center_z"}:
            center_vals[key] = val
        elif key in {"size_x", "size_y", "size_z"}:
            size_vals[key] = val

    if len(center_vals) == 3:
        center = (
            center_vals["center_x"],
            center_vals["center_y"],
            center_vals["center_z"],
        )
    else:
        center = None

    if len(size_vals) == 3:
        size = (size_vals["size_x"], size_vals["size_y"], size_vals["size_z"])
    else:
        size = None

    return center, size


def _mmgbsa_resolve_center_radius(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    stage_dir: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
) -> tuple[Optional[Tuple[float, float, float]], float, Optional[Path]]:
    fallback_radius = float(cfg.get("MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK", 6.0))
    cfg_path = _mmgbsa_find_vina_config(
        run_id, pdb_id, stage_dir, variant_token, ph_label, legacy_mode
    )
    if not cfg_path:
        return None, fallback_radius, None

    center, size = _mmgbsa_parse_vina_config(cfg_path)
    if center is None:
        return None, fallback_radius, cfg_path

    radius = fallback_radius
    if size is not None:
        radius = max(size) / 2.0
    return center, radius, cfg_path


def _mmgbsa_ligand_centroid(
    sdfs: Sequence[Path],
) -> Optional[Tuple[float, float, float]]:
    if not sdfs:
        return None
    try:
        from rdkit import Chem
    except ImportError:
        return None
    coords: list[tuple[float, float, float]] = []
    for sdf_path in sdfs:
        coords.extend(_mmgbsa_sdf_heavy_atom_coords(Chem, sdf_path))
    if not coords:
        return None
    n_atoms = float(len(coords))
    return (
        sum(x for x, _, _ in coords) / n_atoms,
        sum(y for _, y, _ in coords) / n_atoms,
        sum(z for _, _, z in coords) / n_atoms,
    )


def _mmgbsa_sdf_heavy_atom_coords(
    chem_module: Any,
    sdf_path: Path,
) -> list[tuple[float, float, float]]:
    supplier = chem_module.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
    mol = supplier[0] if supplier and len(supplier) else None
    if mol is None or mol.GetNumConformers() == 0:
        return []
    conf = mol.GetConformer()
    coords: list[tuple[float, float, float]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append((float(pos.x), float(pos.y), float(pos.z)))
    return coords


def _mmgbsa_resolve_receptor_pdb(
    cfg: Mapping[str, Any],
    pdb_id: str,
    pdb_file: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
) -> Optional[Path]:
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=pdb_file)
    candidates: List[Path] = []

    if ph_label:
        ensemble_dir = router_ph_ensemble_dir(
            pdb_id, variant=variant_token, legacy=legacy_mode
        )
        tag = f"{pdb_id}_{ph_label}"
        candidates.append(ensemble_dir / f"{tag}.withH.pdb")
        candidates.append(ensemble_dir / f"{tag}.pdb")

    candidates.append(paths.receptor_cleaned_pdb(variant_token))
    candidates.append(paths.receptor_dir(variant_token) / f"{pdb_id}.pdb")

    for cand in candidates:
        if cand.exists():
            return cand
    return None

def _mmgbsa_index_sdfs(sdfs: List[Path], stage_dir_label: str) -> Dict[str, Path]:
    mapping: Dict[str, Path] = {}
    for sdf in sdfs:
        for key in _mmgbsa_sdf_base_candidates(sdf.stem, stage_dir_label):
            if key and key not in mapping:
                mapping[key] = sdf
    return mapping


def _mmgbsa_parse_delta_total(csv_path: Path) -> Optional[float]:
    if not csv_path.exists():
        return None
    try:
        lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except TEXT_FILE_READ_ERRORS as exc:
        logging.getLogger("mmgbsa.pipeline").debug(
            "[mmgbsa.delta_total] action=csv_read_failed path=%s err=%s",
            csv_path,
            exc,
            exc_info=True,
        )
        return None

    in_delta = False
    header = None
    delta_idx = None
    frame_values: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("DELTA Energy Terms"):
            in_delta = True
            header = None
            delta_idx = None
            continue
        if not in_delta:
            continue
        if header is None:
            header = [t.strip() for t in stripped.split(",")]
            try:
                delta_idx = header.index("DELTA TOTAL")
            except ValueError:
                return None
            continue
        values = [t.strip() for t in stripped.split(",")]
        if not header or delta_idx is None:
            return None
        idx = delta_idx
        if idx >= len(values):
            continue
        try:
            frame_values.append(float(values[idx]))
        except FLOAT_COERCE_ERRORS:
            continue
    if not frame_values:
        return None
    return float(frame_values[0])


def _mmgbsa_parse_delta_frames(csv_path: Path) -> Dict[str, object]:
    result = {
        "frame_values": [],
        "n_frames": 0,
        "frame_mean": None,
        "frame_median": None,
        "frame_sd": None,
        "frame_min": None,
        "frame_max": None,
    }
    if not csv_path.exists():
        return result
    try:
        lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except TEXT_FILE_READ_ERRORS as exc:
        logging.getLogger("mmgbsa.pipeline").debug(
            "[mmgbsa.delta_frames] action=csv_read_failed path=%s err=%s",
            csv_path,
            exc,
            exc_info=True,
        )
        return result

    in_delta = False
    header = None
    delta_idx = None
    frames: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("DELTA Energy Terms"):
            in_delta = True
            header = None
            delta_idx = None
            continue
        if not in_delta:
            continue
        if header is None:
            header = [t.strip() for t in stripped.split(",")]
            try:
                delta_idx = header.index("DELTA TOTAL")
            except ValueError:
                break
            continue
        values = [t.strip() for t in stripped.split(",")]
        if not header or delta_idx is None:
            break
        idx = delta_idx
        if idx >= len(values):
            continue
        try:
            frames.append(float(values[idx]))
        except FLOAT_COERCE_ERRORS:
            continue

    if not frames:
        return result

    result["frame_values"] = frames
    result["n_frames"] = len(frames)
    try:
        result["frame_mean"] = float(statistics.mean(frames))
    except STAT_AGG_ERRORS:
        result["frame_mean"] = None
    try:
        result["frame_median"] = float(statistics.median(frames))
    except STAT_AGG_ERRORS:
        result["frame_median"] = None
    try:
        result["frame_sd"] = float(statistics.pstdev(frames))
    except STAT_AGG_ERRORS:
        result["frame_sd"] = None
    try:
        result["frame_min"] = float(min(frames))
        result["frame_max"] = float(max(frames))
    except STAT_AGG_ERRORS:
        pass
    return result


def _mmgbsa_write_replicate_summary(
    summary_path: Path, rows: List[Dict[str, object]], mean: float, sd: float
) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["replicate", "seed", "score", "ok", "notes"]
    tmp_path = summary_path.with_suffix(summary_path.suffix + ".part")
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        writer.writerow({"replicate": "mean", "score": f"{mean:.6g}", "ok": ""})
        writer.writerow({"replicate": "stddev", "score": f"{sd:.6g}", "ok": ""})
    os.replace(tmp_path, summary_path)


def _mmgbsa_five_replicate_runner(
    topo: Dict[str, object],
    out_dir: Path,
    cfg: Mapping[str, Any],
    force: bool,
    md_enabled: bool,
    stage_dir: str,
    ligand_stem: str,
    pdb_id: str,
    variant_dir: str,
    ph_label: str,
    run_id: str,
    logger: logging.Logger,
) -> Dict[str, object]:
    seeds = _mmgbsa_five_rep_seeds(
        run_id, pdb_id, variant_dir, ph_label, stage_dir, ligand_stem
    )
    out_dat_name = str(
        cfg.get("MMGBSA_MMPBSA_OUT_DAT", "FINAL_RESULTS_MMPBSA.dat")
        or "FINAL_RESULTS_MMPBSA.dat"
    )
    out_csv_name = str(
        cfg.get("MMGBSA_MMPBSA_OUT_CSV", "FINAL_RESULTS_MMPBSA.csv")
        or "FINAL_RESULTS_MMPBSA.csv"
    )
    agg_dat = out_dir / out_dat_name
    agg_csv = out_dir / out_csv_name

    if not force and agg_dat.exists() and agg_csv.exists():
        existing_delta = parse_mmpbsa_delta_total(agg_csv)
        if existing_delta is not None:
            return {
                "ok": True,
                "delta_total": existing_delta,
                "results_csv": str(agg_csv),
                "results_dat": str(agg_dat),
                "replicates_ok": 5,
                "replicates_total": 5,
                "summary_path": str(out_dir / "mmgbsa_replicates_summary.csv"),
                "notes": "cached",
            }

    rep_rows: List[Dict[str, object]] = []
    scores: List[float] = []
    summary_path = out_dir / "mmgbsa_replicates_summary.csv"

    if not md_enabled:
        logger.warning(
            "[mmgbsa.replicate] MD_FIVE_REPLICATE enabled but MD disabled; running single MMGBSA and cloning results"
        )
        rep_dir = out_dir / "rep1"
        traj_result = make_mmgbsa_trajectory(
            complex_prmtop=topo["complex_prmtop"],
            complex_inpcrd=topo["complex_inpcrd"],
            out_dir=str(rep_dir),
            cfg=cfg,
            force=force,
            run_cpptraj=_to_bool(cfg.get("MMGBSA_CPPTRAJ_RUN", True)),
        )
        traj_path = Path(traj_result.get("trajout_path") or "")
        mmpbsa_dir = rep_dir / "mmpbsa"
        mmpbsa_dir.mkdir(parents=True, exist_ok=True)
        mm_res = run_mmgbsa(
            complex_prmtop=topo["complex_prmtop"],
            receptor_prmtop=topo["receptor_prmtop"],
            ligand_prmtop=topo["ligand_prmtop"],
            trajectory_path=str(traj_path),
            work_dir=str(mmpbsa_dir),
            cfg=cfg,
            force=force,
            run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
        )
        rep_csv = Path(mm_res.get("out_csv", ""))
        score = parse_mmpbsa_delta_total(rep_csv)
        rep_ok = score is not None and rep_csv.exists()
        rep_rows.append(
            {
                "replicate": 1,
                "seed": seeds[0],
                "score": score if score is not None else "",
                "ok": rep_ok,
                "notes": "",
            }
        )
        if not rep_ok:
            _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=0.0, sd=0.0)
            return {
                "ok": False,
                "notes": "rep1_failed",
                "replicates_ok": 0,
                "replicates_total": 5,
            }
        scores.append(float(score))
        # clone outputs for other reps
        for idx in range(2, 6):
            clone_dir = out_dir / f"rep{idx}" / "mmpbsa"
            clone_dir.mkdir(parents=True, exist_ok=True)
            for src in (
                rep_csv,
                Path(mm_res.get("out_dat", "")),
                Path(mm_res.get("log_path", "")),
                Path(mm_res.get("input_path", "")),
            ):
                if not src:
                    continue
                if not src.exists():
                    continue
                dest = clone_dir / src.name
                if dest.exists() and not force:
                    continue
                shutil.copy2(src, dest)
            rep_rows.append(
                {
                    "replicate": idx,
                    "seed": seeds[idx - 1],
                    "score": score,
                    "ok": True,
                    "notes": "cloned",
                }
            )
            scores.append(float(score))
    else:
        for rep_idx, seed in enumerate(seeds, start=1):
            rep_dir = out_dir / f"rep{rep_idx}"
            md_result = run_implicit_md(
                complex_prmtop=topo["complex_prmtop"],
                complex_inpcrd=topo["complex_inpcrd"],
                out_dir=str(out_dir),
                cfg=cfg,
                replicate_index=rep_idx,
                seed=int(seed),
                force=force,
                run=True,
            )
            traj_path = Path(md_result.get("traj_path") or "")
            rep_ok = False
            score = None
            notes = ""
            if not md_result.get("ok"):
                notes = "md_failed"
            elif not traj_path.exists() or traj_path.stat().st_size == 0:
                notes = "missing_trajectory"
            else:
                mmpbsa_dir = rep_dir / "mmpbsa"
                mmpbsa_dir.mkdir(parents=True, exist_ok=True)
                mm_res = run_mmgbsa(
                    complex_prmtop=topo["complex_prmtop"],
                    receptor_prmtop=topo["receptor_prmtop"],
                    ligand_prmtop=topo["ligand_prmtop"],
                    trajectory_path=str(traj_path),
                    work_dir=str(mmpbsa_dir),
                    cfg=cfg,
                    force=force,
                    run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
                )
                rep_csv = Path(mm_res.get("out_csv", ""))
                score = parse_mmpbsa_delta_total(rep_csv)
                if score is not None and rep_csv.exists():
                    rep_ok = True
                else:
                    notes = "missing_outputs"
            rep_rows.append(
                {
                    "replicate": rep_idx,
                    "seed": seed,
                    "score": score if score is not None else "",
                    "ok": rep_ok,
                    "notes": notes,
                }
            )
            if rep_ok and score is not None:
                scores.append(float(score))

    if len(scores) != 5 or any(not r.get("ok") for r in rep_rows):
        _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=0.0, sd=0.0)
        return {
            "ok": False,
            "notes": "replicate_failed",
            "replicates_ok": len(scores),
            "replicates_total": 5,
            "summary_path": str(summary_path),
        }

    mean_val = float(statistics.mean(scores))
    try:
        sd_val = float(statistics.pstdev(scores))
    except STAT_AGG_ERRORS as exc:
        logging.getLogger("mmgbsa.pipeline").debug(
            "[mmgbsa.replicates] action=pstdev_fallback err=%s n_scores=%s",
            exc,
            len(scores),
            exc_info=True,
        )
        sd_val = 0.0

    _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=mean_val, sd=sd_val)
    agg_paths = write_aggregated_mmpbsa_results(
        work_dir=out_dir, mean_score=mean_val, std_score=sd_val, cfg=cfg, force=force
    )

    logger.info(
        "[mmgbsa.aggregate] n=5 mean=%.6g sd=%.6g out_csv=%s out_dat=%s",
        mean_val,
        sd_val,
        agg_paths.get("out_csv"),
        agg_paths.get("out_dat"),
    )

    return {
        "ok": True,
        "delta_total": mean_val,
        "results_csv": agg_paths.get("out_csv", ""),
        "results_dat": agg_paths.get("out_dat", ""),
        "replicates_ok": 5,
        "replicates_total": 5,
        "summary_path": str(summary_path),
        "notes": "",
        "frame_mean": mean_val,
        "frame_sd": sd_val,
    }


def _mmgbsa_trimmed_mean(values: List[float], trim_fraction: float = 0.1) -> float:
    if not values:
        raise ValueError("no values for trimmed mean")
    if trim_fraction <= 0:
        return float(statistics.mean(values))
    sorted_vals = sorted(values)
    trim_n = int(len(sorted_vals) * trim_fraction)
    max_trim = max(0, (len(sorted_vals) - 1) // 2)
    trim_n = max(0, min(trim_n, max_trim))
    if trim_n > 0:
        trimmed = sorted_vals[trim_n:-trim_n] or sorted_vals
    else:
        trimmed = sorted_vals
    return float(statistics.mean(trimmed))


def _mmgbsa_apply_agg(
    values: List[float],
    method: str,
    logger: Optional[logging.Logger] = None,
    context: str = "mmgbsa",
) -> Optional[float]:
    if not values:
        return None
    method_norm = str(method or "mean").strip().lower()
    try:
        if method_norm == "median":
            return float(statistics.median(values))
        if method_norm == "min":
            return float(min(values))
        if method_norm == "trimmed_mean":
            return _mmgbsa_trimmed_mean(values, trim_fraction=0.1)
        return float(statistics.mean(values))
    except STAT_AGG_ERRORS as exc:
        if logger:
            logger.warning(
                "[mmgbsa.pipeline] action=aggregate_failed context=%s method=%s err=%s",
                context,
                method_norm,
                exc,
                exc_info=True,
            )
    return None


def _mmgbsa_frame_aggregate(
    csv_path: Path, cfg: Mapping[str, Any], logger: logging.Logger
) -> Dict[str, object]:
    parsed = _mmgbsa_parse_delta_frames(csv_path)
    frames: List[float] = list(parsed.get("frame_values") or [])
    total_frames = len(frames)
    try:
        frame_limit = int(cfg.get("MMGBSA_FRAME_LIMIT", 0))
    except INT_COERCE_ERRORS:
        frame_limit = 0
    if frame_limit > 0 and frames:
        frames = frames[:frame_limit]
    frames_used = len(frames)
    frame_method = str(cfg.get("MMGBSA_FRAME_AGG", "mean") or "mean").strip().lower()

    agg_delta = _mmgbsa_apply_agg(frames, frame_method, logger=logger, context="frame")
    qc_payload: Dict[str, object] = {}
    if _to_bool(cfg.get("MMGBSA_QC_ENABLED", True)) and frames:
        try:
            qc_payload = write_frame_qc_artifacts(
                csv_path.parent / "qc",
                frames,
                cfg,
                label=csv_path.parent.name or "mmgbsa",
            )
        except OPTIONAL_QC_STEP_ERRORS as exc:
            logger.warning(
                "[mmgbsa.qc] action=skip reason=frame_qc_failed path=%s err=%s",
                csv_path,
                exc,
                exc_info=True,
            )

    subset_mean = None
    subset_median = None
    subset_sd = None
    try:
        subset_mean = float(statistics.mean(frames))
    except STAT_AGG_ERRORS:
        subset_mean = None
    try:
        subset_median = float(statistics.median(frames))
    except STAT_AGG_ERRORS:
        subset_median = None
    try:
        subset_sd = float(statistics.pstdev(frames))
    except STAT_AGG_ERRORS:
        subset_sd = None

    if frames_used:
        logger.info(
            "[mmgbsa.frame] csv=%s frames_total=%d frames_used=%d limit=%d method=%s delta=%s mean=%s sd=%s",
            csv_path,
            total_frames,
            frames_used,
            frame_limit,
            frame_method,
            f"{agg_delta:.6g}" if agg_delta is not None else "",
            f"{subset_mean:.6g}" if subset_mean is not None else "",
            f"{subset_sd:.6g}" if subset_sd is not None else "",
        )
    else:
        logger.warning(
            "[mmgbsa.frame] csv=%s frames_total=%d frames_used=%d limit=%d method=%s reason=missing_frames",
            csv_path,
            total_frames,
            frames_used,
            frame_limit,
            frame_method,
        )

    return {
        "delta_total": agg_delta,
        "frame_agg_method": frame_method,
        "frames_used": frames_used,
        "frames_total": total_frames,
        "frame_mean": subset_mean,
        "frame_median": subset_median,
        "frame_sd": subset_sd,
        "frame_values": frames,
        "frame_limit": frame_limit,
        "frame_sem": (
            (qc_payload.get("qc") or {}).get("sem")
            if isinstance(qc_payload.get("qc"), dict)
            else None
        ),
        "frame_ci95": (
            (qc_payload.get("qc") or {}).get("ci95")
            if isinstance(qc_payload.get("qc"), dict)
            else None
        ),
        "qc_pass": (
            (qc_payload.get("qc") or {}).get("qc_pass")
            if isinstance(qc_payload.get("qc"), dict)
            else ""
        ),
        "qc_json": str(csv_path.parent / "qc" / "mmgbsa_frame_qc.json") if qc_payload else "",
        "notes": "" if agg_delta is not None else "missing_delta",
    }


def _mmgbsa_summary_frame_fields(frame_meta: Mapping[str, object]) -> Dict[str, object]:
    return {
        "frame_agg_method": frame_meta.get("frame_agg_method", ""),
        "frames_used": frame_meta.get("frames_used", ""),
        "frame_mean": frame_meta.get("frame_mean", ""),
        "frame_sd": frame_meta.get("frame_sd", ""),
        "frames_total": frame_meta.get("frames_total", ""),
        "frame_median": frame_meta.get("frame_median", ""),
        "frame_sem": frame_meta.get("frame_sem", ""),
        "frame_ci95": frame_meta.get("frame_ci95", ""),
        "qc_pass": frame_meta.get("qc_pass", ""),
        "qc_json": frame_meta.get("qc_json", ""),
    }


def _mmgbsa_append_summary(summary_path: Path, row: Dict[str, object]) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    exists = summary_path.exists() and summary_path.stat().st_size > 0
    fieldnames = [
        "stage_dir",
        "ligand_stem",
        "delta_total",
        "results_csv",
        "ok",
        "notes",
        "frame_agg_method",
        "frames_used",
        "frames_total",
        "frame_mean",
        "frame_sd",
        "frame_median",
        "frame_sem",
        "frame_ci95",
        "qc_pass",
        "qc_json",
        "rep_mean",
        "rep_sd",
        "rep_sem",
        "rep_ci95",
        "rep_outliers",
        "rep_agg_method",
        "replicates_ok",
        "replicates_total",
    ]
    with summary_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, extrasaction="ignore", restval=""
        )
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def _mmgbsa_write_pose_aggregate(
    summary_path: Path,
    out_path: Path,
    stage_dir_label: str,
    pose_group_regexes: List[re.Pattern],
    agg_method: str,
    logger: logging.Logger,
) -> None:
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        logger.info(
            "[mmgbsa.pipeline] aggregate=skip reason=missing_summary path=%s",
            summary_path,
        )
        return

    groups: Dict[Tuple[str, str], Dict[str, object]] = {}
    try:
        with summary_path.open(
            "r", encoding="utf-8", errors="ignore", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                stage_dir = str(row.get("stage_dir") or stage_dir_label or "").strip()
                ligand_stem = str(row.get("ligand_stem") or "").strip()
                if not ligand_stem:
                    continue
                ligand_id = _mmgbsa_canonical_ligand_id(
                    ligand_stem, stage_dir, pose_group_regexes, logger
                )
                key = (stage_dir, ligand_id)
                entry = groups.setdefault(key, {"total": 0, "ok": []})
                entry["total"] = int(entry.get("total", 0)) + 1

                ok_val = _to_bool(row.get("ok", False))
                delta_val = None
                if ok_val:
                    raw = str(row.get("delta_total", "")).strip()
                    if raw:
                        try:
                            delta_val = float(raw)
                        except Exception:
                            delta_val = None
                if ok_val and delta_val is not None:
                    entry["ok"].append((delta_val, ligand_stem))
    except Exception as exc:
        logger.warning(
            "[mmgbsa.pipeline] aggregate=skip reason=read_error path=%s err=%s",
            summary_path,
            exc,
        )
        return

    agg_method_norm = str(agg_method or "min").strip().lower()
    if agg_method_norm not in ("min", "median"):
        logger.warning(
            "[mmgbsa.pipeline] aggregate=unknown_method method=%s fallback=min",
            agg_method_norm,
        )
        agg_method_norm = "min"

    rows_out: List[Dict[str, object]] = []
    for (stage_dir, ligand_id), entry in sorted(groups.items()):
        ok_list = list(entry.get("ok", []))
        n_total = int(entry.get("total", 0))
        n_ok = len(ok_list)
        agg_delta = ""
        best_pose = ""
        included: List[str] = []

        if n_ok:
            deltas = [val for val, _ in ok_list]
            if agg_method_norm == "median":
                agg_val = float(statistics.median(deltas))
                best = min(ok_list, key=lambda x: (abs(x[0] - agg_val), x[0], x[1]))
            else:
                agg_val = min(deltas)
                best = min(ok_list, key=lambda x: (x[0], x[1]))
            agg_delta = f"{agg_val:.6g}"
            best_pose = best[1]
            included = [stem for _, stem in ok_list]

        rows_out.append(
            {
                "stage_dir": stage_dir,
                "ligand_id": ligand_id,
                "n_poses_total": n_total,
                "n_poses_ok": n_ok,
                "agg_method": agg_method_norm,
                "agg_delta": agg_delta,
                "best_pose_stem": best_pose,
                "pose_stems_included": ";".join(included),
            }
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    fieldnames = [
        "stage_dir",
        "ligand_id",
        "n_poses_total",
        "n_poses_ok",
        "agg_method",
        "agg_delta",
        "best_pose_stem",
        "pose_stems_included",
    ]
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)
    os.replace(tmp_path, out_path)
    logger.info(
        "[mmgbsa.pipeline] aggregate=pose_delta path=%s rows=%d method=%s",
        out_path,
        len(rows_out),
        agg_method_norm,
    )


def _mmgbsa_md_seed_list(
    cfg: Mapping[str, Any], n_reps: int, logger: logging.Logger
) -> List[int]:
    try:
        base_seed = int(cfg.get("MMGBSA_MD_BASE_SEED", 12345))
    except Exception:
        base_seed = 12345

    mode = (
        str(cfg.get("MMGBSA_MD_SEED_MODE", "increment") or "increment").strip().lower()
    )
    if n_reps < 1:
        n_reps = 1

    seeds: List[int] = []
    if mode == "random":
        rng = random.Random(base_seed)
        for _ in range(n_reps):
            seeds.append(rng.randint(1, 2_147_483_647))
    else:
        for idx in range(n_reps):
            seeds.append(max(1, base_seed + idx))

    logger.info(
        "[mmgbsa.pipeline] md_seeds mode=%s base=%s reps=%d",
        mode,
        base_seed,
        n_reps,
    )
    return seeds


def _mmgbsa_five_rep_seeds(
    run_id: str, pdb_id: str, variant: str, ph_label: str, stage_dir: str, ligand: str
) -> List[int]:
    context = f"{run_id}|{pdb_id}|{variant}|{ph_label}|{stage_dir}|{ligand}"
    base_seed = int(hashlib.md5(context.encode("utf-8")).hexdigest()[:8], 16)
    seeds: List[int] = []
    for idx in range(5):
        seeds.append(max(1, base_seed + idx * 10007))
    return seeds


def _mmgbsa_md_aggregate(deltas: List[float], method: str) -> Optional[float]:
    return _mmgbsa_apply_agg(
        deltas, method, logger=logging.getLogger("mmgbsa.pipeline"), context="replicate"
    )


_MMGBSA_DEPRECATED_MD_KEYS = ("MMGBSA_MD_RUN", "MMGBSA_TRAJ_MODE")
_MMGBSA_DEPRECATED_MD_WARNED = False


def _mmgbsa_effective_md_config(
    cfg: Mapping[str, Any], logger: Optional[logging.Logger]
) -> Dict[str, object]:
    global _MMGBSA_DEPRECATED_MD_WARNED
    deprecated_keys = []
    for key in _MMGBSA_DEPRECATED_MD_KEYS:
        if key in cfg and str(cfg.get(key)).strip() != "":
            deprecated_keys.append(key)

    md_enabled_raw = cfg.get("MMGBSA_MD_ENABLED", None)
    md_enabled_set = "MMGBSA_MD_ENABLED" in cfg and str(md_enabled_raw).strip() != ""
    if md_enabled_set:
        md_enabled = _to_bool(md_enabled_raw)
    else:
        traj_mode = str(cfg.get("MMGBSA_TRAJ_MODE", "") or "").strip().upper()
        md_enabled = traj_mode == "IMPLICIT_MD" or _to_bool(
            cfg.get("MMGBSA_MD_RUN", False)
        )

    if deprecated_keys and not _MMGBSA_DEPRECATED_MD_WARNED and logger is not None:
        logger.warning(
            "[mmgbsa.pipeline] deprecated_keys=%s msg=Deprecated MMGBSA keys detected; please use MMGBSA_MD_ENABLED",
            ",".join(sorted(deprecated_keys)),
        )
        _MMGBSA_DEPRECATED_MD_WARNED = True

    solvent_model = str(cfg.get("MMGBSA_MD_SOLVENT_MODEL", "implicit") or "implicit")
    effective_mode = "ONEFRAME"
    if md_enabled:
        effective_mode = "EXPLICIT_MD" if solvent_model.lower() == "explicit" else "IMPLICIT_MD"
    return {
        "md_enabled": md_enabled,
        "effective_mode": effective_mode,
        "deprecated_keys": deprecated_keys,
    }


def _mmgbsa_apply_env_overrides(
    cfg: Mapping[str, Any], keys: Sequence[str]
) -> Dict[str, Any]:
    updated = dict(cfg)
    for key in keys:
        value = os.environ.get(key)
        if value is None or str(value).strip() == "":
            continue
        updated[key] = value
    return updated


def _mmgbsa_split_cfg_list(value: Any) -> List[str]:
    if value is None:
        return []
    raw = str(value).strip()
    if not raw:
        return []
    return [part.strip() for part in re.split(r"[,;\n]+", raw) if part.strip()]


def _mmgbsa_protocol(cfg: Mapping[str, Any]) -> str:
    token = str(cfg.get("MMGBSA_PROTOCOL", "screening") or "screening").strip().lower()
    if token in {"publication", "publish"}:
        return "production"
    if token in {"screen", "triage", "fast"}:
        return "screening"
    return token


def _mmgbsa_validate_protocol(
    cfg: Mapping[str, Any], logger: logging.Logger, context: str
) -> None:
    protocol = _mmgbsa_protocol(cfg)
    if protocol != "production":
        return
    strict = _to_bool(cfg.get("MMGBSA_PRODUCTION_STRICT", True))
    problems: List[str] = []
    if not _to_bool(cfg.get("MMGBSA_MD_ENABLED", False)):
        problems.append("MMGBSA_MD_ENABLED must be true")
    try:
        reps = int(cfg.get("MMGBSA_MD_NREPLICATES", 1))
    except Exception:
        reps = 1
    if reps < 3:
        problems.append("MMGBSA_MD_NREPLICATES must be at least 3")
    if not _to_bool(cfg.get("MMGBSA_MMPBSA_USE_TRAJ_FRAMES", True)):
        problems.append("MMGBSA_MMPBSA_USE_TRAJ_FRAMES must be true")
    chem_required = _to_bool(cfg.get("MMGBSA_CHEMISTRY_REVIEW_REQUIRED", True))
    chem_status = (
        str(cfg.get("MMGBSA_CHEMISTRY_REVIEW_STATUS", "unreviewed") or "unreviewed")
        .strip()
        .lower()
    )
    if chem_required and chem_status not in {"reviewed", "approved", "ok"}:
        problems.append(
            "chemistry review required; set MMGBSA_CHEMISTRY_REVIEW_STATUS=reviewed "
            "after protonation/tautomer/charge/water/metal review"
        )
    if _to_bool(cfg.get("MMGBSA_CHEMISTRY_SUBGATES_REQUIRED", False)):
        for key in (
            "MMGBSA_PROTONATION_REVIEW_STATUS",
            "MMGBSA_TAUTOMER_REVIEW_STATUS",
            "MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS",
            "MMGBSA_NET_CHARGE_REVIEW_STATUS",
            "MMGBSA_PARAMETER_REVIEW_STATUS",
        ):
            status = str(cfg.get(key, "unreviewed") or "unreviewed").strip().lower()
            if status not in {"reviewed", "approved", "ok"}:
                problems.append(f"{key} must be reviewed")
    if not _to_bool(cfg.get("MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE", False)):
        problems.append(
            "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE must be true for production "
            "after confirming the pose came from curated SDF/MOL2/PDB chemistry, "
            "not PDBQT reconstruction"
        )
    if problems:
        msg = "; ".join(problems)
        logger.error(
            "[mmgbsa.protocol] protocol=production context=%s status=invalid problems=%s",
            context,
            msg,
        )
        if strict:
            raise ValueError(f"MMGBSA production protocol invalid: {msg}")
    logger.info("[mmgbsa.protocol] protocol=production context=%s status=ok", context)


def _mmgbsa_review_status(cfg: Mapping[str, Any], key: str) -> str:
    return str(cfg.get(key, "unreviewed") or "unreviewed").strip().lower()


def _mmgbsa_review_section(cfg: Mapping[str, Any], key: str, source: str) -> Dict[str, object]:
    status = _mmgbsa_review_status(cfg, key)
    reviewed = status in {"reviewed", "approved", "ok"}
    return {
        "reviewed": reviewed,
        "status": status,
        "source": source,
        "notes": "",
    }


def _mmgbsa_ph_value(ph_label: str) -> Optional[float]:
    match = re.search(r"(\d+(?:[_.]\d+)?)", ph_label)
    if not match:
        return None
    try:
        return float(match.group(1).replace("_", "."))
    except Exception:
        return None


def _mmgbsa_write_review_record(
    *,
    cfg: Mapping[str, Any],
    mmgbsa_dir: Path,
    pdb_id: str,
    ligand_stem: str,
    variant_dir: str,
    ph_label: str,
    receptor_result: Mapping[str, object],
    prep_result: Mapping[str, object] | None,
    logger: logging.Logger,
) -> Dict[str, object]:
    review_root = mmgbsa_dir / "review_records"
    record = default_review_record(pdb_id, ligand_stem)
    source = "atlas_config_and_pipeline_metadata"
    record["ligand_protonation"] = _mmgbsa_review_section(
        cfg, "MMGBSA_PROTONATION_REVIEW_STATUS", source
    )
    record["tautomer"] = _mmgbsa_review_section(
        cfg, "MMGBSA_TAUTOMER_REVIEW_STATUS", source
    )
    record["stereochemistry"] = _mmgbsa_review_section(
        cfg, "MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS", source
    )
    record["parameter_review"] = _mmgbsa_review_section(
        cfg, "MMGBSA_PARAMETER_REVIEW_STATUS", source
    )
    if prep_result:
        record["ligand_state_provenance"] = prep_result.get(
            "ligand_state_provenance", {}
        )
        record["suspicious_chemistry"] = prep_result.get("suspicious_chemistry", {})
    record["input_chemistry_authority"] = {
        **_mmgbsa_review_section(cfg, "MMGBSA_CHEMISTRY_REVIEW_STATUS", source),
        "authoritative": _to_bool(cfg.get("MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE", False)),
        "chemistry_source": str(
            cfg.get("MMGBSA_INPUT_CHEMISTRY_SOURCE", "unreviewed_sdf")
            or "unreviewed_sdf"
        ),
    }
    record["net_charge"] = {
        **_mmgbsa_review_section(cfg, "MMGBSA_NET_CHARGE_REVIEW_STATUS", source),
        "value": prep_result.get("net_charge_used") if prep_result else None,
        "charge_method": prep_result.get("charge_method") if prep_result else "",
    }
    water_policy = str(cfg.get("MMGBSA_WATER_POLICY", "ACTIVE_SITE") or "ACTIVE_SITE")
    record["water_policy"] = {
        **_mmgbsa_review_section(cfg, "MMGBSA_CHEMISTRY_REVIEW_STATUS", source),
        "policy": water_policy,
        "retained_waters": [],
        "water_residues_total": receptor_result.get("water_residues_total", 0),
        "water_residues_kept": receptor_result.get("water_residues_kept", 0),
        "water_residues_removed": receptor_result.get("water_residues_removed", 0),
    }
    record["metal_policy"] = {
        **_mmgbsa_review_section(cfg, "MMGBSA_CHEMISTRY_REVIEW_STATUS", source),
        "policy": "keep" if _to_bool(cfg.get("MMGBSA_KEEP_METALS", False)) else "strip",
        "retained_metals": _mmgbsa_split_cfg_list(cfg.get("MMGBSA_METAL_RETAIN_TOKENS")),
        "metal_lines_removed": receptor_result.get("metal_lines_removed", 0),
    }
    record["apo_holo"] = {
        **_mmgbsa_review_section(cfg, "MMGBSA_CHEMISTRY_REVIEW_STATUS", source),
        "value": "HOLO" if str(variant_dir).upper() == "HOLO" else "APO",
    }
    record["ph"] = {
        **_mmgbsa_review_section(cfg, "MMGBSA_CHEMISTRY_REVIEW_STATUS", source),
        "value": _mmgbsa_ph_value(ph_label),
        "label": ph_label,
    }
    record["notes"] = str(cfg.get("MMGBSA_CHEMISTRY_REVIEW_NOTES", "") or "")
    path = review_record_path(review_root, pdb_id, ligand_stem)
    write_review_record(path, record)
    validation = validate_review_record(
        record, production=_mmgbsa_protocol(cfg) == "production"
    )
    logger.info(
        "[mmgbsa.review] pdb=%s ligand=%s production_ready=%s problems=%d path=%s",
        pdb_id,
        ligand_stem,
        validation.get("production_ready"),
        len(validation.get("problems", [])),
        path,
    )
    return {"path": str(path), "record": record, "validation": validation}


def _mmgbsa_frame_qc_ok(cfg: Mapping[str, Any], frame_meta: Mapping[str, object]) -> bool:
    if _mmgbsa_protocol(cfg) != "production":
        return True
    if not _to_bool(cfg.get("MMGBSA_QC_ENABLED", True)):
        return True
    return frame_meta.get("qc_pass") is True


def _maybe_run_mmgbsa_for_pdb(
    cfg: Mapping[str, Any],
    pdb_file: str,
    pdb_id: str,
    variant_token: Optional[str],
    run_id: str,
    test_mode: str,
    legacy_mode: bool,
) -> None:
    logger = logging.getLogger("mmgbsa.pipeline")
    cfg = _mmgbsa_apply_env_overrides(
        cfg,
        [
            "MMGBSA_MAX_LIGANDS",
            "MMGBSA_RERANKED_TOP_PCT",
            "MMGBSA_SELECTED_LIGANDS",
            "MMGBSA_INPUT_PH_LABELS",
            "MMGBSA_PROTOCOL",
            "MMGBSA_PRODUCTION_STRICT",
            "MMGBSA_TRAJ_MODE",
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES",
            "MMGBSA_MD_ENABLED",
            "MMGBSA_MD_SOLVENT_MODEL",
            "MMGBSA_MD_RUN",
            "MMGBSA_MD_ENGINE",
            "MMGBSA_MD_NREPLICATES",
            "MMGBSA_MD_SEED_MODE",
            "MMGBSA_MD_BASE_SEED",
            "MMGBSA_MD_IGB",
            "MMGBSA_MD_SALTCON",
            "MMGBSA_MD_RESTRAIN_PROTEIN_HEAVY",
            "MMGBSA_MD_RESTRAINT_WT",
            "MMGBSA_MD_RESTRAINT_MASK",
            "MMGBSA_MD_DT_PS",
            "MMGBSA_MD_MIN_STEPS",
            "MMGBSA_MD_RESCUE_ENABLED",
            "MMGBSA_MD_RESCUE_STEPS",
            "MMGBSA_MD_HEAT_PS",
            "MMGBSA_MD_EQUIL_PS",
            "MMGBSA_MD_PROD_PS",
            "MMGBSA_MD_TEMP0",
            "MMGBSA_MD_NTT",
            "MMGBSA_MD_GAMMA_LN",
            "MMGBSA_MD_NTC",
            "MMGBSA_MD_NTF",
            "MMGBSA_MD_FRAME_STRIDE_PS",
            "MMGBSA_MD_TRAJ_FORMAT",
            "MMGBSA_MD_TRAJ_NAME",
            "MMGBSA_MD_REP_AGG",
            "MMGBSA_MD_COPY_BEST_REPLICATE",
            "MMGBSA_ANALYSIS_START_PS",
            "MMGBSA_ANALYSIS_END_PS",
            "MMGBSA_ANALYSIS_STARTFRAME",
            "MMGBSA_ANALYSIS_ENDFRAME",
            "MMGBSA_ANALYSIS_INTERVAL",
            "MMGBSA_QC_ENABLED",
            "MMGBSA_QC_RUN_CPPTRAJ",
            "MMGBSA_QC_MIN_FRAMES",
            "MMGBSA_QC_MIN_BLOCKS",
            "MMGBSA_QC_BLOCK_SIZE",
            "MMGBSA_QC_MAX_SEM_KCAL",
            "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL",
            "MMGBSA_PREHEAT_QC_MAX_GMAX",
            "MMGBSA_PREHEAT_QC_MAX_RMS",
            "MMGBSA_PREHEAT_QC_MAX_ABS_ENERGY",
            "MMGBSA_PREHEAT_QC_MAX_ABS_VDW",
            "MMGBSA_PREHEAT_QC_MAX_NONFINITE_MARKERS",
            "MMGBSA_CHEMISTRY_REVIEW_REQUIRED",
            "MMGBSA_CHEMISTRY_REVIEW_STATUS",
            "MMGBSA_CHEMISTRY_REVIEW_NOTES",
            "MMGBSA_CHEMISTRY_SUBGATES_REQUIRED",
            "MMGBSA_PROTONATION_REVIEW_STATUS",
            "MMGBSA_TAUTOMER_REVIEW_STATUS",
            "MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS",
            "MMGBSA_NET_CHARGE_REVIEW_STATUS",
            "MMGBSA_PARAMETER_REVIEW_STATUS",
            "MMGBSA_LIGAND_CHEMISTRY_STRICT",
            "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE",
            "MMGBSA_INPUT_CHEMISTRY_SOURCE",
            "MMGBSA_INPUT_CHEMISTRY_NOTES",
            "MMGBSA_EXPLICIT_TRAJ_PATH",
            "MMGBSA_FRAME_AGG",
            "MMGBSA_FRAME_LIMIT",
        ],
    )
    if not _to_bool(cfg.get("MMGBSA_ENABLED", False)):
        logger.info("[mmgbsa.pipeline] action=skip reason=disabled")
        return
    _mmgbsa_validate_protocol(cfg, logger, context=f"{run_id}:{pdb_id}")

    variant_dir = (variant_token or "legacy").upper()
    post_docked_root = Path(
        str(
            cfg.get("POST_DOCKED_DIR")
            or output_root(Path(str(cfg.get("OVERALL_DIR", "."))), "post_docked")
        )
    ).expanduser()
    if not post_docked_root.is_absolute():
        post_docked_root = Path(str(cfg.get("OVERALL_DIR", "."))) / post_docked_root
    post_root = post_docked_root / run_id / pdb_id
    post_base = post_root if legacy_mode else post_root / variant_dir
    if not post_base.exists():
        logger.info(
            "[mmgbsa.pipeline] action=skip reason=missing_post_docked path=%s",
            post_base,
        )
        return

    try:
        max_ligands_val = int(cfg.get("MMGBSA_MAX_LIGANDS", 1))
    except INT_COERCE_ERRORS:
        max_ligands_val = 1
    if max_ligands_val <= 0:
        max_ligands_val = None

    try:
        top_pct_val = float(cfg.get("MMGBSA_RERANKED_TOP_PCT", 0))
    except FLOAT_COERCE_ERRORS:
        top_pct_val = 0.0
    if top_pct_val <= 0:
        top_pct_val = None

    try:
        poses_per_ligand = int(cfg.get("MMGBSA_POSES_PER_LIGAND", 1))
    except INT_COERCE_ERRORS:
        poses_per_ligand = 1
    if poses_per_ligand <= 0:
        poses_per_ligand = 1

    pose_sort_mode = (
        str(cfg.get("MMGBSA_POSE_SORT_MODE", "name_numeric") or "name_numeric")
        .strip()
        .lower()
    )
    pose_group_regexes = _mmgbsa_pose_group_regexes(cfg, logger)

    agg_enabled = _to_bool(cfg.get("MMGBSA_AGGREGATE_PER_LIGAND", False))
    agg_method = str(cfg.get("MMGBSA_AGG_METHOD", "min") or "min").strip().lower()
    agg_output_name = str(
        cfg.get("MMGBSA_AGG_OUTPUT_CSV", "mmgbsa_pose_aggregate.csv")
        or "mmgbsa_pose_aggregate.csv"
    )

    md_cfg = _mmgbsa_effective_md_config(cfg, logger)
    md_enabled = _to_bool(md_cfg.get("md_enabled", False))
    effective_mode = str(md_cfg.get("effective_mode", "ONEFRAME") or "ONEFRAME")
    five_reps = _to_bool(cfg.get("MD_FIVE_REPLICATE", False))
    try:
        md_reps = int(cfg.get("MMGBSA_MD_NREPLICATES", 1))
    except INT_COERCE_ERRORS:
        md_reps = 1
    if md_reps < 1:
        md_reps = 1
    md_rep_agg = str(cfg.get("MMGBSA_MD_REP_AGG", "mean") or "mean").strip().lower()
    md_copy_best = _to_bool(cfg.get("MMGBSA_MD_COPY_BEST_REPLICATE", False))
    frame_agg_method = (
        str(cfg.get("MMGBSA_FRAME_AGG", "mean") or "mean").strip().lower()
    )
    try:
        frame_limit = int(cfg.get("MMGBSA_FRAME_LIMIT", 0))
    except INT_COERCE_ERRORS:
        frame_limit = 0
    if frame_limit < 0:
        frame_limit = 0

    force = _to_bool(cfg.get("MMGBSA_FORCE", False))
    strict = _to_bool(cfg.get("MMGBSA_STRICT", False))
    stage_dir_name = str(cfg.get("MMGBSA_INPUT_STAGE_DIR", "stage1") or "stage1")
    strip_all_h_for_leap = _to_bool(cfg.get("MMGBSA_TLEAP_STRIP_ALL_H", True))
    map_hoh_to_wat = _to_bool(cfg.get("MMGBSA_TLEAP_MAP_HOH_TO_WAT", True))
    water_model = (
        str(cfg.get("MMGBSA_TLEAP_WATER_MODEL", "tip3p") or "tip3p").strip().lower()
    )
    if water_model not in {"tip3p"}:
        raise ValueError(
            f"MMGBSA_TLEAP_WATER_MODEL supports tip3p only (got {water_model})"
        )
    mmgbsa_logger = ppm._get_logger()

    test_override: Dict[str, object] = {}
    if test_mode != "off" and top_pct_val is None:
        test_stage = stage_dir_name
        candidate_dir = post_base / "pH7_0" / test_stage
        if candidate_dir.is_dir():
            actives = sorted(candidate_dir.glob("actives_final*.sdf"))
            if actives:
                test_override = {
                    "ph_dir": candidate_dir.parent,
                    "sdfs": [actives[0]],
                    "stage_dir": test_stage,
                }
            else:
                sdfs = sorted(candidate_dir.glob("*.sdf"))
                if sdfs:
                    test_override = {
                        "ph_dir": candidate_dir.parent,
                        "sdfs": [sdfs[0]],
                        "stage_dir": test_stage,
                    }

    if test_override:
        ph_dirs = [test_override["ph_dir"]]
    else:
        ph_dirs = sorted([p for p in post_base.iterdir() if p.is_dir()])
        selected_ph_labels = {
            label for label in _mmgbsa_split_cfg_list(cfg.get("MMGBSA_INPUT_PH_LABELS"))
        }
        if selected_ph_labels:
            ph_dirs = [p for p in ph_dirs if p.name in selected_ph_labels]

    if not ph_dirs:
        logger.info(
            "[mmgbsa.pipeline] action=skip reason=no_ph_dirs path=%s", post_base
        )
        return

    amber_prefix = cfg.get("MMGBSA_AMBERTOOLS_PREFIX") or cfg.get("AMBERTOOLS_PREFIX")
    amber_prefix = str(amber_prefix).strip() if amber_prefix else None
    explicit_ligands = _mmgbsa_split_cfg_list(cfg.get("MMGBSA_SELECTED_LIGANDS"))

    for ph_dir in ph_dirs:
        ph_label = ph_dir.name
        if test_override:
            sdfs = list(test_override["sdfs"])
            stage_dir_label = str(test_override["stage_dir"])
            pose_groups = _mmgbsa_group_pose_sdfs(
                sdfs, stage_dir_label, pose_group_regexes, logger
            )
            selected_ligands = sorted(pose_groups.keys())
            logger.info(
                "[mmgbsa.pipeline] selection=poses stage_dir=%s ligand_ids=%d poses_per_ligand=%d poses_selected=%d reranked=%s actives_only=%s",
                stage_dir_label,
                len(selected_ligands),
                poses_per_ligand,
                len(sdfs),
                False,
                False,
            )
        else:
            stage_dir_label = stage_dir_name
            stage_dir = ph_dir / stage_dir_label
            reranked_csv = ph_dir / "consensus_reranked_scorch.csv"
            sdfs, _, _, selected_ligands = _mmgbsa_select_sdfs(
                stage_dir,
                max_ligands_val,
                poses_per_ligand,
                pose_sort_mode,
                pose_group_regexes,
                reranked_csv=reranked_csv,
                top_pct=top_pct_val,
                selected_ligand_ids=explicit_ligands,
                logger=logger,
            )

        if not sdfs:
            logger.info(
                "[mmgbsa.pipeline] action=skip reason=no_sdfs pdb=%s variant=%s ph=%s stage=%s",
                pdb_id,
                variant_dir,
                ph_label,
                stage_dir_label,
            )
            continue

        expected_frames = None
        if md_enabled:
            try:
                prod_ps = float(cfg.get("MMGBSA_MD_PROD_PS", 100.0))
            except Exception:
                prod_ps = 100.0
            try:
                stride_ps = float(cfg.get("MMGBSA_MD_FRAME_STRIDE_PS", 2.0))
            except Exception:
                stride_ps = 2.0
            if stride_ps > 0:
                expected_frames = max(1, int(prod_ps / stride_ps))

        mmpbsa_interval = 1
        for key in ("MMGBSA_GENERAL_INTERVAL", "MMGBSA_MMPBSA_INTERVAL"):
            val = cfg.get(key, None)
            if val is None or str(val).strip() == "":
                continue
            try:
                mmpbsa_interval = int(val)
            except Exception:
                mmpbsa_interval = 1
            break

        try:
            mmpbsa_igb = int(cfg.get("MMGBSA_GB_IGB", 5))
        except Exception:
            mmpbsa_igb = 5
        try:
            mmpbsa_saltcon = float(cfg.get("MMGBSA_GB_SALTCON", 0.150))
        except Exception:
            mmpbsa_saltcon = 0.150
        logger.info(
            "[mmgbsa.pipeline] effective_cfg mode=%s md_enabled=%s expected_frames_per_rep=%s replicates=%d mmpbsa_interval=%d igb=%d saltcon=%.3f",
            effective_mode,
            md_enabled,
            expected_frames if md_enabled else "n/a",
            md_reps,
            mmpbsa_interval,
            mmpbsa_igb,
            mmpbsa_saltcon,
        )

        center, radius, cfg_path = _mmgbsa_resolve_center_radius(
            cfg,
            run_id,
            pdb_id,
            stage_dir_label,
            variant_token,
            ph_label,
            legacy_mode,
        )
        if center is None:
            center = _mmgbsa_ligand_centroid(sdfs)
            if center is None:
                logger.warning(
                    "[mmgbsa.pipeline] action=skip reason=missing_center pdb=%s variant=%s ph=%s cfg=%s",
                    pdb_id,
                    variant_dir,
                    ph_label,
                    cfg_path,
                )
                continue
            logger.info(
                "[mmgbsa.pipeline] center_source=ligand_centroid pdb=%s variant=%s ph=%s center=%.3f,%.3f,%.3f",
                pdb_id,
                variant_dir,
                ph_label,
                center[0],
                center[1],
                center[2],
            )

        receptor_candidates = [
            ph_dir / stage_dir_label / "receptor.pdb",
            ph_dir / "receptor.pdb",
        ]
        receptor_pdb = next(
            (cand for cand in receptor_candidates if cand.exists()), None
        )
        if receptor_pdb is None:
            receptor_pdb = _mmgbsa_resolve_receptor_pdb(
                cfg,
                pdb_id,
                pdb_file,
                variant_token,
                ph_label,
                legacy_mode,
            )

        if receptor_pdb is None:
            logger.warning(
                "[mmgbsa.pipeline] action=skip reason=missing_receptor pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_dir,
                ph_label,
            )
            continue

        receptor_input = receptor_pdb
        if "post_docked" not in receptor_pdb.parts:
            receptor_input_dir = ph_dir / "mmgbsa_receptor_inputs"
            receptor_input_dir.mkdir(parents=True, exist_ok=True)
            receptor_input = receptor_input_dir / receptor_pdb.name
            if not receptor_input.exists() or force:
                try:
                    shutil.copy2(receptor_pdb, receptor_input)
                except Exception as exc:
                    logger.error(
                        "[mmgbsa.pipeline] action=skip reason=receptor_copy_failed pdb=%s variant=%s ph=%s err=%s",
                        pdb_id,
                        variant_dir,
                        ph_label,
                        exc,
                    )
                    if strict:
                        raise
                    continue
            logger.info(
                "[mmgbsa.pipeline] receptor_materialized source=%s dest=%s",
                receptor_pdb,
                receptor_input,
            )

        logger.info(
            "[mmgbsa.pipeline] start pdb=%s variant=%s ph=%s ligands=%d",
            pdb_id,
            variant_dir,
            ph_label,
            len(sdfs),
        )

        try:
            receptor_result = prep_mmgbsa_receptor(
                pdb_path=str(receptor_input),
                runid=run_id,
                center=center,
                radius=radius,
                force=force,
            )
        except Exception as exc:
            logger.error(
                "[mmgbsa.pipeline] action=skip reason=receptor_prep_failed pdb=%s variant=%s ph=%s err=%s",
                pdb_id,
                variant_dir,
                ph_label,
                exc,
            )
            if strict:
                raise
            continue

        receptor_for_leap_path, strip_info = ppm._strip_receptor_h_for_leap(
            Path(receptor_result["output_path"]),
            strip_all_h_for_leap,
        )
        hoh_residue_count = strip_info.get("hoh_residue_count", 0)
        ppm._log_leap_prep(
            mmgbsa_logger,
            "INFO",
            {
                "strip_all_h": strip_all_h_for_leap,
                "input": receptor_result["output_path"],
                "output": str(receptor_for_leap_path),
                "removed_H": strip_info.get("removed_h", 0),
                "hoh_residues": hoh_residue_count,
            },
        )
        ppm._log_leap_prep(
            mmgbsa_logger,
            "INFO",
            {
                "map_hoh_to_wat": map_hoh_to_wat,
                "water_model": water_model,
                "hoh_residues": hoh_residue_count,
            },
            "water_mapping",
        )

        mmgbsa_dir = ph_dir / "mmgbsa"
        summary_path = mmgbsa_dir / "mmgbsa_results_summary.csv"
        work_root = mmgbsa_dir / "work"
        try:
            write_methods_json(
                mmgbsa_dir / "mmgbsa_methods.json",
                cfg,
                run_id=run_id,
                pdb_id=pdb_id,
                variant=variant_dir,
                ph_label=ph_label,
            )
        except Exception as exc:
            logger.warning(
                "[mmgbsa.methods] action=write_failed pdb=%s variant=%s ph=%s err=%s",
                pdb_id,
                variant_dir,
                ph_label,
                exc,
            )

        prep_cfg = dict(cfg)
        source_paths = make_paths(cfg, base_id=pdb_id, pdb_file=pdb_file)
        prep_cfg["MMGBSA_RECEPTOR_PDB_CONTEXT"] = str(
            source_paths.input_pdb_path
            if source_paths.input_pdb_path.exists()
            else receptor_input
        )
        prep_results = prep_mmgbsa_from_sdfs(
            [str(p) for p in sdfs],
            cfg=prep_cfg,
            max_ligands=0,
            force=force,
            amber_prefix=amber_prefix,
        )

        ligand_entries: List[dict] = []
        prep_result_by_ligand: Dict[str, Dict[str, object]] = {}
        for res in prep_results:
            if res.get("error"):
                _mmgbsa_append_summary(
                    summary_path,
                    {
                        "stage_dir": res.get("stage_dir", stage_dir_label),
                        "ligand_stem": Path(res.get("sdf_path", "ligand")).stem,
                        "delta_total": "",
                        "results_csv": "",
                        "ok": False,
                        "notes": res.get("error", "prep_failed"),
                    },
                )
                continue
            ligand_base = Path(res.get("sdf_path", "ligand")).stem
            prep_result_by_ligand[ligand_base] = dict(res)
            ligand_entries.append(
                {
                    "stage_dir": res.get("stage_dir", stage_dir_label),
                    "ligand_base": ligand_base,
                    "mol2_path": res.get("mol2_path"),
                    "frcmod_path": res.get("frcmod_path"),
                }
            )

        if not ligand_entries:
            logger.info(
                "[mmgbsa.pipeline] action=skip reason=ligand_prep_failed pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_dir,
                ph_label,
            )
            continue

        topo_results = write_leap_for_ligands(
            receptor_pdb_path=str(receptor_for_leap_path),
            ligand_mol2_frcmod_pairs=ligand_entries,
            out_dir_base=str(work_root),
            force=force or _to_bool(cfg.get("MMGBSA_TLEAP_FORCE", False)),
            water_model=water_model,
            map_hoh_to_wat=map_hoh_to_wat,
            hoh_residue_count=hoh_residue_count,
        )
        ligand_entry_by_base = {
            str(entry["ligand_base"]): entry for entry in ligand_entries
        }
        for topo in topo_results:
            ligand_base = str(topo.get("ligand_base", ""))
            entry = ligand_entry_by_base.get(ligand_base)
            if entry:
                topo["ligand_mol2"] = entry.get("mol2_path")
                topo["ligand_frcmod"] = entry.get("frcmod_path")

        run_tleap_flag = _to_bool(cfg.get("MMGBSA_TLEAP_ENABLED", True)) and _to_bool(
            cfg.get("MMGBSA_TLEAP_RUN", True)
        )

        for topo in topo_results:
            ligand_stem = topo.get("ligand_base")
            stage_dir = topo.get("stage_dir")
            out_dir = Path(topo.get("output_dir", work_root))
            notes = ""
            ok = False
            results_csv = ""
            delta_total = ""
            summary_frames: Dict[str, object] = {
                "frame_agg_method": frame_agg_method,
                "frames_used": "",
                "frame_mean": "",
                "frame_sd": "",
                "frames_total": "",
                "frame_median": "",
                "frame_sem": "",
                "frame_ci95": "",
                "qc_pass": "",
                "qc_json": "",
            }
            replicates_total = md_reps if md_enabled else 1
            replicates_ok = 0
            rep_agg_method_used = md_rep_agg if md_enabled else "single"
            rep_stats: Dict[str, object] = {}
            rep_outlier_text = ""

            ppm._log_leap(
                mmgbsa_logger,
                "INFO",
                {
                    "stage_dir": stage_dir,
                    "ligand": ligand_stem,
                    "receptor_pdb_used": str(receptor_for_leap_path),
                    "build_leap_path": topo.get("leap_file"),
                    "tleap_log_path": out_dir / "tleap.log",
                    "tleap_run": run_tleap_flag,
                    "skip": topo.get("skip_tleap"),
                },
                "ligand_topology",
            )

            try:
                review_result = _mmgbsa_write_review_record(
                    cfg=cfg,
                    mmgbsa_dir=mmgbsa_dir,
                    pdb_id=pdb_id,
                    ligand_stem=str(ligand_stem),
                    variant_dir=variant_dir,
                    ph_label=ph_label,
                    receptor_result=receptor_result,
                    prep_result=prep_result_by_ligand.get(str(ligand_stem)),
                    logger=logger,
                )
                review_validation = review_result.get("validation", {})
                review_problems = review_validation.get("problems", [])
                if _mmgbsa_protocol(cfg) == "production" and review_problems:
                    notes = "review_record_not_publication_ready"
                    _mmgbsa_append_summary(
                        summary_path,
                        {
                            "stage_dir": stage_dir,
                            "ligand_stem": ligand_stem,
                            "delta_total": "",
                            "results_csv": "",
                            "ok": False,
                            "notes": notes,
                        },
                    )
                    if strict:
                        raise ValueError(
                            "MMGBSA review record invalid: "
                            + "; ".join(str(item) for item in review_problems)
                        )
                    continue

                if run_tleap_flag and not topo.get("skip_tleap"):
                    run_tleap(topo["leap_file"], str(out_dir))

                required = {
                    "complex_prmtop": topo.get("complex_prmtop"),
                    "complex_inpcrd": topo.get("complex_inpcrd"),
                    "receptor_prmtop": topo.get("receptor_prmtop"),
                    "ligand_prmtop": topo.get("ligand_prmtop"),
                }
                missing = []
                for key, value in required.items():
                    if not value:
                        missing.append(key)
                        continue
                    path = Path(value)
                    if not path.exists() or path.stat().st_size == 0:
                        missing.append(key)

                if missing:
                    notes = f"missing_topology:{','.join(missing)}"
                    _mmgbsa_append_summary(
                        summary_path,
                        {
                            "stage_dir": stage_dir,
                            "ligand_stem": ligand_stem,
                            "delta_total": "",
                            "results_csv": "",
                            "ok": False,
                            "notes": notes,
                        },
                    )
                    continue

                if five_reps:
                    rep_result = _mmgbsa_five_replicate_runner(
                        topo=topo,
                        out_dir=out_dir,
                        cfg=cfg,
                        force=force,
                        md_enabled=md_enabled,
                        stage_dir=stage_dir,
                        ligand_stem=ligand_stem,
                        pdb_id=pdb_id,
                        variant_dir=variant_dir,
                        ph_label=ph_label,
                        run_id=run_id,
                        logger=logger,
                    )
                    ok = bool(rep_result.get("ok"))
                    delta_val = rep_result.get("delta_total")
                    if delta_val is not None:
                        delta_total = f"{float(delta_val):.6g}"
                    results_csv = rep_result.get("results_csv", "")
                    notes = rep_result.get("notes", "")
                    summary_frames.update(
                        {
                            "frame_mean": rep_result.get("frame_mean", ""),
                            "frame_sd": rep_result.get("frame_sd", ""),
                        }
                    )
                    replicates_ok = int(rep_result.get("replicates_ok", 0))
                    replicates_total = int(rep_result.get("replicates_total", 5))
                    _mmgbsa_append_summary(
                        summary_path,
                        {
                            "stage_dir": stage_dir,
                            "ligand_stem": ligand_stem,
                            "delta_total": delta_total,
                            "results_csv": results_csv,
                            "ok": ok,
                            "notes": notes or "",
                            "frame_agg_method": summary_frames.get(
                                "frame_agg_method", ""
                            )
                            or "",
                            "frames_used": summary_frames.get("frames_used", "") or "",
                            "frames_total": summary_frames.get("frames_total", "")
                            or "",
                            "frame_mean": summary_frames.get("frame_mean", "") or "",
                            "frame_sd": summary_frames.get("frame_sd", "") or "",
                            "frame_median": summary_frames.get("frame_median", "")
                            or "",
                            "rep_agg_method": rep_agg_method_used,
                            "replicates_ok": replicates_ok,
                            "replicates_total": replicates_total,
                        },
                    )
                    continue

                if not md_enabled:
                    traj_result = make_mmgbsa_trajectory(
                        complex_prmtop=topo["complex_prmtop"],
                        complex_inpcrd=topo["complex_inpcrd"],
                        out_dir=str(out_dir),
                        cfg=cfg,
                        force=force,
                        run_cpptraj=_to_bool(cfg.get("MMGBSA_CPPTRAJ_RUN", True)),
                    )
                    traj_path = Path(traj_result.get("trajout_path") or "")
                    if not traj_path.exists() or traj_path.stat().st_size == 0:
                        default_traj = str(
                            cfg.get("MMGBSA_DEFAULT_TRAJ_NAME", "mdcrd") or "mdcrd"
                        )
                        candidate = out_dir / default_traj
                        if candidate.exists() and candidate.stat().st_size > 0:
                            traj_path = candidate
                        else:
                            notes = "missing_trajectory"
                            _mmgbsa_append_summary(
                                summary_path,
                                {
                                    "stage_dir": stage_dir,
                                    "ligand_stem": ligand_stem,
                                    "delta_total": "",
                                    "results_csv": "",
                                    "ok": False,
                                    "notes": notes,
                                },
                            )
                            continue

                    mmpbsa_result = run_mmgbsa(
                        complex_prmtop=topo["complex_prmtop"],
                        receptor_prmtop=topo["receptor_prmtop"],
                        ligand_prmtop=topo["ligand_prmtop"],
                        trajectory_path=str(traj_path),
                        work_dir=str(out_dir),
                        cfg=cfg,
                        force=force,
                        run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
                    )

                    results_csv = mmpbsa_result.get("out_csv", "")
                    results_dat = mmpbsa_result.get("out_dat", "")
                    if not mmpbsa_result.get("enabled", True):
                        notes = "mmpbsa_disabled"
                    else:
                        if (
                            results_csv
                            and results_dat
                            and Path(results_csv).exists()
                            and Path(results_dat).exists()
                        ):
                            frame_meta = _mmgbsa_frame_aggregate(
                                Path(results_csv), cfg, logger
                            )
                            delta_val = frame_meta.get("delta_total")
                            if delta_val is not None:
                                delta_total = f"{float(delta_val):.6g}"
                                ok = _mmgbsa_frame_qc_ok(cfg, frame_meta)
                                if not ok:
                                    notes = "frame_qc_failed"
                            else:
                                notes = (
                                    frame_meta.get("notes", "missing_delta")
                                    or "missing_delta"
                                )
                            summary_frames.update(_mmgbsa_summary_frame_fields(frame_meta))
                            replicates_ok = 1 if ok else 0
                        else:
                            notes = "missing_outputs"
                else:
                    rep_seeds = _mmgbsa_md_seed_list(cfg, md_reps, logger)
                    solvent_model = (
                        str(cfg.get("MMGBSA_MD_SOLVENT_MODEL", "implicit") or "implicit")
                        .strip()
                        .lower()
                    )
                    rep_results: List[Dict[str, object]] = []
                    for rep_idx, seed in enumerate(rep_seeds, start=1):
                        rep_dir = out_dir / f"rep{rep_idx}"
                        if solvent_model == "explicit":
                            md_result = build_and_run_explicit_workflow(
                                receptor_pdb=str(receptor_for_leap_path),
                                ligand_mol2=str(topo.get("ligand_mol2") or ""),
                                ligand_frcmod=str(topo.get("ligand_frcmod") or ""),
                                dry_complex_prmtop=topo["complex_prmtop"],
                                dry_receptor_prmtop=topo["receptor_prmtop"],
                                dry_ligand_prmtop=topo["ligand_prmtop"],
                                out_dir=rep_dir / "explicit",
                                cfg=cfg,
                                seed=int(seed),
                                force=force,
                                run=True,
                            )
                        else:
                            md_result = run_implicit_md(
                                complex_prmtop=topo["complex_prmtop"],
                                complex_inpcrd=topo["complex_inpcrd"],
                                out_dir=str(out_dir),
                                cfg=cfg,
                                replicate_index=rep_idx,
                                seed=int(seed),
                                force=force,
                                run=True,
                            )
                        traj_path = Path(md_result.get("traj_path") or "")
                        rep_notes = ""
                        rep_ok = False
                        rep_delta = None
                        rep_frames: Dict[str, object] = {
                            "frame_agg_method": frame_agg_method,
                            "frames_used": "",
                            "frames_total": "",
                            "frame_mean": "",
                            "frame_sd": "",
                            "frame_median": "",
                            "frame_sem": "",
                            "frame_ci95": "",
                            "qc_pass": "",
                            "qc_json": "",
                        }
                        rep_csv = ""
                        rep_dat = ""
                        rep_log = ""

                        if not md_result.get("ok"):
                            failed_step = str(md_result.get("step_failed") or "").strip()
                            rep_notes = f"md_failed:{failed_step}" if failed_step else "md_failed"
                        elif not traj_path.exists() or traj_path.stat().st_size == 0:
                            rep_notes = "missing_trajectory"
                        else:
                            try:
                                mmpbsa_result = run_mmgbsa(
                                    complex_prmtop=topo["complex_prmtop"],
                                    receptor_prmtop=topo["receptor_prmtop"],
                                    ligand_prmtop=topo["ligand_prmtop"],
                                    trajectory_path=str(traj_path),
                                    work_dir=str(rep_dir),
                                    cfg=cfg,
                                    force=force,
                                    run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
                                )
                                rep_csv = mmpbsa_result.get("out_csv", "")
                                rep_dat = mmpbsa_result.get("out_dat", "")
                                rep_log = mmpbsa_result.get("log_path", "")
                                if not mmpbsa_result.get("enabled", True):
                                    rep_notes = "mmpbsa_disabled"
                                elif (
                                    rep_csv
                                    and rep_dat
                                    and Path(rep_csv).exists()
                                    and Path(rep_dat).exists()
                                ):
                                    frame_meta = _mmgbsa_frame_aggregate(
                                        Path(rep_csv), cfg, logger
                                    )
                                    delta_val = frame_meta.get("delta_total")
                                    rep_frames.update(_mmgbsa_summary_frame_fields(frame_meta))
                                    if delta_val is not None:
                                        rep_delta = float(delta_val)
                                        rep_ok = _mmgbsa_frame_qc_ok(cfg, frame_meta)
                                        if not rep_ok:
                                            rep_notes = "frame_qc_failed"
                                    else:
                                        rep_notes = (
                                            frame_meta.get("notes", "missing_delta")
                                            or "missing_delta"
                                        )
                                else:
                                    rep_notes = "missing_outputs"
                            except Exception as exc:
                                rep_notes = f"mmpbsa_error:{type(exc).__name__}"
                                logger.warning(
                                    "[mmgbsa.pipeline] md_replicate_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s rep=%d err=%s",
                                    pdb_id,
                                    variant_dir,
                                    ph_label,
                                    stage_dir,
                                    ligand_stem,
                                    rep_idx,
                                    exc,
                                )
                                if strict:
                                    raise

                        rep_results.append(
                            {
                                "replicate": rep_idx,
                                "seed": seed,
                                "ok": rep_ok,
                                "delta_total": rep_delta,
                                "results_csv": rep_csv,
                                "results_dat": rep_dat,
                                "log_path": rep_log,
                                "traj_path": str(traj_path),
                                "work_dir": str(rep_dir),
                                "notes": rep_notes,
                                "frames_used": rep_frames.get("frames_used", ""),
                                "frames_total": rep_frames.get("frames_total", ""),
                                "frame_mean": rep_frames.get("frame_mean", ""),
                                "frame_sd": rep_frames.get("frame_sd", ""),
                                "frame_sem": rep_frames.get("frame_sem", ""),
                                "frame_ci95": rep_frames.get("frame_ci95", ""),
                                "qc_pass": rep_frames.get("qc_pass", ""),
                                "qc_json": rep_frames.get("qc_json", ""),
                                "frame_agg_method": rep_frames.get(
                                    "frame_agg_method", frame_agg_method
                                ),
                                "frame_median": rep_frames.get("frame_median", ""),
                            }
                        )

                    ok_reps = [r for r in rep_results if r.get("ok")]
                    n_ok = len(ok_reps)
                    replicates_ok = n_ok
                    ok_deltas = [
                        float(r["delta_total"])
                        for r in ok_reps
                        if r.get("delta_total") is not None
                    ]
                    agg_delta_val = _mmgbsa_md_aggregate(
                        ok_deltas,
                        md_rep_agg,
                    )
                    rep_stats = summarize_values(ok_deltas)
                    outlier_indexes = detect_outliers(ok_deltas)
                    if outlier_indexes:
                        rep_outlier_text = ";".join(
                            str(ok_reps[idx].get("replicate", idx + 1))
                            for idx in outlier_indexes
                            if idx < len(ok_reps)
                        )
                        for idx in outlier_indexes:
                            if idx < len(ok_reps):
                                ok_reps[idx]["outlier"] = True
                    if agg_delta_val is not None:
                        delta_total = f"{agg_delta_val:.6g}"
                        ok = True
                    else:
                        ok = False

                    best_rep = None
                    if ok_reps:
                        if md_rep_agg == "min":
                            best_rep = min(
                                ok_reps,
                                key=lambda r: (
                                    r.get("delta_total", 0),
                                    r.get("replicate", 0),
                                ),
                            )
                        elif md_rep_agg == "median":
                            best_rep = min(
                                ok_reps,
                                key=lambda r: (
                                    abs(
                                        float(r.get("delta_total", 0))
                                        - float(agg_delta_val or 0)
                                    ),
                                    r.get("replicate", 0),
                                ),
                            )
                        else:
                            best_rep = min(
                                ok_reps,
                                key=lambda r: (
                                    abs(
                                        float(r.get("delta_total", 0))
                                        - float(agg_delta_val or 0)
                                    ),
                                    r.get("replicate", 0),
                                ),
                            )

                    if best_rep:
                        results_csv = str(best_rep.get("results_csv", ""))
                        summary_frames.update(
                            {
                                "frame_agg_method": best_rep.get(
                                    "frame_agg_method", frame_agg_method
                                ),
                                "frames_used": best_rep.get("frames_used", ""),
                                "frame_mean": best_rep.get("frame_mean", ""),
                                "frame_sd": best_rep.get("frame_sd", ""),
                                "frames_total": best_rep.get("frames_total", ""),
                                "frame_median": best_rep.get("frame_median", ""),
                                "frame_sem": best_rep.get("frame_sem", ""),
                                "frame_ci95": best_rep.get("frame_ci95", ""),
                                "qc_pass": best_rep.get("qc_pass", ""),
                                "qc_json": best_rep.get("qc_json", ""),
                            }
                        )
                    failed_notes = sorted(
                        {
                            str(rep.get("notes", "") or "")
                            for rep in rep_results
                            if not rep.get("ok") and str(rep.get("notes", "") or "")
                        }
                    )
                    notes = f"md_reps_ok={n_ok}/{len(rep_results)} rep_agg={md_rep_agg} frame_agg={frame_agg_method}"
                    if failed_notes:
                        notes += f" failures={';'.join(failed_notes)}"

                    replicate_summary = {
                        "ok": ok,
                        "n_reps_total": len(rep_results),
                        "n_reps_ok": n_ok,
                        "agg_method": md_rep_agg,
                        "agg_delta": agg_delta_val,
                        "replicate_stats": rep_stats,
                        "outlier_replicates": rep_outlier_text,
                        "best_replicate": best_rep,
                        "replicates": rep_results,
                    }
                    summary_json = out_dir / "mmgbsa_replicate_summary.json"
                    try:
                        tmp_json = summary_json.with_suffix(".json.part")
                        with tmp_json.open("w", encoding="utf-8") as handle:
                            json.dump(
                                replicate_summary, handle, indent=2, sort_keys=True
                            )
                        os.replace(tmp_json, summary_json)
                    except Exception as exc:
                        logger.warning(
                            "[mmgbsa.pipeline] md_replicate_summary_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s err=%s",
                            pdb_id,
                            variant_dir,
                            ph_label,
                            stage_dir,
                            ligand_stem,
                            exc,
                        )

                    if md_copy_best and best_rep:
                        for filename in (
                            "FINAL_RESULTS_MMPBSA.dat",
                            "FINAL_RESULTS_MMPBSA.csv",
                            "mmpbsa.log",
                            "mmpbsa.in",
                        ):
                            src = Path(best_rep.get("work_dir", "")) / filename
                            dest = out_dir / filename
                            if not src.exists():
                                continue
                            if dest.exists() and not force:
                                continue
                            try:
                                shutil.copy2(src, dest)
                            except Exception as exc:
                                logger.warning(
                                    "[mmgbsa.pipeline] md_copy_best_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s file=%s err=%s",
                                    pdb_id,
                                    variant_dir,
                                    ph_label,
                                    stage_dir,
                                    ligand_stem,
                                    filename,
                                    exc,
                                )

            except Exception as exc:
                notes = f"error:{type(exc).__name__}"
                ok = False
                logger.error(
                    "[mmgbsa.pipeline] ligand_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s err=%s",
                    pdb_id,
                    variant_dir,
                    ph_label,
                    stage_dir,
                    ligand_stem,
                    exc,
                )
                if strict:
                    raise

            _mmgbsa_append_summary(
                summary_path,
                {
                    "stage_dir": stage_dir,
                    "ligand_stem": ligand_stem,
                    "delta_total": delta_total,
                    "results_csv": results_csv,
                    "ok": ok,
                    "notes": notes,
                    "frame_agg_method": summary_frames.get("frame_agg_method", "")
                    or "",
                    "frames_used": summary_frames.get("frames_used", "") or "",
                    "frames_total": summary_frames.get("frames_total", "") or "",
                    "frame_mean": summary_frames.get("frame_mean", "") or "",
                    "frame_sd": summary_frames.get("frame_sd", "") or "",
                    "frame_median": summary_frames.get("frame_median", "") or "",
                    "frame_sem": summary_frames.get("frame_sem", "") or "",
                    "frame_ci95": summary_frames.get("frame_ci95", "") or "",
                    "qc_pass": summary_frames.get("qc_pass", "") or "",
                    "qc_json": summary_frames.get("qc_json", "") or "",
                    "rep_mean": rep_stats.get("mean", "") if rep_stats else "",
                    "rep_sd": rep_stats.get("sd", "") if rep_stats else "",
                    "rep_sem": rep_stats.get("sem", "") if rep_stats else "",
                    "rep_ci95": rep_stats.get("ci95", "") if rep_stats else "",
                    "rep_outliers": rep_outlier_text,
                    "rep_agg_method": rep_agg_method_used,
                    "replicates_ok": replicates_ok,
                    "replicates_total": replicates_total,
                },
            )

        if agg_enabled:
            agg_path = mmgbsa_dir / agg_output_name
            _mmgbsa_write_pose_aggregate(
                summary_path=summary_path,
                out_path=agg_path,
                stage_dir_label=stage_dir_label,
                pose_group_regexes=pose_group_regexes,
                agg_method=agg_method,
                logger=logger,
            )

        logger.info(
            "[mmgbsa.pipeline] done pdb=%s variant=%s ph=%s summary=%s",
            pdb_id,
            variant_dir,
            ph_label,
            summary_path,
        )
