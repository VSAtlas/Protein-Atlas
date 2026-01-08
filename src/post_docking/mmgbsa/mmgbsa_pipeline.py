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
import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
from pathlib import Path

from input_and_export_functions import _to_bool
from post_docking.mmgbsa.prep_for_mmgbsa import prep_mmgbsa_from_sdfs
import post_docking.mmgbsa.protein_prep_mmgbsa as ppm
from post_docking.mmgbsa.protein_prep_mmgbsa import prep_mmgbsa_receptor, run_tleap, write_leap_for_ligands
from post_docking.mmgbsa.mmgbsa_trajectory import make_mmgbsa_trajectory, run_implicit_md
from post_docking.mmgbsa.run_mmgbsa import run_mmgbsa, parse_mmpbsa_delta_total, write_aggregated_mmpbsa_results
from path_router.path_router import make_paths, config_dir as router_config_dir, ph_ensemble_dir as router_ph_ensemble_dir


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
        return None
    txt_files = sorted(cfg_dir.glob("*.txt"))
    return txt_files[0] if txt_files else None


def _mmgbsa_parse_vina_config(cfg_path: Path) -> tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]]]:
    center_vals: dict[str, float] = {}
    size_vals: dict[str, float] = {}
    try:
        lines = cfg_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None, None

    for line in lines:
        if "=" not in line:
            continue
        key, raw_val = line.split("=", 1)
        key = key.strip().lower()
        val_text = raw_val.strip()
        try:
            val = float(val_text)
        except Exception:
            continue
        if key in {"center_x", "center_y", "center_z"}:
            center_vals[key] = val
        elif key in {"size_x", "size_y", "size_z"}:
            size_vals[key] = val

    if len(center_vals) == 3:
        center = (center_vals["center_x"], center_vals["center_y"], center_vals["center_z"])
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
    cfg_path = _mmgbsa_find_vina_config(run_id, pdb_id, stage_dir, variant_token, ph_label, legacy_mode)
    if not cfg_path:
        return None, fallback_radius, None

    center, size = _mmgbsa_parse_vina_config(cfg_path)
    if center is None:
        return None, fallback_radius, cfg_path

    radius = fallback_radius
    if size is not None:
        radius = max(size) / 2.0
    return center, radius, cfg_path


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
        ensemble_dir = router_ph_ensemble_dir(pdb_id, variant=variant_token, legacy=legacy_mode)
        tag = f"{pdb_id}_{ph_label}"
        candidates.append(ensemble_dir / f"{tag}.withH.pdb")
        candidates.append(ensemble_dir / f"{tag}.pdb")

    candidates.append(paths.receptor_cleaned_pdb(variant_token))
    candidates.append(paths.receptor_dir(variant_token) / f"{pdb_id}.pdb")

    for cand in candidates:
        if cand.exists():
            return cand
    return None


def _mmgbsa_normalize_ligand_base(name: str) -> str:
    s = str(name or "").strip()
    if not s:
        return ""
    s = re.sub(r"\.(pdbqt|mol2|sdf)$", "", s, flags=re.IGNORECASE)
    s = s.replace(".sanitized", "")
    return s


def _mmgbsa_sdf_base_candidates(stem: str, stage_dir_label: str) -> List[str]:
    candidates: set[str] = set()
    base = stem.replace(".sanitized", "")
    candidates.add(base)

    if stage_dir_label:
        for sep in ("_", "__"):
            suffix = f"{sep}{stage_dir_label}"
            if base.endswith(suffix):
                candidates.add(base[: -len(suffix)])

    for pattern in (
        r"(_+gnina_stage\d+)$",
        r"(_+ledock_stage\d+)$",
        r"(_+dock6_stage\d+)$",
        r"(_+stage\d+)$",
    ):
        stripped = re.sub(pattern, "", base, flags=re.IGNORECASE)
        if stripped != base:
            candidates.add(stripped)

    cleaned = set()
    for cand in candidates:
        cleaned.add(re.sub(r"_+$", "", cand))
    return sorted(cleaned)


def _mmgbsa_default_pose_group_regexes() -> List[str]:
    return [
        r"(_pose\d+)$",
        r"(_rank\d+)$",
        r"(_conf\d+)$",
        r"(_model\d+)$",
        r"(_p\d+)$",
        r"(_stage\d+)(_pose\d+)$",
        r"(_gnina_stage\d+)(_pose\d+)$",
        r"(_ledock_stage\d+)(_pose\d+)$",
        r"(_dock6_stage\d+)(_pose\d+)$",
    ]


def _mmgbsa_pose_group_regexes(cfg: Mapping[str, Any], logger: logging.Logger) -> List[re.Pattern]:
    raw = str(cfg.get("MMGBSA_POSE_GROUP_REGEXES", "") or "").strip()
    patterns: List[str] = []
    if raw:
        for part in raw.split(","):
            part = part.strip()
            if part:
                patterns.append(part)
    if not patterns:
        patterns = _mmgbsa_default_pose_group_regexes()

    compiled: List[re.Pattern] = []
    for pat in patterns:
        try:
            compiled.append(re.compile(pat, flags=re.IGNORECASE))
        except re.error as exc:
            logger.warning(
                "[mmgbsa.pipeline] selection=pose_regex_invalid regex=%s err=%s",
                pat,
                exc,
            )
    if not compiled:
        compiled = [re.compile(pat, flags=re.IGNORECASE) for pat in _mmgbsa_default_pose_group_regexes()]
    return compiled


def _mmgbsa_pose_candidates(stem: str, stage_dir_label: str, regexes: List[re.Pattern]) -> List[str]:
    base = stem.replace(".sanitized", "")
    seeds: set[str] = {base}
    for regex in regexes:
        stripped = regex.sub("", base)
        if stripped != base:
            seeds.add(stripped)

    candidates: set[str] = set()
    for seed in seeds:
        for cand in _mmgbsa_sdf_base_candidates(seed, stage_dir_label):
            candidates.add(cand)
        for regex in regexes:
            stripped = regex.sub("", seed)
            if stripped != seed:
                for cand in _mmgbsa_sdf_base_candidates(stripped, stage_dir_label):
                    candidates.add(cand)

    cleaned: set[str] = set()
    for cand in candidates:
        cleaned.add(re.sub(r"_+$", "", cand))
    return sorted([c for c in cleaned if c])


def _mmgbsa_canonical_ligand_id(
    stem: str,
    stage_dir_label: str,
    regexes: List[re.Pattern],
    logger: logging.Logger,
) -> str:
    candidates = _mmgbsa_pose_candidates(stem, stage_dir_label, regexes)
    if not candidates:
        return _mmgbsa_normalize_ligand_base(stem)

    min_len = min(len(cand) for cand in candidates)
    shortest = sorted([cand for cand in candidates if len(cand) == min_len])
    chosen = shortest[0]
    if len(shortest) > 1:
        logger.warning(
            "[mmgbsa.pipeline] selection=ligand_id_collision stem=%s candidates=%s chosen=%s",
            stem,
            shortest,
            chosen,
        )
    return chosen


def _mmgbsa_group_pose_sdfs(
    sdfs: List[Path],
    stage_dir_label: str,
    regexes: List[re.Pattern],
    logger: logging.Logger,
) -> Dict[str, List[Path]]:
    grouped: Dict[str, List[Path]] = {}
    for sdf in sdfs:
        ligand_id = _mmgbsa_canonical_ligand_id(sdf.stem, stage_dir_label, regexes, logger)
        if not ligand_id:
            ligand_id = sdf.stem
        grouped.setdefault(ligand_id, []).append(sdf)
    return grouped


def _mmgbsa_pose_numeric_rank(stem: str) -> Optional[int]:
    primary = re.search(r"(?:pose|rank|conf|model|p)(\d+)$", stem, flags=re.IGNORECASE)
    if primary:
        try:
            return int(primary.group(1))
        except Exception:
            return None

    trailing = re.search(r"(\d+)$", stem)
    if trailing:
        try:
            return int(trailing.group(1))
        except Exception:
            return None
    return None


def _mmgbsa_pose_sort_key(path: Path, mode: str) -> Tuple[int, object, str]:
    stem = path.stem.lower()
    if mode == "name_numeric":
        rank = _mmgbsa_pose_numeric_rank(path.stem)
        if rank is not None:
            return (0, rank, stem)
        return (1, stem, stem)
    return (0, stem, stem)


def _mmgbsa_index_sdfs(sdfs: List[Path], stage_dir_label: str) -> Dict[str, Path]:
    mapping: Dict[str, Path] = {}
    for sdf in sdfs:
        for key in _mmgbsa_sdf_base_candidates(sdf.stem, stage_dir_label):
            if key and key not in mapping:
                mapping[key] = sdf
    return mapping


def _mmgbsa_load_reranked_bases(csv_path: Path, logger: logging.Logger) -> List[str]:
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return []

    rows: List[Tuple[Optional[int], int, str]] = []
    try:
        with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle)
            for idx, row in enumerate(reader):
                lig_base = row.get("ligand_base") or _mmgbsa_normalize_ligand_base(row.get("ligand", ""))
                if not lig_base:
                    continue
                rank_val = None
                rank_raw = str(row.get("final_rank", "")).strip()
                if rank_raw:
                    try:
                        rank_val = int(float(rank_raw))
                    except Exception:
                        rank_val = None
                rows.append((rank_val, idx, lig_base))
    except Exception as exc:
        logger.warning(
            "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=csv_read_error path=%s err=%s",
            csv_path,
            exc,
        )
        return []

    if not rows:
        return []

    rows.sort(key=lambda x: (x[0] if x[0] is not None else 1_000_000_000, x[1]))
    seen: set[str] = set()
    ordered: List[str] = []
    for _, _, base in rows:
        if base in seen:
            continue
        seen.add(base)
        ordered.append(base)
    return ordered


def _mmgbsa_select_sdfs(
    stage_dir: Path,
    max_ligands: Optional[int],
    poses_per_ligand: int,
    pose_sort_mode: str,
    pose_group_regexes: List[re.Pattern],
    reranked_csv: Optional[Path] = None,
    top_pct: Optional[float] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[Path], bool, bool, List[str]]:
    sdfs = sorted(stage_dir.glob("*.sdf"))
    if not sdfs:
        return [], False, False, []

    if logger is None:
        logger = logging.getLogger("mmgbsa.pipeline")

    used_reranked = False
    used_actives = False
    stage_label = stage_dir.name
    pose_groups = _mmgbsa_group_pose_sdfs(sdfs, stage_label, pose_group_regexes, logger)
    selected_group_map = pose_groups
    selected_ligands: List[str] = []

    if reranked_csv is not None and top_pct is not None and top_pct > 0:
        if not reranked_csv.exists():
            logger.warning(
                "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=missing_reranked_csv path=%s",
                reranked_csv,
            )
        else:
            ordered_bases = _mmgbsa_load_reranked_bases(reranked_csv, logger)
            if ordered_bases:
                pct_val = float(top_pct)
                if pct_val > 100.0:
                    pct_val = 100.0
                if pct_val <= 0.0:
                    pct_val = 0.0
                ordered_ligands: List[str] = []
                seen: set[str] = set()
                for base in ordered_bases:
                    ligand_id = _mmgbsa_canonical_ligand_id(base, stage_label, pose_group_regexes, logger)
                    if ligand_id in seen or ligand_id not in pose_groups:
                        continue
                    seen.add(ligand_id)
                    ordered_ligands.append(ligand_id)

                total = len(ordered_ligands)
                target = max(1, int(math.ceil(total * pct_val / 100.0))) if total else 0
                if ordered_ligands and target:
                    selected_ligands = ordered_ligands[:target]
                    used_reranked = True
                    logger.info(
                        "[mmgbsa.pipeline] selection=reranked_top_pct pct=%.3g total=%d target=%d selected=%d path=%s",
                        pct_val,
                        total,
                        target,
                        len(selected_ligands),
                        reranked_csv,
                    )
                else:
                    logger.warning(
                        "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=no_matching_sdfs path=%s stage_dir=%s",
                        reranked_csv,
                        stage_dir,
                    )
            else:
                logger.warning(
                    "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=empty_csv path=%s",
                    reranked_csv,
                )

    if not selected_ligands:
        actives = [p for p in sdfs if "actives_final" in p.name]
        if actives:
            used_actives = True
            selected_group_map = _mmgbsa_group_pose_sdfs(actives, stage_label, pose_group_regexes, logger)
            selected_ligands = sorted(selected_group_map.keys())
        else:
            selected_group_map = pose_groups
            selected_ligands = sorted(pose_groups.keys())

    if max_ligands is not None and max_ligands > 0:
        selected_ligands = selected_ligands[:max_ligands]

    selected_sdfs: List[Path] = []
    for ligand_id in selected_ligands:
        group = selected_group_map.get(ligand_id, [])
        if not group:
            continue
        sorted_group = sorted(group, key=lambda p: _mmgbsa_pose_sort_key(p, pose_sort_mode))
        if poses_per_ligand > 0:
            sorted_group = sorted_group[:poses_per_ligand]
        selected_sdfs.extend(sorted_group)

    logger.info(
        "[mmgbsa.pipeline] selection=poses stage_dir=%s ligand_ids=%d poses_per_ligand=%d poses_selected=%d reranked=%s actives_only=%s",
        stage_label,
        len(selected_ligands),
        poses_per_ligand,
        len(selected_sdfs),
        used_reranked,
        used_actives,
    )
    return selected_sdfs, used_reranked, used_actives, selected_ligands


def _mmgbsa_parse_delta_total(csv_path: Path) -> Optional[float]:
    if not csv_path.exists():
        return None
    try:
        lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None

    in_delta = False
    header = None
    frame_values: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("DELTA Energy Terms"):
            in_delta = True
            header = None
            continue
        if not in_delta:
            continue
        if header is None:
            header = [t.strip() for t in stripped.split(",")]
            continue
        values = [t.strip() for t in stripped.split(",")]
        if not header:
            return None
        try:
            idx = header.index("DELTA TOTAL")
        except ValueError:
            return None
        if idx >= len(values):
            continue
        try:
            frame_values.append(float(values[idx]))
        except Exception:
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
    except Exception:
        return result

    in_delta = False
    header = None
    frames: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("DELTA Energy Terms"):
            in_delta = True
            header = None
            continue
        if not in_delta:
            continue
        if header is None:
            header = [t.strip() for t in stripped.split(",")]
            continue
        values = [t.strip() for t in stripped.split(",")]
        if not header:
            break
        try:
            idx = header.index("DELTA TOTAL")
        except ValueError:
            break
        if idx >= len(values):
            continue
        try:
            frames.append(float(values[idx]))
        except Exception:
            continue

    if not frames:
        return result

    result["frame_values"] = frames
    result["n_frames"] = len(frames)
    try:
        result["frame_mean"] = float(statistics.mean(frames))
    except Exception:
        result["frame_mean"] = None
    try:
        result["frame_median"] = float(statistics.median(frames))
    except Exception:
        result["frame_median"] = None
    try:
        result["frame_sd"] = float(statistics.pstdev(frames))
    except Exception:
        result["frame_sd"] = None
    try:
        result["frame_min"] = float(min(frames))
        result["frame_max"] = float(max(frames))
    except Exception:
        pass
    return result


def _mmgbsa_write_replicate_summary(summary_path: Path, rows: List[Dict[str, object]], mean: float, sd: float) -> None:
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
    seeds = _mmgbsa_five_rep_seeds(run_id, pdb_id, variant_dir, ph_label, stage_dir, ligand_stem)
    out_dat_name = str(cfg.get("MMGBSA_MMPBSA_OUT_DAT", "FINAL_RESULTS_MMPBSA.dat") or "FINAL_RESULTS_MMPBSA.dat")
    out_csv_name = str(cfg.get("MMGBSA_MMPBSA_OUT_CSV", "FINAL_RESULTS_MMPBSA.csv") or "FINAL_RESULTS_MMPBSA.csv")
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
        rep_rows.append({"replicate": 1, "seed": seeds[0], "score": score if score is not None else "", "ok": rep_ok, "notes": ""})
        if not rep_ok:
            _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=0.0, sd=0.0)
            return {"ok": False, "notes": "rep1_failed", "replicates_ok": 0, "replicates_total": 5}
        scores.append(float(score))
        # clone outputs for other reps
        for idx in range(2, 6):
            clone_dir = out_dir / f"rep{idx}" / "mmpbsa"
            clone_dir.mkdir(parents=True, exist_ok=True)
            for src in (rep_csv, Path(mm_res.get("out_dat", "")), Path(mm_res.get("log_path", "")), Path(mm_res.get("input_path", ""))):
                if not src:
                    continue
                if not src.exists():
                    continue
                dest = clone_dir / src.name
                if dest.exists() and not force:
                    continue
                shutil.copy2(src, dest)
            rep_rows.append({"replicate": idx, "seed": seeds[idx - 1], "score": score, "ok": True, "notes": "cloned"})
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
            rep_rows.append({"replicate": rep_idx, "seed": seed, "score": score if score is not None else "", "ok": rep_ok, "notes": notes})
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
    except Exception:
        sd_val = 0.0

    _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=mean_val, sd=sd_val)
    agg_paths = write_aggregated_mmpbsa_results(work_dir=out_dir, mean_score=mean_val, std_score=sd_val, cfg=cfg, force=force)

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
    except Exception as exc:
        if logger:
            logger.warning(
                "[mmgbsa.pipeline] aggregate_failed context=%s method=%s err=%s",
                context,
                method_norm,
                exc,
            )
    return None


def _mmgbsa_frame_aggregate(csv_path: Path, cfg: Mapping[str, Any], logger: logging.Logger) -> Dict[str, object]:
    parsed = _mmgbsa_parse_delta_frames(csv_path)
    frames: List[float] = list(parsed.get("frame_values") or [])
    total_frames = len(frames)
    try:
        frame_limit = int(cfg.get("MMGBSA_FRAME_LIMIT", 0))
    except Exception:
        frame_limit = 0
    if frame_limit > 0 and frames:
        frames = frames[:frame_limit]
    frames_used = len(frames)
    frame_method = str(cfg.get("MMGBSA_FRAME_AGG", "mean") or "mean").strip().lower()

    agg_delta = _mmgbsa_apply_agg(frames, frame_method, logger=logger, context="frame")

    subset_mean = None
    subset_median = None
    subset_sd = None
    try:
        subset_mean = float(statistics.mean(frames))
    except Exception:
        subset_mean = None
    try:
        subset_median = float(statistics.median(frames))
    except Exception:
        subset_median = None
    try:
        subset_sd = float(statistics.pstdev(frames))
    except Exception:
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
        "notes": "" if agg_delta is not None else "missing_delta",
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
        "rep_agg_method",
        "replicates_ok",
        "replicates_total",
    ]
    with summary_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", restval="")
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
        logger.info("[mmgbsa.pipeline] aggregate=skip reason=missing_summary path=%s", summary_path)
        return

    groups: Dict[Tuple[str, str], Dict[str, object]] = {}
    try:
        with summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                stage_dir = str(row.get("stage_dir") or stage_dir_label or "").strip()
                ligand_stem = str(row.get("ligand_stem") or "").strip()
                if not ligand_stem:
                    continue
                ligand_id = _mmgbsa_canonical_ligand_id(ligand_stem, stage_dir, pose_group_regexes, logger)
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


def _mmgbsa_md_seed_list(cfg: Mapping[str, Any], n_reps: int, logger: logging.Logger) -> List[int]:
    try:
        base_seed = int(cfg.get("MMGBSA_MD_BASE_SEED", 12345))
    except Exception:
        base_seed = 12345

    mode = str(cfg.get("MMGBSA_MD_SEED_MODE", "increment") or "increment").strip().lower()
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


def _mmgbsa_five_rep_seeds(run_id: str, pdb_id: str, variant: str, ph_label: str, stage_dir: str, ligand: str) -> List[int]:
    context = f"{run_id}|{pdb_id}|{variant}|{ph_label}|{stage_dir}|{ligand}"
    base_seed = int(hashlib.md5(context.encode("utf-8")).hexdigest()[:8], 16)
    seeds: List[int] = []
    for idx in range(5):
        seeds.append(max(1, base_seed + idx * 10007))
    return seeds


def _mmgbsa_md_aggregate(deltas: List[float], method: str) -> Optional[float]:
    return _mmgbsa_apply_agg(deltas, method, logger=logging.getLogger("mmgbsa.pipeline"), context="replicate")


_MMGBSA_DEPRECATED_MD_KEYS = ("MMGBSA_MD_RUN", "MMGBSA_TRAJ_MODE")
_MMGBSA_DEPRECATED_MD_WARNED = False


def _mmgbsa_effective_md_config(cfg: Mapping[str, Any], logger: Optional[logging.Logger]) -> Dict[str, object]:
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
        md_enabled = traj_mode == "IMPLICIT_MD" or _to_bool(cfg.get("MMGBSA_MD_RUN", False))

    if deprecated_keys and not _MMGBSA_DEPRECATED_MD_WARNED and logger is not None:
        logger.warning(
            "[mmgbsa.pipeline] deprecated_keys=%s msg=Deprecated MMGBSA keys detected; please use MMGBSA_MD_ENABLED",
            ",".join(sorted(deprecated_keys)),
        )
        _MMGBSA_DEPRECATED_MD_WARNED = True

    return {
        "md_enabled": md_enabled,
        "effective_mode": "IMPLICIT_MD" if md_enabled else "ONEFRAME",
        "deprecated_keys": deprecated_keys,
    }


def _mmgbsa_apply_env_overrides(cfg: Mapping[str, Any], keys: Sequence[str]) -> Dict[str, Any]:
    updated = dict(cfg)
    for key in keys:
        value = os.environ.get(key)
        if value is None or str(value).strip() == "":
            continue
        updated[key] = value
    return updated


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
            "MMGBSA_TRAJ_MODE",
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES",
            "MMGBSA_MD_ENABLED",
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
            "MMGBSA_FRAME_AGG",
            "MMGBSA_FRAME_LIMIT",
        ],
    )
    if not _to_bool(cfg.get("MMGBSA_ENABLED", False)):
        logger.info("[mmgbsa.pipeline] action=skip reason=disabled")
        return

    variant_dir = (variant_token or "legacy").upper()
    post_root = Path(cfg.get("OVERALL_DIR", ".")) / "post_docked" / run_id / pdb_id
    post_base = post_root if legacy_mode else post_root / variant_dir
    if not post_base.exists():
        logger.info("[mmgbsa.pipeline] action=skip reason=missing_post_docked path=%s", post_base)
        return

    try:
        max_ligands_val = int(cfg.get("MMGBSA_MAX_LIGANDS", 1))
    except Exception:
        max_ligands_val = 1
    if max_ligands_val <= 0:
        max_ligands_val = None

    try:
        top_pct_val = float(cfg.get("MMGBSA_RERANKED_TOP_PCT", 0))
    except Exception:
        top_pct_val = 0.0
    if top_pct_val <= 0:
        top_pct_val = None

    try:
        poses_per_ligand = int(cfg.get("MMGBSA_POSES_PER_LIGAND", 1))
    except Exception:
        poses_per_ligand = 1
    if poses_per_ligand <= 0:
        poses_per_ligand = 1

    pose_sort_mode = str(cfg.get("MMGBSA_POSE_SORT_MODE", "name_numeric") or "name_numeric").strip().lower()
    pose_group_regexes = _mmgbsa_pose_group_regexes(cfg, logger)

    agg_enabled = _to_bool(cfg.get("MMGBSA_AGGREGATE_PER_LIGAND", False))
    agg_method = str(cfg.get("MMGBSA_AGG_METHOD", "min") or "min").strip().lower()
    agg_output_name = str(cfg.get("MMGBSA_AGG_OUTPUT_CSV", "mmgbsa_pose_aggregate.csv") or "mmgbsa_pose_aggregate.csv")

    md_cfg = _mmgbsa_effective_md_config(cfg, logger)
    md_enabled = _to_bool(md_cfg.get("md_enabled", False))
    effective_mode = str(md_cfg.get("effective_mode", "ONEFRAME") or "ONEFRAME")
    five_reps = _to_bool(cfg.get("MD_FIVE_REPLICATE", False))
    try:
        md_reps = int(cfg.get("MMGBSA_MD_NREPLICATES", 1))
    except Exception:
        md_reps = 1
    if md_reps < 1:
        md_reps = 1
    md_rep_agg = str(cfg.get("MMGBSA_MD_REP_AGG", "mean") or "mean").strip().lower()
    md_copy_best = _to_bool(cfg.get("MMGBSA_MD_COPY_BEST_REPLICATE", False))
    frame_agg_method = str(cfg.get("MMGBSA_FRAME_AGG", "mean") or "mean").strip().lower()
    try:
        frame_limit = int(cfg.get("MMGBSA_FRAME_LIMIT", 0))
    except Exception:
        frame_limit = 0
    if frame_limit < 0:
        frame_limit = 0

    force = _to_bool(cfg.get("MMGBSA_FORCE", False))
    strict = _to_bool(cfg.get("MMGBSA_STRICT", False))
    stage_dir_name = str(cfg.get("MMGBSA_INPUT_STAGE_DIR", "stage1") or "stage1")
    strip_all_h_for_leap = _to_bool(cfg.get("MMGBSA_TLEAP_STRIP_ALL_H", True))
    map_hoh_to_wat = _to_bool(cfg.get("MMGBSA_TLEAP_MAP_HOH_TO_WAT", True))
    water_model = str(cfg.get("MMGBSA_TLEAP_WATER_MODEL", "tip3p") or "tip3p").strip().lower()
    if water_model not in {"tip3p"}:
        raise ValueError(f"MMGBSA_TLEAP_WATER_MODEL supports tip3p only (got {water_model})")
    mmgbsa_logger = ppm._get_logger()

    test_override: Dict[str, object] = {}
    if test_mode != "off" and top_pct_val is None:
        test_stage = stage_dir_name
        candidate_dir = post_base / "pH7_0" / test_stage
        if candidate_dir.is_dir():
            actives = sorted(candidate_dir.glob("actives_final*.sdf"))
            if actives:
                test_override = {"ph_dir": candidate_dir.parent, "sdfs": [actives[0]], "stage_dir": test_stage}
            else:
                sdfs = sorted(candidate_dir.glob("*.sdf"))
                if sdfs:
                    test_override = {"ph_dir": candidate_dir.parent, "sdfs": [sdfs[0]], "stage_dir": test_stage}

    if test_override:
        ph_dirs = [test_override["ph_dir"]]
    else:
        ph_dirs = sorted([p for p in post_base.iterdir() if p.is_dir()])

    if not ph_dirs:
        logger.info("[mmgbsa.pipeline] action=skip reason=no_ph_dirs path=%s", post_base)
        return

    amber_prefix = cfg.get("MMGBSA_AMBERTOOLS_PREFIX") or cfg.get("AMBERTOOLS_PREFIX")
    amber_prefix = str(amber_prefix).strip() if amber_prefix else None

    for ph_dir in ph_dirs:
        ph_label = ph_dir.name
        if test_override:
            sdfs = list(test_override["sdfs"])
            stage_dir_label = str(test_override["stage_dir"])
            pose_groups = _mmgbsa_group_pose_sdfs(sdfs, stage_dir_label, pose_group_regexes, logger)
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
            logger.warning(
                "[mmgbsa.pipeline] action=skip reason=missing_center pdb=%s variant=%s ph=%s cfg=%s",
                pdb_id,
                variant_dir,
                ph_label,
                cfg_path,
            )
            continue

        receptor_candidates = [
            ph_dir / stage_dir_label / "receptor.pdb",
            ph_dir / "receptor.pdb",
        ]
        receptor_pdb = next((cand for cand in receptor_candidates if cand.exists()), None)
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
            {"map_hoh_to_wat": map_hoh_to_wat, "water_model": water_model, "hoh_residues": hoh_residue_count},
            "water_mapping",
        )

        mmgbsa_dir = ph_dir / "mmgbsa"
        summary_path = mmgbsa_dir / "mmgbsa_results_summary.csv"
        work_root = mmgbsa_dir / "work"

        prep_results = prep_mmgbsa_from_sdfs(
            [str(p) for p in sdfs],
            cfg=cfg,
            max_ligands=0,
            force=force,
            amber_prefix=amber_prefix,
        )

        ligand_entries: List[dict] = []
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
            ligand_entries.append(
                {
                    "stage_dir": res.get("stage_dir", stage_dir_label),
                    "ligand_base": Path(res.get("sdf_path", "ligand")).stem,
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
            }
            replicates_total = md_reps if md_enabled else 1
            replicates_ok = 0
            rep_agg_method_used = md_rep_agg if md_enabled else "single"

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
                            "frame_agg_method": summary_frames.get("frame_agg_method", "") or "",
                            "frames_used": summary_frames.get("frames_used", "") or "",
                            "frames_total": summary_frames.get("frames_total", "") or "",
                            "frame_mean": summary_frames.get("frame_mean", "") or "",
                            "frame_sd": summary_frames.get("frame_sd", "") or "",
                            "frame_median": summary_frames.get("frame_median", "") or "",
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
                        default_traj = str(cfg.get("MMGBSA_DEFAULT_TRAJ_NAME", "mdcrd") or "mdcrd")
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
                        if results_csv and results_dat and Path(results_csv).exists() and Path(results_dat).exists():
                            frame_meta = _mmgbsa_frame_aggregate(Path(results_csv), cfg, logger)
                            delta_val = frame_meta.get("delta_total")
                            if delta_val is not None:
                                delta_total = f"{float(delta_val):.6g}"
                                ok = True
                            else:
                                notes = frame_meta.get("notes", "missing_delta") or "missing_delta"
                            summary_frames.update(
                                {
                                    "frame_agg_method": frame_meta.get("frame_agg_method", ""),
                                    "frames_used": frame_meta.get("frames_used", ""),
                                    "frame_mean": frame_meta.get("frame_mean", ""),
                                    "frame_sd": frame_meta.get("frame_sd", ""),
                                    "frames_total": frame_meta.get("frames_total", ""),
                                    "frame_median": frame_meta.get("frame_median", ""),
                                }
                            )
                            replicates_ok = 1 if ok else 0
                        else:
                            notes = "missing_outputs"
                else:
                    rep_seeds = _mmgbsa_md_seed_list(cfg, md_reps, logger)
                    rep_results: List[Dict[str, object]] = []
                    for rep_idx, seed in enumerate(rep_seeds, start=1):
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
                        }
                        rep_csv = ""
                        rep_dat = ""
                        rep_log = ""

                        if not md_result.get("ok"):
                            rep_notes = "md_failed"
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
                                elif rep_csv and rep_dat and Path(rep_csv).exists() and Path(rep_dat).exists():
                                    frame_meta = _mmgbsa_frame_aggregate(Path(rep_csv), cfg, logger)
                                    delta_val = frame_meta.get("delta_total")
                                    rep_frames.update(
                                        {
                                            "frame_agg_method": frame_meta.get("frame_agg_method", frame_agg_method),
                                            "frames_used": frame_meta.get("frames_used", ""),
                                            "frames_total": frame_meta.get("frames_total", ""),
                                            "frame_mean": frame_meta.get("frame_mean", ""),
                                            "frame_sd": frame_meta.get("frame_sd", ""),
                                            "frame_median": frame_meta.get("frame_median", ""),
                                        }
                                    )
                                    if delta_val is not None:
                                        rep_delta = float(delta_val)
                                        rep_ok = True
                                    else:
                                        rep_notes = frame_meta.get("notes", "missing_delta") or "missing_delta"
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
                                "frame_agg_method": rep_frames.get("frame_agg_method", frame_agg_method),
                                "frame_median": rep_frames.get("frame_median", ""),
                            }
                        )

                    ok_reps = [r for r in rep_results if r.get("ok")]
                    n_ok = len(ok_reps)
                    replicates_ok = n_ok
                    agg_delta_val = _mmgbsa_md_aggregate(
                        [float(r["delta_total"]) for r in ok_reps if r.get("delta_total") is not None],
                        md_rep_agg,
                    )
                    if agg_delta_val is not None:
                        delta_total = f"{agg_delta_val:.6g}"
                        ok = True
                    else:
                        ok = False

                    best_rep = None
                    if ok_reps:
                        if md_rep_agg == "min":
                            best_rep = min(ok_reps, key=lambda r: (r.get("delta_total", 0), r.get("replicate", 0)))
                        elif md_rep_agg == "median":
                            best_rep = min(
                                ok_reps,
                                key=lambda r: (abs(float(r.get("delta_total", 0)) - float(agg_delta_val or 0)), r.get("replicate", 0)),
                            )
                        else:
                            best_rep = min(
                                ok_reps,
                                key=lambda r: (abs(float(r.get("delta_total", 0)) - float(agg_delta_val or 0)), r.get("replicate", 0)),
                            )

                    if best_rep:
                        results_csv = str(best_rep.get("results_csv", ""))
                        summary_frames.update(
                            {
                                "frame_agg_method": best_rep.get("frame_agg_method", frame_agg_method),
                                "frames_used": best_rep.get("frames_used", ""),
                                "frame_mean": best_rep.get("frame_mean", ""),
                                "frame_sd": best_rep.get("frame_sd", ""),
                                "frames_total": best_rep.get("frames_total", ""),
                                "frame_median": best_rep.get("frame_median", ""),
                            }
                        )
                    notes = f"md_reps_ok={n_ok}/{len(rep_results)} rep_agg={md_rep_agg} frame_agg={frame_agg_method}"

                    replicate_summary = {
                        "ok": ok,
                        "n_reps_total": len(rep_results),
                        "n_reps_ok": n_ok,
                        "agg_method": md_rep_agg,
                        "agg_delta": agg_delta_val,
                        "best_replicate": best_rep,
                        "replicates": rep_results,
                    }
                    summary_json = out_dir / "mmgbsa_replicate_summary.json"
                    try:
                        tmp_json = summary_json.with_suffix(".json.part")
                        with tmp_json.open("w", encoding="utf-8") as handle:
                            json.dump(replicate_summary, handle, indent=2, sort_keys=True)
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
                        for filename in ("FINAL_RESULTS_MMPBSA.dat", "FINAL_RESULTS_MMPBSA.csv", "mmpbsa.log", "mmpbsa.in"):
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
                    "frame_agg_method": summary_frames.get("frame_agg_method", "") or "",
                    "frames_used": summary_frames.get("frames_used", "") or "",
                    "frames_total": summary_frames.get("frames_total", "") or "",
                    "frame_mean": summary_frames.get("frame_mean", "") or "",
                    "frame_sd": summary_frames.get("frame_sd", "") or "",
                    "frame_median": summary_frames.get("frame_median", "") or "",
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
