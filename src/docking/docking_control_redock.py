from __future__ import annotations

import csv
import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from path_router.path_router import Paths, docked_dir, receptor_file, load_ph_tags
from docking.docking_vina import emit_vina_config
from docking.run_vina import run_docking_task
from docking.docking_controls import select_center_via_control_redock
from run_manifest import emit_pocket_detection_event, PocketDetectionEvent
from druggability_orchestrator import decide_engine_policy
from docking.docking_gnina import should_run_gnina_for_target, run_gnina_for_stage
from docking.docking_ledock import (
    should_run_ledock_for_target,
    ensure_ledock_receptor,
    run_ledock_for_stage,
)
from docking.docking_dock6 import should_run_dock6_for_target, run_dock6_for_stage
from prep_docking.prep_for_ledock import ensure_mol2_for_ledock
import prep_docking.prep_dock6 as prep_dock6


def _collect_control_pdbqts(
    paths: Paths, control_stems: List[str], cfg: Dict
) -> List[Path]:
    stems_lower = {s.lower() for s in control_stems if s}
    if not stems_lower:
        return []
    roots = [paths.prepped_ligands_dir]
    out_root = cfg.get("PREPPED_LIGANDS_DIR")
    try:
        if out_root:
            roots.append(Path(out_root) / paths.pdb_id)
    except Exception:
        pass
    hits: List[Path] = []
    for root in roots:
        if not root or not Path(root).exists():
            continue
        for p in Path(root).glob("*.pdbqt"):
            stem0 = p.stem.split("_stage")[0].lower()
            if stem0 in stems_lower:
                hits.append(p)
    return hits


def _control_centers_by_ph(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    *,
    variant_token: Optional[str],
    legacy_mode: bool,
    cleaned_pdb: str,
    fallback_center: Optional[Tuple[float, float, float]],
    fallback_box: Optional[Tuple[float, float, float]],
) -> Tuple[
    Dict[Optional[str], Tuple[float, float, float]],
    Dict[Optional[str], Tuple[float, float, float]],
    Dict[Optional[str], str],
]:
    center_by_ph: Dict[Optional[str], Tuple[float, float, float]] = {}
    box_by_ph: Dict[Optional[str], Tuple[float, float, float]] = {}
    source_by_ph: Dict[Optional[str], str] = {}

    ph_tags: List[Optional[str]]
    if bool(cfg.get("PH_ENSEMBLE", False)):
        ph_tags = [
            str(tag).strip()
            for tag in (load_ph_tags(paths.pdb_id, variant=variant_token) or [])
        ]
        ph_tags = [tag for tag in ph_tags if tag]
        if not ph_tags:
            ph_tags = [None]
    else:
        ph_tags = [None]

    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))

    for ph_label in ph_tags:
        rec_pdbqt = receptor_file(
            paths.pdb_id,
            variant=variant_token,
            ph_tag=ph_label,
            legacy=legacy_mode,
        )
        ph_print = ph_label or "base"
        if not rec_pdbqt.exists():
            logger.warning(
                "[control-redock] receptor_missing pdb=%s variant=%s ph=%s path=%s",
                paths.pdb_id,
                variant_token or "legacy",
                ph_print,
                rec_pdbqt,
            )
            if fallback_center is None or fallback_box is None:
                logger.error(
                    "[control-redock] fallback_missing pdb=%s variant=%s ph=%s",
                    paths.pdb_id,
                    variant_token or "legacy",
                    ph_print,
                )
                continue
            center = tuple(float(x) for x in fallback_center)
            box_size = tuple(float(s) for s in fallback_box)
            source = "receptor_missing"
        else:
            try:
                sel_center, sel_box = select_center_via_control_redock(
                    cfg,
                    paths,
                    str(rec_pdbqt),
                    logger,
                    variant=variant_token,
                    ph_token=ph_label,
                    legacy=legacy_mode,
                )
            except Exception as _e:
                sel_center, sel_box = (None, None)
                logger.exception(
                    "[control-centers] helper errored during ctrl_redock for pdb=%s variant=%s ph=%s: %s",
                    paths.pdb_id,
                    variant_token or "legacy",
                    ph_print,
                    _e,
                )

            if sel_center is not None and sel_box is not None:
                center = tuple(float(x) for x in sel_center)
                box_size = tuple(float(s) for s in sel_box)
                source = "control"
            else:
                if fallback_center is None or fallback_box is None:
                    logger.error(
                        "[control-redock] fallback_missing pdb=%s variant=%s ph=%s",
                        paths.pdb_id,
                        variant_token or "legacy",
                        ph_print,
                    )
                    continue
                center = tuple(float(x) for x in fallback_center)
                box_size = tuple(float(s) for s in fallback_box)
                source = "fallback_no_controls"

        if source == "control":
            side = float(cfg.get("CONTROL_BOX_A", 24.0))
            box_size = (side, side, side)

        box_size = tuple(min(box_cap, float(s)) for s in box_size)
        center_by_ph[ph_label] = center
        box_by_ph[ph_label] = box_size
        source_by_ph[ph_label] = source
        logger.info(
            "[center.final] pdb=%s variant=%s ph=%s source=%s box=%s",
            paths.pdb_id,
            variant_token or "legacy",
            ph_print,
            source,
            box_size,
        )

        try:
            # Note: _cache_active_site is internal to docking.py and not extracted here.
            # If we need to call it, we should probably pass a callback or accept that
            # this side-effect (caching active site for DOCK6 reuse later) is handled elsewhere.
            # However, docking.py's version of _control_centers_by_ph called _cache_active_site.
            # Since _cache_active_site is not in the list of moved definitions, and it modifies 'cfg',
            # we can inline a simplified version or just omit the log if the main side effect is updating 'cfg'.
            # The original code:
            # _cache_active_site(cfg, paths.pdb_id, variant_token, ph_label, center, box_size, logger)
            # This function just updates cfg["_ACTIVE_SITE_CACHE"]. We can replicate the logic here safely.

            cache = cfg.setdefault("_ACTIVE_SITE_CACHE", {})
            cache_key = (
                str(paths.pdb_id).upper(),
                (variant_token or "HOLO"),
                (ph_label or "base"),
            )
            cache[cache_key] = (
                tuple(float(x) for x in center),
                tuple(float(x) for x in box_size),
            )
        except Exception:
            logger.debug(
                "[active-site.cache.store.skip] pdb=%s variant=%s ph=%s",
                paths.pdb_id,
                variant_token,
                ph_print,
            )

        try:
            c_print = tuple(round(float(x), 3) for x in center)
            b_print = tuple(round(float(x), 1) for x in box_size)
            print(
                f"[CENTER] ph={ph_print} source={source} center={c_print} box={b_print}"
            )
        except Exception:
            pass

        try:
            manifest_run_id = cfg.get("RUN_ID")
            logger.debug(
                "[run-manifest.pocket_detection.call] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%r box=%r",
                manifest_run_id,
                paths.pdb_id,
                variant_token or "legacy",
                ph_print,
                source,
                center,
                box_size,
            )
            if manifest_run_id:
                emit_pocket_detection_event(
                    cfg,
                    PocketDetectionEvent(
                        run_id=str(manifest_run_id),
                        pdb_id=paths.pdb_id,
                        variant_label=variant_token or "legacy",
                        ph_tag=ph_label,
                        method=source,
                        center=center,
                        box_size=box_size,
                    ),
                )
            else:
                logger.debug(
                    "[run-manifest.pocket_detection.skip] no RUN_ID for pdb=%s variant=%s ph=%s",
                    paths.pdb_id,
                    variant_token or "legacy",
                    ph_print,
                )
        except Exception:
            logger.warning(
                "[run-manifest.pocket_detection.error] pdb=%s variant=%s ph=%s",
                paths.pdb_id,
                variant_token or "legacy",
                ph_print,
                exc_info=True,
            )

    return center_by_ph, box_by_ph, source_by_ph


def _select_control_engines(
    cfg: Dict,
    use_gnina: bool,
    use_ledock: bool,
    use_dock6: bool,
) -> List[str]:
    engines: List[str] = ["vina"]
    if not bool(cfg.get("CONTROL_CONSENSUS", False)):
        return engines
    if use_gnina:
        engines.append("gnina")
    if use_ledock:
        engines.append("ledock")
    if use_dock6:
        engines.append("dock6")
    return engines


def run_control_docking_multi_engine(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    *,
    variant_token: Optional[str],
    legacy_mode: bool,
    control_lookup: Dict[str, Path],
    control_ligands: List[Path],
    center_by_ph: Dict[Optional[str], Tuple[float, float, float]],
    box_by_ph: Dict[Optional[str], Tuple[float, float, float]],
) -> None:
    control_paths = [Path(p) for p in control_ligands]
    if not control_paths:
        logger.warning(
            "[control-multi] no_control_ligands_found; skipping multi-engine control docking"
        )
        return

    stage_info = {
        "exhaustiveness": int(cfg.get("CTRL_REDOCK_EXHAUSTIVENESS", 64)),
        "num_modes": int(cfg.get("CTRL_REDOCK_NMODES", 9)),
        "energy_range": float(cfg.get("CTRL_REDOCK_ENERGY_RANGE", 6)),
        "verbosity": int(cfg.get("VINA_VERBOSITY", 0)),
    }
    if cfg.get("FAST_MODE"):
        stage_info["exhaustiveness"] = 1
        stage_info["num_modes"] = 1

    ph_labels = list(center_by_ph.keys()) if center_by_ph else [None]
    for ph_label in ph_labels:
        center = center_by_ph.get(ph_label) or center_by_ph.get(None)
        box_size = box_by_ph.get(ph_label) if box_by_ph else None
        if box_size is None:
            box_size = box_by_ph.get(None) if box_by_ph else None
        ph_print = ph_label or "base"
        if center is None or box_size is None:
            logger.warning(
                "[control-multi] skip ph=%s reason=missing_center_or_box", ph_print
            )
            continue

        rec_path = receptor_file(
            paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode
        )
        ctrl_root = (
            docked_dir(
                paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode
            )
            / "ctrl_redock"
        )
        ctrl_root.mkdir(parents=True, exist_ok=True)
        results: List[Dict[str, object]] = []

        ledock_capable = should_run_ledock_for_target(cfg)
        dock6_capable = should_run_dock6_for_target(cfg)
        policy = None
        try:
            policy = decide_engine_policy(
                cfg=cfg,
                pdb_id=paths.pdb_id,
                variant=variant_token or None,
                ph_label=ph_label,
                center=center,
                ledock_enabled=ledock_capable,
                dock6_enabled=dock6_capable,
                logger=logger,
            )
        except Exception as exc:
            logger.debug(
                "[control-multi.policy.error] pdb=%s variant=%s ph=%s reason=%s",
                paths.pdb_id,
                variant_token or "legacy",
                ph_print,
                exc,
            )

        gnina_enabled = bool(should_run_gnina_for_target(None, cfg))
        use_gnina = gnina_enabled
        if policy:
            use_gnina = use_gnina or bool(policy.use_gnina)
        use_ledock = ledock_capable and (bool(policy.use_ledock) if policy else True)
        use_dock6 = dock6_capable and (bool(policy.use_dock6) if policy else True)
        gnina_exe = (cfg.get("GNINA_EXE") or "").strip()
        if use_gnina and not gnina_exe:
            logger.info(
                "[control-gnina] disabled ph=%s reason=missing_gnina_exe", ph_print
            )
            use_gnina = False

        engines = _select_control_engines(cfg, use_gnina, use_ledock, use_dock6)

        logger.info("[control-multi] ph=%s engines=%s", ph_print, ",".join(engines))
        if any(engine in ("ledock", "dock6") for engine in engines):
            try:
                ensure_mol2_for_ledock(cfg, [str(p) for p in control_paths], logger)
                logger.info(
                    "[control-multi.prep] ledock_mol2_prep n=%d", len(control_paths)
                )
            except Exception as exc:
                logger.warning("[control-multi.prep] ledock_mol2_failed reason=%s", exc)

        ledock_receptor: Optional[Path] = None
        if "ledock" in engines:
            try:
                ledock_receptor = ensure_ledock_receptor(
                    cfg, paths.pdb_id, variant_token, ph_label, logger
                )
            except Exception as exc:
                logger.warning(
                    "[control-ledock] prep_failed pdb=%s variant=%s ph=%s reason=%s",
                    paths.pdb_id,
                    variant_token or "legacy",
                    ph_print,
                    exc,
                )

        dock6_ready = True
        if "dock6" in engines:
            try:
                prep_dock6.ensure_dock6_site(
                    cfg=cfg,
                    pdb_id=paths.pdb_id,
                    variant=variant_token or None,
                    ph_label=ph_label,
                    logger=logger,
                )
            except Exception as exc:
                dock6_ready = False
                logger.warning(
                    "[control-dock6.prep.error] pdb=%s variant=%s ph=%s reason=%s",
                    paths.pdb_id,
                    variant_token or "legacy",
                    ph_print,
                    exc,
                )

        # --- Vina control docking
        engine_root = ctrl_root / "vina"
        engine_root.mkdir(parents=True, exist_ok=True)
        vina_exe = cfg.get("VINA_EXE") or cfg.get("VINA_PATH") or "vina"
        vina_stage_name = "ctrl_redock_vina"
        vina_threads = int(cfg.get("THREADS_PER_VINA_CTRL", 8))
        rec_exists = rec_path and Path(rec_path).exists()
        if not rec_exists:
            logger.warning(
                "[control-vina] receptor_missing ph=%s path=%s", ph_print, rec_path
            )
        for lig_path in control_paths:
            lig_base = Path(lig_path).stem.split("_stage")[0]
            pose_path = engine_root / f"{lig_base}_{vina_stage_name}.pdbqt"
            status = "ok"
            score: Optional[float] = None
            try:
                conf_path, out_path = emit_vina_config(
                    cfg,
                    paths.pdb_id,
                    str(rec_path),
                    center,
                    box_size,
                    str(lig_path),
                    vina_stage_name,
                    stage_info,
                    vina_threads,
                    logger=None,
                    variant=variant_token,
                    ph_token=ph_label,
                    legacy=legacy_mode,
                )
                _, score_val = run_docking_task(
                    vina_exe, str(conf_path), str(lig_path), str(out_path)
                )
                score = score_val if isinstance(score_val, (int, float)) else None
                src_pose = Path(out_path)
                if src_pose.exists():
                    try:
                        shutil.copy2(src_pose, pose_path)
                    except Exception:
                        pose_path = src_pose
                status = (
                    "ok"
                    if score is not None
                    else ("missing_inputs" if not rec_exists else "no_score")
                )
            except Exception as exc:
                status = "missing_inputs" if not rec_exists else "failed"
                try:
                    pose_path.write_text(
                        f"vina control redock failed: {exc}\n", encoding="utf-8"
                    )
                except Exception:
                    pass
                logger.warning(
                    "[control-vina] run_failed pdb=%s variant=%s ph=%s lig=%s reason=%s",
                    paths.pdb_id,
                    variant_token or "legacy",
                    ph_print,
                    lig_base,
                    exc,
                )

            results.append(
                {
                    "engine": "vina",
                    "control": lig_base,
                    "ph": ph_print,
                    "status": status,
                    "score": score,
                    "pose_path": str(pose_path),
                }
            )

        # --- GNINA control docking
        if "gnina" in engines:
            engine_root = ctrl_root / "gnina"
            engine_root.mkdir(parents=True, exist_ok=True)
            gn_stage_name = "ctrl_redock"
            gn_stage_info = dict(stage_info)
            gn_stage_info["name"] = gn_stage_name
            gnina_scores: Dict[str, Optional[float]] = {}
            gnina_metrics: Dict[str, Dict[str, Any]] = {}
            if not rec_exists:
                logger.warning(
                    "[control-gnina] skip ph=%s reason=missing_receptor path=%s",
                    ph_print,
                    rec_path,
                )
            else:
                try:
                    gnina_scores, gnina_metrics, _ = run_gnina_for_stage(
                        cfg=cfg,
                        paths=paths,
                        pdb_id=paths.pdb_id,
                        variant=variant_token,
                        ph_label=ph_label,
                        stage_name=gn_stage_name,
                        stage_info=gn_stage_info,
                        ligands=[str(p) for p in control_paths],
                        center=center,
                        box_size=box_size,
                        logger=logger,
                        receptor_pdbqt=str(rec_path),
                        control_lookup=control_lookup,
                    )
                except Exception as exc:
                    logger.warning(
                        "[control-gnina] run_failed pdb=%s variant=%s ph=%s reason=%s",
                        paths.pdb_id,
                        variant_token or "legacy",
                        ph_print,
                        exc,
                    )

            gn_stage_dir = paths.docked_stage_dir(
                variant_token, f"gnina_{gn_stage_name}", ph_label
            )
            for lig_path in control_paths:
                lig_base = Path(lig_path).stem.split("_stage")[0]
                pose_path = (
                    gn_stage_dir / f"{Path(lig_path).stem}_gnina_{gn_stage_name}.pdbqt"
                )
                if pose_path.exists():
                    dest = engine_root / pose_path.name
                    try:
                        shutil.copy2(pose_path, dest)
                        pose_path = dest
                    except Exception:
                        pass
                score = None
                metrics_rec = (
                    gnina_metrics.get(str(lig_path))
                    or gnina_metrics.get(Path(lig_path).as_posix())
                    or {}
                )
                if isinstance(metrics_rec, dict):
                    score = (
                        metrics_rec.get("gnina_primary_score")
                        or metrics_rec.get("cnn_affinity_pK")
                        or metrics_rec.get("minimized_affinity_kcal")
                    )
                if score is None:
                    sc_raw = gnina_scores.get(str(lig_path)) or gnina_scores.get(
                        lig_path
                    )
                    if isinstance(sc_raw, (int, float)):
                        score = sc_raw
                status = (
                    "missing_inputs"
                    if not rec_exists
                    else (
                        "ok" if pose_path.exists() and score is not None else "failed"
                    )
                )
                results.append(
                    {
                        "engine": "gnina",
                        "control": lig_base,
                        "ph": ph_print,
                        "status": status,
                        "score": score,
                        "pose_path": str(pose_path),
                    }
                )

        # --- LeDock control docking
        if "ledock" in engines:
            engine_root = ctrl_root / "ledock"
            engine_root.mkdir(parents=True, exist_ok=True)
            if not ledock_receptor or not Path(ledock_receptor).exists():
                for lig_path in control_paths:
                    lig_base = Path(lig_path).stem.split("_stage")[0]
                    marker = engine_root / f"{lig_base}.dok"
                    try:
                        marker.write_text("missing ledock receptor\n", encoding="utf-8")
                    except Exception:
                        pass
                    results.append(
                        {
                            "engine": "ledock",
                            "control": lig_base,
                            "ph": ph_print,
                            "status": "missing_inputs",
                            "score": None,
                            "pose_path": str(marker),
                        }
                    )
            else:
                try:
                    ledock_scores, ledock_metrics, _ = run_ledock_for_stage(
                        cfg=cfg,
                        paths=paths,
                        pdb_id=paths.pdb_id,
                        variant=variant_token,
                        ph_label=ph_label,
                        stage_name="ctrl_redock_stage1",
                        stage_info={"key": "stage1", **stage_info},
                        ligands=control_paths,
                        center=center,
                        box_size=box_size,
                        logger=logger,
                        receptor_pdb=ledock_receptor,
                        skip_completion=True,
                        output_root=engine_root,
                    )
                except Exception as exc:
                    ledock_scores = {}
                    ledock_metrics = {}
                    logger.warning(
                        "[control-ledock] run_failed pdb=%s variant=%s ph=%s reason=%s",
                        paths.pdb_id,
                        variant_token or "legacy",
                        ph_print,
                        exc,
                    )

                scores_map = {str(k): v for k, v in (ledock_scores or {}).items()}
                metrics_map = {str(k): v for k, v in (ledock_metrics or {}).items()}
                for lig_path in control_paths:
                    lig_base = Path(lig_path).stem.split("_stage")[0]
                    pose_path = engine_root / f"{Path(lig_path).stem}.dok"
                    score = scores_map.get(str(lig_path))
                    metrics_rec = metrics_map.get(str(lig_path), {})
                    if score is None and isinstance(metrics_rec, dict):
                        score = metrics_rec.get("best_score_kcal")
                    status = (
                        "ok"
                        if pose_path.exists() and isinstance(score, (int, float))
                        else "failed"
                    )
                    results.append(
                        {
                            "engine": "ledock",
                            "control": lig_base,
                            "ph": ph_print,
                            "status": status,
                            "score": score,
                            "pose_path": str(pose_path),
                        }
                    )

        # --- DOCK6 control docking
        if "dock6" in engines:
            engine_root = ctrl_root / "dock6"
            engine_root.mkdir(parents=True, exist_ok=True)
            dock6_scores: Dict[Path, Optional[float]] = {}
            dock6_metrics: Dict[Path, Dict[str, Any]] = {}
            if dock6_ready:
                try:
                    dock6_scores, dock6_metrics = run_dock6_for_stage(
                        cfg=cfg,
                        paths=paths,
                        pdb_id=paths.pdb_id,
                        variant=variant_token,
                        ph_label=ph_label,
                        stage_name="ctrl_redock_stage1",
                        stage_info={"key": "stage1", **stage_info},
                        ligands=control_paths,
                        center=center,
                        box_size=box_size,
                        logger=logger,
                        receptor_pdb=None,
                        skip_completion=True,
                        output_root=engine_root,
                    )
                except Exception as exc:
                    logger.warning(
                        "[control-dock6] run_failed pdb=%s variant=%s ph=%s reason=%s",
                        paths.pdb_id,
                        variant_token or "legacy",
                        ph_print,
                        exc,
                    )
            scores_map_d6 = {str(k): v for k, v in (dock6_scores or {}).items()}
            metrics_map_d6 = {str(k): v for k, v in (dock6_metrics or {}).items()}
            ranked_path = engine_root / "dock6_stage1_ranked.mol2"
            for lig_path in control_paths:
                lig_base = Path(lig_path).stem.split("_stage")[0]
                score = scores_map_d6.get(str(lig_path))
                if score is None:
                    rec = metrics_map_d6.get(str(lig_path)) or {}
                    score = rec.get("grid_score")
                pose_path = ranked_path if ranked_path.exists() else engine_root
                status = (
                    "ok"
                    if isinstance(score, (int, float)) and pose_path.exists()
                    else ("missing_inputs" if not dock6_ready else "failed")
                )
                results.append(
                    {
                        "engine": "dock6",
                        "control": lig_base,
                        "ph": ph_print,
                        "status": status,
                        "score": score,
                        "pose_path": str(pose_path),
                    }
                )

        csv_path = ctrl_root / "control_docking_engine_scores.csv"
        try:
            with csv_path.open("w", encoding="utf-8", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=["engine", "control", "score"])
                writer.writeheader()
                for row in results:
                    writer.writerow(
                        {
                            "engine": row.get("engine"),
                            "control": row.get("control"),
                            "score": row.get("score"),
                        }
                    )
            logger.info(
                "[control-multi.csv] pdb=%s variant=%s ph=%s engines=%s path=%s",
                paths.pdb_id,
                variant_token or "legacy",
                ph_print,
                ",".join(engines),
                csv_path,
            )
        except Exception as exc:
            logger.warning(
                "[control-multi.csv.error] pdb=%s variant=%s ph=%s reason=%s",
                paths.pdb_id,
                variant_token or "legacy",
                ph_print,
                exc,
            )
