# docking.py
# Docking orchestration helpers extracted from main.py

from __future__ import annotations

import logging
import os
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from logging_topics import make_protein_logger
from path_router.path_router import Paths, make_paths
from protein_functions import detect_active_site
from run_manifest import update_manifest_for_docking_overall

from docking import docking as _docking
from docking import docking_controls
from docking.docking_control_centering_phase import _phase5b_controls_and_control_redock
from docking.docking_control_redock import (
    _collect_control_pdbqts,
    run_control_docking_multi_engine,
)
from docking.docking_controls import (
    _summarize_ions_file,
    _is_readable_ref,
)
from docking.docking_ligands import (
    select_ligands_for_next,
    _coerce_test_map,
    _resolve_test_mode,
    _read_any_lig,
    compute_rmsd,
    validate_ligand,
)
from docking.library_mode import parse_test_libraries
from docking.docking_ph_ensemble_phase import _phase5_ph_ensemble_global
from docking.docking_receptor_phases import _phase2_to4_receptor_and_center
from docking.docking_stage_runner import RetryManager, run_one_stage, _map_reason_to_category
from docking.docking_subruns import (
    ProteinDockingContext,
    SubrunSpec,
    subruns_for_tokens,
    run_ligand_pipeline_subrun,
)
from docking.docking_utils import norm
from docking.docking_vina import emit_vina_config
from docking.fallback_recenter import RecenterParams
from druggability_evaluation import evaluate_druggability_for_active_site

# Bridge functions for docking helpers
_docking.emit_vina_config = emit_vina_config
_docking.norm = norm
_docking._read_any_lig = _read_any_lig
_docking.validate_ligand = validate_ligand

docking_controls._is_readable_ref = _is_readable_ref
docking_controls.compute_rmsd = compute_rmsd

_ACTIVE_SITE_CACHE_KEY = "_ACTIVE_SITE_CACHE"


def _cache_active_site(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    logger: Optional[logging.Logger] = None,
) -> None:
    """Record active-site center/box for reuse by downstream helpers (e.g., DOCK6 prep)."""
    try:
        cache = cfg.setdefault(_ACTIVE_SITE_CACHE_KEY, {})
        cache_key = (
            str(pdb_id).upper(),
            (variant_token or "HOLO"),
            (ph_label or "base"),
        )
        cache[cache_key] = (
            tuple(float(x) for x in center),
            tuple(float(x) for x in box_size),
        )
        if logger:
            logger.debug(
                "[active-site.cache.store] key=%s center=%s box=%s",
                cache_key,
                cache[cache_key][0],
                cache[cache_key][1],
            )
    except Exception:
        if logger:
            logger.debug("[active-site.cache.skip]")


def get_active_site_center_and_size(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    logger: logging.Logger,
) -> Optional[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    """
    Returns (center, size) for the active site used by Vina/LeDock.

    center: (cx, cy, cz)
    size:   (sx, sy, sz) full box lengths (same values emitted to Vina configs)

    Prefers cached values populated during docking setup; falls back to running
    detect_active_site on the cleaned receptor.
    """
    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    variant_key = variant_token or "HOLO"
    ph_key = (str(ph_label).strip() or "") or "base"

    cache = cfg.get(_ACTIVE_SITE_CACHE_KEY)
    if isinstance(cache, dict):
        cache_key = (str(pdb_id).upper(), variant_key, ph_key)
        hit = cache.get(cache_key)
        if not hit and ph_key != "base":
            hit = cache.get((str(pdb_id).upper(), variant_key, "base"))
        if hit:
            logger.info(
                "[active-site.cache.hit] pdb=%s variant=%s ph=%s center=%s box=%s",
                pdb_id,
                variant_key,
                ph_key,
                hit[0],
                hit[1],
            )
            return hit  # type: ignore[return-value]

    try:
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        cleaned_pdb = paths.receptor_cleaned_pdb(variant_token)
    except Exception as exc:
        logger.warning(
            "[active-site.helper.error] pdb=%s variant=%s ph=%s reason=%s",
            pdb_id,
            variant_key,
            ph_key,
            exc,
        )
        return None

    if not cleaned_pdb or not Path(cleaned_pdb).exists():
        logger.warning(
            "[active-site.helper.skip] reason=missing_cleaned pdb=%s variant=%s ph=%s path=%s",
            pdb_id,
            variant_key,
            ph_key,
            cleaned_pdb,
        )
        return None

    center, box_size, source = detect_active_site(cleaned_pdb)
    if center and box_size:
        center_t = tuple(float(x) for x in center)
        box_t = tuple(float(x) for x in box_size)
        _cache_active_site(cfg, pdb_id, variant_token, ph_key, center_t, box_t, logger)
        logger.info(
            "[active-site.helper.detect] pdb=%s variant=%s ph=%s source=%s center=%s box=%s",
            pdb_id,
            variant_key,
            ph_key,
            source or "activesite",
            center_t,
            box_t,
        )
        try:
            evaluate_druggability_for_active_site(
                cfg,
                pdb_id=pdb_id,
                variant=variant_token,
                ph_label=ph_key,
                center=center_t,
                box_size=box_t,
                logger=logger,
            )
        except Exception as exc:
            logger.warning(
                "[druggability.fpocket.skip] pdb=%s variant=%s ph=%s err=%s",
                pdb_id,
                variant_key,
                ph_key,
                exc,
            )
        return center_t, box_t

    logger.warning(
        "[active-site.helper.miss] pdb=%s variant=%s ph=%s source=%s",
        pdb_id,
        variant_key,
        ph_key,
        source if "source" in locals() else "unknown",
    )
    return None


def _phase0_setup_paths_and_logger(
    cfg: Dict, pdb_file: str
) -> Tuple[Paths, str, logging.Logger]:
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = re.sub(r"(?i)_cleaned$", "", base_id)
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=os.path.basename(pdb_file))

    logger = make_protein_logger(str(paths.docked_pdb_root()), pdb_id, cfg)
    logger.info(f"[paths] base_id={base_id} -> pdb_id={pdb_id}")
    logger.info(f"Processing protein: {pdb_file} (id={pdb_id})")
    return paths, pdb_id, logger


def _phase1_variant_and_ion_context(
    cfg: Dict, paths: Paths, logger: logging.Logger
) -> Tuple[str, Optional[str], str, dict, dict, bool, Path]:
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    variant_label = variant_env or "legacy"
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))

    ion_audit_root: dict = cfg.setdefault("_ION_AUDIT", {})
    pdb_audit: dict = ion_audit_root.setdefault(paths.pdb_id, {})
    clean_audit: dict = pdb_audit.setdefault("clean_counts", {})

    input_summary = _summarize_ions_file(paths.input_pdb_path)
    input_hist = str(input_summary.get("hist", "none"))
    input_error = input_summary.get("error")
    if input_error not in (None, "missing"):
        logger.warning(
            "[ions.input.counts] pdb=%s file=%s action=skip err=%s",
            paths.pdb_id,
            paths.input_pdb_path,
            input_error,
        )
    else:
        logger.info(
            "[ions.input.counts] pdb=%s file=%s present_pdb=%s",
            paths.pdb_id,
            paths.input_pdb_path,
            input_hist,
        )
    pdb_audit["input_counts"] = {
        "hist": input_hist,
        "counts": dict(input_summary.get("counts", {})),
        "metals_present": bool(input_summary.get("metals_present", False)),
        "salts_present": bool(input_summary.get("salts_present", False)),
        "file": str(paths.input_pdb_path),
        "error": input_error,
    }

    cfg["_CURRENT_VARIANT"] = variant_env
    cleaned_target = paths.receptor_cleaned_pdb(variant_token)
    receptor_target = paths.receptor_pdbqt(variant_token, None)
    logger.info(
        "[receptor.path] pdb=%s variant=%s cleaned_pdb=%s exists=%s",
        paths.pdb_id,
        variant_label,
        cleaned_target,
        cleaned_target.exists(),
    )
    logger.info(
        "[receptor.path] pdb=%s variant=%s receptor_pdbqt=%s exists=%s",
        paths.pdb_id,
        variant_label,
        receptor_target,
        receptor_target.exists(),
    )
    return (
        variant_env,
        variant_token,
        variant_label,
        pdb_audit,
        clean_audit,
        legacy_mode,
        receptor_target,
    )


def _phase6_to8_ligands_and_docking(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    pdb_id: str,
    variant_env: str,
    variant_token: Optional[str],
    variant_label: str,
    legacy_mode: bool,
    cleaned_pdb: Optional[str],
    receptor_pdbqt: Optional[str],
    center: Optional[Tuple[float, float, float]],
    box_size: Optional[Tuple[float, float, float]],
    center_by_ph: Optional[Dict[Optional[str], Tuple[float, float, float]]],
    box_by_ph: Optional[Dict[Optional[str], Tuple[float, float, float]]],
    center_source_by_ph: Optional[Dict[Optional[str], str]],
    stages: List[Dict],
    params: RecenterParams,
    control_stems: List[str],
    control_lookup: Dict[str, Path],
) -> None:
    ctx = ProteinDockingContext(
        cfg=cfg,
        paths=paths,
        logger=logger,
        pdb_id=pdb_id,
        variant_env=variant_env,
        variant_token=variant_token,
        variant_label=variant_label,
        legacy_mode=legacy_mode,
        cleaned_pdb=cleaned_pdb,
        receptor_pdbqt=receptor_pdbqt,
        center=center,
        box_size=box_size,
        center_by_ph=center_by_ph,
        box_by_ph=box_by_ph,
        center_source_by_ph=center_source_by_ph,
        stages=stages,
        recenter_params=params,
        control_stems=control_stems,
        control_lookup=control_lookup,
    )

    tokens = parse_test_libraries(cfg)
    subruns: List[SubrunSpec] = subruns_for_tokens(tokens)

    for sub in subruns:
        logger.info(
            "[subrun] tokens=%s run_mode=%r csv_prefix=%r stage_prefix=%r",
            "+".join(tokens),
            sub.run_mode,
            sub.csv_prefix,
            sub.stage_name_prefix,
        )
        run_ligand_pipeline_subrun(ctx, sub)


def process_one_protein(
    cfg: Dict, pdb_file: str, stages: List[Dict], params: RecenterParams
) -> None:
    paths, pdb_id, logger = _phase0_setup_paths_and_logger(cfg, pdb_file)
    (
        variant_env,
        variant_token,
        variant_label,
        pdb_audit,
        clean_audit,
        legacy_mode,
        receptor_target,
    ) = _phase1_variant_and_ion_context(cfg, paths, logger)
    active_ph_label = None

    (
        cleaned_pdb,
        receptor_pdbqt,
        pdb_audit,
        clean_audit,
        center,
        box_size,
        center_source,
        control_stems,
        control_lookup,
    ) = _phase2_to4_receptor_and_center(
        cfg,
        paths,
        logger,
        variant_env,
        variant_token,
        variant_label,
        pdb_audit,
        clean_audit,
        legacy_mode,
        receptor_target,
        active_ph_label,
    )
    if not receptor_pdbqt:
        return

    receptor_pdbqt, active_ph_label = _phase5_ph_ensemble_global(
        cfg,
        paths,
        logger,
        variant_env,
        variant_token,
        legacy_mode,
        cleaned_pdb,
        receptor_pdbqt,
        center,
        box_size,
    )

    (
        center,
        box_size,
        center_by_ph,
        box_by_ph,
        center_source_by_ph,
        control_stems,
        control_lookup,
    ) = _phase5b_controls_and_control_redock(
        cfg,
        paths,
        logger,
        variant_token,
        variant_label,
        legacy_mode,
        cleaned_pdb,
        receptor_pdbqt,
        center,
        box_size,
    )
    if (
        center_by_ph
        and (None not in center_by_ph)
        and center is None
        and box_size is None
    ):
        center = center_by_ph.get(active_ph_label) or center_by_ph.get(None)
        box_size = box_by_ph.get(active_ph_label) if box_by_ph else None
    if not center_by_ph and center is not None and box_size is not None:
        center_by_ph = {active_ph_label: center}
        box_by_ph = {active_ph_label: box_size}
    if (center is None or box_size is None) and not center_by_ph:
        return

    control_ligands = _collect_control_pdbqts(paths, control_stems, cfg)
    if center_by_ph:
        try:
            run_control_docking_multi_engine(
                cfg,
                paths,
                logger,
                variant_token=variant_token,
                legacy_mode=legacy_mode,
                control_lookup=control_lookup,
                control_ligands=control_ligands,
                center_by_ph=center_by_ph,
                box_by_ph=box_by_ph or {},
            )
        except Exception as exc:
            logger.warning(
                "[control-multi.error] pdb=%s variant=%s reason=%s",
                pdb_id,
                variant_label,
                exc,
                exc_info=True,
            )

    run_id = str(cfg.get("RUN_ID", "") or "")
    docking_event = "end"
    docking_error: Optional[str] = None
    docking_start = time.time()
    if run_id:
        try:
            update_manifest_for_docking_overall(
                cfg,
                run_id,
                pdb_id,
                variant_label,
                ph_tag=active_ph_label,
                event="start",
            )
        except Exception:
            logger.warning(
                "[run-manifest.docking-overall] failed to record start pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                active_ph_label if active_ph_label is not None else "base",
                exc_info=True,
            )

    try:
        _phase6_to8_ligands_and_docking(
            cfg,
            paths,
            logger,
            pdb_id,
            variant_env,
            variant_token,
            variant_label,
            legacy_mode,
            cleaned_pdb,
            receptor_pdbqt,
            center,
            box_size,
            center_by_ph,
            box_by_ph,
            center_source_by_ph,
            stages,
            params,
            control_stems,
            control_lookup,
        )
    except Exception as exc:
        docking_event = "fail"
        docking_error = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        docking_elapsed = time.time() - docking_start
        if run_id:
            try:
                update_manifest_for_docking_overall(
                    cfg,
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_tag=active_ph_label,
                    event=docking_event,
                    elapsed_sec=docking_elapsed,
                    error=docking_error,
                )
            except Exception:
                logger.warning(
                    "[run-manifest.docking-overall] failed to record end pdb=%s variant=%s ph=%s event=%s",
                    pdb_id,
                    variant_label,
                    active_ph_label if active_ph_label is not None else "base",
                    docking_event,
                    exc_info=True,
                )


__all__ = [
    "_map_reason_to_category",
    "RetryManager",
    "run_one_stage",
    "process_one_protein",
    "select_ligands_for_next",
    "_coerce_test_map",
    "_resolve_test_mode",
]
