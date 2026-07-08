"""Top-level protein docking orchestration."""

from __future__ import annotations

import time
from typing import Dict, List, Optional

from cli.run_manifest_runtime import (
    update_manifest_for_docking_overall,
    update_manifest_for_prep_stage,
)
from docking import docking as _docking
from docking import docking_controls
from docking.docking_control_centering_phase import _phase5b_controls_and_control_redock
from docking.docking_control_redock import (
    _collect_control_pdbqts,
    run_control_docking_multi_engine,
)
from docking.docking_ligands import validate_ligand
from docking.docking_ph_ensemble_phase import _phase5_ph_ensemble_global
from docking.docking_receptor_phases import _phase2_to4_receptor_and_center
from docking.docking_utils import norm
from docking.docking_vina import emit_vina_config
from docking.fallback_recenter import RecenterParams
from docking.ligand_metrics import _is_readable_ref, _read_any_lig, compute_rmsd
from docking.docking_runtime_state import (
    _COMBO_PREP_BOX_KEY,
    _COMBO_PREP_CENTER_KEY,
    _COMBO_PREP_CONTROL_LOOKUP_KEY,
    _COMBO_PREP_CONTROL_STEMS_KEY,
    _COMBO_PREP_SOURCE_KEY,
    _decode_control_lookup,
    _decode_source_map,
    _decode_vec3_map,
    _encode_ph_key,
    _encode_source_map,
    _encode_vec3_map,
    _signal_distributed_prep_ready,
)
from docking.docking_runtime_phases import (
    _phase0_setup_paths_and_logger,
    _phase1_variant_and_ion_context,
    _phase6_to8_ligands_and_docking,
)

setattr(_docking, "emit_vina_config", emit_vina_config)
setattr(_docking, "norm", norm)
setattr(_docking, "_read_any_lig", _read_any_lig)
setattr(_docking, "validate_ligand", validate_ligand)
setattr(docking_controls, "_is_readable_ref", _is_readable_ref)
setattr(docking_controls, "compute_rmsd", compute_rmsd)


def _restore_control_reuse(cfg: Dict, control_stems, control_lookup):
    center_by_ph = _decode_vec3_map(cfg.get(_COMBO_PREP_CENTER_KEY))
    box_by_ph = _decode_vec3_map(cfg.get(_COMBO_PREP_BOX_KEY))
    center_source_by_ph = _decode_source_map(cfg.get(_COMBO_PREP_SOURCE_KEY))
    reused_control_stems = cfg.get(_COMBO_PREP_CONTROL_STEMS_KEY)
    reused_control_lookup = _decode_control_lookup(cfg.get(_COMBO_PREP_CONTROL_LOOKUP_KEY))
    if isinstance(reused_control_stems, list):
        control_stems = [str(x) for x in reused_control_stems if str(x).strip()]
    if reused_control_lookup:
        control_lookup = reused_control_lookup
    return center_by_ph, box_by_ph, center_source_by_ph, control_stems, control_lookup


def _ensure_control_redock_context(
    cfg: Dict,
    paths,
    logger,
    *,
    pdb_id: str,
    variant_token,
    variant_label: str,
    legacy_mode: bool,
    cleaned_pdb,
    receptor_pdbqt,
    active_ph_label,
    center,
    box_size,
    center_by_ph,
    box_by_ph,
    center_source_by_ph,
    control_stems,
    control_lookup,
):
    skip_control_redock = bool(cfg.get("_SKIP_CONTROL_REDOCK", False))
    if skip_control_redock and center_by_ph and box_by_ph:
        logger.info(
            "[control-redock] skip reason=distributed_prep_reuse pdb=%s variant=%s ph=%s centers=%d",
            pdb_id,
            variant_label,
            _encode_ph_key(active_ph_label),
            len(center_by_ph),
        )
    else:
        if skip_control_redock:
            logger.info(
                "[control-redock] reuse_requested_but_unavailable pdb=%s variant=%s ph=%s; running full control-redock",
                pdb_id,
                variant_label,
                _encode_ph_key(active_ph_label),
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

    if center_by_ph and (None not in center_by_ph) and center is None and box_size is None:
        center = center_by_ph.get(active_ph_label) or center_by_ph.get(None)
        box_size = box_by_ph.get(active_ph_label) if box_by_ph else None
    if not center_by_ph and center is not None and box_size is not None:
        center_by_ph = {active_ph_label: center}
        box_by_ph = {active_ph_label: box_size}
    return (
        center,
        box_size,
        center_by_ph,
        box_by_ph,
        center_source_by_ph,
        control_stems,
        control_lookup,
    )


def _persist_control_context(
    cfg: Dict,
    center_by_ph,
    box_by_ph,
    center_source_by_ph,
    control_stems,
    control_lookup,
) -> None:
    if center_by_ph:
        cfg[_COMBO_PREP_CENTER_KEY] = _encode_vec3_map(center_by_ph)
    if box_by_ph:
        cfg[_COMBO_PREP_BOX_KEY] = _encode_vec3_map(box_by_ph)
    if center_source_by_ph:
        cfg[_COMBO_PREP_SOURCE_KEY] = _encode_source_map(center_source_by_ph)
    if control_stems:
        cfg[_COMBO_PREP_CONTROL_STEMS_KEY] = [str(x) for x in control_stems]
    if control_lookup:
        cfg[_COMBO_PREP_CONTROL_LOOKUP_KEY] = {str(k): str(v) for k, v in control_lookup.items()}


def _active_ph_from_reuse_payload(cfg: Dict, fallback_ph) -> Optional[str]:
    raw_ph = str(fallback_ph or "").strip()
    if raw_ph:
        return raw_ph
    center_by_ph = _decode_vec3_map(cfg.get(_COMBO_PREP_CENTER_KEY))
    if len(center_by_ph) == 1:
        only_key = next(iter(center_by_ph.keys()))
        return only_key
    return None


def _try_run_distributed_prep_reuse_pipeline(
    cfg: Dict,
    paths,
    logger,
    *,
    pdb_id: str,
    variant_token,
    variant_label: str,
    legacy_mode: bool,
    stages: List[Dict],
    params: RecenterParams,
    run_id: str,
) -> bool:
    if not bool(cfg.get("_DISTRIBUTED_PREP_REUSE_ONLY", False)):
        return False

    (
        center_by_ph,
        box_by_ph,
        center_source_by_ph,
        control_stems,
        control_lookup,
    ) = _restore_control_reuse(cfg, [], {})
    if not center_by_ph or not box_by_ph:
        logger.info(
            "[prep.reuse] action=skip reason=missing_center_box pdb=%s variant=%s",
            pdb_id,
            variant_label,
        )
        return False

    active_ph_label = _active_ph_from_reuse_payload(
        cfg,
        cfg.get("_PH_TAG_OVERRIDE"),
    )
    center = center_by_ph.get(active_ph_label) or center_by_ph.get(None)
    box_size = box_by_ph.get(active_ph_label) or box_by_ph.get(None)
    if center is None or box_size is None:
        logger.info(
            "[prep.reuse] action=skip reason=missing_ph_center pdb=%s variant=%s ph=%s centers=%s",
            pdb_id,
            variant_label,
            active_ph_label if active_ph_label is not None else "base",
            ",".join(str(k if k is not None else "base") for k in center_by_ph.keys()),
        )
        return False

    cleaned_pdb = paths.receptor_cleaned_pdb(variant_token)
    receptor_pdbqt = paths.receptor_pdbqt(variant_token, active_ph_label)
    if not receptor_pdbqt.exists():
        logger.info(
            "[prep.reuse] action=skip reason=receptor_missing pdb=%s variant=%s ph=%s receptor=%s",
            pdb_id,
            variant_label,
            active_ph_label if active_ph_label is not None else "base",
            receptor_pdbqt,
        )
        return False

    if run_id:
        try:
            update_manifest_for_prep_stage(
                cfg,
                run_id,
                pdb_id,
                variant_label,
                "completed",
                ph_tag=active_ph_label,
                details={
                    "reuse": "distributed_combo_prep",
                    "cleaned_pdb": str(cleaned_pdb),
                    "receptor_pdbqt": str(receptor_pdbqt),
                },
            )
        except Exception:
            logger.warning(
                "[prep.reuse] manifest update failed pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                active_ph_label if active_ph_label is not None else "base",
                exc_info=True,
            )

    logger.info(
        "[prep.reuse] action=direct_docking pdb=%s variant=%s ph=%s receptor=%s centers=%d",
        pdb_id,
        variant_label,
        active_ph_label if active_ph_label is not None else "base",
        receptor_pdbqt,
        len(center_by_ph),
    )
    _run_docking_pipeline_with_manifest(
        cfg,
        logger,
        run_id=run_id,
        pdb_id=pdb_id,
        variant_label=variant_label,
        active_ph_label=active_ph_label,
        runner=lambda: _phase6_to8_ligands_and_docking(
            cfg,
            paths,
            logger,
            pdb_id,
            str(variant_token or ""),
            variant_token,
            variant_label,
            legacy_mode,
            str(cleaned_pdb) if cleaned_pdb.exists() else None,
            str(receptor_pdbqt),
            center,
            box_size,
            center_by_ph,
            box_by_ph,
            center_source_by_ph,
            stages,
            params,
            control_stems,
            control_lookup,
        ),
    )
    return True


def _maybe_run_control_multi(
    cfg: Dict,
    paths,
    logger,
    *,
    pdb_id: str,
    variant_token,
    variant_label: str,
    legacy_mode: bool,
    active_ph_label,
    center_by_ph,
    box_by_ph,
    control_stems,
    control_lookup,
) -> None:
    if not center_by_ph:
        return
    if bool(cfg.get("_SKIP_CONTROL_DOCKING", False)):
        logger.info(
            "[control-multi] skip reason=distributed_prep_reuse pdb=%s variant=%s ph=%s",
            pdb_id,
            variant_label,
            _encode_ph_key(active_ph_label),
        )
        return

    control_ligands = _collect_control_pdbqts(paths, control_stems, cfg)
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


def _run_docking_pipeline_with_manifest(
    cfg: Dict,
    logger,
    *,
    run_id: str,
    pdb_id: str,
    variant_label: str,
    active_ph_label,
    runner,
) -> None:
    event = "end"
    error: Optional[str] = None
    started = time.time()
    record_overall_manifest = bool(run_id) and not (
        bool(cfg.get("PH_ENSEMBLE", False)) and active_ph_label is None
    )
    if run_id and not record_overall_manifest:
        logger.debug(
            "[run-manifest.docking-overall] skip reason=ph_ensemble_without_single_ph pdb=%s variant=%s",
            pdb_id,
            variant_label,
        )
    if record_overall_manifest:
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
        runner()
    except Exception as exc:
        event = "fail"
        error = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        if record_overall_manifest:
            try:
                update_manifest_for_docking_overall(
                    cfg,
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_tag=active_ph_label,
                    event=event,
                    elapsed_sec=time.time() - started,
                    error=error,
                )
            except Exception:
                logger.warning(
                    "[run-manifest.docking-overall] failed to record end pdb=%s variant=%s ph=%s event=%s",
                    pdb_id,
                    variant_label,
                    active_ph_label if active_ph_label is not None else "base",
                    event,
                    exc_info=True,
                )


def process_one_protein(
    cfg: Dict, pdb_file: str, stages: List[Dict], params: RecenterParams
) -> None:
    paths, pdb_id, logger = _phase0_setup_paths_and_logger(cfg, pdb_file)
    run_id = str(cfg.get("RUN_ID", "") or "")
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
    prep_started = time.time()
    if run_id:
        update_manifest_for_prep_stage(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            "running",
            ph_tag=active_ph_label,
            details={"input_pdb": str(paths.input_pdb_path)},
        )

    if _try_run_distributed_prep_reuse_pipeline(
        cfg,
        paths,
        logger,
        pdb_id=pdb_id,
        variant_token=variant_token,
        variant_label=variant_label,
        legacy_mode=legacy_mode,
        stages=stages,
        params=params,
        run_id=run_id,
    ):
        return

    try:
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
    except Exception as exc:
        if run_id:
            update_manifest_for_prep_stage(
                cfg,
                run_id,
                pdb_id,
                variant_label,
                "failed",
                ph_tag=active_ph_label,
                elapsed_sec=time.time() - prep_started,
                error=f"{type(exc).__name__}: {exc}",
            )
        raise
    if not receptor_pdbqt:
        reason = "receptor_pdbqt_missing"
        if run_id:
            update_manifest_for_prep_stage(
                cfg,
                run_id,
                pdb_id,
                variant_label,
                "failed",
                ph_tag=active_ph_label,
                elapsed_sec=time.time() - prep_started,
                details={"cleaned_pdb": str(cleaned_pdb or "")},
                error=reason,
            )
        raise RuntimeError(reason)
    if run_id:
        update_manifest_for_prep_stage(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            "completed",
            ph_tag=active_ph_label,
            elapsed_sec=time.time() - prep_started,
            details={
                "cleaned_pdb": str(cleaned_pdb or ""),
                "receptor_pdbqt": str(receptor_pdbqt or ""),
                "center_source": str(center_source or "none"),
            },
        )

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
        center_by_ph,
        box_by_ph,
        center_source_by_ph,
        control_stems,
        control_lookup,
    ) = _restore_control_reuse(cfg, control_stems, control_lookup)
    (
        center,
        box_size,
        center_by_ph,
        box_by_ph,
        center_source_by_ph,
        control_stems,
        control_lookup,
    ) = _ensure_control_redock_context(
        cfg,
        paths,
        logger,
        pdb_id=pdb_id,
        variant_token=variant_token,
        variant_label=variant_label,
        legacy_mode=legacy_mode,
        cleaned_pdb=cleaned_pdb,
        receptor_pdbqt=receptor_pdbqt,
        active_ph_label=active_ph_label,
        center=center,
        box_size=box_size,
        center_by_ph=center_by_ph,
        box_by_ph=box_by_ph,
        center_source_by_ph=center_source_by_ph,
        control_stems=control_stems,
        control_lookup=control_lookup,
    )
    _persist_control_context(
        cfg,
        center_by_ph,
        box_by_ph,
        center_source_by_ph,
        control_stems,
        control_lookup,
    )
    _signal_distributed_prep_ready(
        cfg,
        logger,
        pdb_id=pdb_id,
        variant_label=variant_label,
        active_ph_label=active_ph_label,
    )
    if (center is None or box_size is None) and not center_by_ph:
        return

    _maybe_run_control_multi(
        cfg,
        paths,
        logger,
        pdb_id=pdb_id,
        variant_token=variant_token,
        variant_label=variant_label,
        legacy_mode=legacy_mode,
        active_ph_label=active_ph_label,
        center_by_ph=center_by_ph,
        box_by_ph=box_by_ph,
        control_stems=control_stems,
        control_lookup=control_lookup,
    )

    _run_docking_pipeline_with_manifest(
        cfg,
        logger,
        run_id=run_id,
        pdb_id=pdb_id,
        variant_label=variant_label,
        active_ph_label=active_ph_label,
        runner=lambda: _phase6_to8_ligands_and_docking(
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
        ),
    )
