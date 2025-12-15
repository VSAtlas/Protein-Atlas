# docking.py
# Docking orchestration helpers extracted from main.py

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from input_and_export_functions import (
    load_inputs,
    validate_config,
    define_docking_stages,
    write_score_summary_to_csv,
    extract_best_score,
    score_key,
    init_config_run_dir,
)
from docking_vina import emit_vina_config, emit_vina_config as _emit_vina_config_impl
from rdkit.Chem import rdMolAlign
from tqdm import tqdm

import automate_protein_prep as protein_prep
from apo_holo_mode import (
    _debug_normalize_mode_token,
    _record_apo_holo_decision,
    _record_apo_holo_usage,
    _variant_receptor_path,
    delete_variant_trees,
    file_sha1,
    resolve_apo_holo_mode,
)
from checkpoints import (
    checkpoint_invalidate_from,
    checkpoint_mark_done,
    checkpoint_should_skip,
)
from fallback_recenter import (
    BudgetGuard,
    GlobalCenterGuard,
    RecenterParams,
    fallback_recentering_if_empty,
    validate_first_valid_pose,
)
from docking_centering import CenterDecision, CenterSelector
from logging_topics import make_protein_logger
# collapse_sanitized_names lives in main.py and is not used here.
from path_router import Paths, make_paths
from run_manifest import (
    PocketDetectionEvent,
    emit_pocket_detection_event,
    update_manifest_for_docking_overall,
    update_manifest_for_docking_stage,
)
from pose_validation import (
    attempt_fallback_recenter,
    compute_self_rmsd,
    extract_surface_atoms,
    filter_and_rewrite_poses_by_rmsd,
    validate_pose_pdbqt,
    _pose_centroid_from_pdbqt,
)
from prep_ligands import prep_ligands_from_pdb
from protein_functions import detect_active_site
from run_vina import run_docking_task, validate_all_poses
from docking_subruns import ProteinDockingContext, SubrunSpec, subruns_for_test_mode, run_ligand_pipeline_subrun

import docking_controls
from docking_controls import (
    build_control_lookup,
    receptor_sanity_check,
    extract_ligands_to_nolig,
    _ph_values_from_context,
    _ph_ligand_mode,
    _resolve_ph_scope,
    select_center_via_control_redock,
    detect_pocket,
    _summarize_ions_file,
)
from docking_ligands import (
    select_ligands_for_next,
    _count_heavy_atoms_from_pdbqt,
    _dedup_index_roots,
    _lib_roots_for_pdb,
    prepare_and_filter_ligands,
    _coerce_test_map,
    _resolve_test_mode,
    _read_any_lig,
    _is_readable_ref,
    compute_rmsd,
    validate_ligand,
)
from docking_receptor import prepare_receptor
from docking_utils import (
    norm,
    _file_md5,
    _round_tuple,
    _fingerprint_stage,
    final_pose_validation_and_screenshots,
    _write_audit_json,
    _pose_path_for,
    early_recenter_decision,
)


docking_controls._is_readable_ref = _is_readable_ref
docking_controls.compute_rmsd = compute_rmsd

# Bridge functions for docking helpers
import docking as _docking

_docking.emit_vina_config = emit_vina_config
_docking.norm = norm
_docking._read_any_lig = _read_any_lig
_docking.validate_ligand = validate_ligand


def _phase0_setup_paths_and_logger(cfg: Dict, pdb_file: str) -> Tuple[Paths, str, logging.Logger]:
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = re.sub(r'(?i)_cleaned$', '', base_id)
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=os.path.basename(pdb_file))

    logger = make_protein_logger(str(paths.docked_pdb_root()), pdb_id, cfg)
    logger.info(f"[paths] base_id={base_id} -> pdb_id={pdb_id}")
    logger.info(f"Processing protein: {pdb_file} (id={pdb_id})")
    return paths, pdb_id, logger


def _phase1_variant_and_ion_context(cfg: Dict, paths: Paths, logger: logging.Logger) -> Tuple[str, Optional[str], str, dict, dict, bool, Path]:
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
    return variant_env, variant_token, variant_label, pdb_audit, clean_audit, legacy_mode, receptor_target


def _phase2_to4_receptor_and_center(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    variant_env: str,
    variant_token: Optional[str],
    variant_label: str,
    pdb_audit: dict,
    clean_audit: dict,
    legacy_mode: bool,
    receptor_target: Path,
    active_ph_label: Optional[str],
) -> Tuple[
    Optional[str],
    Optional[str],
    dict,
    dict,
    Optional[Tuple[float, float, float]],
    Optional[Tuple[float, float, float]],
    str,
    List[str],
    Dict[str, Path],
]:
    control_stems: List[str] = []
    control_lookup: Dict[str, Path] = {}
    center: Optional[Tuple[float, float, float]] = None
    box_size: Optional[Tuple[float, float, float]] = None
    center_source = "none"

    logger.info("[ph.debug] calling prepare_receptor; PH_ENSEMBLE=%s", cfg.get("PH_ENSEMBLE", False))
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    provenance = getattr(prepare_receptor, "last_provenance", None)
    if provenance is None:
        try:
            provenance = getattr(protein_prep, "get_clean_provenance", lambda: "unknown")()
        except Exception:
            provenance = "unknown"
    logger.info("[receptor.clean.provenance] created_by=%s", provenance)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("Skipping protein due to prep failure.")
        return cleaned_pdb, None, pdb_audit, clean_audit, center, box_size, center_source, control_stems, control_lookup

    cleaned_hist = "none"
    if cleaned_pdb:
        clean_summary = _summarize_ions_file(cleaned_pdb)
        cleaned_hist = str(clean_summary.get("hist", "none"))
        clean_error = clean_summary.get("error")
        if clean_error not in (None, "missing"):
            logger.warning(
                "[ions.clean.counts] pdb=%s variant=%s action=skip err=%s",
                paths.pdb_id,
                variant_label,
                clean_error,
            )
        else:
            logger.info(
                "[ions.clean.counts] pdb=%s variant=%s file=%s present_pdb=%s",
                paths.pdb_id,
                variant_label,
                cleaned_pdb,
                cleaned_hist,
            )
        clean_audit[variant_label] = {
            "hist": cleaned_hist,
            "counts": dict(clean_summary.get("counts", {})),
            "metals_present": bool(clean_summary.get("metals_present", False)),
            "salts_present": bool(clean_summary.get("salts_present", False)),
            "file": str(cleaned_pdb),
            "error": clean_error,
        }
    else:
        clean_audit[variant_label] = {
            "hist": "missing",
            "counts": {},
            "metals_present": False,
            "salts_present": False,
            "file": "",
            "error": "missing",
        }

    try:
        probe_map = protein_prep.get_ion_probe_map(paths.pdb_id)  # type: ignore[attr-defined]
    except Exception as exc:
        logger.warning("[ions.summary] pdb=%s variant=%s action=skip err=%s", paths.pdb_id, variant_label, exc)
    else:
        before_counts = probe_map.get("strip_nonstandard:before", {})
        final_counts = probe_map.get("receptor_write", {})
        metals_before = sum(before_counts.values())
        metals_kept = sum(final_counts.values())
        metals_stripped = max(0, metals_before - metals_kept)
        logger.info("[ions.summary] pdb=%s variant=%s metals_kept=%d metals_stripped=%d", paths.pdb_id, variant_label, metals_kept, metals_stripped)

    if cleaned_pdb:
        if variant_env == "HOLO":
            skip_reason = "variant"
        elif not variant_env:
            skip_reason = "legacy"
        else:
            skip_reason = "disabled"
        logger.info(
            "[ions.prep-early] pdb=%s variant=%s action=skip reason=%s file=%s",
            paths.pdb_id,
            variant_label,
            skip_reason,
            cleaned_pdb,
        )

    _lig_count, control_stems = extract_ligands_to_nolig(paths, logger)
    try:
        from prep_ligands import prep_ligands_from_pdb
        prep_ligands_from_pdb(
            ligand_output_dir=paths.ligand_output_dir,
            ligands_mol2_dir=paths.ligands_mol2_dir,
            prepped_ligands_dir=paths.prepped_ligands_dir,
        )
        logger.info("[Controls] Prepped extracted controls ahead of redock.")
    except Exception as e:
        logger.warning(f"[Controls] Prepping extracted controls failed: {e}")

    ctrl_pdbqts: list[Path] = []
    for root in {paths.prepped_ligands_dir, Path(cfg["OUTPUT_LIGANDS_DIR"])}:
        if root.exists():
            ctrl_pdbqts.extend(root.glob("*.pdbqt"))

    logger.info(f"[Controls] Prepped control PDBQTs found (union): {len(ctrl_pdbqts)}")
    for p in ctrl_pdbqts[:10]:
        logger.info(f"[Controls]   {p.name}")

    control_lookup = build_control_lookup(paths)

    try:
        sel_center, sel_box = select_center_via_control_redock(
            cfg,
            paths,
            receptor_pdbqt,
            logger,
            variant=variant_token,
            ph_token=active_ph_label,
            legacy=legacy_mode,
        )
    except Exception as _e:
        sel_center, sel_box = (None, None)
        logger.exception(
            "[control-centers] helper errored during ctrl_redock for pdb=%s variant=%s: %s",
            paths.pdb_id,
            variant_label,
            _e,
        )
    if sel_center is not None:
        center, box_size, center_source = sel_center, sel_box, "control"
        logger.info(f"[control-redock] Using control-derived center {center} with box {box_size}")
    else:
        c2, b2, src = detect_active_site(cleaned_pdb)
        if c2:
            box_size = tuple(min(28.0, float(s)) for s in b2)
            center = c2
            center_source = src or "activesite"
            logger.info("[active-site] Using center %s with box %s source=%s", center, box_size, center_source)
        else:
            logger.error("Active-site detection failed (no usable controls, active-site returned None).")
            return cleaned_pdb, receptor_pdbqt, pdb_audit, clean_audit, center, box_size, center_source, control_stems, control_lookup
    if center is None:
        return cleaned_pdb, receptor_pdbqt, pdb_audit, clean_audit, center, box_size, center_source, control_stems, control_lookup

    try:
        import automate_protein_prep as _auto_prep_mod
    except Exception as ions_err:
        logger.warning("[ions] pocket_refine_skip err=%s", ions_err)
    else:
        if cleaned_pdb and variant_env == "HOLO":
            logger.info(
                "[ions.pocket-pass] pdb=%s variant=%s action=refine_with_center file=%s",
                paths.pdb_id,
                variant_label,
                cleaned_pdb,
            )
            try:
                _auto_prep_mod._maybe_strip_ions(
                    Path(cleaned_pdb),
                    cfg=cfg,
                    variant=variant_token,
                    pocket_center=center,
                )
            except Exception as pocket_err:
                logger.warning(
                    "[ions] pocket_refine_skip err=%s",
                    pocket_err,
                )

    if center_source == "control":
        side = float(cfg.get("CONTROL_BOX_A", 24.0))
        box_size = (side, side, side)

    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
    box_size = tuple(min(box_cap, float(s)) for s in box_size)
    logger.info(f"Initial box clamped to {box_size} (cap={box_cap} A)")

    try:
        c_print = tuple(round(float(x), 3) for x in center)
        b_print = tuple(round(float(x), 1) for x in box_size)
        print(f"[CENTER] source={center_source} center={c_print} box={b_print}")
    except Exception:
        pass

    # --- manifest: record pocket detection metadata (best-effort) ---
    try:
        manifest_run_id = cfg.get("RUN_ID")
        logger.debug(
            "[run-manifest.pocket_detection.call] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%r box=%r",
            manifest_run_id,
            paths.pdb_id,
            variant_label,
            active_ph_label,
            center_source,
            center,
            box_size,
        )
        if manifest_run_id:
            emit_pocket_detection_event(
                cfg,
                PocketDetectionEvent(
                    run_id=str(manifest_run_id),
                    pdb_id=paths.pdb_id,
                    variant_label=variant_label,
                    ph_tag=active_ph_label,
                    method=center_source,
                    center=center,
                    box_size=box_size,
                ),
            )
        else:
            logger.debug(
                "[run-manifest.pocket_detection.skip] no RUN_ID for pdb=%s variant=%s ph=%s",
                paths.pdb_id,
                variant_label,
                active_ph_label,
            )
    except Exception:
        logger.warning(
            "[run-manifest.pocket_detection.error] pdb=%s variant=%s ph=%s",
            paths.pdb_id,
            variant_label,
            active_ph_label,
            exc_info=True,
        )

    metals_added = 0
    cofactors_added = 0
    regen = False
    try:
        logger.info("[holo.restore.call] invoking for pdb=%s", paths.pdb_id)
        metals_added, cofactors_added, regen = protein_prep._holo_restore_from_input_if_needed(
            pdb_id=paths.pdb_id,
            cleaned_pdb=cleaned_pdb,
            output_pdbqt=receptor_pdbqt,
            config=cfg,
            center=center,
            box_size=box_size,
        )

    except Exception as _restore_err:
        logger.warning("[holo.restore] action=skip reason=%s", _restore_err)

    if variant_env == "HOLO":
        logger.info(
            "[holo.restore.summary] pdb=%s variant=%s metals_added=%d cofactors_added=%d regen=%s",
            paths.pdb_id,
            variant_env,
            metals_added,
            cofactors_added,
            regen,
        )

        if regen:
            logger.warning(
                "[holo.restore.regen] pdb=%s variant=%s regen=True; rebuilding receptor PDBQT from %s -> %s",
                paths.pdb_id,
                variant_env,
                cleaned_pdb,
                receptor_target,
            )
            try:
                ok_after = protein_prep.run_prepare_receptor(
                    input_pdb=cleaned_pdb,
                    output_pdbqt=str(receptor_target),
                    cfg=cfg,
                )
                if not ok_after or not receptor_target.exists():
                    logger.warning(
                        "[holo.restore.regen] status=failed pdb=%s; keeping previous receptor PDBQT=%s",
                        paths.pdb_id,
                        receptor_pdbqt,
                    )
                else:
                    receptor_pdbqt = str(receptor_target)
                    logger.info(
                        "[holo.restore.regen] status=ok pdb=%s receptor_pdbqt=%s",
                        paths.pdb_id,
                        receptor_pdbqt,
                    )
            except Exception as regen_err:
                logger.warning(
                    "[holo.restore.regen] status=error pdb=%s err=%s; keeping previous receptor PDBQT=%s",
                    paths.pdb_id,
                    regen_err,
                    receptor_pdbqt,
                )
        else:
            logger.info(
                "[holo.restore.regen] pdb=%s variant=%s regen=False; skipping receptor PDBQT rebuild",
                paths.pdb_id,
                variant_env,
            )

        if receptor_pdbqt:
            try:
                protein_prep.run_metal_site_audit(
                    pdb_id=paths.pdb_id,
                    router_paths=paths,
                    input_pdb_path=str(paths.input_pdb_path),
                    receptor_pdb_path=cleaned_pdb,
                    receptor_pdbqt_path=receptor_pdbqt,
                    center=center,
                    variant_label=variant_label,
                    ph_label=active_ph_label,
                )
            except Exception as audit_err:
                logger.warning(
                    "[holo.metal_audit] action=skip pdb=%s reason=%s",
                    paths.pdb_id,
                    audit_err,
                )

    resolved_mode = (str(cfg.get("_RESOLVED_APO_HOLO_MODE")) or "").strip().lower() or "legacy"
    if variant_env == "HOLO" and resolved_mode == "apo_vs_holo":
        apo_clean = _variant_receptor_path(paths.pdb_id, "APO", cfg)
        holo_clean = cleaned_pdb or _variant_receptor_path(paths.pdb_id, "HOLO", cfg)
        apo_path = Path(apo_clean) if apo_clean else None
        holo_path = Path(holo_clean) if holo_clean else None
        apo_exists = apo_path.exists() if apo_path else False
        holo_exists = holo_path.exists() if holo_path else False

        if not apo_exists or not holo_exists:
            logger.warning(
                "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=continue reason=missing_paths apo=%s holo=%s",
                paths.pdb_id,
                apo_clean,
                holo_clean,
            )
            _record_apo_holo_decision(cfg, paths.pdb_id, "HOLO", "missing_paths")
        else:
            logger.info(
                "[apo-vs-holo] compare.preflight apo=%s exists=%s holo=%s exists=%s",
                norm(apo_path), ("T" if apo_exists else "F"),
                norm(holo_path), ("T" if holo_exists else "F"),
            )
            try:
                apo_sha = file_sha1(str(apo_path))
                holo_sha = file_sha1(str(holo_path))
            except Exception as hash_err:
                logger.warning(
                    "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=skip reason=sha_error err=%s",
                    paths.pdb_id,
                    hash_err,
                )
                _record_apo_holo_decision(cfg, paths.pdb_id, "HOLO", "sha_error")
            else:
                if apo_sha == holo_sha:
                    logger.info(
                        "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=skip reason=identical apo_sha=%s holo_sha=%s",
                        paths.pdb_id,
                        apo_sha,
                        holo_sha,
                    )
                    audit_root = cfg.get("_ION_AUDIT", {})
                    pdb_entry = audit_root.get(paths.pdb_id) or audit_root.get(paths.pdb_id)
                    warn_needed = False
                    if isinstance(pdb_entry, dict):
                        input_info = pdb_entry.get("input_counts", {})
                        clean_map = pdb_entry.get("clean_counts", {}) or {}
                        holo_info = clean_map.get("HOLO") or clean_map.get(variant_label) or {}
                        if input_info.get("metals_present") or input_info.get("salts_present"):
                            warn_needed = True
                        if holo_info.get("metals_present") or holo_info.get("salts_present"):
                            warn_needed = True
                    if warn_needed:
                        logger.warning(
                            "[apo-vs-holo] unexpected_identical_after_ion_policy pdb=%s apo_sha=%s holo_sha=%s",
                            paths.pdb_id,
                            apo_sha,
                            holo_sha,
                        )
                    if receptor_pdbqt:
                        _record_apo_holo_usage(cfg, paths.pdb_id, variant_token, None, receptor_pdbqt)
                    _record_apo_holo_decision(cfg, paths.pdb_id, "HOLO", "skipped_preflight")
                    try:
                        delete_variant_trees(paths.pdb_id, "HOLO", cfg)
                    except Exception as cleanup_err:
                        logger.warning(
                            "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=cleanup_warn err=%s",
                            paths.pdb_id,
                            cleanup_err,
                        )
                    return cleaned_pdb, None, pdb_audit, clean_audit, center, box_size, center_source, control_stems, control_lookup
                else:
                    logger.info(
                        "[apo-vs-holo] pdb_id=%s variant=HOLO stage=preflight action=continue reason=not_identical apo_sha=%s holo_sha=%s",
                        paths.pdb_id,
                        apo_sha,
                        holo_sha,
                    )
                    _record_apo_holo_decision(cfg, paths.pdb_id, "HOLO", "not_identical")

    return cleaned_pdb, receptor_pdbqt, pdb_audit, clean_audit, center, box_size, center_source, control_stems, control_lookup


def _phase5_ph_ensemble_global(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    variant_env: str,
    variant_token: Optional[str],
    legacy_mode: bool,
    cleaned_pdb: Optional[str],
    receptor_pdbqt: Optional[str],
    center: Optional[Tuple[float, float, float]],
    box_size: Optional[Tuple[float, float, float]],
) -> Tuple[Optional[str], Optional[str]]:
    active_ph_label = None
    if bool(cfg.get("PH_ENSEMBLE", False)):
        try:
            from context_ph import select_ph_values_for_protonation
            from ph_ensemble import build_ph_ensemble

            raw_vals = select_ph_values_for_protonation(str(paths.input_pdb_path))
            logger.info("[ph.ctx.raw] path=%s values=%s", str(paths.input_pdb_path),
                        ",".join(f"{v:.2f}" for v in (raw_vals or [])))

            ph_values = sorted({max(3.0, min(10.5, round(float(x), 1))) for x in (raw_vals or [])})
            if not ph_values:
                logger.warning("[ph.ctx.fallback] context list empty -> using [7.0]")
                ph_values = [7.0]

            logger.info("[ph.list] n=%d values=%s", len(ph_values),
                        ",".join(f"{v:.1f}" for v in ph_values))

            manifest_path = build_ph_ensemble(
                pdb_id=paths.pdb_id,
                cleaned_receptor_pdb=str(Path(cleaned_pdb)),
                out_dir=str(Path(cfg["OUTPUT_DIR"])),
                center=(0.0, 0.0, 0.0),
                radius=1_000_000.0,
                ph_values=ph_values,
                member_index_start=0,
                variant=variant_token,
                legacy=legacy_mode,
            )
            logger.info("[ph_ensemble.manifest] path=%s", manifest_path)

        except Exception as e:
            logger.warning("[ph_ensemble.skip] error=%s", e)
    return receptor_pdbqt, active_ph_label


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
        stages=stages,
        recenter_params=params,
        control_stems=control_stems,
        control_lookup=control_lookup,
    )

    test_mode = _resolve_test_mode(cfg)
    subruns: List[SubrunSpec] = subruns_for_test_mode(test_mode)

    for sub in subruns:
        logger.info(
            "[subrun] mode=%s run_mode=%r csv_prefix=%r stage_prefix=%r",
            test_mode,
            sub.run_mode,
            sub.csv_prefix,
            sub.stage_name_prefix,
        )
        run_ligand_pipeline_subrun(ctx, sub)


def process_one_protein(cfg: Dict, pdb_file: str, stages: List[Dict], params: RecenterParams) -> None:
    paths, pdb_id, logger = _phase0_setup_paths_and_logger(cfg, pdb_file)
    variant_env, variant_token, variant_label, pdb_audit, clean_audit, legacy_mode, receptor_target = _phase1_variant_and_ion_context(cfg, paths, logger)
    active_ph_label = None

    cleaned_pdb, receptor_pdbqt, pdb_audit, clean_audit, center, box_size, _center_source, control_stems, control_lookup = _phase2_to4_receptor_and_center(
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
    if not receptor_pdbqt or center is None or box_size is None:
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
    ph_label: Optional[str] = None,
) -> Tuple[
    Dict[str, float],
    List[str],
    List[float],
    Dict[str, str],
    Dict[str, Tuple[Optional[float], str]]
]:
    from sys import stdout as _stdout  # for tqdm
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

    def _validate_with_rmsd_gate(
        lig_path: str,
        lig_name: str,
        out_pdbqt_path: str,
        score_val: float
    ) -> Tuple[bool, Optional[str]]:
        """
        Controls: crystal redock RMSD is a hard gate.
        Non-controls: self-RMSD is logged upstream; do not gate here.
        """
        base = Path(lig_path).stem.split("_stage")[0]
        crystal_ref = control_lookup.get(base)

        if crystal_ref:
            best_pdb = _best_pose_pdb_from_pdbqt(out_pdbqt_path, obabel_path=cfg.get("OPENBABEL_PATH"))
            if not best_pdb:
                logger.warning(f"{lig_name} | unable to extract best pose PDB for redock RMSD.")
                return False, "no_best_pose_for_rmsd"

            # --- AUDIT: control redock (compute RMSD just for logging) ---
            rmsd_val = compute_rmsd(str(crystal_ref), best_pdb)

            ok = validate_ligand(
                ligand_name=lig_name,
                docked_path=best_pdb,
                crystal_path=str(crystal_ref),
                rmsd_thresh=float(cfg.get("CONTROL_RMSD_MAX_ANG", 2.0)),
                self_rmsd=None,
                logger=logger
            )

            logger.info("[control-redock] lig=%s rmsd=%.2f A score=%.2f",
                        lig_name, (rmsd_val if rmsd_val is not None else float('nan')), float(score_val))

            if ok:
                logger.info(f"{lig_name} | {stage['name']} score: {score_val:.2f} kcal/mol (redock-RMSD PASS)")
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
        nonlocal scores, validated_ligands, all_distances, raw_docked_ligands, invalids, processed

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

        if _handle_near_miss_retry(
            lig=lig,
            lig_name=lig_name,
            out_path=out_path,
            score=score,
            result=result,
            guard=guard,
            pbar=pbar,
        ):
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

    def _handle_near_miss_retry(
        lig: str,
        lig_name: str,
        out_path: str,
        score: float,
        result: Dict[str, Any],
        guard: BudgetGuard,
        pbar: "tqdm",
    ) -> bool:
        nonlocal processed
        if not bool(cfg.get("RETRY_NEAR_MISS", True)):
            return False

        try:
            if bool(cfg.get("LOG_SELF_RMSD", True)):
                try:
                    self_rmsd_val = compute_self_rmsd(out_path)
                except Exception as _e:
                    self_rmsd_val = None
                    logger.warning(f"self-RMSD failed for {lig_name}: {_e}")
                try:
                    logger.info(f"[self-rmsd] lig={lig_name} rmsd={self_rmsd_val}")
                except Exception:
                    pass

            reason = result.get("reason", "") or ""
            near_miss = (
                ("clash" in reason) or
                (result.get("distance_to_surface") or 0.0) < 6.5 or
                (result.get("distance_to_centroid") or 0.0) < 4.0
            )

            if not near_miss:
                return False

            stage_retry = dict(stage)
            stage_retry["name"] = f"{stage['name']}_retry"
            stage_retry["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
            if cfg.get("FAST_MODE"):
                stage_retry["exhaustiveness"] = 1

            try:
                base_seed = int(stage_retry.get("seed", 0)) if "seed" in stage_retry else 0
            except Exception:
                base_seed = 0
            stage_retry["seed"] = base_seed + 137
            stage_retry["num_modes"] = int(cfg.get("NEAR_MISS_NUM_MODES", stage_retry.get("num_modes", 4)))
            stage_retry["energy_range"] = float(cfg.get("NEAR_MISS_ENERGY_RANGE", stage_retry.get("energy_range", 4.0)))

            if bool(cfg.get("NEAR_MISS_RECENTER", True)):
                try:
                    cent = _pose_centroid_from_pdbqt(str(out_path))
                except Exception as _e:
                    logger.warning(f"[near-miss] failed to compute centroid for {lig_name}: {_e}")
                    cent = None

                if cent and isinstance(cent, (list, tuple)) and len(cent) == 3:
                    try:
                        center_nm = tuple(float(x) for x in cent)
                        logger.info(f"[near-miss] recentering on best pose centroid {center_nm}")
                        try:
                            max_box = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                        except Exception:
                            max_box = 28.0
                        box_nm = tuple(min(max_box, s + 1.0) for s in box_size)
                    except Exception:
                        center_nm = center
                        box_nm = box_size
                else:
                    center_nm = center
                    box_nm = box_size
            else:
                center_nm = center
                box_nm = box_size

            conf_path2, out_path2 = emit_vina_config(
                cfg,
                pdb_id,
                receptor_pdbqt,
                center_nm,
                box_nm,
                lig,
                stage_retry["name"],
                stage_retry,
                threads_per_vina,
                logger,
                variant=variant_token,
                ph_token=ph_label,
                legacy=legacy_mode,
            )
            logger.info("[cfg.emit] %s -> %s", os.path.basename(lig), conf_path2)
            try:
                Path(conf_path2).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
            except Exception:
                raise RuntimeError(
                    f"Refusing to launch Vina with config outside current RUN_DIR: {conf_path2}")
            logger.info(f"[vina.call] config={conf_path2}")

            if guard.expired():
                invalids[lig] = (float(score), "budget_exceeded")
                processed += 1
                if (processed % 25 == 0) or (processed == len(futures)):
                    pbar.set_postfix(ok=len(scores), inv=len(invalids))
                pbar.update(1)
                return True

            _, score2 = run_docking_task(cfg["VINA_EXE"], conf_path2, lig, out_path2)

            try:
                filter_and_rewrite_poses_by_rmsd(
                    out_path2,
                    rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                    max_models=int(cfg.get("RMSD_MAX_MODELS", 3))
                )
            except Exception:
                pass

            result2 = validate_first_valid_pose(
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

            if result2.get("valid", False) and score2 is not None:
                ok2, reason2 = _validate_with_rmsd_gate(lig, lig_name, out_path2, float(score2))
                if ok2:
                    scores[lig] = float(score2)
                    validated_ligands.append(lig)
                    raw_docked_ligands[lig] = norm(out_path2)
                    logger.info(f"{lig_name} | {stage_retry['name']} score: {score2:.2f} kcal/mol (rescued)")
                    processed += 1
                    if (processed % 25 == 0) or (processed == len(futures)):
                        pbar.set_postfix(ok=len(scores), inv=len(invalids))
                    pbar.update(1)
                    return True
                else:
                    invalids[lig] = (float(score2), reason2 or "rmsd_fail")
                    processed += 1
                    if (processed % 25 == 0) or (processed == len(futures)):
                        pbar.set_postfix(ok=len(scores), inv=len(invalids))
                    pbar.update(1)
                    return True
        except Exception as _e:
            logger.warning(f"Retry path failed for {lig_name}: {_e}")
        return False

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
        err_cat = _map_reason_to_category(result.get("reason", ""))
        attempt = 0
        retained_invalid = True

        while attempt < retry_mgr.max_retries:
            if guard.expired():
                invalids[lig] = (float(score), "budget_exceeded")
                break

            recipe = retry_mgr.apply(stage, err_cat, attempt)
            if not recipe:
                break

            stage_retry2 = dict(stage)
            stage_retry2["name"] = f"{stage['name']}_r{attempt + 1}"
            stage_retry2["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))

            if "exhaustiveness" in recipe:
                stage_retry2["exhaustiveness"] = recipe["exhaustiveness"]
            if "num_modes" in recipe:
                stage_retry2["num_modes"] = recipe["num_modes"]
            if cfg.get("FAST_MODE"):
                stage_retry2["exhaustiveness"] = 1
            if recipe.get("seed_jitter", False):
                try:
                    base_seed = int(stage_retry2.get("seed", 0)) if "seed" in stage_retry2 else 0
                except Exception:
                    base_seed = 0
                stage_retry2["seed"] = base_seed + (attempt + 1) * 137

            retry_center = center
            retry_box = box_size
            try:
                if recipe.get("recenter", False):
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
        return True

    # ---- scheduling & submission ----
    default_budget_seconds = float(
        cfg.get("MAX_RETRY_SECONDS_PER_LIGAND", cfg.get("BENCH_MAX_SECONDS", 300.0))
    )

    submit_queue = []
    guards_for_ligand: Dict[str, BudgetGuard] = {}

    for lig in ligands:
        guard = (budget_guards.get(lig) if budget_guards else None)
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
                    Path(conf_path).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
                except Exception:
                    raise RuntimeError(f"Refusing to launch Vina with config outside current RUN_DIR: {conf_path}")

                logger.info(f"[vina.call] config={conf_path}")

                lig_n, out_n = norm(lig), norm(out_path)
                raw_docked_ligands[lig_n] = out_n
                futures[pool.submit(run_docking_task, cfg["VINA_EXE"], conf_path, lig, out_path)] = (lig_n, out_n)

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
                    guard = guards_for_ligand.get(lig) or BudgetGuard(default_budget_seconds)

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


__all__ = [
    "_map_reason_to_category",
    "RetryManager",
    "run_one_stage",
    "process_one_protein",
    "select_ligands_for_next",
    "_coerce_test_map",
    "_resolve_test_mode",
]
