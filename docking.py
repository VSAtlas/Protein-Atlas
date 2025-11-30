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
    emit_vina_config,
    emit_vina_config as _emit_vina_config_impl,
    record_score,
    score_key,
    _to_bool,
    init_config_run_dir,
)
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
from path_router import (
    RouterPaths,
    Paths,
    docked_dir,
    make_paths,
    receptor_file,
)
from single_ligand_index import (
    _ensure_single_ligand_index,
    _resolve_single_ligand,
)
from ph_ensemble_docking import (
    enumerate_ligands_for_ph_context,
    init_ph_tags_and_manifest,
    prewarm_ph_ligand_microstates,
)
from run_manifest import (
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
    update_manifest_for_protein_success,
)
from pose_validation import (
    attempt_fallback_recenter,
    compute_self_rmsd,
    extract_surface_atoms,
    filter_and_rewrite_poses_by_rmsd,
    validate_pose_pdbqt,
    _pose_centroid_from_pdbqt,
)
from prep_ligands import (
    enumerate_ligands_for_docking,
    prep_ligands_from_pdb,
)
from protein_functions import detect_active_site
from record_data import record_le, write_scores_csv
from run_vina import run_docking_task, validate_all_poses

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
    _coerce_test_map,
    _resolve_test_mode,
    _dedup_index_roots,
    _lib_roots_for_pdb,
    prepare_and_filter_ligands,
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
        logger.debug(f"[control-centers] helper errored: {_e}")
    if sel_center is not None:
        center, box_size, center_source = sel_center, sel_box, "control"
        logger.info(f"[control-redock] Using control-derived center {center} with box {box_size}")
    else:
        c2, b2 = detect_active_site(cleaned_pdb)
        if c2:
            box_size = tuple(min(28.0, float(s)) for s in b2)
            center = c2
            center_source = "p2rank"
            logger.info(f"[P2Rank] Using P2Rank center {center} with box {box_size}")
        else:
            logger.error("Active-site detection failed (no usable controls, P2Rank returned None).")
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


def _run_ligand_pipeline_subrun(
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
    *,
    run_mode: Optional[str],
    csv_prefix: str,
    stage_name_prefix: str,
) -> None:
    """
    Run the existing ligands + multi-stage docking pipeline once,
    but parameterized by:
      - run_mode: None | "dud" | "fda"
      - csv_prefix: "" or "dud_"
      - stage_name_prefix: "" or "dud_"
    Variant and pH behavior must remain unchanged: only the final stage
    component gets the prefix.
    """
    cfg.setdefault("_EFFECTIVE_SINGLE_LIGAND", "")
    single_ligand_hit: Optional[Path] = None
    cfg.pop("_SINGLE_RESOLVED_PATH", None)
    if cfg["_EFFECTIVE_SINGLE_LIGAND"]:
        _ensure_single_ligand_index(cfg, paths, logger)
        cfg.setdefault("paths", {})
        cfg["paths"]["prepped_ligands_dir"] = str(paths.prepped_ligands_dir)

        hit = _resolve_single_ligand(cfg["_EFFECTIVE_SINGLE_LIGAND"], pdb_id, cfg, logger)
        if hit:
            cfg["_SINGLE_RESOLVED_PATH"] = str(hit)
            single_ligand_hit = hit
        else:
            selector_token = cfg["_EFFECTIVE_SINGLE_LIGAND"]
            suggestions: list[str] = []
            try:
                import difflib

                fda_map = _load_fda_name_map(cfg, logger)
                suggestions = difflib.get_close_matches(
                    selector_token,
                    list(fda_map.keys()),
                    n=5,
                    cutoff=0.7,
                )
            except Exception:
                suggestions = []
            if suggestions:
                logger.error("[single.miss.suggest] did_you_mean=%s", ", ".join(suggestions))

            allow_flag = os.environ.get("ALLOW_FDA_FALLBACK")
            if allow_flag is None:
                allow_flag = cfg.get("ALLOW_FDA_FALLBACK", False)
            if not _to_bool(allow_flag):
                logger.error(
                    "[single.block] selector '%s' not found in fda_library via FDA_MAPPING_CSV; aborting instead of fallback.",
                    selector_token,
                )
                raise SystemExit(2)
            logger.warning(
                "[single.block] selector '%s' not found; ALLOW_FDA_FALLBACK enabled, continuing with fallback flow.",
                selector_token,
            )
            cfg["_EFFECTIVE_SINGLE_LIGAND"] = ""

    if cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
        if not single_ligand_hit:
            return
        ligands = [str(single_ligand_hit)]
        ha = _count_heavy_atoms_from_pdbqt(single_ligand_hit)
        heavy_atom_counts = {str(single_ligand_hit): ha}
        pains_flags = {}
        logger.info(f"[single] Active ? docking only: {single_ligand_hit.name} (heavy={ha})")

        if cfg.get("PH_LIGAND_MODE", "").lower() == "context_window" and cfg.get("PH_ENSEMBLE_IN_PREP"):
            try:
                from prep_ligands import enumerate_ligands_for_docking

                ph_values = [6.0, 8.0]
                if "_PH_CONTEXT_VALUES" in cfg:
                    ph_values = cfg["_PH_CONTEXT_VALUES"]
                ligand_window = sorted(
                    {round(p, 1) for ph in ph_values for p in (float(ph) - 1.0, float(ph), float(ph) + 1.0)}
                )
                logger.info(f"[single.ph_ligand] Using ligand window {ligand_window}")

                ph_override_single = "dud" if run_mode == "dud" else None
                _, noncontrol_roots_single = _lib_roots_for_pdb(
                    cfg,
                    paths.pdb_id.upper(),
                    paths,
                    logger,
                    test_mode_override=ph_override_single,
                )
                ph_root_path = noncontrol_roots_single[0] if noncontrol_roots_single else None
                ph_root_cfg = str(ph_root_path) if ph_root_path else ""

                logger.info(
                    "[single.ph_ligand.bridge] ph_root_cfg=%s ph_root_path=%s exists=%s",
                    ph_root_cfg,
                    str(ph_root_path) if ph_root_path is not None else "",
                    ph_root_path.exists() if ph_root_path is not None else False,
                )
                if ph_root_path is not None and ph_root_path.exists():
                    enumerate_ligands_for_docking(
                        requested_ph_values=ligand_window,
                        root_dir=ph_root_path,
                        microstate_dedup=True,
                        force=False,
                    )
                else:
                    logger.info(
                        "[single.ph_ligand.bridge.skip] no valid ph_ligand_root; "
                        "skipping microstate priming for single-ligand mode"
                    )
            except Exception as e:
                logger.warning(f"[single.ph_ligand.skip] Could not run PH-ligand window for single mode: {e}")
    else:
        ligands, heavy_atom_counts, pains_flags = prepare_and_filter_ligands(
            cfg,
            paths,
            logger,
            run_mode=run_mode,
        )

    if not cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
        ctrl_stems_lower = {s.lower() for s in control_stems}

        prepped_control_pdbqts = []
        scan_roots = [paths.prepped_ligands_dir]
        if cfg.get("OUTPUT_LIGANDS_DIR"):
            try:
                out_root = Path(cfg["OUTPUT_LIGANDS_DIR"])
                if out_root.exists():
                    scan_roots.append(out_root)
            except Exception:
                pass

        for root in scan_roots:
            if root and root.exists():
                for p in root.glob("*.pdbqt"):
                    stem0 = p.stem.split("_stage")[0].lower()
                    if stem0 in ctrl_stems_lower:
                        prepped_control_pdbqts.append(p)

        lig_set = {norm(x) for x in ligands}
        missing_controls = [p for p in prepped_control_pdbqts if norm(p) not in lig_set]

        if missing_controls:
            logger.info(f"[Controls] Adding {len(missing_controls)} prepared control(s) to Stage1.")
            ligands = [str(p) for p in missing_controls] + ligands
            for p in missing_controls:
                try:
                    heavy_atom_counts.setdefault(str(p), _count_heavy_atoms_from_pdbqt(p))
                except Exception:
                    heavy_atom_counts.setdefault(str(p), 0)

    def _norm_dedupe(seq):
        seen = set()
        out = []
        for p in seq:
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                out.append(pn)
        return out

    ligands = _norm_dedupe(ligands)

    base_ligands = ligands[:]
    base_heavy_atoms = dict(heavy_atom_counts)
    base_pains_flags = dict(pains_flags)
    base_center = tuple(center)
    base_box = tuple(box_size)

    ph_log = logging.getLogger("ph_ensemble")
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))
    plan_only = os.environ.get("A2_PLAN_ONLY") == "1"

    manifest_run_id = cfg.get("RUN_ID")
    library_for_manifest = None
    try:
        mode_for_manifest = _resolve_test_mode(cfg)
        if mode_for_manifest != "off":
            lib_map = cfg.get("_TEST_LIBRARY_CANONICAL", {}) or {}
            library_for_manifest = lib_map.get(paths.pdb_id.upper())
        if not library_for_manifest:
            library_for_manifest = cfg.get("LIBRARY_SUBDIR_DEFAULT")
    except Exception:
        library_for_manifest = cfg.get("LIBRARY_SUBDIR_DEFAULT")

    ph_tags = init_ph_tags_and_manifest(cfg, paths.pdb_id, variant_token, legacy_mode)
    if ph_enabled and not ph_tags:
        logger.info(
            "[subrun.ph] run_mode=%s ph_enabled=True but no ph_tags; skipping PH run",
            run_mode or "None",
        )
        return

    ph_test_mode_override = "dud" if run_mode == "dud" else None
    _ctrl_roots_ph, noncontrol_roots_ph = _lib_roots_for_pdb(
        cfg,
        paths.pdb_id.upper(),
        paths,
        logger,
        test_mode_override=ph_test_mode_override,
    )
    ph_ligand_root = noncontrol_roots_ph[0] if noncontrol_roots_ph else None

    logger.info(
        "[subrun.ph] run_mode=%s ph_tags=%s ph_ligand_root=%s override=%s",
        run_mode or "None",
        ",".join(ph_tags) if ph_tags else "(none)",
        str(ph_ligand_root) if ph_ligand_root else "(none)",
        ph_test_mode_override or "(none)",
    )
    prewarm_ph_ligand_microstates(cfg, ph_tags, ph_ligand_root)

    stages_for_run = [
        {**stage, "name": f"{stage_name_prefix}{stage['name']}"}
        for stage in stages
    ]

    for ph_label in ph_tags:
        ph_start_ts = time.time()
        try:
            update_manifest_for_protein_start(
                cfg,
                manifest_run_id or "",
                paths.pdb_id,
                variant_label,
                library_for_manifest,
                ph_tag=ph_label,
            )
        except Exception:
            ph_log.warning(
                "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=start",
                paths.pdb_id,
                variant_label,
                ph_label if ph_label else "base",
                exc_info=True,
            )

        rec_path = receptor_file(paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode)
        out_root = docked_dir(paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode)
        ph_print = ph_label or "(none)"
        rec_exists = rec_path.exists()
        logger.info(
            "[router] pdb=%s variant=%s ph=%s receptor_file=%s docked_dir=%s exists=%s",
            paths.pdb_id,
            variant_label,
            ph_print,
            str(rec_path),
            str(out_root),
            rec_exists,
        )

        if plan_only:
            print(
                f"pdb={paths.pdb_id} variant={variant_label} ph={ph_print} "
                f"receptor_file={rec_path} docked_dir={out_root} exists={rec_exists}"
            )
            continue

        _record_apo_holo_usage(cfg, paths.pdb_id, variant_token, ph_label, rec_path)

        if not rec_exists:
            ph_log.warning(
                "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s receptor_missing=%s",
                paths.pdb_id,
                variant_label,
                ph_print,
                str(rec_path),
            )
            try:
                update_manifest_for_protein_failure(
                    cfg,
                    manifest_run_id or "",
                    paths.pdb_id,
                    variant_label,
                    rec_path,
                    ph_tag=ph_label,
                )
            except Exception:
                ph_log.warning(
                    "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=receptor-missing",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    exc_info=True,
                )
            continue

        try:
    
            ligands = base_ligands[:]
            heavy_atom_counts = dict(base_heavy_atoms)
            pains_flags = dict(base_pains_flags)
            center = tuple(base_center)
            box_size = tuple(base_box)
            receptor_pdbqt = str(rec_path)
    
            enumerated = enumerate_ligands_for_ph_context(
                cfg=cfg,
                pdb_id=paths.pdb_id,
                ph_label=ph_label,
                ph_ligand_root=ph_ligand_root,
            )
    
            if enumerated:
                ligands = [str(p) for p in enumerated]
                heavy_atom_counts = {
                    str(p): _count_heavy_atoms_from_pdbqt(p) for p in enumerated
                }
                pains_flags = {
                    k: base_pains_flags.get(
                        k,
                        base_pains_flags.get(Path(k).stem, False),
                    )
                    for k in ligands
                }
            else:
                ligands = base_ligands[:]
                heavy_atom_counts = dict(base_heavy_atoms)
                pains_flags = dict(base_pains_flags)
    
            ctrl_stems_lower = {s.lower() for s in control_stems}
            ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
            min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))
    
            def _is_control_path(p: str) -> bool:
                stem = Path(p).stem.split("_stage")[0]
                if stem.upper() in ctrl_blacklist:
                    return False
                if stem.lower() not in ctrl_stems_lower:
                    return False
                ha = heavy_atom_counts.get(p)
                return (ha is None) or (ha >= min_ha)
    
            ctrls = [p for p in ligands if _is_control_path(p)]
            non_ctrls = [p for p in ligands if not _is_control_path(p)]
            if ctrls:
                ligands = ctrls + non_ctrls
                logger.info(
                    f"[Controls] Front-loading {len(ctrls)} controls. "
                    f"First wave: {[Path(x).name for x in ligands[:int(cfg.get('MAX_PARALLEL_JOBS', 1))]]}"
                )
    
            present_ctrls = [
                Path(l).stem.split("_stage")[0].lower()
                for l in ligands
                if Path(l).stem.split("_stage")[0].lower() in ctrl_stems_lower
            ]
    
            if not present_ctrls:
                logger.warning(
                    "[Controls] No control ligands present in Stage1 ligand list -- "
                    "self-RMSD/locking will not be possible. (Check prep errors above.)"
                )
            if not ligands:
                logger.warning("No valid ligands after filtering; skipping protein.")
                continue
    
            selector = CenterSelector(cfg, logger, control_stems, heavy_atom_counts, center)
            guard = GlobalCenterGuard(
                max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2))
            )
    
            stage1_original = ligands[:]
    
            control_stems_lower = {s.lower() for s in control_stems}
            forced_extracted_for_stage3 = {
                lig for lig in stage1_original
                if Path(lig).stem.split("_stage")[0].lower() in control_stems_lower
            }
            logger.info(
                f"[Force-carry] Extracted ligands earmarked for Stage3: {len(forced_extracted_for_stage3)}"
            )
    
            score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
            validated_ligands_last: List[str] = []
            recenter_attempts = 0
            docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()
    
            retry_mgr = RetryManager()
    
            i = 0
            while i < len(stages_for_run):
                guard.reset_stage()
                stage = stages_for_run[i]
    
                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                    if checkpoint_should_skip(
                        cfg,
                        paths.pdb_id,
                        stage["name"],
                        fp,
                        ph_label=ph_label,
                        variant=variant_env or None,
                    ):
                        logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                        i += 1
                        continue
    
                if not ligands:
                    logger.warning(f"No ligands to dock at {stage['name']}; stopping for this protein.")
                    break
    
                logger.info(f"Starting {stage['name']} with {len(ligands)} ligands...")
                if ph_label:
                    stage_dir = paths.docked_stage_dir(variant_env or None, stage["name"], ph_label)
                    ph_log.info(
                        "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s stage=%s receptor=%s out=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label,
                        stage["name"],
                        receptor_pdbqt,
                        str(stage_dir),
                    )
    
                if i == 0 and ctrls and non_ctrls:
                    logger.info(
                        f"Stage1 two-wave: {len(ctrls)} controls first, then {len(non_ctrls)} others."
                    )
    
                    s1, v1, d1, rd1, inv1 = run_one_stage(
                        cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                        ctrls, logger, retry_mgr, control_lookup, ph_label=ph_label
                    )
    
                    try:
                        dec = selector.consider_switch(stage['name'], s1, v1, rd1, receptor_pdbqt, center, guard)
                        if dec.promoted and dec.new_center is not None:
                            old = center
                            center = dec.new_center
                            guard.mark_switch()
                            logger.info(
                                f"[CENTER] Switched before library run: {old} -> {center} ({dec.reason}) [global switch]"
                            )
                    except Exception as e:
                        logger.warning(f"CenterSelector (controls-only) failed gracefully: {e}")
    
                    lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", -6.0))
                    lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                    lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))
                    qualified_controls = []
    
                    for lig in v1:
                        stem = Path(lig).stem.split("_stage")[0].lower()
                        ha = heavy_atom_counts.get(lig)
                        if stem in control_stems_lower and (ha is None or ha >= min_ha):
                            sc = s1.get(lig)
                            if sc is not None and np.isfinite(sc) and sc <= lock_score_max:
                                pose_path = rd1.get(lig)
                                c = CenterSelector._pdbqt_centroid(pose_path) if pose_path else None
                                if c is not None and np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                                    qualified_controls.append(lig)
    
                    if len(qualified_controls) >= lock_min_hits and not guard.locked:
                        guard.lock()
                        logger.info(
                            "[CONTROL-LOCK] Early lock from controls-only wave "
                            f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} A); "
                            "future center switches disabled."
                        )
    
                    s2, v2, d2, rd2, inv2 = run_one_stage(
                        cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                        non_ctrls, logger, retry_mgr, control_lookup, ph_label=ph_label
                    )
    
                    scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
                    raw_docked = {**rd1, **rd2}
                    invalids = {**inv1, **inv2}
                else:
                    scores, validated, distances, raw_docked, invalids = run_one_stage(
                        cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                        ligands, logger, retry_mgr, control_lookup, ph_label=ph_label
                    )
    
                validated_ligands_last = validated
    
                def _is_control(lig: str) -> bool:
                    stem = Path(lig).stem.split("_stage")[0].lower()
                    if stem.upper() in {s.strip().upper() for s in cfg.get("CONTROL_BLACKLIST", "").split(",") if s.strip()}:
                        return False
                    ha = heavy_atom_counts.get(lig)
                    if ha is not None and ha < int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10)):
                        return False
                    return stem in {s.lower() for s in control_stems}
    
                control_anchor_hit = any(_is_control(lig) for lig in validated)
    
                lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", float("inf")))
                lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))
    
                qualified_controls = []
                for lig in validated:
                    if not _is_control(lig):
                        continue
                    sc = scores.get(lig)
                    if sc is None or not np.isfinite(sc):
                        continue
                    if sc > lock_score_max:
                        continue
                    pose_path = raw_docked.get(lig)
                    if not pose_path:
                        continue
                    c = CenterSelector._pdbqt_centroid(pose_path)
                    if c is None:
                        continue
                    if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                        qualified_controls.append(lig)
    
                if len(qualified_controls) >= lock_min_hits and not guard.locked:
                    guard.lock()
                    logger.info(
                        "[CONTROL-LOCK] Control(s) validated with strong confidence "
                        f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Ang); "
                        "center is now anchored; future center switches are disabled."
                    )
    
                try:
                    processed = {norm(x) for x in ligands}
                    valid_set = {norm(x) for x in scores.keys()}
                    invalid_set = {norm(x) for x in invalids.keys()}
                    both = valid_set & invalid_set
                    missing = processed - (valid_set | invalid_set)
                    if both or missing:
                        logger.error(
                            f"Invariant violation at {stage['name']}: both={len(both)}, missing={len(missing)}"
                        )
                        if both:
                            logger.error(
                                "Ligands marked both valid & invalid: "
                                + ", ".join(os.path.basename(x) for x in list(both)[:10])
                            )
                        if missing:
                            logger.error(
                                "Ligands missing from results: "
                                + ", ".join(os.path.basename(x) for x in list(missing)[:10])
                            )
                except Exception as _e:
                    logger.warning(f"Invariant check failed: {_e}")
    
                for lig, sc in scores.items():
                    record_score(score_history, stage['name'], lig, sc, True)
                    record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)
                for lig, (sc, reason) in invalids.items():
                    record_score(score_history, stage['name'], lig, sc, False, reason=reason)
                    record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)
    
                promoted_this_stage = False
                try:
                    decision = selector.consider_switch(
                        stage['name'], scores, validated, raw_docked, receptor_pdbqt, center, guard
                    )
                    if decision.promoted and decision.new_center is not None:
                        old = center
                        center = decision.new_center
                        promoted_this_stage = True
                        guard.mark_switch()
                        if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                            checkpoint_invalidate_from(
                                cfg,
                                paths.pdb_id,
                                stages_for_run,
                                start_index=i,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
                        logger.info(
                            f"[CENTER] Switched from {old} -> {center} ({decision.reason}, "
                            f"SwitchScore={decision.switchscore:.2f}) [global switch]"
                        )
                except Exception as e:
                    logger.warning(f"CenterSelector failed gracefully: {e}")
    
                if ph_label:
                    ph_log.info(
                        "[ph_ensemble.dock.scores] pdb_id=%s variant=%s ph=%s stage=%s valid=%d invalid=%d",
                        paths.pdb_id,
                        variant_label,
                        ph_label,
                        stage["name"],
                        len(scores),
                        len(invalids),
                    )
    
                if not promoted_this_stage:
                    restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                        i, scores, distances, box_size, center, stage1_original, recenter_attempts, params,
                        cfg, paths.pdb_id, receptor_pdbqt, logger, raw_docked, guard, control_anchor_hit
                    )
                    if restart:
                        ligands = redo_ligands
                        if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                            checkpoint_invalidate_from(
                                cfg,
                                paths.pdb_id,
                                stages_for_run,
                                start_index=0,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
                        i = 0
                        continue
    
                try:
                    if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                        med = (
                            float(np.median([d for d in distances if isinstance(d, (int, float))]))
                            if distances
                            else None
                        )
                        if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                            dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                            min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                            new_box = tuple(max(min_box, s - dec) for s in box_size)
                            if new_box != box_size:
                                logger.info(
                                    f"Adaptive shrink: median dist {med:.2f} A -> box {box_size} -> {new_box}"
                                )
                                box_size = new_box
                except Exception as _e:
                    logger.warning(f"Adaptive shrink skipped: {_e}")
    
                if i < len(stages_for_run) - 1:
                    if not scores:
                        restart, center, box_size, redo_ligands = fallback_recentering_if_empty(
                            cfg, paths.pdb_id, stage['name'], scores, raw_docked,
                            receptor_pdbqt, center, box_size, stage1_original, logger, guard, control_anchor_hit
                        )
                        if restart:
                            ligands = redo_ligands
                            if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                                checkpoint_invalidate_from(
                                    cfg,
                                    paths.pdb_id,
                                    stages_for_run,
                                    start_index=0,
                                    ph_label=ph_label,
                                    variant=variant_env or None,
                                )
                            i = 0
                            continue
                        else:
                            break
    
                    use_stage1_base = (docking_mode == "polypharmacology" and i == 1)
    
                    rescue = []
                    if i < len(stages_for_run) - 1:
                        for lig, (sc, reason) in invalids.items():
                            if sc is not None and "self_rmsd_" in str(reason).lower() and sc <= float(
                                    cfg.get("RESCUE_SELF_RMSD_SCORE_MAX", -8.0)):
                                rescue.append((sc, lig))
                        rescue = [lig for _, lig in sorted(rescue)[:int(cfg.get("RESCUE_SELF_RMSD_TOP_N", 10))]]
    
                    selected = select_ligands_for_next(
                        docking_mode,
                        i,
                        stages_for_run,
                        scores,
                        logger,
                        base_pool_n=(len(stage1_original) if use_stage1_base else None),
                        force_include=(forced_extracted_for_stage3 if use_stage1_base else None)
                    )
    
                    if rescue:
                        sel_set = set(selected)
                        rescue_unique = [r for r in rescue if r not in sel_set]
                        ligands = rescue_unique + selected
                    else:
                        ligands = selected
                    if not ligands:
                        logger.warning(f"No ligands selected for {stages_for_run[i + 1]['name']}; stopping.")
                        break
    
                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    try:
                        fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                        checkpoint_mark_done(
                            cfg,
                            paths.pdb_id,
                            stage["name"],
                            fp,
                            ph_label=ph_label,
                            variant=variant_env or None,
                        )
                    except Exception:
                        pass
    
                i += 1
    
            final_pose_validation_and_screenshots(
                cfg, paths.pdb_id, stages_for_run, receptor_pdbqt, center, validated_ligands_last,
                score_history, cleaned_pdb, docking_mode, logger, ph_label
            )
    
            csv_path = write_scores_csv(
                cfg,
                paths.pdb_id,
                score_history,
                ph_label=ph_label,
                variant=variant_env or None,
                csv_prefix=csv_prefix,
            )
            logger.info(
                "[Scores] ph_label=%s summary=%s",
                ph_label if ph_label else "base",
                csv_path,
            )
    
            try:
                summary = {
                    "pdb_id": paths.pdb_id,
                    "center": tuple(map(float, center)) if center else None,
                    "box_size": tuple(map(float, box_size)) if box_size else None,
                    "n_ligands_stage1": len(stage1_original),
                    "n_valid_last_stage": len(validated_ligands_last),
                    "switch_history": getattr(selector, "switch_history", []),
                    "global_switches": guard.global_switches,
                    "stages": [s["name"] for s in stages_for_run],
                    "ph_label": ph_label,
                }
                _write_audit_json(cfg, paths.pdb_id, summary, ph_label=ph_label, variant=variant_env or None)
            except Exception as _e:
                logger.warning(f"Audit JSON write failed: {_e}")
    
            try:
                update_manifest_for_protein_success(
                    cfg,
                    manifest_run_id or "",
                    paths.pdb_id,
                    variant_label,
                    time.time() - ph_start_ts,
                    ph_tag=ph_label,
                )
            except Exception:
                ph_log.warning(
                    "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=success",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    exc_info=True,
                )
        except Exception:
            try:
                update_manifest_for_protein_failure(
                    cfg,
                    manifest_run_id or "",
                    paths.pdb_id,
                    variant_label,
                    None,
                    ph_tag=ph_label,
                )
            except Exception:
                ph_log.warning(
                    "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=fail",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    exc_info=True,
                )
            raise


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
    test_mode = _resolve_test_mode(cfg)

    if test_mode != "fda+dud":

        logger.info(
            "[subrun] mode=%s run_mode=None csv_prefix='' stage_prefix='' (single subrun)",
            test_mode,
        )
        _run_ligand_pipeline_subrun(
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
                run_mode=None,
                csv_prefix="",
                stage_name_prefix="",
        )
        return

    logger.info(
        "[subrun] mode=fda+dud -> running DUD-only subrun then FDA-only subrun "
        "(stage_prefix='dud_' for DUD only)"
    )


    # fda+dud: 1) DUD-only sub-run, 2) FDA-only sub-run

    # DUD: use run_mode="dud", stage and CSV prefixes "dud_"
    _run_ligand_pipeline_subrun(
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
        run_mode="dud",
        csv_prefix="dud_",
        stage_name_prefix="dud_",
    )

    # FDA: run_mode="fda", no prefixes (standard behavior)
    _run_ligand_pipeline_subrun(
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
        run_mode="fda",
        csv_prefix="",
        stage_name_prefix="",
    )


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
    max_workers = int(cfg["MAX_PARALLEL_JOBS"])
    variant_env = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

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

    return scores, validated_ligands, all_distances, raw_docked_ligands, invalids


__all__ = [
    "_map_reason_to_category",
    "RetryManager",
    "run_one_stage",
    "process_one_protein",
    "_coerce_test_map",
    "_resolve_test_mode",
]
