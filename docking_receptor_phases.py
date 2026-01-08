# docking_receptor_phases.py

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import automate_protein_prep as protein_prep
from apo_holo_mode import (
    _record_apo_holo_decision,
    _record_apo_holo_usage,
    _variant_receptor_path,
    delete_variant_trees,
    file_sha1,
)
from docking_centering import CenterSelector
from docking_controls import _summarize_ions_file
from docking_receptor import prepare_receptor
from docking_utils import norm
from path_router import Paths
from protein_functions import detect_active_site


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
    """
    Orchestrates receptor preparation, active site detection, ion audit/strip,
    holo-restore (rebuilding if needed), and APO vs HOLO preflight checks.
    """
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
        probe_map = protein_prep.get_ion_probe_map(paths.pdb_id)
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

    c2, b2, src = detect_active_site(cleaned_pdb)
    if c2:
        box_size = tuple(min(28.0, float(s)) for s in b2)
        center = c2
        center_source = src or "activesite"
        logger.info("[active-site] Using center %s with box %s source=%s", center, box_size, center_source)
    else:
        logger.error(
            "[active-site] detection_failed pdb=%s variant=%s action=continue",
            paths.pdb_id,
            variant_label,
        )

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

    if box_size is not None:
        box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
        box_size = tuple(min(box_cap, float(s)) for s in box_size)
        logger.info(f"Initial box clamped to {box_size} (cap={box_cap} A)")

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
