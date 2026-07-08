"""Phase helpers for docking runtime orchestration."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, cast

from cli.runtime_logging import make_protein_logger
from path_router.path_router import Paths, make_paths
from docking.docking_controls import _summarize_ions_file
from docking.docking_subruns import (
    ProteinDockingContext,
    SubrunSpec,
    subruns_for_tokens,
    run_ligand_pipeline_subrun,
)
from docking.fallback_recenter import RecenterParams
from docking.library_mode import parse_test_libraries


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


def _record_input_ion_audit(paths: Paths, pdb_audit: dict, logger: logging.Logger) -> None:
    input_summary = cast(dict[str, Any], _summarize_ions_file(paths.input_pdb_path))
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


def _log_receptor_targets(
    paths: Paths,
    logger: logging.Logger,
    variant_label: str,
    cleaned_target: Path,
    receptor_target: Path,
) -> None:
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

    _record_input_ion_audit(paths, pdb_audit, logger)

    cfg["_CURRENT_VARIANT"] = variant_env
    cleaned_target = paths.receptor_cleaned_pdb(variant_token)
    receptor_target = paths.receptor_pdbqt(variant_token, None)
    _log_receptor_targets(paths, logger, variant_label, cleaned_target, receptor_target)
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
    chunk_run_mode = str(cfg.get("_CHUNK_RUN_MODE") or "").strip().lower()
    if cfg.get("_CHUNK_LIGAND_KEYS") and chunk_run_mode:
        subruns = [
            sub
            for sub in subruns
            if str(sub.run_mode or "fda").strip().lower() == chunk_run_mode
        ]
        if not subruns:
            raise RuntimeError(
                "chunk_run_mode_not_in_test_libraries:"
                f"run_mode={chunk_run_mode}:tokens={'+'.join(tokens)}"
            )

    for sub in subruns:
        logger.info(
            "[subrun] tokens=%s run_mode=%r csv_prefix=%r stage_prefix=%r",
            "+".join(tokens),
            sub.run_mode,
            sub.csv_prefix,
            sub.stage_name_prefix,
        )
        run_ligand_pipeline_subrun(ctx, sub)
