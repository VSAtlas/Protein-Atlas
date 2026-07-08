from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from docking.docking_subrun_selection import _use_ledock, _use_dock6
from docking.druggability_policy import decide_engine_policy
from path_router.path_router import Paths, docked_dir
from cli.run_manifest_runtime import update_manifest_for_druggability_and_engine_plan


def maybe_handle_no_library_docking(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    variant_env: str,
    variant_token: Optional[str],
    variant_label: str,
    legacy_mode: bool,
    ph_label: str,
    run_mode: Optional[str],
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    stage1_original: List[str],
) -> bool:
    """
    Checks NO_LIBRARY_DOCKING flag. If true, runs plan-only logic (druggability, engine policy),
    writes preview file, and returns True (indicating the caller should skip actual docking).
    """
    if not bool(cfg.get("NO_LIBRARY_DOCKING", False)):
        return False

    run_id_token = str(cfg.get("RUN_ID") or "")
    if run_id_token:
        try:
            policy = decide_engine_policy(
                cfg=cfg,
                pdb_id=paths.pdb_id,
                variant=variant_env or None,
                ph_label=ph_label,
                center=center,
                ledock_enabled=_use_ledock(cfg),
                dock6_enabled=_use_dock6(cfg),
                logger=logger,
            )
            update_manifest_for_druggability_and_engine_plan(
                cfg=cfg,
                run_id=run_id_token,
                pdb_id=paths.pdb_id,
                variant_label=variant_env or None,
                ph_tag=ph_label,
                tier=policy.tier,
                use_gnina=policy.use_gnina,
                use_ledock=policy.use_ledock,
                use_dock6=policy.use_dock6,
            )
        except Exception:
            logger.warning(
                "[run-manifest.druggability-plan.skip] run_id=%s pdb=%s variant=%s ph=%s",
                run_id_token,
                paths.pdb_id,
                variant_env,
                ph_label,
                exc_info=True,
            )

    # How many ligands would be docked for this PH / mode?
    preview_n = int(cfg.get("NO_DOCKING_PREVIEW_N", 10))
    preview_names = ", ".join(
        Path(ligand_path).name for ligand_path in stage1_original[:preview_n]
    )

    logger.info(
        "[no-docking-planned] pdb=%s variant=%s ph=%s mode=%s n_stage1=%d preview=[%s]",
        paths.pdb_id,
        variant_label,
        ph_label if ph_label else "base",
        run_mode or "(unspecified)",
        len(stage1_original),
        preview_names,
    )

    # Also drop a simple text file with the full Stage1 ligand list
    try:
        out_dir = docked_dir(
            paths.pdb_id,
            variant=variant_token,
            ph_tag=ph_label,
            legacy=legacy_mode,
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        suffix = run_mode or "run"
        out_txt = out_dir / f"planned_ligands_{suffix}.txt"

        with out_txt.open("w") as fh:
            for lig in stage1_original:
                fh.write(f"{lig}\n")

        logger.info(
            "[no-docking-planned] Wrote planned ligands to %s",
            out_txt,
        )
    except Exception as e:
        logger.warning(
            "[no-docking-planned] Failed to write planned ligand list: %s",
            e,
        )

    # Skip all docking stages for this PH context.
    return True
