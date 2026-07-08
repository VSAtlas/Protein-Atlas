from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Mapping

from cli.distributed_chunk_planner import (
    _chunk_ligand_key,
    _load_post_scored_ligand_keys,
    _load_scored_ligand_keys_from_summary,
    _normalize_ph_tag_token,
    _resolve_combo_docking_summary_csv,
    _resolve_combo_output_dir,
    _resolve_combo_post_consensus_csv,
    _resolve_global_scheduler_plan,
    _update_combo_coverage_snapshot,
    _verify_chunk_combo_outputs,
)
import post_docking.mmgbsa.mmgbsa_pipeline as _mmgbsa_pipeline
from post_docking.mmgbsa.mmgbsa_pipeline import (
    make_mmgbsa_trajectory,
    run_implicit_md,
    run_mmgbsa,
)


def run_mmgbsa_five_replicate_with_hooks(
    *,
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
    run_implicit_md_hook,
    run_mmgbsa_hook,
    make_mmgbsa_trajectory_hook,
) -> Dict[str, object]:
    # Compatibility wrapper only.
    # No new logic here.
    # Use src/post_docking/mmgbsa and src/cli owner modules for new code.
    _mmgbsa_pipeline.run_implicit_md = run_implicit_md_hook
    _mmgbsa_pipeline.run_mmgbsa = run_mmgbsa_hook
    _mmgbsa_pipeline.make_mmgbsa_trajectory = make_mmgbsa_trajectory_hook
    return _mmgbsa_pipeline._mmgbsa_five_replicate_runner(
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


__all__ = [
    "_chunk_ligand_key",
    "_load_post_scored_ligand_keys",
    "_load_scored_ligand_keys_from_summary",
    "_normalize_ph_tag_token",
    "_resolve_combo_docking_summary_csv",
    "_resolve_combo_output_dir",
    "_resolve_combo_post_consensus_csv",
    "_resolve_global_scheduler_plan",
    "_update_combo_coverage_snapshot",
    "_verify_chunk_combo_outputs",
    "make_mmgbsa_trajectory",
    "run_implicit_md",
    "run_mmgbsa",
    "run_mmgbsa_five_replicate_with_hooks",
]

