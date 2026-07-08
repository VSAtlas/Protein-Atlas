"""
Best-effort run manifest writer for Atlas runs.

All helpers must remain non-fatal: any error (missing PyYAML, IO, git)
is logged at WARNING level and silently skipped so scientific behavior
is unchanged.
"""

from __future__ import annotations

from cli.run_manifest_distributed import (
    reduce_distributed_protein_states,
    reconcile_distributed_manifest_state,
)
from cli.run_manifest_io import (
    get_manifest_paths,
    load_manifest as _load_manifest,
    protein_key as _protein_key,
)
from cli.run_manifest_support import (
    PocketDetectionEvent,
    compute_config_hash,
    refresh_summary as _refresh_summary,
)
from cli.run_manifest_updates_lifecycle import (
    finalize_run_manifest,
    init_run_manifest,
    load_run_manifest,
    update_manifest_for_config_hash,
    update_manifest_for_run_config,
    update_manifest_for_scheduled_proteins,
)
from cli.run_manifest_updates_protein import (
    update_manifest_for_druggability_and_engine_plan,
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
    update_manifest_for_protein_success,
)
from cli.run_manifest_updates_stages import (
    apply_pocket_detection_events,
    emit_pocket_detection_event,
    update_manifest_for_docking_overall,
    update_manifest_for_docking_stage,
    update_manifest_for_pocket_detection,
    update_manifest_for_postprocessing,
    update_manifest_for_prep_stage,
    update_manifest_for_stage,
)

__all__ = [
    "PocketDetectionEvent",
    "apply_pocket_detection_events",
    "compute_config_hash",
    "emit_pocket_detection_event",
    "finalize_run_manifest",
    "get_manifest_paths",
    "init_run_manifest",
    "load_run_manifest",
    "reconcile_distributed_manifest_state",
    "reduce_distributed_protein_states",
    "update_manifest_for_config_hash",
    "update_manifest_for_docking_overall",
    "update_manifest_for_docking_stage",
    "update_manifest_for_druggability_and_engine_plan",
    "update_manifest_for_pocket_detection",
    "update_manifest_for_postprocessing",
    "update_manifest_for_prep_stage",
    "update_manifest_for_protein_failure",
    "update_manifest_for_protein_start",
    "update_manifest_for_protein_success",
    "update_manifest_for_run_config",
    "update_manifest_for_scheduled_proteins",
    "update_manifest_for_stage",
    "_load_manifest",
    "_protein_key",
    "_refresh_summary",
]

