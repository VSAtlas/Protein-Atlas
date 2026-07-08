"""Thin orchestration facade; implementation split into orchestrate_*.py helpers."""

from analysis.dud_eval_core.orchestrate_cli import main
from analysis.dud_eval_core.orchestrate_discovery import (
    CONTROL_CENTERS_RE,
    CONTROL_REDOCK_RE,
    FULL_RUN_MIN_LIGANDS,
    _CONTROL_PATTERNS_LOGGED,
    _build_control_records,
    _candidate_protein_logs,
    _compute_analysis_root,
    _format_run_label,
    _load_default_cfg,
    _manifest_library_for_target,
    _manifest_protein_entries,
    _manifest_proteins_by_pdb,
    _resolve_run_label,
    select_default_run_id,
)
from analysis.dud_eval_core.orchestrate_metrics import (
    DUD_E_PROTEIN_CLASSES,
    _emit_class_aggregate_plots,
)

__all__ = [
    "CONTROL_CENTERS_RE",
    "CONTROL_REDOCK_RE",
    "DUD_E_PROTEIN_CLASSES",
    "FULL_RUN_MIN_LIGANDS",
    "_CONTROL_PATTERNS_LOGGED",
    "_build_control_records",
    "_candidate_protein_logs",
    "_compute_analysis_root",
    "_emit_class_aggregate_plots",
    "_format_run_label",
    "_load_default_cfg",
    "_manifest_library_for_target",
    "_manifest_protein_entries",
    "_manifest_proteins_by_pdb",
    "_resolve_run_label",
    "main",
    "select_default_run_id",
]
