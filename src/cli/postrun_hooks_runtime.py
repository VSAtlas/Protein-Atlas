# -*- coding: utf-8 -*-
from __future__ import annotations

import subprocess  # noqa: F401  # tests monkeypatch postrun_hooks.subprocess

from cli.postrun_hooks_common import _is_no_library_docking, _write_json_atomic
from cli.postrun_hooks_integrity import _maybe_run_throughput_integrity
from cli.postrun_hooks_report import (
    _collect_sacct_cpu_metrics,
    _manifest_indicates_dud_run,
    _maybe_run_dud_eval,
    _maybe_run_master_schema_export,
    _maybe_run_report_generation,
    _maybe_run_static_heatmaps,
    _write_run_efficiency_report,
)
from cli.postrun_hooks_retention import (
    _as_bool,
    _invoke_artifact_retention,
    _maybe_run_artifact_retention,
    _maybe_run_artifact_retention_for_combo,
    _maybe_run_artifact_retention_for_combos,
    _maybe_run_artifact_retention_for_pdb,
    _normalize_artifact_retention_mode,
    _to_globs,
)
from cli.postrun_hooks_scorch import (
    _log_rescore_verification,
    _maybe_run_scorch_rescore,
    _maybe_run_scorch_rescore_for_pdb,
    _run_provisional_scorch_rescore_for_pdb,
)
from cli.postrun_hooks_scorch_queue import (
    ScorchMixedQueueService,
    _start_scorch_mixed_queue_service,
)
from cli.postrun_hooks_support import (
    _resolve_scorch_execution_mode,
    _scorch_cmd_prefix,
    _scorch_subprocess_env,
)
from cli.run_manifest_runtime import update_manifest_for_postprocessing

__all__ = [
    "ScorchMixedQueueService",
    "_as_bool",
    "_collect_sacct_cpu_metrics",
    "_invoke_artifact_retention",
    "_is_no_library_docking",
    "_log_rescore_verification",
    "_manifest_indicates_dud_run",
    "_maybe_run_artifact_retention",
    "_maybe_run_artifact_retention_for_combo",
    "_maybe_run_artifact_retention_for_combos",
    "_maybe_run_artifact_retention_for_pdb",
    "_maybe_run_dud_eval",
    "_maybe_run_master_schema_export",
    "_maybe_run_report_generation",
    "_maybe_run_scorch_rescore",
    "_maybe_run_scorch_rescore_for_pdb",
    "_maybe_run_static_heatmaps",
    "_maybe_run_throughput_integrity",
    "_normalize_artifact_retention_mode",
    "_resolve_scorch_execution_mode",
    "_run_provisional_scorch_rescore_for_pdb",
    "_scorch_cmd_prefix",
    "_scorch_subprocess_env",
    "_start_scorch_mixed_queue_service",
    "_to_globs",
    "update_manifest_for_postprocessing",
    "_write_json_atomic",
    "_write_run_efficiency_report",
]

