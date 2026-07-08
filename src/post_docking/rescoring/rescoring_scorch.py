from __future__ import annotations

import logging
import sys
from typing import Callable, Optional

from pathlib import Path

from post_docking.rescoring import rescoring_scorch_support as _scorch_support_mod
from post_docking.rescoring import scorch_runtime_state as _scorch_runtime_state_mod
from post_docking.rescoring import scorch_stage_io as _scorch_stage_io_mod
from post_docking.rescoring import scorch_main as _scorch_main_mod

from post_docking.rescoring.scorch_scoring_bridge import (  # noqa: F401
    _aggregate_combo,
    _annotate_scorch_file,
    _collect_best_pose_per_base,
    _collect_combo_from_rel,
    _collect_stage_pdbqts,
    _combo_matches_filters,
    _compute_best_composites,
    _control_base_from_path,
    _decoy_engine_stage_dirs,
    _decoy_engine_stage_dirs_all,
    _decoy_engine_stage_dirs_legacy,
    _decoy_post_stage_dirs,
    _decoy_prefixes_from_test_mode,
    _decoy_stage_dirs,
    _discover_mode_dirs,
    _dock6_stage_dirs,
    _done_sentinel_path,
    _filter_combos,
    _filter_done_combos,
    _filter_done_combos_with_stats,
    _gnina_stage_dirs,
    _launch_pose_bust_async,
    _ledock_stage_dirs,
    _load_consensus_top_bases,
    _load_control_bases,
    _log_no_combos_after_filters,
    _log_score_csv_coverage,
    _mark_combo_done,
    _normalized_filter_token,
    _parse_test_mode_tokens,
    _pose_base_from_path,
    _posebusters_missing,
    _post_stage_dirs,
    _prep_missing,
    _resolve_decoy_prefix_override,
    _run_pose_bust,
    _run_prep_for_scorch,
    _score_csv_for_spec,
    _score_stage,
    _scorch_command,
    _select_top_bases_from_score_csv,
    _stage_priority,
    _vina_stage_dirs,
    annotate_scorch_t_scores,
    annotate_scorch_z_scores,
    discover_combos,
    stage_dir_candidates,
)

_preflight = _scorch_runtime_state_mod._preflight
_resolve_roots = _scorch_runtime_state_mod._resolve_roots
configure_logging = _scorch_runtime_state_mod.configure_logging
parse_args = _scorch_runtime_state_mod.parse_args
shutil = _scorch_runtime_state_mod.shutil


def main(*args, **kwargs):
    _scorch_main_mod._preflight = _preflight
    return _scorch_main_mod.main(*args, **kwargs)


_RUNTIME_STATE_EXPORTS = {
    "SCORCH_SCRIPT",
    "SCORCH_ENV",
    "SCORCH_ENV_PREFIX",
    "SCORCH_ROOT",
    "SCORCH_USE_MICROMAMBA",
    "SCORCH_PYTHON",
    "SCORCH_TOP_FRACTION_DEFAULT",
    "SCORCH_TOP_FRACTION_KEY",
    "SCORCH_DONE_DIR",
    "SCORCH_DONE_SENTINEL",
}

def __getattr__(name: str):
    if name in _RUNTIME_STATE_EXPORTS:
        return getattr(_scorch_runtime_state_mod, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

_find_consensus_csv_import: Optional[Callable[[Path], Optional[Path]]]
_rerank_consensus_with_scorch_import: Optional[Callable[..., bool]]
try:
    from post_docking.rescoring.rescore_reranker import (
        find_consensus_csv as _find_consensus_csv_loaded,
        rerank_consensus_with_scorch as _rerank_consensus_with_scorch_loaded,
    )

    _find_consensus_csv_import = _find_consensus_csv_loaded
    _rerank_consensus_with_scorch_import = _rerank_consensus_with_scorch_loaded
except Exception:  # pragma: no cover - optional dependency
    _find_consensus_csv_import = None
    _rerank_consensus_with_scorch_import = None

find_consensus_csv_fn = _find_consensus_csv_import
rerank_consensus_with_scorch_fn = _rerank_consensus_with_scorch_import

COMPONENT = "[scorch-rescore]"
DECOY_PREFIX_KEY = _scorch_support_mod.DECOY_PREFIX_KEY
DUD_PREFIX_KEY = _scorch_support_mod.DUD_PREFIX_KEY
SCORCH_CHUNK_SIZE = _scorch_support_mod.SCORCH_CHUNK_SIZE
SCORCH_CHUNK_MIN = _scorch_support_mod.SCORCH_CHUNK_MIN
SCORCH_TAIL_SPLIT_TRIGGER = _scorch_support_mod.SCORCH_TAIL_SPLIT_TRIGGER
SCORCH_RUNTIME_RECHUNK_STALE_SEC = (
    _scorch_support_mod.SCORCH_RUNTIME_RECHUNK_STALE_SEC
)
SCORCH_HEDGE_TIMEOUT_SEC = _scorch_support_mod.SCORCH_HEDGE_TIMEOUT_SEC
SCORCH_HEDGE_MAX_INFLIGHT = _scorch_support_mod.SCORCH_HEDGE_MAX_INFLIGHT
SCORCH_HEDGE_MIN_ALLOWED = _scorch_support_mod.SCORCH_HEDGE_MIN_ALLOWED
SCORCH_PARALLEL_PROFILE_DEFAULT = _scorch_support_mod.SCORCH_PARALLEL_PROFILE_DEFAULT
SCORCH_PARALLEL_PROFILE_ENV = _scorch_support_mod.SCORCH_PARALLEL_PROFILE_ENV
_normalize_parallel_profile = _scorch_support_mod._normalize_parallel_profile
_util_target_for_profile = _scorch_support_mod._util_target_for_profile
_allocate_scorch_parallelism = _scorch_support_mod._allocate_scorch_parallelism
_decoy_prefix_value = _scorch_support_mod._decoy_prefix_value
_normalize_decoy_prefix = _scorch_support_mod._normalize_decoy_prefix
_resolve_decoy_prefix_from_config = _scorch_support_mod._resolve_decoy_prefix_from_config
_parse_test_mode_value = _scorch_support_mod._parse_test_mode_value
discover_stage3_roots = _scorch_support_mod.discover_stage3_roots
_prefer_existing = _scorch_support_mod._prefer_existing
_scheduler_runtime_snapshot = _scorch_support_mod._scheduler_runtime_snapshot
_adaptive_chunk_size = _scorch_support_mod._adaptive_chunk_size
_tail_split_allowed_chunks = _scorch_support_mod._tail_split_allowed_chunks
_chunk_allowed_bases = _scorch_support_mod._chunk_allowed_bases
_task_allowed_count = _scorch_support_mod._task_allowed_count
_task_estimate_seconds = _scorch_support_mod._task_estimate_seconds
_elastic_scorch_threads = _scorch_support_mod._elastic_scorch_threads
_split_task_for_rechunk = _scorch_support_mod._split_task_for_rechunk
DECOY_PREFIX_VALUE = _scorch_support_mod.DECOY_PREFIX_VALUE
_materialize_inputs = _scorch_stage_io_mod.materialize_inputs
_annotate_csv = _scorch_stage_io_mod.annotate_csv


def _set_decoy_prefix(value: str, logger: Optional[logging.Logger] = None) -> str:
    global DECOY_PREFIX_VALUE
    DECOY_PREFIX_VALUE = _scorch_support_mod._set_decoy_prefix(value, logger)
    return DECOY_PREFIX_VALUE


if __name__ == "__main__":
    sys.exit(main())
