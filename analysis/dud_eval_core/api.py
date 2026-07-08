#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ruff: noqa: F401

"""
DUD/DUD-E style evaluator for Atlas docking outputs (filename-labeled actives/decoys).
"""

from analysis.dud_eval_core.discovery import (
    _LIB_INDEX_CACHE,
    _dedup,
    _get_library_index,
    _infer_library_name_via_manifest,
    _is_probable_run_id_dirname,
    _is_under,
    _normalize_pdb_ids,
    _read_pdb_header_lines,
    _resolve_reranked_scorch_path,
    _resolve_scan_roots,
    choose_target_name,
    derive_target_name,
    extract_compnd_molecules,
    extract_uniprot_from_dbref,
    infer_library_name,
    make_target_key,
)
from analysis.dud_eval_core.engine import (
    emit_consensus_summary,
    emit_reranked_scorch_summary,
    evaluate_target,
    evaluate_target_consensus,
    evaluate_target_post_docked_reranked_scorch,
)
from analysis.dud_eval_core.labels import (
    EXT_RE,
    MICROSTATE_TAIL_RE,
    POSE_TAIL_RE,
    TOKEN_RE,
    _collapse_best_scores,
    _make_placeholder_row,
    compute_decoy_stats_from_long_csv,
    parse_name_and_label,
)
from analysis.dud_eval_core.log import (
    LEVEL_ABBREV,
    LOG_LEVELS,
    _BACKEND_LOGGED,
    _current_log_level,
    dbg,
    set_log_level,
)
from analysis.dud_eval_core.metrics import bedroc, ef_at_fractions, log_auc_from_roc, pr_auc
from analysis.dud_eval_core.orchestrate import (
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
    main,
    select_default_run_id,
)
from analysis.dud_eval_core.reporting import (
    _compute_metrics_and_plots,
    _df_to_pretty_text,
    _fmt_counts_and_round,
    _write_pretty_summary,
    _write_pretty_table_noformat,
)
from analysis.dud_eval_core.schema import (
    CONSENSUS_SCORE_CANDIDATES,
    LIGFILE_CANDIDATES,
    RUN_ID_COL_CANDIDATES,
    SCORE_CANDIDATES,
    VALID_COL_CANDIDATES,
    VALID_TRUE_STRINGS,
    _filter_consensus_no_data,
    _normalize_run_id_token,
    filter_df_by_run_id,
    guess_col,
    guess_consensus_ligfile_col,
    guess_consensus_score_col,
    guess_ligfile_col,
    guess_reranked_scorch_score_col,
    guess_score_col,
    parse_valid_mask,
    read_reranked_scorch_csv,
    resolve_valid_col,
)
from analysis.dud_eval_core.types import (
    CONSENSUS_CSV_BASENAME,
    CSV_BASENAMES,
    RERANKED_SCORCH_BASENAME,
    TargetEvaluation,
    TargetSpec,
)
