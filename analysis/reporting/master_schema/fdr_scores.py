# -*- coding: utf-8 -*-
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

from analysis.reporting.value_utils import as_float
from analysis.reporting.master_schema.constants import DECOY_PREFIX_DEFAULT, MIN_DECOYS_FOR_FDR
from analysis.reporting.master_schema.decoy_helpers import _row_is_decoy
from analysis.reporting.master_schema.io_utils import _read_csv_rows

def _candidate_fdr_score_fields() -> Tuple[str, ...]:
    return (
        "final_score",
        "z_vs_decoys_blend",
        "z_selected",
        "z_vs_decoys_consensus",
        "z_vs_decoys_consensus_pre",
        "consensus_score",
        "consensus_score_pre",
    )


def _choose_fdr_score_field(
    decoy_rows: List[Dict[str, Any]], min_decoys: int = MIN_DECOYS_FOR_FDR
) -> str:
    counts: Dict[str, int] = {}
    for field in _candidate_fdr_score_fields():
        vals = []
        for r in decoy_rows:
            v = as_float(r.get(field))
            if v is None or not math.isfinite(v):
                continue
            vals.append(v)
        counts[field] = len(vals)
    preferred_field = None
    preferred_count = -1
    for field, count in counts.items():
        if count >= min_decoys and count > preferred_count:
            preferred_field = field
            preferred_count = count
    if preferred_field:
        return preferred_field

    if counts:
        fallback_field = max(counts.items(), key=lambda kv: kv[1])[0]
        if counts[fallback_field] > 0:
            return fallback_field

    return "final_score"


def _get_decoy_scores_for_combo(
    combo_dir: Path,
    fallback_rows: List[Dict[str, Any]],
    min_decoys: int = 1,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
) -> Tuple[str, List[float], bool, Dict[str, List[float]], str]:
    primary_path = combo_dir / f"{decoy_prefix}_consensus_reranked_scorch.csv"
    legacy_path = combo_dir / "dud_consensus_reranked_scorch.csv"
    fallback_used = False
    decoy_rows: List[Dict[str, Any]] = []
    source = ""
    if primary_path.exists() and primary_path.stat().st_size > 0:
        decoy_rows = _read_csv_rows(primary_path)
        source = "post_scorch_decoy"
    if not decoy_rows and legacy_path.exists() and legacy_path.stat().st_size > 0:
        decoy_rows = _read_csv_rows(legacy_path)
        source = "post_scorch_legacy_decoy"
    if not decoy_rows:
        fallback_used = True
        decoy_rows = fallback_rows
        source = "master_rows_fallback"
    decoys_only = [r for r in decoy_rows if _row_is_decoy(r, decoy_prefix)]
    if not decoys_only:
        fallback_used = True
        decoys_only = [r for r in fallback_rows if _row_is_decoy(r, decoy_prefix)]
        source = "master_rows_fallback"
    scores_by_field: Dict[str, List[float]] = {}
    for field_name in _candidate_fdr_score_fields():
        vals: List[float] = []
        for r in decoys_only:
            val = as_float(r.get(field_name))
            if val is None or not math.isfinite(val):
                continue
            vals.append(val)
        scores_by_field[field_name] = vals

    field = _choose_fdr_score_field(decoys_only, min_decoys)
    scores = scores_by_field.get(field, [])
    return field, scores, fallback_used, scores_by_field, source


def _docking_decoy_score_paths(
    docking_combo_dir: Path, decoy_prefix: str
) -> List[Path]:
    prefixes: List[str] = []
    for prefix in (decoy_prefix, "fda_dud", "dud"):
        if prefix and prefix not in prefixes:
            prefixes.append(prefix)
    paths: List[Path] = []
    for prefix in prefixes:
        paths.append(docking_combo_dir / f"{prefix}_consensus_docking_scores.csv")
    return paths


def _get_docking_decoy_scores_for_combo(
    docking_combo_dir: Path,
    min_decoys: int = 1,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
) -> Tuple[str, List[float], Dict[str, List[float]], str]:
    decoy_rows: List[Dict[str, Any]] = []
    source_path = ""
    for path in _docking_decoy_score_paths(docking_combo_dir, decoy_prefix):
        if not path.exists() or path.stat().st_size <= 0:
            continue
        rows = _read_csv_rows(path)
        decoys_only = [r for r in rows if _row_is_decoy(r, decoy_prefix)]
        if decoys_only:
            decoy_rows = decoys_only
            source_path = str(path)
            break
    scores_by_field: Dict[str, List[float]] = {}
    for field_name in _candidate_fdr_score_fields():
        vals: List[float] = []
        for r in decoy_rows:
            val = as_float(r.get(field_name))
            if val is None or not math.isfinite(val):
                continue
            vals.append(val)
        scores_by_field[field_name] = vals
    field = _choose_fdr_score_field(decoy_rows, min_decoys)
    scores = scores_by_field.get(field, [])
    return field, scores, scores_by_field, source_path
