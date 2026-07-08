from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional

import pandas as pd

from analysis.dud_eval_core.log import dbg

SCORE_CANDIDATES = [
    "score",
    "docking_score",
    "vina_score",
    "affinity",
    "pose_score",
    "dockscore",
    "energy",
    "gnina_score",
    "cnnscore",
    "cnn_affinity",
]

LIGFILE_CANDIDATES = [
    "ligand_file",
    "ligand_path",
    "ligand",
    "ligand_name",
    "ligand",
    "pose_file",
    "pose_path",
    "output_ligand",
    "output_ligand_path",
    "file",
    "filepath",
    "filename",
    "pdbqt_path",
    "pdbqt",
    "out_path",
]

VALID_COL_CANDIDATES = [
    "valid",
    "pose_valid",
    "is_valid",
    "passed_pose_validation",
    "pose_is_valid",
    "valid_pose",
    "pose_ok",
    "ok",
    "passed",
    "pass",
]

VALID_TRUE_STRINGS = {"true", "t", "yes", "y", "pass", "passed", "ok", "valid"}

CONSENSUS_SCORE_CANDIDATES = ["consensus_score", "score", "consensus", "final_score"]

RUN_ID_COL_CANDIDATES = ["run_id", "runid", "run"]

def guess_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    lower = {c.lower(): c for c in df.columns}
    for key in candidates:
        if key in lower:
            return lower[key]
    for frag in candidates:
        for c in df.columns:
            if frag in c.lower():
                return c
    return None

def guess_score_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    col = guess_col(df, SCORE_CANDIDATES)
    if col:
        return col
    # last resort: first numeric column that smells like score
    for c in df.columns:
        if df[c].dtype.kind in "fi" and any(
            k in c.lower() for k in ("score", "affin", "dock", "energy")
        ):
            return c
    raise ValueError("Could not find score column. Use --score-col.")

def guess_consensus_score_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    lower = {c.lower(): c for c in df.columns}
    for key in ("consensus_score", "score"):
        if key in lower:
            return lower[key]
    for cand in CONSENSUS_SCORE_CANDIDATES:
        if cand in lower:
            return lower[cand]
    col = guess_col(df, CONSENSUS_SCORE_CANDIDATES)
    if col:
        return col
    return guess_score_col(df, override=None)

def guess_reranked_scorch_score_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    lower = {c.lower(): c for c in df.columns}

    def _has_numeric(col: str) -> bool:
        if col not in df.columns:
            return False
        try:
            ser = pd.to_numeric(df[col], errors="coerce")
            return ser.notna().any()
        except Exception:
            return False

    if "final_score" in df.columns and _has_numeric("final_score"):
        return "final_score"
    if "final_score" in lower and _has_numeric(lower["final_score"]):
        return lower["final_score"]

    for cand in (
        "final_rank",
        "scorch_composite",
        "SCORCH_score_used",
        "consensus_score",
    ):
        if cand in df.columns:
            return cand
    for cand in ("final_rank", "scorch_score_used", "consensus_score"):
        if cand in lower:
            return lower[cand]
    return guess_score_col(df, None)

def read_reranked_scorch_csv(path: Path) -> pd.DataFrame:
    """
    Read a reranked SCORCH CSV, tolerating optional metadata preamble lines.

    The file may start with comment-style metadata (e.g., '# run_id=...').
    We first attempt a normal read; if required columns are missing, re-read
    after skipping metadata lines.
    """

    def _detect_preamble(lines: List[str]) -> int:
        count = 0
        meta_prefixes = ("run_id", "pdb_id", "variant", "ph", "ph_label")
        for line in lines:
            stripped = line.strip()
            if not stripped:
                count += 1
                continue
            if stripped.startswith("#"):
                count += 1
                continue
            low = stripped.lower()
            if any(
                low.startswith(f"{p}=") or low.startswith(f"{p}:")
                for p in meta_prefixes
            ):
                count += 1
                continue
            break
        return count

    try:
        return pd.read_csv(path)
    except Exception:
        pass

    try:
        with path.open("r", encoding="utf-8") as fh:
            head: List[str] = []
            for _ in range(5):
                try:
                    head.append(next(fh))
                except StopIteration:
                    break
    except Exception:
        return pd.read_csv(path)

    skip = _detect_preamble(head)
    if skip <= 0:
        return pd.read_csv(path)
    return pd.read_csv(path, skiprows=skip)

def guess_ligfile_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    col = guess_col(df, LIGFILE_CANDIDATES)
    if col:
        return col
    # fallback: first object column
    for c in df.columns:
        if df[c].dtype == object:
            return c
    raise ValueError("Could not find ligand filename/path column. Use --lig-col.")

def guess_consensus_ligfile_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    lower = {c.lower(): c for c in df.columns}
    for key in ("ligand", "ligand_file"):
        if key in lower:
            return lower[key]
    return guess_ligfile_col(df, override=None)

def resolve_valid_col(
    df: pd.DataFrame, override: Optional[str], enabled: bool
) -> Optional[str]:
    if not enabled:
        return None
    if override:
        if override in df.columns:
            return override
        raise ValueError(
            f"Valid-only requested but --valid-col '{override}' is missing. "
            "Pass --valid-col with an existing column name."
        )
    for name in VALID_COL_CANDIDATES:
        if name in df.columns:
            return name
    raise ValueError(
        "Valid-only requested but no validity column found. Use --valid-col."
    )

def parse_valid_mask(series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        numeric = pd.to_numeric(series, errors="coerce")
        return numeric.notna() & (numeric != 0)
    text = series.fillna("").astype(str).str.strip()
    numeric = pd.to_numeric(text, errors="coerce")
    numeric_mask = numeric.notna() & (numeric != 0)
    token_mask = text.str.lower().isin(VALID_TRUE_STRINGS)
    return numeric_mask | token_mask

def _normalize_run_id_token(val: str) -> str:
    """Normalize run_id tokens for comparison."""
    if val is None:
        return ""
    text = str(val).strip().strip("'\"")
    return re.sub(r"[^0-9A-Za-z]+", "", text).lower()

def filter_df_by_run_id(
    df: pd.DataFrame,
    cli_run_id: Optional[str],
    *,
    pdb_id: str,
    csv_path: Path,
    strict: bool,
    layout_label: str = "legacy",
) -> tuple[pd.DataFrame, str]:
    """
    Filter a dataframe by run_id with robust normalization.

    Returns (filtered_df, status_reason).
    status_reason is "ok" on match, or a mismatch/fallback label otherwise.
    """
    if not cli_run_id:
        return df, "ok"

    run_col = None
    for cand in RUN_ID_COL_CANDIDATES:
        if cand in df.columns:
            run_col = cand
            break
    if run_col is None:
        dbg("WARN", "run.filter", f"pdb={pdb_id} path={csv_path} reason=no_run_id_col")
        return df, "ok_no_run_id_col"

    target_norm = _normalize_run_id_token(cli_run_id)
    series_raw = df[run_col].astype(str).str.strip().str.strip("'\"")
    series_norm = series_raw.str.replace(r"[^0-9A-Za-z]+", "", regex=True).str.lower()
    mask = series_norm == target_norm
    kept = int(mask.sum())
    total = len(df)
    if kept > 0:
        dbg(
            "INFO",
            "run.filter",
            f"pdb={pdb_id} path={csv_path} run_id={cli_run_id} kept={kept}/{total} col={run_col} layout={layout_label}",
        )
        return df.loc[mask].copy(), "ok"

    value_counts = series_norm.value_counts(dropna=False)
    sample = value_counts.head(5).to_dict()
    dbg(
        "WARN",
        "run.filter",
        f"pdb={pdb_id} path={csv_path} run_id={cli_run_id} matched=0/{total} layout={layout_label} sample={sample}",
    )

    if strict:
        dbg(
            "ERROR",
            "run.filter",
            f"pdb={pdb_id} path={csv_path} run_id={cli_run_id} matched=0 action=abort layout={layout_label}",
        )
        return df.loc[mask].copy(), "run_id_mismatch_new_layout"

    if len(value_counts) == 1:
        dbg(
            "WARN",
            "run.filter",
            f"pdb={pdb_id} layout=legacy_fallback action=use_singleton run_id_val={value_counts.index[0]} rows={total}",
        )
        return df.copy(), "run_id_mismatch_legacy_fallback"

    top_value = value_counts.index[0]
    top_rows = int(value_counts.iloc[0])
    dbg(
        "WARN",
        "run.filter",
        f"pdb={pdb_id} layout=legacy_fallback action=use_top run_id_val={top_value} rows={top_rows}/{total}",
    )
    top_mask = series_norm == top_value
    return df.loc[top_mask].copy(), "run_id_mismatch_legacy_fallback"

def _filter_consensus_no_data(
    df: pd.DataFrame, score_col: str
) -> tuple[pd.DataFrame, int]:
    """
    Drop rows that clearly have no engine data.

    Priority:
      1) n_engines_with_data > 0
      2) Boolean-like validity column (valid/ok/has_data)
      3) All engine columns zero (excluding consensus_score + obvious metadata)
    """
    start = len(df)
    if start == 0:
        return df, 0

    if "n_engines_with_data" in df.columns:
        counts = pd.to_numeric(df["n_engines_with_data"], errors="coerce").fillna(0)
        mask = counts > 0
        filtered = df.loc[mask].copy()
        dropped = start - len(filtered)
        dbg(
            "INFO",
            "consensus.filter",
            f"method=n_engines_with_data kept={len(filtered)}/{start}",
        )
        return filtered, dropped

    bool_candidates = [
        c
        for c in df.columns
        if c.lower() in {"valid", "ok", "has_data", "usable", "is_valid", "available"}
    ]
    for col in bool_candidates:
        try:
            mask = parse_valid_mask(df[col])
        except Exception:
            continue
        filtered = df.loc[mask].copy()
        dropped = start - len(filtered)
        dbg(
            "INFO",
            "consensus.filter",
            f"method=bool_col col={col} kept={len(filtered)}/{start}",
        )
        return filtered, dropped

    numeric_cols = list(df.select_dtypes(include=["number", "bool"]).columns)
    exclude_tokens = {"rank", "percent", "perc", "stage", "pose", "order"}
    engine_cols: List[str] = []
    for col in numeric_cols:
        low = col.lower()
        if low == score_col.lower():
            continue
        if low in {"n_engines_with_data", "run_id", "is_active"}:
            continue
        if any(tok in low for tok in exclude_tokens):
            continue
        engine_cols.append(col)

    if engine_cols:
        engines = df[engine_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        consensus_scores = pd.to_numeric(df[score_col], errors="coerce").fillna(0)
        mask_zero = (engines == 0).all(axis=1) & (consensus_scores == 0)
        filtered = df.loc[~mask_zero].copy()
        dropped = int(mask_zero.sum())
        if dropped > 0:
            dbg(
                "INFO",
                "consensus.filter",
                f"method=engine_zero engine_cols={engine_cols} dropped={dropped} kept={len(filtered)}/{start}",
            )
        else:
            dbg(
                "DEBUG",
                "consensus.filter",
                f"method=engine_zero engine_cols={engine_cols} dropped=0 kept={len(filtered)}/{start}",
            )
        return filtered, dropped

    dbg("WARN", "consensus.filter", "method=engine_zero reason=no_numeric_candidates")
    return df, 0
