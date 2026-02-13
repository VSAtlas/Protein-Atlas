from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd

from analysis.dud_eval_schema import guess_ligfile_col, guess_score_col

TOKEN_RE = re.compile(r"(?<![A-Za-z0-9])(active|decoy)s?(?![A-Za-z0-9])", re.IGNORECASE)

POSE_TAIL_RE = re.compile(
    r"(?:_pose\d+|_mode\d+|_conf\d+|_rank\d+|_cluster\d+|_p\d+)$", re.IGNORECASE
)

MICROSTATE_TAIL_RE = re.compile(r"(?:__ms|_ms)_[A-Za-z0-9]+$", re.IGNORECASE)

EXT_RE = re.compile(r"\.(pdbqt|sdf|mol2)(?:\.gz)?$", re.IGNORECASE)

def parse_name_and_label(path_str: str) -> Tuple[str, Optional[int]]:
    """Return (ligand_root_id, is_active) from a path or filename.
    is_active: 1 for active, 0 for decoy, None if not found."""
    base = os.path.basename(str(path_str))
    stem = EXT_RE.sub("", base)  # drop .pdbqt/.sdf/.mol2(.gz)
    stem = POSE_TAIL_RE.sub("", stem)  # drop trailing _pose1 etc.
    stem = MICROSTATE_TAIL_RE.sub("", stem)  # drop trailing microstate token
    lig_id = stem  # normalized ligand identifier (pose+microstate collapsed)
    m = TOKEN_RE.search(base)
    if not m:
        return lig_id, None
    tok = m.group(1).lower()
    return lig_id, 1 if tok == "active" else 0

def compute_decoy_stats_from_long_csv(
    csv_path: Path | str,
    lig_col_override: Optional[str] = None,
    score_col_override: Optional[str] = None,
) -> tuple[float, float, int]:
    """
    Compute (mean, std, n) of best decoy scores from a docking_score_long.csv-style file.

    - csv_path: path to the CSV file (per-POSE rows).
    - lig_col_override / score_col_override: optional explicit column names.

    Decoy/active labels are derived from the ligand filename using parse_name_and_label,
    and best scores are defined as the minimum score per ligand.

    Returns (mean_decoy_score, std_decoy_score, n_decoys).
    If no decoys are found or the file is missing, returns (nan, nan, 0).
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return float("nan"), float("nan"), 0

    df = pd.read_csv(csv_path)

    lig_col = guess_ligfile_col(df, lig_col_override)
    score_col = guess_score_col(df, score_col_override)

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    lig_series = df[lig_col].astype(str)
    parsed = lig_series.apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])

    best = df.groupby("lig_id", as_index=False).agg(
        best_score=(score_col, "min"),
        is_active=("is_active", "max"),
    )

    decoys = best[best["is_active"] == 0]
    n_decoys = int(len(decoys))
    if n_decoys == 0:
        return float("nan"), float("nan"), 0

    mu = float(decoys["best_score"].mean())
    sigma = float(decoys["best_score"].std(ddof=1))
    return mu, sigma, n_decoys

def _make_placeholder_row(
    pdb_id: str,
    variant: Optional[str],
    ph_tag: Optional[str],
    run_id: Optional[str],
    bedroc_alpha: float,
    status_reason: str,
) -> pd.Series:
    base_row = {
        "run_id": str(run_id) if run_id else "(none)",
        "pdb_id": pdb_id,
        "variant": variant or "",
        "pH": ph_tag or "",
        "N": 0,
        "n_actives": 0,
        "actives_fraction": float("nan"),
        "ROC_AUC": float("nan"),
        "PR_AUC": float("nan"),
        "logAUC": float("nan"),
        "logAUC_adj": float("nan"),
        f"BEDROC_alpha_{bedroc_alpha:g}": float("nan"),
        "EF@1%": float("nan"),
        "EF@2%": float("nan"),
        "EF@5%": float("nan"),
        "EF@10%": float("nan"),
        "status_reason": status_reason,
    }
    return pd.Series(base_row)

def _collapse_best_scores(
    df: pd.DataFrame, score_col: str, *, best_is_min: bool = True
) -> tuple[pd.DataFrame, int]:
    """
    Collapse to one row per lig_id using either min or max score.

    Returns (best_df, drop_nan_best) where drop_nan_best counts ligands removed
    due to NaN best_score.
    """
    agg_fn = "min" if best_is_min else "max"
    best = df.groupby("lig_id", as_index=False).agg(
        best_score=(score_col, agg_fn),
        is_active=("is_active", "max"),  # any active -> active
    )
    before_best = len(best)
    best = best.dropna(subset=["best_score"])
    drop_nan_best = before_best - len(best)
    return best, drop_nan_best
