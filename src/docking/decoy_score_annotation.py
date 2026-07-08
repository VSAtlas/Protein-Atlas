from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd  # type: ignore[import-untyped]
from path_router import make_paths  # type: ignore[import-not-found]

def annotate_fda_long_csv_with_z_scores_vs_decoys(
    cfg: Dict[str, Any],
    pdb_id: str,
    ph_label: Optional[str] = None,
    *,
    csv_prefix: str = "",
    decoy_csv_prefix: str = "dud_",
    logger=None,
) -> Optional[str]:
    """
    Post-processing helper used by docking.py when TEST_MODE_ENABLE includes DUD
    and we are running the FDA subrun.

    It reads <decoy_prefix>docking_score_long.csv to compute mean/std of the
    DUD-role background scores, then annotates <csv_prefix>docking_score_long.csv
    with a canonical z_vs_decoys column and legacy t_vs_decoys alias. Filename
    active/decoy labels are honored when present. If the DUD-role CSV has no
    labels at all, it is treated as an unlabeled decoy background.

    Returns the path to the updated FDA long CSV, or None if skipped.
    """
    try:
        from analysis.dud_eval import (
            compute_decoy_stats_from_long_csv,
            guess_ligfile_col,
            guess_score_col,
        )
    except Exception as e:
        if logger:
            logger.warning("[z-score.skip] pdb_id=%s reason=import_error %s", pdb_id, e)
        return None

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    ph_token = (ph_label or "").strip() or None
    variant_root = Path(paths.docked_variant_root(var, ph_token))

    dud_csv = variant_root / f"{decoy_csv_prefix}docking_score_long.csv"
    fda_csv = variant_root / f"{csv_prefix}docking_score_long.csv"

    if not dud_csv.exists() or not fda_csv.exists():
        if logger:
            logger.info(
                "[z-score.skip] pdb_id=%s ph=%s reason=missing_csv dud=%s fda=%s",
                pdb_id,
                ph_label or "base",
                str(dud_csv),
                str(fda_csv),
            )
        return None

    mu, sigma, n_decoys = compute_decoy_stats_from_long_csv(
        dud_csv, assume_unlabeled_decoys=True
    )
    if (
        not n_decoys
        or not math.isfinite(mu)
        or not math.isfinite(sigma)
        or sigma == 0.0
    ):
        if logger:
            logger.info(
                "[z-score.skip] pdb_id=%s ph=%s reason=degenerate_stats n=%s mu=%s sigma=%s",
                pdb_id,
                ph_label or "base",
                n_decoys,
                mu,
                sigma,
            )
        return None

    df = pd.read_csv(fda_csv)

    lig_col = guess_ligfile_col(df, None)
    score_col = guess_score_col(df, None)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    # Best score per ligand
    best = df.groupby(lig_col, as_index=False).agg(best_score=(score_col, "min"))
    best["z_vs_decoys"] = (mu - best["best_score"]) / sigma
    z_map = dict(zip(best[lig_col], best["z_vs_decoys"]))

    df["z_vs_decoys"] = df[lig_col].map(z_map)
    # Legacy compatibility alias.
    df["t_vs_decoys"] = df["z_vs_decoys"]

    df.to_csv(fda_csv, index=False)

    if logger:
        logger.info(
            "[z-score.ok] pdb_id=%s ph=%s n_decoys=%s mean=%.3f std=%.3f out=%s",
            pdb_id,
            ph_label or "base",
            n_decoys,
            mu,
            sigma,
            str(fda_csv),
        )

    return str(fda_csv)


def annotate_fda_long_csv_with_t_scores_vs_decoys(
    cfg: Dict[str, Any],
    pdb_id: str,
    ph_label: Optional[str] = None,
    *,
    csv_prefix: str = "",
    decoy_csv_prefix: str = "dud_",
    logger=None,
) -> Optional[str]:
    """Compatibility wrapper for legacy call sites."""
    return annotate_fda_long_csv_with_z_scores_vs_decoys(
        cfg,
        pdb_id,
        ph_label=ph_label,
        csv_prefix=csv_prefix,
        decoy_csv_prefix=decoy_csv_prefix,
        logger=logger,
    )


# -------------------------
# Common path builder per protein
# -------------------------
# >>> BUILD_PATHS SHIM START
