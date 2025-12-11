from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional

logger = logging.getLogger(__name__)

EASY_THRESHOLD = 0.80
MEDIUM_THRESHOLD = 0.60


@dataclass
class TargetDifficulty:
    pdb_id: str
    roc_auc: float
    difficulty: str  # "easy" | "medium" | "hard" | "degenerate"
    N: int
    n_actives: int


def bucket_from_auc(roc_auc: float) -> str:
    """
    Map ROC_AUC to a coarse difficulty bucket.
    """
    try:
        val = float(roc_auc)
    except Exception:
        return "degenerate"
    if math.isnan(val):
        return "degenerate"
    if val < MEDIUM_THRESHOLD:
        return "hard"
    if val < EASY_THRESHOLD:
        return "medium"
    return "easy"


def difficulty_from_metrics(pdb_id: str, metrics: Optional[Mapping[str, object]]) -> TargetDifficulty:
    """
    Translate a dud_eval TargetEvaluation.metrics Series (or mapping) into TargetDifficulty.
    """
    roc_auc = float("nan")
    N = 0
    n_actives = 0
    if metrics is not None:
        try:
            roc_auc = float(metrics.get("ROC_AUC", float("nan")))  # type: ignore[arg-type]
        except Exception:
            roc_auc = float("nan")
        try:
            N = int(metrics.get("N", 0))  # type: ignore[arg-type]
        except Exception:
            N = 0
        try:
            n_actives = int(metrics.get("n_actives", 0))  # type: ignore[arg-type]
        except Exception:
            n_actives = 0

    bucket = bucket_from_auc(roc_auc)
    return TargetDifficulty(
        pdb_id=pdb_id,
        roc_auc=roc_auc,
        difficulty=bucket,
        N=N,
        n_actives=n_actives,
    )


def evaluate_difficulty_for_target(
    pdb_id: str,
    csv_path: Path,
    analysis_root: Path,
    lig_col: str = "ligand",
    score_col: str = "score",
    bedroc_alpha: float = 20.0,
    logauc_lambda: float = 1e-3,
    run_id: Optional[str] = None,
) -> TargetDifficulty:
    """
    Use dud_eval.evaluate_target to compute ROC_AUC and map it to a difficulty bucket.
    """
    from dud_eval import evaluate_target  # late import to avoid heavy deps unless needed

    out_dir = analysis_root / pdb_id
    csv_path = Path(csv_path)

    evaluated = None
    try:
        evaluated = evaluate_target(
            pdb_id=pdb_id,
            csv_path=csv_path,
            out_dir=out_dir,
            lig_col_cli=lig_col,
            score_col_cli=score_col,
            bedroc_alpha=bedroc_alpha,
            logauc_lambda=logauc_lambda,
            run_id=run_id,
        )
    except Exception as exc:
        logger.warning(
            "[difficulty] pdb=%s csv=%s action=skip reason=%s",
            pdb_id,
            csv_path,
            exc,
        )
        return TargetDifficulty(
            pdb_id=pdb_id,
            roc_auc=float("nan"),
            difficulty="degenerate",
            N=0,
            n_actives=0,
        )

    metrics = evaluated.metrics if evaluated is not None else None
    if metrics is None:
        logger.warning(
            "[difficulty] pdb=%s csv=%s action=skip reason=no_metrics",
            pdb_id,
            csv_path,
        )
    td = difficulty_from_metrics(pdb_id, metrics)
    logger.info(
        "[difficulty] pdb=%s roc_auc=%.3f difficulty=%s N=%d n_actives=%d",
        td.pdb_id,
        td.roc_auc,
        td.difficulty,
        td.N,
        td.n_actives,
    )
    return td


def evaluate_difficulty_for_targets(
    targets: Dict[str, Path],
    analysis_root: Path,
    lig_col: str = "ligand",
    score_col: str = "score",
    bedroc_alpha: float = 20.0,
    logauc_lambda: float = 1e-3,
    run_id: Optional[str] = None,
) -> Dict[str, TargetDifficulty]:
    """
    Compute difficulty for a set of targets keyed by pdb_id.
    """
    results: Dict[str, TargetDifficulty] = {}
    for pdb_id, csv_path in targets.items():
        results[pdb_id] = evaluate_difficulty_for_target(
            pdb_id=pdb_id,
            csv_path=csv_path,
            analysis_root=analysis_root,
            lig_col=lig_col,
            score_col=score_col,
            bedroc_alpha=bedroc_alpha,
            logauc_lambda=logauc_lambda,
            run_id=run_id,
        )
    return results
