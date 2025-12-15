from __future__ import annotations

import json
import logging
import math
import os
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


def get_or_compute_target_difficulty(
    cfg: Mapping[str, object],
    pdb_id: str,
    csv_path: Path,
    analysis_root: Path,
    run_id: Optional[str] = None,
) -> TargetDifficulty:
    """
    Convenience helper for docking-time difficulty gating.

    - Caches difficulty on disk under analysis_root / f"{pdb_id}.json".
    - If cached, loads it instead of recomputing.
    - Otherwise, calls evaluate_difficulty_for_target(...).
    - If csv_path missing or dud_eval fails, returns a degenerate difficulty.

    Env override:
      FORCE_TARGET_DIFFICULTY=<easy|medium|hard|degenerate> forces the bucket.
      Used for tests; no CLI surface area change.
    """
    force = (os.environ.get("FORCE_TARGET_DIFFICULTY") or "").strip().lower()
    if force in {"easy", "medium", "hard", "degenerate"}:
        return TargetDifficulty(
            pdb_id=pdb_id,
            roc_auc=float("nan"),
            difficulty=force,
            N=0,
            n_actives=0,
        )

    analysis_root = Path(analysis_root)
    analysis_root.mkdir(parents=True, exist_ok=True)
    cache_path = analysis_root / f"{pdb_id}.json"

    if cache_path.exists():
        try:
            data = json.loads(cache_path.read_text(encoding="utf-8"))
            return TargetDifficulty(
                pdb_id=str(data.get("pdb_id", pdb_id)),
                roc_auc=float(data.get("roc_auc", float("nan"))),
                difficulty=str(data.get("difficulty", "degenerate")),
                N=int(data.get("N", 0)),
                n_actives=int(data.get("n_actives", 0)),
            )
        except Exception:
            pass

    csv_path = Path(csv_path)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        td = TargetDifficulty(
            pdb_id=pdb_id,
            roc_auc=float("nan"),
            difficulty="degenerate",
            N=0,
            n_actives=0,
        )
    else:
        td = evaluate_difficulty_for_target(
            pdb_id=pdb_id,
            csv_path=csv_path,
            analysis_root=analysis_root,
            lig_col="ligand",
            score_col="score",
            bedroc_alpha=float(cfg.get("DUD_EVAL_BEDROC_ALPHA", 20.0)),  # type: ignore[arg-type]
            logauc_lambda=float(cfg.get("DUD_EVAL_LOGAUC_LAMBDA", 1e-3)),  # type: ignore[arg-type]
            run_id=run_id,
        )

    try:
        cache_path.write_text(
            json.dumps(
                {
                    "pdb_id": td.pdb_id,
                    "roc_auc": td.roc_auc,
                    "difficulty": td.difficulty,
                    "N": td.N,
                    "n_actives": td.n_actives,
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
    except Exception:
        pass

    return td
