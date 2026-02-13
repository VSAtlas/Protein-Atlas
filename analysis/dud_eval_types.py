from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Set

@dataclass
class TargetSpec:
    target_key: str
    pdb_id: str
    variant: Optional[str]
    ph_tag: Optional[str]
    csv_path: Optional[Path]
    source: str  # "manifest" | "scan" | "legacy"

@dataclass
class TargetEvaluation:
    metrics: Optional[Any]
    ligand_basenames: Set[str]
    has_run_id_column: bool
    run_ids: Set[str]
    status_reason: str = "ok"

CSV_BASENAMES = (
    "dud_docking_score_long.csv",  # new DUD runs (preferred)
    "docking_score_long.csv",  # legacy name (fallback)
)

CONSENSUS_CSV_BASENAME = "consensus_docking_scores.csv"

RERANKED_SCORCH_BASENAME = "consensus_reranked_scorch.csv"
