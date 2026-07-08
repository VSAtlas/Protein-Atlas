from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

try:
    from typing import TypeAlias
except ImportError:  # Python 3.7 in scorch-env
    TypeAlias = object


@dataclass(frozen=True)
class StageSpec:
    source: str  # vina, gnina, ledock, dock6
    stage_dir: str
    output_name: str


@dataclass(frozen=True)
class MaterializeResult:
    input_dir: Optional[Path]
    materialized_count: int
    failed_count: int
    failed_examples: Tuple[str, ...]


@dataclass(frozen=True)
class AnnotateResult:
    ok: bool
    degraded: bool
    reason: str = "ok"


@dataclass(frozen=True)
class SelectionResult:
    allowed_bases: Set[str]
    n_pool: int
    k: int
    controls_total: int
    selected_stage_by_base: Dict[str, str]
    selected_score_by_base: Dict[str, float]


@dataclass(frozen=True)
class PoseSelectionResult:
    ligands: List[Path]
    stage_counts: Dict[int, int]
    total_candidates: int
    available_bases: Set[str]
    rescored_stage_by_base: Dict[str, str]
    stage_fallback_reason_by_base: Dict[str, str]


ScorchTask: TypeAlias = Tuple[
    StageSpec,
    Tuple[str, str, str],
    Path,
    Set[str],
    Set[str],
    str,
    Optional[Sequence[str]],
    Optional[Path],
    Dict[str, str],
    Dict[str, float],
    Optional[str],
]
