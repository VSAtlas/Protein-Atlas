"""FDR utilities for Atlas report exports."""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple


@dataclass(frozen=True)
class ScoreFdrResult:
    row_metrics: Dict[int, Dict[str, float]]
    summary: Dict[str, float]


STATISTICAL_FORMULAS = (
    {
        "title": "Empirical decoy p-value",
        "formula": "p = (1 + #{decoy_score >= observed_score}) / (1 + n_decoys)",
        "notes": "Implemented by compute_score_fdr for higher-is-better scores; lower-is-better callers reverse the extreme-score inequality.",
        "source": "analysis.reporting.fdr_methods.compute_score_fdr",
    },
    {
        "title": "Benjamini-Hochberg q-value",
        "formula": "q_(i) = min_{j >= i} min(1, p_(j) * m / j)",
        "notes": "Implemented by _bh_q_values on sorted empirical p-values within each target/PDB condition.",
        "source": "analysis.reporting.fdr_methods._bh_q_values",
    },
    {
        "title": "Minimum resolvable BH q",
        "formula": "min_possible_q = min(1, m / (1 + n_decoys))",
        "notes": "Computed from the empirical p-value floor. If this exceeds 0.10, that target cannot produce q <= 0.10 under empirical-p/BH correction.",
        "source": "analysis.reporting.fdr_methods.compute_score_fdr",
    },
    {
        "title": "Target-decoy competition q",
        "formula": "FDR_k = (D_{score >= score_k} * n_tested / n_decoys) / k; q_k = min_{j >= k} min(1, FDR_j)",
        "notes": "The report also emits a finite-sample version using D_{score >= score_k} + 1.",
        "source": "analysis.reporting.fdr_methods.compute_score_fdr",
    },
    {
        "title": "SCORCH conditional FDR",
        "formula": "SCORCH conditional q = target-decoy competition q on rows with SCORCH_score_used",
        "notes": "This is a top-tranche conditional null because SCORCH is not run for every FDA/DUD row.",
        "source": "analysis.reporting.master_schema_export",
    },
)


def _bh_q_values(indexed_p_values: List[Tuple[int, float]]) -> Dict[int, float]:
    if not indexed_p_values:
        return {}
    ordered = sorted(indexed_p_values, key=lambda item: item[1])
    m = len(ordered)
    raw: List[Tuple[int, float]] = []
    for rank, (idx, p_value) in enumerate(ordered, start=1):
        raw.append((idx, p_value * m / rank))
    q_values: Dict[int, float] = {}
    running = 1.0
    for idx, q_value in reversed(raw):
        running = min(running, q_value, 1.0)
        q_values[idx] = running
    return q_values


def _monotone_rank_q(raw_values: List[Tuple[int, float]]) -> Dict[int, float]:
    q_values: Dict[int, float] = {}
    running = 1.0
    for idx, value in reversed(raw_values):
        running = min(running, value, 1.0)
        q_values[idx] = running
    return q_values


def compute_score_fdr(
    tested_scores: Iterable[Tuple[int, float]],
    decoy_scores: Iterable[float],
    *,
    higher_is_better: bool = True,
) -> ScoreFdrResult:
    """Compute empirical/BH and target-decoy competition FDR for one ranked list.

    The target-decoy estimate is reported both without and with a +1 finite-sample
    pseudocount. A decoy-size correction is applied because Atlas FDA and DUD
    retained score sets can be unequal after conditional rescoring.
    """

    finite_decoys = [
        float(score)
        for score in decoy_scores
        if score is not None and math.isfinite(float(score))
    ]
    finite_tested = [
        (idx, float(score))
        for idx, score in tested_scores
        if score is not None and math.isfinite(float(score))
    ]
    n_decoys = len(finite_decoys)
    n_tested = len(finite_tested)
    row_metrics: Dict[int, Dict[str, float]] = {}
    if n_decoys == 0 or n_tested == 0:
        return ScoreFdrResult(
            row_metrics={},
            summary={
                "n_decoys": float(n_decoys),
                "n_tested": float(n_tested),
                "best_bh_q": 1.0,
                "best_decoy_competition_q": 1.0,
                "best_decoy_competition_q_plus1": 1.0,
                "min_possible_bh_q": 1.0,
            },
        )

    decoys_sorted = sorted(finite_decoys)
    tested_ranked = sorted(finite_tested, key=lambda item: item[1], reverse=higher_is_better)
    p_values: List[Tuple[int, float]] = []
    competition_raw: List[Tuple[int, float]] = []
    competition_raw_plus1: List[Tuple[int, float]] = []
    decoy_size_factor = n_tested / max(1, n_decoys)

    for rank, (idx, score) in enumerate(tested_ranked, start=1):
        if higher_is_better:
            n_decoy_extreme = n_decoys - bisect.bisect_left(decoys_sorted, score)
        else:
            n_decoy_extreme = bisect.bisect_right(decoys_sorted, score)
        empirical_p = (1.0 + n_decoy_extreme) / (1.0 + n_decoys)
        p_values.append((idx, empirical_p))
        competition_raw.append((idx, min(1.0, (n_decoy_extreme * decoy_size_factor) / rank)))
        competition_raw_plus1.append(
            (idx, min(1.0, ((n_decoy_extreme + 1.0) * decoy_size_factor) / rank))
        )
        row_metrics[idx] = {
            "p_empirical": empirical_p,
            "n_decoy_extreme": float(n_decoy_extreme),
            "rank": float(rank),
        }

    bh_q = _bh_q_values(p_values)
    competition_q = _monotone_rank_q(competition_raw)
    competition_q_plus1 = _monotone_rank_q(competition_raw_plus1)
    for idx, metrics in row_metrics.items():
        metrics["q_bh"] = bh_q.get(idx, 1.0)
        metrics["decoy_competition_q"] = competition_q.get(idx, 1.0)
        metrics["decoy_competition_q_plus1"] = competition_q_plus1.get(idx, 1.0)

    min_possible_bh_q = min(1.0, (1.0 / (1.0 + n_decoys)) * n_tested)
    summary = {
        "n_decoys": float(n_decoys),
        "n_tested": float(n_tested),
        "best_bh_q": min((m.get("q_bh", 1.0) for m in row_metrics.values()), default=1.0),
        "best_decoy_competition_q": min(
            (m.get("decoy_competition_q", 1.0) for m in row_metrics.values()),
            default=1.0,
        ),
        "best_decoy_competition_q_plus1": min(
            (m.get("decoy_competition_q_plus1", 1.0) for m in row_metrics.values()),
            default=1.0,
        ),
        "min_possible_bh_q": min_possible_bh_q,
        "n_decoy_competition_hits_q05": float(
            sum(1 for m in row_metrics.values() if m.get("decoy_competition_q", 1.0) <= 0.05)
        ),
        "n_decoy_competition_hits_q10": float(
            sum(1 for m in row_metrics.values() if m.get("decoy_competition_q", 1.0) <= 0.10)
        ),
        "n_decoy_competition_hits_q05_plus1": float(
            sum(
                1
                for m in row_metrics.values()
                if m.get("decoy_competition_q_plus1", 1.0) <= 0.05
            )
        ),
        "n_decoy_competition_hits_q10_plus1": float(
            sum(
                1
                for m in row_metrics.values()
                if m.get("decoy_competition_q_plus1", 1.0) <= 0.10
            )
        ),
    }
    return ScoreFdrResult(row_metrics=row_metrics, summary=summary)
