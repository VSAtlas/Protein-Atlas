"""Collect formula metadata from the modules that implement statistics."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List

from analysis.reporting import fdr_methods
from analysis import statistics
from analysis.controls import matched_controls
from analysis.external import spd


REPORT_ONLY_FORMULAS = (
    {
        "title": "Bootstrap confidence interval",
        "formula": "CI = percentile(metric_bootstrap, [2.5%, 97.5%])",
        "notes": "Rows are resampled with replacement, the metric is recomputed, and the percentile interval is displayed when bootstrap outputs are available.",
        "source": "analysis.run_enrichment.run_enrichment",
    },
    {
        "title": "Ablation reranking",
        "formula": "delta_metric = metric_full_score - metric_ablated_score",
        "notes": "Each ablation removes one score component, recomputes the ranking, then recomputes enrichment/AUPRC metrics.",
        "source": "analysis.run_enrichment.run_ablation_suite",
    },
)


def _formula_rows(rows: Iterable[Dict[str, Any]]) -> List[Dict[str, str]]:
    return [
        {
            "title": str(row.get("title", "")).strip(),
            "formula": str(row.get("formula", "")).strip(),
            "notes": str(row.get("notes", "")).strip(),
            "source": str(row.get("source", "")).strip(),
        }
        for row in rows
    ]


def collect_statistical_formulas() -> List[Dict[str, str]]:
    """Return formula metadata exported by the statistical implementation modules."""

    formulas: List[Dict[str, str]] = []
    for module in (fdr_methods, statistics, matched_controls, spd):
        formulas.extend(_formula_rows(getattr(module, "STATISTICAL_FORMULAS", ())))
    formulas.extend(_formula_rows(REPORT_ONLY_FORMULAS))
    return formulas
