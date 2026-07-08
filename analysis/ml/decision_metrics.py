from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics


def _binary(frame: pd.DataFrame, label_col: str) -> pd.Series:
    return pd.to_numeric(frame[label_col], errors="coerce").fillna(0).astype(int)


def decision_metrics(
    df: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    group_col: str = "target_id",
    top_n: int = 20,
    hit_n: int = 10,
) -> dict[str, float]:
    work = df[[col for col in [label_col, score_col, group_col, "is_control", "known_control"] if col in df.columns]].dropna(
        subset=[label_col, score_col]
    )
    if work.empty:
        return {
            "precision@20_per_target_mean": 0.0,
            "hits@10_per_target_mean": 0.0,
            "recall_of_known_controls": 0.0,
        }
    labels = _binary(work, label_col)
    work = work.assign(_label=labels)
    grouped = work.groupby(group_col, dropna=False) if group_col in work.columns else [("all", work)]
    precision_values: list[float] = []
    hits_values: list[float] = []
    for _key, group in grouped:
        ranked = group.sort_values(score_col, ascending=False)
        top_precision = ranked.head(top_n)
        top_hits = ranked.head(hit_n)
        if len(top_precision):
            precision_values.append(float(top_precision["_label"].mean()))
        if len(top_hits):
            hits_values.append(float(top_hits["_label"].sum()))
    control_cols = [col for col in ["is_control", "known_control"] if col in work.columns]
    recall_controls = 0.0
    if control_cols:
        control_mask = pd.Series(False, index=work.index)
        for col in control_cols:
            control_mask |= work[col].fillna(False).astype(str).str.lower().isin({"1", "true", "yes"})
        controls = work.loc[control_mask]
        if len(controls):
            ranked = work.sort_values(score_col, ascending=False)
            k = max(1, math.ceil(len(work) * 0.10))
            top_ids = set(ranked.head(k).index)
            recall_controls = len(top_ids & set(controls.index)) / len(controls)
    return {
        "precision@20_per_target_mean": float(pd.Series(precision_values).mean()) if precision_values else 0.0,
        "hits@10_per_target_mean": float(pd.Series(hits_values).mean()) if hits_values else 0.0,
        "recall_of_known_controls": float(recall_controls),
    }


def score_metric_row(
    df: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    method: str,
) -> dict[str, Any]:
    work = df.dropna(subset=[label_col, score_col]).copy()
    if work.empty or pd.to_numeric(work[label_col], errors="coerce").nunique() < 2:
        return {"method": method, "score_col": score_col, "n": int(len(work)), "status": "insufficient_labels"}
    labels = pd.to_numeric(work[label_col], errors="coerce").astype(int).tolist()
    scores = pd.to_numeric(work[score_col], errors="coerce").tolist()
    row: dict[str, Any] = {
        "method": method,
        "score_col": score_col,
        "n": int(len(work)),
        "positive_rate": float(pd.Series(labels).mean()),
        "status": "ok",
    }
    row.update(calibration_metrics(scores, labels, scores))
    row.update(decision_metrics(work.assign(**{score_col: scores}), label_col=label_col, score_col=score_col))
    return row


def write_decision_metrics(
    pred: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    out_path: str | Path,
    method: str = "trained_model",
) -> pd.DataFrame:
    table = pd.DataFrame([score_metric_row(pred, label_col=label_col, score_col=score_col, method=method)])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_path, index=False)
    return table


def random_baseline_rows(
    df: pd.DataFrame,
    *,
    label_col: str,
    n_repeats: int = 20,
    seed: int = 42,
) -> list[dict[str, Any]]:
    rng = random.Random(seed)
    rows: list[dict[str, Any]] = []
    work = df.dropna(subset=[label_col]).copy()
    for idx in range(n_repeats):
        score_col = f"_random_score_{idx}"
        work[score_col] = [rng.random() for _ in range(len(work))]
        row = score_metric_row(work, label_col=label_col, score_col=score_col, method="random")
        row["repeat"] = idx
        rows.append(row)
    return rows


def group_topk_recovery(
    df: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    group_cols: list[str] | None = None,
    top_ks: list[int] | None = None,
) -> pd.DataFrame:
    """Return top-K precision/enrichment within audit groups.

    These rows are decision metrics, not calibration claims. They are useful
    for sparse ADR/tissue/mechanism tasks where the practical endpoint is a
    short follow-up list inside a target, source, family, or scaffold stratum.
    """

    groups = group_cols or [
        "target_id",
        "target_family",
        "protein_class",
        "label_source",
        "source_family",
        "external_source_family",
        "external_upstream_source",
        "scaffold_key",
        "chemical_cluster",
        "ligand_chemotype",
        "assay_type",
        "endpoint_type",
        "activity_type",
        "applicability_domain",
        "chemical_fingerprint_ad",
        "target_family_ad",
    ]
    ks = top_ks or [10, 20, 50]
    if label_col not in df.columns or score_col not in df.columns:
        return pd.DataFrame()
    labels = pd.to_numeric(df[label_col], errors="coerce")
    scores = pd.to_numeric(df[score_col], errors="coerce")
    valid = labels.isin([0, 1]) & scores.notna()
    work = df.loc[valid].copy()
    if work.empty:
        return pd.DataFrame()
    work["_label"] = labels.loc[valid].astype(int)
    work["_score"] = scores.loc[valid].astype(float)
    rows: list[dict[str, Any]] = []
    for group_col in groups:
        if group_col not in work.columns:
            continue
        for group_value, group in work.groupby(group_col, dropna=False):
            n = int(len(group))
            positives = int(group["_label"].sum())
            negatives = int(n - positives)
            prevalence = positives / n if n else 0.0
            ranked = group.sort_values("_score", ascending=False)
            for top_k in ks:
                k = min(int(top_k), n)
                if k <= 0:
                    continue
                top = ranked.head(k)
                precision = float(top["_label"].mean())
                hit_count = int(top["_label"].sum())
                rows.append(
                    {
                        "group_col": group_col,
                        "group_value": group_value,
                        "top_k": int(top_k),
                        "effective_k": k,
                        "n": n,
                        "n_positive": positives,
                        "n_negative": negatives,
                        "prevalence": prevalence,
                        "hit_count_at_k": hit_count,
                        "precision_at_k": precision,
                        "enrichment_at_k": precision / prevalence if prevalence > 0 else None,
                    }
                )
    return pd.DataFrame(rows)


def write_group_topk_recovery(
    pred: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    out_path: str | Path,
    group_cols: list[str] | None = None,
    top_ks: list[int] | None = None,
) -> pd.DataFrame:
    table = group_topk_recovery(
        pred,
        label_col=label_col,
        score_col=score_col,
        group_cols=group_cols,
        top_ks=top_ks,
    )
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_path, index=False)
    return table
