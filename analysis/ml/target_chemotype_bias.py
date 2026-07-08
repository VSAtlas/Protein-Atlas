from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any
import warnings

import pandas as pd
from pandas.errors import PerformanceWarning

from analysis.ml.labels import binary_label_series
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier_core import (
    _fit_design_matrix,
    _fit_model_with_optional_weights,
    _model,
    _model_probabilities,
    _transform_design_matrix,
)
from analysis.statistics import auroc, average_precision, brier_score, ranked_binary_metrics


DEFAULT_GROUP_COLS = [
    "target_id",
    "target_family",
    "protein_class",
    "ligand_chemotype",
    "scaffold_key",
    "chemical_cluster",
]

DEFAULT_PAIRWISE_GROUPS = [
    ("target_family", "ligand_chemotype"),
    ("target_family", "scaffold_key"),
    ("protein_class", "ligand_chemotype"),
    ("target_id", "ligand_chemotype"),
]

DEFAULT_SHORTCUT_FEATURE_SETS: dict[str, list[str]] = {
    "target_id_only": ["target_id"],
    "target_context": ["target_family", "protein_class"],
    "chemotype_only": ["ligand_chemotype"],
    "scaffold_only": ["scaffold_key"],
    "chemical_cluster_only": ["chemical_cluster"],
    "target_context_plus_chemotype": ["target_family", "protein_class", "ligand_chemotype"],
    "target_context_plus_scaffold": ["target_family", "protein_class", "scaffold_key"],
}


def _available(cols: list[str], df: pd.DataFrame) -> list[str]:
    return [col for col in cols if col in df.columns and df[col].notna().any()]


def _topk(labels: list[int], scores: list[float], *, k: int) -> dict[str, float | int]:
    if not labels:
        return {"precision_at_K": math.nan, "enrichment_at_K": math.nan, "K": int(k)}
    pairs = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    top_n = min(int(k), len(pairs))
    prevalence = sum(labels) / len(labels) if labels else 0.0
    precision = sum(label for _score, label in pairs[:top_n]) / top_n if top_n else 0.0
    return {
        "precision_at_K": float(precision),
        "enrichment_at_K": float(precision / prevalence) if prevalence else 0.0,
        "K": int(top_n),
    }


def _format_group_key(frame: pd.DataFrame, cols: list[str]) -> pd.Series:
    if not cols:
        return pd.Series(["all"] * len(frame), index=frame.index)
    return frame[cols].fillna("__missing__").astype(str).agg("|".join, axis=1)


def _group_prevalence(
    data: pd.DataFrame,
    *,
    label_col: str,
    group_cols: list[str],
    min_group_n: int,
    min_positive: int,
    min_negative: int,
    dominance_threshold: float,
) -> pd.DataFrame:
    available = _available(group_cols, data)
    if not available:
        return pd.DataFrame()
    global_prev = float(data[label_col].mean()) if len(data) else math.nan
    total_pos = int(data[label_col].sum())
    total_neg = int((1 - data[label_col]).sum())
    frames: list[pd.DataFrame] = []
    for group_col in available:
        work = data[[group_col, label_col]].copy()
        work[group_col] = work[group_col].fillna("__missing__").astype(str)
        grouped = (
            work.groupby(group_col, dropna=False)[label_col]
            .agg(n="count", positives="sum")
            .reset_index()
            .rename(columns={group_col: "group_value"})
        )
        grouped["group_col"] = group_col
        grouped["negatives"] = grouped["n"] - grouped["positives"]
        grouped["positive_rate"] = grouped["positives"] / grouped["n"]
        grouped["prevalence_lift"] = grouped["positive_rate"] / global_prev if global_prev else math.nan
        grouped["positive_share"] = grouped["positives"] / total_pos if total_pos else 0.0
        grouped["negative_share"] = grouped["negatives"] / total_neg if total_neg else 0.0
        grouped["is_all_positive"] = grouped["negatives"].eq(0)
        grouped["is_all_negative"] = grouped["positives"].eq(0)
        grouped["too_small_for_fold"] = (
            grouped["n"].lt(min_group_n)
            | grouped["positives"].lt(min_positive)
            | grouped["negatives"].lt(min_negative)
        )
        grouped["dominates_positives"] = grouped["positive_share"].ge(dominance_threshold)
        grouped["dominates_negatives"] = grouped["negative_share"].ge(dominance_threshold)
        grouped["risk_score"] = (
            grouped["is_all_positive"].astype(int)
            + grouped["is_all_negative"].astype(int)
            + grouped["too_small_for_fold"].astype(int)
            + grouped["dominates_positives"].astype(int)
            + grouped["dominates_negatives"].astype(int)
        )
        frames.append(grouped)
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["risk_score", "n"], ascending=[False, False], kind="stable")


def _pairwise_prevalence(
    data: pd.DataFrame,
    *,
    label_col: str,
    group_pairs: list[tuple[str, str]],
    min_group_n: int,
    min_positive: int,
    min_negative: int,
    dominance_threshold: float,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for left, right in group_pairs:
        available = _available([left, right], data)
        if len(available) != 2:
            continue
        work = data[[left, right, label_col]].copy()
        work["group_key"] = _format_group_key(work, [left, right])
        summary = _group_prevalence(
            work,
            label_col=label_col,
            group_cols=["group_key"],
            min_group_n=min_group_n,
            min_positive=min_positive,
            min_negative=min_negative,
            dominance_threshold=dominance_threshold,
        )
        if summary.empty:
            continue
        summary.insert(0, "left_group_col", left)
        summary.insert(1, "right_group_col", right)
        frames.append(summary)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _shortcut_row(
    data: pd.DataFrame,
    *,
    label_col: str,
    feature_set_name: str,
    features: list[str],
    split: str,
    model_type: str,
    seed: int,
    top_k: int,
    model_n_jobs: int,
) -> dict[str, Any]:
    available = _available(features, data)
    if not available:
        return {
            "feature_set": feature_set_name,
            "features_used": "",
            "split_method": split,
            "model_type": model_type,
            "status": "skipped_no_available_features",
        }
    try:
        train_idx, test_idx = make_split(data, split_mode=split, seed=seed)
        train = data.loc[train_idx].copy()
        test = data.loc[test_idx].copy()
        if train[label_col].nunique() < 2 or test[label_col].nunique() < 2:
            return {
                "feature_set": feature_set_name,
                "features_used": ";".join(available),
                "split_method": split,
                "model_type": model_type,
                "status": "insufficient_labels",
                "n_train": int(len(train)),
                "n_test": int(len(test)),
                "train_positives": int(train[label_col].sum()),
                "test_positives": int(test[label_col].sum()),
            }
        x_train, preprocessing = _fit_design_matrix(train, available)
        x_test = _transform_design_matrix(test, preprocessing)
        params = {"n_jobs": model_n_jobs, "thread_count": model_n_jobs} if model_n_jobs > 0 else None
        clf = _model(model_type, seed, class_weight="balanced", model_params=params)
        clf = _fit_model_with_optional_weights(clf, x_train, train[label_col])
        scores = _model_probabilities(clf, x_test)
        labels = pd.to_numeric(test[label_col], errors="coerce").astype(int).tolist()
        split_summary = split_overlap_summary(train, test, split)
    except Exception as exc:  # noqa: BLE001 - audit records failures for incomplete metadata.
        return {
            "feature_set": feature_set_name,
            "features_used": ";".join(available),
            "split_method": split,
            "model_type": model_type,
            "status": "failed",
            "error": str(exc),
        }
    row: dict[str, Any] = {
        "feature_set": feature_set_name,
        "features_used": ";".join(available),
        "split_method": split,
        "model_type": model_type,
        "status": "ok",
        "number_of_rows": int(len(data)),
        "number_of_positives": int(data[label_col].sum()),
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "train_positives": int(train[label_col].sum()),
        "test_positives": int(test[label_col].sum()),
        "test_positive_rate": float(pd.Series(labels).mean()) if labels else math.nan,
        "split_passes_holdout": bool(split_summary.get("passes_holdout", False)),
    }
    row.update(ranked_binary_metrics(scores, labels, [0.01, 0.05, 0.10]))
    row.update(_topk(labels, scores, k=top_k))
    row["Brier_score"] = float(brier_score(scores, labels))
    row["AUROC"] = float(auroc(scores, labels))
    row["PR_AUC"] = float(average_precision(scores, labels))
    return row


def audit_target_chemotype_bias(
    dataset_path: str | Path,
    *,
    label_col: str,
    out_dir: str | Path,
    group_cols: list[str] | None = None,
    splits: list[str] | None = None,
    model_type: str = "logistic_regression",
    seed: int = 42,
    top_k: int = 20,
    min_group_n: int = 10,
    min_positive: int = 5,
    min_negative: int = 5,
    dominance_threshold: float = 0.25,
    model_n_jobs: int = 4,
) -> dict[str, Any]:
    warnings.filterwarnings("ignore", category=PerformanceWarning)
    frame = pd.read_csv(dataset_path, low_memory=False)
    observed = binary_label_series(frame[label_col])
    data = frame.loc[observed.notna()].copy()
    data[label_col] = observed.loc[observed.notna()].astype(int)

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    group_cols = group_cols or DEFAULT_GROUP_COLS
    splits = splits or ["target_holdout", "target_family_holdout", "chemical_cluster_holdout"]

    prevalence = _group_prevalence(
        data,
        label_col=label_col,
        group_cols=group_cols,
        min_group_n=min_group_n,
        min_positive=min_positive,
        min_negative=min_negative,
        dominance_threshold=dominance_threshold,
    )
    pairwise = _pairwise_prevalence(
        data,
        label_col=label_col,
        group_pairs=DEFAULT_PAIRWISE_GROUPS,
        min_group_n=min_group_n,
        min_positive=min_positive,
        min_negative=min_negative,
        dominance_threshold=dominance_threshold,
    )

    shortcut_rows: list[dict[str, Any]] = []
    for split in splits:
        for name, features in DEFAULT_SHORTCUT_FEATURE_SETS.items():
            shortcut_rows.append(
                _shortcut_row(
                    data,
                    label_col=label_col,
                    feature_set_name=name,
                    features=features,
                    split=split,
                    model_type=model_type,
                    seed=seed,
                    top_k=top_k,
                    model_n_jobs=model_n_jobs,
                )
            )
    shortcut = pd.DataFrame(shortcut_rows)

    prevalence.to_csv(out / "target_chemotype_group_prevalence.csv", index=False)
    pairwise.to_csv(out / "target_chemotype_pairwise_prevalence.csv", index=False)
    shortcut.to_csv(out / "target_chemotype_shortcut_models.csv", index=False)
    high_risk = prevalence.loc[prevalence["risk_score"].gt(0)].copy() if not prevalence.empty else pd.DataFrame()
    high_risk.to_csv(out / "target_chemotype_high_risk_groups.csv", index=False)

    ok_shortcut = shortcut.loc[shortcut.get("status", pd.Series(dtype=str)).eq("ok")].copy()
    best_shortcut: dict[str, Any] | None = None
    if not ok_shortcut.empty:
        best = ok_shortcut.sort_values(["PR_AUC", "AUROC"], ascending=[False, False]).iloc[0]
        best_shortcut = {key: (value.item() if hasattr(value, "item") else value) for key, value in best.to_dict().items()}
    findings: list[str] = []
    if not high_risk.empty:
        findings.append(f"{len(high_risk)} target/chemotype groups have shortcut-risk flags.")
    if best_shortcut and float(best_shortcut.get("PR_AUC", 0.0) or 0.0) >= 0.10:
        findings.append(
            "Shortcut-only metadata predicts labels above a low PR-AUC threshold; "
            "treat target/chemotype metadata as potential assay-context bias."
        )
    manifest: dict[str, Any] = {
        "dataset_path": str(dataset_path),
        "label_col": label_col,
        "n_labelable": int(len(data)),
        "n_positive": int(data[label_col].sum()),
        "n_negative": int((1 - data[label_col]).sum()),
        "positive_rate": float(data[label_col].mean()) if len(data) else math.nan,
        "group_cols": _available(group_cols, data),
        "splits": splits,
        "shortcut_model_type": model_type,
        "best_shortcut_model": best_shortcut,
        "findings": findings,
        "outputs": {
            "group_prevalence": str(out / "target_chemotype_group_prevalence.csv"),
            "pairwise_prevalence": str(out / "target_chemotype_pairwise_prevalence.csv"),
            "shortcut_models": str(out / "target_chemotype_shortcut_models.csv"),
            "high_risk_groups": str(out / "target_chemotype_high_risk_groups.csv"),
        },
    }
    (out / "target_chemotype_bias_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
