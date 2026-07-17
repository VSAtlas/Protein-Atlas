from __future__ import annotations

from collections.abc import Iterable, Mapping
import json
import math
import random
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.train_classifier_core import (
    _fit_design_matrix,
    _fit_model_with_optional_weights,
    _model,
    _model_probabilities,
    _sample_weights,
    _transform_design_matrix,
)
from analysis.statistics import permutation_p_value


def _present_features(
    df: pd.DataFrame,
    *,
    label_col: str,
    feature_set: str,
    exclude_features: Iterable[str] | None = None,
    strict_feature_set: bool = False,
) -> tuple[list[str], list[str], list[str]]:
    excluded = effective_exclude_features(
        label_col,
        list(exclude_features) if exclude_features is not None else None,
    )
    requested = [feature for feature in get_feature_set(feature_set) if feature not in excluded]
    missing = [feature for feature in requested if feature not in df.columns]
    present = [feature for feature in requested if feature in df.columns]
    all_missing = [feature for feature in present if df[feature].notna().sum() == 0]
    if strict_feature_set and (missing or all_missing):
        parts = []
        if missing:
            parts.append(f"missing columns: {', '.join(missing[:25])}")
        if all_missing:
            parts.append(f"all-missing columns: {', '.join(all_missing[:25])}")
        raise ValueError(f"strict feature-set check failed for {feature_set!r}: {'; '.join(parts)}")
    features = [feature for feature in present if feature not in all_missing]
    assert_no_leakage(features)
    return features, missing, all_missing


def _binary_metrics(labels: pd.Series, scores: pd.Series) -> dict[str, float]:
    from sklearn.metrics import average_precision_score, roc_auc_score

    y = pd.to_numeric(labels, errors="coerce").astype(int)
    probability = pd.to_numeric(scores, errors="coerce").astype(float)
    calibration = pd.DataFrame({"label": y, "score": probability.clip(0.0, 1.0)})
    calibration["bin"] = pd.cut(
        calibration["score"],
        bins=[idx / 10 for idx in range(11)],
        include_lowest=True,
    )
    ece = 0.0
    for _, group in calibration.groupby("bin", observed=False):
        if group.empty:
            continue
        ece += (len(group) / len(calibration)) * abs(float(group["score"].mean() - group["label"].mean()))
    brier = float(((probability - y) ** 2).mean())
    prevalence = float(y.mean())
    brier_baseline = prevalence * (1.0 - prevalence)
    return {
        "AUROC": float(roc_auc_score(y, probability)),
        "AUPRC": float(average_precision_score(y, probability)),
        "Brier": brier,
        "Brier_prevalence_baseline": brier_baseline,
        "Brier_skill_score": (
            1.0 - brier / brier_baseline if brier_baseline > 0 else float("nan")
        ),
        "ECE": float(ece),
    }


def _topk_metrics(labels: pd.Series, scores: pd.Series, *, k: int) -> dict[str, float | int]:
    work = pd.DataFrame({"label": labels.astype(int), "score": scores.astype(float)}).dropna()
    if work.empty:
        return {"precision_at_K": math.nan, "enrichment_at_K": math.nan, "K": int(k)}
    top_n = min(max(1, int(k)), len(work))
    ranked = work.sort_values("score", ascending=False)
    precision = float(ranked.head(top_n)["label"].mean())
    prevalence = float(work["label"].mean())
    return {
        "precision_at_K": precision,
        "enrichment_at_K": precision / prevalence if prevalence else 0.0,
        "K": int(top_n),
    }


def _metric_row(
    pred: pd.DataFrame,
    *,
    label_col: str,
    score_col: str = "ml_prediction_score",
    top_k: int = 20,
) -> dict[str, Any]:
    valid = pred[[label_col, score_col]].copy()
    valid[label_col] = pd.to_numeric(valid[label_col], errors="coerce")
    valid[score_col] = pd.to_numeric(valid[score_col], errors="coerce")
    valid = valid.dropna()
    if valid.empty or valid[label_col].nunique() < 2:
        return {
            "n_test": int(len(valid)),
            "n_test_positive": int(valid[label_col].eq(1).sum()),
            "n_test_negative": int(valid[label_col].eq(0).sum()),
            "AUROC": math.nan,
            "AUPRC": math.nan,
            "Brier": math.nan,
            "Brier_prevalence_baseline": math.nan,
            "Brier_skill_score": math.nan,
            "ECE": math.nan,
            **_topk_metrics(valid[label_col], valid[score_col], k=top_k),
        }
    metrics = _binary_metrics(valid[label_col], valid[score_col])
    return {
        "n_test": int(len(valid)),
        "n_test_positive": int(valid[label_col].eq(1).sum()),
        "n_test_negative": int(valid[label_col].eq(0).sum()),
        "AUROC": metrics.get("AUROC"),
        "AUPRC": metrics.get("AUPRC"),
        "Brier": metrics.get("Brier"),
        "Brier_prevalence_baseline": metrics.get("Brier_prevalence_baseline"),
        "Brier_skill_score": metrics.get("Brier_skill_score"),
        "ECE": metrics.get("ECE"),
        **_topk_metrics(valid[label_col], valid[score_col], k=top_k),
    }


def _bootstrap_ci(
    pred: pd.DataFrame,
    *,
    label_col: str,
    score_col: str = "ml_prediction_score",
    n_bootstraps: int = 200,
    seed: int = 42,
    top_k: int = 20,
    cluster_col: str | None = None,
) -> pd.DataFrame:
    if n_bootstraps <= 0 or pred.empty:
        return pd.DataFrame(columns=["metric", "ci_low", "ci_high", "n_bootstraps"])
    rng = random.Random(seed)
    values: dict[str, list[float]] = {}
    cluster_frames: dict[str, pd.DataFrame] = {}
    if cluster_col and cluster_col in pred.columns:
        cluster_values = pred[cluster_col].fillna("missing").astype(str)
        cluster_frames = {
            str(value): pred.loc[index].copy()
            for value, index in cluster_values.groupby(cluster_values).groups.items()
        }
    clusters = list(cluster_frames)
    n = len(pred)
    for _idx in range(n_bootstraps):
        if clusters:
            sampled_clusters = [clusters[rng.randrange(len(clusters))] for _ in clusters]
            sample = pd.concat(
                [cluster_frames[value] for value in sampled_clusters],
                ignore_index=True,
            )
        else:
            idx = [rng.randrange(n) for _ in range(n)]
            sample = pred.iloc[idx]
        row = _metric_row(sample, label_col=label_col, score_col=score_col, top_k=top_k)
        if row.get("n_test_positive", 0) <= 0 or row.get("n_test_negative", 0) <= 0:
            continue
        for metric in ["AUROC", "AUPRC", "Brier", "ECE", "precision_at_K", "enrichment_at_K"]:
            value = row.get(metric)
            if value is not None and pd.notna(value):
                values.setdefault(metric, []).append(float(value))
    rows = []
    for metric, metric_values in values.items():
        series = pd.Series(metric_values)
        rows.append(
            {
                "metric": metric,
                "ci_low": float(series.quantile(0.025)),
                "ci_high": float(series.quantile(0.975)),
                "n_bootstraps": int(len(series)),
            }
        )
    return pd.DataFrame(rows)


def _permutation_rows(
    pred: pd.DataFrame,
    *,
    label_col: str,
    score_col: str = "ml_prediction_score",
    n_permutations: int = 1000,
    seed: int = 42,
) -> list[dict[str, Any]]:
    if n_permutations <= 0 or pred.empty:
        return []
    labels = pd.to_numeric(pred[label_col], errors="coerce")
    scores = pd.to_numeric(pred[score_col], errors="coerce")
    valid = labels.notna() & scores.notna()
    labels = labels.loc[valid].astype(int)
    scores = scores.loc[valid].astype(float)
    if labels.nunique() < 2:
        return []
    rows = []
    for metric in ["AUPRC", "AUROC", "EF@1%", "EF@5%", "EF@10%"]:
        rows.append(
            {
                "metric": metric,
                "permutation_p": permutation_p_value(
                    scores.tolist(),
                    labels.tolist(),
                    metric,
                    n_permutations=n_permutations,
                    seed=seed,
                ),
                "n_permutations": int(n_permutations),
            }
        )
    return rows


def _fit_predict_fold(
    data: pd.DataFrame,
    *,
    train_idx: pd.Index,
    test_idx: pd.Index,
    features: list[str],
    label_col: str,
    model_type: str,
    seed: int,
    class_weight: str | None,
    pu_mode: str,
    model_params: Mapping[str, object] | None = None,
    compute_applicability_domain: bool = False,
) -> pd.DataFrame:
    train = data.loc[train_idx].copy()
    test = data.loc[test_idx].copy()
    x_train, preprocessing = _fit_design_matrix(train, features)
    x_test = _transform_design_matrix(test, preprocessing)
    model = _model(model_type, seed, class_weight=class_weight, model_params=model_params)
    weights = _sample_weights(train, label_col, pu_mode=pu_mode)
    _fit_model_with_optional_weights(model, x_train, train[label_col], weights)
    probs = _model_probabilities(model, x_test)
    meta_cols = [
        "drug_id",
        "target_id",
        "pdb_id",
        "target_family",
        "protein_class",
        "label_source",
        "source_family",
        "upstream_source",
        "source_objective",
        "assay_type",
        "endpoint_type",
        "spd_activity_relation",
        "scaffold_key",
        "chemical_cluster",
        "butina_cluster",
        "ligand_chemotype",
        "mechanism_panel",
        "mechanism_panel_match_basis",
        "canonical_pair_key",
        label_col,
    ]
    pred = test[[col for col in meta_cols if col in test.columns]].copy()
    pred["ml_prediction_score"] = [float(value) for value in probs]
    pred["_source_index"] = test.index.astype(str)
    if compute_applicability_domain and not x_train.empty and not x_test.empty:
        from sklearn.neighbors import NearestNeighbors

        means = x_train.mean(axis=0)
        stds = x_train.std(axis=0).replace(0, 1.0).fillna(1.0)
        train_scaled = ((x_train - means) / stds).fillna(0.0)
        test_scaled = ((x_test - means) / stds).fillna(0.0)
        train_neighbors = 2 if len(train_scaled) > 1 else 1
        train_nn = NearestNeighbors(n_neighbors=train_neighbors)
        train_nn.fit(train_scaled)
        train_distance = train_nn.kneighbors(train_scaled, return_distance=True)[0][:, -1]
        threshold = float(pd.Series(train_distance).quantile(0.95))
        test_nn = NearestNeighbors(n_neighbors=1)
        test_nn.fit(train_scaled)
        test_distance = test_nn.kneighbors(test_scaled, return_distance=True)[0][:, 0]
        denominator = threshold if threshold > 0 else 1.0
        pred["applicability_distance"] = test_distance
        pred["applicability_threshold_train_q95"] = threshold
        pred["applicability_distance_ratio"] = test_distance / denominator
        pred["applicability_domain"] = [
            "in_domain" if float(distance) <= threshold else "out_of_domain"
            for distance in test_distance
        ]
    return pred


def _split_valid(
    data: pd.DataFrame,
    train_idx: pd.Index,
    test_idx: pd.Index,
    *,
    label_col: str,
    min_test_positives: int,
    min_test_negatives: int,
) -> tuple[bool, dict[str, Any]]:
    train_labels = data.loc[train_idx, label_col]
    test_labels = data.loc[test_idx, label_col]
    summary = {
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "n_train_positive": int(train_labels.eq(1).sum()),
        "n_train_negative": int(train_labels.eq(0).sum()),
        "n_test_positive": int(test_labels.eq(1).sum()),
        "n_test_negative": int(test_labels.eq(0).sum()),
    }
    ok = (
        summary["n_train_positive"] > 0
        and summary["n_train_negative"] > 0
        and summary["n_test_positive"] >= min_test_positives
        and summary["n_test_negative"] >= min_test_negatives
    )
    return ok, summary


def _stratified_group_folds(
    data: pd.DataFrame,
    *,
    label_col: str,
    group_col: str,
    n_splits: int,
    seed: int,
) -> list[tuple[pd.Index, pd.Index]]:
    groups = data[group_col].fillna("missing").astype(str)
    n_groups = groups.nunique()
    if n_groups < 2:
        return []
    split_count = min(max(2, int(n_splits)), int(n_groups))
    try:
        from sklearn.model_selection import StratifiedGroupKFold

        splitter = StratifiedGroupKFold(n_splits=split_count, shuffle=True, random_state=seed)
        return [
            (data.index[train_pos], data.index[test_pos])
            for train_pos, test_pos in splitter.split(data, data[label_col].astype(int), groups)
        ]
    except Exception:
        rng = random.Random(seed)
        values = list(groups.drop_duplicates())
        rng.shuffle(values)
        chunks = [values[idx::split_count] for idx in range(split_count)]
        folds = []
        for chunk in chunks:
            test_mask = groups.isin(chunk)
            folds.append((data.index[~test_mask], data.index[test_mask]))
        return folds


def _leave_one_group_folds(data: pd.DataFrame, *, group_col: str) -> list[tuple[str, pd.Index, pd.Index]]:
    values = data[group_col].fillna("missing").astype(str)
    out = []
    for value in sorted(values.unique()):
        test_mask = values.eq(value)
        out.append((value, data.index[~test_mask], data.index[test_mask]))
    return out



def _stratified_target_holdout_by_family_folds(
    data: pd.DataFrame,
    *,
    target_col: str = "target_id",
    family_col: str = "target_family",
    test_fraction: float = 0.2,
    repeats: int = 10,
    seed: int = 42,
) -> list[tuple[str, int, pd.Index, pd.Index]]:
    if target_col not in data.columns or family_col not in data.columns:
        return []
    families = data[family_col].fillna("missing").astype(str)
    targets = data[target_col].fillna("missing").astype(str)
    folds: list[tuple[str, int, pd.Index, pd.Index]] = []
    for repeat in range(int(repeats)):
        rng = random.Random(seed + repeat)
        heldout_targets: set[str] = set()
        skipped_families: list[str] = []
        for family in sorted(families.unique()):
            family_targets = sorted(targets.loc[families.eq(family)].unique().tolist())
            if len(family_targets) < 2:
                skipped_families.append(family)
                continue
            rng.shuffle(family_targets)
            n_holdout = max(1, int(round(len(family_targets) * float(test_fraction))))
            n_holdout = min(n_holdout, len(family_targets) - 1)
            heldout_targets.update(family_targets[:n_holdout])
        if not heldout_targets:
            continue
        test_mask = targets.isin(heldout_targets)
        heldout = ";".join(sorted(heldout_targets)[:50])
        if len(heldout_targets) > 50:
            heldout = f"{heldout};...(+{len(heldout_targets) - 50} more)"
        if skipped_families:
            heldout = f"{heldout} | single_target_families_kept_in_train={','.join(skipped_families)}"
        folds.append((heldout, repeat, data.index[~test_mask], data.index[test_mask]))
    return folds


def _summarize_fold_metrics(fold_metrics: pd.DataFrame) -> pd.DataFrame:
    if fold_metrics.empty:
        return pd.DataFrame()
    metrics = [
        "AUROC",
        "AUPRC",
        "Brier",
        "Brier_prevalence_baseline",
        "Brier_skill_score",
        "ECE",
        "precision_at_K",
        "enrichment_at_K",
    ]
    rows = []
    group_cols = ["model_type", "evaluation", "group_col"]
    for key, group in fold_metrics.groupby(group_cols, dropna=False):
        base = {col: value for col, value in zip(group_cols, key, strict=False)}
        base.update(
            {
                "n_folds": int(len(group)),
                "median_n_test": float(group["n_test"].median()),
                "median_n_test_positive": float(group["n_test_positive"].median()),
                "median_n_test_negative": float(group["n_test_negative"].median()),
            }
        )
        for metric in metrics:
            vals = pd.to_numeric(group[metric], errors="coerce").dropna()
            base[f"{metric}_median"] = float(vals.median()) if not vals.empty else math.nan
            base[f"{metric}_q025"] = float(vals.quantile(0.025)) if not vals.empty else math.nan
            base[f"{metric}_q975"] = float(vals.quantile(0.975)) if not vals.empty else math.nan
        rows.append(base)
    return pd.DataFrame(rows)


def _mean_oof_predictions(predictions: pd.DataFrame, *, label_col: str) -> pd.DataFrame:
    if predictions.empty:
        return predictions
    meta_cols = [
        col
        for col in predictions.columns
        if col
        not in {
            "ml_prediction_score",
            "repeat",
            "fold",
            "evaluation",
            "group_col",
            "heldout_group",
        }
    ]
    agg: dict[str, Any] = {"ml_prediction_score": "mean"}
    mean_cols = {
        "applicability_distance",
        "applicability_threshold_train_q95",
        "applicability_distance_ratio",
    }
    for col in meta_cols:
        agg[col] = "mean" if col in mean_cols else "first"
    out = predictions.groupby("_source_index", dropna=False).agg(agg).reset_index(drop=True)
    if "applicability_distance_ratio" in out.columns:
        out["applicability_domain"] = out["applicability_distance_ratio"].map(
            lambda value: "in_domain" if float(value) <= 1.0 else "out_of_domain"
        )
    if label_col in out.columns:
        out[label_col] = pd.to_numeric(out[label_col], errors="coerce")
    return out


def _write_group_feasibility(
    data: pd.DataFrame,
    *,
    label_col: str,
    group_cols: Iterable[str],
    out_path: Path,
) -> pd.DataFrame:
    rows = []
    for group_col in group_cols:
        if group_col not in data.columns:
            rows.append({"group_col": group_col, "group_value": "missing_column", "n": 0, "n_positive": 0, "n_negative": 0})
            continue
        for value, group in data.groupby(group_col, dropna=False):
            rows.append(
                {
                    "group_col": group_col,
                    "group_value": value,
                    "n": int(len(group)),
                    "n_positive": int(group[label_col].eq(1).sum()),
                    "n_negative": int(group[label_col].eq(0).sum()),
                    "positive_rate": float(group[label_col].mean()) if len(group) else math.nan,
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(out_path / "group_feasibility.csv", index=False)
    return frame


def run_grouped_cv_stability(
    dataset_path: str | Path,
    *,
    label_col: str,
    feature_set: str,
    out_dir: str | Path,
    model_types: list[str] | None = None,
    group_cols: list[str] | None = None,
    leave_one_group_cols: list[str] | None = None,
    n_splits: int = 5,
    repeats: int = 20,
    min_test_positives: int = 10,
    min_test_negatives: int = 10,
    top_k: int = 20,
    n_bootstraps: int = 200,
    n_permutations: int = 1000,
    seed: int = 42,
    class_weight: str | None = "balanced",
    pu_mode: str = "standard_binary",
    exclude_features: list[str] | None = None,
    strict_feature_set: bool = False,
    model_params: Mapping[str, object] | None = None,
    compute_applicability_domain: bool = False,
    run_stratified_target_holdout_by_family: bool = False,
    stratified_target_holdout_repeats: int | None = None,
    stratified_target_holdout_fraction: float = 0.2,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    model_types = model_types if model_types is not None else ["logistic_regression", "random_forest", "lightgbm"]
    group_cols = group_cols if group_cols is not None else ["target_id", "target_family"]
    leave_one_group_cols = leave_one_group_cols if leave_one_group_cols is not None else ["target_family"]
    raw = pd.read_csv(dataset_path, low_memory=False)
    raw["_atlas_observed_label"] = binary_label_series(raw[label_col])
    features, missing_features, all_missing_features = _present_features(
        raw,
        label_col=label_col,
        feature_set=feature_set,
        exclude_features=exclude_features,
        strict_feature_set=strict_feature_set,
    )
    data = raw.loc[raw["_atlas_observed_label"].notna()].copy()
    data[label_col] = data["_atlas_observed_label"].astype(int)
    _write_group_feasibility(data, label_col=label_col, group_cols=[*group_cols, *leave_one_group_cols], out_path=out)

    fold_metric_rows: list[dict[str, Any]] = []
    skipped_rows: list[dict[str, Any]] = []
    prediction_frames: list[pd.DataFrame] = []
    pooled_rows: list[dict[str, Any]] = []
    bootstrap_frames: list[pd.DataFrame] = []
    permutation_rows: list[dict[str, Any]] = []
    applicability_cache: dict[tuple[str, str, int, int], pd.DataFrame] = {}

    def run_fold(
        *,
        model_type: str,
        evaluation: str,
        group_col: str,
        fold: int,
        repeat: int,
        train_idx: pd.Index,
        test_idx: pd.Index,
        heldout_group: str,
    ) -> None:
        valid, split_counts = _split_valid(
            data,
            train_idx,
            test_idx,
            label_col=label_col,
            min_test_positives=min_test_positives,
            min_test_negatives=min_test_negatives,
        )
        if not valid:
            skipped_rows.append(
                {
                    "model_type": model_type,
                    "evaluation": evaluation,
                    "group_col": group_col,
                    "repeat": int(repeat),
                    "fold": int(fold),
                    "heldout_group": heldout_group,
                    "reason": "insufficient train/test class counts",
                    **split_counts,
                }
            )
            return
        applicability_key = (evaluation, group_col, int(repeat), int(fold))
        cached_applicability = applicability_cache.get(applicability_key)
        pred = _fit_predict_fold(
            data,
            train_idx=train_idx,
            test_idx=test_idx,
            features=features,
            label_col=label_col,
            model_type=model_type,
            seed=seed + repeat * 1000 + fold,
            class_weight=class_weight,
            pu_mode=pu_mode,
            model_params=model_params,
            compute_applicability_domain=(
                compute_applicability_domain and cached_applicability is None
            ),
        )
        applicability_cols = [
            "_source_index",
            "applicability_distance",
            "applicability_threshold_train_q95",
            "applicability_distance_ratio",
            "applicability_domain",
        ]
        if compute_applicability_domain:
            if cached_applicability is None:
                present = [col for col in applicability_cols if col in pred.columns]
                applicability_cache[applicability_key] = pred[present].copy()
            else:
                pred = pred.merge(cached_applicability, on="_source_index", how="left")
        pred["model_type"] = model_type
        pred["evaluation"] = evaluation
        pred["group_col"] = group_col
        pred["repeat"] = int(repeat)
        pred["fold"] = int(fold)
        pred["heldout_group"] = heldout_group
        prediction_frames.append(pred)
        row = _metric_row(pred, label_col=label_col, top_k=top_k)
        row.update(
            {
                "model_type": model_type,
                "evaluation": evaluation,
                "group_col": group_col,
                "repeat": int(repeat),
                "fold": int(fold),
                "heldout_group": heldout_group,
                **split_counts,
            }
        )
        fold_metric_rows.append(row)

    for model_type in model_types:
        for group_col in group_cols:
            if group_col not in data.columns:
                skipped_rows.append(
                    {
                        "model_type": model_type,
                        "evaluation": "repeated_grouped_cv",
                        "group_col": group_col,
                        "reason": "missing group column",
                    }
                )
                continue
            for repeat in range(int(repeats)):
                folds = _stratified_group_folds(
                    data,
                    label_col=label_col,
                    group_col=group_col,
                    n_splits=n_splits,
                    seed=seed + repeat,
                )
                for fold, (train_idx, test_idx) in enumerate(folds):
                    heldout = ";".join(sorted(data.loc[test_idx, group_col].fillna("missing").astype(str).unique()))
                    run_fold(
                        model_type=model_type,
                        evaluation="repeated_grouped_cv",
                        group_col=group_col,
                        repeat=repeat,
                        fold=fold,
                        train_idx=train_idx,
                        test_idx=test_idx,
                        heldout_group=heldout,
                    )
        for group_col in leave_one_group_cols:
            if group_col not in data.columns:
                skipped_rows.append(
                    {
                        "model_type": model_type,
                        "evaluation": "leave_one_group_out",
                        "group_col": group_col,
                        "reason": "missing group column",
                    }
                )
                continue
            for fold, (heldout, train_idx, test_idx) in enumerate(_leave_one_group_folds(data, group_col=group_col)):
                run_fold(
                    model_type=model_type,
                    evaluation="leave_one_group_out",
                    group_col=group_col,
                    repeat=0,
                    fold=fold,
                    train_idx=train_idx,
                    test_idx=test_idx,
                    heldout_group=heldout,
                )

        if run_stratified_target_holdout_by_family:
            family_folds = _stratified_target_holdout_by_family_folds(
                data,
                target_col="target_id",
                family_col="target_family",
                test_fraction=stratified_target_holdout_fraction,
                repeats=stratified_target_holdout_repeats if stratified_target_holdout_repeats is not None else repeats,
                seed=seed,
            )
            if not family_folds:
                skipped_rows.append(
                    {
                        "model_type": model_type,
                        "evaluation": "stratified_target_holdout_by_family",
                        "group_col": "target_id_within_target_family",
                        "reason": "missing target_id/target_family columns or no families with multiple targets",
                    }
                )
            for fold, (heldout, repeat_idx, train_idx, test_idx) in enumerate(family_folds):
                run_fold(
                    model_type=model_type,
                    evaluation="stratified_target_holdout_by_family",
                    group_col="target_id_within_target_family",
                    repeat=repeat_idx,
                    fold=fold,
                    train_idx=train_idx,
                    test_idx=test_idx,
                    heldout_group=heldout,
                )

    predictions = pd.concat(prediction_frames, ignore_index=True) if prediction_frames else pd.DataFrame()
    predictions.to_csv(out / "grouped_cv_oof_predictions.csv", index=False)
    fold_metrics = pd.DataFrame(fold_metric_rows)
    fold_metrics.to_csv(out / "grouped_cv_fold_metrics.csv", index=False)
    pd.DataFrame(skipped_rows).to_csv(out / "grouped_cv_skipped_folds.csv", index=False)
    summary = _summarize_fold_metrics(fold_metrics)
    summary.to_csv(out / "grouped_cv_stability_summary.csv", index=False)

    if not predictions.empty:
        for key, group in predictions.groupby(["model_type", "evaluation", "group_col"], dropna=False):
            model_type, evaluation, group_col = key
            mean_pred = _mean_oof_predictions(group, label_col=label_col)
            pred_path = out / f"pooled_oof_predictions__{model_type}__{evaluation}__{group_col}.csv"
            mean_pred.to_csv(pred_path, index=False)
            pooled = _metric_row(mean_pred, label_col=label_col, top_k=top_k)
            pooled.update(
                {
                    "model_type": model_type,
                    "evaluation": evaluation,
                    "group_col": group_col,
                    "prediction_file": str(pred_path),
                }
            )
            pooled_rows.append(pooled)
            boot = _bootstrap_ci(
                mean_pred,
                label_col=label_col,
                n_bootstraps=n_bootstraps,
                seed=seed,
                top_k=top_k,
                cluster_col=str(group_col) if str(group_col) in mean_pred.columns else None,
            )
            if not boot.empty:
                boot["bootstrap_unit"] = str(group_col) if str(group_col) in mean_pred.columns else "row"
                boot["model_type"] = model_type
                boot["evaluation"] = evaluation
                boot["group_col"] = group_col
                bootstrap_frames.append(boot)
            for row in _permutation_rows(
                mean_pred,
                label_col=label_col,
                n_permutations=n_permutations,
                seed=seed,
            ):
                row.update({"model_type": model_type, "evaluation": evaluation, "group_col": group_col})
                permutation_rows.append(row)
    pooled_frame = pd.DataFrame(pooled_rows)
    pooled_frame.to_csv(out / "pooled_oof_metrics.csv", index=False)
    bootstrap_frame = pd.concat(bootstrap_frames, ignore_index=True) if bootstrap_frames else pd.DataFrame()
    bootstrap_frame.to_csv(out / "pooled_oof_bootstrap_ci.csv", index=False)
    pd.DataFrame(permutation_rows).to_csv(out / "pooled_oof_permutation_p.csv", index=False)

    manifest = {
        "dataset_path": str(dataset_path),
        "label_col": label_col,
        "feature_set": feature_set,
        "features": features,
        "missing_features": missing_features,
        "all_missing_features": all_missing_features,
        "model_types": model_types,
        "group_cols": group_cols,
        "leave_one_group_cols": leave_one_group_cols,
        "n_splits": int(n_splits),
        "repeats": int(repeats),
        "min_test_positives": int(min_test_positives),
        "min_test_negatives": int(min_test_negatives),
        "top_k": int(top_k),
        "n_labelable": int(len(data)),
        "n_positive": int(data[label_col].eq(1).sum()),
        "n_negative": int(data[label_col].eq(0).sum()),
        "n_folds_run": int(len(fold_metrics)),
        "n_folds_skipped": int(len(skipped_rows)),
        "compute_applicability_domain": bool(compute_applicability_domain),
        "policy": (
            "Grouped CV holds out all rows for a target/group together. Metrics are reported from held-out "
            "predictions only; pooled OOF metrics average repeated predictions per row before scoring."
        ),
    }
    (out / "grouped_cv_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "manifest": manifest,
        "summary": summary,
        "pooled_metrics": pooled_frame,
        "fold_metrics": fold_metrics,
        "predictions": predictions,
    }
