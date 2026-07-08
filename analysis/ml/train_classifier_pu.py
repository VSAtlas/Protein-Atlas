from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.pu_correction import (
    DEFAULT_PU_STRATA,
    correct_elkan_noto_probabilities,
    estimate_elkan_noto_c,
    write_pu_selection_diagnostics,
    write_sar_selection_diagnostics,
)
from analysis.ml.train_classifier_core import (
    _design_matrix,
    _fit_model_with_optional_weights,
    _model,
    _model_probabilities,
    _sample_weights,
)


def _available_strata(*frames: pd.DataFrame) -> list[str]:
    cols = set().union(*(set(frame.columns) for frame in frames if frame is not None))
    return [col for col in DEFAULT_PU_STRATA if col in cols]


def _stratum_key(frame: pd.DataFrame, strata_cols: list[str]) -> pd.Series:
    if not strata_cols:
        return pd.Series("__all__", index=frame.index)
    values = frame[strata_cols].astype("object").where(pd.notna(frame[strata_cols]), "missing").astype(str)
    key = values[strata_cols[0]]
    for col in strata_cols[1:]:
        key = key.str.cat(values[col], sep="|")
    return key


def _sample_unlabeled_indices(
    *,
    positives: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    n_sample: int,
    rng: random.Random,
    stratified: bool,
) -> list[Any]:
    if not stratified:
        return rng.sample(list(unlabeled_pool.index), n_sample)
    strata_cols = _available_strata(positives, unlabeled_pool)
    if not strata_cols:
        return rng.sample(list(unlabeled_pool.index), n_sample)
    pos_keys = _stratum_key(positives, strata_cols)
    unlabeled_keys = _stratum_key(unlabeled_pool, strata_cols)
    desired = pos_keys.value_counts(normalize=True)
    selected: list[Any] = []
    available_by_key = {
        key: list(unlabeled_keys[unlabeled_keys.eq(key)].index)
        for key in unlabeled_keys.dropna().unique()
    }
    for key, frac in desired.items():
        pool = available_by_key.get(key, [])
        if not pool:
            continue
        take = min(len(pool), max(1, int(round(n_sample * float(frac)))))
        selected.extend(rng.sample(pool, take))
    if len(selected) < n_sample:
        remaining = [idx for idx in unlabeled_pool.index if idx not in set(selected)]
        if remaining:
            selected.extend(rng.sample(remaining, min(len(remaining), n_sample - len(selected))))
    if len(selected) > n_sample:
        selected = rng.sample(selected, n_sample)
    return selected


def _fit_bagging_pu(
    *,
    model_type: str,
    seed: int,
    class_weight: str | None,
    model_train: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    features: list[str],
    label_col: str,
    x_test: pd.DataFrame,
    out_path: Path,
    n_bags: int = 50,
    unlabeled_ratio: float = 1.0,
    stratified: bool = False,
) -> tuple[Any, list[float], dict[str, object]]:
    labels = pd.to_numeric(model_train[label_col], errors="coerce")
    positives = model_train.loc[labels.eq(1)].copy()
    measured_negatives = model_train.loc[labels.eq(0)].copy()
    if positives.empty or unlabeled_pool.empty:
        raise ValueError("bagging_pu requires training positives and fold-local unlabeled rows")
    rng = random.Random(seed)
    n_sample = min(len(unlabeled_pool), max(1, int(round(len(positives) * float(unlabeled_ratio)))))
    bag_probs: list[list[float]] = []
    bag_rows: list[dict[str, object]] = []
    last_model: Any | None = None
    for bag_idx in range(max(1, int(n_bags))):
        sampled_idx = _sample_unlabeled_indices(
            positives=positives,
            unlabeled_pool=unlabeled_pool,
            n_sample=n_sample,
            rng=rng,
            stratified=stratified,
        )
        sampled = unlabeled_pool.loc[sampled_idx].copy()
        sampled[label_col] = 0
        sampled["_sample_weight"] = 0.25
        sampled["negative_evidence_type"] = "pu_temporary_unlabeled_negative"
        train_frame = pd.concat([positives, measured_negatives, sampled], ignore_index=False)
        x_bag = _design_matrix(train_frame, features)
        x_bag_test = x_test.reindex(columns=x_bag.columns, fill_value=0)
        model = _model(model_type, seed + bag_idx, class_weight=class_weight)
        weights = _sample_weights(train_frame, label_col, pu_mode="positive_unlabeled_weighted")
        _fit_model_with_optional_weights(
            model,
            x_bag,
            train_frame[label_col].astype(int),
            weights.reindex(train_frame.index) if weights is not None else None,
        )
        bag_probs.append(_model_probabilities(model, x_bag_test))
        last_model = model
        bag_rows.append(
            {
                "bag": bag_idx,
                "seed": seed + bag_idx,
                "n_positive": int(len(positives)),
                "n_measured_negative": int(len(measured_negatives)),
                "n_unlabeled_sampled_as_temporary_negative": int(len(sampled)),
                "sampled_unlabeled_index_hash": str(abs(hash(tuple(sorted(str(i) for i in sampled_idx))))),
            }
        )
    frame = pd.DataFrame(bag_probs)
    probs = frame.mean(axis=0).astype(float).tolist()
    variance = frame.var(axis=0).fillna(0.0).astype(float).tolist()
    pd.DataFrame(bag_rows).to_csv(out_path / "bagging_pu_members.csv", index=False)
    manifest: dict[str, object] = {
        "status": "ok",
        "mode": "stratified_bagging_pu" if stratified else "bagging_pu",
        "n_bags": int(max(1, int(n_bags))),
        "n_train_positive": int(len(positives)),
        "n_train_measured_negative": int(len(measured_negatives)),
        "n_unlabeled_pool_fold_local": int(len(unlabeled_pool)),
        "n_unlabeled_sampled_per_bag": int(n_sample),
        "stratified_sampling": bool(stratified),
        "strata_columns_used": _available_strata(positives, unlabeled_pool) if stratified else [],
        "unlabeled_sample_weight": 0.25,
        "mean_prediction_variance": float(pd.Series(variance).mean()) if variance else None,
        "policy": "unlabeled rows are sampled as temporary negatives inside the training fold only",
    }
    (out_path / "pu_training_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return last_model or _model(model_type, seed, class_weight=class_weight), probs, manifest


def _selection_design_matrix(
    positive_frame: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    features: list[str],
) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    strata_cols = _available_strata(positive_frame, unlabeled_pool)
    cols = [col for col in [*features, *strata_cols] if col in positive_frame.columns or col in unlabeled_pool.columns]
    pos = positive_frame.reindex(columns=cols).copy()
    unlabeled = unlabeled_pool.reindex(columns=cols).copy()
    pos["_pu_observed_positive_label"] = 1
    unlabeled["_pu_observed_positive_label"] = 0
    frame = pd.concat([pos, unlabeled], ignore_index=True)
    x = _design_matrix(frame, cols)
    y = frame["_pu_observed_positive_label"].astype(int)
    return x, y, cols


def _fit_propensity_weighted_pu(
    *,
    model_type: str,
    seed: int,
    class_weight: str | None,
    model_train: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    features: list[str],
    label_col: str,
    x_test: pd.DataFrame,
    out_path: Path,
) -> tuple[Any, list[float], dict[str, object]]:
    labels = pd.to_numeric(model_train[label_col], errors="coerce")
    positives = model_train.loc[labels.eq(1)].copy()
    measured_negatives = model_train.loc[labels.eq(0)].copy()
    if positives.empty or measured_negatives.empty or unlabeled_pool.empty:
        raise ValueError("propensity_weighted_pu requires positives, measured negatives, and fold-local unlabeled rows")
    try:
        from sklearn.linear_model import LogisticRegression
    except ImportError as exc:
        raise ValueError("propensity_weighted_pu requires scikit-learn") from exc
    x_sel, y_sel, selection_cols = _selection_design_matrix(positives, unlabeled_pool, features)
    propensity_model = LogisticRegression(max_iter=1000, class_weight="balanced")
    propensity_model.fit(x_sel, y_sel)
    x_pos_sel = _design_matrix(positives, selection_cols).reindex(columns=x_sel.columns, fill_value=0)
    propensity = pd.Series(propensity_model.predict_proba(x_pos_sel)[:, 1], index=positives.index).clip(0.05, 1.0)
    inv = (1.0 / propensity).clip(upper=10.0)
    inv = inv / inv.mean()
    train_frame = model_train.copy()
    weights = pd.Series(1.0, index=train_frame.index)
    weights.loc[positives.index] = inv
    existing = _sample_weights(train_frame, label_col, pu_mode="positive_unlabeled_weighted")
    if existing is not None:
        weights = weights * existing.reindex(train_frame.index).fillna(1.0)
    model = _model(model_type, seed, class_weight=class_weight)
    x_train_local = _design_matrix(train_frame, features)
    x_test_local = x_test.reindex(columns=x_train_local.columns, fill_value=0)
    _fit_model_with_optional_weights(model, x_train_local, train_frame[label_col].astype(int), weights)
    probs = _model_probabilities(model, x_test_local)
    write_pu_selection_diagnostics(
        train_frame=model_train,
        unlabeled_pool=unlabeled_pool,
        label_col=label_col,
        out_path=out_path,
    )
    pd.DataFrame(
        {
            "row_index": [str(idx) for idx in positives.index],
            "propensity_p_labeled_given_features": propensity.astype(float).tolist(),
            "positive_inverse_propensity_weight": inv.astype(float).tolist(),
        }
    ).to_csv(out_path / "pu_propensity_weights.csv", index=False)
    manifest = {
        "status": "ok",
        "mode": "propensity_weighted_pu",
        "n_train_positive": int(len(positives)),
        "n_train_measured_negative": int(len(measured_negatives)),
        "n_unlabeled_pool_fold_local": int(len(unlabeled_pool)),
        "selection_columns": selection_cols,
        "mean_positive_propensity": float(propensity.mean()),
        "median_positive_propensity": float(propensity.median()),
        "max_inverse_propensity_weight": float(inv.max()),
        "policy": "Known positives are inverse-propensity weighted; measured inactives remain measured negatives and unlabeled rows remain unlabeled.",
    }
    (out_path / "pu_training_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return model, probs, manifest


def _fit_pulsnar_style_pu(
    *,
    model_type: str,
    seed: int,
    class_weight: str | None,
    model_train: pd.DataFrame,
    validation: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    features: list[str],
    label_col: str,
    x_test: pd.DataFrame,
    out_path: Path,
) -> tuple[Any, list[float], dict[str, object]]:
    labels = pd.to_numeric(model_train[label_col], errors="coerce")
    positives = model_train.loc[labels.eq(1)].copy()
    measured_negatives = model_train.loc[labels.eq(0)].copy()
    if len(positives) < 10 or unlabeled_pool.empty:
        raise ValueError("pulsnar_style requires at least 10 positives and fold-local unlabeled rows")
    try:
        from sklearn.cluster import KMeans
    except ImportError as exc:
        raise ValueError("pulsnar_style requires scikit-learn clustering") from exc
    x_pos = _design_matrix(positives, features)
    n_clusters = max(1, min(5, len(positives) // 20))
    clusterer = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    clusters = clusterer.fit_predict(x_pos)
    cluster_rows: list[dict[str, object]] = []
    cluster_probs: list[list[float]] = []
    last_model: Any | None = None
    for cluster_id in range(n_clusters):
        cluster_positive = positives.iloc[[idx for idx, val in enumerate(clusters) if int(val) == cluster_id]].copy()
        if cluster_positive.empty:
            continue
        pu_train = pd.concat([cluster_positive, unlabeled_pool], ignore_index=False).copy()
        pu_train["_pu_observed_positive_label"] = 0
        pu_train.loc[cluster_positive.index, "_pu_observed_positive_label"] = 1
        x_pu = _design_matrix(pu_train, features)
        x_pu_test = x_test.reindex(columns=x_pu.columns, fill_value=0)
        model = _model(model_type, seed + cluster_id, class_weight=class_weight)
        model.fit(x_pu, pu_train["_pu_observed_positive_label"].astype(int))
        raw = _model_probabilities(model, x_pu_test)
        validation_positive = validation.loc[pd.to_numeric(validation.get(label_col), errors="coerce").eq(1)].copy()
        c_source = "training_cluster_positives"
        c_frame = cluster_positive
        if not validation_positive.empty:
            c_source = "heldout_validation_positives_all_clusters"
            c_frame = validation_positive
        x_c = _design_matrix(c_frame, features).reindex(columns=x_pu.columns, fill_value=0)
        c_scores = c_frame.copy()
        c_scores["_cluster_observed_positive_probability"] = _model_probabilities(model, x_c)
        c_manifest = estimate_elkan_noto_c(
            positive_score_frame=c_scores,
            score_col="_cluster_observed_positive_probability",
            out_path=None,
        )
        corrected = correct_elkan_noto_probabilities(raw, c_estimate=c_manifest.get("global_c"))
        cluster_probs.append(corrected)
        last_model = model
        cluster_rows.append(
            {
                "cluster_id": cluster_id,
                "n_cluster_positive": int(len(cluster_positive)),
                "n_unlabeled_pool_fold_local": int(len(unlabeled_pool)),
                "c_source": c_source,
                "cluster_c": c_manifest.get("global_c"),
            }
        )
    if not cluster_probs:
        raise ValueError("pulsnar_style failed to train any positive clusters")
    probs = pd.DataFrame(cluster_probs).max(axis=0).clip(0.0, 1.0).astype(float).tolist()
    pd.DataFrame(cluster_rows).to_csv(out_path / "pu_pulsnar_clusters.csv", index=False)
    write_pu_selection_diagnostics(
        train_frame=model_train,
        unlabeled_pool=unlabeled_pool,
        label_col=label_col,
        out_path=out_path,
    )
    sar_manifest = write_sar_selection_diagnostics(
        positive_frame=positives,
        unlabeled_pool=unlabeled_pool,
        out_path=out_path,
        seed=seed,
    )
    manifest = {
        "status": "ok",
        "mode": "pulsnar_style",
        "n_positive_clusters": int(len(cluster_rows)),
        "n_train_positive": int(len(positives)),
        "n_train_measured_negative_excluded_from_pu_denominator": int(len(measured_negatives)),
        "n_unlabeled_pool_fold_local": int(len(unlabeled_pool)),
        "sar_warning": sar_manifest.get("sar_warning"),
        "policy": "PULSNAR-style sensitivity: positives are clustered, one positive-vs-unlabeled model is trained per cluster, and corrected probabilities are aggregated by max.",
    }
    (out_path / "pu_training_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return last_model or _model(model_type, seed, class_weight=class_weight), probs, manifest


def _fit_elkan_noto_pu(
    *,
    model_type: str,
    seed: int,
    class_weight: str | None,
    model_train: pd.DataFrame,
    validation: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    features: list[str],
    label_col: str,
    x_test: pd.DataFrame,
    out_path: Path,
) -> tuple[Any, list[float], dict[str, object]]:
    labels = pd.to_numeric(model_train[label_col], errors="coerce")
    positives = model_train.loc[labels.eq(1)].copy()
    measured_negatives = model_train.loc[labels.eq(0)].copy()
    if positives.empty or unlabeled_pool.empty:
        raise ValueError("elkan_noto PU requires training positives and fold-local unlabeled rows")

    pu_train = pd.concat([positives, unlabeled_pool], ignore_index=False).copy()
    pu_train["_pu_observed_positive_label"] = 0
    pu_train.loc[positives.index, "_pu_observed_positive_label"] = 1
    x_pu = _design_matrix(pu_train, features)
    x_pu_test = x_test.reindex(columns=x_pu.columns, fill_value=0)
    model = _model(model_type, seed, class_weight=class_weight)
    model.fit(x_pu, pu_train["_pu_observed_positive_label"].astype(int))
    raw_probs = _model_probabilities(model, x_pu_test)

    validation_positive = validation.loc[pd.to_numeric(validation.get(label_col), errors="coerce").eq(1)].copy()
    c_score_source = "heldout_validation_positives"
    if validation_positive.empty:
        validation_positive = positives.copy()
        c_score_source = "training_positives_fallback"
    x_c = _design_matrix(validation_positive, features).reindex(columns=x_pu.columns, fill_value=0)
    positive_scores = validation_positive.copy()
    positive_scores["_pu_observed_positive_probability"] = _model_probabilities(model, x_c)
    c_manifest = estimate_elkan_noto_c(
        positive_score_frame=positive_scores,
        score_col="_pu_observed_positive_probability",
        out_path=out_path,
    )
    sar_manifest = write_sar_selection_diagnostics(
        positive_frame=positives,
        unlabeled_pool=unlabeled_pool,
        out_path=out_path,
        seed=seed,
    )
    global_c = c_manifest.get("global_c")
    corrected_probs = correct_elkan_noto_probabilities(raw_probs, c_estimate=float(global_c) if global_c is not None else None)
    write_pu_selection_diagnostics(
        train_frame=model_train,
        unlabeled_pool=unlabeled_pool,
        label_col=label_col,
        out_path=out_path,
    )
    manifest = {
        "status": "ok",
        "mode": "elkan_noto",
        "n_train_positive": int(len(positives)),
        "n_train_measured_negative_excluded_from_pu_denominator": int(len(measured_negatives)),
        "n_unlabeled_pool_fold_local": int(len(unlabeled_pool)),
        "c_score_source": c_score_source,
        "global_c": global_c,
        "scar_warning": c_manifest.get("scar_warning"),
        "sar_warning": sar_manifest.get("sar_warning"),
        "policy": (
            "Elkan-Noto trains observed-positive vs fold-local unlabeled examples. "
            "Measured inactives are tracked separately and are not mixed into the unlabeled denominator."
        ),
    }
    (out_path / "pu_training_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return model, corrected_probs, manifest
