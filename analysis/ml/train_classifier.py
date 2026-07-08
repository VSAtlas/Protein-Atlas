from __future__ import annotations

from collections.abc import Mapping
import json
import pickle
from pathlib import Path

import pandas as pd

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.applicability_domain import append_abstention_flags, append_chemical_fingerprint_ad, append_target_family_ad
from analysis.ml.artifact_manifest import write_artifact_manifest
from analysis.ml.baseline_panel import write_standard_baseline_panel
from analysis.ml.cards import write_model_card
from analysis.ml.claim_readiness import write_model_claim_readiness
from analysis.ml.conformal import append_split_conformal_sets
from analysis.ml.dataset_manifest import write_dataset_version_manifest
from analysis.ml.decision_metrics import write_decision_metrics, write_group_topk_recovery
from analysis.ml.experiment_tracking import track_model_run
from analysis.ml.explain import write_feature_importance
from analysis.ml.feature_sets import RETAINED_CONTEXT_AUDIT_COLUMNS, effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.model_run_ledger import write_model_run_record
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.split_manifest import load_locked_split_manifest, write_split_manifest
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier_core import (
    _design_matrix,
    _fit_design_matrix,
    _model,
    _fit_model_with_optional_weights,
    _model_probabilities,
    _sample_weights,
    combine_sample_weights,
    group_reweighting_weights,
    _select_model,
    _transform_design_matrix,
    write_preprocessing_artifacts,
)
from analysis.ml.train_classifier_outputs import (
    _append_applicability_domain,
    _bootstrap_metric_ci,
    _write_grouped_calibration_outputs,
    _write_reliability_outputs,
    _write_source_transfer_diagnostics,
    write_subgroup_metrics,
)
from analysis.ml.train_classifier_pu import (
    _fit_bagging_pu,
    _fit_elkan_noto_pu,
    _fit_propensity_weighted_pu,
    _fit_pulsnar_style_pu,
)
from analysis.ml.train_classifier_splits import (
    _custom_split,
    _exclude_unlabeled_holdout_entities,
    _fixed_validation_split,
    _temporal_split,
)

__all__ = ["train_ml_model", "_design_matrix", "_model"]


def train_ml_model(
    dataset_path: str | Path,
    label_col: str,
    feature_set: str,
    model_type: str,
    split_mode: str,
    out_dir: str | Path,
    seed: int = 42,
    split_column: str | None = None,
    train_values: list[str] | None = None,
    test_values: list[str] | None = None,
    temporal_year_col: str | None = None,
    temporal_cutoff_year: int | None = None,
    n_bootstraps: int = 200,
    pu_mode: str = "standard_binary",
    class_weight: str | None = "balanced",
    exclude_features: list[str] | None = None,
    allow_label_definition_features: bool = False,
    allow_partial_rescoring_features: bool = False,
    validation_fold_col: str | None = None,
    validation_fold_value: str | None = None,
    validation_fraction: float = 0.15,
    nested_model_selection: bool = False,
    temporal_max_missing_fraction: float = 0.20,
    pu_bags: int = 50,
    pu_unlabeled_ratio: float = 1.0,
    claim_mode: str = "exploratory",
    repo_root: str | Path | None = None,
    model_params: Mapping[str, object] | None = None,
    hpo_metadata: dict[str, object] | None = None,
    split_manifest: str | Path | None = None,
    dataset_provenance: dict[str, object] | None = None,
    calibration_method: str = "none",
    strict_feature_set: bool = False,
    group_reweight_cols: list[str] | None = None,
    group_reweight_include_label: bool = True,
    group_reweight_max_factor: float = 5.0,
) -> dict[str, object]:
    df = pd.read_csv(dataset_path, low_memory=False)
    df["_atlas_observed_label"] = binary_label_series(df[label_col])
    excluded = effective_exclude_features(
        label_col,
        exclude_features,
        allow_label_definition_features=allow_label_definition_features,
    )
    requested_features = get_feature_set(feature_set)
    effective_requested_features = [feature for feature in requested_features if feature not in excluded]
    retained_excluded_features = sorted(feature for feature in excluded if feature in df.columns)
    missing_features = [feature for feature in effective_requested_features if feature not in df.columns]
    present_requested_features = [feature for feature in effective_requested_features if feature in df.columns]
    all_missing_features = [feature for feature in present_requested_features if df[feature].notna().sum() == 0]
    if strict_feature_set and (missing_features or all_missing_features):
        parts = []
        if missing_features:
            preview = ", ".join(missing_features[:20])
            suffix = "" if len(missing_features) <= 20 else f" ... +{len(missing_features) - 20} more"
            parts.append(f"{len(missing_features)} missing columns: {preview}{suffix}")
        if all_missing_features:
            preview = ", ".join(all_missing_features[:20])
            suffix = "" if len(all_missing_features) <= 20 else f" ... +{len(all_missing_features) - 20} more"
            parts.append(f"{len(all_missing_features)} all-missing columns: {preview}{suffix}")
        raise ValueError(f"strict feature-set check failed for {feature_set!r}: " + "; ".join(parts))
    features = [feature for feature in effective_requested_features if feature in df.columns]
    retained_context_not_trained = sorted(
        col for col in RETAINED_CONTEXT_AUDIT_COLUMNS if col in df.columns and col not in features
    )
    assert_no_leakage(
        features,
        allow_label_definition_features=allow_label_definition_features,
        allow_partial_rescoring_features=allow_partial_rescoring_features,
    )
    data = df.loc[df["_atlas_observed_label"].notna()].copy()
    data[label_col] = data["_atlas_observed_label"]
    data = data.dropna(subset=[label_col]).copy()
    data[label_col] = data[label_col].astype(int)
    locked_validation_idx = pd.Index([])
    locked_split_metadata: dict[str, object] | None = None
    if split_manifest is not None:
        locked = load_locked_split_manifest(data, split_manifest, label_col=label_col, split_mode=split_mode)
        train_idx = locked["train_idx"]  # type: ignore[assignment]
        test_idx = locked["test_idx"]  # type: ignore[assignment]
        locked_validation_idx = locked["validation_idx"]  # type: ignore[assignment]
        split_summary = dict(locked.get("split_summary") or {})
        locked_split_metadata = dict(locked.get("metadata") or {})
        split_summary.setdefault("split_mode", split_mode)
        split_summary.setdefault("locked_split_manifest", str(split_manifest))
        split_summary.setdefault("passes_holdout", True)
        split_summary.setdefault("overlaps", {})
    else:
        temporal = _temporal_split(
            data,
            temporal_year_col,
            temporal_cutoff_year,
            max_missing_fraction=temporal_max_missing_fraction,
        )
        custom = _custom_split(data, split_column, train_values, test_values)
        split_summary: dict[str, object] | None
        if temporal is not None:
            train_idx, test_idx, split_summary = temporal
        elif custom is not None:
            train_idx, test_idx, split_summary = custom
        else:
            train_idx, test_idx = make_split(data, split_mode=split_mode, seed=seed)
            split_summary = None
    train = data.loc[train_idx]
    test = data.loc[test_idx]
    if split_summary is None:
        split_summary = split_overlap_summary(train, test, split_mode)
    if split_summary["overlaps"] and not split_summary["passes_holdout"]:
        raise ValueError(f"holdout split leakage detected: {split_summary}")
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    if model_params is not None and pu_mode in {
        "bagging_pu",
        "stratified_bagging_pu",
        "propensity_weighted_pu",
        "elkan_noto",
        "pulsnar_style",
    }:
        raise ValueError("external model_params are supported only for supervised or weighted binary PU modes")
    dataset_manifest = write_dataset_version_manifest(
        dataset_path=dataset_path,
        frame=df,
        out_path=out_path / "dataset_version_manifest.json",
        label_col=label_col,
        feature_set=feature_set,
        features=features,
        requested_features=effective_requested_features,
        missing_features=missing_features,
        all_missing_features=all_missing_features,
        exclude_features=sorted(excluded),
        retained_excluded_features=retained_excluded_features,
        retained_context_not_trained=retained_context_not_trained,
        provenance={"stage": "train_ml_model", "split_mode": split_mode, **dict(dataset_provenance or {})},
    )
    unlabeled_pool = _exclude_unlabeled_holdout_entities(
        df.loc[df["_atlas_observed_label"].isna()].copy(),
        test,
        split_mode,
    )
    if len(locked_validation_idx):
        model_train = train.copy()
        validation = data.loc[locked_validation_idx].copy()
        validation_idx = locked_validation_idx
    else:
        model_train, validation, validation_idx = _fixed_validation_split(
            train,
            validation_fold_col=validation_fold_col,
            validation_fold_value=validation_fold_value,
            seed=seed,
            validation_fraction=validation_fraction if nested_model_selection or validation_fold_col or hpo_metadata is not None else 0.0,
        )
    x_train, preprocessing = _fit_design_matrix(model_train, features)
    x_test = _transform_design_matrix(test, preprocessing)
    x_validation = _transform_design_matrix(validation, preprocessing) if not validation.empty else pd.DataFrame()
    preprocessing_manifest = write_preprocessing_artifacts(
        out_path=out_path,
        preprocessing=preprocessing,
        x_train=x_train,
    )
    selected_model_name = model_type
    if nested_model_selection or validation_fold_col:
        selected_model_name, model, validation_metrics = _select_model(
            model_type,
            seed,
            class_weight,
            x_train,
            model_train[label_col],
            x_validation,
            validation[label_col] if not validation.empty else pd.Series(dtype=int),
            model_params=model_params,
        )
    else:
        model = _model(model_type, seed, class_weight=class_weight, model_params=model_params)
        validation_metrics = pd.DataFrame()
    group_weights, group_reweight_manifest = group_reweighting_weights(
        model_train,
        label_col,
        group_cols=group_reweight_cols,
        include_label=group_reweight_include_label,
        max_factor=group_reweight_max_factor,
    )
    (out_path / "group_reweighting_manifest.json").write_text(
        json.dumps(group_reweight_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    if group_weights is not None:
        pd.DataFrame({"row_index": group_weights.index.astype(str), "group_reweight": group_weights.astype(float).values}).to_csv(
            out_path / "group_reweighting_weights.csv",
            index=False,
        )

    if pu_mode in {"bagging_pu", "stratified_bagging_pu"}:
        model, probs, pu_manifest = _fit_bagging_pu(
            model_type=selected_model_name if selected_model_name in {"logistic_regression", "random_forest", "xgboost", "gradient_boosted_trees"} else model_type,
            seed=seed,
            class_weight=class_weight,
            model_train=model_train,
            unlabeled_pool=unlabeled_pool,
            features=features,
            label_col=label_col,
            x_test=x_test,
            out_path=out_path,
            n_bags=pu_bags,
            unlabeled_ratio=pu_unlabeled_ratio,
            stratified=pu_mode == "stratified_bagging_pu",
        )
        validation_metrics = pd.concat(
            [validation_metrics, pd.DataFrame([{"candidate": pu_mode, **pu_manifest}])],
            ignore_index=True,
        )
    elif pu_mode == "propensity_weighted_pu":
        model, probs, pu_manifest = _fit_propensity_weighted_pu(
            model_type=selected_model_name if selected_model_name in {"logistic_regression", "random_forest", "xgboost", "gradient_boosted_trees"} else model_type,
            seed=seed,
            class_weight=class_weight,
            model_train=model_train,
            unlabeled_pool=unlabeled_pool,
            features=features,
            label_col=label_col,
            x_test=x_test,
            out_path=out_path,
        )
        validation_metrics = pd.concat(
            [validation_metrics, pd.DataFrame([{"candidate": "propensity_weighted_pu", **pu_manifest}])],
            ignore_index=True,
        )
    elif pu_mode == "elkan_noto":
        model, probs, pu_manifest = _fit_elkan_noto_pu(
            model_type=selected_model_name if selected_model_name in {"logistic_regression", "random_forest", "xgboost", "gradient_boosted_trees"} else model_type,
            seed=seed,
            class_weight=class_weight,
            model_train=model_train,
            validation=validation,
            unlabeled_pool=unlabeled_pool,
            features=features,
            label_col=label_col,
            x_test=x_test,
            out_path=out_path,
        )
        validation_metrics = pd.concat(
            [validation_metrics, pd.DataFrame([{"candidate": "elkan_noto", **pu_manifest}])],
            ignore_index=True,
        )
    elif pu_mode == "pulsnar_style":
        model, probs, pu_manifest = _fit_pulsnar_style_pu(
            model_type=selected_model_name if selected_model_name in {"logistic_regression", "random_forest", "xgboost", "gradient_boosted_trees"} else model_type,
            seed=seed,
            class_weight=class_weight,
            model_train=model_train,
            validation=validation,
            unlabeled_pool=unlabeled_pool,
            features=features,
            label_col=label_col,
            x_test=x_test,
            out_path=out_path,
        )
        validation_metrics = pd.concat(
            [validation_metrics, pd.DataFrame([{"candidate": "pulsnar_style", **pu_manifest}])],
            ignore_index=True,
        )
    else:
        base_weight = _sample_weights(model_train, label_col, pu_mode=pu_mode)
        sample_weight = combine_sample_weights(base_weight, group_weights)
        _fit_model_with_optional_weights(model, x_train, model_train[label_col], sample_weight)
        probs = _model_probabilities(model, x_test)
    if pu_mode not in {"bagging_pu", "stratified_bagging_pu", "propensity_weighted_pu", "elkan_noto", "pulsnar_style"}:
        (out_path / "pu_training_manifest.json").write_text(
            json.dumps(
                {
                    "status": "not_applicable",
                    "mode": pu_mode,
                    "n_unlabeled_pool_fold_local": int(len(unlabeled_pool)),
                    "policy": "no fold-local PU bagging requested",
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
    raw_probs = [float(v) for v in probs]
    calibration_name = str(calibration_method or "none").strip().lower()
    calibration_manifest: dict[str, object] = {
        "method": calibration_name,
        "status": "not_requested" if calibration_name in {"", "none", "raw"} else "skipped",
        "reason": None,
    }
    if calibration_name in {"sigmoid", "isotonic"}:
        if pu_mode in {"bagging_pu", "stratified_bagging_pu", "propensity_weighted_pu", "elkan_noto", "pulsnar_style"}:
            calibration_manifest.update({"status": "skipped", "reason": "calibration currently applies only to direct supervised model outputs"})
        elif validation.empty or x_validation.empty:
            calibration_manifest.update({"status": "skipped", "reason": "no validation/calibration split available"})
        elif pd.to_numeric(validation[label_col], errors="coerce").nunique() < 2:
            calibration_manifest.update({"status": "skipped", "reason": "calibration split does not contain both classes"})
        else:
            validation_raw = _model_probabilities(model, x_validation)
            if calibration_name == "sigmoid":
                from sklearn.linear_model import LogisticRegression

                calibrator = LogisticRegression(max_iter=1000)
                calibrator.fit(pd.DataFrame({"score": validation_raw}), validation[label_col].astype(int))
                probs = calibrator.predict_proba(pd.DataFrame({"score": raw_probs}))[:, 1].astype(float).tolist()
            else:
                from sklearn.isotonic import IsotonicRegression

                calibrator = IsotonicRegression(out_of_bounds="clip")
                calibrator.fit(validation_raw, validation[label_col].astype(int))
                probs = [float(v) for v in calibrator.predict(raw_probs)]
            with (out_path / "probability_calibrator.pkl").open("wb") as handle:
                pickle.dump(calibrator, handle)
            calibration_manifest.update(
                {
                    "status": "fit",
                    "n_calibration": int(len(validation)),
                    "calibration_role": (
                        "validation_model_selection_and_calibration"
                        if nested_model_selection or validation_fold_col or hpo_metadata is not None
                        else "pure_calibration"
                    ),
                }
            )
    elif calibration_name not in {"", "none", "raw"}:
        raise ValueError(f"unsupported calibration_method: {calibration_method}")
    (out_path / "probability_calibration_manifest.json").write_text(
        json.dumps(calibration_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    prediction_metadata = [
        "drug_id",
        "target_id",
        "pdb_id",
        "scaffold_key",
        "ligand_chemotype",
        "label_source",
        "source_family",
        "upstream_source",
        "external_source_family",
        "external_upstream_source",
        "external_evidence_sources",
        "external_parent_sources",
        "mechanism_label_source",
        "assay_type",
        "endpoint_type",
        "activity_type",
        "target_family",
        "protein_family",
        "chemical_cluster",
        "ligand_cluster",
        "smiles",
        "canonical_smiles",
        "database_release_year",
        "activity_publication_year",
        "activity_publication_year_source",
        label_col,
    ]
    pred = test[[col for col in prediction_metadata if col in test.columns]].copy()
    if probs != raw_probs:
        pred["raw_ml_prediction_score"] = raw_probs
    pred["ml_prediction_score"] = probs
    conformal_status = "skipped"
    if pu_mode not in {"bagging_pu", "stratified_bagging_pu", "propensity_weighted_pu", "elkan_noto", "pulsnar_style"} and not validation.empty and not x_validation.empty:
        pred = append_split_conformal_sets(
            pred,
            model=model,
            x_calibration=x_validation,
            y_calibration=validation[label_col],
            alpha=0.10,
            label_col=label_col,
            out_path=out_path,
            calibration_role=(
                "validation_model_selection_and_calibration"
                if nested_model_selection or validation_fold_col or hpo_metadata is not None
                else "pure_calibration"
            ),
        )
        try:
            summary = json.loads((out_path / "conformal_summary.json").read_text(encoding="utf-8"))
            conformal_status = str(summary.get("status", "skipped"))
        except Exception:
            conformal_status = "skipped"
    else:
        (out_path / "conformal_summary.json").write_text(
            json.dumps(
                {
                    "status": "skipped",
                    "reason": "no_validation_or_calibration_split",
                    "alpha": 0.10,
                    "calibration_role": "none",
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
    pred = _append_applicability_domain(pred, x_train, x_test, label_col, out_path=out_path)
    pred = append_chemical_fingerprint_ad(pred, model_train, test, label_col=label_col, out_path=out_path)
    pred = append_target_family_ad(pred, model_train, test, label_col=label_col, out_path=out_path)
    pred = append_abstention_flags(pred, out_path=out_path)
    metrics: dict[str, object] = dict(calibration_metrics([float(v) for v in probs], test[label_col].tolist()))
    metrics["selected_model"] = selected_model_name
    metrics["label_definition_leakage_policy"] = (
        "allowed_for_explicit_sensitivity_model"
        if allow_label_definition_features
        else "auto_excluded_for_label_definition"
    )
    pred.to_csv(out_path / "model_predictions.csv", index=False)
    pd.DataFrame([{"metric": key, "value": value} for key, value in metrics.items()]).to_csv(out_path / "model_metrics.csv", index=False)
    hpo_trials = validation_metrics.copy()
    if not hpo_trials.empty:
        hpo_trials.to_csv(out_path / "model_selection_validation_metrics.csv", index=False)
    hpo_trials.to_csv(out_path / "hpo_trials.csv", index=False)
    try:
        effective_model_params: dict[str, object] = dict(model.get_params(deep=False))
    except Exception:
        effective_model_params = dict(model_params or {})
    if hpo_metadata is not None:
        hpo_n_trials = hpo_metadata.get("n_trials", 0)
        if not isinstance(hpo_n_trials, int | float | str):
            hpo_n_trials = 0
        hpo_manifest = {
            "status": "run",
            "backend": "external",
            "model_type": model_type,
            "model_params": dict(model_params or effective_model_params),
            "calibration_method": calibration_name,
            "strict_feature_set": strict_feature_set,
            "all_missing_features": all_missing_features,
            "selected_model": selected_model_name,
            "selection_metric": hpo_metadata.get("selection_metric", "AUPRC"),
            "n_trials": int(hpo_n_trials or 0),
            "trials_path": str(hpo_metadata.get("trials_path", out_path / "hpo_trials.csv")),
            "best_params": dict(model_params or {}),
            "effective_model_params": effective_model_params,
            "policy": "external HPO selected hyperparameters on a fixed validation split; outer holdout is reserved for final scoring",
        }
        hpo_manifest.update(hpo_metadata)
        hpo_best_params: Mapping[str, object] = model_params or {}
        if model_params is None:
            raw_best_params = hpo_manifest.get("best_params", {})
            if isinstance(raw_best_params, Mapping):
                hpo_best_params = {str(key): value for key, value in raw_best_params.items()}
        hpo_manifest["best_params"] = dict(hpo_best_params)
    else:
        hpo_manifest = {
            "status": "run" if not hpo_trials.empty else "not_requested",
            "model_type": model_type,
            "selected_model": selected_model_name,
            "selection_metric": "AUPRC",
            "n_trials": int(len(hpo_trials)),
            "trials_path": str(out_path / "hpo_trials.csv"),
            "policy": "local validation split over built-in candidates; use atlas ml hpo for Optuna-backed sweeps",
            "effective_model_params": effective_model_params,
        }
    (out_path / "hpo_manifest.json").write_text(
        json.dumps(hpo_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _bootstrap_metric_ci(pred, label_col, n_bootstraps=n_bootstraps, seed=seed).to_csv(
        out_path / "model_metric_bootstrap_ci.csv",
        index=False,
    )
    _write_reliability_outputs(pred, label_col, out_path)
    _write_grouped_calibration_outputs(pred, label_col, out_path)
    subgroup_metrics = write_subgroup_metrics(pred, label_col, out_path)
    write_decision_metrics(pred, label_col=label_col, score_col="ml_prediction_score", out_path=out_path / "model_decision_metrics.csv")
    write_group_topk_recovery(
        pred,
        label_col=label_col,
        score_col="ml_prediction_score",
        out_path=out_path / "model_group_topk_recovery.csv",
    )
    baseline_frame = test.copy()
    baseline_frame["ml_prediction_score"] = probs
    write_standard_baseline_panel(baseline_frame, label_col=label_col, out_path=out_path / "standard_baseline_panel.csv", seed=seed)
    _write_source_transfer_diagnostics(train, test, pred, label_col, split_column, out_path)
    write_feature_importance(model, list(x_train.columns), out_path / "feature_importance.csv")
    write_split_manifest(
        data,
        model_train.index,
        test.index,
        out_path,
        label_col=label_col,
        split_mode=split_mode,
        split_summary=split_summary,
        validation_idx=validation_idx,
    )
    with (out_path / "split_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(split_summary, handle, indent=2, sort_keys=True)
    readiness = write_model_claim_readiness(
        out_path=out_path,
        pred=pred,
        label_col=label_col,
        metrics=metrics,
        split_summary=split_summary,
        features=features,
        validation_rows=len(validation),
        conformal_status=conformal_status,
        baseline_path=out_path / "standard_baseline_panel.csv",
        claim_mode=claim_mode,
        dataset_provenance=dataset_provenance,
    )
    model_card = write_model_card(
        out_path=out_path,
        model_type=model_type,
        selected_model=selected_model_name,
        label_col=label_col,
        feature_set=feature_set,
        features=features,
        split_mode=split_mode,
        dataset_manifest=dataset_manifest,
        metrics=metrics,
        claim_readiness=readiness,
        pu_mode=pu_mode,
        hpo_manifest=hpo_manifest,
        subgroup_metrics_path=out_path / "model_subgroup_metrics.csv",
    )
    with (out_path / "trained_model.pkl").open("wb") as handle:
        pickle.dump(model, handle)
    artifact_manifest = write_artifact_manifest(out_path)
    tracking_record = track_model_run(
        model_dir=out_path,
        repo_root=repo_root or Path.cwd(),
        params={
            "label_col": label_col,
            "feature_set": feature_set,
            "model_type": model_type,
            "selected_model": selected_model_name,
            "split_mode": split_mode,
            "pu_mode": pu_mode,
            "seed": seed,
            "claim_mode": claim_mode,
            "model_params": dict(model_params or effective_model_params),
            "calibration_method": calibration_name,
            "strict_feature_set": strict_feature_set,
            "group_reweight_cols": group_reweight_cols or [],
            "group_reweight_status": group_reweight_manifest.get("status"),
            "retained_excluded_features": retained_excluded_features,
            "retained_context_not_trained": retained_context_not_trained,
            "split_manifest": str(split_manifest) if split_manifest is not None else None,
            "split_manifest_hash": (locked_split_metadata or {}).get("split_manifest_hash") if locked_split_metadata else None,
            "sampled": bool((dataset_provenance or {}).get("sampled")),
        },
        metrics=metrics,
        dataset_manifest=dataset_manifest,
        artifact_manifest=artifact_manifest,
        claim_readiness=readiness,
    )
    model_run_record = write_model_run_record(
        out_path,
        repo_root=repo_root or Path.cwd(),
        notes=str((dataset_provenance or {}).get("notes", "")),
        extra_config={
            "run_id": (dataset_provenance or {}).get("run_id"),
            "label": label_col,
            "feature_set": feature_set,
            "split": split_mode,
            "pu_strategy": pu_mode,
            "model": model_type,
            "calibration": calibration_name,
            "strict_feature_set": strict_feature_set,
            "group_reweight_cols": group_reweight_cols or [],
            "group_reweight_status": group_reweight_manifest.get("status"),
            "retained_excluded_features": retained_excluded_features,
            "retained_context_not_trained": retained_context_not_trained,
            "seed": seed,
        },
    )
    return {
        "metrics": metrics,
        "n_train": len(train),
        "n_test": len(test),
        "features": features,
        "requested_features": effective_requested_features,
        "missing_features": missing_features,
        "all_missing_features": all_missing_features,
        "retained_excluded_features": retained_excluded_features,
        "retained_context_not_trained": retained_context_not_trained,
        "group_reweighting": group_reweight_manifest,
        "claim_readiness": readiness,
        "dataset_manifest": dataset_manifest,
        "model_card": model_card,
        "artifact_manifest": artifact_manifest,
        "tracking_record": tracking_record,
        "model_run_record": model_run_record,
        "subgroup_metrics_rows": int(len(subgroup_metrics)),
        "preprocessing_artifact": preprocessing_manifest,
        "split_manifest_input": str(split_manifest) if split_manifest is not None else None,
    }
