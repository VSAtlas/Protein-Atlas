from __future__ import annotations

import json
import math
import pickle
from pathlib import Path
from typing import Any

import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.split_manifest import write_split_manifest
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier_core import _design_matrix


def _first_numeric(df: pd.DataFrame, names: list[str]) -> pd.Series:
    for name in names:
        if name in df.columns:
            return pd.to_numeric(df[name], errors="coerce")
    return pd.Series(pd.NA, index=df.index, dtype="Float64")


def _normal_cdf(values: pd.Series) -> pd.Series:
    return values.map(lambda z: 0.5 * (1.0 + math.erf(float(z) / math.sqrt(2.0))))


def _regressor(model_type: str, seed: int) -> Any:
    if model_type == "ridge":
        return Ridge(alpha=1.0, random_state=seed)
    if model_type == "random_forest":
        from sklearn.ensemble import RandomForestRegressor

        return RandomForestRegressor(
            n_estimators=300,
            min_samples_leaf=3,
            random_state=seed,
            n_jobs=1,
        )
    raise ValueError(f"unsupported regressor: {model_type}")


def _regression_metrics(observed_log: pd.Series, predicted_log: pd.Series, prefix: str) -> dict[str, float]:
    residual = observed_log - predicted_log
    rmse = math.sqrt(float((residual**2).mean())) if len(residual) else float("nan")
    metrics = {
        f"{prefix}_log10_mae": float(residual.abs().mean()),
        f"{prefix}_log10_rmse": float(rmse),
        f"{prefix}_fold_error_median": float((10 ** residual.abs()).median()),
    }
    if len(observed_log) > 1 and observed_log.nunique(dropna=True) > 1:
        metrics[f"{prefix}_r2"] = float(r2_score(observed_log, predicted_log))
    return metrics


def _unique_drug_table(frame: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    cols = ["drug_id", "_log10_free_cmax_um", *features]
    cols = [col for col in cols if col in frame.columns]
    work = frame[cols].dropna(subset=["drug_id", "_log10_free_cmax_um"]).copy()
    if work.empty:
        return work
    numeric = work.select_dtypes(exclude="object").columns.difference(["_log10_free_cmax_um"])
    categorical = [col for col in work.columns if col not in set(numeric) | {"_log10_free_cmax_um"}]
    agg: dict[str, str] = {col: "median" for col in numeric}
    agg.update({col: "first" for col in categorical})
    agg["_log10_free_cmax_um"] = "median"
    return work.groupby("drug_id", as_index=False).agg(agg)


def run_spd_estimated_pk_exposure_model(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    label_col: str = "spd_exposure_label",
    potency_feature_set: str = "spd_potency_physchem",
    pk_feature_set: str = "ligand_physchem_descriptors",
    split_mode: str = "drug_holdout",
    model_type: str = "ridge",
    seed: int = 42,
    residual_sigma_floor: float = 0.25,
) -> dict[str, Any]:
    """Predict SPD exposure relevance without using observed free Cmax as an input.

    The model has two non-memorizing stages:
    1. pair-level potency regression: non-PK features -> log10(AC50 uM)
    2. drug-level exposure regression: RDKit descriptors -> log10(free Cmax uM)

    The final exposure probability combines both predictions as
    P(predicted AC50 / predicted free Cmax <= 10). Observed free Cmax is a
    training target for the exposure regressor, not a classifier feature.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(dataset_path, low_memory=False)
    ac50 = _first_numeric(df, ["spd_ac50_uM", "ac50_um", "ac50_uM"])
    if ac50.isna().all() and "ac50_nM" in df.columns:
        ac50 = pd.to_numeric(df["ac50_nM"], errors="coerce") / 1000.0
    free_cmax = _first_numeric(df, ["free_cmax_um", "free_cmax_uM", "free_cmax"])
    label = binary_label_series(df[label_col]) if label_col in df.columns else pd.Series(pd.NA, index=df.index)

    data = df.copy()
    data["_observed_ac50_um"] = ac50.astype(float)
    data["_observed_free_cmax_um"] = free_cmax.astype(float)
    data["_observed_exposure_label"] = label
    usable = data["_observed_ac50_um"].gt(0) & data["_observed_free_cmax_um"].gt(0) & data["_observed_exposure_label"].notna()
    data = data.loc[usable].copy()
    if len(data) < 20:
        raise ValueError(f"estimated-PK SPD exposure model requires at least 20 usable rows; found {len(data)}")
    data["_log10_ac50_um"] = data["_observed_ac50_um"].map(math.log10)
    data["_log10_free_cmax_um"] = data["_observed_free_cmax_um"].map(math.log10)

    excluded = effective_exclude_features(label_col, None, allow_label_definition_features=False)
    potency_features = [feature for feature in get_feature_set(potency_feature_set) if feature in data.columns and feature not in excluded]
    pk_features = [feature for feature in get_feature_set(pk_feature_set) if feature in data.columns and feature not in excluded]
    assert_no_leakage(potency_features, allow_label_definition_features=False)
    assert_no_leakage(pk_features, allow_label_definition_features=False)
    if not potency_features:
        raise ValueError(f"potency feature set {potency_feature_set!r} has no usable columns")
    if not pk_features:
        raise ValueError(f"PK feature set {pk_feature_set!r} has no usable descriptor columns")

    train_idx, test_idx = make_split(data, split_mode=split_mode, seed=seed)
    train = data.loc[train_idx].copy()
    test = data.loc[test_idx].copy()
    split_summary = split_overlap_summary(train, test, split_mode)
    if split_summary["overlaps"] and not split_summary["passes_holdout"]:
        raise ValueError(f"holdout split leakage detected: {split_summary}")

    potency_model = _regressor(model_type, seed)
    x_pot_train = _design_matrix(train, potency_features)
    x_pot_test = _design_matrix(test, potency_features).reindex(columns=x_pot_train.columns, fill_value=0)
    potency_model.fit(x_pot_train, train["_log10_ac50_um"])
    train_log_ac50_pred = pd.Series(potency_model.predict(x_pot_train), index=train.index)
    test_log_ac50_pred = pd.Series(potency_model.predict(x_pot_test), index=test.index)
    ac50_sigma = float((train["_log10_ac50_um"] - train_log_ac50_pred).std(ddof=1))
    if not math.isfinite(ac50_sigma) or ac50_sigma < residual_sigma_floor:
        ac50_sigma = residual_sigma_floor

    pk_train = _unique_drug_table(train, pk_features)
    pk_test = _unique_drug_table(test, pk_features)
    if len(pk_train) < 5 or len(pk_test) < 2:
        raise ValueError(f"free-Cmax regressor needs >=5 train and >=2 test drugs; found train={len(pk_train)} test={len(pk_test)}")
    pk_model = _regressor(model_type, seed)
    x_pk_train = _design_matrix(pk_train, pk_features)
    x_pk_test = _design_matrix(pk_test, pk_features).reindex(columns=x_pk_train.columns, fill_value=0)
    pk_model.fit(x_pk_train, pk_train["_log10_free_cmax_um"])
    pk_train_pred = pd.Series(pk_model.predict(x_pk_train), index=pk_train.index)
    pk_test_pred = pd.Series(pk_model.predict(x_pk_test), index=pk_test.index)
    cmax_sigma = float((pk_train["_log10_free_cmax_um"] - pk_train_pred).std(ddof=1))
    if not math.isfinite(cmax_sigma) or cmax_sigma < residual_sigma_floor:
        cmax_sigma = residual_sigma_floor
    pk_pred_by_drug = dict(zip(pk_test["drug_id"].astype(str), pk_test_pred.astype(float), strict=False))
    test_log_cmax_pred = test["drug_id"].astype(str).map(pk_pred_by_drug)
    test = test.loc[test_log_cmax_pred.notna()].copy()
    test_log_ac50_pred = test_log_ac50_pred.loc[test.index]
    test_log_cmax_pred = test_log_cmax_pred.loc[test.index].astype(float)

    sigma = math.sqrt(ac50_sigma**2 + cmax_sigma**2)
    exposure_z = (1.0 + test_log_cmax_pred - test_log_ac50_pred) / sigma
    estimated_exposure_probability = _normal_cdf(exposure_z).clip(1e-9, 1.0 - 1e-9)
    observed_labels = test["_observed_exposure_label"].astype(int).tolist()
    metrics = calibration_metrics(estimated_exposure_probability.astype(float).tolist(), observed_labels)
    metrics.update(_regression_metrics(test["_log10_ac50_um"], test_log_ac50_pred, "ac50"))
    pk_test_observed = pk_test.set_index("drug_id").loc[list(pk_pred_by_drug.keys()), "_log10_free_cmax_um"]
    pk_test_pred_aligned = pd.Series(pk_pred_by_drug).loc[pk_test_observed.index]
    metrics.update(_regression_metrics(pk_test_observed.astype(float), pk_test_pred_aligned.astype(float), "free_cmax"))
    metrics.update(
        {
            "model_type": model_type,
            "potency_feature_set": potency_feature_set,
            "pk_feature_set": pk_feature_set,
            "split_mode": split_mode,
            "n_train_pairs": int(len(train)),
            "n_test_pairs": int(len(test)),
            "n_train_drugs_for_pk": int(len(pk_train)),
            "n_test_drugs_for_pk": int(len(pk_test)),
            "ac50_residual_sigma_log10": ac50_sigma,
            "free_cmax_residual_sigma_log10": cmax_sigma,
            "policy": "descriptor_estimated_free_cmax_plus_non_pk_potency",
        }
    )

    pred_cols = [col for col in ["drug_id", "target_id", "pdb_id", "ligand_chemotype", "scaffold_key", "target_family"] if col in test.columns]
    pred = test[pred_cols].copy()
    pred["observed_ac50_um"] = test["_observed_ac50_um"].astype(float)
    pred["predicted_ac50_um"] = (10 ** test_log_ac50_pred).astype(float)
    pred["observed_free_cmax_um"] = test["_observed_free_cmax_um"].astype(float)
    pred["predicted_free_cmax_um"] = (10 ** test_log_cmax_pred).astype(float)
    pred["estimated_exposure_margin"] = pred["predicted_ac50_um"] / pred["predicted_free_cmax_um"]
    pred["estimated_exposure_probability"] = estimated_exposure_probability.astype(float)
    pred[label_col] = test["_observed_exposure_label"].astype(int)
    pred.to_csv(out / "estimated_pk_exposure_predictions.csv", index=False)
    pd.DataFrame([{"metric": key, "value": value} for key, value in metrics.items()]).to_csv(
        out / "estimated_pk_exposure_metrics.csv",
        index=False,
    )
    with (out / "potency_regressor.pkl").open("wb") as fh:
        pickle.dump(potency_model, fh)
    with (out / "free_cmax_regressor.pkl").open("wb") as fh:
        pickle.dump(pk_model, fh)
    write_split_manifest(
        data,
        train_idx,
        test.index,
        out / "split_manifest",
        label_col=label_col,
        split_mode=split_mode,
        split_summary=split_summary,
    )
    manifest: dict[str, Any] = {
        "dataset": str(dataset_path),
        "out_dir": str(out),
        "label_col": label_col,
        "potency_feature_set": potency_feature_set,
        "pk_feature_set": pk_feature_set,
        "potency_features": potency_features,
        "pk_features": pk_features,
        "model_type": model_type,
        "split_summary": split_summary,
        "metrics": metrics,
        "policy": [
            "Use RDKit descriptors to predict free Cmax at drug level; observed free Cmax is never an input feature.",
            "Use Atlas/chemistry/target features to predict SPD AC50 at pair level.",
            "Combine predicted AC50 and predicted free Cmax to estimate exposure relevance.",
            "Drug-holdout is required for claim-grade free-Cmax evaluation because free Cmax is drug-level.",
        ],
        "outputs": {
            "predictions": str(out / "estimated_pk_exposure_predictions.csv"),
            "metrics": str(out / "estimated_pk_exposure_metrics.csv"),
            "potency_model": str(out / "potency_regressor.pkl"),
            "free_cmax_model": str(out / "free_cmax_regressor.pkl"),
            "split_manifest": str(out / "split_manifest"),
        },
    }
    if split_mode != "drug_holdout":
        manifest["warning"] = "Use drug_holdout for claim-grade estimated-PK exposure results; free Cmax repeats across target rows."
    (out / "estimated_pk_exposure_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
