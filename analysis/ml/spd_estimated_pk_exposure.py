from __future__ import annotations

from datetime import datetime, timezone
import json
import math
import pickle
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.model_run_ledger import LEDGER_COLUMNS
from analysis.ml.split_manifest import dataframe_content_hash, write_split_manifest
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier_core import _fit_design_matrix, _transform_design_matrix


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
    if model_type in {"lightgbm", "lgbm"}:
        try:
            from lightgbm import LGBMRegressor
        except ImportError as exc:
            raise RuntimeError("lightgbm regressor requires the optional lightgbm package") from exc
        return LGBMRegressor(
            random_state=seed,
            n_estimators=200,
            max_depth=3,
            num_leaves=7,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_lambda=1.0,
            min_child_samples=10,
            verbosity=-1,
            n_jobs=4,
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


def _partition_profile(name: str, frame: pd.DataFrame) -> dict[str, Any]:
    labels = frame["_observed_exposure_label"]
    relation = frame["_activity_relation"]
    return {
        "partition": name,
        "n_rows": int(len(frame)),
        "n_positive": int(labels.eq(1).sum()),
        "n_negative": int(labels.eq(0).sum()),
        "n_unknown": int(labels.isna().sum()),
        "positive_rate_labeled": (
            float(labels.dropna().astype(int).mean()) if labels.notna().any() else None
        ),
        "n_exact_ac50": int(relation.eq("=").sum()),
        "n_right_censored_ac50": int(relation.str.startswith(">").sum()),
        "n_drugs": int(frame["drug_id"].nunique()) if "drug_id" in frame.columns else None,
        "n_targets": (
            int(frame["target_id"].nunique()) if "target_id" in frame.columns else None
        ),
        "n_pdbs": int(frame["pdb_id"].nunique()) if "pdb_id" in frame.columns else None,
    }


def _write_exposure_data_profile(
    *,
    all_data: pd.DataFrame,
    usable: pd.DataFrame,
    train: pd.DataFrame,
    test: pd.DataFrame,
    potency_train: pd.DataFrame,
    pk_train: pd.DataFrame,
    pk_test: pd.DataFrame,
    out: Path,
) -> dict[str, Any]:
    partitions = pd.DataFrame(
        [
            _partition_profile("full_table", all_data),
            _partition_profile("model_usable", usable),
            _partition_profile("train", train),
            _partition_profile("test", test),
            _partition_profile("potency_train", potency_train),
        ]
    )
    partitions.to_csv(out / "exposure_data_partitions.csv", index=False)

    family_rows: list[dict[str, Any]] = []
    if "target_family" in usable.columns:
        for family, group in usable.groupby("target_family", dropna=False):
            row = _partition_profile(str(family), group)
            row["target_family"] = family
            family_rows.append(row)
    pd.DataFrame(family_rows).to_csv(
        out / "exposure_target_family_profile.csv",
        index=False,
    )

    source_rows: list[dict[str, Any]] = []
    for source_col in ("label_source", "source_family", "upstream_source"):
        if source_col not in usable.columns:
            continue
        for source, group in usable.groupby(source_col, dropna=False):
            row = _partition_profile(str(source), group)
            row["source_col"] = source_col
            row["source_value"] = source
            source_rows.append(row)
    pd.DataFrame(source_rows).to_csv(
        out / "exposure_source_profile.csv",
        index=False,
    )

    censored = all_data["_activity_relation"].str.startswith(">")
    ac50 = all_data["_observed_ac50_um"]
    free_cmax = all_data["_observed_free_cmax_um"]
    free_values = (
        all_data.loc[free_cmax.gt(0), ["drug_id", "_observed_free_cmax_um"]]
        .drop_duplicates()
        .copy()
    )
    usable_free_values = (
        usable.loc[
            usable["_observed_free_cmax_um"].gt(0),
            ["drug_id", "_observed_free_cmax_um"],
        ]
        .drop_duplicates()
        .copy()
    )
    conflicts = (
        free_values.groupby("drug_id")["_observed_free_cmax_um"].nunique().gt(1)
        if not free_values.empty
        else pd.Series(dtype=bool)
    )
    profile = {
        "partitions": partitions.to_dict("records"),
        "train_to_test_prevalence_ratio": (
            float(train["_observed_exposure_label"].mean())
            / float(test["_observed_exposure_label"].mean())
            if float(test["_observed_exposure_label"].mean()) > 0
            else None
        ),
        "n_censored_bound_ge_10_um_full_table": int(
            (censored & ac50.ge(10.0)).sum()
        ),
        "n_censored_bound_ge_30_um_full_table": int(
            (censored & ac50.ge(30.0)).sum()
        ),
        "n_censored_bound_ge_10_um_model_usable": int(
            (
                usable["_activity_relation"].str.startswith(">")
                & usable["_observed_ac50_um"].ge(10.0)
            ).sum()
        ),
        "n_censored_bound_ge_30_um_model_usable": int(
            (
                usable["_activity_relation"].str.startswith(">")
                & usable["_observed_ac50_um"].ge(30.0)
            ).sum()
        ),
        "n_censored_bound_eq_10_um_model_usable": int(
            (
                usable["_activity_relation"].str.startswith(">")
                & usable["_observed_ac50_um"].eq(10.0)
            ).sum()
        ),
        "n_censored_bound_eq_30_um_model_usable": int(
            (
                usable["_activity_relation"].str.startswith(">")
                & usable["_observed_ac50_um"].eq(30.0)
            ).sum()
        ),
        "n_unknown_with_ac50_and_free_cmax": int(
            (
                all_data["_observed_exposure_label"].isna()
                & ac50.gt(0)
                & free_cmax.gt(0)
            ).sum()
        ),
        "n_missing_free_cmax": int((~free_cmax.gt(0)).sum()),
        "n_distinct_pk_drugs_full_table": int(
            free_values["drug_id"].nunique()
        ),
        "n_distinct_free_cmax_values_full_table": int(
            free_values["_observed_free_cmax_um"].nunique()
        ),
        "n_distinct_pk_drugs_model_usable": int(
            usable_free_values["drug_id"].nunique()
        ),
        "n_distinct_free_cmax_values_model_usable": int(
            usable_free_values["_observed_free_cmax_um"].nunique()
        ),
        "n_drugs_with_conflicting_free_cmax": int(conflicts.sum()),
        "free_cmax_min_um": (
            float(free_values["_observed_free_cmax_um"].min())
            if not free_values.empty
            else None
        ),
        "free_cmax_max_um": (
            float(free_values["_observed_free_cmax_um"].max())
            if not free_values.empty
            else None
        ),
        "n_pk_train_drugs": int(len(pk_train)),
        "n_pk_test_drugs": int(len(pk_test)),
        "n_multi_assay_usable_pairs": (
            int(pd.to_numeric(usable["spd_assay_count"], errors="coerce").gt(1).sum())
            if "spd_assay_count" in usable.columns
            else None
        ),
        "n_multi_assay_exact_pairs": (
            int(
                pd.to_numeric(
                    usable.loc[
                        usable["_activity_relation"].eq("="),
                        "spd_assay_count",
                    ],
                    errors="coerce",
                )
                .gt(1)
                .sum()
            )
            if "spd_assay_count" in usable.columns
            else None
        ),
        "policy_flags": [
            "Right-censored AC50 bounds are excluded from potency fitting.",
            "Free Cmax is a drug-level target repeated across target rows.",
            "The usable population is selected by availability of AC50, free Cmax, and a resolved exposure label.",
            "Multiple-assay aggregation can favor the strongest measured AC50.",
        ],
    }
    path = out / "exposure_data_profile.json"
    path.write_text(
        json.dumps(profile, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return {
        "summary": str(path),
        "partitions": str(out / "exposure_data_partitions.csv"),
        "target_families": str(out / "exposure_target_family_profile.csv"),
        "sources": str(out / "exposure_source_profile.csv"),
    }


def run_spd_estimated_pk_exposure_model(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    label_col: str = "spd_exposure_label",
    potency_feature_set: str = "spd_binding_pair_final_full_no_qed",
    pk_feature_set: str = "ligand_physchem_descriptors_no_qed",
    split_mode: str = "drug_holdout",
    model_type: str = "ridge",
    seed: int = 42,
    residual_sigma_floor: float = 0.25,
    censored_policy: str = "exclude",
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
    data["_activity_relation"] = data.get(
        "spd_activity_relation",
        pd.Series("=", index=data.index),
    ).fillna("").astype(str).str.strip()
    all_data = data.copy()
    usable = (
        data["_observed_ac50_um"].gt(0)
        & data["_observed_free_cmax_um"].gt(0)
        & data["_observed_exposure_label"].notna()
    )
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

    if censored_policy not in {"exclude", "bound"}:
        raise ValueError("censored_policy must be 'exclude' or 'bound'")
    potency_train = train.loc[train["_activity_relation"].eq("=")].copy() if censored_policy == "exclude" else train
    if len(potency_train) < 20:
        raise ValueError(
            "potency regression requires at least 20 training rows under the "
            f"{censored_policy!r} censoring policy; found {len(potency_train)}"
        )

    potency_model = _regressor(model_type, seed)
    x_pot_train, potency_preprocessing = _fit_design_matrix(potency_train, potency_features)
    x_pot_test = _transform_design_matrix(test, potency_preprocessing)
    potency_model.fit(x_pot_train, potency_train["_log10_ac50_um"])
    train_log_ac50_pred = pd.Series(potency_model.predict(x_pot_train), index=potency_train.index)
    test_log_ac50_pred = pd.Series(potency_model.predict(x_pot_test), index=test.index)
    ac50_sigma = float(
        (potency_train["_log10_ac50_um"] - train_log_ac50_pred).std(ddof=1)
    )
    if not math.isfinite(ac50_sigma) or ac50_sigma < residual_sigma_floor:
        ac50_sigma = residual_sigma_floor

    pk_train = _unique_drug_table(train, pk_features)
    pk_test = _unique_drug_table(test, pk_features)
    data_profile_outputs = _write_exposure_data_profile(
        all_data=all_data,
        usable=data,
        train=train,
        test=test,
        potency_train=potency_train,
        pk_train=pk_train,
        pk_test=pk_test,
        out=out,
    )
    if len(pk_train) < 5 or len(pk_test) < 2:
        raise ValueError(f"free-Cmax regressor needs >=5 train and >=2 test drugs; found train={len(pk_train)} test={len(pk_test)}")
    pk_model = _regressor(model_type, seed)
    x_pk_train, pk_preprocessing = _fit_design_matrix(pk_train, pk_features)
    x_pk_test = _transform_design_matrix(pk_test, pk_preprocessing)
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
    exact_test = test["_activity_relation"].eq("=")
    if exact_test.any():
        metrics.update(
            _regression_metrics(
                test.loc[exact_test, "_log10_ac50_um"],
                test_log_ac50_pred.loc[exact_test],
                "ac50_exact",
            )
        )
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
            "n_potency_train_rows": int(len(potency_train)),
            "n_exact_test_rows": int(exact_test.sum()),
            "censored_policy": censored_policy,
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
    pred["observed_ac50_relation"] = test["_activity_relation"].astype(str)
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
    top_n = min(20, len(pred))
    top_precision = float(pred.nlargest(top_n, "estimated_exposure_probability")[label_col].mean())
    prevalence = float(pred[label_col].mean())
    repo_root = Path(__file__).resolve().parents[2]
    git_result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )
    model_run_record = {
        "run_id": out.name,
        "date": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_result.stdout.strip() if git_result.returncode == 0 else "unknown",
        "dataset_version": dataframe_content_hash(df),
        "dataset_path": str(dataset_path),
        "label_used": label_col,
        "feature_set": f"{potency_feature_set} + {pk_feature_set}",
        "excluded_columns": "observed free Cmax and label-definition fields",
        "split_method": split_mode,
        "PU_strategy": "standard_binary",
        "model_type": f"two_stage_{model_type}",
        "selected_model": False,
        "hyperparameters": json.dumps({"censored_policy": censored_policy}, sort_keys=True),
        "calibration_method": "none",
        "number_of_rows": int(len(data)),
        "number_of_positives": int(data["_observed_exposure_label"].eq(1).sum()),
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "AUROC": metrics.get("AUROC"),
        "PR_AUC": metrics.get("AUPRC"),
        "precision_at_K": top_precision,
        "precision_K": top_n,
        "enrichment_at_K": top_precision / prevalence if prevalence else 0.0,
        "enrichment_K": top_n,
        "Brier_score": metrics.get("Brier"),
        "notes": "Explicit predicted AC50 / predicted free-Cmax margin; censored AC50 policy recorded in hyperparameters.",
        "model_dir": str(out),
    }
    (out / "model_run_record.json").write_text(
        json.dumps(model_run_record, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    pd.DataFrame([model_run_record]).reindex(columns=LEDGER_COLUMNS).to_csv(
        out / "model_run_record.csv",
        index=False,
    )
    with (out / "potency_regressor.pkl").open("wb") as fh:
        pickle.dump(potency_model, fh)
    with (out / "free_cmax_regressor.pkl").open("wb") as fh:
        pickle.dump(pk_model, fh)
    (out / "potency_preprocessing.json").write_text(
        json.dumps(potency_preprocessing, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (out / "free_cmax_preprocessing.json").write_text(
        json.dumps(pk_preprocessing, indent=2, sort_keys=True),
        encoding="utf-8",
    )
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
            "model_run_record": str(out / "model_run_record.json"),
            "potency_model": str(out / "potency_regressor.pkl"),
            "free_cmax_model": str(out / "free_cmax_regressor.pkl"),
            "potency_preprocessing": str(out / "potency_preprocessing.json"),
            "free_cmax_preprocessing": str(out / "free_cmax_preprocessing.json"),
            "split_manifest": str(out / "split_manifest"),
            "data_profile": data_profile_outputs,
        },
    }
    if split_mode != "drug_holdout":
        manifest["warning"] = "Use drug_holdout for claim-grade estimated-PK exposure results; free Cmax repeats across target rows."
    (out / "estimated_pk_exposure_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
