from __future__ import annotations

from collections.abc import Mapping
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics

def _coerce_int(value: object, default: int) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        return int(value)
    return default


def _coerce_float(value: object, default: float) -> float:
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, int | float):
        return float(value)
    if isinstance(value, str):
        return float(value)
    return default


def _optional_int(value: object, default: int) -> int | None:
    if value is None or (isinstance(value, str) and value.lower() == "none"):
        return None
    return _coerce_int(value, default)


def _max_features(value: object) -> object:
    if value is None or (isinstance(value, str) and value.lower() == "none"):
        return None
    return value


def _model(
    model_type: str,
    seed: int,
    class_weight: str | None = "balanced",
    model_params: Mapping[str, object] | None = None,
) -> Any:
    params = dict(model_params or {})
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression, SGDClassifier
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.svm import SVC

    if model_type == "logistic_regression":
        return LogisticRegression(
            max_iter=_coerce_int(params.get("max_iter"), 5000),
            class_weight=class_weight,
            C=_coerce_float(params.get("C"), 1.0),
            n_jobs=_coerce_int(params.get("n_jobs"), 1),
        )
    if model_type in {"sgd_logistic", "fast_logistic"}:
        return SGDClassifier(
            loss="log_loss",
            max_iter=_coerce_int(params.get("max_iter"), 1000),
            tol=_coerce_float(params.get("tol"), 1e-3),
            alpha=_coerce_float(params.get("alpha"), 1e-4),
            penalty=str(params.get("penalty", "l2")),
            class_weight=class_weight,
            random_state=seed,
            n_jobs=_coerce_int(params.get("n_jobs"), 1),
        )
    if model_type in {"elastic_net", "elastic_net_logistic", "elastic_net_logistic_regression", "logistic_elastic_net"}:
        return LogisticRegression(
            max_iter=_coerce_int(params.get("max_iter"), 5000),
            class_weight=class_weight,
            C=_coerce_float(params.get("C"), 1.0),
            penalty="elasticnet",
            solver="saga",
            l1_ratio=_coerce_float(params.get("l1_ratio"), 0.5),
            random_state=seed,
            n_jobs=_coerce_int(params.get("n_jobs"), 1),
        )
    if model_type == "random_forest":
        max_depth = params.get("max_depth")
        max_features = params.get("max_features", "sqrt")
        return RandomForestClassifier(
            n_estimators=_coerce_int(params.get("n_estimators"), 200),
            random_state=seed,
            class_weight=class_weight,
            max_depth=_optional_int(max_depth, 3),
            min_samples_leaf=_coerce_int(params.get("min_samples_leaf"), 1),
            min_samples_split=_coerce_int(params.get("min_samples_split"), 2),
            max_features=_max_features(max_features),
            n_jobs=_coerce_int(params.get("n_jobs"), 4),
        )
    if model_type in {"svm", "svc"}:
        return SVC(
            probability=True,
            random_state=seed,
            class_weight=class_weight,
            C=_coerce_float(params.get("C"), 1.0),
            kernel=str(params.get("kernel", "rbf")),
            gamma=str(params.get("gamma", "scale")),
        )
    if model_type in {"knn", "k_nearest_neighbors"}:
        return KNeighborsClassifier(
            n_neighbors=_coerce_int(params.get("n_neighbors"), 5),
            weights=str(params.get("weights", "distance")),
            metric=str(params.get("metric", "minkowski")),
            n_jobs=_coerce_int(params.get("n_jobs"), 4),
        )
    if model_type in {"xgboost", "gradient_boosted_trees"}:
        try:
            from xgboost import XGBClassifier

            return XGBClassifier(
                random_state=seed,
                eval_metric="logloss",
                n_estimators=_coerce_int(params.get("n_estimators"), 200),
                max_depth=_coerce_int(params.get("max_depth"), 3),
                learning_rate=_coerce_float(params.get("learning_rate"), 0.1),
                subsample=_coerce_float(params.get("subsample"), 1.0),
                colsample_bytree=_coerce_float(params.get("colsample_bytree"), 1.0),
                reg_lambda=_coerce_float(params.get("reg_lambda"), 1.0),
                min_child_weight=_coerce_float(params.get("min_child_weight"), 1.0),
                n_jobs=_coerce_int(params.get("n_jobs"), 4),
            )
        except ImportError:
            return GradientBoostingClassifier(
                random_state=seed,
                n_estimators=_coerce_int(params.get("n_estimators"), 100),
                max_depth=_coerce_int(params.get("max_depth"), 3),
                learning_rate=_coerce_float(params.get("learning_rate"), 0.1),
                subsample=_coerce_float(params.get("subsample"), 1.0),
            )
    if model_type in {"lightgbm", "shallow_lightgbm", "lgbm"}:
        try:
            from lightgbm import LGBMClassifier
        except ImportError as exc:
            raise RuntimeError("lightgbm model_type requires the optional lightgbm package") from exc
        return LGBMClassifier(
            random_state=seed,
            class_weight=class_weight,
            n_estimators=_coerce_int(params.get("n_estimators"), 200),
            max_depth=_coerce_int(params.get("max_depth"), 3),
            num_leaves=_coerce_int(params.get("num_leaves"), 7),
            learning_rate=_coerce_float(params.get("learning_rate"), 0.05),
            subsample=_coerce_float(params.get("subsample"), 0.8),
            colsample_bytree=_coerce_float(params.get("colsample_bytree"), 0.8),
            reg_lambda=_coerce_float(params.get("reg_lambda"), 1.0),
            min_child_samples=_coerce_int(params.get("min_child_samples"), 10),
            verbosity=-1,
            n_jobs=_coerce_int(params.get("n_jobs"), 4),
        )
    if model_type in {"catboost", "catboost_shallow"}:
        try:
            from catboost import CatBoostClassifier
        except ImportError as exc:
            raise RuntimeError("catboost model_type requires the optional catboost package") from exc
        auto_class_weights = "Balanced" if class_weight == "balanced" else None
        return CatBoostClassifier(
            iterations=_coerce_int(params.get("iterations"), 300),
            depth=_coerce_int(params.get("depth"), 3),
            learning_rate=_coerce_float(params.get("learning_rate"), 0.05),
            l2_leaf_reg=_coerce_float(params.get("l2_leaf_reg"), 3.0),
            loss_function="Logloss",
            eval_metric="AUC",
            random_seed=seed,
            verbose=False,
            allow_writing_files=False,
            auto_class_weights=auto_class_weights,
            thread_count=_coerce_int(params.get("thread_count", params.get("n_jobs")), 4),
        )
    if model_type in {"explainable_boosting_machine", "ebm", "interpret_ebm"}:
        try:
            from interpret.glassbox import ExplainableBoostingClassifier
        except ImportError as exc:
            raise RuntimeError("EBM model_type requires the optional interpret package") from exc
        return ExplainableBoostingClassifier(
            random_state=seed,
            max_rounds=_coerce_int(params.get("max_rounds"), 5000),
            learning_rate=_coerce_float(params.get("learning_rate"), 0.01),
            interactions=_coerce_int(params.get("interactions"), 10),
            max_bins=_coerce_int(params.get("max_bins"), 256),
            n_jobs=_coerce_int(params.get("n_jobs"), 4),
        )
    raise ValueError(f"unsupported model_type: {model_type}")


def _missing_raw_features(df: pd.DataFrame, features: list[str]) -> list[str]:
    return [feature for feature in features if feature not in df.columns]


def _raw_feature_kinds(df: pd.DataFrame, features: list[str]) -> tuple[list[str], list[str]]:
    numeric_features: list[str] = []
    categorical_features: list[str] = []
    for col in features:
        if df[col].dtype == object:
            categorical_features.append(col)
        else:
            numeric_features.append(col)
    return numeric_features, categorical_features


def _design_matrix(df: pd.DataFrame, features: list[str], *, add_missing_indicators: bool = True) -> pd.DataFrame:
    missing = _missing_raw_features(df, features)
    if missing:
        raise ValueError(f"missing raw feature columns: {', '.join(missing)}")
    x = df[features].copy()
    indicator_cols: dict[str, pd.Series] = {}
    for col in x.columns:
        if x[col].dtype == object:
            x[col] = x[col].fillna("missing")
        else:
            x[col] = pd.to_numeric(x[col], errors="coerce")
            if add_missing_indicators and x[col].isna().any():
                indicator_cols[f"{col}__missing"] = x[col].isna().astype(float)
            median = x[col].median()
            x[col] = x[col].fillna(0.0 if pd.isna(median) else median)
    if indicator_cols:
        x = pd.concat([x, pd.DataFrame(indicator_cols, index=x.index)], axis=1)
    return pd.get_dummies(x, dummy_na=True).astype(float)


def _fit_design_matrix(
    df: pd.DataFrame,
    features: list[str],
    *,
    add_missing_indicators: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    missing = _missing_raw_features(df, features)
    if missing:
        raise ValueError(f"missing raw feature columns: {', '.join(missing)}")
    numeric_features, categorical_features = _raw_feature_kinds(df, features)
    numeric_medians: dict[str, float | None] = {}
    missing_indicator_features: list[str] = []
    categorical_levels: dict[str, list[str]] = {}
    for col in numeric_features:
        values = pd.to_numeric(df[col], errors="coerce")
        median = values.median()
        numeric_medians[col] = None if pd.isna(median) else float(median)
        if add_missing_indicators and values.isna().any():
            missing_indicator_features.append(col)
    for col in categorical_features:
        categorical_levels[col] = sorted(df[col].fillna("missing").astype(str).unique().tolist())
    x = _design_matrix(df, features, add_missing_indicators=add_missing_indicators)
    means = x.mean(axis=0).fillna(0.0)
    stds = x.std(axis=0).replace(0, 1.0).fillna(1.0)
    scaled = ((x - means) / stds).fillna(0.0)
    centroid_distance = (scaled.pow(2).sum(axis=1) ** 0.5) if not scaled.empty else pd.Series(dtype=float)
    metadata: dict[str, Any] = {
        "version": 2,
        "raw_features": list(features),
        "numeric_features": numeric_features,
        "categorical_features": categorical_features,
        "categorical_levels": categorical_levels,
        "numeric_medians": numeric_medians,
        "missing_indicator_features": missing_indicator_features,
        "add_missing_indicators": bool(add_missing_indicators),
        "design_columns": list(x.columns),
        "n_fit_rows": int(len(df)),
        "applicability_domain": {
            "method": "scaled_training_centroid_distance",
            "center": {str(k): float(v) for k, v in means.items()},
            "scale": {str(k): float(v) for k, v in stds.items()},
            "distance_threshold_train_q95": float(centroid_distance.quantile(0.95)) if len(centroid_distance) else None,
            "n_reference_rows": int(len(x)),
        },
    }
    return x, metadata


def _transform_design_matrix(
    df: pd.DataFrame,
    preprocessing: Mapping[str, Any],
    *,
    strict_raw_features: bool = True,
    strict_categories: bool = False,
) -> pd.DataFrame:
    features = [str(feature) for feature in preprocessing.get("raw_features", [])]
    if not features:
        raise ValueError("preprocessing artifact has no raw_features")
    missing = _missing_raw_features(df, features)
    if missing and strict_raw_features:
        raise ValueError(f"missing raw feature columns: {', '.join(missing)}")
    if missing:
        df = df.copy()
        for col in missing:
            df[col] = pd.NA
    numeric_features = [str(col) for col in preprocessing.get("numeric_features", [])]
    categorical_features = [str(col) for col in preprocessing.get("categorical_features", [])]
    categorical_levels = {str(k): {str(v) for v in values} for k, values in dict(preprocessing.get("categorical_levels", {})).items()}
    numeric_medians = dict(preprocessing.get("numeric_medians", {}))
    missing_indicator_features = {str(col) for col in preprocessing.get("missing_indicator_features", [])}
    x = df[features].copy()
    unseen: dict[str, list[str]] = {}
    indicator_cols: dict[str, pd.Series] = {}
    for col in categorical_features:
        if col in x.columns:
            filled = x[col].fillna("missing").astype(str)
            known = categorical_levels.get(col, set())
            extra = sorted(set(filled.unique().tolist()) - known) if known else []
            if extra:
                unseen[col] = extra
                if not strict_categories:
                    filled = filled.where(filled.isin(known), "missing")
            x[col] = filled
    if unseen and strict_categories:
        formatted = "; ".join(f"{col}={values[:5]}" for col, values in unseen.items())
        raise ValueError(f"unseen categorical levels for fitted preprocessing: {formatted}")
    for col in numeric_features:
        if col not in x.columns:
            continue
        values = pd.to_numeric(x[col], errors="coerce")
        if col in missing_indicator_features:
            indicator_cols[f"{col}__missing"] = values.isna().astype(float)
        median = numeric_medians.get(col)
        fill_value = 0.0 if median is None or pd.isna(median) else float(median)
        x[col] = values.fillna(fill_value)
    if indicator_cols:
        x = pd.concat([x, pd.DataFrame(indicator_cols, index=x.index)], axis=1)
    transformed = pd.get_dummies(x, dummy_na=True).astype(float)
    design_columns = [str(col) for col in preprocessing.get("design_columns", [])]
    if not design_columns:
        raise ValueError("preprocessing artifact has no design_columns")
    extra_design = sorted(set(transformed.columns) - set(design_columns))
    if extra_design and strict_categories:
        raise ValueError(f"unexpected transformed design columns: {extra_design[:10]}")
    aligned = pd.DataFrame(index=transformed.index)
    for col in design_columns:
        aligned[col] = transformed[col] if col in transformed.columns else 0.0
    return aligned.astype(float)


def write_preprocessing_artifacts(
    *,
    out_path: str | Path,
    preprocessing: Mapping[str, Any],
    x_train: pd.DataFrame,
) -> dict[str, Any]:
    out = Path(out_path)
    out.mkdir(parents=True, exist_ok=True)
    manifest = dict(preprocessing)
    manifest["training_design_matrix"] = str(out / "training_design_matrix.csv")
    (out / "model_preprocessing.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    x_train.to_csv(out / "training_design_matrix.csv", index=False)
    return manifest


def load_preprocessing_artifact(model_dir: str | Path) -> dict[str, Any]:
    path = Path(model_dir) / "model_preprocessing.json"
    if not path.exists():
        raise ValueError(f"missing preprocessing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _clip_weight(value: object, default: float = 1.0) -> float:
    try:
        parsed = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default
    if pd.isna(parsed):
        return default
    return min(1.0, max(0.05, parsed))


def _sample_weights(frame: pd.DataFrame, label_col: str, *, pu_mode: str) -> pd.Series | None:
    weights = pd.Series(1.0, index=frame.index)
    used = False
    if pu_mode in {"positive_unlabeled_weighted", "case_control_matched"} and "_sample_weight" in frame.columns:
        weights = pd.to_numeric(frame["_sample_weight"], errors="coerce").fillna(1.0).clip(lower=0.05, upper=1.0)
        used = True
    labels = pd.to_numeric(frame[label_col], errors="coerce")
    if "negative_confidence" in frame.columns:
        confidence = frame["negative_confidence"].map(_clip_weight)
        neg_mask = labels.eq(0)
        if neg_mask.any():
            weights.loc[neg_mask] = weights.loc[neg_mask] * confidence.loc[neg_mask]
            used = True
    if "negative_evidence_type" in frame.columns:
        evidence_weight = {
            "fold_specific_reliable_negative": 0.35,
            "reliable_negative": 0.35,
            "pu_temporary_unlabeled_negative": 0.25,
            "faers_nonsignal": 0.25,
            "offsides_nonsignal": 0.25,
            "omop_negative_control": 0.75,
            "measured_inactive": 1.0,
            "measured_or_reliable_negative": 0.75,
        }
        mapped = frame["negative_evidence_type"].fillna("").astype(str).map(evidence_weight).dropna()
        if not mapped.empty:
            weights.loc[mapped.index] = weights.loc[mapped.index] * mapped
            used = True
    return weights if used else None


def group_reweighting_weights(
    frame: pd.DataFrame,
    label_col: str,
    *,
    group_cols: list[str] | None,
    include_label: bool = True,
    max_factor: float = 5.0,
) -> tuple[pd.Series | None, dict[str, object]]:
    """Return capped inverse-frequency weights for source/family strata.

    This is a practical source/family reweighting control: common strata are
    downweighted and rare strata are upweighted. It does not relabel rows or turn
    unknowns into negatives.
    """

    available = [col for col in (group_cols or []) if col in frame.columns]
    manifest: dict[str, object] = {
        "requested_group_cols": list(group_cols or []),
        "available_group_cols": available,
        "include_label": bool(include_label),
        "max_factor": float(max_factor),
        "status": "not_requested" if not group_cols else "skipped",
    }
    if not available:
        manifest["reason"] = "no requested group columns are present"
        return None, manifest
    work = frame[available].astype("object").where(pd.notna(frame[available]), "missing").astype(str)
    if include_label and label_col in frame.columns:
        work["__label__"] = pd.to_numeric(frame[label_col], errors="coerce").fillna("unlabeled").astype(str)
    key = work[available[0]]
    for col in work.columns[1:]:
        key = key.str.cat(work[col], sep="|")
    counts = key.value_counts(dropna=False)
    if counts.empty:
        manifest["reason"] = "empty group counts"
        return None, manifest
    cap = max(1.0, float(max_factor))
    raw = len(frame) / (len(counts) * key.map(counts).astype(float))
    weights = raw.clip(lower=1.0 / cap, upper=cap)
    weights = weights / weights.mean() if weights.mean() else weights
    weights = weights.clip(lower=1.0 / cap, upper=cap)
    manifest.update(
        {
            "status": "ok",
            "n_groups": int(len(counts)),
            "min_group_size": int(counts.min()),
            "max_group_size": int(counts.max()),
            "min_weight": float(weights.min()),
            "max_weight": float(weights.max()),
            "mean_weight": float(weights.mean()),
            "policy": "capped inverse-frequency source/family weights; rows remain in their original label state",
        }
    )
    return weights.astype(float), manifest


def combine_sample_weights(*weights: pd.Series | None) -> pd.Series | None:
    present = [weight for weight in weights if weight is not None]
    if not present:
        return None
    out = pd.Series(1.0, index=present[0].index)
    for weight in present:
        out = out * weight.reindex(out.index).fillna(1.0)
    return out


def _exclude_unlabeled_holdout_entities(unlabeled: pd.DataFrame, test: pd.DataFrame, split_mode: str) -> pd.DataFrame:
    if unlabeled.empty:
        return unlabeled
    fields_by_mode = {
        "drug_holdout": ["drug_id"],
        "target_holdout": ["target_id"],
        "scaffold_holdout": ["scaffold_key", "ligand_chemotype"],
        "chemical_cluster_holdout": ["chemical_cluster", "scaffold_key", "ligand_chemotype"],
        "target_family_holdout": ["target_family", "protein_class"],
        "protein_class_holdout": ["protein_class"],
    }
    out = unlabeled.copy()
    for field in fields_by_mode.get(split_mode, []):
        if field not in out.columns or field not in test.columns:
            continue
        held_out = set(test[field].dropna().astype(str))
        if held_out:
            return out[~out[field].fillna("").astype(str).isin(held_out)].copy()
    return out



def _fit_model_with_optional_weights(
    model: Any,
    x: pd.DataFrame,
    y: pd.Series,
    sample_weight: pd.Series | None = None,
) -> Any:
    if sample_weight is None:
        model.fit(x, y)
        return model
    try:
        model.fit(x, y, sample_weight=sample_weight)
    except TypeError:
        model.fit(x, y)
    return model

def _model_probabilities(model: Any, x: pd.DataFrame) -> list[float]:
    if hasattr(model, "predict_proba"):
        return [float(v) for v in model.predict_proba(x)[:, 1]]
    raw = model.decision_function(x)
    return [float(1 / (1 + pow(2.718281828, -float(v)))) for v in raw]


def _candidate_models(
    model_type: str,
    seed: int,
    class_weight: str | None,
    model_params: Mapping[str, object] | None = None,
) -> list[tuple[str, Any]]:
    params = dict(model_params or {})
    if model_type == "logistic_regression":
        candidates: list[tuple[str, Any]] = []
        for c_value in [0.1, 1.0, 10.0]:
            candidates.append(
                (
                    f"logistic_regression_C{c_value:g}",
                    _model(
                        "logistic_regression",
                        seed,
                        class_weight=class_weight,
                        model_params={**params, "C": c_value},
                    ),
                )
            )
        return candidates
    if model_type in {"elastic_net", "elastic_net_logistic", "elastic_net_logistic_regression", "logistic_elastic_net"}:
        candidates = []
        for c_value in [0.1, 1.0, 10.0]:
            for l1_ratio in [0.2, 0.5, 0.8]:
                name = f"elastic_net_logistic_C{c_value:g}_l1{l1_ratio:g}"
                candidates.append(
                    (
                        name,
                        _model(
                            "elastic_net_logistic",
                            seed,
                            class_weight=class_weight,
                            model_params={**params, "C": c_value, "l1_ratio": l1_ratio},
                        ),
                    )
                )
        return candidates
    return [(model_type, _model(model_type, seed, class_weight=class_weight, model_params=model_params))]


def _select_model(
    model_type: str,
    seed: int,
    class_weight: str | None,
    x_train: pd.DataFrame,
    y_train: pd.Series,
    x_val: pd.DataFrame,
    y_val: pd.Series,
    model_params: Mapping[str, object] | None = None,
) -> tuple[str, Any, pd.DataFrame]:
    if x_val.empty or y_val.nunique() < 2:
        return (
            model_type,
            _model(model_type, seed, class_weight=class_weight, model_params=model_params),
            pd.DataFrame(),
        )
    rows: list[dict[str, object]] = []
    best_name = model_type
    best_model: Any | None = None
    best_score = -1.0
    for name, model in _candidate_models(model_type, seed, class_weight, model_params=model_params):
        model.fit(x_train, y_train)
        probs = model.predict_proba(x_val)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(x_val)
        metrics = calibration_metrics([float(v) for v in probs], y_val.astype(int).tolist())
        score = float(metrics.get("AUPRC", 0.0))
        rows.append({"candidate": name, **metrics})
        if score > best_score:
            best_score = score
            best_name = name
            best_model = model
    return best_name, best_model or _model(model_type, seed, class_weight=class_weight), pd.DataFrame(rows)
