from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics


DEFAULT_PU_STRATA = [
    "target_family",
    "protein_class",
    "source_family",
    "label_source",
    "upstream_source",
    "assay_type",
    "endpoint_type",
    "activity_type",
    "ligand_chemotype",
    "chemical_cluster",
    "scaffold_key",
    "drug_status",
    "drug_class",
    "approval_status",
    "safety_panel_target",
    "source_objective",
]


def _coerce_probability(values: list[float] | pd.Series) -> pd.Series:
    return pd.to_numeric(pd.Series(values), errors="coerce").clip(lower=0.0, upper=1.0)


def _estimate_c_from_positive_scores(scores: pd.Series, *, floor: float = 0.05) -> dict[str, Any]:
    valid = pd.to_numeric(scores, errors="coerce").dropna().clip(lower=0.0, upper=1.0)
    if valid.empty:
        return {
            "status": "insufficient_positive_scores",
            "c_estimate": None,
            "n_positive_scores": 0,
        }
    c_value = max(float(valid.mean()), floor)
    return {
        "status": "ok",
        "c_estimate": c_value,
        "n_positive_scores": int(len(valid)),
        "positive_score_mean": float(valid.mean()),
        "positive_score_median": float(valid.median()),
        "positive_score_q10": float(valid.quantile(0.10)),
        "positive_score_q90": float(valid.quantile(0.90)),
        "floor": floor,
    }


def correct_elkan_noto_probabilities(
    raw_probs: list[float],
    *,
    c_estimate: float | None,
) -> list[float]:
    raw = _coerce_probability(raw_probs)
    if c_estimate is None or pd.isna(c_estimate) or float(c_estimate) <= 0:
        return raw.astype(float).tolist()
    return (raw / float(c_estimate)).clip(lower=0.0, upper=1.0).astype(float).tolist()


def estimate_elkan_noto_c(
    *,
    positive_score_frame: pd.DataFrame,
    score_col: str,
    strata_cols: list[str] | None = None,
    out_path: str | Path | None = None,
    min_positives_per_stratum: int = 5,
    variation_ratio_warn: float = 2.0,
) -> dict[str, Any]:
    strata_cols = strata_cols or DEFAULT_PU_STRATA
    scores = pd.to_numeric(positive_score_frame.get(score_col), errors="coerce")
    global_estimate = _estimate_c_from_positive_scores(scores)
    rows: list[dict[str, Any]] = []
    available_strata = [col for col in strata_cols if col in positive_score_frame.columns]
    for col in available_strata:
        for value, group in positive_score_frame.groupby(col, dropna=False):
            estimate = _estimate_c_from_positive_scores(group[score_col])
            if estimate.get("status") != "ok":
                continue
            rows.append(
                {
                    "stratum_col": col,
                    "stratum_value": str(value),
                    **estimate,
                    "usable_for_group_correction": bool(
                        int(estimate.get("n_positive_scores") or 0) >= min_positives_per_stratum
                    ),
                }
            )
    strata = pd.DataFrame(rows)
    usable = strata.loc[strata.get("usable_for_group_correction", pd.Series(dtype=bool)).astype(bool)].copy()
    c_values = pd.to_numeric(usable.get("c_estimate"), errors="coerce").dropna()
    scar_warning = None
    variation_ratio = None
    if len(c_values) >= 2 and float(c_values.min()) > 0:
        variation_ratio = float(c_values.max() / c_values.min())
        if variation_ratio >= variation_ratio_warn:
            scar_warning = (
                "Elkan-Noto c varies strongly across strata; SCAR is likely violated. "
                "Treat global correction as a sensitivity baseline and prefer bagged or SAR-aware PU diagnostics."
            )
    manifest: dict[str, Any] = {
        "method": "elkan_noto",
        "status": global_estimate.get("status"),
        "global_c": global_estimate.get("c_estimate"),
        "global_estimate": global_estimate,
        "strata_columns_evaluated": available_strata,
        "n_strata_estimates": int(len(strata)),
        "n_usable_strata_estimates": int(len(usable)),
        "min_positives_per_stratum": int(min_positives_per_stratum),
        "c_variation_ratio": variation_ratio,
        "scar_warning": scar_warning,
        "assumption": (
            "Elkan-Noto correction assumes labeled positives are selected completely at random "
            "from all positives. Atlas reports strata-specific c to diagnose likely SAR/SNAR bias."
        ),
    }
    if out_path is not None:
        out = Path(out_path)
        out.mkdir(parents=True, exist_ok=True)
        strata.to_csv(out / "pu_elkan_noto_c_by_stratum.csv", index=False)
        (out / "pu_elkan_noto_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    return manifest


def write_pu_selection_diagnostics(
    *,
    train_frame: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    label_col: str,
    out_path: str | Path,
    strata_cols: list[str] | None = None,
) -> dict[str, Any]:
    strata_cols = strata_cols or DEFAULT_PU_STRATA
    out = Path(out_path)
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    positive = train_frame[pd.to_numeric(train_frame[label_col], errors="coerce").eq(1)].copy()
    measured_negative = train_frame[pd.to_numeric(train_frame[label_col], errors="coerce").eq(0)].copy()
    pools = [
        ("known_positive", positive),
        ("measured_negative", measured_negative),
        ("unlabeled_pool_fold_local", unlabeled_pool),
    ]
    for col in [candidate for candidate in strata_cols if candidate in train_frame.columns or candidate in unlabeled_pool.columns]:
        for role, frame in pools:
            if col not in frame.columns:
                continue
            counts = frame[col].fillna("missing").astype(str).value_counts(dropna=False)
            for value, count in counts.items():
                rows.append(
                    {
                        "stratum_col": col,
                        "stratum_value": value,
                        "role": role,
                        "n": int(count),
                    }
                )
    table = pd.DataFrame(rows)
    table.to_csv(out / "pu_selection_strata_balance.csv", index=False)
    manifest = {
        "status": "ok",
        "n_train_positive": int(len(positive)),
        "n_train_measured_negative": int(len(measured_negative)),
        "n_unlabeled_pool_fold_local": int(len(unlabeled_pool)),
        "strata_columns_evaluated": sorted(set(table["stratum_col"].tolist())) if not table.empty else [],
        "policy": "Measured inactives and unlabeled examples are counted separately for PU diagnostics.",
    }
    (out / "pu_selection_diagnostics_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest


def write_sar_selection_diagnostics(
    *,
    positive_frame: pd.DataFrame,
    unlabeled_pool: pd.DataFrame,
    out_path: str | Path,
    strata_cols: list[str] | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    strata_cols = strata_cols or DEFAULT_PU_STRATA
    out = Path(out_path)
    out.mkdir(parents=True, exist_ok=True)
    available = [col for col in strata_cols if col in positive_frame.columns or col in unlabeled_pool.columns]
    manifest: dict[str, Any]
    if not available or positive_frame.empty or unlabeled_pool.empty:
        manifest = {
            "status": "skipped",
            "reason": "missing_strata_or_empty_positive_unlabeled_pool",
            "strata_columns_evaluated": available,
        }
        (out / "pu_sar_selection_diagnostics.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return manifest
    pos = positive_frame.reindex(columns=available).copy()
    unlabeled = unlabeled_pool.reindex(columns=available).copy()
    pos["_observed_positive"] = 1
    unlabeled["_observed_positive"] = 0
    frame = pd.concat([pos, unlabeled], ignore_index=True)
    if frame["_observed_positive"].nunique() < 2 or len(frame) < 20:
        manifest = {
            "status": "skipped",
            "reason": "insufficient_rows_or_classes",
            "n_rows": int(len(frame)),
            "n_positive": int(frame["_observed_positive"].sum()),
            "n_unlabeled": int((frame["_observed_positive"] == 0).sum()),
            "strata_columns_evaluated": available,
        }
        (out / "pu_sar_selection_diagnostics.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return manifest
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import train_test_split
    except ImportError:
        manifest = {
            "status": "skipped",
            "reason": "sklearn_unavailable",
            "strata_columns_evaluated": available,
        }
        (out / "pu_sar_selection_diagnostics.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return manifest
    x = pd.get_dummies(frame[available].fillna("missing").astype(str), dummy_na=True).astype(float)
    y = frame["_observed_positive"].astype(int)
    try:
        x_train, x_test, y_train, y_test = train_test_split(
            x,
            y,
            test_size=0.3,
            random_state=seed,
            stratify=y,
        )
    except ValueError:
        x_train, x_test, y_train, y_test = train_test_split(
            x,
            y,
            test_size=0.3,
            random_state=seed,
        )
    model = LogisticRegression(max_iter=1000, class_weight="balanced")
    model.fit(x_train, y_train)
    scores = model.predict_proba(x_test)[:, 1]
    metrics = calibration_metrics([float(v) for v in scores], y_test.astype(int).tolist())
    auroc = metrics.get("AUROC")
    sar_warning = None
    if auroc is not None and pd.notna(auroc) and float(auroc) >= 0.70:
        sar_warning = (
            "Observed positives are separable from unlabeled examples using provenance/strata fields; "
            "SCAR is likely violated and SAR-aware or bagged PU should be preferred."
        )
    manifest = {
        "status": "ok",
        "n_rows": int(len(frame)),
        "n_positive": int(y.sum()),
        "n_unlabeled": int((y == 0).sum()),
        "strata_columns_evaluated": available,
        "metrics": metrics,
        "sar_warning": sar_warning,
        "interpretation": "This model predicts whether a true/observed positive became labeled using only strata/provenance fields.",
    }
    (out / "pu_sar_selection_diagnostics.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
