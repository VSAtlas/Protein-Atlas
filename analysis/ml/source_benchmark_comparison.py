from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.source_benchmark_tables import build_source_benchmark_tables
from analysis.ml.source_calibration import run_source_calibration
from analysis.ml.source_transfer_audit import audit_source_transfer
from analysis.ml.train_classifier import train_ml_model

SPD_LABEL_DEFINING_FEATURES = ["free_cmax_um", "free_cmax_uM", "cmax_um", "fraction_unbound_plasma"]


def _metric_rows(model_dir: Path, analysis: str) -> list[dict[str, object]]:
    metrics_path = model_dir / "model_metrics.csv"
    if not metrics_path.exists():
        return []
    metrics = pd.read_csv(metrics_path)
    rows = []
    for row in metrics.to_dict(orient="records"):
        rows.append(
            {
                "analysis": analysis,
                "metric": row.get("metric"),
                "value": row.get("value"),
                "model_dir": str(model_dir),
            }
        )
    return rows


def _prediction_context_rows(model_dir: Path, analysis: str, label_col: str) -> list[dict[str, object]]:
    pred_path = model_dir / "model_predictions.csv"
    if not pred_path.exists():
        return []
    pred = pd.read_csv(pred_path, low_memory=False)
    if label_col not in pred.columns or "ml_prediction_score" not in pred.columns:
        return []
    labels = pd.to_numeric(pred[label_col], errors="coerce")
    scores = pd.to_numeric(pred["ml_prediction_score"], errors="coerce")
    valid = labels.notna() & scores.notna()
    labels = labels.loc[valid].astype(int)
    scores = scores.loc[valid]
    if labels.empty:
        return []
    prevalence = float(labels.mean())
    mean_prediction = float(scores.mean())
    median_prediction = float(scores.median())
    rows: list[dict[str, object]] = [
        {"analysis": analysis, "metric": "test_positive_rate", "value": prevalence, "model_dir": str(model_dir)},
        {"analysis": analysis, "metric": "test_rows", "value": int(labels.shape[0]), "model_dir": str(model_dir)},
        {"analysis": analysis, "metric": "mean_prediction", "value": mean_prediction, "model_dir": str(model_dir)},
        {"analysis": analysis, "metric": "median_prediction", "value": median_prediction, "model_dir": str(model_dir)},
    ]
    metric_path = model_dir / "model_metrics.csv"
    if metric_path.exists() and prevalence > 0:
        metrics = pd.read_csv(metric_path)
        auprc = pd.to_numeric(metrics.loc[metrics["metric"].eq("AUPRC"), "value"], errors="coerce")
        if not auprc.empty and pd.notna(auprc.iloc[0]):
            rows.append(
                {
                    "analysis": analysis,
                    "metric": "AUPRC_over_prevalence",
                    "value": float(auprc.iloc[0]) / prevalence,
                    "model_dir": str(model_dir),
                }
            )
    return rows


def _calibration_metric_rows(model_dir: Path, analysis: str) -> list[dict[str, object]]:
    metrics_path = model_dir / "source_calibration_metrics.csv"
    if not metrics_path.exists():
        return []
    metrics = pd.read_csv(metrics_path)
    rows = []
    for row in metrics.to_dict(orient="records"):
        rows.append(
            {
                "analysis": f"{analysis}:{row.get('method')}",
                "metric": row.get("metric"),
                "value": row.get("value"),
                "model_dir": str(model_dir),
            }
        )
    return rows


def _calibration_context_rows(model_dir: Path, analysis: str, label_col: str) -> list[dict[str, object]]:
    pred_path = model_dir / "source_calibrated_predictions.csv"
    if not pred_path.exists():
        return []
    pred = pd.read_csv(pred_path, low_memory=False)
    if label_col not in pred.columns:
        return []
    labels = pd.to_numeric(pred[label_col], errors="coerce").dropna()
    if labels.empty:
        return []
    rows: list[dict[str, object]] = [
        {
            "analysis": f"{analysis}:raw_transfer",
            "metric": "test_positive_rate",
            "value": float(labels.mean()),
            "model_dir": str(model_dir),
        },
        {
            "analysis": f"{analysis}:toxcast_hts_logistic_calibrated",
            "metric": "test_positive_rate",
            "value": float(labels.mean()),
            "model_dir": str(model_dir),
        },
    ]
    metrics_path = model_dir / "source_calibration_metrics.csv"
    if metrics_path.exists() and labels.mean() > 0:
        metrics = pd.read_csv(metrics_path)
        for method, group in metrics.groupby("method"):
            auprc = pd.to_numeric(group.loc[group["metric"].eq("AUPRC"), "value"], errors="coerce")
            if not auprc.empty and pd.notna(auprc.iloc[0]):
                rows.append(
                    {
                        "analysis": f"{analysis}:{method}",
                        "metric": "AUPRC_over_prevalence",
                        "value": float(auprc.iloc[0]) / float(labels.mean()),
                        "model_dir": str(model_dir),
                    }
                )
    return rows


def _copy_label_balance(table_path: Path, label_col: str, out_path: Path) -> dict[str, object]:
    if not table_path.exists():
        return {"path": str(table_path), "status": "missing"}
    df = pd.read_csv(table_path, low_memory=False)
    labels = pd.to_numeric(df.get(label_col), errors="coerce")
    summary = {
        "path": str(table_path),
        "rows": int(len(df)),
        "label_col": label_col,
        "labelable_rows": int(labels.notna().sum()),
        "positive_rows": int(labels.eq(1).sum()),
        "negative_rows": int(labels.eq(0).sum()),
        "positive_rate": float(labels.dropna().mean()) if labels.notna().any() else None,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(out_path, index=False)
    return summary


def _write_metric_pivot(rows: list[dict[str, object]], out_path: Path) -> pd.DataFrame:
    metrics = pd.DataFrame(rows)
    if metrics.empty:
        metrics.to_csv(out_path, index=False)
        return metrics
    metrics.to_csv(out_path, index=False)
    pivot = metrics.pivot_table(index="analysis", columns="metric", values="value", aggfunc="first").reset_index()
    pivot.to_csv(out_path.with_name("ml_source_benchmark_metric_matrix.csv"), index=False)
    return pivot


def _append_run_metrics(
    metric_rows: list[dict[str, object]],
    run_dir: Path,
    analysis: str,
    label_col: str,
) -> None:
    metric_rows.extend(_metric_rows(run_dir, analysis))
    metric_rows.extend(_prediction_context_rows(run_dir, analysis, label_col))


def _write_calibration_plot(metrics_dir: Path, run_dirs: dict[str, Path]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    wrote_any = False
    for name, run_dir in run_dirs.items():
        table_path = run_dir / "model_reliability_table.csv"
        if not table_path.exists():
            continue
        table = pd.read_csv(table_path)
        if {"mean_predicted_probability", "observed_positive_rate"}.issubset(table.columns):
            ax.plot(
                table["mean_predicted_probability"],
                table["observed_positive_rate"],
                marker="o",
                linewidth=1.2,
                label=name,
            )
            wrote_any = True
    if not wrote_any:
        plt.close(fig)
        return
    ax.plot([0, 1], [0, 1], color="0.5", linestyle="--", linewidth=1)
    ax.set_xlabel("Mean predicted probability")
    ax.set_ylabel("Observed positive rate")
    ax.set_title("Source-stratified ML calibration")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(metrics_dir / "source_stratified_calibration.png", dpi=160)
    plt.close(fig)


def run_source_benchmark_comparison(
    bioactivity_source_path: str | Path,
    out_dir: str | Path,
    *,
    spd_table_path: str | Path | None = None,
    feature_set: str = "pilot_binding_only",
    model_type: str = "logistic_regression",
    split_mode: str = "drug_holdout",
    seed: int = 42,
    n_bootstraps: int = 100,
    class_weight: str | None = "balanced",
) -> dict[str, Any]:
    out = Path(out_dir)
    tables_dir = out / "tables"
    models_dir = out / "models"
    audit_dir = out / "audits"
    out.mkdir(parents=True, exist_ok=True)
    tables = build_source_benchmark_tables(bioactivity_source_path, tables_dir, spd_table_path=spd_table_path)
    curated_path = tables_dir / "ml_curated_bioactivity_table.csv"
    toxcast_path = tables_dir / "ml_toxcast_hts_table.csv"
    spd_path = tables_dir / "ml_spd_exposure_table.csv"
    transfer_path = tables_dir / "ml_source_transfer_table.csv"
    matched_path = tables_dir / "ml_matched_source_transfer_table.csv"

    balances = {
        "curated": _copy_label_balance(curated_path, "curated_bioactivity_label", out / "ml_curated_bioactivity_label_balance.csv"),
        "toxcast": _copy_label_balance(toxcast_path, "toxcast_hts_label", out / "ml_toxcast_hts_label_balance.csv"),
        "spd": _copy_label_balance(spd_path, "spd_exposure_label", out / "ml_spd_exposure_label_balance.csv"),
    }
    metric_rows: list[dict[str, object]] = []
    run_dirs: dict[str, Path] = {}

    train_ml_model(
        curated_path,
        "curated_bioactivity_label",
        feature_set,
        model_type,
        split_mode,
        models_dir / "within_curated_bioactivity",
        seed,
        n_bootstraps=n_bootstraps,
        class_weight=class_weight,
    )
    run_dirs["within_curated"] = models_dir / "within_curated_bioactivity"
    _append_run_metrics(metric_rows, run_dirs["within_curated"], "within_source_curated_bioactivity", "curated_bioactivity_label")

    train_ml_model(
        toxcast_path,
        "toxcast_hts_label",
        feature_set,
        model_type,
        split_mode,
        models_dir / "within_toxcast_hts",
        seed,
        n_bootstraps=n_bootstraps,
        class_weight=class_weight,
    )
    run_dirs["within_toxcast"] = models_dir / "within_toxcast_hts"
    _append_run_metrics(metric_rows, run_dirs["within_toxcast"], "within_source_toxcast_hts", "toxcast_hts_label")

    if spd_path.exists():
        try:
            train_ml_model(
                spd_path,
                "spd_exposure_label",
                feature_set,
                model_type,
                split_mode,
                models_dir / "within_spd_exposure",
                seed,
                n_bootstraps=n_bootstraps,
                class_weight=class_weight,
                exclude_features=SPD_LABEL_DEFINING_FEATURES,
            )
            run_dirs["within_spd"] = models_dir / "within_spd_exposure"
            _append_run_metrics(metric_rows, run_dirs["within_spd"], "within_source_spd_exposure", "spd_exposure_label")
        except ValueError as exc:
            metric_rows.append({"analysis": "within_source_spd_exposure", "metric": "skipped", "value": str(exc), "model_dir": ""})

    raw_transfer_dir = models_dir / "curated_to_toxcast_raw_transfer"
    train_ml_model(
        transfer_path,
        "source_specific_activity_label",
        feature_set,
        model_type,
        split_mode,
        raw_transfer_dir,
        seed,
        split_column="source_objective",
        train_values=["curated_bioactivity"],
        test_values=["toxcast_hts"],
        n_bootstraps=n_bootstraps,
        class_weight=class_weight,
    )
    run_dirs["raw_transfer"] = raw_transfer_dir
    _append_run_metrics(metric_rows, raw_transfer_dir, "curated_to_toxcast_raw_transfer", "source_specific_activity_label")
    audit_source_transfer(
        transfer_path,
        "source_specific_activity_label",
        audit_dir / "curated_to_toxcast_raw_transfer",
        train_sources=["curated_bioactivity"],
        test_sources=["toxcast_hts"],
        predictions_path=raw_transfer_dir / "model_predictions.csv",
        feature_cols=["atlas_score", "consensus_score", "free_cmax_um", "cmax_um", "fraction_unbound_plasma"],
    )

    calibration_dir = models_dir / "curated_to_toxcast_source_calibrated"
    run_source_calibration(
        transfer_path,
        "source_specific_activity_label",
        feature_set,
        model_type,
        ["curated_bioactivity"],
        "toxcast_hts",
        calibration_dir,
        seed=seed,
        class_weight=class_weight,
    )
    metric_rows.extend(_calibration_metric_rows(calibration_dir, "curated_to_toxcast_source_calibration"))
    metric_rows.extend(
        _calibration_context_rows(
            calibration_dir,
            "curated_to_toxcast_source_calibration",
            "source_specific_activity_label",
        )
    )

    matched_dir = models_dir / "curated_to_toxcast_matched_transfer"
    if matched_path.exists():
        train_ml_model(
            matched_path,
            "source_specific_activity_label",
            feature_set,
            model_type,
            split_mode,
            matched_dir,
            seed,
            split_column="source_objective",
            train_values=["curated_bioactivity"],
            test_values=["toxcast_hts"],
            n_bootstraps=n_bootstraps,
            class_weight=class_weight,
        )
        run_dirs["matched_transfer"] = matched_dir
        _append_run_metrics(metric_rows, matched_dir, "curated_to_toxcast_matched_transfer", "source_specific_activity_label")
        audit_source_transfer(
            matched_path,
            "source_specific_activity_label",
            audit_dir / "curated_to_toxcast_matched_transfer",
            train_sources=["curated_bioactivity"],
            test_sources=["toxcast_hts"],
            predictions_path=matched_dir / "model_predictions.csv",
            feature_cols=["atlas_score", "consensus_score", "free_cmax_um", "cmax_um", "fraction_unbound_plasma"],
        )

    metric_matrix = _write_metric_pivot(metric_rows, out / "ml_source_benchmark_comparison.csv")
    _write_calibration_plot(out, run_dirs)
    manifest: dict[str, Any] = {
        "bioactivity_source_path": str(bioactivity_source_path),
        "spd_table_path": str(spd_table_path) if spd_table_path else "",
        "feature_set": feature_set,
        "model_type": model_type,
        "split_mode": split_mode,
        "table_manifest": tables,
        "label_balances": balances,
        "metric_matrix_rows": int(len(metric_matrix)),
        "interpretation": [
            "Curated ChEMBL/Papyrus, ToxCast HTS, and SPD exposure relevance are separate objectives.",
            "Raw curated-to-ToxCast transfer is a domain-transfer stress test, not a direct publishable performance claim.",
            "Source-calibrated transfer uses held-out ToxCast calibration data and must be reported separately.",
            "Matched-domain transfer narrows the comparison to overlapping target/chemotype/score ranges.",
        ],
    }
    (out / "ml_source_benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
