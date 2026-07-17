from __future__ import annotations

from collections.abc import Sequence
import json
import math
import pickle
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_GROUP_COLUMNS = [
    "target_family",
    "protein_class",
    "source_family",
    "label_source",
    "external_source_family",
    "upstream_source",
    "scaffold_key",
    "chemical_cluster",
    "ligand_chemotype",
]
DEFAULT_TARGET_COLUMNS = ["target_id", "pdb_id"]
DEFAULT_TOP_KS = [10, 20, 50, 100]
DEFAULT_FAILURE_GROUP_COLUMNS = [
    "target_family",
    "protein_class",
    "target_id",
    "pdb_id",
    "scaffold_key",
    "chemical_cluster",
    "ligand_chemotype",
    "label_source",
    "source_family",
    "external_source_family",
]
DEFAULT_SOURCE_COLUMNS = ["label_source", "source_family", "upstream_source", "external_source_family"]
DEFAULT_SCAFFOLD_COLUMNS = ["scaffold_key", "chemical_cluster", "ligand_chemotype"]
ID_OR_PROVENANCE_TOKENS = [
    "drug_id",
    "target_id",
    "pdb_id",
    "scaffold",
    "source",
    "target_family",
    "protein_class",
    "cluster",
    "chemotype",
]


def render_model_figure_suite(
    *,
    model_dir: str | Path | None = None,
    predictions: str | Path | None = None,
    dataset: str | Path | None = None,
    label_col: str | None = None,
    score_col: str = "ml_prediction_score",
    out_dir: str | Path | None = None,
    group_cols: Sequence[str] | None = None,
    target_cols: Sequence[str] | None = None,
    top_ks: Sequence[int] | None = None,
    max_groups: int = 20,
    min_group_size: int = 5,
    shap_sample_rows: int = 250,
    shap_background_rows: int = 100,
    random_state: int = 42,
    skip_shap: bool = False,
    force_generic_shap: bool = False,
    max_generic_shap_features: int = 200,
    threshold: float = 0.5,
    top_n_errors: int = 50,
    title_prefix: str | None = None,
) -> dict[str, Any]:
    """Write publication-oriented ML diagnostic figures and compact findings.

    The primary input is a trained Atlas model directory. Passing only
    ``model_predictions.csv`` still generates all prediction-based figures; SHAP
    additionally needs ``trained_model.pkl`` and ``training_design_matrix.csv``.
    """

    model_path = Path(model_dir).resolve() if model_dir else None
    pred_path = _resolve_predictions(model_path, predictions)
    output_dir = _resolve_out_dir(model_path, out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pred = pd.read_csv(pred_path, low_memory=False)
    label = label_col or detect_binary_label(pred)
    if label not in pred.columns:
        raise ValueError(f"label column not found in predictions: {label}")
    if score_col not in pred.columns:
        raise ValueError(f"score column not found in predictions: {score_col}")

    label_frame = _load_label_frame(dataset, pred, label)
    prefix = title_prefix or _default_title_prefix(model_path, pred_path)
    selected_groups = [col for col in (group_cols or DEFAULT_GROUP_COLUMNS) if col in label_frame.columns]
    selected_targets = [col for col in (target_cols or DEFAULT_TARGET_COLUMNS) if col in pred.columns]
    ks = _normalize_top_ks(top_ks or DEFAULT_TOP_KS, len(pred))

    outputs: dict[str, str] = {}
    skipped: dict[str, str] = {}
    findings: list[str] = []

    _configure_matplotlib()
    label_distribution = write_label_distribution_outputs(
        label_frame,
        label_col=label,
        group_cols=selected_groups,
        out_dir=output_dir,
        max_groups=max_groups,
        title_prefix=prefix,
    )
    outputs.update(label_distribution.outputs)
    skipped.update(label_distribution.skipped)

    target_rates = write_positive_rate_outputs(
        pred,
        label_col=label,
        group_cols=selected_targets,
        out_dir=output_dir,
        max_groups=max_groups,
        min_group_size=min_group_size,
        title_prefix=prefix,
    )
    outputs.update(target_rates.outputs)
    skipped.update(target_rates.skipped)

    pr_outputs = write_precision_recall_outputs(
        pred,
        label_col=label,
        score_col=score_col,
        out_dir=output_dir,
        title_prefix=prefix,
    )
    outputs.update(pr_outputs.outputs)
    skipped.update(pr_outputs.skipped)

    topk_outputs = write_topk_outputs(
        pred,
        label_col=label,
        score_col=score_col,
        top_ks=ks,
        out_dir=output_dir,
        title_prefix=prefix,
    )
    outputs.update(topk_outputs.outputs)
    skipped.update(topk_outputs.skipped)

    calibration_outputs = write_calibration_outputs(
        pred,
        label_col=label,
        score_col=score_col,
        out_dir=output_dir,
        title_prefix=prefix,
    )
    outputs.update(calibration_outputs.outputs)
    skipped.update(calibration_outputs.skipped)

    ad_outputs = write_applicability_domain_outputs(
        pred,
        label_col=label,
        score_col=score_col,
        out_dir=output_dir,
        title_prefix=prefix,
    )
    outputs.update(ad_outputs.outputs)
    skipped.update(ad_outputs.skipped)

    failure_outputs = write_failure_analysis_outputs(
        pred,
        label_col=label,
        score_col=score_col,
        out_dir=output_dir,
        threshold=threshold,
        top_n=top_n_errors,
        title_prefix=prefix,
    )
    outputs.update(failure_outputs.outputs)
    skipped.update(failure_outputs.skipped)

    feature_outputs = write_feature_importance_outputs(
        model_path=model_path,
        out_dir=output_dir,
        title_prefix=prefix,
        max_features=max_groups,
    )
    outputs.update(feature_outputs.outputs)
    skipped.update(feature_outputs.skipped)

    if skip_shap:
        skipped["shap"] = "skipped by --skip-shap"
    else:
        shap_outputs = write_shap_outputs(
            model_path=model_path,
            out_dir=output_dir,
            sample_rows=shap_sample_rows,
            background_rows=shap_background_rows,
            random_state=random_state,
            force_generic=force_generic_shap,
            max_generic_features=max_generic_shap_features,
            max_features=max_groups,
            title_prefix=prefix,
        )
        outputs.update(shap_outputs.outputs)
        skipped.update(shap_outputs.skipped)

    findings.extend(
        summarize_figure_findings(
            pred=pred,
            label_col=label,
            score_col=score_col,
            label_distribution=label_distribution.tables.get("label_distribution"),
            positive_rates=target_rates.tables.get("positive_rate"),
            topk=topk_outputs.tables.get("topk"),
            calibration=calibration_outputs.tables.get("calibration"),
            ad_summary=ad_outputs.tables.get("ad_summary"),
            failure_summary=failure_outputs.tables.get("failure_summary"),
            source_clustering=failure_outputs.tables.get("source_clustering"),
            scaffold_dominance=failure_outputs.tables.get("scaffold_dominance"),
            feature_importance=feature_outputs.tables.get("feature_importance"),
            shap_summary=outputs.get("shap_summary_csv"),
        )
    )

    manifest: dict[str, Any] = {
        "model_dir": str(model_path) if model_path else None,
        "predictions": str(pred_path),
        "dataset": str(Path(dataset).resolve()) if dataset else None,
        "out_dir": str(output_dir),
        "label_col": label,
        "score_col": score_col,
        "n_prediction_rows": int(len(pred)),
        "n_labeled_prediction_rows": int(_valid_binary_frame(pred, label, score_col).shape[0]),
        "group_cols": selected_groups,
        "target_cols": selected_targets,
        "top_ks": ks,
        "threshold": float(threshold),
        "top_n_errors": int(top_n_errors),
        "outputs": outputs,
        "skipped": skipped,
        "findings": findings,
    }
    manifest_path = output_dir / "figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    outputs["manifest"] = str(manifest_path)
    findings_path = write_findings_markdown(output_dir / "figure_findings.md", manifest)
    outputs["findings_markdown"] = str(findings_path)
    manifest["outputs"] = outputs
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return manifest


class FigureResult:
    def __init__(
        self,
        *,
        outputs: dict[str, str] | None = None,
        skipped: dict[str, str] | None = None,
        tables: dict[str, pd.DataFrame] | None = None,
    ) -> None:
        self.outputs = outputs or {}
        self.skipped = skipped or {}
        self.tables = tables or {}


def detect_binary_label(frame: pd.DataFrame) -> str:
    preferred = [
        "spd_binding_label",
        "spd_exposure_label",
        "spd_exposure_relevant",
        "combined_activity_ml_label",
        "external_four_state_ml_label",
        "mechanism_ml_label_clean",
        "mechanism_ml_label",
        "tissue_site_label",
    ]
    for col in preferred:
        if col in frame.columns and _is_binary_label_series(frame[col]):
            return col
    candidates = []
    for col in frame.columns:
        lowered = str(col).lower()
        if not lowered.endswith("label"):
            continue
        if "source" in lowered or "status" in lowered:
            continue
        if _is_binary_label_series(frame[col]):
            candidates.append(col)
    if candidates:
        return candidates[0]
    raise ValueError("could not auto-detect a binary label column; pass --label")


def write_label_distribution_outputs(
    frame: pd.DataFrame,
    *,
    label_col: str,
    group_cols: Sequence[str],
    out_dir: Path,
    max_groups: int,
    title_prefix: str,
) -> FigureResult:
    rows = []
    for group_col in group_cols:
        if group_col not in frame.columns:
            continue
        work = frame[[group_col, label_col]].copy()
        labels = pd.to_numeric(work[label_col], errors="coerce")
        valid = labels.isin([0, 1])
        work = work.loc[valid].copy()
        if work.empty:
            continue
        work["_label"] = labels.loc[valid].astype(int)
        work["_group"] = _string_group(work[group_col])
        grouped = work.groupby("_group", dropna=False)["_label"]
        for value, values in grouped:
            n = int(values.shape[0])
            n_positive = int(values.sum())
            rows.append(
                {
                    "group_col": group_col,
                    "group_value": value,
                    "n": n,
                    "n_positive": n_positive,
                    "n_negative": int(n - n_positive),
                    "positive_rate": float(n_positive / n) if n else np.nan,
                }
            )
    table = pd.DataFrame(rows)
    if table.empty:
        return FigureResult(skipped={"label_distribution": "no requested group columns with binary labels"}, tables={"label_distribution": table})
    csv_path = out_dir / "label_distribution_by_group.csv"
    table.sort_values(["group_col", "n"], ascending=[True, False]).to_csv(csv_path, index=False)
    fig_path = out_dir / "label_distribution_by_group.png"
    _plot_label_distribution(table, fig_path, max_groups=max_groups, title_prefix=title_prefix)
    return FigureResult(
        outputs={"label_distribution_csv": str(csv_path), "label_distribution_png": str(fig_path)},
        tables={"label_distribution": table},
    )


def write_positive_rate_outputs(
    frame: pd.DataFrame,
    *,
    label_col: str,
    group_cols: Sequence[str],
    out_dir: Path,
    max_groups: int,
    min_group_size: int,
    title_prefix: str,
) -> FigureResult:
    rows = []
    for group_col in group_cols:
        if group_col not in frame.columns:
            continue
        work = frame[[group_col, label_col]].copy()
        labels = pd.to_numeric(work[label_col], errors="coerce")
        valid = labels.isin([0, 1])
        work = work.loc[valid].copy()
        if work.empty:
            continue
        work["_label"] = labels.loc[valid].astype(int)
        work["_group"] = _string_group(work[group_col])
        for value, group in work.groupby("_group", dropna=False):
            n = int(group.shape[0])
            positives = int(group["_label"].sum())
            rows.append(
                {
                    "group_col": group_col,
                    "group_value": value,
                    "n": n,
                    "n_positive": positives,
                    "positive_rate": float(positives / n) if n else np.nan,
                }
            )
    table = pd.DataFrame(rows)
    if table.empty:
        return FigureResult(skipped={"positive_rate": "target/PDB columns absent or unlabeled"}, tables={"positive_rate": table})
    csv_path = out_dir / "positive_rate_by_target_pdb.csv"
    table.sort_values(["group_col", "n"], ascending=[True, False]).to_csv(csv_path, index=False)
    fig_path = out_dir / "positive_rate_by_target_pdb.png"
    _plot_positive_rates(table, fig_path, max_groups=max_groups, min_group_size=min_group_size, title_prefix=title_prefix)
    return FigureResult(
        outputs={"positive_rate_csv": str(csv_path), "positive_rate_png": str(fig_path)},
        tables={"positive_rate": table},
    )


def write_precision_recall_outputs(
    frame: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    out_dir: Path,
    title_prefix: str,
) -> FigureResult:
    work = _valid_binary_frame(frame, label_col, score_col)
    if work.empty or work[label_col].nunique() < 2:
        return FigureResult(skipped={"precision_recall": "requires at least one positive and one negative label"})
    from sklearn.metrics import average_precision_score, precision_recall_curve

    y = work[label_col].astype(int).to_numpy()
    scores = work[score_col].astype(float).to_numpy()
    precision, recall, thresholds = precision_recall_curve(y, scores)
    table = pd.DataFrame(
        {
            "precision": precision,
            "recall": recall,
            "threshold": list(thresholds) + [np.nan],
        }
    )
    ap = float(average_precision_score(y, scores))
    prevalence = float(y.mean())
    csv_path = out_dir / "precision_recall_curve.csv"
    table.to_csv(csv_path, index=False)
    fig_path = out_dir / "precision_recall_curve.png"
    _plot_precision_recall(table, fig_path, ap=ap, prevalence=prevalence, title_prefix=title_prefix)
    return FigureResult(
        outputs={"precision_recall_csv": str(csv_path), "precision_recall_png": str(fig_path)},
        tables={"precision_recall": table},
    )


def write_topk_outputs(
    frame: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    top_ks: Sequence[int],
    out_dir: Path,
    title_prefix: str,
) -> FigureResult:
    work = _valid_binary_frame(frame, label_col, score_col)
    if work.empty:
        return FigureResult(skipped={"topk": "no rows with binary labels and scores"})
    ranked = work.sort_values(score_col, ascending=False).reset_index(drop=True)
    effective_ks = _clip_and_dedupe_top_ks(top_ks, len(ranked))
    prevalence = float(ranked[label_col].mean()) if len(ranked) else np.nan
    rows = []
    for k in effective_ks:
        top = ranked.head(k)
        hits = int(top[label_col].sum())
        precision = float(hits / k)
        rows.append(
            {
                "top_k": k,
                "effective_k": k,
                "hit_count_at_k": hits,
                "precision_at_k": precision,
                "recall_at_k": float(hits / ranked[label_col].sum()) if ranked[label_col].sum() else np.nan,
                "prevalence": prevalence,
                "enrichment_at_k": float(precision / prevalence) if prevalence else np.nan,
            }
        )
    table = pd.DataFrame(rows)
    csv_path = out_dir / "topk_recovery.csv"
    table.to_csv(csv_path, index=False)
    fig_path = out_dir / "topk_recovery.png"
    _plot_topk(table, fig_path, title_prefix=title_prefix)
    return FigureResult(outputs={"topk_csv": str(csv_path), "topk_png": str(fig_path)}, tables={"topk": table})


def write_calibration_outputs(
    frame: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    out_dir: Path,
    title_prefix: str,
    bins: int = 10,
) -> FigureResult:
    work = _valid_binary_frame(frame, label_col, score_col)
    if work.empty:
        return FigureResult(skipped={"calibration": "no rows with binary labels and scores"})
    table = _calibration_table(work, label_col=label_col, score_col=score_col, bins=bins)
    csv_path = out_dir / "calibration_curve.csv"
    table.to_csv(csv_path, index=False)
    fig_path = out_dir / "calibration_curve.png"
    _plot_calibration(table, work[score_col].astype(float), fig_path, title_prefix=title_prefix)
    return FigureResult(outputs={"calibration_csv": str(csv_path), "calibration_png": str(fig_path)}, tables={"calibration": table})


def write_applicability_domain_outputs(
    frame: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    out_dir: Path,
    title_prefix: str,
) -> FigureResult:
    if "applicability_distance" not in frame.columns:
        return FigureResult(skipped={"applicability_domain": "applicability_distance column absent"})
    work = _valid_binary_frame(frame, label_col, score_col, extra_cols=["applicability_distance", "applicability_domain"])
    if work.empty:
        return FigureResult(skipped={"applicability_domain": "no labeled rows with applicability distances"})
    work["applicability_distance"] = pd.to_numeric(work["applicability_distance"], errors="coerce")
    work = work.dropna(subset=["applicability_distance"])
    if work.empty:
        return FigureResult(skipped={"applicability_domain": "all applicability distances are missing"})
    summary = _applicability_domain_summary(work, label_col=label_col, score_col=score_col)
    csv_path = out_dir / "applicability_domain_summary.csv"
    summary.to_csv(csv_path, index=False)
    fig_path = out_dir / "applicability_domain_plot.png"
    _plot_applicability_domain(work, label_col=label_col, score_col=score_col, fig_path=fig_path, title_prefix=title_prefix)
    return FigureResult(
        outputs={"applicability_domain_csv": str(csv_path), "applicability_domain_png": str(fig_path)},
        tables={"ad_summary": summary},
    )


def write_failure_analysis_outputs(
    frame: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    out_dir: Path,
    threshold: float,
    top_n: int,
    title_prefix: str,
) -> FigureResult:
    work = _failure_frame(frame, label_col=label_col, score_col=score_col, threshold=threshold)
    if work.empty:
        return FigureResult(skipped={"failure_analysis": "no rows with binary labels and scores"})

    outputs: dict[str, str] = {}
    tables: dict[str, pd.DataFrame] = {}
    context_cols = _context_columns(work, label_col, score_col)

    false_positives = _rank_threshold_errors(
        work,
        error_type="false_positive",
        score_col=score_col,
        ascending=False,
        top_n=top_n,
        context_cols=context_cols,
    )
    false_negatives = _rank_threshold_errors(
        work,
        error_type="false_negative",
        score_col=score_col,
        ascending=True,
        top_n=top_n,
        context_cols=context_cols,
    )
    fp_path = out_dir / "top_false_positive_candidates.csv"
    fn_path = out_dir / "top_false_negative_candidates.csv"
    false_positives.to_csv(fp_path, index=False)
    false_negatives.to_csv(fn_path, index=False)
    outputs["top_false_positives_csv"] = str(fp_path)
    outputs["top_false_negatives_csv"] = str(fn_path)
    tables["top_false_positives"] = false_positives
    tables["top_false_negatives"] = false_negatives

    failure_summary = _group_failure_summary(
        work,
        label_col=label_col,
        score_col=score_col,
        group_cols=[col for col in DEFAULT_FAILURE_GROUP_COLUMNS if col in work.columns],
    )
    if not failure_summary.empty:
        failure_path = out_dir / "failure_summary_by_group.csv"
        failure_summary.to_csv(failure_path, index=False)
        failure_png = out_dir / "failure_summary_by_group.png"
        _plot_failure_summary(failure_summary, failure_png, title_prefix=title_prefix)
        outputs["failure_summary_csv"] = str(failure_path)
        outputs["failure_summary_png"] = str(failure_png)
        tables["failure_summary"] = failure_summary

    scaffold_dominance = _dominance_summary(
        work,
        label_col=label_col,
        score_col=score_col,
        group_cols=[col for col in DEFAULT_SCAFFOLD_COLUMNS if col in work.columns],
    )
    if not scaffold_dominance.empty:
        scaffold_path = out_dir / "scaffold_dominance.csv"
        scaffold_dominance.to_csv(scaffold_path, index=False)
        scaffold_png = out_dir / "scaffold_dominance.png"
        _plot_dominance(scaffold_dominance, scaffold_png, title_prefix=title_prefix, title="scaffold dominance")
        outputs["scaffold_dominance_csv"] = str(scaffold_path)
        outputs["scaffold_dominance_png"] = str(scaffold_png)
        tables["scaffold_dominance"] = scaffold_dominance

    source_clustering = _dominance_summary(
        work,
        label_col=label_col,
        score_col=score_col,
        group_cols=[col for col in DEFAULT_SOURCE_COLUMNS if col in work.columns],
    )
    if not source_clustering.empty:
        source_path = out_dir / "positive_source_clustering.csv"
        source_clustering.to_csv(source_path, index=False)
        source_png = out_dir / "positive_source_clustering.png"
        _plot_dominance(source_clustering, source_png, title_prefix=title_prefix, title="positive clustering by source")
        outputs["positive_source_clustering_csv"] = str(source_path)
        outputs["positive_source_clustering_png"] = str(source_png)
        tables["source_clustering"] = source_clustering

    return FigureResult(outputs=outputs, tables=tables)


def _rank_threshold_errors(
    frame: pd.DataFrame,
    *,
    error_type: str,
    score_col: str,
    ascending: bool,
    top_n: int,
    context_cols: Sequence[str],
) -> pd.DataFrame:
    """Filter to one threshold error class before ranking its candidates."""

    return (
        frame.loc[frame["error_type"].eq(error_type)]
        .sort_values(score_col, ascending=ascending, kind="stable")
        .head(max(0, int(top_n)))[list(context_cols)]
        .copy()
    )


def write_feature_importance_outputs(
    *,
    model_path: Path | None,
    out_dir: Path,
    title_prefix: str,
    max_features: int,
) -> FigureResult:
    if model_path is None:
        return FigureResult(skipped={"feature_importance": "requires --model-dir"})
    path = model_path / "feature_importance.csv"
    if not path.exists():
        return FigureResult(skipped={"feature_importance": f"not found: {path}"})
    table = pd.read_csv(path)
    if table.empty or "feature" not in table.columns or "importance" not in table.columns:
        return FigureResult(skipped={"feature_importance": f"invalid feature importance table: {path}"})
    copied_csv = out_dir / "feature_importance_top.csv"
    out = table.copy()
    out["abs_importance"] = pd.to_numeric(out["importance"], errors="coerce").abs()
    out.sort_values("abs_importance", ascending=False).head(max_features).to_csv(copied_csv, index=False)
    fig_path = out_dir / "feature_importance_top.png"
    _plot_feature_importance(out, fig_path, max_features=max_features, title_prefix=title_prefix)
    return FigureResult(
        outputs={"feature_importance_csv": str(copied_csv), "feature_importance_png": str(fig_path)},
        tables={"feature_importance": out},
    )


def write_shap_outputs(
    *,
    model_path: Path | None,
    out_dir: Path,
    sample_rows: int,
    background_rows: int,
    random_state: int,
    force_generic: bool,
    max_generic_features: int,
    max_features: int,
    title_prefix: str,
) -> FigureResult:
    if model_path is None:
        return FigureResult(skipped={"shap": "requires --model-dir"})
    model_file = model_path / "trained_model.pkl"
    design_file = model_path / "training_design_matrix.csv"
    if not model_file.exists() or not design_file.exists():
        return FigureResult(skipped={"shap": "requires trained_model.pkl and training_design_matrix.csv"})
    try:
        import shap  # type: ignore[import-not-found]
    except ImportError:
        return FigureResult(skipped={"shap": "optional package shap is not installed"})

    try:
        with model_file.open("rb") as handle:
            model = pickle.load(handle)
    except Exception as exc:
        return FigureResult(skipped={"shap": f"could not load trained model: {exc}"})

    try:
        x_all = pd.read_csv(design_file, low_memory=False)
    except Exception as exc:
        return FigureResult(skipped={"shap": f"could not load training design matrix: {exc}"})
    if x_all.empty:
        return FigureResult(skipped={"shap": "training design matrix is empty"})
    x_all = x_all.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    x_sample = _sample_rows(x_all, sample_rows, random_state=random_state)
    background = _sample_rows(x_all, background_rows, random_state=random_state + 1)

    try:
        values = _compute_shap_values(
            shap_module=shap,
            model=model,
            x_sample=x_sample,
            background=background,
            force_generic=force_generic,
            max_generic_features=max_generic_features,
        )
    except Exception as exc:
        return FigureResult(skipped={"shap": f"SHAP computation failed: {exc}"})

    if values.shape[0] != x_sample.shape[0] or values.shape[1] != x_sample.shape[1]:
        return FigureResult(skipped={"shap": f"unexpected SHAP shape {values.shape} for design matrix {x_sample.shape}"})
    summary = pd.DataFrame(
        {
            "feature": list(x_sample.columns),
            "mean_abs_shap": np.abs(values).mean(axis=0),
            "mean_shap": values.mean(axis=0),
        }
    ).sort_values("mean_abs_shap", ascending=False)
    summary_csv = out_dir / "shap_summary.csv"
    summary.to_csv(summary_csv, index=False)
    bar_png = out_dir / "shap_summary.png"
    _plot_shap_summary(summary, bar_png, max_features=max_features, title_prefix=title_prefix)
    outputs = {"shap_summary_csv": str(summary_csv), "shap_summary_png": str(bar_png)}

    beeswarm_png = out_dir / "shap_beeswarm.png"
    if _plot_shap_beeswarm(values, x_sample, beeswarm_png, max_features=max_features):
        outputs["shap_beeswarm_png"] = str(beeswarm_png)
    return FigureResult(outputs=outputs, tables={"shap_summary": summary})


def summarize_figure_findings(
    *,
    pred: pd.DataFrame,
    label_col: str,
    score_col: str,
    label_distribution: pd.DataFrame | None,
    positive_rates: pd.DataFrame | None,
    topk: pd.DataFrame | None,
    calibration: pd.DataFrame | None,
    ad_summary: pd.DataFrame | None,
    failure_summary: pd.DataFrame | None,
    source_clustering: pd.DataFrame | None,
    scaffold_dominance: pd.DataFrame | None,
    feature_importance: pd.DataFrame | None,
    shap_summary: str | None,
) -> list[str]:
    findings: list[str] = []
    work = _valid_binary_frame(pred, label_col, score_col)
    if work.empty:
        return ["No rows had binary labels plus prediction scores; figure interpretation is limited."]
    n = int(len(work))
    positives = int(work[label_col].sum())
    prevalence = positives / n if n else math.nan
    findings.append(f"Prediction set has {n} labeled rows, {positives} positives, and prevalence {prevalence:.3f}.")

    try:
        from sklearn.metrics import average_precision_score, roc_auc_score

        if work[label_col].nunique() >= 2:
            ap = float(average_precision_score(work[label_col].astype(int), work[score_col].astype(float)))
            auroc = float(roc_auc_score(work[label_col].astype(int), work[score_col].astype(float)))
            findings.append(f"PR curve AP/AUPRC is {ap:.3f} versus baseline prevalence {prevalence:.3f}; AUROC is {auroc:.3f}.")
            if prevalence and ap < prevalence * 2:
                findings.append("AUPRC is less than 2x prevalence, so rank ordering may be weak for short-list triage.")
    except Exception:
        pass

    predicted_at_half = int((work[score_col].astype(float) >= 0.5).sum())
    if predicted_at_half == 0:
        findings.append("No rows score above 0.5; use rank-based thresholds or calibrate a decision threshold instead of using 0.5.")

    if topk is not None and not topk.empty:
        first = topk.sort_values("effective_k").iloc[0]
        findings.append(
            "Top-K recovery starts at "
            f"K={int(first['effective_k'])}: {int(first['hit_count_at_k'])} hits, "
            f"precision {float(first['precision_at_k']):.3f}, "
            f"enrichment {float(first['enrichment_at_k']):.2f}x."
        )

    if calibration is not None and not calibration.empty:
        weighted_abs = (
            (calibration["n"].astype(float) * (calibration["observed_positive_rate"] - calibration["mean_predicted_probability"]).abs()).sum()
            / max(1.0, calibration["n"].astype(float).sum())
        )
        findings.append(f"Mean absolute calibration error from plotted bins is {weighted_abs:.3f}.")
        if weighted_abs > 0.05:
            findings.append("Calibration error is above 0.05; probability values should be treated cautiously even if ranking is useful.")

    if ad_summary is not None and not ad_summary.empty and "applicability_domain" in ad_summary.columns:
        out_rows = ad_summary[ad_summary["applicability_domain"].astype(str).str.contains("out", case=False, na=False)]
        if not out_rows.empty:
            frac = float(out_rows["n"].sum() / ad_summary["n"].sum()) if ad_summary["n"].sum() else 0.0
            findings.append(f"Applicability-domain summary marks {frac:.1%} of labeled predictions as out-of-domain.")
            if frac > 0.10:
                findings.append("Out-of-domain fraction exceeds 10%; review AD-stratified performance before presenting global metrics.")

    if label_distribution is not None and not label_distribution.empty:
        for group_col in ["source_family", "label_source", "target_family", "scaffold_key", "chemical_cluster"]:
            subset = label_distribution[label_distribution["group_col"].eq(group_col)]
            if subset.empty:
                continue
            dominant = subset.sort_values("n", ascending=False).iloc[0]
            share = float(dominant["n"] / subset["n"].sum()) if subset["n"].sum() else 0.0
            if share > 0.50:
                findings.append(
                    f"{group_col} is concentrated: {dominant['group_value']} accounts for {share:.1%} of rows."
                )
            if group_col in {"scaffold_key", "chemical_cluster"}:
                small = int((subset["n"] < 3).sum())
                if small:
                    findings.append(f"{group_col} is sparse: {small} groups have fewer than 3 labeled rows.")

    if positive_rates is not None and not positive_rates.empty:
        reliable = positive_rates[positive_rates["n"] >= 5]
        if not reliable.empty:
            high = reliable.sort_values("positive_rate", ascending=False).head(1).iloc[0]
            findings.append(
                f"Highest positive-rate group with n>=5 is {high['group_col']}={high['group_value']} "
                f"at {float(high['positive_rate']):.3f} over n={int(high['n'])}."
            )

    if failure_summary is not None and not failure_summary.empty:
        for group_col in ["target_family", "protein_class", "scaffold_key", "chemical_cluster", "label_source", "source_family"]:
            subset = failure_summary[failure_summary["group_col"].eq(group_col)]
            subset = subset[subset["n"] >= 5]
            if subset.empty:
                continue
            worst = subset.sort_values(["error_pressure", "n"], ascending=[False, False]).iloc[0]
            findings.append(
                f"Worst {group_col} failure-pressure group is {worst['group_value']} "
                f"(n={int(worst['n'])}, error_pressure={float(worst['error_pressure']):.3f})."
            )
            break

    if source_clustering is not None and not source_clustering.empty:
        source = source_clustering[source_clustering["group_col"].isin(DEFAULT_SOURCE_COLUMNS)]
        if not source.empty:
            dominant = source.sort_values("positive_share", ascending=False).iloc[0]
            findings.append(
                f"Positive labels are most concentrated in {dominant['group_col']}={dominant['group_value']} "
                f"with {float(dominant['positive_share']):.1%} of positives."
            )
            if float(dominant["positive_share"]) >= 0.50:
                findings.append("At least half of positives come from one source group; source-transfer claims need extra care.")

    if scaffold_dominance is not None and not scaffold_dominance.empty:
        scaffold = scaffold_dominance.sort_values("positive_share", ascending=False).iloc[0]
        if float(scaffold["positive_share"]) >= 0.25:
            findings.append(
                f"Scaffold/chemotype positives are concentrated: {scaffold['group_col']}={scaffold['group_value']} "
                f"contains {float(scaffold['positive_share']):.1%} of positives."
            )

    importance = _top_importance_frame(feature_importance, shap_summary)
    if importance is not None and not importance.empty:
        top_features = importance["feature"].astype(str).head(10).tolist()
        findings.append("Top explanatory features include: " + ", ".join(top_features[:5]) + ".")
        risky = [feature for feature in top_features if any(token in feature.lower() for token in ID_OR_PROVENANCE_TOKENS)]
        if risky:
            findings.append(
                "Feature explanation is dominated by provenance/group features: "
                + ", ".join(risky[:5])
                + ". Treat this as a possible source, scaffold, or family shortcut to audit."
            )

    return findings


def write_findings_markdown(path: Path, manifest: dict[str, Any]) -> Path:
    lines = [
        "# ML Figure Findings",
        "",
        f"Predictions: `{manifest['predictions']}`",
        f"Label: `{manifest['label_col']}`",
        f"Score: `{manifest['score_col']}`",
        f"Rows: `{manifest['n_prediction_rows']}`",
        "",
        "## Findings",
        "",
    ]
    for finding in manifest.get("findings", []):
        lines.append(f"- {finding}")
    lines.extend(["", "## Outputs", ""])
    for key, value in sorted(dict(manifest.get("outputs", {})).items()):
        lines.append(f"- `{key}`: `{value}`")
    if manifest.get("skipped"):
        lines.extend(["", "## Skipped", ""])
        for key, value in sorted(dict(manifest.get("skipped", {})).items()):
            lines.append(f"- `{key}`: {value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _resolve_predictions(model_path: Path | None, predictions: str | Path | None) -> Path:
    if predictions:
        path = Path(predictions).resolve()
    elif model_path is not None:
        path = model_path / "model_predictions.csv"
    else:
        raise ValueError("provide --model-dir or --predictions")
    if not path.exists():
        raise FileNotFoundError(f"model predictions not found: {path}")
    return path


def _resolve_out_dir(model_path: Path | None, out_dir: str | Path | None) -> Path:
    if out_dir:
        return Path(out_dir).resolve()
    if model_path is not None:
        return model_path / "figures"
    return Path("figures").resolve()


def _default_title_prefix(model_path: Path | None, pred_path: Path) -> str:
    if model_path is not None:
        parts = model_path.parts
        return "/".join(parts[-4:-1]) if len(parts) >= 4 else model_path.name
    return pred_path.parent.name


def _load_label_frame(dataset: str | Path | None, pred: pd.DataFrame, label_col: str) -> pd.DataFrame:
    if not dataset:
        return pred
    path = Path(dataset)
    if not path.exists():
        raise FileNotFoundError(f"dataset not found: {path}")
    if path.suffix == ".parquet":
        frame = pd.read_parquet(path)
    else:
        frame = pd.read_csv(path, low_memory=False)
    if label_col not in frame.columns:
        return pred
    return frame


def _is_binary_label_series(series: pd.Series) -> bool:
    values = pd.to_numeric(series, errors="coerce").dropna().unique()
    if len(values) == 0:
        return False
    return set(values).issubset({0, 1, 0.0, 1.0})


def _string_group(series: pd.Series) -> pd.Series:
    return series.astype("object").where(pd.notna(series), "missing").astype(str).replace({"": "missing"})


def _valid_binary_frame(
    frame: pd.DataFrame,
    label_col: str,
    score_col: str,
    extra_cols: Sequence[str] | None = None,
) -> pd.DataFrame:
    cols = [label_col, score_col, *(extra_cols or [])]
    cols = [col for col in cols if col in frame.columns]
    work = frame[cols].copy()
    labels = pd.to_numeric(work[label_col], errors="coerce")
    scores = pd.to_numeric(work[score_col], errors="coerce")
    valid = labels.isin([0, 1]) & scores.notna()
    work = work.loc[valid].copy()
    work[label_col] = labels.loc[valid].astype(int)
    work[score_col] = scores.loc[valid].astype(float)
    return work


def _normalize_top_ks(values: Sequence[int], n_rows: int) -> list[int]:
    dynamic = [max(1, int(math.ceil(n_rows * frac))) for frac in (0.01, 0.05, 0.10)]
    return _clip_and_dedupe_top_ks([*values, *dynamic], n_rows)


def _clip_and_dedupe_top_ks(values: Sequence[int], n_rows: int) -> list[int]:
    if n_rows <= 0:
        return []
    return sorted({min(int(value), n_rows) for value in values if int(value) > 0})


def _configure_matplotlib() -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)


def _plot_label_distribution(table: pd.DataFrame, fig_path: Path, *, max_groups: int, title_prefix: str) -> None:
    import matplotlib.pyplot as plt

    groups = [col for col in ["target_family", "source_family", "label_source", "scaffold_key", "chemical_cluster"] if col in set(table["group_col"])]
    if not groups:
        groups = sorted(table["group_col"].unique().tolist())[:3]
    fig, axes = plt.subplots(len(groups), 1, figsize=(10, max(3.2, 3.0 * len(groups))), squeeze=False)
    for ax, group_col in zip(axes.ravel(), groups):
        subset = table[table["group_col"].eq(group_col)].sort_values("n", ascending=False).head(max_groups).iloc[::-1]
        y = np.arange(len(subset))
        ax.barh(y, subset["n_negative"], color="#d8dee9", label="negative")
        ax.barh(y, subset["n_positive"], left=subset["n_negative"], color="#bf616a", label="positive")
        ax.set_yticks(y, [_short_label(v) for v in subset["group_value"]])
        ax.set_xlabel("labeled rows")
        ax.set_title(group_col)
        ax.legend(loc="lower right", frameon=False)
    fig.suptitle(f"{title_prefix}: label distribution by group")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _plot_positive_rates(
    table: pd.DataFrame,
    fig_path: Path,
    *,
    max_groups: int,
    min_group_size: int,
    title_prefix: str,
) -> None:
    import matplotlib.pyplot as plt

    groups = sorted(table["group_col"].unique().tolist())
    fig, axes = plt.subplots(len(groups), 1, figsize=(10, max(3.2, 3.0 * len(groups))), squeeze=False)
    for ax, group_col in zip(axes.ravel(), groups):
        subset = table[(table["group_col"].eq(group_col)) & (table["n"] >= min_group_size)].copy()
        if subset.empty:
            subset = table[table["group_col"].eq(group_col)].copy()
        subset = subset.sort_values(["positive_rate", "n"], ascending=[False, False]).head(max_groups).iloc[::-1]
        y = np.arange(len(subset))
        sizes = np.maximum(20, np.sqrt(subset["n"].astype(float)) * 12)
        ax.scatter(subset["positive_rate"], y, s=sizes, color="#5e81ac", alpha=0.85)
        for idx, row in enumerate(subset.itertuples(index=False)):
            ax.text(float(row.positive_rate) + 0.01, idx, f"n={int(row.n)}", va="center", fontsize=8)
        ax.set_xlim(0, min(1.05, max(0.15, float(subset["positive_rate"].max()) + 0.1 if not subset.empty else 1.0)))
        ax.set_yticks(y, [_short_label(v) for v in subset["group_value"]])
        ax.set_xlabel("positive rate")
        ax.set_title(f"{group_col} positive rate")
    fig.suptitle(f"{title_prefix}: positive rate by target/PDB")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _plot_precision_recall(table: pd.DataFrame, fig_path: Path, *, ap: float, prevalence: float, title_prefix: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    ax.plot(table["recall"], table["precision"], color="#5e81ac", linewidth=2)
    ax.axhline(prevalence, color="#a3a3a3", linestyle="--", linewidth=1, label=f"baseline={prevalence:.3f}")
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_title(f"{title_prefix}: precision-recall (AP={ap:.3f})")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _plot_topk(table: pd.DataFrame, fig_path: Path, *, title_prefix: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(table["effective_k"], table["precision_at_k"], marker="o", label="precision@K", color="#5e81ac")
    ax.plot(table["effective_k"], table["recall_at_k"], marker="o", label="recall@K", color="#a3be8c")
    ax.axhline(float(table["prevalence"].iloc[0]), color="#a3a3a3", linestyle="--", linewidth=1, label="prevalence")
    ax.set_xscale("log")
    ax.set_xlabel("K")
    ax.set_ylabel("rate")
    ax.set_ylim(0, min(1.05, max(0.2, float(table[["precision_at_k", "recall_at_k"]].max().max()) + 0.1)))
    ax.set_title(f"{title_prefix}: top-K recovery")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _calibration_table(work: pd.DataFrame, *, label_col: str, score_col: str, bins: int) -> pd.DataFrame:
    tmp = work[[label_col, score_col]].copy()
    try:
        tmp["_bin"] = pd.qcut(tmp[score_col].rank(method="first"), q=min(bins, len(tmp)), duplicates="drop")
    except ValueError:
        tmp["_bin"] = pd.cut(tmp[score_col], bins=bins, include_lowest=True)
    rows = []
    for idx, (_bin, group) in enumerate(tmp.groupby("_bin", observed=False)):
        rows.append(
            {
                "bin_index": idx,
                "n": int(len(group)),
                "mean_predicted_probability": float(group[score_col].mean()),
                "observed_positive_rate": float(group[label_col].mean()),
                "score_min": float(group[score_col].min()),
                "score_max": float(group[score_col].max()),
            }
        )
    return pd.DataFrame(rows)


def _plot_calibration(table: pd.DataFrame, scores: pd.Series, fig_path: Path, *, title_prefix: str) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), gridspec_kw={"width_ratios": [2, 1]})
    axes[0].plot([0, 1], [0, 1], color="#a3a3a3", linestyle="--", linewidth=1)
    axes[0].plot(table["mean_predicted_probability"], table["observed_positive_rate"], marker="o", color="#5e81ac")
    axes[0].set_xlabel("mean predicted probability")
    axes[0].set_ylabel("observed positive rate")
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].set_title("calibration")
    axes[1].hist(scores, bins=20, color="#d8dee9", edgecolor="#4c566a")
    axes[1].set_xlabel("score")
    axes[1].set_ylabel("rows")
    axes[1].set_title("score histogram")
    fig.suptitle(f"{title_prefix}: calibration")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _applicability_domain_summary(work: pd.DataFrame, *, label_col: str, score_col: str) -> pd.DataFrame:
    if "applicability_domain" not in work.columns:
        work = work.copy()
        threshold = float(work["applicability_distance"].quantile(0.95))
        work["applicability_domain"] = np.where(work["applicability_distance"] <= threshold, "in_domain", "out_of_domain")
    rows = []
    for domain, group in work.groupby("applicability_domain", dropna=False):
        rows.append(
            {
                "applicability_domain": domain,
                "n": int(len(group)),
                "positive_rate": float(group[label_col].mean()),
                "mean_prediction": float(group[score_col].mean()),
                "median_distance": float(group["applicability_distance"].median()),
                "distance_q95": float(group["applicability_distance"].quantile(0.95)),
            }
        )
    return pd.DataFrame(rows)


def _plot_applicability_domain(
    work: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    fig_path: Path,
    title_prefix: str,
) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    colors = np.where(work[label_col].astype(int).eq(1), "#bf616a", "#5e81ac")
    axes[0].scatter(work["applicability_distance"], work[score_col], c=colors, alpha=0.65, s=18, linewidth=0)
    axes[0].set_xlabel("applicability distance")
    axes[0].set_ylabel("prediction score")
    axes[0].set_title("score vs distance")
    if "applicability_domain" in work.columns:
        summary = work.groupby("applicability_domain", dropna=False)[label_col].agg(["size", "mean"]).reset_index()
        axes[1].bar(summary["applicability_domain"].astype(str), summary["mean"], color="#a3be8c")
        for idx, row in enumerate(summary.itertuples(index=False)):
            axes[1].text(idx, float(row.mean) + 0.01, f"n={int(row.size)}", ha="center", fontsize=8)
        axes[1].set_ylim(0, min(1.05, max(0.15, float(summary["mean"].max()) + 0.1)))
        axes[1].set_xlabel("domain")
        axes[1].set_ylabel("positive rate")
        axes[1].set_title("label rate by domain")
    else:
        axes[1].hist(work["applicability_distance"], bins=20, color="#d8dee9", edgecolor="#4c566a")
        axes[1].set_title("distance histogram")
    fig.suptitle(f"{title_prefix}: applicability domain")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _plot_feature_importance(table: pd.DataFrame, fig_path: Path, *, max_features: int, title_prefix: str) -> None:
    import matplotlib.pyplot as plt

    subset = table.sort_values("abs_importance", ascending=False).head(max_features).iloc[::-1]
    colors = np.where(pd.to_numeric(subset["importance"], errors="coerce") >= 0, "#5e81ac", "#bf616a")
    fig, ax = plt.subplots(figsize=(9, max(4, 0.28 * len(subset) + 1.5)))
    ax.barh(np.arange(len(subset)), subset["importance"].astype(float), color=colors)
    ax.set_yticks(np.arange(len(subset)), [_short_label(v, length=52) for v in subset["feature"]])
    ax.axvline(0, color="#4c566a", linewidth=0.8)
    ax.set_xlabel("importance")
    ax.set_title(f"{title_prefix}: feature importance")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _compute_shap_values(
    *,
    shap_module: Any,
    model: Any,
    x_sample: pd.DataFrame,
    background: pd.DataFrame,
    force_generic: bool,
    max_generic_features: int,
) -> np.ndarray:
    model_name = type(model).__name__.lower()
    if not force_generic and any(token in model_name for token in ["xgb", "forest", "lgbm", "gradientboost", "catboost"]):
        explainer = shap_module.TreeExplainer(model)
        raw = explainer.shap_values(x_sample)
    elif not force_generic and hasattr(model, "coef_"):
        explainer = shap_module.LinearExplainer(model, background)
        raw = explainer.shap_values(x_sample)
    else:
        if x_sample.shape[1] > max_generic_features and not force_generic:
            raise RuntimeError(
                f"generic SHAP skipped for {x_sample.shape[1]} features; use --force-generic-shap to run it"
            )

        columns = list(x_sample.columns)

        def predict_fn(values: np.ndarray) -> np.ndarray:
            frame = pd.DataFrame(values, columns=columns)
            if hasattr(model, "predict_proba"):
                return np.asarray(model.predict_proba(frame))[:, 1]
            raw_scores = np.asarray(model.decision_function(frame), dtype=float)
            return 1.0 / (1.0 + np.exp(-raw_scores))

        explainer = shap_module.Explainer(predict_fn, background, feature_names=columns)
        raw = explainer(x_sample)
    return _coerce_shap_array(raw)


def _coerce_shap_array(raw: Any) -> np.ndarray:
    if hasattr(raw, "values"):
        values = np.asarray(raw.values)
    elif isinstance(raw, list):
        values = np.asarray(raw[1] if len(raw) > 1 else raw[0])
    else:
        values = np.asarray(raw)
    if values.ndim == 3:
        values = values[:, :, 1] if values.shape[2] > 1 else values[:, :, 0]
    if values.ndim == 1:
        values = values.reshape(1, -1)
    return values.astype(float)


def _plot_shap_summary(summary: pd.DataFrame, fig_path: Path, *, max_features: int, title_prefix: str) -> None:
    import matplotlib.pyplot as plt

    subset = summary.sort_values("mean_abs_shap", ascending=False).head(max_features).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, max(4, 0.28 * len(subset) + 1.5)))
    ax.barh(np.arange(len(subset)), subset["mean_abs_shap"], color="#5e81ac")
    ax.set_yticks(np.arange(len(subset)), [_short_label(v, length=52) for v in subset["feature"]])
    ax.set_xlabel("mean absolute SHAP value")
    ax.set_title(f"{title_prefix}: SHAP summary")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _plot_shap_beeswarm(values: np.ndarray, x_sample: pd.DataFrame, fig_path: Path, *, max_features: int) -> bool:
    try:
        import matplotlib.pyplot as plt
        import shap  # type: ignore[import-not-found]

        shap.summary_plot(values, x_sample, show=False, max_display=max_features)
        plt.tight_layout()
        plt.savefig(fig_path, dpi=180, bbox_inches="tight")
        plt.close()
        return True
    except Exception:
        return False


def _sample_rows(frame: pd.DataFrame, n_rows: int, *, random_state: int) -> pd.DataFrame:
    if n_rows <= 0 or len(frame) <= n_rows:
        return frame.copy()
    return frame.sample(n=n_rows, random_state=random_state).copy()


def _failure_frame(frame: pd.DataFrame, *, label_col: str, score_col: str, threshold: float) -> pd.DataFrame:
    context = [
        col
        for col in [
            "drug_id",
            "target_id",
            "pdb_id",
            "target_family",
            "protein_class",
            "scaffold_key",
            "chemical_cluster",
            "ligand_chemotype",
            "label_source",
            "source_family",
            "upstream_source",
            "external_source_family",
            "applicability_domain",
            "applicability_distance",
        ]
        if col in frame.columns
    ]
    work = _valid_binary_frame(frame, label_col, score_col, extra_cols=context)
    if work.empty:
        return work
    work["predicted_label"] = (work[score_col].astype(float) >= float(threshold)).astype(int)
    work["error_type"] = np.select(
        [work[label_col].eq(0) & work["predicted_label"].eq(1), work[label_col].eq(1) & work["predicted_label"].eq(0)],
        ["false_positive", "false_negative"],
        default="correct",
    )
    work["is_threshold_error"] = work["error_type"].ne("correct")
    work["false_positive_pressure"] = np.where(work[label_col].eq(0), work[score_col].astype(float), 0.0)
    work["false_negative_pressure"] = np.where(work[label_col].eq(1), 1.0 - work[score_col].astype(float), 0.0)
    work["error_pressure"] = np.where(
        work[label_col].eq(1),
        work["false_negative_pressure"],
        work["false_positive_pressure"],
    )
    work["absolute_probability_error"] = (work[score_col].astype(float) - work[label_col].astype(float)).abs()
    return work


def _context_columns(frame: pd.DataFrame, label_col: str, score_col: str) -> list[str]:
    preferred = [
        "drug_id",
        "target_id",
        "pdb_id",
        "target_family",
        "protein_class",
        "scaffold_key",
        "chemical_cluster",
        "ligand_chemotype",
        "label_source",
        "source_family",
        "upstream_source",
        "external_source_family",
        "applicability_domain",
        "applicability_distance",
        label_col,
        score_col,
        "predicted_label",
        "error_type",
        "false_positive_pressure",
        "false_negative_pressure",
        "error_pressure",
        "absolute_probability_error",
    ]
    return [col for col in preferred if col in frame.columns]


def _group_failure_summary(
    frame: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    group_cols: Sequence[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group_col in group_cols:
        if group_col not in frame.columns:
            continue
        work = frame.copy()
        work["_group"] = _string_group(work[group_col])
        for value, group in work.groupby("_group", dropna=False):
            labels = group[label_col].astype(int)
            positives = int(labels.sum())
            negatives = int(len(group) - positives)
            fp_count = int(group["error_type"].eq("false_positive").sum())
            fn_count = int(group["error_type"].eq("false_negative").sum())
            rows.append(
                {
                    "group_col": group_col,
                    "group_value": value,
                    "n": int(len(group)),
                    "n_positive": positives,
                    "n_negative": negatives,
                    "positive_rate": float(labels.mean()) if len(group) else np.nan,
                    "mean_prediction": float(group[score_col].mean()),
                    "false_positive_count": fp_count,
                    "false_negative_count": fn_count,
                    "threshold_error_rate": float(group["is_threshold_error"].mean()),
                    "false_positive_rate_within_negatives": float(fp_count / negatives) if negatives else np.nan,
                    "false_negative_rate_within_positives": float(fn_count / positives) if positives else np.nan,
                    "mean_false_positive_pressure": float(group.loc[labels.eq(0), "false_positive_pressure"].mean()) if negatives else np.nan,
                    "mean_false_negative_pressure": float(group.loc[labels.eq(1), "false_negative_pressure"].mean()) if positives else np.nan,
                    "mean_abs_probability_error": float(group["absolute_probability_error"].mean()),
                    "error_pressure": float(group["error_pressure"].mean()),
                }
            )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["error_pressure", "n"], ascending=[False, False], kind="stable")


def _dominance_summary(
    frame: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    group_cols: Sequence[str],
) -> pd.DataFrame:
    total = len(frame)
    total_pos = int(frame[label_col].sum())
    total_neg = int(total - total_pos)
    score_cutoff = float(frame[score_col].quantile(0.90)) if total else np.nan
    high_score = frame[score_col].astype(float) >= score_cutoff
    high_score_neg_total = int((high_score & frame[label_col].eq(0)).sum())
    rows: list[dict[str, object]] = []
    for group_col in group_cols:
        if group_col not in frame.columns:
            continue
        work = frame.copy()
        work["_group"] = _string_group(work[group_col])
        for value, group in work.groupby("_group", dropna=False):
            labels = group[label_col].astype(int)
            n = int(len(group))
            positives = int(labels.sum())
            negatives = int(n - positives)
            group_high_neg = int(((group[score_col].astype(float) >= score_cutoff) & labels.eq(0)).sum())
            rows.append(
                {
                    "group_col": group_col,
                    "group_value": value,
                    "n": n,
                    "row_share": float(n / total) if total else np.nan,
                    "n_positive": positives,
                    "n_negative": negatives,
                    "positive_rate": float(positives / n) if n else np.nan,
                    "positive_share": float(positives / total_pos) if total_pos else 0.0,
                    "negative_share": float(negatives / total_neg) if total_neg else 0.0,
                    "mean_prediction": float(group[score_col].mean()),
                    "top_decile_negative_share": float(group_high_neg / high_score_neg_total) if high_score_neg_total else 0.0,
                    "mean_error_pressure": float(group["error_pressure"].mean()),
                }
            )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    return out.sort_values(["positive_share", "n"], ascending=[False, False], kind="stable")


def _plot_failure_summary(table: pd.DataFrame, fig_path: Path, *, title_prefix: str, max_groups: int = 12) -> None:
    import matplotlib.pyplot as plt

    group_cols = [col for col in ["target_family", "scaffold_key", "chemical_cluster", "label_source", "source_family"] if col in set(table["group_col"])]
    if not group_cols:
        group_cols = sorted(table["group_col"].unique().tolist())[:3]
    fig, axes = plt.subplots(len(group_cols), 1, figsize=(10, max(3.2, 2.8 * len(group_cols))), squeeze=False)
    for ax, group_col in zip(axes.ravel(), group_cols):
        subset = table[table["group_col"].eq(group_col)].sort_values(["error_pressure", "n"], ascending=[False, False]).head(max_groups).iloc[::-1]
        y = np.arange(len(subset))
        ax.barh(y, subset["error_pressure"].astype(float), color="#bf616a")
        ax.set_yticks(y, [_short_label(v) for v in subset["group_value"]])
        ax.set_xlabel("mean failure pressure")
        ax.set_title(group_col)
    fig.suptitle(f"{title_prefix}: failure pressure by group")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _plot_dominance(table: pd.DataFrame, fig_path: Path, *, title_prefix: str, title: str, max_groups: int = 15) -> None:
    import matplotlib.pyplot as plt

    subset = table.sort_values(["positive_share", "n"], ascending=[False, False]).head(max_groups).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, max(4, 0.3 * len(subset) + 1.5)))
    y = np.arange(len(subset))
    ax.barh(y, subset["row_share"].astype(float), color="#d8dee9", label="row share")
    ax.barh(y, subset["positive_share"].astype(float), color="#bf616a", alpha=0.82, label="positive share")
    ax.set_yticks(y, [f"{row.group_col}={_short_label(row.group_value, length=30)}" for row in subset.itertuples(index=False)])
    ax.set_xlabel("share")
    ax.set_title(f"{title_prefix}: {title}")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)


def _top_importance_frame(feature_importance: pd.DataFrame | None, shap_summary: str | None) -> pd.DataFrame | None:
    if shap_summary and Path(shap_summary).exists():
        try:
            shap_frame = pd.read_csv(shap_summary)
            if "mean_abs_shap" in shap_frame.columns:
                return shap_frame.sort_values("mean_abs_shap", ascending=False)
        except Exception:
            pass
    if feature_importance is None or feature_importance.empty:
        return None
    if "abs_importance" not in feature_importance.columns and "importance" in feature_importance.columns:
        feature_importance = feature_importance.copy()
        feature_importance["abs_importance"] = pd.to_numeric(feature_importance["importance"], errors="coerce").abs()
    return feature_importance.sort_values("abs_importance", ascending=False)


def _short_label(value: object, *, length: int = 38) -> str:
    text = str(value)
    if len(text) <= length:
        return text
    return text[: max(1, length - 3)] + "..."
