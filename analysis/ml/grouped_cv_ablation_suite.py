from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import pandas as pd

from analysis.ml.dataset_manifest import write_dataset_version_manifest
from analysis.ml.feature_metadata import enrich_ml_feature_metadata
from analysis.ml.feature_sets import get_feature_set
from analysis.ml.grouped_cv_stability import _binary_metrics, run_grouped_cv_stability
from analysis.ml.labels import binary_label_series
from analysis.ml.model_run_ledger import LEDGER_COLUMNS
from analysis.ml.score_scale_audit import write_score_scale_audit
from analysis.ml.spd_estimated_pk_exposure import run_spd_estimated_pk_exposure_model


DEFAULT_OBJECTIVES = {
    "spd_binding": "spd_binding_label",
    "addon_aware_activity": "combined_activity_ml_label",
}
DEFAULT_ABLATIONS = {
    "rdkit_only": "ligand_physchem_descriptors_no_qed",
    "shortcut_reduced": "ligand_physchem_shortcut_reduced",
    "pair_score_only": "spd_binding_pair_final_scores_only",
    "full": "spd_binding_pair_final_full_no_qed",
}
DEFAULT_GROUPS = ["drug_id", "chemical_cluster", "target_id", "target_family"]
FAILURE_GROUPS = [
    "drug_id",
    "target_id",
    "target_family",
    "label_source",
    "source_family",
    "scaffold_key",
    "chemical_cluster",
    "butina_cluster",
    "applicability_domain",
]


def _git_commit(repo_root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )
    value = result.stdout.strip()
    return value if result.returncode == 0 and value else "unknown"


def _valid_predictions(pred: pd.DataFrame, label_col: str) -> pd.DataFrame:
    work = pred.copy()
    work[label_col] = pd.to_numeric(work[label_col], errors="coerce")
    work["ml_prediction_score"] = pd.to_numeric(work["ml_prediction_score"], errors="coerce")
    return work.dropna(subset=[label_col, "ml_prediction_score"]).copy()


def _rank_error_cases(
    pred: pd.DataFrame,
    *,
    label_col: str,
    top_k: int | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    work = _valid_predictions(pred, label_col)
    work = work.sort_values("ml_prediction_score", ascending=False, kind="mergesort").reset_index(drop=True)
    work["score_rank"] = range(1, len(work) + 1)
    n_positive = int(work[label_col].eq(1).sum())
    selected_k = n_positive if top_k is None else int(top_k)
    selected_k = min(max(selected_k, 1), len(work)) if len(work) else 0
    work["rank_predicted_positive"] = work["score_rank"].le(selected_k)
    positive = work[label_col].eq(1)
    called = work["rank_predicted_positive"]
    work["error_type"] = "true_negative"
    work.loc[positive & called, "error_type"] = "true_positive"
    work.loc[~positive & called, "error_type"] = "false_positive"
    work.loc[positive & ~called, "error_type"] = "false_negative"
    cutoff = (
        float(work.loc[work["score_rank"].eq(selected_k), "ml_prediction_score"].iloc[0])
        if selected_k and len(work)
        else math.nan
    )
    summary = {
        "decision_policy": "prevalence_matched_top_k" if top_k is None else "explicit_top_k",
        "K": int(selected_k),
        "n": int(len(work)),
        "n_positive": n_positive,
        "score_cutoff_at_K": cutoff,
        "true_positive": int(work["error_type"].eq("true_positive").sum()),
        "false_positive": int(work["error_type"].eq("false_positive").sum()),
        "false_negative": int(work["error_type"].eq("false_negative").sum()),
        "true_negative": int(work["error_type"].eq("true_negative").sum()),
        "warning": "This is a rank operating point, not a 0.5 probability threshold.",
    }
    return work, summary


def _group_metric_rows(cases: pd.DataFrame, *, label_col: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    total_positive = int(cases[label_col].eq(1).sum())
    for group_col in FAILURE_GROUPS:
        if group_col not in cases.columns:
            continue
        for value, group in cases.groupby(group_col, dropna=False):
            labels = group[label_col].astype(int)
            scores = group["ml_prediction_score"].astype(float)
            metrics = _binary_metrics(labels, scores) if labels.nunique() > 1 else {}
            n_positive = int(labels.eq(1).sum())
            n_negative = int(labels.eq(0).sum())
            rows.append(
                {
                    "summary_group_col": group_col,
                    "summary_group_value": value,
                    "n": int(len(group)),
                    "n_positive": n_positive,
                    "n_negative": n_negative,
                    "positive_rate": float(labels.mean()) if len(group) else math.nan,
                    "positive_share": n_positive / total_positive if total_positive else 0.0,
                    "mean_score": float(scores.mean()),
                    "median_score": float(scores.median()),
                    "n_true_positive": int(group["error_type"].eq("true_positive").sum()),
                    "n_false_positive": int(group["error_type"].eq("false_positive").sum()),
                    "n_false_negative": int(group["error_type"].eq("false_negative").sum()),
                    "n_true_negative": int(group["error_type"].eq("true_negative").sum()),
                    "AUROC": metrics.get("AUROC"),
                    "AUPRC": metrics.get("AUPRC"),
                    "Brier": metrics.get("Brier"),
                    "ECE": metrics.get("ECE"),
                }
            )
    return pd.DataFrame(rows)


def _write_failure_outputs(
    pred: pd.DataFrame,
    *,
    label_col: str,
    out_dir: Path,
    top_k: int | None,
    top_errors: int = 100,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cases, summary = _rank_error_cases(pred, label_col=label_col, top_k=top_k)
    cases.to_csv(out_dir / "ranked_error_cases.csv", index=False)
    cases.loc[cases["error_type"].eq("false_positive")].head(top_errors).to_csv(
        out_dir / "top_false_positives.csv",
        index=False,
    )
    cases.loc[cases["error_type"].eq("false_negative")].sort_values(
        "ml_prediction_score",
        ascending=True,
    ).head(top_errors).to_csv(out_dir / "top_false_negatives.csv", index=False)
    group_summary = _group_metric_rows(cases, label_col=label_col)
    group_summary.to_csv(out_dir / "group_failure_summary.csv", index=False)
    (out_dir / "rank_decision_manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return group_summary, summary


def _write_score_audit(df: pd.DataFrame, out_dir: Path) -> dict[str, Any]:
    score_cols = [
        col
        for col in [
            "atlas_score",
            "consensus_score",
            "consensus_score_raw",
            "consensus_z_score",
            "z_selected",
            "SCORCH_score_used",
            "final_score",
            "banana_score",
            "banana_binding_probability",
            "banana_score_normalized",
            "banana_atlas_blend_score",
            "binding_expert_score",
            "atlas_binding_prior",
        ]
        if col in df.columns
    ]
    numeric = df[score_cols].apply(pd.to_numeric, errors="coerce")
    coverage = pd.DataFrame(
        [
            {
                "score": col,
                "n_present": int(numeric[col].notna().sum()),
                "coverage_fraction": float(numeric[col].notna().mean()),
                "min": numeric[col].min(),
                "max": numeric[col].max(),
                "mean": numeric[col].mean(),
                "std": numeric[col].std(),
            }
            for col in score_cols
        ]
    )
    coverage.to_csv(out_dir / "score_coverage.csv", index=False)

    correlation_rows: list[dict[str, Any]] = []
    for left_idx, left in enumerate(score_cols):
        for right in score_cols[left_idx + 1 :]:
            pair = numeric[[left, right]].dropna()
            if len(pair) < 3:
                continue
            correlation_rows.append(
                {
                    "score_a": left,
                    "score_b": right,
                    "n": int(len(pair)),
                    "pearson": float(pair[left].corr(pair[right], method="pearson")),
                    "spearman": float(pair[left].corr(pair[right], method="spearman")),
                }
            )
    correlations = pd.DataFrame(correlation_rows)
    correlations.to_csv(out_dir / "score_correlations.csv", index=False)

    source_rows: list[dict[str, Any]] = []
    for source_col in [
        "z_selected_source",
        "consensus_z_score_source",
        "atlas_binding_prior_source",
        "binding_expert_source",
        "banana_score_feature_source",
    ]:
        if source_col not in df.columns:
            continue
        counts = df[source_col].fillna("missing").astype(str).value_counts(dropna=False)
        for value, count in counts.items():
            source_rows.append(
                {
                    "source_col": source_col,
                    "source_value": value,
                    "n": int(count),
                    "fraction": float(count / len(df)) if len(df) else math.nan,
                }
            )
    pd.DataFrame(source_rows).to_csv(out_dir / "score_source_counts.csv", index=False)

    within_pdb: list[dict[str, Any]] = []
    consensus_col = (
        "consensus_z_score"
        if "consensus_z_score" in df.columns
        else "consensus_score"
    )
    required = {"pdb_id", consensus_col, "banana_score"}
    if required.issubset(df.columns):
        work = df[list(required)].copy().rename(columns={consensus_col: "consensus_feature"})
        work["consensus_feature"] = pd.to_numeric(
            work["consensus_feature"],
            errors="coerce",
        )
        work["banana_score"] = pd.to_numeric(work["banana_score"], errors="coerce")
        for pdb_id, group in work.dropna().groupby("pdb_id"):
            if (
                len(group) < 10
                or group["consensus_feature"].nunique() < 2
                or group["banana_score"].nunique() < 2
            ):
                continue
            within_pdb.append(
                {
                    "pdb_id": pdb_id,
                    "n": int(len(group)),
                    "consensus_feature": consensus_col,
                    "pearson": float(
                        group["consensus_feature"].corr(group["banana_score"])
                    ),
                    "spearman": float(
                        group["consensus_feature"].corr(
                            group["banana_score"],
                            method="spearman",
                        )
                    ),
                }
            )
    within_frame = pd.DataFrame(within_pdb)
    within_frame.to_csv(out_dir / "consensus_banana_within_pdb.csv", index=False)

    atlas_equals_z = None
    if {"atlas_score", "z_selected"}.issubset(numeric.columns):
        comparable = numeric[["atlas_score", "z_selected"]].dropna()
        atlas_equals_z = int(comparable["atlas_score"].eq(comparable["z_selected"]).sum())
    findings = {
        "n_rows": int(len(df)),
        "atlas_equals_z_selected_rows": atlas_equals_z,
        "z_selected_source_counts": (
            {str(key): int(value) for key, value in df["z_selected_source"].fillna("missing").value_counts().items()}
            if "z_selected_source" in df.columns
            else {}
        ),
        "within_pdb_consensus_banana_spearman_median": (
            float(within_frame["spearman"].median()) if not within_frame.empty else None
        ),
        "interpretation": [
            "consensus_score_raw is a deterministic within-PDB docking percentile, not an absolute potency.",
            "consensus_z_score uses an explicitly marked DUD-E consensus null and is the homogeneous model input.",
            "atlas_score, z_selected, and atlas_binding_prior alias consensus_z_score after repair; raw values remain in provenance columns.",
            "BANANA is an independent learned pocket-ligand logit; one-to-one agreement with docking is not expected.",
        ],
    }
    (out_dir / "score_provenance_findings.json").write_text(
        json.dumps(findings, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return findings


def _write_clozapine_audit(df: pd.DataFrame, out_dir: Path, label_col: str = "spd_binding_label") -> dict[str, Any]:
    if "drug_id" not in df.columns or label_col not in df.columns:
        return {"status": "skipped", "reason": "drug_id_or_label_missing"}
    work = df[["drug_id", label_col]].copy()
    work[label_col] = binary_label_series(work[label_col])
    work = work.loc[work[label_col].notna()].copy()
    grouped = (
        work.groupby("drug_id", dropna=False)[label_col]
        .agg(n="size", n_positive="sum", positive_rate="mean")
        .reset_index()
        .sort_values(["positive_rate", "n"], ascending=[False, False])
        .reset_index(drop=True)
    )
    grouped["positive_rate_rank"] = range(1, len(grouped) + 1)
    grouped.to_csv(out_dir / "drug_positive_rate_comparison.csv", index=False)
    clozapine = grouped.loc[grouped["drug_id"].astype(str).str.casefold().eq("clozapine")]
    result: dict[str, Any]
    if clozapine.empty:
        result = {"status": "skipped", "reason": "clozapine_not_found"}
    else:
        row = clozapine.iloc[0]
        result = {
            "status": "ok",
            "n": int(row["n"]),
            "n_positive": int(row["n_positive"]),
            "positive_rate": float(row["positive_rate"]),
            "positive_rate_rank": int(row["positive_rate_rank"]),
            "n_drugs": int(len(grouped)),
            "percentile": float(grouped["positive_rate"].rank(pct=True).loc[clozapine.index[0]]),
            "overall_positive_rate": float(work[label_col].astype(int).mean()),
            "median_drug_positive_rate": float(grouped["positive_rate"].median()),
        }
    (out_dir / "clozapine_positive_rate_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return result


def _write_label_source_profiles(
    df: pd.DataFrame,
    objectives: Mapping[str, str],
    out_dir: Path,
) -> None:
    rows: list[dict[str, Any]] = []
    for objective, label_col in objectives.items():
        if label_col not in df.columns:
            continue
        labels = binary_label_series(df[label_col])
        for group_col in ["label_source", "source_family", "scaffold_key", "target_family"]:
            if group_col not in df.columns:
                continue
            work = pd.DataFrame({group_col: df[group_col], "label": labels}).dropna(subset=["label"])
            total_positive = int(work["label"].eq(1).sum())
            for value, group in work.groupby(group_col, dropna=False):
                n_positive = int(group["label"].eq(1).sum())
                rows.append(
                    {
                        "objective": objective,
                        "label_col": label_col,
                        "group_col": group_col,
                        "group_value": value,
                        "n": int(len(group)),
                        "n_positive": n_positive,
                        "n_negative": int(group["label"].eq(0).sum()),
                        "positive_rate": float(group["label"].astype(int).mean()),
                        "positive_share": n_positive / total_positive if total_positive else 0.0,
                    }
                )
    pd.DataFrame(rows).to_csv(out_dir / "label_source_scaffold_family_profile.csv", index=False)


def _architecture_rows() -> pd.DataFrame:
    rdkit = ", ".join(get_feature_set("ligand_physchem_descriptors"))
    pair = ", ".join(get_feature_set("spd_binding_consensus_z_banana_scores_only"))
    full = ", ".join(get_feature_set("spd_binding_nonleaky_consensus_z"))
    rows: list[dict[str, str]] = [
        {
            "Expert / Model": "Decoy-standardized docking-consensus feature",
            "What it predicts": "Relative within-PDB multi-engine docking rank; not AC50 or probability",
            "Label used": "None (physics/rank-based docking)",
            "Inputs used": (
                "consensus_score_raw, pdb_id, comparison/reference run ID, "
                "explicitly marked DUD-E consensus_score values -> consensus_z_score"
            ),
        },
        {
            "Expert / Model": "BANANA binding prior",
            "What it predicts": "Learned pocket-ligand binding logit; uncalibrated",
            "Label used": "BANANA/BigBind training labels external to SPD",
            "Inputs used": "banana_smiles, pocket_pdb",
        },
    ]
    ablations = {
        "RDKit-only": rdkit,
        "pair-score-only": pair,
        "full": full,
    }
    for model_type in ("logistic regression", "LightGBM"):
        for ablation, inputs in ablations.items():
            rows.append(
                {
                    "Expert / Model": f"SPD binding {model_type} ({ablation})",
                    "What it predicts": "Held-out rank score for SPD binary binding",
                    "Label used": "spd_binding_label",
                    "Inputs used": inputs,
                }
            )
    for model_type in ("logistic regression", "LightGBM"):
        for ablation, inputs in ablations.items():
            rows.append(
                {
                    "Expert / Model": f"Addon-aware {model_type} ({ablation})",
                    "What it predicts": "Combined SPD plus external activity evidence",
                    "Label used": "combined_activity_ml_label",
                    "Inputs used": inputs,
                }
            )
    rows.extend(
        [
            {
                "Expert / Model": "Potency regressor",
                "What it predicts": "log10(AC50 uM) for exact-measurement pairs",
                "Label used": "spd_ac50_uM where spd_activity_relation == '='",
                "Inputs used": full,
            },
            {
                "Expert / Model": "Free-Cmax regressor",
                "What it predicts": "log10(free Cmax uM) at drug level",
                "Label used": "free_cmax_um",
                "Inputs used": rdkit,
            },
            {
                "Expert / Model": "Exposure-margin combiner",
                "What it predicts": (
                    "Uncalibrated Gaussian-residual estimate of P(predicted AC50 / "
                    "predicted free Cmax <= 10)"
                ),
                "Label used": "spd_exposure_relevant for held-out evaluation",
                "Inputs used": (
                    "predicted_log10_ac50_um, predicted_log10_free_cmax_um, "
                    "potency_residual_sigma, free_cmax_residual_sigma"
                ),
            },
        ]
    )
    return pd.DataFrame(rows)


def _build_data_priorities(group_summaries: pd.DataFrame, selected: pd.DataFrame) -> pd.DataFrame:
    if group_summaries.empty or selected.empty:
        return pd.DataFrame()
    keys = ["objective", "cv_group_col", "model_type", "ablation"]
    chosen = selected[keys].drop_duplicates()
    merged = group_summaries.merge(chosen, on=keys, how="inner")
    merged = merged.loc[
        merged["summary_group_col"].isin(
            ["drug_id", "target_id", "target_family", "scaffold_key", "chemical_cluster"]
        )
    ].copy()
    merged["false_negative_rate"] = merged["n_false_negative"] / merged["n_positive"].clip(lower=1)
    merged["false_positive_rate"] = merged["n_false_positive"] / merged["n_negative"].clip(lower=1)
    merged["positive_gap_priority"] = merged["false_negative_rate"] * merged["n_positive"].pow(0.5)
    merged["negative_control_priority"] = merged["false_positive_rate"] * merged["n_negative"].pow(0.5)
    merged["data_addition_priority"] = merged[["positive_gap_priority", "negative_control_priority"]].max(axis=1)
    merged["recommended_addition"] = "targeted positive measurements"
    merged.loc[
        merged["negative_control_priority"].gt(merged["positive_gap_priority"])
        | merged["n_negative"].eq(0),
        "recommended_addition",
    ] = "matched inactive or censored-negative controls"
    return merged.sort_values("data_addition_priority", ascending=False)


def _write_ledger(
    pooled: pd.DataFrame,
    *,
    dataset_path: Path,
    dataset_version: str,
    objectives: Mapping[str, str],
    out_dir: Path,
    run_id: str,
    git_commit: str,
    model_n_jobs: int,
) -> pd.DataFrame:
    now = datetime.now(timezone.utc).isoformat()
    rows: list[dict[str, Any]] = []
    for row in pooled.itertuples(index=False):
        label_col = objectives[str(row.objective)]
        record = {
            "run_id": run_id,
            "date": now,
            "git_commit": git_commit,
            "dataset_version": dataset_version,
            "dataset_path": str(dataset_path),
            "label_used": label_col,
            "feature_set": row.feature_set,
            "excluded_columns": "label-definition and leakage columns",
            "split_method": f"pooled_grouped_oof:{row.group_col}",
            "PU_strategy": "standard_binary",
            "model_type": row.model_type,
            "selected_model": False,
            "hyperparameters": json.dumps({"n_jobs": model_n_jobs}, sort_keys=True),
            "calibration_method": "none",
            "number_of_rows": int(row.n_test),
            "number_of_positives": int(row.n_test_positive),
            "n_train": None,
            "n_test": int(row.n_test),
            "AUROC": row.AUROC,
            "PR_AUC": row.AUPRC,
            "precision_at_K": row.precision_at_K,
            "precision_K": row.K,
            "enrichment_at_K": row.enrichment_at_K,
            "enrichment_K": row.K,
            "Brier_score": row.Brier,
            "notes": "Held-out pooled OOF metric; probabilities are not calibrated.",
            "model_dir": row.run_dir,
            "objective": row.objective,
            "ablation": row.ablation,
            "cv_group_col": row.group_col,
        }
        rows.append(record)
    columns = [*LEDGER_COLUMNS, "objective", "ablation", "cv_group_col"]
    ledger = pd.DataFrame(rows).reindex(columns=columns)
    ledger.to_csv(out_dir / "ml_model_run_ledger.csv", index=False)
    return ledger


def run_grouped_cv_ablation_suite(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    objectives: Mapping[str, str] | None = None,
    ablations: Mapping[str, str] | None = None,
    model_types: list[str] | None = None,
    group_cols: list[str] | None = None,
    n_splits: int = 5,
    repeats: int = 1,
    min_test_positives: int = 10,
    min_test_negatives: int = 10,
    top_k: int = 20,
    error_top_k: int | None = None,
    n_bootstraps: int = 200,
    n_permutations: int = 0,
    seed: int = 42,
    model_n_jobs: int = 4,
    compute_applicability_domain: bool = True,
    run_exposure: bool = True,
    exposure_models: list[str] | None = None,
    refresh_feature_metadata: bool = True,
    run_dir: str | Path | None = None,
    extra_feature_tables: list[str | Path] | None = None,
) -> dict[str, Any]:
    source_dataset = Path(dataset_path)
    dataset = source_dataset
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    objectives = dict(objectives or DEFAULT_OBJECTIVES)
    ablations = dict(ablations or DEFAULT_ABLATIONS)
    model_types = list(model_types or ["logistic_regression", "lightgbm"])
    group_cols = list(group_cols or DEFAULT_GROUPS)
    exposure_models = list(exposure_models or ["ridge", "lightgbm"])

    df = pd.read_csv(source_dataset, low_memory=False)
    feature_refresh_summary: dict[str, Any] | None = None
    if refresh_feature_metadata:
        repo_root = Path(__file__).resolve().parents[2]
        refreshed, feature_refresh_summary = enrich_ml_feature_metadata(
            df,
            chemical_cluster="auto",
            target_family="auto",
            source_lineage="auto",
            run_dir=Path(run_dir) if run_dir is not None else None,
            repo_root=repo_root,
            extra_feature_tables=[Path(value) for value in (extra_feature_tables or [])],
        )
        feature_dir = out / "_feature_metadata"
        feature_dir.mkdir(parents=True, exist_ok=True)
        dataset = feature_dir / source_dataset.name
        refreshed.to_csv(dataset, index=False)
        (feature_dir / "feature_metadata_summary.json").write_text(
            json.dumps(feature_refresh_summary, indent=2, sort_keys=True, default=str),
            encoding="utf-8",
        )
        df = refreshed
    missing_labels = [label for label in objectives.values() if label not in df.columns]
    if missing_labels:
        raise ValueError(f"dataset missing objective labels: {', '.join(missing_labels)}")
    for feature_set in ablations.values():
        missing = [feature for feature in get_feature_set(feature_set) if feature not in df.columns]
        if missing:
            raise ValueError(f"feature set {feature_set!r} missing columns: {', '.join(missing)}")

    dataset_manifest = write_dataset_version_manifest(
        dataset_path=dataset,
        frame=df,
        out_path=out / "dataset_version_manifest.json",
        label_col=next(iter(objectives.values())),
        feature_set="grouped_cv_ablation_suite",
        provenance={"objectives": objectives, "ablations": ablations},
    )
    score_findings = _write_score_audit(df, out)
    selected_score_features = sorted(
        {
            feature
            for feature_set in ablations.values()
            for feature in get_feature_set(feature_set)
        }
    )
    score_scale = write_score_scale_audit(
        df,
        out / "score_scale",
        feature_names=selected_score_features,
    )
    clozapine = _write_clozapine_audit(df, out)
    _write_label_source_profiles(df, objectives, out)
    architecture = _architecture_rows()
    architecture.to_csv(out / "model_architecture_table.csv", index=False)

    pooled_frames: list[pd.DataFrame] = []
    fold_frames: list[pd.DataFrame] = []
    bootstrap_frames: list[pd.DataFrame] = []
    prediction_index_rows: list[dict[str, Any]] = []
    group_summary_frames: list[pd.DataFrame] = []
    decision_rows: list[dict[str, Any]] = []

    for objective, label_col in objectives.items():
        for ablation, feature_set in ablations.items():
            run_dir = out / "grouped_cv" / objective / ablation
            result = run_grouped_cv_stability(
                dataset,
                label_col=label_col,
                feature_set=feature_set,
                out_dir=run_dir,
                model_types=model_types,
                group_cols=group_cols,
                leave_one_group_cols=[],
                n_splits=n_splits,
                repeats=repeats,
                min_test_positives=min_test_positives,
                min_test_negatives=min_test_negatives,
                top_k=top_k,
                n_bootstraps=n_bootstraps,
                n_permutations=n_permutations,
                seed=seed,
                class_weight="balanced",
                pu_mode="standard_binary",
                strict_feature_set=True,
                model_params={"n_jobs": model_n_jobs},
                compute_applicability_domain=compute_applicability_domain,
            )
            pooled = result["pooled_metrics"].copy()
            pooled["objective"] = objective
            pooled["label_col"] = label_col
            pooled["ablation"] = ablation
            pooled["feature_set"] = feature_set
            pooled["run_dir"] = str(run_dir)
            pooled_frames.append(pooled)
            folds = result["fold_metrics"].copy()
            folds["objective"] = objective
            folds["label_col"] = label_col
            folds["ablation"] = ablation
            folds["feature_set"] = feature_set
            fold_frames.append(folds)
            bootstrap_path = run_dir / "pooled_oof_bootstrap_ci.csv"
            if bootstrap_path.exists():
                bootstrap = pd.read_csv(bootstrap_path)
                bootstrap["objective"] = objective
                bootstrap["label_col"] = label_col
                bootstrap["ablation"] = ablation
                bootstrap["feature_set"] = feature_set
                bootstrap_frames.append(bootstrap)

            for metric_row in pooled.itertuples(index=False):
                pred_path = Path(str(metric_row.prediction_file))
                pred = pd.read_csv(pred_path, low_memory=False)
                failure_dir = run_dir / "failure_analysis" / str(metric_row.model_type) / str(metric_row.group_col)
                group_summary, decision = _write_failure_outputs(
                    pred,
                    label_col=label_col,
                    out_dir=failure_dir,
                    top_k=error_top_k,
                )
                if not group_summary.empty:
                    group_summary["objective"] = objective
                    group_summary["label_col"] = label_col
                    group_summary["ablation"] = ablation
                    group_summary["feature_set"] = feature_set
                    group_summary["model_type"] = metric_row.model_type
                    group_summary["cv_group_col"] = metric_row.group_col
                    group_summaries_path = failure_dir / "group_failure_summary.csv"
                    group_summary["source_file"] = str(group_summaries_path)
                    group_summary_frames.append(group_summary)
                decision_rows.append(
                    {
                        "objective": objective,
                        "label_col": label_col,
                        "ablation": ablation,
                        "feature_set": feature_set,
                        "model_type": metric_row.model_type,
                        "cv_group_col": metric_row.group_col,
                        **decision,
                    }
                )
                prediction_index_rows.append(
                    {
                        "objective": objective,
                        "label_col": label_col,
                        "ablation": ablation,
                        "feature_set": feature_set,
                        "model_type": metric_row.model_type,
                        "cv_group_col": metric_row.group_col,
                        "prediction_file": str(pred_path),
                        "failure_dir": str(failure_dir),
                    }
                )

    pooled_all = pd.concat(pooled_frames, ignore_index=True) if pooled_frames else pd.DataFrame()
    pooled_all.to_csv(out / "grouped_cv_ablation_pooled_metrics.csv", index=False)
    folds_all = pd.concat(fold_frames, ignore_index=True) if fold_frames else pd.DataFrame()
    folds_all.to_csv(out / "grouped_cv_ablation_fold_metrics.csv", index=False)
    bootstrap_all = pd.concat(bootstrap_frames, ignore_index=True) if bootstrap_frames else pd.DataFrame()
    bootstrap_all.to_csv(out / "grouped_cv_ablation_bootstrap_ci.csv", index=False)
    pd.DataFrame(prediction_index_rows).to_csv(out / "pooled_oof_prediction_index.csv", index=False)
    pd.DataFrame(decision_rows).to_csv(out / "rank_decision_summary.csv", index=False)
    group_all = pd.concat(group_summary_frames, ignore_index=True) if group_summary_frames else pd.DataFrame()
    group_all.to_csv(out / "all_group_failure_summaries.csv", index=False)

    full = pooled_all.loc[pooled_all["ablation"].eq("full")].copy()
    selected = (
        full.sort_values("AUPRC", ascending=False)
        .groupby(["objective", "group_col"], as_index=False, dropna=False)
        .head(1)
        .rename(columns={"group_col": "cv_group_col"})
    )
    selected.to_csv(out / "diagnostic_best_full_model_by_group.csv", index=False)
    priorities = _build_data_priorities(group_all, selected)
    priorities.to_csv(out / "data_addition_priorities.csv", index=False)

    exposure_rows: list[dict[str, Any]] = []
    if run_exposure and "spd_exposure_relevant" in df.columns:
        for model_type in exposure_models:
            exposure_dir = out / "exposure_margin" / model_type
            exposure_manifest = run_spd_estimated_pk_exposure_model(
                dataset,
                exposure_dir,
                label_col="spd_exposure_relevant",
                potency_feature_set="spd_potency_consensus_z",
                pk_feature_set="ligand_physchem_descriptors",
                split_mode="drug_holdout",
                model_type=model_type,
                seed=seed,
                censored_policy="exclude",
            )
            exposure_predictions = pd.read_csv(
                exposure_dir / "estimated_pk_exposure_predictions.csv",
                low_memory=False,
            )
            exposure_label = pd.to_numeric(
                exposure_predictions["spd_exposure_relevant"],
                errors="coerce",
            )
            exposure_rows.append(
                {
                    "model_type": model_type,
                    "out_dir": str(exposure_dir),
                    "n_test_positive": int(exposure_label.eq(1).sum()),
                    "n_test_negative": int(exposure_label.eq(0).sum()),
                    **exposure_manifest["metrics"],
                }
            )
    pd.DataFrame(exposure_rows).to_csv(out / "exposure_margin_model_comparison.csv", index=False)

    repo_root = Path(__file__).resolve().parents[2]
    git_commit = _git_commit(repo_root)
    ledger = _write_ledger(
        pooled_all,
        dataset_path=dataset,
        dataset_version=str(dataset_manifest["dataset_version_id"]),
        objectives=objectives,
        out_dir=out,
        run_id=out.name,
        git_commit=git_commit,
        model_n_jobs=model_n_jobs,
    )
    if exposure_rows:
        exposure_ledger_rows = []
        for row in exposure_rows:
            exposure_ledger_rows.append(
                {
                    "run_id": out.name,
                    "date": datetime.now(timezone.utc).isoformat(),
                    "git_commit": git_commit,
                    "dataset_version": dataset_manifest["dataset_version_id"],
                    "dataset_path": str(dataset),
                    "label_used": "spd_exposure_relevant",
                    "feature_set": "spd_potency_consensus_z + ligand_physchem_descriptors",
                    "excluded_columns": "observed free Cmax and label-definition fields",
                    "split_method": "drug_holdout",
                    "PU_strategy": "standard_binary",
                    "model_type": f"two_stage_{row['model_type']}",
                    "selected_model": False,
                    "hyperparameters": json.dumps({"censored_policy": "exclude"}),
                    "calibration_method": "none",
                    "number_of_rows": int(row["n_train_pairs"] + row["n_test_pairs"]),
                    "number_of_positives": None,
                    "n_train": int(row["n_train_pairs"]),
                    "n_test": int(row["n_test_pairs"]),
                    "AUROC": row.get("AUROC"),
                    "PR_AUC": row.get("AUPRC"),
                    "precision_at_K": None,
                    "precision_K": None,
                    "enrichment_at_K": None,
                    "enrichment_K": None,
                    "Brier_score": row.get("Brier"),
                    "notes": "Explicit predicted AC50 / predicted free-Cmax margin; right-censored AC50 excluded from potency fitting.",
                    "model_dir": row["out_dir"],
                    "objective": "exposure_margin",
                    "ablation": "atlas_prior_plus_predicted_pk",
                    "cv_group_col": "drug_id",
                }
            )
        ledger = pd.concat([ledger, pd.DataFrame(exposure_ledger_rows)], ignore_index=True, sort=False)
        ledger.to_csv(out / "ml_model_run_ledger.csv", index=False)
    readiness = {
        "status": "exploratory_not_claim_ready",
        "blockers": [
            "OOF probabilities are not cross-fitted calibrated probabilities; use ranking and Top-K metrics.",
            "The sparse external addon contributes positives but no matched negatives, so source transfer is not established.",
            "Model selection from this same OOF table is diagnostic only; freeze a model before a final untouched holdout.",
            "Consensus-z and BANANA show real within-PDB disagreement and require separate ablation reporting.",
        ],
        "strengths": [
            "All reported classifier scores are pooled predictions for rows excluded from the corresponding training fold.",
            "Drug, chemical-cluster, target, and target-family grouping axes are evaluated separately.",
            "False positives and false negatives use an explicit rank operating point rather than threshold 0.5.",
            "Potency regression excludes right-censored AC50 bounds by default.",
        ],
    }
    (out / "suite_claim_readiness.json").write_text(
        json.dumps(readiness, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    manifest = {
        "source_dataset": str(source_dataset),
        "dataset": str(dataset),
        "dataset_version": dataset_manifest["dataset_version_id"],
        "feature_refresh_summary": feature_refresh_summary,
        "out_dir": str(out),
        "objectives": objectives,
        "ablations": ablations,
        "model_types": model_types,
        "group_cols": group_cols,
        "n_splits": n_splits,
        "repeats": repeats,
        "compute_applicability_domain": compute_applicability_domain,
        "score_findings": score_findings,
        "score_scale": score_scale,
        "clozapine": clozapine,
        "n_pooled_metric_rows": int(len(pooled_all)),
        "n_ledger_rows": int(len(ledger)),
        "outputs": {
            "pooled_metrics": str(out / "grouped_cv_ablation_pooled_metrics.csv"),
            "bootstrap_ci": str(out / "grouped_cv_ablation_bootstrap_ci.csv"),
            "prediction_index": str(out / "pooled_oof_prediction_index.csv"),
            "failure_summaries": str(out / "all_group_failure_summaries.csv"),
            "data_priorities": str(out / "data_addition_priorities.csv"),
            "architecture": str(out / "model_architecture_table.csv"),
            "score_scale": score_scale["outputs"],
            "ledger": str(out / "ml_model_run_ledger.csv"),
        },
    }
    (out / "grouped_cv_ablation_suite_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return manifest
