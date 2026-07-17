from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.dataset_manifest import write_dataset_version_manifest
from analysis.ml.feature_sets import PHYSICHEM_DESCRIPTOR_FEATURES
from analysis.ml.grouped_cv_stability import run_grouped_cv_stability
from analysis.ml.labels import binary_label_series
from analysis.ml.model_run_ledger import LEDGER_COLUMNS


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


def _association_tables(
    frame: pd.DataFrame,
    *,
    label_col: str,
    descriptors: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = frame.copy()
    work[label_col] = binary_label_series(work[label_col])
    work = work.loc[work[label_col].notna()].copy()
    numeric = work[descriptors].apply(pd.to_numeric, errors="coerce")
    correlation = numeric.corr(method="spearman")
    correlation.index.name = "descriptor"
    correlation = correlation.reset_index()

    rows: list[dict[str, Any]] = []
    drug_cols = ["drug_id", label_col, *descriptors]
    drug_level = (
        work[drug_cols]
        .groupby("drug_id", dropna=False)
        .agg(
            {
                label_col: ["mean", "sum", "count"],
                **{descriptor: "median" for descriptor in descriptors},
            }
        )
    )
    drug_level.columns = [
        "drug_positive_rate",
        "drug_n_positive",
        "drug_n_rows",
        *descriptors,
    ]
    for descriptor in descriptors:
        row_values = pd.to_numeric(work[descriptor], errors="coerce")
        drug_values = pd.to_numeric(drug_level[descriptor], errors="coerce")
        rows.append(
            {
                "descriptor": descriptor,
                "row_label_spearman": float(
                    row_values.corr(work[label_col].astype(float), method="spearman")
                ),
                "drug_positive_rate_spearman": float(
                    drug_values.corr(
                        drug_level["drug_positive_rate"],
                        method="spearman",
                    )
                ),
                "drug_n_positive_spearman": float(
                    drug_values.corr(
                        drug_level["drug_n_positive"],
                        method="spearman",
                    )
                ),
                "n_labeled_rows": int(row_values.notna().sum()),
                "n_drugs": int(drug_values.notna().sum()),
            }
        )
    return pd.DataFrame(rows), correlation


def _ledger_rows(
    pooled: pd.DataFrame,
    *,
    dataset: Path,
    dataset_version: str,
    label_col: str,
    feature_set: str,
    out: Path,
    model_n_jobs: int,
) -> pd.DataFrame:
    now = datetime.now(timezone.utc).isoformat()
    commit = _git_commit(Path(__file__).resolve().parents[2])
    rows: list[dict[str, Any]] = []
    for row in pooled.itertuples(index=False):
        excluded = "" if row.excluded_descriptor == "none" else row.excluded_descriptor
        rows.append(
            {
                "run_id": out.name,
                "date": now,
                "git_commit": commit,
                "dataset_version": dataset_version,
                "dataset_path": str(dataset),
                "label_used": label_col,
                "feature_set": feature_set,
                "excluded_columns": excluded,
                "split_method": f"pooled_grouped_oof:{row.group_col}",
                "PU_strategy": "standard_binary",
                "model_type": row.model_type,
                "selected_model": False,
                "hyperparameters": json.dumps(
                    {"n_jobs": model_n_jobs, "class_weight": "balanced"},
                    sort_keys=True,
                ),
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
                "notes": (
                    "Descriptor leave-one-feature-out shortcut audit; "
                    "held-out pooled OOF probabilities are uncalibrated."
                ),
                "model_dir": row.run_dir,
            }
        )
    return pd.DataFrame(rows).reindex(columns=LEDGER_COLUMNS)


def run_descriptor_shortcut_audit(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    label_col: str = "spd_binding_label",
    feature_set: str = "spd_binding_nonleaky_consensus_z",
    model_type: str = "lightgbm",
    group_cols: list[str] | None = None,
    n_splits: int = 5,
    n_bootstraps: int = 100,
    seed: int = 42,
    model_n_jobs: int = 8,
) -> dict[str, Any]:
    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    groups = group_cols or [
        "drug_id",
        "chemical_cluster",
        "target_id",
        "target_family",
    ]
    frame = pd.read_csv(dataset, low_memory=False)
    descriptors = [
        descriptor
        for descriptor in PHYSICHEM_DESCRIPTOR_FEATURES
        if descriptor in frame.columns
    ]
    missing = sorted(set(PHYSICHEM_DESCRIPTOR_FEATURES) - set(descriptors))
    if missing:
        raise ValueError(f"missing RDKit descriptors: {', '.join(missing)}")

    version = write_dataset_version_manifest(
        dataset_path=dataset,
        frame=frame,
        out_path=out / "dataset_version_manifest.json",
        label_col=label_col,
        feature_set=feature_set,
        provenance={
            "audit": "descriptor_leave_one_feature_out",
            "group_cols": groups,
        },
    )
    associations, descriptor_correlation = _association_tables(
        frame,
        label_col=label_col,
        descriptors=descriptors,
    )
    associations.to_csv(out / "descriptor_label_associations.csv", index=False)
    descriptor_correlation.to_csv(
        out / "descriptor_spearman_correlation.csv",
        index=False,
    )

    pooled_frames: list[pd.DataFrame] = []
    conditions: list[tuple[str, list[str]]] = [("none", [])]
    conditions.extend((descriptor, [descriptor]) for descriptor in descriptors)
    for condition, excluded in conditions:
        run_dir = out / "grouped_cv" / f"without_{condition}"
        result = run_grouped_cv_stability(
            dataset,
            label_col=label_col,
            feature_set=feature_set,
            out_dir=run_dir,
            model_types=[model_type],
            group_cols=groups,
            leave_one_group_cols=[],
            n_splits=n_splits,
            repeats=1,
            min_test_positives=10,
            min_test_negatives=10,
            top_k=20,
            n_bootstraps=n_bootstraps,
            n_permutations=0,
            seed=seed,
            class_weight="balanced",
            pu_mode="standard_binary",
            exclude_features=excluded,
            strict_feature_set=True,
            model_params={"n_jobs": model_n_jobs},
            compute_applicability_domain=False,
        )
        pooled = result["pooled_metrics"].copy()
        pooled["excluded_descriptor"] = condition
        pooled["run_dir"] = str(run_dir)
        pooled_frames.append(pooled)

    pooled_all = pd.concat(pooled_frames, ignore_index=True)
    pooled_all.to_csv(out / "descriptor_leave_one_out_metrics.csv", index=False)
    baseline = pooled_all.loc[
        pooled_all["excluded_descriptor"].eq("none"),
        ["model_type", "evaluation", "group_col", "AUPRC", "AUROC", "Brier", "ECE"],
    ].rename(
        columns={
            "AUPRC": "baseline_AUPRC",
            "AUROC": "baseline_AUROC",
            "Brier": "baseline_Brier",
            "ECE": "baseline_ECE",
        }
    )
    deltas = pooled_all.loc[~pooled_all["excluded_descriptor"].eq("none")].merge(
        baseline,
        on=["model_type", "evaluation", "group_col"],
        how="left",
        validate="many_to_one",
    )
    for metric in ("AUPRC", "AUROC", "Brier", "ECE"):
        deltas[f"removal_delta_{metric}"] = (
            deltas[metric] - deltas[f"baseline_{metric}"]
        )
    deltas.to_csv(out / "descriptor_leave_one_out_deltas.csv", index=False)

    pivot = deltas.pivot_table(
        index="excluded_descriptor",
        columns="group_col",
        values="removal_delta_AUPRC",
        aggfunc="mean",
    ).reset_index()
    pivot.columns.name = None
    ranking = associations.rename(
        columns={"descriptor": "excluded_descriptor"}
    ).merge(pivot, on="excluded_descriptor", how="left")
    transfer_cols = [
        col
        for col in ("drug_id", "chemical_cluster", "target_id", "target_family")
        if col in ranking.columns
    ]
    positive_transfer = ranking[transfer_cols].clip(lower=0.0)
    ranking["mean_transfer_gain_when_removed"] = positive_transfer.mean(axis=1)
    ranking["max_transfer_gain_when_removed"] = positive_transfer.max(axis=1)
    ranking["shortcut_rank_score"] = (
        ranking["mean_transfer_gain_when_removed"]
        + ranking["drug_positive_rate_spearman"].abs() * 0.01
    )
    ranking = ranking.sort_values(
        ["shortcut_rank_score", "max_transfer_gain_when_removed"],
        ascending=False,
    )
    ranking.to_csv(out / "descriptor_shortcut_ranking.csv", index=False)

    ledger = _ledger_rows(
        pooled_all,
        dataset=dataset,
        dataset_version=str(version["dataset_version_id"]),
        label_col=label_col,
        feature_set=feature_set,
        out=out,
        model_n_jobs=model_n_jobs,
    )
    ledger.to_csv(out / "ml_model_run_ledger.csv", index=False)
    manifest = {
        "dataset": str(dataset),
        "dataset_version": version["dataset_version_id"],
        "out_dir": str(out),
        "label_col": label_col,
        "feature_set": feature_set,
        "model_type": model_type,
        "group_cols": groups,
        "descriptors": descriptors,
        "n_splits": n_splits,
        "n_bootstraps": n_bootstraps,
        "interpretation": (
            "Positive removal_delta_AUPRC means transfer improved when the "
            "descriptor was removed. Correlated descriptors can substitute for "
            "one another, so leave-one-out deltas are lower bounds on shortcut reliance."
        ),
        "outputs": {
            "metrics": str(out / "descriptor_leave_one_out_metrics.csv"),
            "deltas": str(out / "descriptor_leave_one_out_deltas.csv"),
            "associations": str(out / "descriptor_label_associations.csv"),
            "correlations": str(out / "descriptor_spearman_correlation.csv"),
            "ranking": str(out / "descriptor_shortcut_ranking.csv"),
            "ledger": str(out / "ml_model_run_ledger.csv"),
        },
    }
    (out / "descriptor_shortcut_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
