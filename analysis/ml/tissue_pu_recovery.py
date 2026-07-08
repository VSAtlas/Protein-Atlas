from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.dataset_manifest import write_dataset_version_manifest
from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.splits import (
    CHEMICAL_CLUSTER_FIELDS,
    SOURCE_HOLDOUT_FIELDS,
    TARGET_FAMILY_FIELDS,
    make_split,
    split_overlap_summary,
)
from analysis.ml.train_classifier_core import _design_matrix
from analysis.ml.train_classifier_pu import _fit_bagging_pu


DEFAULT_CASE_COLUMNS = [
    "drug_id",
    "drug_name",
    "ligand_base",
    "target_id",
    "target_gene",
    "gene_symbol",
    "pdb_id",
    "adr_site_group",
    "adr_term",
    "label_source",
    "source_family",
    "protein_class",
    "target_family",
    "ligand_chemotype",
    "scaffold_key",
    "site_relevance_score",
    "site_relevance_score_by_adr",
    "expression_presence_score",
    "site_specificity_score",
    "expression_concordance_score",
]


def _git_commit(repo_root: str | Path | None) -> str:
    if repo_root is None:
        return "unknown"
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(repo_root),
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return "unknown"
    value = result.stdout.strip()
    return value if result.returncode == 0 and value else "unknown"


def _first_available_field(df: pd.DataFrame, fields: list[str]) -> str | None:
    for field in fields:
        if field in df.columns and df[field].notna().any():
            return field
    return None


def _context_fields(df: pd.DataFrame, split_mode: str) -> list[str]:
    fields_by_mode = {
        "drug_holdout": ["drug_id"],
        "target_holdout": ["target_id"],
        "scaffold_holdout": ["scaffold_key", "ligand_chemotype"],
        "chemical_cluster_holdout": CHEMICAL_CLUSTER_FIELDS,
        "target_family_holdout": TARGET_FAMILY_FIELDS,
        "protein_class_holdout": ["protein_class"],
        "source_holdout": SOURCE_HOLDOUT_FIELDS,
    }
    return [field for field in fields_by_mode.get(split_mode, []) if field in df.columns]


def _matches_context(pool: pd.DataFrame, heldout: pd.DataFrame, fields: list[str]) -> pd.Series:
    if pool.empty or heldout.empty or not fields:
        return pd.Series(False, index=pool.index)
    mask = pd.Series(False, index=pool.index)
    for field in fields:
        if field not in pool.columns or field not in heldout.columns:
            continue
        values = set(heldout[field].dropna().astype(str))
        if values:
            mask |= pool[field].fillna("").astype(str).isin(values)
    return mask


def _random_split_background(pool: pd.DataFrame, *, seed: int, test_fraction: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    if pool.empty:
        return pool.copy(), pool.copy()
    n_eval = max(1, int(round(len(pool) * float(test_fraction))))
    eval_idx = pool.sample(n=min(n_eval, len(pool)), random_state=seed).index
    return pool.drop(index=eval_idx).copy(), pool.loc[eval_idx].copy()


def _bounded_eval_background(
    frame: pd.DataFrame,
    *,
    max_eval_background: int,
    seed: int,
) -> pd.DataFrame:
    if max_eval_background <= 0 or len(frame) <= max_eval_background:
        return frame
    return frame.sample(n=max_eval_background, random_state=seed).copy()


def _split_tissue_pu_frames(
    df: pd.DataFrame,
    *,
    label_col: str,
    split_mode: str,
    seed: int,
    test_fraction: float,
    max_eval_background: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    labels = binary_label_series(df[label_col])
    work = df.copy()
    work["_tissue_pu_label_state"] = labels
    positives = work.loc[labels.eq(1)].copy()
    negatives = work.loc[labels.eq(0)].copy()
    unknown = work.loc[labels.isna()].copy()
    if positives.empty:
        raise ValueError(f"{label_col} has no positive rows for tissue PU recovery")
    train_pos_idx, test_pos_idx = make_split(positives, split_mode, seed=seed, test_fraction=test_fraction)
    train_pos = positives.loc[train_pos_idx].copy()
    eval_pos = positives.loc[test_pos_idx].copy()
    if train_pos.empty or eval_pos.empty:
        raise ValueError(
            f"{split_mode} produced empty positive train/test for tissue PU recovery: "
            f"train={len(train_pos)} test={len(eval_pos)}"
        )
    fields = _context_fields(work, split_mode)
    background_status = "heldout_context"
    if split_mode == "random" or not fields:
        train_unknown, eval_unknown = _random_split_background(
            unknown,
            seed=seed,
            test_fraction=test_fraction,
        )
        train_neg, eval_neg = _random_split_background(
            negatives,
            seed=seed,
            test_fraction=test_fraction,
        )
        background_status = "random_background_split" if split_mode == "random" else "fallback_random_unknown"
    else:
        unknown_eval_mask = _matches_context(unknown, eval_pos, fields)
        neg_eval_mask = _matches_context(negatives, eval_pos, fields)
        eval_unknown = unknown.loc[unknown_eval_mask].copy()
        train_unknown = unknown.loc[~unknown_eval_mask].copy()
        eval_neg = negatives.loc[neg_eval_mask].copy()
        train_neg = negatives.loc[~neg_eval_mask].copy()
        if eval_unknown.empty and not unknown.empty:
            train_unknown, eval_unknown = _random_split_background(
                unknown,
                seed=seed,
                test_fraction=test_fraction,
            )
            background_status = "fallback_random_unknown"
    eval_unknown = _bounded_eval_background(
        eval_unknown,
        max_eval_background=max_eval_background,
        seed=seed,
    )
    eval_frame = pd.concat([eval_pos, eval_unknown, eval_neg], ignore_index=False)
    eval_status = pd.Series("unknown_background", index=eval_frame.index, dtype="object")
    eval_status.loc[eval_pos.index] = "heldout_positive"
    if not eval_neg.empty:
        eval_status.loc[eval_neg.index] = "heldout_measured_negative"
    eval_frame["tissue_pu_eval_status"] = eval_status
    eval_frame["tissue_pu_recovery_label"] = eval_status.eq("heldout_positive").astype(int)
    train_labeled = pd.concat([train_pos, train_neg], ignore_index=False)
    manifest = {
        "split_mode": split_mode,
        "test_fraction": float(test_fraction),
        "context_fields": fields,
        "background_selection_status": background_status,
        "n_total": int(len(work)),
        "n_positive_total": int(len(positives)),
        "n_measured_negative_total": int(len(negatives)),
        "n_unknown_total": int(len(unknown)),
        "n_train_positive": int(len(train_pos)),
        "n_train_measured_negative": int(len(train_neg)),
        "n_train_unlabeled_pool": int(len(train_unknown)),
        "n_eval_positive": int(len(eval_pos)),
        "n_eval_unknown_background": int(len(eval_unknown)),
        "n_eval_measured_negative": int(len(eval_neg)),
    }
    return train_labeled, train_unknown, eval_frame, manifest


def _write_tissue_split_manifest(
    *,
    train_labeled: pd.DataFrame,
    train_unknown: pd.DataFrame,
    eval_frame: pd.DataFrame,
    out_path: Path,
    label_col: str,
    split_mode: str,
) -> None:
    id_cols = [
        col
        for col in [
            "canonical_pair_key",
            "drug_id",
            "target_id",
            "pdb_id",
            "scaffold_key",
            "chemical_cluster",
            "target_family",
            "protein_class",
            "label_source",
            "source_family",
            label_col,
        ]
        if col in train_labeled.columns or col in train_unknown.columns or col in eval_frame.columns
    ]
    pieces: list[pd.DataFrame] = []
    for split_name, frame in [
        ("train_labeled", train_labeled),
        ("train_unlabeled_pool", train_unknown),
        ("evaluation", eval_frame),
    ]:
        table = frame.reindex(columns=id_cols).copy()
        table["_row_index"] = frame.index.astype(str)
        table["split"] = split_name
        table["fold"] = split_mode
        if "tissue_pu_eval_status" in frame.columns:
            table["tissue_pu_eval_status"] = frame["tissue_pu_eval_status"]
        pieces.append(table)
    manifest = pd.concat(pieces, ignore_index=True)
    manifest.to_csv(out_path, index=False)


def _topk_recovery_rows(
    pred: pd.DataFrame,
    *,
    score_col: str,
    top_k: list[int],
    group_cols: list[str] | None = None,
) -> list[dict[str, Any]]:
    group_cols = list(group_cols or [])
    rows: list[dict[str, Any]] = []

    def add_rows(scope: str, value: object, group: pd.DataFrame) -> None:
        work = group.dropna(subset=[score_col]).copy()
        if work.empty:
            return
        work["tissue_pu_recovery_label"] = pd.to_numeric(
            work["tissue_pu_recovery_label"],
            errors="coerce",
        ).fillna(0)
        n_pos = int(work["tissue_pu_recovery_label"].sum())
        if n_pos <= 0:
            return
        prevalence = n_pos / len(work)
        ranked = work.sort_values(score_col, ascending=False)
        for k in top_k:
            top_n = min(max(1, int(k)), len(ranked))
            top = ranked.head(top_n)
            hits = int(top["tissue_pu_recovery_label"].sum())
            fraction = hits / top_n if top_n else 0.0
            rows.append(
                {
                    "scope": scope,
                    "scope_value": value,
                    "score_col": score_col,
                    "top_k": int(top_n),
                    "n_eval": int(len(work)),
                    "n_eval_positive": n_pos,
                    "background_positive_rate": float(prevalence),
                    "positive_hits_at_k": hits,
                    "positive_fraction_at_k": float(fraction),
                    "heldout_positive_recovery_at_k": float(hits / n_pos),
                    "enrichment_over_background": float(fraction / prevalence) if prevalence else None,
                    "unknown_background_at_k": int(top["tissue_pu_eval_status"].eq("unknown_background").sum()),
                    "measured_negative_at_k": int(top["tissue_pu_eval_status"].eq("heldout_measured_negative").sum()),
                }
            )

    add_rows("all", "all", pred)
    for group_col in group_cols:
        if group_col not in pred.columns:
            continue
        for value, group in pred.groupby(group_col, dropna=False):
            add_rows(group_col, value, group)
    return rows


def _case_studies(pred: pd.DataFrame, *, score_col: str, rows_per_group: int) -> pd.DataFrame:
    keep = [col for col in DEFAULT_CASE_COLUMNS if col in pred.columns]
    keep.extend([score_col, "tissue_pu_rank", "tissue_pu_eval_status", "tissue_pu_recovery_label"])
    pieces: list[pd.DataFrame] = []
    ranked = pred.sort_values(score_col, ascending=False).copy()
    for status in ["heldout_positive", "unknown_background", "heldout_measured_negative"]:
        subset = ranked.loc[ranked["tissue_pu_eval_status"].eq(status)]
        if subset.empty:
            continue
        table = subset.head(rows_per_group).reindex(columns=keep).copy()
        table.insert(0, "case_study_type", f"top_ranked_{status}")
        pieces.append(table)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=["case_study_type", *keep])


def run_tissue_pu_recovery_report(
    dataset_path: str | Path,
    *,
    label_col: str = "tissue_site_label",
    feature_set: str = "spd_tissue_site_nonleaky",
    model_type: str = "logistic_regression",
    split_mode: str = "target_holdout",
    out_dir: str | Path,
    seed: int = 42,
    pu_mode: str = "stratified_bagging_pu",
    pu_bags: int = 50,
    pu_unlabeled_ratio: float = 1.0,
    top_k: list[int] | None = None,
    test_fraction: float = 0.2,
    max_eval_background: int = 10000,
    class_weight: str | None = "balanced",
    allow_tissue_expression_features: bool = True,
    strict_feature_set: bool = False,
    case_study_rows: int = 25,
    repo_root: str | Path | None = ".",
) -> dict[str, Any]:
    if pu_mode not in {"bagging_pu", "stratified_bagging_pu"}:
        raise ValueError("tissue PU recovery supports bagging_pu or stratified_bagging_pu")
    top_k = sorted({int(k) for k in (top_k or [10, 20, 50]) if int(k) > 0})
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(dataset_path, low_memory=False)
    if label_col not in df.columns:
        raise ValueError(f"missing tissue label column: {label_col}")
    requested = get_feature_set(feature_set)
    excluded = effective_exclude_features(
        label_col,
        None,
        allow_label_definition_features=allow_tissue_expression_features,
    )
    features = [feature for feature in requested if feature in df.columns and feature not in excluded]
    missing = [feature for feature in requested if feature not in df.columns]
    if strict_feature_set and missing:
        raise ValueError(f"missing requested tissue feature columns: {', '.join(missing)}")
    if not features:
        raise ValueError(f"no usable features for tissue PU recovery feature_set={feature_set}")
    assert_no_leakage(features)
    train_labeled, train_unknown, eval_frame, split_manifest = _split_tissue_pu_frames(
        df,
        label_col=label_col,
        split_mode=split_mode,
        seed=seed,
        test_fraction=test_fraction,
        max_eval_background=max_eval_background,
    )
    if train_unknown.empty:
        raise ValueError("tissue PU recovery requires fold-local unlabeled/background rows")
    train_labeled = train_labeled.copy()
    train_labeled[label_col] = binary_label_series(train_labeled[label_col])
    x_eval = _design_matrix(eval_frame, features)
    _, probabilities, pu_manifest = _fit_bagging_pu(
        model_type=model_type,
        seed=seed,
        class_weight=class_weight,
        model_train=train_labeled,
        unlabeled_pool=train_unknown,
        features=features,
        label_col=label_col,
        x_test=x_eval,
        out_path=out_path,
        n_bags=pu_bags,
        unlabeled_ratio=pu_unlabeled_ratio,
        stratified=pu_mode == "stratified_bagging_pu",
    )
    pred = eval_frame.copy()
    pred["tissue_pu_score"] = probabilities
    pred = pred.sort_values("tissue_pu_score", ascending=False)
    pred["tissue_pu_rank"] = range(1, len(pred) + 1)
    pred.to_csv(out_path / "tissue_pu_predictions.csv", index=False)
    topk_rows = _topk_recovery_rows(
        pred,
        score_col="tissue_pu_score",
        top_k=top_k,
        group_cols=["adr_site_group", "target_id", "pdb_id", "target_family", "protein_class"],
    )
    topk = pd.DataFrame(topk_rows)
    topk.to_csv(out_path / "tissue_pu_topk_recovery.csv", index=False)
    cases = _case_studies(pred, score_col="tissue_pu_score", rows_per_group=case_study_rows)
    cases.to_csv(out_path / "tissue_pu_case_studies.csv", index=False)
    _write_tissue_split_manifest(
        train_labeled=train_labeled,
        train_unknown=train_unknown,
        eval_frame=eval_frame,
        out_path=out_path / "tissue_pu_split_manifest.csv",
        label_col=label_col,
        split_mode=split_mode,
    )
    split_summary = split_overlap_summary(train_labeled, eval_frame.loc[pred.index], split_mode)
    dataset_manifest = write_dataset_version_manifest(
        dataset_path=dataset_path,
        frame=df,
        out_path=out_path / "dataset_version_manifest.json",
        label_col=label_col,
        feature_set=feature_set,
        features=features,
        requested_features=requested,
        missing_features=missing,
        exclude_features=sorted(excluded),
        provenance={
            "task": "tissue_pu_recovery",
            "allow_tissue_expression_features": bool(allow_tissue_expression_features),
        },
    )
    primary_k = 20 if 20 in top_k else top_k[0]
    primary_row = {}
    if not topk.empty:
        match = topk.loc[topk["scope"].eq("all") & topk["top_k"].eq(min(primary_k, len(pred)))]
        if not match.empty:
            primary_row = match.iloc[0].to_dict()
    record = {
        "run_id": None,
        "date": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(repo_root),
        "dataset_version": dataset_manifest.get("dataset_version_id"),
        "dataset_path": str(dataset_path),
        "label_used": label_col,
        "feature_set": feature_set,
        "excluded_columns": ";".join(sorted(excluded)),
        "split_method": split_mode,
        "PU_strategy": pu_mode,
        "model_type": model_type,
        "selected_model": model_type,
        "hyperparameters": json.dumps(
            {
                "pu_bags": int(pu_bags),
                "pu_unlabeled_ratio": float(pu_unlabeled_ratio),
                "class_weight": class_weight,
            },
            sort_keys=True,
        ),
        "calibration_method": "none",
        "number_of_rows": int(len(df)),
        "number_of_positives": int(binary_label_series(df[label_col]).eq(1).sum()),
        "n_train": int(len(train_labeled) + len(train_unknown)),
        "n_test": int(len(eval_frame)),
        "AUROC": None,
        "PR_AUC": None,
        "precision_at_K": primary_row.get("positive_fraction_at_k"),
        "precision_K": primary_row.get("top_k", primary_k),
        "enrichment_at_K": primary_row.get("enrichment_over_background"),
        "enrichment_K": primary_row.get("top_k", primary_k),
        "Brier_score": None,
        "notes": (
            "Tissue PU recovery report: positives are held out and ranked against unknown/background rows. "
            "Ordinary AUROC/AUPRC are not primary claims without independent tissue negatives."
        ),
        "model_dir": str(out_path),
        "task": "tissue_pu_recovery",
        "heldout_positive_recovery_at_K": primary_row.get("heldout_positive_recovery_at_k"),
    }
    (out_path / "model_run_record.json").write_text(json.dumps(record, indent=2, sort_keys=True, default=str), encoding="utf-8")
    pd.DataFrame([record]).to_csv(out_path / "model_run_record.csv", index=False)
    manifest = {
        "status": "ok",
        "dataset_path": str(dataset_path),
        "label_col": label_col,
        "feature_set": feature_set,
        "features": features,
        "missing_requested_features": missing,
        "split": split_manifest,
        "split_overlap_summary": split_summary,
        "pu_training": pu_manifest,
        "topk_recovery": str(out_path / "tissue_pu_topk_recovery.csv"),
        "predictions": str(out_path / "tissue_pu_predictions.csv"),
        "case_studies": str(out_path / "tissue_pu_case_studies.csv"),
        "policy": {
            "unlabeled_handling": "unknown rows are sampled only as temporary fold-local background inside PU bags",
            "primary_metrics": [
                "top-K positive recovery",
                "enrichment over background",
                "held-out-positive recovery",
                "case studies",
            ],
            "not_primary_metrics": ["AUROC", "AUPRC"],
            "expression_feature_policy": (
                "Expression features are allowed here only when tissue_site_label is externally sourced; "
                "rule-derived labels must be named weak/rule labels and interpreted as evidence-layer checks."
            ),
        },
    }
    (out_path / "tissue_pu_recovery_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return {
        "status": "ok",
        "manifest": manifest,
        "record": record,
        "topk_rows": topk_rows,
        "n_predictions": int(len(pred)),
    }
