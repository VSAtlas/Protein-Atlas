from __future__ import annotations

import gc
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.audit_utils import (
    OVERLAP_FIELDS,
    duplicate_label_conflicts,
    ensure_canonical_pair_key,
    feature_missingness,
    grouped_label_balance,
    high_risk_column_scan,
    label_evidence_composition,
    label_balance,
    load_table,
    overlap_row,
    source_holdout,
    source_label_balance,
)
from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.leakage_checks import find_leaky_features
from analysis.ml.split_manifest import write_split_manifest
from analysis.ml.splits import make_split, split_overlap_summary


DEFAULT_SPLITS = [
    "random",
    "drug_holdout",
    "target_holdout",
    "scaffold_holdout",
    "chemical_cluster_holdout",
    "target_family_holdout",
    "temporal_holdout",
]

AUDIT_CONTEXT_COLUMNS = {
    "canonical_pair_key",
    "drug_id",
    "target_id",
    "pdb_id",
    "ligand_chemotype",
    "scaffold_key",
    "chemical_cluster",
    "chemical_cluster_id",
    "ecfp_cluster",
    "butina_cluster",
    "umap_cluster",
    "scaffold_cluster",
    "ligand_cluster",
    "target_family",
    "protein_family",
    "target_class",
    "protein_class",
    "dedup_drug_key",
    "dedup_target_key",
    "label_source",
    "source_family",
    "upstream_source",
    "assay_type",
    "assay_mode",
    "endpoint_type",
    "activity_type",
    "standard_type",
    "relation_domain",
    "source_objective",
    "benchmark_only",
    "_sample_weight",
    "negative_evidence_type",
    "label_publication_year",
    "evidence_publication_year",
    "activity_publication_year",
    "database_release_year",
    "smiles",
    "canonical_smiles",
    "ligand_smiles",
}


def audit_ml_leakage_overlap(
    dataset_path: str | Path,
    label_col: str,
    out_dir: str | Path,
    *,
    feature_set: str | None = None,
    exclude_features: list[str] | None = None,
    split_modes: list[str] | None = None,
    source_col: str = "label_source",
    seed: int = 42,
    test_fraction: float = 0.2,
    allow_partial_rescoring_features: bool = False,
    allow_label_definition_features: bool = False,
) -> dict[str, Any]:
    df = load_table(dataset_path)
    all_columns = list(df.columns)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    features: list[str] = []
    if feature_set:
        excluded = effective_exclude_features(
            label_col,
            exclude_features,
            allow_label_definition_features=allow_label_definition_features,
        )
        features = [
            feature
            for feature in get_feature_set(feature_set)
            if feature in df.columns and feature not in excluded
        ]
    projection = [
        column
        for column in all_columns
        if column in ({label_col, source_col, *features} | AUDIT_CONTEXT_COLUMNS)
    ]
    data = ensure_canonical_pair_key(df.loc[:, projection].copy())
    del df
    gc.collect()
    leaky_features = find_leaky_features(
        features,
        allow_label_definition_features=allow_label_definition_features,
        allow_partial_rescoring_features=allow_partial_rescoring_features,
    )

    labeled = data.dropna(subset=[label_col]).copy() if label_col in data.columns else data.iloc[0:0].copy()
    split_rows: list[dict[str, Any]] = []
    overlap_rows: list[dict[str, Any]] = []
    label_rows: list[dict[str, Any]] = []
    modes = split_modes or list(DEFAULT_SPLITS)
    for mode in modes:
        try:
            if mode == "source_holdout":
                train_idx, test_idx, split_summary = source_holdout(labeled, source_col, seed, test_fraction)
            else:
                train_idx, test_idx = make_split(labeled, split_mode=mode, seed=seed, test_fraction=test_fraction)
                split_summary = split_overlap_summary(labeled.loc[train_idx], labeled.loc[test_idx], mode)
            train = labeled.loc[train_idx]
            test = labeled.loc[test_idx]
            split_rows.append(
                {
                    "split_mode": mode,
                    "status": "ok",
                    "n_train": int(len(train)),
                    "n_test": int(len(test)),
                    "passes_holdout": bool(split_summary.get("passes_holdout", False)),
                    "held_out_sources": ";".join(split_summary.get("held_out_sources", [])),
                    "error": "",
                }
            )
            for field in OVERLAP_FIELDS:
                overlap_rows.append(overlap_row(train, test, field, mode))
            for part_name, part in (("train", train), ("test", test)):
                balance = label_balance(part, label_col)
                label_rows.append({"split_mode": mode, "partition": part_name, **balance})
            write_split_manifest(
                labeled,
                train_idx,
                test_idx,
                out / "split_manifests" / mode,
                label_col=label_col,
                split_mode=mode,
                split_summary=split_summary,
            )
        except Exception as exc:
            message = str(exc)
            status = "skipped" if "requires" in message else "error"
            split_rows.append(
                {
                    "split_mode": mode,
                    "status": status,
                    "n_train": 0,
                    "n_test": 0,
                    "passes_holdout": status == "skipped",
                    "held_out_sources": "",
                    "error": message,
                }
            )

    splits = pd.DataFrame(split_rows)
    overlaps = pd.DataFrame(overlap_rows)
    labels = pd.DataFrame(label_rows)
    splits.to_csv(out / "split_overlap_summary.csv", index=False)
    overlaps.to_csv(out / "field_overlap_by_split.csv", index=False)
    labels.to_csv(out / "label_balance_by_split.csv", index=False)
    feature_missingness(data, features).to_csv(out / "feature_missingness.csv", index=False)
    source_label_balance(data, label_col, source_col).to_csv(out / "source_label_balance.csv", index=False)
    grouped_label_balance(
        data,
        label_col,
        ["source_family", "upstream_source", "assay_type", "endpoint_type", "activity_type", "standard_type"],
    ).to_csv(out / "assay_context_label_balance.csv", index=False)
    label_evidence_composition(data, label_col).to_csv(out / "label_evidence_composition.csv", index=False)
    high_risk_column_scan(all_columns, label_col, features).to_csv(
        out / "high_risk_columns.csv", index=False
    )

    manifest: dict[str, Any] = {
        "dataset_path": str(dataset_path),
        "label_col": label_col,
        "n_rows": int(len(data)),
        "input_column_count": int(len(all_columns)),
        "audit_projection_columns": list(data.columns),
        "label_balance": label_balance(data, label_col),
        "duplicate_label_conflicts": duplicate_label_conflicts(data, label_col),
        "feature_set": feature_set,
        "features_checked": features,
        "leaky_features": leaky_features,
        "allow_partial_rescoring_features": allow_partial_rescoring_features,
        "allow_label_definition_features": allow_label_definition_features,
        "split_modes": modes,
        "split_results": split_rows,
        "all_requested_splits_passed": bool(
            len(splits) > 0
            and splits["status"].isin(["ok", "skipped"]).all()
            and splits["passes_holdout"].fillna(False).all()
        ),
        "outputs": {
            "split_summary": str(out / "split_overlap_summary.csv"),
            "field_overlap": str(out / "field_overlap_by_split.csv"),
            "label_balance": str(out / "label_balance_by_split.csv"),
            "feature_missingness": str(out / "feature_missingness.csv"),
            "source_label_balance": str(out / "source_label_balance.csv"),
            "assay_context_label_balance": str(out / "assay_context_label_balance.csv"),
            "label_evidence_composition": str(out / "label_evidence_composition.csv"),
            "high_risk_columns": str(out / "high_risk_columns.csv"),
            "split_manifests": str(out / "split_manifests"),
            "manifest": str(out / "leakage_overlap_audit.json"),
        },
    }
    (out / "leakage_overlap_audit.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
