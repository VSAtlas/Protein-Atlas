from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


def _add_issue(report: dict[str, Any], severity: str, message: str, *, dataset: str | None = None) -> None:
    entry = {"severity": severity, "message": message}
    if dataset:
        entry["dataset"] = dataset
    report.setdefault("issues", []).append(entry)
    report.setdefault("blockers" if severity == "blocker" else "warnings", []).append(message if not dataset else f"{dataset}: {message}")


def _metadata_issue_severity(claim_mode: str) -> str:
    return "blocker" if claim_mode == "publication" else "warning"


def dataset_preflight_checks(
    dataset: str | Path,
    *,
    label_col: str | None,
    feature_set: str | None,
    claim_mode: str = "exploratory",
) -> dict[str, Any]:
    path = Path(dataset)
    report: dict[str, Any] = {
        "dataset": str(path),
        "label_col": label_col,
        "feature_set": feature_set,
        "status": "ok",
        "checks": {},
        "warnings": [],
        "blockers": [],
        "issues": [],
    }
    if not path.exists():
        _add_issue(report, "blocker", f"dataset does not exist: {path}")
        report["status"] = "blocked"
        return report
    try:
        frame = pd.read_csv(path, low_memory=False)
    except Exception as exc:
        _add_issue(report, "blocker", f"could not read dataset: {exc}")
        report["status"] = "blocked"
        return report
    report["checks"]["rows"] = int(len(frame))
    report["checks"]["columns"] = int(len(frame.columns))
    if frame.empty:
        _add_issue(report, "blocker", "dataset has zero rows")
    if label_col:
        if label_col not in frame.columns:
            _add_issue(report, "blocker", f"label column missing: {label_col}")
        else:
            labels = pd.to_numeric(frame[label_col], errors="coerce")
            n_pos = int((labels == 1).sum())
            n_neg = int((labels == 0).sum())
            n_unlabeled = int(labels.isna().sum())
            invalid = int((labels.notna() & ~labels.isin([0, 1])).sum())
            report["checks"].update({"n_positive": n_pos, "n_negative": n_neg, "n_unlabeled": n_unlabeled, "n_invalid_labels": invalid})
            if invalid:
                _add_issue(report, "blocker", f"non-binary labels present: {invalid}")
            if n_pos < 20:
                _add_issue(report, "warning", f"few positives for claim-grade training: {n_pos} < 20")
            if n_neg < 20:
                _add_issue(report, "warning", f"few negatives for claim-grade training: {n_neg} < 20")
            if n_unlabeled and not any(col in frame.columns for col in ["negative_evidence_type", "negative_confidence", "_sample_weight"]):
                _add_issue(report, "warning", "unlabeled rows are present without PU/negative provenance columns")
    provenance_cols = [col for col in ["label_source", "source_family", "upstream_source", "assay_type", "endpoint_type"] if col in frame.columns]
    report["checks"]["source_columns_present"] = provenance_cols
    if not provenance_cols:
        _add_issue(report, _metadata_issue_severity(claim_mode), "label provenance columns are missing")
    scaffold_cols = [col for col in ["scaffold_key", "chemical_cluster", "ligand_chemotype", "canonical_smiles", "smiles"] if col in frame.columns]
    report["checks"]["chemical_holdout_columns_present"] = scaffold_cols
    if not scaffold_cols:
        _add_issue(report, _metadata_issue_severity(claim_mode), "scaffold/chemical-cluster metadata is missing")
    family_cols = [col for col in ["target_family", "protein_class", "protein_family", "target_class"] if col in frame.columns]
    report["checks"]["target_family_columns_present"] = family_cols
    if not family_cols:
        _add_issue(report, _metadata_issue_severity(claim_mode), "target-family/protein-class metadata is missing")
    year_cols = [col for col in ["activity_publication_year", "label_publication_year", "database_release_year", "evidence_publication_year"] if col in frame.columns]
    report["checks"]["temporal_columns_present"] = year_cols
    if not year_cols:
        _add_issue(report, _metadata_issue_severity(claim_mode), "temporal year metadata is missing")
    if "label_source" in frame.columns and not frame.empty:
        top_fraction = float(frame["label_source"].fillna("missing").value_counts(normalize=True).max())
        report["checks"]["largest_label_source_fraction"] = top_fraction
        if top_fraction > 0.90:
            _add_issue(report, "warning", f"label source imbalance is high: largest source fraction={top_fraction:.3f}")
    if label_col and feature_set:
        try:
            from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
            from analysis.ml.leakage_checks import find_leaky_features
            from analysis.ml.score_scale_audit import audit_score_scale_frame

            excluded = effective_exclude_features(label_col, None)
            features = [feature for feature in get_feature_set(feature_set) if feature in frame.columns and feature not in excluded]
            leaky = find_leaky_features(features)
            report["checks"]["features_available"] = len(features)
            if not features:
                _add_issue(report, "blocker", f"feature set has no available non-excluded features: {feature_set}")
            if "final_score" in features:
                final_score = pd.to_numeric(frame["final_score"], errors="coerce")
                eligible = pd.Series(True, index=frame.index)
                if label_col and label_col in frame.columns:
                    eligible = pd.to_numeric(frame[label_col], errors="coerce").isin([0, 1])
                denominator = int(eligible.sum())
                covered = int((eligible & final_score.notna()).sum())
                coverage = float(covered / denominator) if denominator else 0.0
                report["checks"]["final_score_coverage"] = {
                    "eligible_rows": denominator,
                    "covered_rows": covered,
                    "fraction": coverage,
                    "source_column_present": "final_score_source" in frame.columns,
                }
                if "final_score_source" in frame.columns:
                    report["checks"]["final_score_source_counts"] = (
                        frame.loc[eligible, "final_score_source"]
                        .fillna("missing")
                        .astype(str)
                        .value_counts(dropna=False)
                        .to_dict()
                    )
                if "ml_blend_mode" in frame.columns:
                    mode_counts = (
                        frame.loc[eligible & final_score.notna(), "ml_blend_mode"]
                        .fillna("missing")
                        .astype(str)
                        .loc[lambda values: values.str.strip().ne("")]
                        .value_counts()
                        .to_dict()
                    )
                    report["checks"]["final_score_blend_mode_counts"] = mode_counts
                    if len(mode_counts) > 1:
                        severity = "blocker" if claim_mode == "publication" else "warning"
                        _add_issue(report, severity, "final_score mixes SCORCH/CNN availability modes; use mode-stratified DUD nulls")
                if coverage < 0.80:
                    severity = "blocker" if claim_mode == "publication" else "warning"
                    _add_issue(report, severity, f"final_score coverage is incomplete: {coverage:.3f} < 0.800")
            if leaky:
                _add_issue(report, "blocker", f"leaky predictive features selected: {', '.join(leaky)}")
            score_scale, _, _ = audit_score_scale_frame(frame, feature_names=features)
            report["checks"]["score_scale"] = score_scale
            if score_scale["status"] == "failed":
                _add_issue(
                    report,
                    "blocker",
                    "selected score features mix raw consensus percentiles with z-score semantics",
                )
            elif score_scale["status"] == "warning":
                _add_issue(
                    report,
                    "warning",
                    "score provenance is incomplete; inspect score-scale audit before interpretation",
                )
        except Exception as exc:
            _add_issue(report, "warning", f"could not run feature leakage check: {exc}")
    if report["blockers"]:
        report["status"] = "blocked"
    elif report["warnings"]:
        report["status"] = "ready_with_warnings"
    return report


def discover_default_expert_datasets(run_dir: str | Path) -> list[dict[str, Any]]:
    from analysis.cli.run_ml_training_pass import EXPERT_DEFAULTS

    root = Path(run_dir)
    rows: list[dict[str, Any]] = []
    for expert, defaults in EXPERT_DEFAULTS.items():
        candidates = [root / rel for rel in defaults.get("candidates", [])]
        found = next((path for path in candidates if path.exists()), None)
        rows.append(
            {
                "expert": expert,
                "dataset": str(found) if found else None,
                "label": defaults.get("label"),
                "feature_set": defaults.get("feature_set"),
                "split": defaults.get("split"),
                "searched": [str(path) for path in candidates],
            }
        )
    return rows


def build_ml_doctor_report(
    *,
    run_id: str,
    run_dir: str | Path,
    dataset: str | Path | None = None,
    label_col: str | None = None,
    feature_set: str | None = None,
    claim_mode: str = "exploratory",
    strict: bool = False,
) -> dict[str, Any]:
    root = Path(run_dir)
    report: dict[str, Any] = {
        "run_id": str(run_id),
        "run_dir": str(root),
        "claim_mode": claim_mode,
        "status": "ok",
        "blockers": [],
        "warnings": [],
        "checks": {},
        "datasets": [],
    }
    required_modules = [
        "analysis.cli.run_ml_training_pass",
        "analysis.cli.run_ml_audit_suite",
        "analysis.cli.run_ml_source_pu_suite",
        "analysis.cli.refresh_ml_feature_metadata",
    ]
    for module_name in required_modules:
        try:
            importlib.import_module(module_name)
            report["checks"][module_name] = "ok"
        except Exception as exc:
            report["checks"][module_name] = "failed"
            report["blockers"].append(f"cannot import {module_name}: {exc}")
    if not root.exists():
        report["warnings"].append(f"run directory does not exist yet: {root}")
    if dataset:
        dataset_report = dataset_preflight_checks(dataset, label_col=label_col, feature_set=feature_set, claim_mode=claim_mode)
        report["datasets"].append(dataset_report)
    else:
        found_default_dataset = False
        for item in discover_default_expert_datasets(root):
            if not item.get("dataset"):
                report["warnings"].append(f"default dataset for expert {item['expert']!r} was not found under {root}")
                report["datasets"].append({**item, "status": "missing"})
                continue
            found_default_dataset = True
            dataset_report = dataset_preflight_checks(
                item["dataset"],
                label_col=str(item.get("label") or ""),
                feature_set=str(item.get("feature_set") or ""),
                claim_mode=claim_mode,
            )
            dataset_report.update({"expert": item["expert"], "split": item.get("split")})
            report["datasets"].append(dataset_report)
        if not found_default_dataset:
            report["blockers"].append(f"no default ML datasets were found under {root}")
    for dataset_report in report["datasets"]:
        for blocker in dataset_report.get("blockers", []):
            report["blockers"].append(f"{dataset_report.get('expert') or dataset_report.get('dataset')}: {blocker}")
        for warning in dataset_report.get("warnings", []):
            report["warnings"].append(f"{dataset_report.get('expert') or dataset_report.get('dataset')}: {warning}")
    if strict and report["warnings"]:
        report["blockers"].extend(report["warnings"])
        report["warnings"] = []
    if report["blockers"]:
        report["status"] = "blocked"
    elif report["warnings"]:
        report["status"] = "ready_with_warnings"
    return report


def write_doctor_report(report: dict[str, Any], out_path: str | Path) -> Path:
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return path
