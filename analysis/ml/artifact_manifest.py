from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


KNOWN_MODEL_ARTIFACTS = [
    "dataset_version_manifest.json",
    "model_predictions.csv",
    "model_metrics.csv",
    "model_metric_bootstrap_ci.csv",
    "model_reliability_table.csv",
    "model_grouped_calibration.csv",
    "model_subgroup_metrics.csv",
    "model_decision_metrics.csv",
    "model_group_topk_recovery.csv",
    "model_abstention_summary.csv",
    "standard_baseline_panel.csv",
    "feature_importance.csv",
    "split_manifest.csv",
    "split_manifest.json",
    "split_summary.json",
    "hpo_trials.csv",
    "hpo_manifest.json",
    "conformal_summary.json",
    "conformal_predictions.csv",
    "model_applicability_domain_summary.csv",
    "chemical_fingerprint_applicability_domain_summary.csv",
    "target_family_applicability_domain_summary.csv",
    "pu_training_manifest.json",
    "bagging_pu_members.csv",
    "pu_propensity_weights.csv",
    "pu_elkan_noto_manifest.json",
    "model_claim_readiness.json",
    "model_claim_readiness.csv",
    "model_card.json",
    "model_card.md",
    "model_preprocessing.json",
    "training_design_matrix.csv",
    "trained_model.pkl",
]


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_known_model_artifacts(model_dir: str | Path) -> list[dict[str, Any]]:
    root = Path(model_dir)
    artifacts: list[dict[str, Any]] = []
    for name in KNOWN_MODEL_ARTIFACTS:
        path = root / name
        if not path.exists() or not path.is_file():
            continue
        artifacts.append(
            {
                "path": str(path),
                "name": name,
                "bytes": int(path.stat().st_size),
                "sha256": file_sha256(path),
            }
        )
    return artifacts


def write_artifact_manifest(model_dir: str | Path, artifacts: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    root = Path(model_dir)
    manifest = {
        "model_dir": str(root),
        "artifacts": artifacts if artifacts is not None else collect_known_model_artifacts(root),
    }
    (root / "model_artifacts.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
