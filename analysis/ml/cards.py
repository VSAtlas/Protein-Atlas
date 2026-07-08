from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_model_card(
    *,
    out_path: str | Path,
    model_type: str,
    selected_model: str,
    label_col: str,
    feature_set: str,
    features: list[str],
    split_mode: str,
    dataset_manifest: dict[str, Any],
    metrics: dict[str, Any],
    claim_readiness: dict[str, Any],
    pu_mode: str,
    hpo_manifest: dict[str, Any],
    subgroup_metrics_path: str | Path,
) -> dict[str, Any]:
    out = Path(out_path)
    card: dict[str, Any] = {
        "model_type": model_type,
        "selected_model": selected_model,
        "label_col": label_col,
        "feature_set": feature_set,
        "features": list(features),
        "retained_excluded_features": list(dataset_manifest.get("retained_excluded_features") or []),
        "retained_excluded_policy": dataset_manifest.get("retained_excluded_policy"),
        "retained_context_not_trained": list(dataset_manifest.get("retained_context_not_trained") or []),
        "retained_context_not_trained_policy": dataset_manifest.get("retained_context_not_trained_policy"),
        "split_mode": split_mode,
        "dataset_version_id": dataset_manifest.get("dataset_version_id"),
        "dataset_path": dataset_manifest.get("dataset_path"),
        "metrics": metrics,
        "claim_readiness": claim_readiness,
        "pu_mode": pu_mode,
        "hpo": hpo_manifest,
        "subgroup_metrics": str(subgroup_metrics_path),
        "limitations": [
            "Report AUROC together with AUPRC, prevalence, calibration, enrichment, and applicability domain.",
            "Use random splits only as smoke tests; drug, target, scaffold, source, temporal, or external holdouts are required for stronger claims.",
            "Unknown labels, PU background, weak negatives, and benchmark decoys must remain distinct in interpretation.",
        ],
    }
    (out / "model_card.json").write_text(json.dumps(card, indent=2, sort_keys=True, default=str), encoding="utf-8")
    lines = [
        f"# Atlas ML Model Card: {selected_model}",
        "",
        f"- Label: `{label_col}`",
        f"- Feature set: `{feature_set}`",
        f"- Split: `{split_mode}`",
        f"- Retained but explicitly excluded columns: `{len(card.get('retained_excluded_features') or [])}`",
        f"- Retained context columns not trained: `{len(card.get('retained_context_not_trained') or [])}`",
        f"- Dataset version: `{card.get('dataset_version_id')}`",
        f"- Claim status: `{claim_readiness.get('overall_status')}`",
        "",
        "## Required Caveats",
        "",
        *[f"- {item}" for item in card["limitations"]],
    ]
    (out / "model_card.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return card
