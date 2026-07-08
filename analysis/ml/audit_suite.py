from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from analysis.ml.audit_data_card import write_ml_data_card
from analysis.ml.claim_readiness import write_audit_claim_readiness
from analysis.ml.data_gap_report import run_data_gap_report
from analysis.ml.dataset_independence_audit import audit_dataset_independence
from analysis.ml.decoy_bias_audit import audit_decoy_bias
from analysis.ml.leakage_overlap_audit import DEFAULT_SPLITS, audit_ml_leakage_overlap
from analysis.ml.source_transfer_audit import audit_source_transfer


def run_ml_audit_suite(
    dataset_path: str | Path,
    label_col: str,
    out_dir: str | Path,
    *,
    feature_set: str | None = None,
    exclude_features: list[str] | None = None,
    splits: list[str] | None = None,
    source_col: str = "label_source",
    seed: int = 42,
    test_fraction: float = 0.2,
    allow_partial_rescoring_features: bool = False,
    allow_label_definition_features: bool = False,
    candidates: list[str] | None = None,
    strict_source_overlap: bool = False,
    source_transfer_train: list[str] | None = None,
    source_transfer_test: list[str] | None = None,
    predictions_path: str | Path | None = None,
    claim_mode: str = "exploratory",
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    data_card = write_ml_data_card(dataset_path, label_col, feature_set, exclude_features, out)
    leakage = audit_ml_leakage_overlap(
        dataset_path,
        label_col,
        out / "leakage",
        feature_set=feature_set,
        exclude_features=exclude_features,
        split_modes=splits or DEFAULT_SPLITS,
        source_col=source_col,
        seed=seed,
        test_fraction=test_fraction,
        allow_partial_rescoring_features=allow_partial_rescoring_features,
        allow_label_definition_features=allow_label_definition_features,
    )
    outputs: dict[str, Any] = {
        "data_card": str(out / "ml_data_card.json"),
        "leakage": leakage["outputs"],
    }
    decoy_bias = audit_decoy_bias(dataset_path, out / "decoy_bias", seed=seed)
    outputs["decoy_bias"] = decoy_bias.get("outputs", str(out / "decoy_bias" / "decoy_bias_summary.json"))
    data_gap = run_data_gap_report(dataset_path, out / "data_gap", labels=[label_col])
    outputs["data_gap"] = data_gap.get("outputs", str(out / "data_gap" / "data_gap_manifest.json"))
    independence: dict[str, Any] | None = None
    if candidates:
        independence = audit_dataset_independence(
            dataset_path,
            candidates,
            out / "independence",
            strict_pair_overlap=True,
            strict_source_overlap=strict_source_overlap,
        )
        outputs["independence"] = independence["outputs"]
    source_transfer: dict[str, Any] | None = None
    if source_transfer_train and source_transfer_test:
        source_transfer = audit_source_transfer(
            dataset_path,
            label_col,
            out / "source_transfer",
            train_sources=source_transfer_train,
            test_sources=source_transfer_test,
            predictions_path=predictions_path,
            feature_cols=data_card.get("features_available") or None,
        )
        outputs["source_transfer"] = str(out / "source_transfer" / "source_transfer_audit_manifest.json")
    readiness = write_audit_claim_readiness(
        out_path=out,
        data_card=data_card,
        leakage_manifest=leakage,
        source_transfer_findings=(source_transfer.get("findings") or []) if source_transfer else [],
        independence_passed=independence.get("all_candidates_passed") if independence else None,
        decoy_bias=decoy_bias,
        claim_mode=claim_mode,
    )
    outputs["claim_readiness"] = str(out / "ml_audit_claim_readiness.json")
    manifest: dict[str, Any] = {
        "dataset_path": str(dataset_path),
        "label_col": label_col,
        "feature_set": feature_set,
        "claim_mode": claim_mode,
        "allow_label_definition_features": allow_label_definition_features,
        "data_card": data_card,
        "leakage_all_requested_splits_passed": leakage.get("all_requested_splits_passed", False),
        "independence_all_candidates_passed": independence.get("all_candidates_passed") if independence else None,
        "source_transfer_findings": source_transfer.get("findings") if source_transfer else [],
        "decoy_bias": decoy_bias,
        "data_gap": data_gap,
        "claim_readiness": readiness,
        "outputs": outputs,
    }
    (out / "ml_audit_suite_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
