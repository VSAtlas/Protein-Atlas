from __future__ import annotations

import gc
import json
from pathlib import Path
from typing import Any

from analysis.ml.audit_data_card import write_ml_data_card
from analysis.ml.audit_utils import load_table
from analysis.ml.binding_label_contract_audit import audit_binding_label_contract
from analysis.ml.claim_readiness import write_audit_claim_readiness
from analysis.ml.data_gap_report import run_data_gap_report
from analysis.ml.dataset_independence_audit import audit_dataset_independence
from analysis.ml.decoy_bias_audit import audit_decoy_bias
from analysis.ml.descriptor_constellation_audit import (
    run_descriptor_constellation_audit,
)
from analysis.ml.leakage_overlap_audit import DEFAULT_SPLITS, audit_ml_leakage_overlap
from analysis.ml.labels import training_eligibility_mask
from analysis.ml.score_scale_audit import write_score_scale_audit
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
    raw = load_table(Path(dataset_path))
    eligibility, eligibility_summary = training_eligibility_mask(raw)
    effective_dataset = Path(dataset_path)
    if eligibility_summary["excluded_rows"]:
        effective_dataset = out / "training_eligible_audit_input.csv"
        raw.loc[eligibility].to_csv(effective_dataset, index=False)
    (out / "training_eligibility_manifest.json").write_text(
        json.dumps(eligibility_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    # Wide model-ready tables contain hundreds of object columns. Keeping the
    # initial frame alive while each component reloads the table can double the
    # peak memory of the audit suite and trigger an OOM kill.
    del raw, eligibility
    gc.collect()
    data_card = write_ml_data_card(
        effective_dataset, label_col, feature_set, exclude_features, out
    )
    binding_label_contract = audit_binding_label_contract(
        effective_dataset, out / "binding_label_contract"
    )
    score_scale = write_score_scale_audit(
        effective_dataset,
        out / "score_scale",
        feature_names=data_card.get("features_available") or [],
    )
    try:
        ligand_identity = run_descriptor_constellation_audit(
            effective_dataset,
            out / "ligand_identity",
            label_col=label_col,
        )
    except ValueError as exc:
        ligand_identity = {"status": "skipped", "reason": str(exc)}
        identity_dir = out / "ligand_identity"
        identity_dir.mkdir(parents=True, exist_ok=True)
        (identity_dir / "descriptor_constellation_audit_manifest.json").write_text(
            json.dumps(ligand_identity, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    leakage = audit_ml_leakage_overlap(
        effective_dataset,
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
        "binding_label_contract": binding_label_contract.get("outputs"),
        "score_scale": score_scale["outputs"],
        "ligand_identity": str(
            out
            / "ligand_identity"
            / "descriptor_constellation_audit_manifest.json"
        ),
        "leakage": leakage["outputs"],
    }
    decoy_bias = audit_decoy_bias(dataset_path, out / "decoy_bias", seed=seed)
    outputs["decoy_bias"] = decoy_bias.get("outputs", str(out / "decoy_bias" / "decoy_bias_summary.json"))
    data_gap = run_data_gap_report(
        effective_dataset, out / "data_gap", labels=[label_col]
    )
    outputs["data_gap"] = data_gap.get("outputs", str(out / "data_gap" / "data_gap_manifest.json"))
    independence: dict[str, Any] | None = None
    if candidates:
        independence = audit_dataset_independence(
            effective_dataset,
            candidates,
            out / "independence",
            strict_pair_overlap=True,
            strict_source_overlap=strict_source_overlap,
        )
        outputs["independence"] = independence["outputs"]
    source_transfer: dict[str, Any] | None = None
    if source_transfer_train and source_transfer_test:
        source_transfer = audit_source_transfer(
            effective_dataset,
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
        score_scale=score_scale,
        ligand_identity=ligand_identity,
        binding_label_contract=binding_label_contract,
        claim_mode=claim_mode,
    )
    outputs["claim_readiness"] = str(out / "ml_audit_claim_readiness.json")
    manifest: dict[str, Any] = {
        "dataset_path": str(dataset_path),
        "effective_dataset_path": str(effective_dataset),
        "training_eligibility": eligibility_summary,
        "label_col": label_col,
        "feature_set": feature_set,
        "claim_mode": claim_mode,
        "allow_label_definition_features": allow_label_definition_features,
        "data_card": data_card,
        "binding_label_contract": binding_label_contract,
        "leakage_all_requested_splits_passed": leakage.get("all_requested_splits_passed", False),
        "independence_all_candidates_passed": independence.get("all_candidates_passed") if independence else None,
        "source_transfer_findings": source_transfer.get("findings") if source_transfer else [],
        "decoy_bias": decoy_bias,
        "score_scale": score_scale,
        "ligand_identity": ligand_identity,
        "data_gap": data_gap,
        "claim_readiness": readiness,
        "outputs": outputs,
    }
    outputs["training_eligibility"] = str(
        out / "training_eligibility_manifest.json"
    )
    (out / "ml_audit_suite_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
