from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _metric_value(metrics: dict[str, Any], name: str) -> float | None:
    value = metrics.get(name)
    try:
        return float(value) if value is not None and pd.notna(value) else None
    except (TypeError, ValueError):
        return None


def _status(blockers: list[str], warnings: list[str]) -> str:
    if blockers:
        return "blocked"
    if warnings:
        return "ready_with_caveats"
    return "ready"


def write_model_claim_readiness(
    *,
    out_path: str | Path,
    pred: pd.DataFrame,
    label_col: str,
    metrics: dict[str, Any],
    split_summary: dict[str, Any],
    features: list[str],
    validation_rows: int,
    conformal_status: str,
    baseline_path: str | Path,
    claim_mode: str = "exploratory",
    dataset_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    out = Path(out_path)
    mode = str(claim_mode or "exploratory").strip().lower()
    publication = mode == "publication"
    blockers: list[str] = []
    warnings: list[str] = []
    labels = pd.to_numeric(pred[label_col], errors="coerce") if label_col in pred.columns else pd.Series(dtype=float)
    n_pos = int((labels == 1).sum()) if len(labels) else 0
    n_neg = int((labels == 0).sum()) if len(labels) else 0
    if n_pos < 20:
        blockers.append(f"test positives below publication-grade minimum: {n_pos} < 20")
    if n_neg < 20:
        blockers.append(f"test negatives below publication-grade minimum: {n_neg} < 20")
    if not bool(split_summary.get("passes_holdout", False)):
        blockers.append("requested outer split did not pass holdout overlap checks")
    split_mode = str(split_summary.get("split_mode") or "")
    if publication and split_mode in {"", "random"}:
        blockers.append("publication claim requires a non-random outer holdout split")
    if conformal_status != "ok":
        target = blockers if publication else warnings
        target.append("conformal prediction-set coverage unavailable for this model run")
    if validation_rows <= 0:
        target = blockers if publication else warnings
        target.append("no validation/calibration split was held out for model selection or uncertainty calibration")
    if not Path(baseline_path).exists():
        target = blockers if publication else warnings
        target.append("standard baseline panel missing")

    ece = _metric_value(metrics, "ECE")
    brier = _metric_value(metrics, "Brier")
    if ece is not None and ece > 0.10:
        target = blockers if publication else warnings
        target.append(f"calibration ECE is high: {ece:.3f}")
    if brier is not None and brier > 0.25:
        target = blockers if publication else warnings
        target.append(f"Brier score is high: {brier:.3f}")

    if "source_family" not in pred.columns and "label_source" not in pred.columns:
        target = blockers if publication else warnings
        target.append("source lineage columns absent from prediction output")
    if "chemical_cluster" not in pred.columns and "scaffold_key" not in pred.columns:
        target = blockers if publication else warnings
        target.append("chemical cluster/scaffold metadata absent from prediction output")
    if "target_family" not in pred.columns and "protein_class" not in pred.columns:
        target = blockers if publication else warnings
        target.append("target family/protein class metadata absent from prediction output")
    if any(col in features for col in {"drug_id", "target_id", "pdb_id", "literature_supported_label"}):
        blockers.append("raw identifiers or evaluation labels are present in predictive features")
    sampled = bool((dataset_provenance or {}).get("sampled"))
    if sampled:
        target = blockers if publication else warnings
        target.append("model was trained on a deterministic sampled/smoke dataset, not the full model-ready table")

    manifest = {
        "claim_mode": mode,
        "overall_status": _status(blockers, warnings),
        "blockers": blockers,
        "warnings": warnings,
        "n_test": int(len(pred)),
        "n_test_positive": n_pos,
        "n_test_negative": n_neg,
        "positive_rate": float((labels == 1).mean()) if len(labels.dropna()) else None,
        "sampled": sampled,
        "dataset_provenance": dataset_provenance or {},
        "split_summary": split_summary,
        "validation_rows": int(validation_rows),
        "conformal_status": conformal_status,
        "primary_metrics": {
            "AUROC": _metric_value(metrics, "AUROC"),
            "AUPRC": _metric_value(metrics, "AUPRC"),
            "Brier": brier,
            "ECE": ece,
            "EF@1%": _metric_value(metrics, "EF@1%"),
            "EF@5%": _metric_value(metrics, "EF@5%"),
            "EF@10%": _metric_value(metrics, "EF@10%"),
        },
        "required_claim_language": [
            "Report AUROC with AUPRC, prevalence, calibration, early enrichment, and applicability domain.",
            "Use random split only as a smoke test; primary claims should use drug/target/cluster/source/temporal or external validation as appropriate.",
            "Unknown labels remain unknown; weak negatives, PU background, and benchmark decoys must be reported separately.",
        ],
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "model_claim_readiness.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    rows = [
        {"severity": "blocker", "finding": finding}
        for finding in blockers
    ] + [
        {"severity": "warning", "finding": finding}
        for finding in warnings
    ]
    if not rows:
        rows = [{"severity": "ok", "finding": "no claim-readiness blockers or warnings detected"}]
    pd.DataFrame(rows).to_csv(out / "model_claim_readiness.csv", index=False)
    return manifest


def write_audit_claim_readiness(
    *,
    out_path: str | Path,
    data_card: dict[str, Any],
    leakage_manifest: dict[str, Any],
    source_transfer_findings: list[dict[str, Any]] | list[Any],
    independence_passed: bool | None,
    decoy_bias: dict[str, Any] | None = None,
    score_scale: dict[str, Any] | None = None,
    ligand_identity: dict[str, Any] | None = None,
    binding_label_contract: dict[str, Any] | None = None,
    claim_mode: str = "exploratory",
) -> dict[str, Any]:
    out = Path(out_path)
    mode = str(claim_mode or "exploratory").strip().lower()
    publication = mode == "publication"
    blockers: list[str] = []
    warnings: list[str] = []
    if leakage_manifest.get("leaky_features"):
        blockers.append(f"leaky features detected: {', '.join(leakage_manifest.get('leaky_features', []))}")
    if not leakage_manifest.get("all_requested_splits_passed", False):
        blockers.append("one or more requested leakage/holdout splits failed")
    if independence_passed is False:
        blockers.append("candidate benchmark independence audit failed")
    if data_card.get("n_positive", 0) < 50:
        warnings.append("few observed positives; report confidence intervals and avoid broad claims")
    if data_card.get("n_negative", 0) < 50:
        warnings.append("few observed negatives; prefer PU/matched-background sensitivity framing")
    if data_card.get("n_chemical_clusters") in {None, 0, 1}:
        target = blockers if publication else warnings
        target.append("chemical cluster metadata is missing or not diverse enough for OOD claims")
    if data_card.get("n_target_families") in {None, 0, 1}:
        target = blockers if publication else warnings
        target.append("target-family metadata is missing or not diverse enough for target-family OOD claims")
    if not data_card.get("has_temporal_metadata"):
        target = blockers if publication else warnings
        target.append("temporal metadata absent; prospective/temporal claims are blocked for this table")
    split_results = leakage_manifest.get("split_results") or []
    skipped_required = [
        row.get("split_mode")
        for row in split_results
        if row.get("status") != "ok" and row.get("split_mode") != "random"
    ]
    if publication and skipped_required:
        blockers.append(f"publication claim requires non-skipped OOD/temporal/source splits: {sorted(set(skipped_required))}")
    if source_transfer_findings:
        target = blockers if publication else warnings
        target.append("source-transfer audit emitted findings; inspect source_transfer outputs before external-validation claims")
    if decoy_bias and decoy_bias.get("flags"):
        warnings.append("decoy-bias audit emitted flags; decoy benchmark claims need caveats or stronger controls")
    if score_scale and score_scale.get("status") == "failed":
        blockers.append(
            "selected score features mix raw consensus percentiles with decoy-standardized z-scores"
        )
    elif score_scale and score_scale.get("status") == "warning":
        warnings.append("score-scale audit emitted warnings; inspect score_scale outputs")
    if ligand_identity and ligand_identity.get("status") == "failed":
        target = blockers if publication else warnings
        target.append(
            "ligand identity audit found "
            f"{ligand_identity.get('n_drugs_with_multiple_practical_signatures', 0)} "
            "drug IDs mapped to multiple RDKit descriptor constellations"
        )
    if binding_label_contract and binding_label_contract.get("status") == "failed":
        blockers.append(
            "SPD binding labels contradict the censor-aware AC50 relation policy or lack "
            "the required binding-policy version"
        )
    manifest = {
        "claim_mode": mode,
        "overall_status": _status(blockers, warnings),
        "blockers": blockers,
        "warnings": warnings,
        "data_card_summary": {
            "n_rows": data_card.get("n_rows"),
            "n_positive": data_card.get("n_positive"),
            "n_negative": data_card.get("n_negative"),
            "positive_rate": data_card.get("positive_rate"),
            "n_chemical_clusters": data_card.get("n_chemical_clusters"),
            "n_target_families": data_card.get("n_target_families"),
            "n_source_families": data_card.get("n_source_families"),
            "temporal_metadata_columns": data_card.get("temporal_metadata_columns"),
        },
        "decoy_bias_status": decoy_bias.get("status") if decoy_bias else None,
        "decoy_bias_flags": decoy_bias.get("flags") if decoy_bias else [],
        "score_scale_status": score_scale.get("status") if score_scale else None,
        "score_scale_flags": score_scale.get("flags") if score_scale else [],
        "ligand_identity_status": (
            ligand_identity.get("status") if ligand_identity else None
        ),
        "ligand_identity_conflicted_drugs": (
            ligand_identity.get("n_drugs_with_multiple_practical_signatures")
            if ligand_identity
            else None
        ),
        "binding_label_contract_status": (
            binding_label_contract.get("status") if binding_label_contract else None
        ),
        "binding_label_contradictions": (
            binding_label_contract.get("n_label_contradictions")
            if binding_label_contract
            else None
        ),
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "ml_audit_claim_readiness.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    rows = [{"severity": "blocker", "finding": value} for value in blockers]
    rows += [{"severity": "warning", "finding": value} for value in warnings]
    if not rows:
        rows = [{"severity": "ok", "finding": "audit passed claim-readiness checks"}]
    pd.DataFrame(rows).to_csv(out / "ml_audit_claim_readiness.csv", index=False)
    return manifest
