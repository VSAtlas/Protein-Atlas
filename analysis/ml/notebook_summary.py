from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"status": "unreadable", "error": str(exc), "path": str(path)}


def _manifest_rows(root: Path) -> list[dict[str, Any]]:
    patterns = [
        "training_pass/ml_training_pass_manifest.json",
        "audit_suite/ml_audit_suite_manifest.json",
        "source_pu_suite/ml_source_pu_suite_manifest.json",
        "optuna_sweep/ml_optuna_sweep_manifest.json",
        "optuna_sweep/hpo/optuna/optuna_sweep_manifest.json",
        "*/ml_optuna_sweep_manifest.json",
        "*/hpo/optuna/optuna_sweep_manifest.json",
        "chemprop/chemprop_training_manifest.json",
        "*/chemprop_training_manifest.json",
        "reinvent_generation/reinvent_generation_manifest.json",
        "*/reinvent_generation_manifest.json",
        "*_benchmark/*_benchmark_manifest.json",
        "*/*_benchmark_manifest.json",
        "tdc/*.manifest.json",
    ]
    rows: list[dict[str, Any]] = []
    seen_paths: set[Path] = set()
    for pattern in patterns:
        for path in sorted(root.glob(pattern)):
            if path in seen_paths:
                continue
            seen_paths.add(path)
            payload = _read_json(path)
            rows.append(
                {
                    "kind": path.stem,
                    "path": str(path),
                    "status": payload.get("status") or payload.get("overall_status"),
                    "adapter": payload.get("adapter"),
                    "claim_status": (payload.get("claim_readiness") or {}).get("overall_status"),
                    "rows": payload.get("rows") or payload.get("staged_rows") or payload.get("n_trials"),
                    "dataset": payload.get("dataset") or payload.get("dataset_path"),
                    "manifest": payload,
                }
            )
    return rows


def load_ml_run_summary(run_id_or_ml_root: str | Path, *, repo_root: str | Path = ".") -> pd.DataFrame:
    """Load run-scoped Atlas ML manifests into a notebook-friendly table."""
    candidate = Path(run_id_or_ml_root)
    if candidate.exists():
        root = candidate
    else:
        root = Path(repo_root) / "outputs" / "data" / str(run_id_or_ml_root) / "ml"
    rows = _manifest_rows(root)
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=["kind", "path", "status", "adapter", "claim_status", "rows", "dataset"])
    return frame


def load_model_registry(model_root: str | Path) -> pd.DataFrame:
    """Load an Atlas model registry index or experiment JSONL into a DataFrame."""
    root = Path(model_root)
    index = root / "model_registry_index.csv"
    if index.exists():
        return pd.read_csv(index)
    jsonl = root / "experiment_runs.jsonl"
    if jsonl.exists():
        rows = [json.loads(line) for line in jsonl.read_text(encoding="utf-8").splitlines() if line.strip()]
        return pd.DataFrame(rows)
    return pd.DataFrame()


def display_ml_run_summary(run_id_or_ml_root: str | Path, *, repo_root: str | Path = ".") -> pd.DataFrame:
    """Display and return the Atlas ML manifest summary when IPython is available."""
    frame = load_ml_run_summary(run_id_or_ml_root, repo_root=repo_root)
    try:
        from IPython.display import display
    except Exception:
        return frame
    display(frame.drop(columns=["manifest"], errors="ignore"))
    return frame


def load_ml_leaderboard(run_id_or_ml_root: str | Path, *, repo_root: str | Path = ".") -> pd.DataFrame:
    """Load the run-scoped Atlas ML leaderboard, building from manifests if needed."""
    candidate = Path(run_id_or_ml_root)
    if candidate.exists():
        root = candidate
    else:
        root = Path(repo_root) / "outputs" / "data" / str(run_id_or_ml_root) / "ml"
    csv_path = root / "leaderboard" / "model_leaderboard.csv"
    if csv_path.exists():
        return pd.read_csv(csv_path)
    try:
        from analysis.ml.leaderboard import collect_ml_leaderboard

        return collect_ml_leaderboard(root)
    except Exception:
        return pd.DataFrame(columns=["kind", "expert", "model_dir", "AUROC", "AUPRC", "claim_status", "artifact_path"])


def display_ml_leaderboard(run_id_or_ml_root: str | Path, *, repo_root: str | Path = ".") -> pd.DataFrame:
    """Display and return the Atlas ML leaderboard when IPython is available."""
    frame = load_ml_leaderboard(run_id_or_ml_root, repo_root=repo_root)
    try:
        from IPython.display import display
    except Exception:
        return frame
    preferred = [
        "kind", "expert", "selected_model", "feature_set", "split_mode", "sampled",
        "claim_status", "AUROC", "AUPRC", "Brier", "ECE", "mlflow_run_id", "model_dir",
    ]
    display(frame[[col for col in preferred if col in frame.columns]] if not frame.empty else frame)
    return frame

def expert_feature_architecture() -> pd.DataFrame:
    """Return the current Atlas expert feature architecture for notebook review."""

    rows = [
        {
            "expert": "Clean Binding expert",
            "predicts": "Can this drug-target pair show measured activity?",
            "label": "spd_binding_label",
            "feature_set": "spd_binding_nonleaky",
            "input_policy": "Raw binding/structure/class/chemistry features; aggregate binding priors excluded.",
        },
        {
            "expert": "Frozen Binding prior / baseline",
            "predicts": "Rule-based structural binding prior from Atlas/BANANA/SCORCH-style scores.",
            "label": "None",
            "feature_set": "frozen_binding_prior_baseline",
            "input_policy": "Outputs banana_atlas_blend_score and binding_expert_score; use as a baseline or compact downstream prior, not as clean Binding expert input.",
        },
        {
            "expert": "Direct SPD Exposure expert",
            "predicts": "Is the interaction exposure-relevant under the SPD definition?",
            "label": "spd_exposure_label",
            "feature_set": "spd_exposure_nonleaky",
            "input_policy": "Uses frozen binding prior plus structure/class/chemistry; excludes AC50, Cmax, free Cmax, and exposure margin.",
        },
        {
            "expert": "Estimated-PK Exposure model",
            "predicts": "Estimates potency and free exposure separately, then derives exposure relevance.",
            "label": "spd_ac50_uM; free_cmax_um; final spd_exposure_label",
            "feature_set": "explicit opt-in via --run-estimated-pk-exposure",
            "input_policy": "Deferred v1 path; predictions must be out-of-fold before downstream use.",
        },
        {
            "expert": "Tissue expression / site scorer",
            "predicts": "Is the target expressed or biologically present in the ADR-relevant tissue/site?",
            "label": "None by default; optional independent tissue_site_label",
            "feature_set": "tissue_site_score_only",
            "input_policy": "Expression/tissue features only; site_relevance_score is the deterministic composite output.",
        },
        {
            "expert": "Pair-level Tissue/Site relevance expert",
            "predicts": "Is this drug-target pair relevant to the ADR tissue/site after combining binding, exposure, and expression?",
            "label": "independent tissue_site_label if available",
            "feature_set": "spd_tissue_site_nonleaky",
            "input_policy": "Supports binding_expert_score or out-of-fold binding_probability_clean plus exposure and tissue relevance features.",
        },
        {
            "expert": "Clean Mechanism expert",
            "predicts": "Is there plausible ADR mechanism support without direct graph/safety evidence?",
            "label": "mechanism_ml_label_clean preferred; mechanism_ml_label with caveats",
            "feature_set": "spd_mechanism_nonleaky",
            "input_policy": "Uses binding/exposure/tissue priors plus structure/class/chemistry; excludes graph/pathway/target-safety evidence.",
        },
        {
            "expert": "Mechanism evidence / sensitivity model",
            "predicts": "How much known graph/pathway/drug-ADR/target-ADR/safety evidence supports the mechanism?",
            "label": "mechanism_ml_label",
            "feature_set": "spd_mechanism_evidence_sensitivity",
            "input_policy": "Evidence-augmented sensitivity model; requires explicit --allow-label-definition-features.",
        },
        {
            "expert": "Graph/KGE mechanism model candidate",
            "predicts": "Can sparse drug-target-ADR/pathway graph structure predict missing mechanism-support edges?",
            "label": "held-out graph edges or panel-specific mechanism labels",
            "feature_set": "mechanism graph nodes/edges, not a tabular feature set",
            "input_policy": "Optional PyKEEN/PyG/DGL path. Keep separate from tabular experts until edge-split leakage, relation leakage, and negative-edge sampling are audited.",
        },
        {
            "expert": "Hypergraph mechanism model candidate",
            "predicts": "Can drug-target-ADR/tissue/pathway tuples be modeled as hyperedges instead of pairwise labels?",
            "label": "held-out mechanism hyperedges or panel-specific mechanism labels",
            "feature_set": "hyperedge incidence graph, not a tabular feature set",
            "input_policy": "Optional PyG/DGL path. Best for sparse multi-entity mechanism tuples; requires hyperedge split manifests and careful negative hyperedge sampling.",
        },
        {
            "expert": "KAN / DeepADR-style model candidate",
            "predicts": "Can a nonlinear tabular model improve sparse ADR mechanism prediction over logistic/RF/XGBoost?",
            "label": "same label as selected tabular expert",
            "feature_set": "same nonleaky tabular features as the selected expert",
            "input_policy": "Experimental opt-in candidate. Use nested validation/HPO and compare against simple baselines before making claims.",
        },
        {
            "expert": "Expanded tabular model candidates",
            "predicts": "Can regularized linear, additive boosting, or shallow tree boosting improve the selected expert?",
            "label": "same label as selected tabular expert",
            "feature_set": "same nonleaky tabular features as the selected expert",
            "input_policy": "Available model_type values include elastic_net_logistic, explainable_boosting_machine/ebm, catboost, and shallow_lightgbm/lightgbm. Compare under identical split manifests before selecting a claim model.",
        },
    ]
    return pd.DataFrame(rows)


def display_expert_feature_architecture() -> pd.DataFrame:
    """Display the current expert feature architecture when IPython is available."""

    frame = expert_feature_architecture()
    try:
        from IPython.display import display
    except Exception:
        return frame
    display(frame)
    return frame
