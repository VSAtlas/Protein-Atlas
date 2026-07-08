from __future__ import annotations

import json

import pandas as pd
import pytest

from analysis.ml.external_eval import run_external_model_eval
from analysis.ml.train_classifier import train_ml_model


def _training_rows(n: int = 80) -> list[dict[str, object]]:
    rows = []
    for idx in range(n):
        rows.append(
            {
                "drug_id": f"D{idx % 20:02d}",
                "target_id": f"T{idx % 8:02d}",
                "label_source": "chembl" if idx % 3 else "papyrus",
                "source_family": "curated",
                "scaffold_key": f"S{idx % 10:02d}",
                "chemical_cluster": f"C{idx % 9:02d}",
                "target_family": f"F{idx % 4:02d}",
                "protein_class": f"P{idx % 3:02d}",
                "activity_publication_year": 2010 + (idx % 10),
                "atlas_score": float(idx % 11) / 10.0,
                "consensus_score": float((idx * 3) % 13) / 10.0,
                "y": int(idx % 2 == 0),
                "validation_fold": "validation" if idx % 7 == 0 else "train",
            }
        )
    return rows


def test_train_ml_model_writes_publishable_artifact_contract(tmp_path):
    dataset = tmp_path / "model_ready.csv"
    pd.DataFrame(_training_rows()).to_csv(dataset, index=False)
    out = tmp_path / "model"

    result = train_ml_model(
        dataset,
        "y",
        "pilot_nonleaky",
        "logistic_regression",
        "random",
        out,
        seed=7,
        n_bootstraps=2,
        validation_fold_col="validation_fold",
        validation_fold_value="validation",
        nested_model_selection=True,
        repo_root=tmp_path,
    )

    required = [
        "dataset_version_manifest.json",
        "model_predictions.csv",
        "model_metrics.csv",
        "model_metric_bootstrap_ci.csv",
        "model_reliability_table.csv",
        "model_grouped_calibration.csv",
        "model_subgroup_metrics.csv",
        "model_decision_metrics.csv",
        "standard_baseline_panel.csv",
        "feature_importance.csv",
        "split_manifest.csv",
        "split_manifest.json",
        "split_summary.json",
        "hpo_trials.csv",
        "hpo_manifest.json",
        "conformal_summary.json",
        "pu_training_manifest.json",
        "model_claim_readiness.json",
        "model_card.json",
        "model_artifacts.json",
        "model_preprocessing.json",
        "training_design_matrix.csv",
        "registry.json",
        "trained_model.pkl",
    ]
    for name in required:
        assert (out / name).exists(), name

    assert result["features"] == ["atlas_score", "consensus_score"]
    pred = pd.read_csv(out / "model_predictions.csv")
    assert {"drug_id", "target_id", "label_source", "ml_prediction_score", "y"}.issubset(pred.columns)
    preprocessing = json.loads((out / "model_preprocessing.json").read_text())
    assert preprocessing["raw_features"] == ["atlas_score", "consensus_score"]
    assert preprocessing["design_columns"] == ["atlas_score", "consensus_score"]
    artifact_manifest = json.loads((out / "model_artifacts.json").read_text())
    assert any(item["name"] == "model_card.json" for item in artifact_manifest["artifacts"])
    assert any(item["name"] == "model_preprocessing.json" for item in artifact_manifest["artifacts"])



def test_external_eval_replays_training_preprocessing_and_writes_uncertainty_outputs(tmp_path):
    dataset = tmp_path / "model_ready.csv"
    pd.DataFrame(_training_rows()).to_csv(dataset, index=False)
    model_dir = tmp_path / "model"
    train_ml_model(
        dataset,
        "y",
        "pilot_nonleaky",
        "logistic_regression",
        "random",
        model_dir,
        seed=11,
        n_bootstraps=1,
        validation_fold_col="validation_fold",
        validation_fold_value="validation",
        nested_model_selection=True,
        repo_root=tmp_path,
    )

    external_rows = _training_rows(24)
    external_rows[0]["atlas_score"] = None
    external = tmp_path / "external.csv"
    pd.DataFrame(external_rows).to_csv(external, index=False)
    out = tmp_path / "external_eval"

    manifest = run_external_model_eval(
        model_dir=model_dir,
        dataset=external,
        label_col="y",
        name="contract_external",
        out_dir=out,
    )

    assert manifest["status"] == "evaluated"
    assert manifest["calibration_status"] == "ok"
    assert manifest["conformal_status"] == "ok"
    assert manifest["applicability_domain_status"] == "ok"
    pred = pd.read_csv(out / "external_predictions.csv")
    assert {"ml_prediction_score", "conformal_prediction_set", "applicability_domain"}.issubset(pred.columns)
    assert (out / "external_calibration_summary.json").exists()


def test_external_eval_blocks_missing_raw_feature_instead_of_zero_filling(tmp_path):
    dataset = tmp_path / "model_ready.csv"
    pd.DataFrame(_training_rows()).to_csv(dataset, index=False)
    model_dir = tmp_path / "model"
    train_ml_model(
        dataset,
        "y",
        "pilot_nonleaky",
        "logistic_regression",
        "random",
        model_dir,
        seed=13,
        n_bootstraps=1,
        validation_fold_col="validation_fold",
        validation_fold_value="validation",
        nested_model_selection=True,
        repo_root=tmp_path,
    )

    external = tmp_path / "external_missing.csv"
    pd.DataFrame(_training_rows(20)).drop(columns=["consensus_score"]).to_csv(external, index=False)
    out = tmp_path / "blocked_external_eval"

    with pytest.raises(ValueError, match="missing raw feature columns"):
        run_external_model_eval(
            model_dir=model_dir,
            dataset=external,
            label_col="y",
            name="missing_feature_external",
            out_dir=out,
        )

    blocked = json.loads((out / "external_eval_manifest.json").read_text())
    assert blocked["status"] == "blocked"
    assert blocked["reason"] == "missing_raw_features"
    assert blocked["details"]["missing_raw_features"] == ["consensus_score"]
