from __future__ import annotations

import pandas as pd
import pytest

from analysis.ml.source_transfer_audit import audit_source_transfer


def test_source_transfer_audit_writes_balance_shift_prediction_and_findings(tmp_path):
    rows = []
    for idx in range(50):
        source = "chembl;papyrus" if idx < 30 else "toxcast"
        label = 1 if idx < 24 else 0
        if source == "toxcast":
            label = 1 if idx == 31 else 0
        rows.append(
            {
                "drug_id": f"D{idx % 20:02d}",
                "target_id": f"T{idx % 8:02d}",
                "label_source": source,
                "atlas_score": float(idx % 5),
                "consensus_score": float(idx % 7),
                "y": label,
            }
        )
    dataset = tmp_path / "source.csv"
    pred = tmp_path / "pred.csv"
    frame = pd.DataFrame(rows)
    frame.to_csv(dataset, index=False)
    frame.assign(ml_prediction_score=0.8).to_csv(pred, index=False)

    manifest = audit_source_transfer(
        dataset,
        "y",
        tmp_path / "audit",
        train_sources=["chembl", "papyrus"],
        test_sources=["toxcast"],
        predictions_path=pred,
    )

    assert "large_label_prevalence_shift" in manifest["findings"]
    assert manifest["prediction_audit"]["status"] == "written"
    assert (tmp_path / "audit" / "source_transfer_audit_manifest.json").exists()
    assert (tmp_path / "audit" / "prediction_calibration_by_source.csv").exists()


def test_source_transfer_audit_requires_label_source(tmp_path):
    dataset = tmp_path / "bad.csv"
    pd.DataFrame({"drug_id": ["D1"], "y": [1]}).to_csv(dataset, index=False)
    with pytest.raises(ValueError, match="label_source"):
        audit_source_transfer(dataset, "y", tmp_path / "audit", train_sources=["a"], test_sources=["b"])

