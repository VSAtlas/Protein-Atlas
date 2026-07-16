from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.spl_pk_context import load_adjudicated_spl_pk_context


def test_maximum_dose_is_not_loaded_as_primary_pk_context(tmp_path: Path) -> None:
    clearance = tmp_path / "clearance.csv"
    bioavailability = tmp_path / "bioavailability.csv"
    pd.DataFrame(
        [
            {
                "adjudication_status": "accepted",
                "drug_id": "example",
                "normalized_value": 2.5,
                "normalized_unit": "L/h",
                "semantic_measurement_context": "systemic_clearance",
                "route_evidence_text": "intravenous",
                "population_evidence_text": "adult subjects",
                "canonical_candidate_id": "candidate-1",
            }
        ]
    ).to_csv(clearance, index=False)
    pd.DataFrame(columns=["adjudication_status"]).to_csv(bioavailability, index=False)

    context = load_adjudicated_spl_pk_context(
        clearance_path=clearance,
        bioavailability_path=bioavailability,
    )

    assert len(context) == 1
    assert context.iloc[0]["clearance_value"] == 2.5
    assert pd.isna(context.iloc[0]["dose_value"])
    assert context.iloc[0]["measurement_context"] == "systemic_clearance"
    assert not bool(context.iloc[0]["training_allowed"])
    assert context.iloc[0]["source_name"] == "DailyMed_openFDA_SPL"
    assert "review" not in str(context.iloc[0]["extraction_method"])


def test_absolute_bioavailability_uses_extravascular_test_route(tmp_path: Path) -> None:
    clearance = tmp_path / "clearance.csv"
    bioavailability = tmp_path / "bioavailability.csv"
    pd.DataFrame(columns=["adjudication_status"]).to_csv(clearance, index=False)
    pd.DataFrame(
        [
            {
                "adjudication_status": "accepted",
                "drug_id": "example",
                "normalized_value": 50.0,
                "normalized_unit": "%",
                "semantic_measurement_context": "absolute_bioavailability",
                "route_evidence_text": "intravenous",
                "extravascular_route": "oral",
                "canonical_candidate_id": "candidate-2",
                "review_decision": "approved",
                "review_application_status": "approved_for_context_training",
                "review_approved": True,
                "review_decision_source": "local/review_decisions.csv",
                "review_decision_file_sha256": "a" * 64,
            }
        ]
    ).to_csv(bioavailability, index=False)

    context = load_adjudicated_spl_pk_context(
        clearance_path=clearance,
        bioavailability_path=bioavailability,
    )

    assert context.iloc[0]["route"] == "oral"
    assert context.iloc[0]["reference_route"] == "intravenous"
    assert bool(context.iloc[0]["training_allowed"])
    assert "explicit_review_decision" in str(context.iloc[0]["extraction_method"])
