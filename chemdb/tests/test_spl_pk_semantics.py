from __future__ import annotations

from analysis.external.spl_pk_semantics import semantic_evidence_for_text


def _row(candidate_type: str, drug_id: str, text: str, value: str) -> dict[str, object]:
    start = text.index(value)
    return {
        "candidate_type": candidate_type,
        "drug_id": drug_id,
        "match_start": start,
        "match_end": start + len(value),
        "raw_value": value.split()[0],
        "section": "pharmacokinetics",
    }


def test_creatinine_clearance_cutoff_is_not_drug_clearance() -> None:
    text = (
        "In adult patients, meropenem plasma clearance correlates with creatinine "
        "clearance 50 mL/min or less after intravenous administration."
    )
    evidence = semantic_evidence_for_text(
        text,
        _row("clearance", "meropenem", text, "50 mL/min"),
    )
    assert evidence["semantic_context_status"] == "not_verified"
    assert "subject_renal_function_not_drug_clearance" in str(
        evidence["semantic_context_reasons"]
    )


def test_direct_human_clearance_context_is_verified() -> None:
    text = (
        "Following intravenous administration in healthy adult subjects, eliglustat "
        "total body clearance was 88 L/h."
    )
    evidence = semantic_evidence_for_text(
        text,
        _row("clearance", "eliglustat", text, "88 L/h"),
    )
    assert evidence["semantic_context_status"] == "verified"
    assert evidence["adjudicated_analyte"] == "eliglustat"


def test_absolute_bioavailability_requires_same_analyte() -> None:
    text = (
        "In healthy subjects, the absolute bioavailability of citalopram is 80% "
        "relative to an intravenous dose after oral administration."
    )
    evidence = semantic_evidence_for_text(
        text,
        _row("absolute_bioavailability", "escitalopram", text, "80%"),
    )
    assert evidence["semantic_context_status"] == "not_verified"
    assert "parent_analyte_not_in_value_clause" in str(
        evidence["semantic_context_reasons"]
    )


def test_maximum_tolerated_dose_stays_out_of_sensitivity_context() -> None:
    text = (
        "In adult subjects, the maximum tolerated dose of imatinib was 800 mg "
        "during dose-escalation after oral administration."
    )
    row = _row("maximum_labeled_dose", "imatinib", text, "800 mg")
    row["section"] = "dosage_and_administration"
    evidence = semantic_evidence_for_text(text, row)
    assert evidence["semantic_context_status"] == "not_verified"
    assert "highest_studied_or_mtd_not_recommended_maximum" in str(
        evidence["semantic_context_reasons"]
    )


def test_stratified_absolute_bioavailability_requires_manual_context() -> None:
    text = (
        "In healthy adult subjects after oral and intravenous administration, "
        "the absolute bioavailability of tramadol was 73% in males and 79% in females."
    )
    evidence = semantic_evidence_for_text(
        text,
        _row("absolute_bioavailability", "tramadol", text, "73%"),
    )

    assert evidence["semantic_context_status"] == "not_verified"
    assert "multiple_bioavailability_values_in_clause" in str(
        evidence["semantic_context_reasons"]
    )


def test_clearance_route_in_neighboring_clause_is_not_evidence() -> None:
    text = (
        "The drug was administered intravenously. In healthy adult subjects, "
        "example total body clearance was 10 L/h."
    )
    evidence = semantic_evidence_for_text(
        text,
        _row("clearance", "example", text, "10 L/h"),
    )

    assert evidence["semantic_context_status"] == "not_verified"
    assert "administration_route_not_in_value_clause" in str(
        evidence["semantic_context_reasons"]
    )
    assert evidence["route_evidence_text"] == ""
    assert evidence["route_evidence_json"] == "[]"


def test_bioavailability_routes_in_neighboring_clause_are_not_evidence() -> None:
    text = (
        "Healthy subjects received oral and intravenous doses. In healthy adult "
        "subjects, the absolute bioavailability of example was 50%."
    )
    evidence = semantic_evidence_for_text(
        text,
        _row("absolute_bioavailability", "example", text, "50%"),
    )

    reasons = str(evidence["semantic_context_reasons"])
    assert evidence["semantic_context_status"] == "not_verified"
    assert "intravenous_reference_not_in_value_clause" in reasons
    assert "extravascular_test_route_not_in_value_clause" in reasons
    assert evidence["route_evidence_json"] == "[]"


def test_maximum_dose_requires_recommended_adult_same_clause() -> None:
    text = "The maximum daily dose of example after oral administration is 800 mg."
    row = _row("maximum_labeled_dose", "example", text, "800 mg")
    row["section"] = "dosage_and_administration"

    evidence = semantic_evidence_for_text(text, row)

    reasons = str(evidence["semantic_context_reasons"])
    assert evidence["semantic_context_status"] == "not_verified"
    assert "explicit_adult_context_not_in_value_clause" in reasons
    assert "explicit_recommended_dose_context_not_in_value_clause" in reasons


def test_decimal_value_does_not_split_semantic_clause() -> None:
    text = (
        "In healthy adult subjects after oral administration, the absolute "
        "bioavailability of ibrutinib was 2.9% (90% CI: 2.7% to 3.1%)."
    )
    evidence = semantic_evidence_for_text(
        text,
        _row("absolute_bioavailability", "ibrutinib", text, "2.9%"),
    )

    assert "2.9%" in str(evidence["semantic_clause"])
    assert "90% CI" in str(evidence["semantic_clause"])
