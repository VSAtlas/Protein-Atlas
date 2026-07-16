from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from analysis.external.spl_pk_adjudication import (
    ACCEPTED_BIOAVAILABILITY_NAME,
    REVIEW_TEMPLATE_NAME,
    _maximum_dose_machine_rejection_reviewable,
    _structured_review_errors,
    adjudicate_spl_pk_candidates,
    adjudicate_spl_pk_rows,
)
from analysis.external.spl_pk_context import load_adjudicated_spl_pk_context


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def test_review_can_resolve_ambiguous_absolute_bioavailability(
    tmp_path: Path,
) -> None:
    text = (
        "After oral administration, the absolute bioavailability of example "
        "is approximately 73%."
    )
    cache_file = tmp_path / "set-1.json"
    cache_file.write_text(
        json.dumps(
            {
                "records": [
                    {
                        "set_id": "set-1",
                        "version": "1",
                        "pharmacokinetics": [text],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    start = text.index("73%")
    row = {
        "candidate_id": "bioavailability-1",
        "candidate_type": "absolute_bioavailability",
        "drug_id": "example",
        "cache_file": str(cache_file),
        "cache_file_sha256": hashlib.sha256(cache_file.read_bytes()).hexdigest(),
        "source_name": "DailyMed_openFDA_SPL_cache",
        "source_record_id": "set-1",
        "source_set_id": "set-1",
        "source_document_id": "set-1",
        "spl_version": "1",
        "source_url": "https://example.invalid/set-1",
        "section": "pharmacokinetics",
        "section_item_index": "0",
        "section_text_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "match_start": str(start),
        "match_end": str(start + len("73%")),
        "source_excerpt": text,
        "raw_value": "73",
        "raw_unit": "%",
        "normalized_value": "73",
        "normalized_unit": "%",
        "normalization_status": "normalized_unambiguous",
        "measurement_context": "absolute_bioavailability",
        "analyte": "example",
        "analyte_status": "context_label_match",
        "label_generic_names": "example",
        "label_substance_names": "example",
        "route": "oral",
        "population": "",
        "formulation": "",
        "steady_state": "",
        "regimen": "",
        "exclusion_reason": "",
        "extraction_method": "spl_pk_candidate_review_v1",
    }
    candidates = tmp_path / "candidates.csv"
    with candidates.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)

    first_out = tmp_path / "first"
    adjudicate_spl_pk_candidates(candidates=candidates, out_dir=first_out)
    columns, template = _read_csv(first_out / REVIEW_TEMPLATE_NAME)
    assert len(template) == 1
    assert template[0]["machine_adjudication_status"] == "ambiguous"
    assert template[0]["eligible_for_approval"] == "true"

    template[0].update(
        {
            "review_decision": "approved",
            "reviewer": "source-text-reviewer",
            "reviewed_at": "2026-07-15",
            "review_notes": (
                "Official SPL explicitly reports parent-drug absolute oral "
                "bioavailability; local IV comparator text is unnecessary."
            ),
            "reviewed_value": "73",
            "reviewed_unit": "%",
            "reviewed_value_qualifier": "central_estimate",
            "reviewed_analyte": "example",
            "reviewed_endpoint_semantics": "absolute_bioavailability",
            "reviewed_reference_basis": "explicit_absolute_bioavailability",
            "reviewed_extravascular_route": "oral",
            "reviewed_formulation": "tablet",
            "reviewed_population": "human SPL label context",
        }
    )
    decisions = tmp_path / "review_decisions.csv"
    with decisions.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(template)

    reviewed_out = tmp_path / "reviewed"
    manifest = adjudicate_spl_pk_candidates(
        candidates=candidates,
        out_dir=reviewed_out,
        review_decisions=decisions,
    )
    _, accepted = _read_csv(reviewed_out / ACCEPTED_BIOAVAILABILITY_NAME)
    assert len(accepted) == 1
    assert accepted[0]["machine_adjudication_status"] == "ambiguous"
    assert accepted[0]["review_approved"] == "true"
    assert accepted[0]["training_allowed"] == "true"
    assert accepted[0]["normalized_value"] == "73"
    assert accepted[0]["endpoint_semantics"] == (
        "absolute_bioavailability:oral_vs_intravenous"
    )
    assert manifest["training_allowed_context_rows"] == 1

    normalized = load_adjudicated_spl_pk_context(
        clearance_path=tmp_path / "missing_clearance.csv",
        bioavailability_path=reviewed_out / ACCEPTED_BIOAVAILABILITY_NAME,
    )
    assert len(normalized) == 1
    assert normalized.iloc[0]["route"] == "oral"
    assert normalized.iloc[0]["formulation"] == "tablet"
    assert normalized.iloc[0]["population"] == "human SPL label context"
    assert bool(normalized.iloc[0]["training_allowed"])


def test_absolute_value_with_relative_comparator_is_reviewable(
    tmp_path: Path,
) -> None:
    text = (
        "The mean absolute bioavailability of example extended-release "
        "tablets is 90%, and the relative bioavailability is about 100%."
    )
    cache_file = tmp_path / "set-relative.json"
    cache_file.write_text(
        json.dumps(
            {
                "records": [
                    {
                        "set_id": "set-relative",
                        "version": "1",
                        "pharmacokinetics": [text],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    start = text.index("90%")
    row = {
        "candidate_id": "bioavailability-relative",
        "candidate_type": "absolute_bioavailability",
        "drug_id": "example",
        "cache_file": str(cache_file),
        "cache_file_sha256": hashlib.sha256(cache_file.read_bytes()).hexdigest(),
        "source_name": "DailyMed_openFDA_SPL_cache",
        "source_record_id": "set-relative",
        "source_set_id": "set-relative",
        "source_document_id": "set-relative",
        "spl_version": "1",
        "source_url": "https://example.invalid/set-relative",
        "section": "pharmacokinetics",
        "section_item_index": "0",
        "section_text_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "match_start": str(start),
        "match_end": str(start + len("90%")),
        "source_excerpt": text,
        "raw_value": "90",
        "raw_unit": "%",
        "normalized_value": "90",
        "normalized_unit": "%",
        "normalization_status": "normalized_unambiguous",
        "measurement_context": "absolute_bioavailability",
        "analyte": "example",
        "analyte_status": "context_label_match",
        "label_generic_names": "example",
        "label_substance_names": "example",
        "route": "oral",
        "population": "healthy subjects",
        "formulation": "extended-release tablet",
        "steady_state": "",
        "regimen": "",
        "exclusion_reason": "",
        "extraction_method": "spl_pk_candidate_review_v1",
    }
    candidates = tmp_path / "relative_candidates.csv"
    with candidates.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)

    out_dir = tmp_path / "relative"
    adjudicate_spl_pk_candidates(candidates=candidates, out_dir=out_dir)
    _, template = _read_csv(out_dir / REVIEW_TEMPLATE_NAME)
    assert len(template) == 1
    assert template[0]["machine_adjudication_status"] == "ambiguous"
    assert template[0]["eligible_for_approval"] == "true"


def test_mixed_adult_pediatric_maximum_is_reviewable_not_rejected() -> None:
    rows, _ = adjudicate_spl_pk_rows(
        [
            {
                "candidate_id": "maximum-mixed",
                "candidate_type": "maximum_labeled_dose",
                "drug_id": "example",
                "section": "dosage_and_administration",
                "normalized_value": "10",
                "normalized_unit": "mg",
                "normalization_status": "normalized_unambiguous",
                "measurement_context": "maximum_daily_dose",
                "semantic_measurement_context": "maximum_daily_dose",
                "semantic_context_status": "verified",
                "semantic_clause": (
                    "For adults, the maximum recommended daily oral dose is "
                    "10 mg. Pediatric patients use a lower dose."
                ),
                "route": "oral",
                "source_excerpt": "Human dosage and administration.",
            }
        ]
    )

    assert rows[0]["adjudication_status"] == "ambiguous"
    assert (
        "mixed_adult_pediatric_maximum_requires_structured_review"
        in rows[0]["adjudication_reasons"]
    )
    assert "pediatric_maximum_dose_not_allowed" not in rows[0]["adjudication_reasons"]


def test_only_parser_limited_maximum_rejections_are_reviewable() -> None:
    base = {
        "candidate_type": "maximum_labeled_dose",
        "machine_adjudication_reasons": (
            "maximum_dose_interval_unspecified;"
            "same_clause_semantics_not_verified:multiple_dose_values_in_clause;"
            "parent_analyte_not_in_value_clause"
        ),
    }
    assert _maximum_dose_machine_rejection_reviewable(base)

    base["machine_adjudication_reasons"] += ";pediatric_maximum_dose_not_allowed"
    assert not _maximum_dose_machine_rejection_reviewable(base)


def test_ambiguous_maximum_dose_review_is_sensitivity_only() -> None:
    row = {
        "candidate_type": "maximum_labeled_dose",
        "drug_id": "example",
        "reviewed_value": "10",
        "reviewed_unit": "mg",
        "reviewed_value_qualifier": "recommended_maximum",
        "reviewed_analyte": "example",
        "reviewed_endpoint_semantics": "maximum_daily_dose",
        "reviewed_reference_basis": "explicit_recommended_adult_maximum",
        "reviewed_extravascular_route": "oral",
        "reviewed_population": "adult patients",
    }
    assert _structured_review_errors(row) == []

    row["reviewed_population"] = "pediatric patients"
    assert "reviewed_population_not_explicitly_adult" in _structured_review_errors(row)


def test_ambiguous_bioavailability_review_requires_structured_fields() -> None:
    row = {
        "candidate_type": "absolute_bioavailability",
        "drug_id": "example",
        "reviewed_value": "73",
        "reviewed_unit": "%",
        "reviewed_value_qualifier": "central_estimate",
        "reviewed_analyte": "example",
        "reviewed_endpoint_semantics": "absolute_bioavailability",
        "reviewed_reference_basis": "explicit_absolute_bioavailability",
        "reviewed_extravascular_route": "",
    }
    assert "missing_reviewed_extravascular_route" in _structured_review_errors(row)

    row["reviewed_extravascular_route"] = "oral"
    assert _structured_review_errors(row) == []
