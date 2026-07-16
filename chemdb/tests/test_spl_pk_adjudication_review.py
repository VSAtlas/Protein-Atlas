from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from analysis.external.spl_pk_adjudication import (
    ACCEPTED_CLEARANCE_NAME,
    ALL_ADJUDICATED_NAME,
    MAX_DOSE_SENSITIVITY_NAME,
    REVIEW_TEMPLATE_NAME,
    adjudicate_spl_pk_candidates,
)


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _candidate(
    *,
    candidate_id: str,
    candidate_type: str,
    text: str,
    value_text: str,
    normalized_value: str,
    normalized_unit: str,
    measurement_context: str,
    section: str,
    cache_file: Path,
    cache_sha256: str,
) -> dict[str, str]:
    start = text.index(value_text)
    return {
        "candidate_id": candidate_id,
        "candidate_type": candidate_type,
        "drug_id": "example",
        "cache_file": str(cache_file),
        "cache_file_sha256": cache_sha256,
        "source_name": "DailyMed_openFDA_SPL_cache",
        "source_record_id": "set-1",
        "source_set_id": "set-1",
        "source_document_id": "set-1",
        "spl_version": "1",
        "source_url": "https://example.invalid/set-1",
        "section": section,
        "section_item_index": "0",
        "section_text_sha256": _sha256_text(text),
        "match_start": str(start),
        "match_end": str(start + len(value_text)),
        "source_excerpt": text,
        "raw_value": normalized_value,
        "raw_unit": normalized_unit,
        "normalized_value": normalized_value,
        "normalized_unit": normalized_unit,
        "normalization_status": "normalized_unambiguous",
        "measurement_context": measurement_context,
        "analyte": "example",
        "analyte_status": "context_label_match",
        "label_generic_names": "example",
        "label_substance_names": "example",
        "route": "",
        "population": "",
        "formulation": "",
        "steady_state": "",
        "regimen": "",
        "exclusion_reason": "",
        "extraction_method": "spl_pk_candidate_review_v1",
    }


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def test_review_template_controls_training_and_preserves_sensitivity(
    tmp_path: Path,
) -> None:
    clearance_text = (
        "Following intravenous administration in healthy adult subjects, "
        "example total body clearance was 10 L/h."
    )
    maximum_text = (
        "For adults, the maximum recommended daily dose of example after oral "
        "administration is 800 mg."
    )
    cache_file = tmp_path / "set-1.json"
    cache_file.write_text(
        json.dumps(
            {
                "records": [
                    {
                        "set_id": "set-1",
                        "version": "1",
                        "pharmacokinetics": [clearance_text],
                        "dosage_and_administration": [maximum_text],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    cache_sha256 = hashlib.sha256(cache_file.read_bytes()).hexdigest()
    rows = [
        _candidate(
            candidate_id="clearance-1",
            candidate_type="clearance",
            text=clearance_text,
            value_text="10 L/h",
            normalized_value="10",
            normalized_unit="L/h",
            measurement_context="total_body_clearance",
            section="pharmacokinetics",
            cache_file=cache_file,
            cache_sha256=cache_sha256,
        ),
        _candidate(
            candidate_id="maximum-1",
            candidate_type="maximum_labeled_dose",
            text=maximum_text,
            value_text="800 mg",
            normalized_value="800",
            normalized_unit="mg",
            measurement_context="maximum_daily_dose",
            section="dosage_and_administration",
            cache_file=cache_file,
            cache_sha256=cache_sha256,
        ),
    ]
    candidates = tmp_path / "candidates.csv"
    with candidates.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    first_out = tmp_path / "first"
    first_manifest = adjudicate_spl_pk_candidates(
        candidates=candidates,
        out_dir=first_out,
    )

    template_columns, template_rows = _read_csv(first_out / REVIEW_TEMPLATE_NAME)
    assert len(template_rows) == 2
    assert {
        "canonical_candidate_id",
        "review_decision",
        "reviewer",
        "reviewed_at",
        "review_notes",
        "source_lineage_json",
    }.issubset(template_columns)
    _, unreviewed = _read_csv(first_out / ALL_ADJUDICATED_NAME)
    assert all(row["training_allowed"] == "false" for row in unreviewed)
    assert first_manifest["training_allowed_context_rows"] == 0

    for row in template_rows:
        row["review_decision"] = "approved"
        row["reviewer"] = "local-reviewer"
        row["reviewed_at"] = "2026-07-14"
        row["review_notes"] = "source clause checked"
    decisions = tmp_path / "review_decisions.csv"
    with decisions.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=template_columns)
        writer.writeheader()
        writer.writerows(template_rows)

    reviewed_out = tmp_path / "reviewed"
    reviewed_manifest = adjudicate_spl_pk_candidates(
        candidates=candidates,
        out_dir=reviewed_out,
        review_decisions=decisions,
    )

    _, clearance_rows = _read_csv(reviewed_out / ACCEPTED_CLEARANCE_NAME)
    _, maximum_rows = _read_csv(reviewed_out / MAX_DOSE_SENSITIVITY_NAME)
    assert clearance_rows[0]["training_allowed"] == "true"
    assert (
        clearance_rows[0]["review_decision_file_sha256"]
        == hashlib.sha256(decisions.read_bytes()).hexdigest()
    )
    assert maximum_rows[0]["review_approved"] == "true"
    assert maximum_rows[0]["sensitivity_only"] == "true"
    assert maximum_rows[0]["training_allowed"] == "false"
    assert maximum_rows[0]["review_application_status"] == "approved_sensitivity_only"
    assert reviewed_manifest["training_allowed_context_rows"] == 1
