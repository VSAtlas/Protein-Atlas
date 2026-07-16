from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.pk_context import PK_CONTEXT_COLUMNS, finalize_context


def _text(frame: pd.DataFrame, column: str) -> pd.Series:
    return (
        frame.get(column, pd.Series("", index=frame.index))
        .fillna("")
        .astype(str)
        .str.strip()
    )


def _accepted(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    if not source.is_file() or source.stat().st_size == 0:
        return pd.DataFrame()
    frame = pd.read_csv(source, low_memory=False)
    if "adjudication_status" not in frame.columns:
        return pd.DataFrame()
    accepted = _text(frame, "adjudication_status").eq("accepted") | (
        _text(frame, "review_approved").str.casefold().isin({"1", "true", "yes"})
        & _text(frame, "review_application_status").eq("approved_for_context_training")
    )
    not_review_rejected = ~_text(frame, "review_decision").str.casefold().eq("rejected")
    return frame.loc[accepted & not_review_rejected].copy()


def _explicit_review_approval(frame: pd.DataFrame) -> pd.Series:
    return (
        _text(frame, "canonical_candidate_id").ne("")
        & _text(frame, "review_decision").str.casefold().eq("approved")
        & _text(frame, "review_application_status").eq("approved_for_context_training")
        & _text(frame, "review_approved").str.casefold().isin({"1", "true", "yes"})
        & _text(frame, "review_decision_source").ne("")
        & _text(frame, "review_decision_file_sha256")
        .str.casefold()
        .str.fullmatch(r"[0-9a-f]{64}")
    )


def load_adjudicated_spl_pk_context(
    *,
    clearance_path: str | Path,
    bioavailability_path: str | Path,
) -> pd.DataFrame:
    """Load only same-clause-adjudicated SPL PK contexts.

    Maximum recommended doses are intentionally excluded here. They remain a
    separately named sensitivity artifact and never replace a Cmax-matched dose.
    """

    parts: list[pd.DataFrame] = []
    for measurement, path in (
        ("clearance", clearance_path),
        ("absolute_bioavailability", bioavailability_path),
    ):
        raw = _accepted(path)
        if raw.empty:
            continue
        out = pd.DataFrame(
            {
                column: pd.Series(pd.NA, index=raw.index, dtype="object")
                for column in PK_CONTEXT_COLUMNS
            }
        )
        review_approved = _explicit_review_approval(raw)
        canonical_id = _text(raw, "canonical_candidate_id")
        candidate_id = _text(raw, "candidate_id")
        lineage = _text(raw, "source_lineage_json")
        review_reference = (
            "review_decision_file="
            + _text(raw, "review_decision_source")
            + "; review_decision_sha256="
            + _text(raw, "review_decision_file_sha256")
        )

        out["drug_id"] = raw["drug_id"]
        out["drug_name"] = raw["drug_id"]
        out["source_name"] = "DailyMed_openFDA_SPL"
        out["source_version"] = _text(raw, "spl_versions").where(
            _text(raw, "spl_versions").ne(""),
            _text(raw, "spl_version"),
        )
        out["source_record_id"] = canonical_id.where(canonical_id.ne(""), candidate_id)
        out["source_url"] = _text(raw, "source_url")
        out["study_id"] = _text(raw, "source_record_id")
        out["reference"] = lineage
        out.loc[review_approved, "reference"] = (
            lineage.where(lineage.eq(""), lineage + "; ") + review_reference
        ).loc[review_approved]
        out["population"] = _text(raw, "reviewed_population").where(
            _text(raw, "reviewed_population").ne(""),
            _text(raw, "population_evidence_text"),
        )
        out["species"] = "Homo sapiens"
        out["route"] = _text(raw, "route_evidence_text")
        out["regimen"] = _text(raw, "regimen")
        out["formulation"] = _text(raw, "reviewed_formulation").where(
            _text(raw, "reviewed_formulation").ne(""),
            _text(raw, "formulation"),
        )
        out["parent_or_metabolite"] = "parent"
        out["measurement_context"] = _text(raw, "semantic_measurement_context")
        out["extraction_method"] = "same_clause_machine_semantic_adjudication"
        out.loc[review_approved, "extraction_method"] = (
            "same_clause_machine_semantic_adjudication_with_explicit_review_decision"
        )
        out["source_confidence"] = "medium"
        out.loc[review_approved, "source_confidence"] = "high"
        review_state = _text(raw, "review_application_status").where(
            _text(raw, "review_application_status").ne(""),
            "review_not_provided",
        )
        out["context_status"] = (
            "machine_semantic_accepted; "
            + review_state
            + "; context_specific_not_universal_drug_pk; "
            + _text(raw, "semantic_measurement_context")
        )
        out["license_note"] = "public FDA labeling; retain SPL set/version provenance"
        out["training_allowed"] = review_approved.astype(bool)
        if measurement == "clearance":
            out["clearance_value"] = raw["normalized_value"]
            out["clearance_unit"] = raw["normalized_unit"]
        else:
            out["bioavailability_value"] = raw["normalized_value"]
            out["bioavailability_unit"] = raw["normalized_unit"]
            out["route"] = _text(raw, "reviewed_extravascular_route").where(
                _text(raw, "reviewed_extravascular_route").ne(""),
                _text(raw, "extravascular_route"),
            )
            out["reference_route"] = "intravenous"
        parts.append(out)
    if not parts:
        return pd.DataFrame(columns=["pk_context_id", *PK_CONTEXT_COLUMNS])
    combined = pd.DataFrame(
        {
            column: pd.concat(
                [part[column] for part in parts],
                ignore_index=True,
            )
            for column in PK_CONTEXT_COLUMNS
        }
    )
    return finalize_context(combined)


__all__ = ["load_adjudicated_spl_pk_context"]
