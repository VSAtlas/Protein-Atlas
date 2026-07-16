from __future__ import annotations

import json

import pandas as pd
import pytest

from analysis.cli.merge_spd_addon_pipeline import discover_latest_scored_addon
from analysis.external.openfda_pk import _openfda_eligible_mask
from analysis.ml.merge_pk_context import _detect_backfill_columns
from analysis.external.pk_context import (
    combine_pk_context,
    finalize_context,
    join_endpoint_pk_context,
    join_representative_pk_context,
    load_existing_phase1_pk_context,
    load_flat_pk_context,
    load_ncats_frdb_pk_context,
    load_reviewed_openfda_pk_context,
    select_endpoint_pk_context,
    select_representative_pk_context,
    write_pk_context_outputs,
)


def test_ncats_frdb_parser_preserves_context_and_converts_units(tmp_path) -> None:
    source = tmp_path / "frdb-pk.tsv"
    pd.DataFrame(
        [
            {
                "id": 1,
                "compound_id": 3,
                "pk_application_pt": "AMITRIPTYLINE",
                "pk_analyte_pt": "AMITRIPTYLINE",
                "pk_analyte_mw": 277.4033,
                "pk_application_unii": "1806D8D52K",
                "pk_analyte_unii": "1806D8D52K",
                "pk_age_group": "ADULT",
                "pk_health_status": "HEALTHY",
                "pk_species": "Homo sapiens",
                "pk_dose_value": 75,
                "pk_dose_units": "mg",
                "pk_dose_type": "RECOMMENDED",
                "pk_routes": "Oral",
                "pk_experiment_type": "SINGLE",
                "pk_cmax_value": 15.3,
                "pk_cmax_units": "ng/mL",
                "pk_funbound_value": 7.7,
                "pk_source_uri": "https://pubmed.ncbi.nlm.nih.gov/10383563",
            },
            {
                "id": 2,
                "compound_id": 4,
                "pk_application_pt": "EXAMPLE",
                "pk_analyte_pt": "EXAMPLE",
                "pk_analyte_mw": 100,
                "pk_cmax_value": 5,
                "pk_cmax_units": "ng/mL/kg",
            },
        ]
    ).to_csv(source, sep="\t", index=False)

    context = load_ncats_frdb_pk_context(source)

    first = context.loc[context["source_record_id"].astype(str).eq("1")].iloc[0]
    assert first["drug_id"] == "frdb:3"
    assert first["dose_value"] == 75
    assert first["route"] == "Oral"
    assert first["regimen"] == "SINGLE"
    assert first["fraction_unbound_plasma"] == pytest.approx(0.077)
    assert first["cmax_um"] == pytest.approx(15.3 / 277.4033)
    assert first["free_cmax_um"] == pytest.approx(15.3 / 277.4033 * 0.077)

    unsupported = context.loc[context["source_record_id"].astype(str).eq("2")].iloc[0]
    assert pd.isna(unsupported["cmax_um"])
    assert "cmax_unit_not_safely_convertible" in unsupported["context_status"]


def test_reviewed_openfda_loader_excludes_unaccepted_contexts(tmp_path) -> None:
    source = tmp_path / "reviewed.csv"
    pd.DataFrame(
        [
            {
                "drug_id": "selegiline",
                "source_record_id": "accepted",
                "source_url": "https://dailymed.nlm.nih.gov/example",
                "spl_version": "1",
                "acceptable_for_model_training": True,
                "adjudication_status": "accept_model_context",
                "adjudicated_dose_value": 10,
                "adjudicated_dose_unit": "mg",
                "adjudicated_route": "oral",
                "adjudicated_regimen": "single dose",
                "adjudicated_context": "parent drug",
                "cmax_values_raw": "1.0",
                "cmax_units_raw": "ng/mL",
                "cmax_converted_um_candidates": "0.0053",
            },
            {
                "drug_id": "unsafe",
                "source_record_id": "rejected",
                "acceptable_for_model_training": False,
                "adjudication_status": "reject",
                "cmax_converted_um_candidates": "99",
            },
        ]
    ).to_csv(source, index=False)

    context = load_reviewed_openfda_pk_context(source)

    assert context["source_record_id"].tolist() == ["accepted"]
    assert context["cmax_um"].iloc[0] == pytest.approx(0.0053)
    assert context["source_confidence"].iloc[0] == "high"


def test_latest_addon_discovery_skips_newer_incomplete_directory(tmp_path) -> None:
    complete = tmp_path / "target_positive_addon_20260713_090333"
    reference = complete / "reference_vina_full"
    reference.mkdir(parents=True)
    selected = complete / "selected.csv"
    scores = complete / "scores.csv"
    pd.DataFrame([{"pdb_id": "1ABC", "ligand_base": "rdk_0000001"}]).to_csv(
        selected, index=False
    )
    pd.DataFrame([{"pdb_id": "1ABC", "ligand_base": "rdk_0000001"}]).to_csv(
        scores, index=False
    )
    (reference / "reference_vina_manifest.json").write_text(
        json.dumps(
            {
                "selected_pairs": str(selected),
                "outputs": {"scores": str(scores)},
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "target_positive_addon_20260713_093646").mkdir()

    chosen, manifest = discover_latest_scored_addon(tmp_path)

    assert chosen == complete.resolve()
    assert manifest["mode"] == "latest_completed_score_ready"
    assert manifest["rejected_newer_candidates"]


def test_pk_join_accepts_explicitly_missing_inchikey(tmp_path) -> None:
    source = tmp_path / "reviewed.csv"
    pd.DataFrame(
        [
            {
                "drug_id": "selegiline",
                "source_record_id": "accepted",
                "acceptable_for_model_training": True,
                "adjudication_status": "accept_model_context",
                "cmax_converted_um_candidates": "0.0053",
            }
        ]
    ).to_csv(source, index=False)
    context = load_reviewed_openfda_pk_context(source)
    representative = select_representative_pk_context(context)
    model = pd.DataFrame(
        [
            {
                "drug_id": "rdk_0000001",
                "generic_name": "selegiline",
                "pk_context_route": "stale",
            }
        ]
    )

    joined = join_representative_pk_context(model, representative)

    assert joined.columns.is_unique
    assert joined["pk_context_join_status"].iloc[0] == "exact_name"
    assert joined["pk_context_cmax_um"].iloc[0] == pytest.approx(0.0053)


def test_representative_pk_prefers_measured_external_over_empty_placeholder(
    tmp_path,
) -> None:
    model = pd.DataFrame(
        [
            {
                "drug_id": "rdk_0000001",
                "generic_name": "selegiline",
                "cmax_um": pd.NA,
                "fraction_unbound_plasma": pd.NA,
                "free_cmax_um": pd.NA,
            }
        ]
    )
    source = tmp_path / "reviewed.csv"
    pd.DataFrame(
        [
            {
                "drug_id": "selegiline",
                "source_record_id": "accepted",
                "acceptable_for_model_training": True,
                "adjudication_status": "accept_model_context",
                "cmax_converted_um_candidates": "0.0053",
                "adjudicated_dose_value": 10,
                "adjudicated_dose_unit": "mg",
            }
        ]
    ).to_csv(source, index=False)
    context = combine_pk_context(
        [
            load_existing_phase1_pk_context(model),
            load_reviewed_openfda_pk_context(source),
        ]
    )

    representative = select_representative_pk_context(context)

    assert representative["source_name"].iloc[0] == "DailyMed_openFDA_SPL"
    assert representative["cmax_um"].iloc[0] == pytest.approx(0.0053)


def test_openfda_query_filter_excludes_non_fda_and_probe_rows() -> None:
    frame = pd.DataFrame(
        [
            {
                "drug_id": "approved",
                "canonical_identity_regulatory_status": "drugcentral_fda_approved",
                "probe_sensitivity_only": False,
            },
            {
                "drug_id": "not-approved",
                "canonical_identity_regulatory_status": "fda_not_verified",
                "probe_sensitivity_only": False,
            },
            {
                "drug_id": "probe",
                "canonical_identity_regulatory_status": "drugcentral_fda_approved",
                "probe_sensitivity_only": True,
            },
        ]
    )

    eligible, policy = _openfda_eligible_mask(frame)

    assert policy == "canonical_fda_status_or_primary_claim"
    assert frame.loc[eligible, "drug_id"].tolist() == ["approved"]


def test_endpoint_join_preserves_clearance_separately_from_cmax() -> None:
    context = finalize_context(
        pd.DataFrame(
            [
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL_semantic_review",
                    "source_record_id": "clearance",
                    "measurement_context": "systemic_clearance",
                    "clearance_value": 2.5,
                    "clearance_unit": "L/h",
                    "route": "intravenous",
                    "source_confidence": "high",
                    "training_allowed": True,
                },
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL",
                    "source_record_id": "unreviewed",
                    "measurement_context": "plasma_clearance",
                    "clearance_value": 99.0,
                    "clearance_unit": "L/h",
                    "source_confidence": "low",
                    "training_allowed": False,
                },
            ]
        )
    )

    selected = select_endpoint_pk_context(
        context,
        value_column="clearance_value",
        unit_column="clearance_unit",
    )
    joined = join_endpoint_pk_context(
        pd.DataFrame([{"drug_id": "rdk_1", "generic_name": "example"}]),
        selected,
        prefix="pk_clearance_",
    )

    assert len(selected) == 1
    assert selected["clearance_value"].iloc[0] == pytest.approx(2.5)
    assert joined["pk_clearance_clearance_value"].iloc[0] == pytest.approx(2.5)
    assert joined["pk_clearance_join_status"].iloc[0] == "exact_name"


def test_endpoint_selection_partitions_clearance_semantics() -> None:
    context = finalize_context(
        pd.DataFrame(
            [
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL_semantic_review",
                    "source_record_id": "renal",
                    "measurement_context": "renal_clearance",
                    "clearance_value": 1.0,
                    "clearance_unit": "L/h",
                    "route": "oral",
                    "source_confidence": "high",
                    "training_allowed": True,
                },
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL_semantic_review",
                    "source_record_id": "systemic",
                    "measurement_context": "systemic_clearance",
                    "clearance_value": 2.0,
                    "clearance_unit": "L/h",
                    "route": "intravenous",
                    "source_confidence": "high",
                    "training_allowed": True,
                },
            ]
        )
    )

    selected = select_endpoint_pk_context(
        context,
        value_column="clearance_value",
        unit_column="clearance_unit",
    )

    assert set(selected["measurement_context"]) == {
        "renal_clearance",
        "systemic_clearance",
    }


def test_endpoint_selection_excludes_conflicting_same_context() -> None:
    context = finalize_context(
        pd.DataFrame(
            [
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL_semantic_review",
                    "source_record_id": "one",
                    "measurement_context": "renal_clearance",
                    "clearance_value": 1.0,
                    "clearance_unit": "L/h",
                    "route": "oral",
                    "population": "healthy adults",
                    "source_confidence": "high",
                    "training_allowed": True,
                },
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL_semantic_review",
                    "source_record_id": "two",
                    "measurement_context": "renal_clearance",
                    "clearance_value": 3.0,
                    "clearance_unit": "L/h",
                    "route": "oral",
                    "population": "healthy adults",
                    "source_confidence": "high",
                    "training_allowed": True,
                },
            ]
        )
    )

    selected = select_endpoint_pk_context(
        context,
        value_column="clearance_value",
        unit_column="clearance_unit",
    )

    assert selected.empty


def test_pk_merge_auto_detects_endpoint_specific_columns() -> None:
    frame = pd.DataFrame(
        {
            "drug_id": ["example"],
            "pk_context_cmax_um": [1.0],
            "pk_clearance_renal_clearance_value": [2.0],
            "pk_bioavailability_absolute_bioavailability_value": [50.0],
            "unrelated": [3.0],
        }
    )

    selected = _detect_backfill_columns(frame, None)

    assert "pk_context_cmax_um" in selected
    assert "pk_clearance_renal_clearance_value" in selected
    assert "pk_bioavailability_absolute_bioavailability_value" in selected
    assert "unrelated" not in selected


def test_reviewed_cmax_does_not_fall_back_to_unadjudicated_dose(tmp_path) -> None:
    source = tmp_path / "reviewed_without_bound_dose.csv"
    pd.DataFrame(
        [
            {
                "drug_id": "example",
                "source_record_id": "accepted",
                "acceptable_for_model_training": True,
                "adjudication_status": "accept_model_context",
                "dose_value": 999,
                "dose_unit": "mg",
                "cmax_converted_um_candidates": "0.5",
            }
        ]
    ).to_csv(source, index=False)

    context = load_reviewed_openfda_pk_context(source)

    assert context["cmax_um"].iloc[0] == pytest.approx(0.5)
    assert pd.isna(context["dose_value"].iloc[0])
    assert pd.isna(context["dose_context_type"].iloc[0])


def test_flat_loader_quarantines_unspecified_bioavailability(tmp_path) -> None:
    source = tmp_path / "generic_pk.csv"
    pd.DataFrame(
        [
            {
                "drug_name": "example",
                "bioavailability": 80,
                "bioavailability_unit": "%",
                "training_allowed": True,
            }
        ]
    ).to_csv(source, index=False)

    context = load_flat_pk_context(
        source,
        source_name="example_source",
    )

    assert pd.isna(context["bioavailability_value"].iloc[0])
    assert "generic_bioavailability_quarantined" in str(
        context["context_status"].iloc[0]
    )


def test_pk_writer_emits_separate_clearance_endpoint_columns(tmp_path) -> None:
    context = finalize_context(
        pd.DataFrame(
            [
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL_semantic_review",
                    "source_record_id": "renal",
                    "measurement_context": "renal_clearance",
                    "clearance_value": 1.0,
                    "clearance_unit": "L/h",
                    "route": "oral",
                    "source_confidence": "high",
                    "training_allowed": True,
                },
                {
                    "drug_name": "example",
                    "source_name": "DailyMed_openFDA_SPL_semantic_review",
                    "source_record_id": "systemic",
                    "measurement_context": "systemic_clearance",
                    "clearance_value": 2.0,
                    "clearance_unit": "L/h",
                    "route": "intravenous",
                    "source_confidence": "high",
                    "training_allowed": True,
                },
            ]
        )
    )
    model = pd.DataFrame([{"drug_id": "rdk_1", "generic_name": "example"}])

    manifest = write_pk_context_outputs(
        context=context,
        model_table=model,
        out_dir=tmp_path,
    )
    enriched = pd.read_csv(tmp_path / "AtlasSPD_phase1_pk_enriched.csv")

    assert manifest["validation"]["status"] == "pass"
    assert enriched["pk_clearance_renal_clearance_value"].iloc[0] == pytest.approx(1.0)
    assert enriched["pk_clearance_systemic_clearance_value"].iloc[0] == pytest.approx(
        2.0
    )
    assert "pk_clearance_clearance_value" not in enriched


def test_primary_selector_quarantines_explicit_training_denial() -> None:
    context = finalize_context(
        pd.DataFrame(
            [
                {
                    "drug_name": "example",
                    "source_name": "restricted_source",
                    "source_record_id": "restricted",
                    "cmax_um": 10.0,
                    "dose_value": 100.0,
                    "dose_unit": "mg",
                    "dose_context_type": "cmax_study_matched",
                    "source_confidence": "high",
                    "training_allowed": False,
                }
            ]
        )
    )

    selected = select_representative_pk_context(context)

    assert pd.isna(selected["cmax_um"].iloc[0])
    assert pd.isna(selected["dose_value"].iloc[0])
    assert "training_not_allowed" in str(selected["context_status"].iloc[0])
