from __future__ import annotations

import json

import pandas as pd
import pytest

from analysis.cli.merge_spd_addon_pipeline import discover_latest_scored_addon
from analysis.external.openfda_pk import _openfda_eligible_mask
from analysis.external.pk_context import (
    combine_pk_context,
    join_representative_pk_context,
    load_existing_phase1_pk_context,
    load_ncats_frdb_pk_context,
    load_reviewed_openfda_pk_context,
    select_representative_pk_context,
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

    unsupported = context.loc[
        context["source_record_id"].astype(str).eq("2")
    ].iloc[0]
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
