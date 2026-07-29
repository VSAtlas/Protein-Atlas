from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from analysis.cli.audit_spd_interval_parser import (
    BEHAVIOR_COLUMNS,
    BOUNDARY_COLUMNS,
    RELATION_COLUMNS,
    _boundary_cases,
    _compare_materialized_expert,
    _compare_pair_outputs,
    _relation_summary,
)
from analysis.cli.verify_spd_boundary_correction import (
    EVIDENCE_COLUMNS,
    _decimal,
    _production_exposure_label,
    _qualifier_is_bound,
)
from analysis.external.spd import aggregate_spd_assays
from analysis.ml.spd_four_expert_tables import (
    BINDING_LABEL_POLICY_VERSION,
    _binding_label,
)
from analysis.ml.spd_label_enrichment import _read_spd_panel


def _synthetic_observations() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "drug_id": ["drug_a", "drug_b"],
            "target_id": ["TARGET_A", "TARGET_B"],
            "assay_id": ["a1", "b1"],
            "assay_name": ["assay a", "assay b"],
            "assay_type": ["biochemical", "target"],
            "ac50_nM": [500.0, 20_000.0],
            "free_cmax_nM": [100.0, 1_000.0],
            "total_cmax_nM": [200.0, 2_000.0],
            "activity_relation": ["=", ">"],
            "source": ["SPD", "SPD"],
            "spd_target_protein_class": ["class_a", "class_b"],
            "spd_drugcentral_struct_id": ["1", "2"],
            "spd_inchikey": ["AAAA", "BBBB"],
            "_source_row_id": ["101", "102"],
            "_ac50_uM": [0.5, 20.0],
            "_free_cmax_uM": [0.1, 1.0],
        }
    )


def test_boundary_cases_are_independently_expected_and_pass() -> None:
    cases = _boundary_cases()

    assert cases.columns.tolist() == BOUNDARY_COLUMNS
    assert bool(cases["passed"].all())
    assert (
        cases.loc[
            cases["case_id"].eq("binding_open_inactive"), "actual_exposure_label"
        ].iloc[0]
        == 0
    )
    assert pd.isna(
        cases.loc[
            cases["case_id"].eq("binding_closed_inactive"),
            "actual_exposure_label",
        ].iloc[0]
    )


def test_relation_summary_keeps_zero_categories() -> None:
    frame = pd.DataFrame({"relation": ["=", ">"], "value": [1.0, 2.0]})

    summary = pd.DataFrame(
        _relation_summary(
            frame,
            stage="synthetic",
            relation_col="relation",
            activity_col="value",
        ),
        columns=RELATION_COLUMNS,
    )

    counts = summary.groupby("relation_category")["row_count"].sum().to_dict()
    assert counts["="] == 1
    assert counts[">"] == 1
    assert counts["<"] == 0
    assert counts["<="] == 0
    assert counts[">="] == 0
    assert counts["blank"] == 0
    assert counts["null"] == 0
    assert counts["unexpected"] == 0


def test_three_way_pair_comparison_matches_synthetic_materialization(
    tmp_path: Path,
) -> None:
    observations = _synthetic_observations()
    materialized = aggregate_spd_assays(observations)
    pair_path = tmp_path / "pairs.csv"
    materialized.to_csv(pair_path, index=False)
    selection_path = tmp_path / "selection.csv"
    pd.DataFrame(
        {
            "code_state": ["committed_head", "committed_head"],
            "selection_stage": [
                "raw_spd_pair_aggregation",
                "raw_spd_pair_aggregation",
            ],
            "drug_id": ["drug_a", "drug_b"],
            "target_id": ["TARGET_A", "TARGET_B"],
            "current_row_id": ["101", "102"],
            "current_binding_label": ["positive", "negative"],
            "current_exposure_label": ["positive", "negative"],
            "current_exposure_status": [
                "labeled_relevant",
                "labeled_censored_not_relevant",
            ],
        }
    ).to_csv(selection_path, index=False)

    pair, summary = _compare_pair_outputs(observations, pair_path, selection_path)

    assert len(pair) == 2
    assert summary["materialized_audit_vs_legacy_binding_mismatches"] == 0
    assert summary["materialized_panel_vs_legacy_exposure_mismatches"] == 0
    assert summary["legacy_vs_canonical_binding_mismatches"] == 0
    assert summary["legacy_vs_canonical_exposure_mismatches"] == 0
    assert summary["materialized_vs_replayed_selected_row_id_mismatches"] == 0


def test_behavior_schema_is_stable() -> None:
    assert BEHAVIOR_COLUMNS[:7] == [
        "source_row_id",
        "source_excel_row",
        "drug_identity",
        "target_identity",
        "raw_relation",
        "raw_activity_value",
        "activity_unit",
    ]
    assert BEHAVIOR_COLUMNS[14] == "production_selected"
    assert "binding_changed" in BEHAVIOR_COLUMNS
    assert "exposure_changed" in BEHAVIOR_COLUMNS


def test_materialized_expert_comparison_separates_legacy_and_canonical_delta(
    tmp_path: Path,
) -> None:
    path = tmp_path / "expert.csv"
    pd.DataFrame(
        {
            "drug_id": ["drug_a", "drug_b"],
            "target_id": ["TARGET_A", "TARGET_B"],
            "pdb_id": ["1ABC", "2DEF"],
            "spd_ac50_uM": [0.5, 10.0],
            "spd_activity_relation": ["=", ">"],
            "free_cmax_um": [0.1, 1.0],
            "spd_binding_label": [1, 0],
            "spd_exposure_label": [1, pd.NA],
        }
    ).to_csv(path, index=False)

    comparison = _compare_materialized_expert(path)

    assert comparison["materialized_vs_legacy_binding_mismatches"] == 0
    assert comparison["materialized_vs_legacy_exposure_mismatches"] == 0
    assert comparison["materialized_vs_canonical_binding_mismatches"] == 0
    assert comparison["materialized_vs_canonical_exposure_mismatches"] == 1
    assert comparison["changed_rows"][0]["pdb_id"] == "2DEF"


def _one_observation(
    *,
    relation: object,
    ac50_nm: object,
    free_cmax_nm: object = 1_000.0,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "drug_id": ["drug"],
            "target_id": ["TARGET"],
            "assay_id": ["assay"],
            "assay_name": ["assay"],
            "assay_type": ["biochemical"],
            "ac50_nM": [ac50_nm],
            "free_cmax_nM": [free_cmax_nm],
            "total_cmax_nM": [free_cmax_nm],
            "activity_relation": [relation],
            "source": ["SPD"],
        }
    )


def test_production_open_lower_bound_at_margin_ten_is_negative() -> None:
    pair = aggregate_spd_assays(_one_observation(relation=">", ac50_nm=10_000.0)).iloc[
        0
    ]

    assert pair["spd_label_status"] == "labeled_censored_not_relevant"
    assert pair["spd_missing_reason"] == ""
    assert bool(pair["spd_exposure_relevant"]) is False
    assert pd.isna(pair["spd_exposure_weak"])
    assert pd.isna(pair["spd_exposure_unlikely"])


def test_production_closed_lower_bound_at_margin_ten_remains_unknown() -> None:
    pair = aggregate_spd_assays(_one_observation(relation=">=", ac50_nm=10_000.0)).iloc[
        0
    ]

    assert pair["spd_label_status"] == (
        "unknown_censored_ac50_gt_crosses_relevant_threshold"
    )
    assert pd.isna(pair["spd_exposure_relevant"])


def test_production_censored_weak_boundary_statuses() -> None:
    open_pair = aggregate_spd_assays(
        _one_observation(relation=">", ac50_nm=100_000.0)
    ).iloc[0]
    closed_pair = aggregate_spd_assays(
        _one_observation(relation=">=", ac50_nm=100_000.0)
    ).iloc[0]

    assert open_pair["spd_label_status"] == "labeled_censored_unlikely"
    assert bool(open_pair["spd_exposure_relevant"]) is False
    assert pd.isna(open_pair["spd_exposure_weak"])
    assert bool(open_pair["spd_exposure_unlikely"]) is True
    assert closed_pair["spd_label_status"] == "labeled_censored_not_relevant"
    assert bool(closed_pair["spd_exposure_relevant"]) is False
    assert pd.isna(closed_pair["spd_exposure_weak"])
    assert pd.isna(closed_pair["spd_exposure_unlikely"])


def test_production_decimal_weak_boundary_statuses_are_exact() -> None:
    open_pair = aggregate_spd_assays(
        _one_observation(relation=">", ac50_nm=7.0, free_cmax_nm=0.07)
    ).iloc[0]
    closed_pair = aggregate_spd_assays(
        _one_observation(relation=">=", ac50_nm=7.0, free_cmax_nm=0.07)
    ).iloc[0]
    exact_pair = aggregate_spd_assays(
        _one_observation(relation="=", ac50_nm=7.0, free_cmax_nm=0.07)
    ).iloc[0]

    assert open_pair["spd_label_status"] == "labeled_censored_unlikely"
    assert closed_pair["spd_label_status"] == "labeled_censored_not_relevant"
    assert exact_pair["spd_label_status"] == "labeled_weak"


def test_materialized_panel_definitive_open_interval_is_negative(
    tmp_path: Path,
) -> None:
    panel_path = tmp_path / "panel.csv"
    pd.DataFrame(
        {
            "drug_id": ["drug"],
            "target_id": ["TARGET"],
            "ac50_nM": [0.7],
            "free_cmax_nM": [0.07],
            "spd_activity_relation": [">"],
            "spd_label_status": [""],
        }
    ).to_csv(panel_path, index=False)

    panel = _read_spd_panel(panel_path)

    assert panel.loc[0, "spd_exposure_label"] == 0


@pytest.mark.parametrize(
    ("relation", "ac50_nm"),
    [
        (">=", 0.7),
        (">", 0.5),
    ],
)
def test_materialized_panel_canonical_unknown_does_not_fall_back_to_margin(
    tmp_path: Path,
    relation: str,
    ac50_nm: float,
) -> None:
    panel_path = tmp_path / f"panel_{relation.replace('>', 'gt')}.csv"
    pd.DataFrame(
        {
            "drug_id": ["drug"],
            "target_id": ["TARGET"],
            "ac50_nM": [ac50_nm],
            "free_cmax_nM": [0.07],
            "exposure_margin": [5.0],
            "spd_activity_relation": [relation],
            "spd_label_status": [""],
        }
    ).to_csv(panel_path, index=False)

    panel = _read_spd_panel(panel_path)

    assert pd.isna(panel.loc[0, "spd_exposure_label"])


@pytest.mark.parametrize("relation", ["", "approximately"])
def test_materialized_panel_invalid_relation_cannot_use_scalar_margin(
    tmp_path: Path,
    relation: str,
) -> None:
    panel_path = tmp_path / "invalid_relation_panel.csv"
    pd.DataFrame(
        {
            "drug_id": ["drug"],
            "target_id": ["TARGET"],
            "ac50_nM": [0.7],
            "free_cmax_nM": [0.07],
            "exposure_margin": [5.0],
            "spd_activity_relation": [relation],
            "spd_label_status": [""],
        }
    ).to_csv(panel_path, index=False)

    panel = _read_spd_panel(panel_path)

    assert pd.isna(panel.loc[0, "spd_exposure_label"])


@pytest.mark.parametrize(
    ("status", "relation", "free_cmax_nm", "expected"),
    [
        ("", "=", None, None),
        ("", "=", 0.0, None),
        ("", "=", -0.07, None),
        ("", "=", 0.07, 1),
        ("labeled_relevant", "", 0.07, 1),
        ("labeled_censored_not_relevant", "", 0.07, 0),
    ],
)
def test_materialized_panel_fail_closed_and_direct_status_precedence(
    tmp_path: Path,
    status: str,
    relation: str,
    free_cmax_nm: float | None,
    expected: int | None,
) -> None:
    panel_path = tmp_path / "panel_precedence.csv"
    pd.DataFrame(
        {
            "drug_id": ["drug"],
            "target_id": ["TARGET"],
            "ac50_nM": [0.7],
            "free_cmax_nM": [free_cmax_nm],
            "exposure_margin": [5.0],
            "spd_activity_relation": [relation],
            "spd_label_status": [status],
        }
    ).to_csv(panel_path, index=False)

    panel = _read_spd_panel(panel_path)
    label = panel.loc[0, "spd_exposure_label"]

    if expected is None:
        assert pd.isna(label)
    else:
        assert label == expected


def test_binding_adapter_preserves_direct_precedence_contract() -> None:
    source = pd.DataFrame(
        {
            "spd_binding_label": [0, 1, 1, 0],
            "spd_ac50_uM": [0.5, 5.0, 0.5, 0.5],
            "spd_activity_relation": ["=", "=", "", "\u2264"],
        }
    )

    actual = _binding_label(source, active_um=1.0, inactive_um=10.0)

    assert actual.iloc[0] == 1
    assert pd.isna(actual.iloc[1])
    assert actual.iloc[2] == 1
    assert actual.iloc[3] == 0
    assert BINDING_LABEL_POLICY_VERSION == "spd_binding_interval_aware_v2"


def test_label_attachment_preserves_two_same_target_pdb_observations() -> None:
    source = pd.DataFrame(
        {
            "drug_id": ["drug", "drug"],
            "target_id": ["DRD1", "DRD1"],
            "pdb_id": ["9I52", "9LLG"],
            "structure_quality": [0.25, 0.95],
            "spd_ac50_uM": [0.5, 0.5],
            "spd_activity_relation": ["=", "="],
        }
    )

    source["spd_binding_label"] = _binding_label(
        source, active_um=1.0, inactive_um=10.0
    )

    assert len(source) == 2
    assert source["pdb_id"].tolist() == ["9I52", "9LLG"]
    assert source["structure_quality"].tolist() == [0.25, 0.95]
    assert source["spd_binding_label"].tolist() == [1, 1]


def test_boundary_verifier_uses_decimal_and_explicit_ppb_bounds() -> None:
    assert _decimal("30") / _decimal("3") == _decimal("10")
    assert _decimal("1") / _decimal("0.1") == _decimal("10")
    assert _qualifier_is_bound("=") is False
    assert _qualifier_is_bound(">=") is True


def test_production_status_projection_preserves_weak_and_censored_negatives() -> None:
    assert _production_exposure_label("labeled_relevant") == "positive"
    assert _production_exposure_label("labeled_weak") == "negative"
    assert _production_exposure_label("labeled_censored_not_relevant") == "negative"
    assert (
        _production_exposure_label(
            "unknown_censored_ac50_gt_crosses_relevant_threshold"
        )
        == "unknown"
    )


def test_boundary_provenance_schema_keeps_required_raw_and_pk_fields() -> None:
    assert {
        "source_row_id",
        "activity_relation_xml_value",
        "activity_ac50_xml_value",
        "cmax_source",
        "ppb_source",
        "ppb_qualifier",
        "free_cmax_xml_value",
        "lower_margin_decimal",
        "stored_value_contract_verdict",
        "biological_provenance_verdict",
    }.issubset(EVIDENCE_COLUMNS)
