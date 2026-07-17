from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from analysis.ml.feature_sets import FEATURE_SETS
from analysis.ml.grouped_cv_stability import _metric_row, _summarize_fold_metrics
from analysis.ml.ligand_descriptors import (
    DESCRIPTOR_COLUMNS,
    DESCRIPTOR_GROUP_COLUMNS,
    _mapping_smiles,
    add_ligand_physchem_descriptors,
)
from analysis.ml.pair_interaction_features import (
    PairKey,
    PairProvenanceError,
    _require_unique_keys,
    _verified_hash,
    sha256_file,
)
from analysis.ml.pair_residual_binding import _outer_folds
from analysis.ml.pk_mechanistic_features import (
    PK_CONTEXT_FREE_CMAX_TARGET,
    SPD_FREE_CMAX_TARGET,
    _alignment,
)
from analysis.ml.spd_censor_relabel import relabel_spd_binding_censor_aware
from analysis.ml.spd_exposure_grouped_oof import _make_folds
from analysis.ml.spd_four_expert_tables import _binding_label


def test_descriptor_registry_is_complete_deterministic_and_no_qed_primary() -> None:
    expected_columns = [
        column
        for columns in DESCRIPTOR_GROUP_COLUMNS.values()
        for column in columns
    ]
    assert DESCRIPTOR_COLUMNS == expected_columns
    assert len(DESCRIPTOR_COLUMNS) == len(set(DESCRIPTOR_COLUMNS))
    assert list(DESCRIPTOR_GROUP_COLUMNS) == [
        "physchem",
        "functional_counts",
        "branch_counts",
        "topological_shape",
    ]

    mapped_smiles, _source = _mapping_smiles(
        pd.Series({"rdk_id": "rdk_000001", "drug_id": "collision"}),
        {
            "base:rdk_000001": ("CCO", "exact"),
            "name:collision": ("c1ccccc1", "name"),
        },
    )
    assert mapped_smiles == "CCO"

    source = pd.DataFrame(
        {
            "smiles": ["CCO", "CC(=O)Oc1ccccc1C(=O)O"],
            "rdkit_mol_wt": [123.0, pd.NA],
        }
    )
    described, summary = add_ligand_physchem_descriptors(source)
    assert described.loc[0, "rdkit_mol_wt"] == 123.0
    assert described[DESCRIPTOR_COLUMNS].notna().all().all()
    assert summary["ligand_descriptor_columns"] == DESCRIPTOR_COLUMNS
    assert list(summary["ligand_descriptor_groups"]) == list(
        DESCRIPTOR_GROUP_COLUMNS
    )

    for name in (
        "ligand_physchem_descriptors_no_qed",
        "ligand_full_interpretable_no_qed",
        "ligand_primary_functional_branch_no_qed",
        "spd_binding_pair_final_full_no_qed",
        "spd_binding_pair_final_full_structure_no_qed",
        "spd_binding_pair_final_primary_functional_branch_no_qed",
        "spd_exposure_descriptor_functional_branch_no_qed",
    ):
        assert "rdkit_qed" not in FEATURE_SETS[name]


def test_double_cold_outer_folds_exclude_both_axes() -> None:
    frame = pd.DataFrame(
        [
            {
                "drug_id": f"D{drug}",
                "target_id": f"T{target}",
                "label": (drug + target) % 2,
            }
            for drug in range(6)
            for target in range(6)
        ]
    )
    folds, row_cells = _outer_folds(
        frame,
        label_col="label",
        group_columns=("drug_id", "target_id"),
        n_splits=3,
        seed=17,
    )

    all_positions = set(range(len(frame)))
    assert len(row_cells) == len(frame)
    assert row_cells.nunique() == 9
    assert len(folds) == 9
    for fold in folds:
        train_positions = set(map(int, fold["train_pos"]))
        test_positions = set(map(int, fold["test_pos"]))
        embargo_positions = set(map(int, fold["embargo_pos"]))
        assert test_positions
        assert train_positions.isdisjoint(test_positions)
        assert train_positions.isdisjoint(embargo_positions)
        assert test_positions.isdisjoint(embargo_positions)
        assert train_positions | test_positions | embargo_positions == all_positions

        train = frame.iloc[sorted(train_positions)]
        test = frame.iloc[sorted(test_positions)]
        assert set(train["drug_id"]).isdisjoint(test["drug_id"])
        assert set(train["target_id"]).isdisjoint(test["target_id"])


def test_grouped_exposure_folds_are_deterministic_and_drug_disjoint() -> None:
    data = pd.DataFrame(
        [
            {"drug_id": f"D{drug}", "_observed_exposure_label": label}
            for drug in range(8)
            for label in (0, 1)
        ]
    )
    first, first_strategy, first_count = _make_folds(
        data,
        group_col="drug_id",
        n_splits=4,
        seed=23,
    )
    second, second_strategy, second_count = _make_folds(
        data,
        group_col="drug_id",
        n_splits=4,
        seed=23,
    )
    assert (first_strategy, first_count) == (second_strategy, second_count)
    assert [list(test) for _, test in first] == [list(test) for _, test in second]

    assigned: set[int] = set()
    for train_index, test_index in first:
        assert set(data.loc[train_index, "drug_id"]).isdisjoint(
            data.loc[test_index, "drug_id"]
        )
        assert assigned.isdisjoint(map(int, test_index))
        assigned.update(map(int, test_index))
    assert assigned == set(map(int, data.index))


def test_censor_aware_binding_thresholds_do_not_invent_labels() -> None:
    source = pd.DataFrame(
        {
            "spd_binding_label": [pd.NA] * 6,
            "spd_ac50_uM": [0.5, 1.0, 10.0, 20.0, 5.0, 12.0],
            "spd_activity_relation": ["=", "<=", ">=", "<=", ">", "="],
        }
    )
    labels = _binding_label(source, active_um=1.0, inactive_um=10.0)
    assert labels.astype("string").fillna("unknown").tolist() == [
        "1",
        "1",
        "0",
        "unknown",
        "unknown",
        "0",
    ]


def test_pair_hash_and_duplicate_guards_fail_closed(tmp_path: Path) -> None:
    key = PairKey(("drug_id", "target_id"), ("D1", "T1"))
    with pytest.raises(PairProvenanceError, match="duplicate canonical pair key"):
        _require_unique_keys([key, key])

    pose = tmp_path / "pose.pdbqt"
    pose.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    expected = sha256_file(pose)
    assert _verified_hash(pose, expected, {}, "pose") == expected
    with pytest.raises(PairProvenanceError, match="SHA-256 mismatch"):
        _verified_hash(pose, "0" * 64, {}, "pose")


def test_pk_endpoint_alignment_rejects_cross_scenario_inputs() -> None:
    frame = pd.DataFrame(
        {
            "pk_mech_dose_normalized_value": [10.0],
            "free_cmax_um": [0.25],
            "source_name": ["SPD"],
            "pk_context_id": ["spd-cmax-1"],
            "measurement_context": ["observed_cmax"],
            "dose_context_type": ["cmax_study_matched"],
            "training_allowed": [True],
            "pk_context_free_cmax_um": [0.30],
            "pk_context_source_name": ["PKDB"],
            "pk_context_context_id": ["pk-context-1"],
            "pk_context_study_id": ["study-1"],
            "pk_context_dose_context_type": ["cmax_study_matched"],
            "pk_context_training_allowed": [True],
        }
    )

    status, _reason, allowed = _alignment(
        frame,
        0,
        field="dose",
        output_column="pk_mech_dose_normalized_value",
        selected_source="spd_dose_mg",
        target_endpoint=SPD_FREE_CMAX_TARGET,
    )
    assert (status, allowed) == ("aligned_same_cmax_scenario", True)

    status, _reason, allowed = _alignment(
        frame,
        0,
        field="dose",
        output_column="pk_mech_dose_normalized_value",
        selected_source="pk_context_dose_normalized_value",
        target_endpoint=SPD_FREE_CMAX_TARGET,
    )
    assert status == "mismatched_external_context_for_spd_target"
    assert allowed is False

    status, _reason, allowed = _alignment(
        frame,
        0,
        field="dose",
        output_column="pk_mech_dose_normalized_value",
        selected_source="pk_context_dose_normalized_value",
        target_endpoint=PK_CONTEXT_FREE_CMAX_TARGET,
    )
    assert (status, allowed) == ("aligned_same_cmax_scenario", True)


def test_grouped_cv_summary_includes_brier_baseline_and_skill() -> None:
    prediction = pd.DataFrame(
        {"label": [0, 1, 0, 1], "ml_prediction_score": [0.1, 0.9, 0.2, 0.8]}
    )
    row = _metric_row(prediction, label_col="label")
    assert row["Brier_prevalence_baseline"] == pytest.approx(0.25)
    assert row["Brier_skill_score"] > 0
    row.update(
        {"model_type": "logistic", "evaluation": "drug", "group_col": "drug_id"}
    )
    summary = _summarize_fold_metrics(pd.DataFrame([row]))
    assert summary.loc[0, "Brier_prevalence_baseline_median"] == pytest.approx(
        0.25
    )
    assert summary.loc[0, "Brier_skill_score_median"] > 0

    one_class = _metric_row(prediction.iloc[[0, 2]], label_col="label")
    assert "Brier_prevalence_baseline" in one_class
    assert "Brier_skill_score" in one_class


def test_censor_relabel_accepts_missing_optional_combined_label(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.csv"
    output = tmp_path / "relabelled.csv"
    pd.DataFrame(
        {
            "spd_binding_label": [pd.NA, pd.NA],
            "spd_ac50_uM": [0.5, 20.0],
            "spd_activity_relation": ["<=", ">="],
        }
    ).to_csv(source, index=False)

    manifest = relabel_spd_binding_censor_aware(source, output)

    assert manifest["combined_activity_recomputed"] is False
    assert manifest["n_combined_labels_changed"] == 0
    assert output.is_file()
