from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.spd import SPD_LABEL_POLICY_VERSION
from analysis.ml.spd_four_expert_tables import build_spd_four_expert_tables
from analysis.ml.spd_label_enrichment import enrich_spd_labels_for_run_master


def _write_synthetic_inputs(tmp_path: Path) -> tuple[pd.DataFrame, Path, dict[str, Path]]:
    source = pd.DataFrame(
        {
            "row_token": ["nine", "six", "eight", "other"],
            "pdb_id": ["9I52", "6HUJ", "8HCQ", "1ABC"],
            "ligand_base": ["drug_a", "drug_b", "drug_c", "drug_d"],
            "drug_id": ["drug_a", "drug_b", "drug_c", "drug_d"],
            "target_id": ["ADRB2", "GABRA1", "EDNRA", "OTHER"],
            "target_gene": ["ADRB2", "GABRA1", "EDNRA", "OTHER"],
            "target_uniprot": ["P07550", "P14867", "P25101", "P00000"],
            "protein_class": ["legacy_beta", "legacy_gaba", "legacy_ednr", "other"],
            "target_family": [
                "legacy_beta_family",
                "gaba_family",
                "ednr_family",
                "other_family",
            ],
            "spd_ac50_uM": [20.0, 0.4, 20.0, 0.2],
            "spd_activity_relation": ["=", "=", "=", "="],
            "spd_binding_label": [0, 1, 0, 1],
            "spd_exposure_label": [0, 1, 0, 1],
            "exposure_margin": [200.0, 4.0, 200.0, 2.0],
            "tissue_site_label": [1, 1, 0, 1],
            "mechanism_ml_label": [1, 0, 1, 1],
            "spd_label_policy_version": [SPD_LABEL_POLICY_VERSION] * 4,
            # Historical model-ready inputs already contain projections from
            # earlier ligand/target-map joins. Strict refresh must replace
            # these for 9I52 without creating duplicate column names.
            "_mapped_drug_id": ["drug_a", "drug_b", "drug_c", "drug_d"],
            "_mapped_target_id": ["ADRB2", "GABRA1", "EDNRA", "OTHER"],
            "generic_name_fda_map": [
                "drug_a",
                "drug_b",
                "drug_c",
                "drug_d",
            ],
            "target_gene_target_map": [
                "ADRB2",
                "GABRA1",
                "EDNRA",
                "OTHER",
            ],
            "protein_class_target_map": [
                "legacy_beta",
                "legacy_gaba",
                "legacy_ednr",
                "other",
            ],
        }
    )
    source_path = tmp_path / "source.csv"
    source.to_csv(source_path, index=False)
    pd.DataFrame(
        {
            "ligand_base": ["drug_a", "drug_b", "drug_c", "drug_d"],
            "generic_name": ["drug_a", "drug_b", "drug_c", "drug_d"],
        }
    ).to_csv(tmp_path / "ligands.csv", index=False)
    pd.DataFrame(
        {
            "pdb_id": ["9I52", "6HUJ", "8HCQ", "1ABC"],
            "gene": ["ADRB2", "GABRA1", "EDNRA", "OTHER"],
        }
    ).to_csv(tmp_path / "targets.csv", index=False)
    pd.DataFrame(
        {
            "gene": ["ADRB2", "DRD1", "GABRA1", "EDNRA", "OTHER"],
            "protein_class": ["beta", "dopamine", "gaba", "ednr", "other"],
        }
    ).to_csv(tmp_path / "metadata.csv", index=False)
    pd.DataFrame(
        {
            "drug_id": ["drug_a"],
            "target_id": ["DRD1"],
            "assay_id": ["a1"],
            "assay_name": ["synthetic"],
            "ac50_nM": [500.0],
            "free_cmax_nM": [100.0],
            "total_cmax_nM": [200.0],
            "source": ["synthetic"],
            "assay_count": [1],
            "source_assay_ids": ["a1"],
            "exposure_margin": [5.0],
            "spd_label_status": ["labeled_relevant"],
            "spd_exposure_relevant": [True],
            "spd_activity_relation": ["="],
            "spd_label_policy_version": [SPD_LABEL_POLICY_VERSION],
        }
    ).to_csv(tmp_path / "panel.csv", index=False)
    paths = {
        "spd_panel_path": tmp_path / "panel.csv",
        "target_map_path": tmp_path / "targets.csv",
        "ligand_map_path": tmp_path / "ligands.csv",
        "target_metadata_path": tmp_path / "metadata.csv",
    }
    return source, source_path, paths


def test_receptor_mapping_modes_preserve_or_selectively_rebuild(tmp_path: Path) -> None:
    source, _, paths = _write_synthetic_inputs(tmp_path)

    legacy, legacy_summary = enrich_spd_labels_for_run_master(
        source.copy(),
        spd_panel_path=paths["spd_panel_path"],
        target_map_path=paths["target_map_path"],
        ligand_map_path=paths["ligand_map_path"],
        target_metadata_path=paths["target_metadata_path"],
    )
    pd.testing.assert_frame_equal(legacy, source)
    assert legacy_summary["status"] == "skipped"

    strict, strict_summary = enrich_spd_labels_for_run_master(
        source.copy(),
        receptor_mapping_mode="strict",
        spd_panel_path=paths["spd_panel_path"],
        target_map_path=paths["target_map_path"],
        ligand_map_path=paths["ligand_map_path"],
        target_metadata_path=paths["target_metadata_path"],
    )
    assert strict["row_token"].tolist() == source["row_token"].tolist()
    assert len(strict) == len(source)
    assert not strict.columns.duplicated().any()
    nine = strict.loc[strict["pdb_id"].eq("9I52")].iloc[0]
    assert (nine["target_id"], nine["target_gene"], nine["target_uniprot"]) == (
        "DRD1",
        "DRD1",
        "P21728",
    )
    assert float(nine["spd_ac50_uM"]) == 0.5
    assert int(nine["spd_exposure_label"]) == 1

    unchanged_columns = [
        "target_id",
        "target_gene",
        "target_uniprot",
        "spd_ac50_uM",
        "spd_binding_label",
        "spd_exposure_label",
        "exposure_margin",
        "tissue_site_label",
        "mechanism_ml_label",
    ]
    for pdb_id in ["6HUJ", "8HCQ", "1ABC"]:
        before = source.loc[source["pdb_id"].eq(pdb_id), unchanged_columns].reset_index(
            drop=True
        )
        after = strict.loc[strict["pdb_id"].eq(pdb_id), unchanged_columns].reset_index(
            drop=True
        )
        pd.testing.assert_frame_equal(after, before, check_dtype=False)
    unresolved = strict["pdb_id"].isin(["6HUJ", "8HCQ"])
    assert not strict.loc[
        unresolved, "spd_receptor_mapping_strict_eligible"
    ].astype(bool).any()
    assert strict_summary["strict_target_derived_rebuilt_rows"] == 1
    assert strict_summary["rows"] == len(source)


def test_strict_excludes_unresolved_receptors_from_all_model_ready_outputs(
    tmp_path: Path,
) -> None:
    source, source_path, paths = _write_synthetic_inputs(tmp_path)
    out_dir = tmp_path / "strict_out"

    manifest = build_spd_four_expert_tables(
        source_path,
        out_dir,
        receptor_mapping_mode="strict",
        spd_panel_path=paths["spd_panel_path"],
        target_map_path=paths["target_map_path"],
        ligand_map_path=paths["ligand_map_path"],
        target_metadata_path=paths["target_metadata_path"],
    )

    expected_pdb_order = source["pdb_id"].tolist()
    for expert in ["binding", "exposure", "tissue", "mechanism"]:
        raw = pd.read_csv(out_dir / f"ml_spd_{expert}_table.csv")
        assert raw["pdb_id"].tolist() == expected_pdb_order
        assert len(raw) == len(source)
        assert raw["pdb_id"].isin(["6HUJ", "8HCQ"]).sum() == 2
        unresolved_raw = raw.loc[raw["pdb_id"].isin(["6HUJ", "8HCQ"])]
        assert unresolved_raw["target_gene"].tolist() == ["GABRA1", "EDNRA"]
        if expert in {"tissue", "mechanism"}:
            label_col = (
                "tissue_site_label"
                if expert == "tissue"
                else "mechanism_ml_label"
            )
            assert pd.isna(
                raw.loc[raw["pdb_id"].eq("9I52"), label_col].iloc[0]
            )

        model_ready_path = (
            out_dir / "model_ready" / f"ml_spd_{expert}_model_ready.csv"
        )
        assert model_ready_path.exists()
        model_ready = pd.read_csv(model_ready_path)
        assert not model_ready["pdb_id"].isin(["6HUJ", "8HCQ"]).any()
        assert manifest["experts"][expert]["strict_ineligible_rows"] == 2
        assert (
            manifest["experts"][expert]["strict_ineligible_labelable_rows_excluded"]
            == 2
        )


def test_default_and_explicit_legacy_builds_are_equivalent(tmp_path: Path) -> None:
    _, source_path, paths = _write_synthetic_inputs(tmp_path)
    default_out = tmp_path / "legacy_default"
    explicit_out = tmp_path / "legacy_explicit"

    build_spd_four_expert_tables(source_path, default_out, **paths)
    build_spd_four_expert_tables(
        source_path,
        explicit_out,
        receptor_mapping_mode="legacy",
        **paths,
    )

    for expert in ["binding", "exposure", "tissue", "mechanism"]:
        default_raw = pd.read_csv(default_out / f"ml_spd_{expert}_table.csv")
        explicit_raw = pd.read_csv(explicit_out / f"ml_spd_{expert}_table.csv")
        pd.testing.assert_frame_equal(default_raw, explicit_raw)

        default_model_ready = (
            default_out / "model_ready" / f"ml_spd_{expert}_model_ready.csv"
        )
        explicit_model_ready = (
            explicit_out / "model_ready" / f"ml_spd_{expert}_model_ready.csv"
        )
        assert default_model_ready.exists() == explicit_model_ready.exists()
        if default_model_ready.exists():
            pd.testing.assert_frame_equal(
                pd.read_csv(default_model_ready),
                pd.read_csv(explicit_model_ready),
            )
