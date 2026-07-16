from __future__ import annotations

import pandas as pd

from analysis.external.pkdb_recovery import (
    build_canonical_pk_identity_table,
    build_drug_alias_map,
    collect_endpoint_inventory_rows,
    parse_study_sources,
    study_tsv_urls,
    to_pk_context,
)


def _model_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "drug_id": "RDK0001",
                "generic_name": "glimepiride",
                "display_name": "Glimepiride",
                "inchikey": "WGZKDVQMWBSCNE-UHFFFAOYSA-N",
                "rdkit_mol_wt": 490.617,
            }
        ]
    )


def _parse_outputs(
    model_table: pd.DataFrame,
    output_rows: list[dict[str, str]],
    *,
    reference: object = "12345",
    source_rights: pd.DataFrame | None = None,
) -> pd.DataFrame:
    intervention = pd.DataFrame(
        [
            {
                "name": "DOSE",
                "route": "oral",
                "application": "single dose",
                "measurement_type": "dosing",
                "form": "tablet",
                "substance": "glimepiride",
                "value": "3",
                "unit": "mg",
            }
        ]
    )
    groups = pd.DataFrame(
        [
            {
                "name": "adults",
                "parent": "",
                "measurement_type": "species",
                "choice": "human",
            }
        ]
    )
    output = pd.DataFrame(
        [
            {
                "group": "adults",
                "intervention": "DOSE",
                "substance": "glimepiride",
                "tissue": "plasma",
                "method": "HPLC",
                **row,
            }
            for row in output_rows
        ]
    )
    return parse_study_sources(
        study={
            "sid": "PKDB00001",
            "name": "Example",
            "date": "2026-01-01",
            "reference": reference,
            "licence": "open",
        },
        source_frames=[
            ("https://pk-db.com/media/data/interventions.tsv", "a" * 64, intervention),
            ("https://pk-db.com/media/data/groups.tsv", "b" * 64, groups),
            ("https://pk-db.com/media/data/outputs.tsv", "c" * 64, output),
        ],
        aliases=build_drug_alias_map(model_table),
        source_rights=source_rights,
    )


def test_canonical_identity_table_replaces_evidence_join_conflicts() -> None:
    model = pd.DataFrame(
        [
            {
                "drug_id": "example",
                "generic_name": "wrong one",
                "display_name": "wrong one",
                "ligand_base": "rdk_0000001",
                "inchikey": "WRONG-ONE-A",
                "rdkit_mol_wt": 123.4,
            },
            {
                "drug_id": "example",
                "generic_name": "wrong two",
                "display_name": "wrong two",
                "ligand_base": "rdk_0000002",
                "inchikey": "WRONG-TWO-B",
                "rdkit_mol_wt": 123.4,
            },
        ]
    )
    mapping = pd.DataFrame(
        [
            {
                "rdk_id": "rdk_0000001",
                "preferred_identity": "Example drug",
                "identity_parent_inchikey": "PARENT-INCHIKEY-A",
                "identity_structure_validated": True,
            },
            {
                "rdk_id": "rdk_0000002",
                "preferred_identity": "Example drug",
                "identity_parent_inchikey": "PARENT-INCHIKEY-A",
                "identity_structure_validated": True,
            },
            {
                "rdk_id": "rdk_9999999",
                "preferred_identity": "Out of universe",
                "identity_parent_inchikey": "OUTSIDE-INCHIKEY-A",
                "identity_structure_validated": True,
            },
        ]
    )

    canonical = build_canonical_pk_identity_table(model, mapping)
    assert len(canonical) == 2
    assert canonical["drug_id"].eq("PARENT-INCHIKEY-A").all()
    assert canonical["inchikey"].eq("PARENT-INCHIKEY-A").all()
    aliases = build_drug_alias_map(canonical)
    assert aliases["example drug"].inchikey == "PARENT-INCHIKEY-A"
    assert aliases["example drug"].rdkit_mol_wt == 123.4


def test_study_tsv_urls_only_accepts_official_source_files() -> None:
    detail = {
        "files": [
            "https://pk-db.com/media/data/.Study_Tab.tsv",
            "/media/data/.Study_Relative.tsv",
            "https://pk-db.com/media/data/Study.pdf",
            "https://example.org/not-accepted.tsv",
        ],
        "outputset": {
            "outputs": [{"source": "https://pk-db.com/media/data/.Study_Output.tsv"}]
        },
    }
    assert study_tsv_urls(detail) == [
        "https://pk-db.com/media/data/.Study_Output.tsv",
        "https://pk-db.com/media/data/.Study_Relative.tsv",
        "https://pk-db.com/media/data/.Study_Tab.tsv",
    ]


def test_parse_pkdb_sources_preserves_endpoint_and_intervention_context() -> None:
    aliases = build_drug_alias_map(_model_table())
    intervention = pd.DataFrame(
        [
            {
                "study": "Example",
                "name": "GLI3",
                "time": "0",
                "time_end": "NA",
                "time_unit": "hr",
                "route": "oral",
                "application": "single dose",
                "measurement_type": "dosing",
                "form": "tablet",
                "substance": "glimepiride",
                "value": "3",
                "unit": "mg",
            },
            {
                "study": "Example",
                "name": "GLI_IV",
                "time": "0",
                "time_end": "NA",
                "time_unit": "hr",
                "route": "intravenous",
                "application": "single dose",
                "measurement_type": "dosing",
                "form": "injection",
                "substance": "glimepiride",
                "value": "1",
                "unit": "mg",
            },
        ]
    )
    groups = pd.DataFrame(
        [
            {
                "study": "Example",
                "name": "adults",
                "parent": "",
                "measurement_type": "species",
                "choice": "homo sapiens",
            },
            {
                "study": "Example",
                "name": "adults",
                "parent": "",
                "measurement_type": "healthy",
                "choice": "Y",
            },
        ]
    )
    output = pd.DataFrame(
        [
            {
                "study": "Example",
                "group": "adults",
                "intervention": "GLI3",
                "measurement_type": "oral clearance",
                "substance": "glimepiride",
                "tissue": "plasma",
                "method": "HPLC",
                "mean": "41.6",
                "unit": "ml/min",
            },
            {
                "study": "Example",
                "group": "adults",
                "intervention": "GLI3",
                "measurement_type": "absolute bioavailability",
                "substance": "glimepiride",
                "tissue": "plasma",
                "method": "HPLC",
                "mean": "56",
                "unit": "%",
            },
            {
                "study": "Example",
                "group": "adults",
                "intervention": "GLI3",
                "measurement_type": "clearance",
                "substance": "creatinine",
                "mean": "77.7",
                "unit": "ml/min",
            },
        ]
    )
    study = {
        "sid": "PKDB00001",
        "name": "Example",
        "date": "2026-01-01",
        "reference": "12345",
        "licence": "open",
    }
    source_frames = [
        ("https://pk-db.com/media/data/interventions.tsv", "a" * 64, intervention),
        ("https://pk-db.com/media/data/groups.tsv", "b" * 64, groups),
        ("https://pk-db.com/media/data/outputs.tsv", "c" * 64, output),
    ]
    recovered = parse_study_sources(
        study=study,
        source_frames=source_frames,
        aliases=aliases,
    )
    inventory = collect_endpoint_inventory_rows(
        study=study,
        source_frames=source_frames,
        aliases=aliases,
    )

    assert len(inventory) == 2
    assert inventory["has_numeric_value"].all()
    assert inventory["species"].eq("Homo sapiens").all()
    assert len(recovered) == 2
    clearance = recovered.loc[
        recovered["measurement_context"].eq("apparent_oral_clearance")
    ].iloc[0]
    assert clearance["normalized_value"] == 2.496
    assert clearance["normalized_unit"] == "L/h"
    assert clearance["dose_value"] == 3.0
    assert clearance["route"] == "oral"
    assert bool(clearance["semantic_ready"])
    assert bool(clearance["endpoint_evidence_ready"])
    assert not bool(clearance["license_allows_ml_training"])
    assert not bool(clearance["training_allowed"])
    assert not bool(clearance["model_ready"])
    assert "source_level_training_rights_not_provided" in clearance["context_status"]

    bioavailability = recovered.loc[
        recovered["measurement_context"].eq("absolute_bioavailability")
    ].iloc[0]
    assert bioavailability["normalized_value"] == 56.0
    assert bioavailability["reference_route"] == ""
    assert not bool(bioavailability["endpoint_evidence_ready"])
    assert not bool(bioavailability["model_ready"])
    assert (
        "absolute_bioavailability_endpoint_linked_iv_comparator_not_verified"
        in bioavailability["context_status"]
    )

    context = to_pk_context(recovered)
    assert len(context) == 2
    assert context["clearance_value"].notna().sum() == 1
    assert context["bioavailability_value"].notna().sum() == 1
    assert not context["training_allowed"].fillna(False).astype(bool).any()
    assert context["missing_reason"].fillna("").ne("").all()


def test_cmax_and_surface_area_clearance_normalization_keep_provenance() -> None:
    recovered = _parse_outputs(
        _model_table(),
        [
            {"measurement_type": "Cmax", "mean": "250", "unit": "nmol/L"},
            {"measurement_type": "Cmax", "mean": "490.617", "unit": "ng/mL"},
            {
                "measurement_type": "oral clearance",
                "mean": "1.5",
                "unit": "L/min/m^2",
            },
        ],
        reference={
            "pmid": "14691614",
            "doi": "10.1007/s00228-003-0714-5",
            "sid": "14691614",
            "date": "2004-04-30",
        },
    )

    assert len(recovered) == 3
    nmolar = recovered.loc[recovered["raw_unit"].eq("nmol/L")].iloc[0]
    assert nmolar["normalized_value"] == 0.25
    assert nmolar["normalized_unit"] == "uM"
    assert nmolar["normalization_method"] == "molar_concentration_unit_scale"

    mass = recovered.loc[recovered["raw_unit"].eq("ng/mL")].iloc[0]
    assert mass["normalized_value"] == 1.0
    assert mass["normalized_unit"] == "uM"
    assert mass["rdkit_mol_wt"] == 490.617
    assert mass["normalization_method"] == "mass_concentration_via_rdkit_mol_wt"

    clearance = recovered.loc[recovered["raw_unit"].eq("L/min/m^2")].iloc[0]
    assert clearance["normalized_value"] == 90.0
    assert clearance["normalized_unit"] == "L/h/m2"
    assert clearance["normalization_method"] == "clearance_unit_scale"

    assert recovered["study_reference"].eq("PMID:14691614").all()
    assert recovered["study_reference_type"].eq("pmid").all()
    assert recovered["study_reference_id"].eq("14691614").all()
    assert recovered["publication_year"].eq(2004).all()
    assert not recovered["source_rights_allow_ml_training"].astype(bool).any()
    assert not recovered["training_allowed"].astype(bool).any()
    assert not recovered["model_ready"].astype(bool).any()


def test_structured_reference_matches_compact_source_rights_reference() -> None:
    recovered = _parse_outputs(
        _model_table(),
        [{"measurement_type": "Cmax", "mean": "250", "unit": "nmol/L"}],
        reference={
            "pmid": "14691614",
            "date": "2004-04-30",
        },
        source_rights=pd.DataFrame(
            [
                {
                    "study_sid": "PKDB00001",
                    "study_reference": "PMID:14691614",
                    "source_file_url": "",
                    "training_allowed": True,
                    "rights_reference": "https://example.org/license",
                    "rights_basis": "explicit source license",
                }
            ]
        ),
    )

    row = recovered.iloc[0]
    assert bool(row["source_rights_allow_ml_training"])
    assert bool(row["training_allowed"])
    assert bool(row["model_ready"])
    context = to_pk_context(recovered)
    assert context.iloc[0]["source_version"] == "2004"


def test_mass_cmax_requires_model_table_rdkit_mol_wt() -> None:
    table = _model_table().rename(columns={"rdkit_mol_wt": "molecular_weight"})
    recovered = _parse_outputs(
        table,
        [{"measurement_type": "Cmax", "mean": "490.617", "unit": "ng/mL"}],
    )
    row = recovered.iloc[0]

    assert pd.isna(row["rdkit_mol_wt"])
    assert pd.isna(row["normalized_value"])
    assert row["normalized_unit"] == ""
    assert row["normalization_method"] == "mass_concentration_requires_rdkit_mol_wt"
    assert "mass_concentration_requires_rdkit_mol_wt" in row["context_status"]
    assert not bool(row["semantic_ready"])
    assert not bool(row["model_ready"])


def test_alias_resolution_is_deterministic_for_duplicate_model_rows() -> None:
    table = pd.DataFrame(
        [
            {
                "drug_id": "RDK0001",
                "generic_name": "zeta",
                "display_name": "shared name",
                "inchikey": "AAAA-BBBB-C",
                "rdkit_mol_wt": 123.4,
            },
            {
                "drug_id": "RDK0001",
                "generic_name": "alpha",
                "display_name": "shared name",
                "inchikey": "AAAA-BBBB-C",
                "rdkit_mol_wt": 123.4,
            },
        ]
    )

    aliases = build_drug_alias_map(table)
    reversed_aliases = build_drug_alias_map(table.iloc[::-1].reset_index(drop=True))

    assert aliases == reversed_aliases
    assert aliases["shared name"].drug_name == "alpha"
    assert aliases["shared name"].inchikey == "AAAA-BBBB-C"
    assert aliases["shared name"].rdkit_mol_wt == 123.4


def test_same_drug_conflicting_inchikeys_reject_all_aliases() -> None:
    table = pd.DataFrame(
        [
            {
                "drug_id": "RDK0001",
                "generic_name": "alpha",
                "display_name": "shared name",
                "inchikey": "AAAA-BBBB-C",
            },
            {
                "drug_id": "RDK0001",
                "generic_name": "zeta",
                "display_name": "shared name",
                "inchikey": "XXXX-YYYY-Z",
            },
        ]
    )

    aliases = build_drug_alias_map(table)

    assert {"alpha", "zeta", "shared name"}.isdisjoint(aliases)


def test_ambiguous_aliases_are_not_used() -> None:
    table = pd.DataFrame(
        [
            {"drug_id": "RDK0001", "generic_name": "shared name"},
            {"drug_id": "RDK0002", "generic_name": "shared name"},
        ]
    )
    assert "shared name" not in build_drug_alias_map(table)
