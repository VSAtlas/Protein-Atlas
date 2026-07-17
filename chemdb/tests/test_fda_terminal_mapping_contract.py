from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from analysis.external.pilot_bioactivity_sources import build_pilot_pair_table
from analysis.ml.spd_label_enrichment import _read_fda_ligand_map
from analysis.reporting.fda_name_map import mapping_file_provenance_path
from analysis.reporting.ligand_side_effect_cache import build_ligand_query
from prep_ligands.fda_mapping_identity import (
    authoritative_drugcentral_id,
    clean_mapping_value,
)
from prep_ligands.fda_mapping_promotion import (
    _validate_mapping_dispositions,
    _validate_verifier_checks,
)


CANONICAL_MAPPING_SHA256 = (
    "4697ca38ed842cbc5d673b3e97c7c5fe91479eb3997091400a2a5d4d626ff483"
)
REQUIRED_VERIFIER_CHECKS = (
    "mapping_row_and_rdk_uniqueness",
    "prepared_path_sha_parity",
    "prepared_connectivity_and_score_reuse_partition",
    "deserpidine_regression_sentinel",
    "terminal_disposition_partition",
    "zero_nonterminal_usable_rows",
    "fda_rows_have_structure_and_approval_evidence",
    "fda_canonical_parent_readiness",
    "unsafe_parent_collapse_repair_quarantine",
    "named_library_manifest_counts_and_hashes",
    "fda_source_record_partition",
    "pubchem_manifest_batch_hashes_and_offline_coverage",
)


def _passing_verifier_checks(*, mapping_rows: int) -> list[dict[str, object]]:
    return [
        {
            "name": name,
            "passed": True,
            "metrics": (
                {"mapping_rows": mapping_rows}
                if name == "mapping_row_and_rdk_uniqueness"
                else {}
            ),
        }
        for name in REQUIRED_VERIFIER_CHECKS
    ]


def test_canonical_fda_mapping_bytes_and_shape_are_pinned() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    mapping = repo_root / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv"
    assert hashlib.sha256(mapping.read_bytes()).hexdigest() == CANONICAL_MAPPING_SHA256

    with mapping.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    assert len(rows) == 4_730
    assert len(fields) == 76
    assert len({row["rdk_id"] for row in rows}) == len(rows)
    assert len({row["mapping_row_number"] for row in rows}) == len(rows)


def test_mapping_values_and_authoritative_drugcentral_ids_fail_closed() -> None:
    assert clean_mapping_value(pd.NA) == ""
    assert (
        authoritative_drugcentral_id(
            {
                "identity_resolution_status": "terminal_structure_resolved",
                "resolved_drugcentral_id": "222",
                "drugcentral_id": "111",
            }
        )
        == "222"
    )
    assert (
        authoritative_drugcentral_id(
            {
                "identity_resolution_status": "terminal_structure_resolved",
                "resolved_drugcentral_id": pd.NA,
                "drugcentral_id": "111",
            }
        )
        == ""
    )
    assert authoritative_drugcentral_id({"drugcentral_id": "333"}) == "333"


def test_terminal_side_effect_query_suppresses_stale_auxiliary_ids(
    tmp_path: Path,
) -> None:
    terminal = {
        "identity_resolution_status": "terminal_structure_resolved",
        "resolved_preferred_name": "proguanil",
        "rxnorm_rxcui": "891788",
        "pubchem_unii_list": "SOI2LOH54Z",
        "generic_name": "zinc oxide",
    }
    query = build_ligand_query(
        tmp_path,
        {"ligand_base": "rdk_0000353", "ligand_display": "rdk_0000353"},
        mapping_rows={"rdk 0000353": terminal},
    )
    assert query.query_names == ("proguanil",)
    assert query.rxcuis == ()
    assert query.uniis == ()

    legacy = {
        "generic_name": "zinc oxide",
        "rxnorm_rxcui": "891788",
        "pubchem_unii_list": "SOI2LOH54Z",
    }
    legacy_query = build_ligand_query(
        tmp_path,
        {"ligand_base": "rdk_0000353", "ligand_display": "rdk_0000353"},
        mapping_rows={"rdk 0000353": legacy},
    )
    assert legacy_query.rxcuis == ("891788",)
    assert legacy_query.uniis == ("SOI2LOH54Z",)


def test_promotion_requires_exact_independent_check_contract() -> None:
    checks = _passing_verifier_checks(mapping_rows=4_730)
    _validate_verifier_checks(
        {"checks": list(reversed(checks))},
        expected_rows=4_730,
    )

    invalid_reports: list[list[dict[str, object]]] = [
        checks[:1],
        checks[:-1],
        [*checks[:-1], dict(checks[0])],
        [
            *checks[:-1],
            {"name": "forged_replacement", "passed": True, "metrics": {}},
        ],
        [
            {**checks[0], "name": f" {checks[0]['name']}"},
            *checks[1:],
        ],
        [
            {**checks[0], "passed": False},
            *checks[1:],
        ],
        [
            {
                **checks[0],
                "metrics": {"mapping_rows": 1},
            },
            *checks[1:],
        ],
    ]
    for invalid in invalid_reports:
        with pytest.raises(ValueError):
            _validate_verifier_checks({"checks": invalid}, expected_rows=4_730)

    with pytest.raises(ValueError, match="unsupported terminal dispositions"):
        _validate_mapping_dispositions(
            [{"terminal_disposition": "unreviewed_custom_value"}]
        )


def test_mapping_provenance_path_never_leaks_external_absolute_root(
    tmp_path: Path,
) -> None:
    repo_root = tmp_path / "repo"
    inside = repo_root / "chemdb" / "data" / "mapping.csv"
    outside = tmp_path / "private" / "mapping.csv"
    assert (
        mapping_file_provenance_path(inside, repo_root)
        == "chemdb/data/mapping.csv"
    )
    assert mapping_file_provenance_path(outside, repo_root) == "external/mapping.csv"
    assert mapping_file_provenance_path(None, repo_root) == ""


def test_spd_mapping_consumer_prefers_only_authorized_resolved_ids(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    pd.DataFrame(
        [
            {
                "path": "library/rdk_0000001.pdbqt",
                "generic_name": "resolved",
                "drugcentral_id": "111",
                "resolved_drugcentral_id": "222",
                "identity_resolution_status": "terminal_structure_resolved",
            },
            {
                "path": "library/rdk_0000002.pdbqt",
                "generic_name": "legacy",
                "drugcentral_id": "333",
                "resolved_drugcentral_id": "",
                "identity_resolution_status": "",
            },
            {
                "path": "library/rdk_0000003.pdbqt",
                "generic_name": "blank-resolved",
                "drugcentral_id": "444",
                "resolved_drugcentral_id": "",
                "identity_resolution_status": "terminal_structure_resolved",
            },
        ]
    ).to_csv(mapping, index=False)

    consumed = _read_fda_ligand_map(mapping).set_index("ligand_base")
    assert consumed.loc["rdk_0000001", "drugcentral_id"] == "222"
    assert consumed.loc["rdk_0000002", "drugcentral_id"] == "333"
    assert consumed.loc["rdk_0000003", "drugcentral_id"] == ""


def test_pilot_pair_table_retains_authoritative_drugcentral_id(
    tmp_path: Path,
) -> None:
    master = tmp_path / "master.csv"
    mapping = tmp_path / "mapping.csv"
    pdb_dir = tmp_path / "pdbs"
    output = tmp_path / "pilot.csv"
    pdb_dir.mkdir()
    pd.DataFrame(
        [
            {
                "pdb_id": "1ABC",
                "ligand_base": "rdk_0000001",
                "is_decoy": False,
                "z_selected": -1.0,
                "consensus_score": 0.8,
            }
        ]
    ).to_csv(master, index=False)
    pd.DataFrame(
        [
            {
                "ligand_base": "rdk_0000001",
                "display_name": "resolved",
                "generic_name": "legacy",
                "inchikey": "ABCDEFGHIJKLMN-AAAA",
                "smiles": "CCO",
                "drugcentral_id": "111",
                "resolved_drugcentral_id": "222",
                "identity_resolution_status": "terminal_structure_resolved",
            }
        ]
    ).to_csv(mapping, index=False)
    (pdb_dir / "1ABC.pdb").write_text(
        "COMPND    MOLECULE: TEST TARGET;\n"
        "SOURCE    GENE: TEST1;\n"
        "DBREF  1ABC A    1   10  UNP    P12345\n"
        "ATOM\n",
        encoding="utf-8",
    )

    result = build_pilot_pair_table(master, mapping, pdb_dir, output)

    assert result.loc[0, "drugcentral_id"] == "222"
    assert output.is_file()
