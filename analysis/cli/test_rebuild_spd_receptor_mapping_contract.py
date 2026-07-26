from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from analysis.cli import rebuild_spd_receptor_mapping_contract as contract_cli


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_input(path: Path) -> None:
    pd.DataFrame(
        {
            "pdb_id": ["9I52", "6HUJ", "8HCQ", "1ABC"],
            "ligand_base": ["rdk_1", "rdk_2", "rdk_3", "rdk_4"],
            "drug_id": ["a", "b", "c", "d"],
            "target_id": ["OLD1", "OLD2", "OLD3", "KEEP"],
            "target_gene": ["OLD1", "OLD2", "OLD3", "KEEP"],
            "target_uniprot": ["U1", "U2", "U3", "U4"],
        }
    ).to_csv(path, index=False)


def _write_required_inputs(root: Path) -> dict[str, Path]:
    paths = {
        "spd_table": root / "source.csv",
        "spd_panel": root / "panel.csv",
        "target_map": root / "target_map.csv",
        "ligand_map": root / "ligand_map.csv",
        "contract": root / "contract.csv",
    }
    _write_input(paths["spd_table"])
    for name, path in paths.items():
        if name not in {"spd_table", "contract"}:
            pd.DataFrame({"value": [name]}).to_csv(path, index=False)
    review = root / "review.csv"
    pd.DataFrame(
        {
            "pdb_id": ["9I52"],
            "author_chain_id": ["R"],
            "candidate_primary_accession": ["P21728"],
            "candidate_gene_symbols": ["DRD1"],
            "chain_level_classification": ["probable_true_target_site_segment"],
            "docking_site_evidence_json": [
                json.dumps(
                    {
                        "segment_contacts": [
                            {
                                "chain_id": "R",
                                "reported_accession": "P21728",
                                "ligand_atoms_within_5A": 4,
                            }
                        ]
                    }
                )
            ],
        }
    ).to_csv(review, index=False)
    mappings = []
    for pdb_id in ("9I52", "6HUJ", "8HCQ"):
        local = root / f"{pdb_id}.pdb"
        prepared = root / f"{pdb_id}.pdbqt"
        local.write_text(f"HEADER {pdb_id}\n", encoding="utf-8")
        chain = "R" if pdb_id == "9I52" else "A"
        prepared.write_text(
            f"ATOM      1  C   GLY {chain}   1       0.000   0.000   0.000\n",
            encoding="utf-8",
        )
        mappings.append(
            {
                "pdb_id": pdb_id,
                "selected_chain": "R" if pdb_id == "9I52" else None,
                "site_chain": "R" if pdb_id == "9I52" else None,
                "revised_target_id": "DRD1" if pdb_id == "9I52" else None,
                "revised_target_gene": "DRD1" if pdb_id == "9I52" else None,
                "revised_target_uniprot": "P21728" if pdb_id == "9I52" else None,
                "evidence": [
                    {
                        "kind": "audit_mapping_review",
                        "path": review.name,
                        "sha256": _sha256(review),
                    },
                    {
                        "kind": "local_structure",
                        "path": local.name,
                        "sha256": _sha256(local),
                    },
                    {
                        "kind": "prepared_receptor",
                        "path": prepared.name,
                        "sha256": _sha256(prepared),
                    },
                ],
            }
        )
    paths["contract"].write_text(
        json.dumps({"mappings": mappings}),
        encoding="utf-8",
    )
    paths["evidence_root"] = root
    return paths


def _fake_builder(
    source_path: str | Path,
    out_dir: str | Path,
    *,
    receptor_mapping_mode: str,
    **_kwargs: object,
) -> dict[str, object]:
    source = pd.read_csv(source_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame = source.copy()
    frame["spd_binding_label"] = [1, 0, pd.NA, 1]
    frame["spd_exposure_label"] = [1, 0, 0, 1]
    frame["tissue_site_label"] = pd.NA
    frame["mechanism_ml_label"] = pd.NA
    frame["spd_receptor_mapping_legacy_target_id"] = source["target_id"]
    frame["spd_receptor_mapping_legacy_target_gene"] = source["target_gene"]
    frame["spd_receptor_mapping_legacy_target_uniprot"] = source["target_uniprot"]
    frame["spd_receptor_mapping_revised_target_id"] = source["target_id"]
    frame["spd_receptor_mapping_revised_target_gene"] = source["target_gene"]
    frame["spd_receptor_mapping_revised_target_uniprot"] = source["target_uniprot"]
    frame["spd_receptor_mapping_policy_version"] = f"{receptor_mapping_mode}_v1"
    frame["spd_receptor_mapping_status"] = "legacy_preserved"
    frame["spd_receptor_mapping_contract_source"] = "contract.csv"
    frame["spd_receptor_mapping_contract_sha256"] = "abc123"
    frame["spd_receptor_mapping_strict_eligible"] = True
    frame["spd_receptor_mapping_changed"] = False
    frame["spd_receptor_mapping_exclusion_reason"] = ""
    frame["spd_receptor_mapping_requires_target_derived_rebuild"] = False
    frame["spd_receptor_mapping_selected_chain"] = "A"
    frame["spd_receptor_mapping_site_chain"] = "A"
    if receptor_mapping_mode == "strict":
        affected = frame["pdb_id"].isin(contract_cli.EXPECTED_AFFECTED_PDB_SCOPE)
        frame.loc[affected, "target_id"] = ["NEW1", "NEW2", "NEW3"]
        frame.loc[affected, "target_gene"] = ["NEW1", "NEW2", "NEW3"]
        frame.loc[affected, "target_uniprot"] = ["N1", "N2", "N3"]
        frame.loc[affected, "spd_receptor_mapping_revised_target_id"] = [
            "NEW1",
            "NEW2",
            "NEW3",
        ]
        frame.loc[affected, "spd_receptor_mapping_revised_target_gene"] = [
            "NEW1",
            "NEW2",
            "NEW3",
        ]
        frame.loc[affected, "spd_receptor_mapping_revised_target_uniprot"] = [
            "N1",
            "N2",
            "N3",
        ]
        frame.loc[affected, "spd_receptor_mapping_status"] = "strict_contract_applied"
        frame.loc[affected, "spd_receptor_mapping_changed"] = True
        frame.loc[affected, "spd_receptor_mapping_requires_target_derived_rebuild"] = (
            True
        )
        unresolved = frame["pdb_id"].isin(["6HUJ", "8HCQ"])
        frame.loc[unresolved, "spd_receptor_mapping_strict_eligible"] = False
        frame.loc[
            unresolved,
            "spd_receptor_mapping_exclusion_reason",
        ] = "contract_ineligible"
        frame.loc[frame["pdb_id"].eq("9I52"), "spd_binding_label"] = 0
    for filename in contract_cli.EXPERT_TABLES.values():
        frame.to_csv(out / filename, index=False)
    model_ready = out / "model_ready"
    model_ready.mkdir(parents=True, exist_ok=True)
    label_columns = dict(
        zip(
            contract_cli.EXPERT_TABLES.values(),
            contract_cli.LABEL_COLUMNS,
            strict=True,
        )
    )
    for filename, label_column in label_columns.items():
        eligible = frame["spd_receptor_mapping_strict_eligible"].astype(bool)
        label = pd.to_numeric(frame[label_column], errors="coerce")
        selected = frame.loc[label.isin([0, 1]) & eligible].copy()
        if not selected.empty:
            selected.to_csv(
                model_ready / filename.replace("_table.csv", "_model_ready.csv"),
                index=False,
            )
    manifest = {"mode": receptor_mapping_mode, "rows": len(frame)}
    (out / "spd_four_expert_tables_manifest.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    return manifest


def test_comparison_writes_contract_artifacts_and_preserves_inputs(
    tmp_path: Path,
) -> None:
    paths = _write_required_inputs(tmp_path)
    before = {name: _sha256(path) for name, path in paths.items() if path.is_file()}

    manifest = contract_cli.run_receptor_mapping_comparison(
        **paths,
        out_dir=tmp_path / "out",
        expected_affected_counts={"9I52": 1, "6HUJ": 1, "8HCQ": 1},
        builder=_fake_builder,
    )

    assert manifest["status"] == "pass"
    assert all(isinstance(value, bool) for value in manifest["checks"].values())
    assert set(manifest["code_sha256"]) == {
        "orchestration_cli",
        "four_expert_builder",
        "label_enrichment",
        "receptor_mapping_module",
        "model_ready_deduplication",
        "mechanism_projection",
        "tissue_projection",
    }
    assert all(manifest["code_sha256"].values())
    assert manifest["affected_counts"] == {"6HUJ": 1, "8HCQ": 1, "9I52": 1}
    assert manifest["strict_eligibility"]["ineligible_rows"] == 2
    assert manifest["label_counts"]["legacy"]["spd_binding_label"] == {
        "positive": 2,
        "negative": 1,
        "excluded": 0,
        "unknown": 1,
    }
    assert manifest["label_counts"]["strict"]["spd_binding_label"]["positive"] == 1
    assert all(
        _sha256(path) == before[name] for name, path in paths.items() if path.is_file()
    )
    comparison = pd.read_csv(tmp_path / "out/receptor_mapping_before_after.csv")
    assert comparison.columns.tolist() == [
        "pdb_id",
        "legacy_target",
        "revised_status",
        "rows_affected",
        "binding_label_changes",
        "exposure_label_changes",
        "strict_eligibility",
    ]
    assert set(comparison["pdb_id"]) == {"9I52", "6HUJ", "8HCQ"}
    row_deltas = pd.read_csv(tmp_path / "out/receptor_mapping_row_deltas.csv")
    assert row_deltas["source_row_number"].tolist() == [0, 1, 2]
    assert manifest["model_ready_pdb_counts"]["strict"]["binding"] == {
        "1ABC": 1,
        "9I52": 1,
    }
    assert (tmp_path / "out/receptor_mapping_before_after.md").is_file()
    assert (tmp_path / "out/mapping_provenance_manifest.json").is_file()
    assert (tmp_path / "out/legacy/ml_spd_binding_table.csv").is_file()
    assert (tmp_path / "out/strict/ml_spd_binding_table.csv").is_file()


def test_unrelated_pdb_change_fails_closed_and_records_manifest(
    tmp_path: Path,
) -> None:
    paths = _write_required_inputs(tmp_path)

    def builder(
        source_path: str | Path,
        out_dir: str | Path,
        *,
        receptor_mapping_mode: str,
        **kwargs: object,
    ) -> dict[str, object]:
        manifest = _fake_builder(
            source_path,
            out_dir,
            receptor_mapping_mode=receptor_mapping_mode,
            **kwargs,
        )
        if receptor_mapping_mode == "strict":
            out = Path(out_dir)
            for filename in contract_cli.EXPERT_TABLES.values():
                frame = pd.read_csv(out / filename)
                frame.loc[frame["pdb_id"].eq("1ABC"), "target_id"] = "UNRELATED"
                frame.to_csv(out / filename, index=False)
        return manifest

    with pytest.raises(
        contract_cli.ContractViolation,
        match="unrelated_pdb_rows_unchanged",
    ):
        contract_cli.run_receptor_mapping_comparison(
            **paths,
            out_dir=tmp_path / "out",
            expected_affected_counts={"9I52": 1, "6HUJ": 1, "8HCQ": 1},
            builder=builder,
        )

    manifest = json.loads(
        (tmp_path / "out/mapping_provenance_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["status"] == "fail"
    assert manifest["checks"]["unrelated_pdb_rows_unchanged"] is False
    assert manifest["unrelated_changed_columns"] == ["target_id"]


def test_row_multiplication_fails_before_comparison(
    tmp_path: Path,
) -> None:
    paths = _write_required_inputs(tmp_path)

    def builder(
        source_path: str | Path,
        out_dir: str | Path,
        *,
        receptor_mapping_mode: str,
        **kwargs: object,
    ) -> dict[str, object]:
        manifest = _fake_builder(
            source_path,
            out_dir,
            receptor_mapping_mode=receptor_mapping_mode,
            **kwargs,
        )
        if receptor_mapping_mode == "strict":
            out = Path(out_dir)
            for filename in contract_cli.EXPERT_TABLES.values():
                frame = pd.read_csv(out / filename)
                pd.concat([frame, frame.iloc[[0]]], ignore_index=True).to_csv(
                    out / filename,
                    index=False,
                )
        return manifest

    with pytest.raises(contract_cli.ContractViolation, match="multiplied or dropped"):
        contract_cli.run_receptor_mapping_comparison(
            **paths,
            out_dir=tmp_path / "out",
            expected_affected_counts={"9I52": 1, "6HUJ": 1, "8HCQ": 1},
            builder=builder,
        )


def test_evidence_hash_mismatch_fails_before_builder(
    tmp_path: Path,
) -> None:
    paths = _write_required_inputs(tmp_path)
    (tmp_path / "9I52.pdbqt").write_text(
        "ATOM      1  C   GLY A   1       0.000   0.000   0.000\n",
        encoding="utf-8",
    )
    called = False

    def builder(*_args: object, **_kwargs: object) -> dict[str, object]:
        nonlocal called
        called = True
        return {}

    with pytest.raises(contract_cli.ContractViolation, match="SHA-256 mismatch"):
        contract_cli.run_receptor_mapping_comparison(
            **paths,
            out_dir=tmp_path / "out",
            expected_affected_counts={"9I52": 1, "6HUJ": 1, "8HCQ": 1},
            builder=builder,
        )
    assert called is False
    assert not (tmp_path / "out").exists()


def test_semantic_comparison_accepts_numeric_serialization_only() -> None:
    left = pd.Series([0.0, 1.0, pd.NA, "same"])
    right = pd.Series(["0.0", "1", pd.NA, "same"])
    assert contract_cli._semantic_equal(left, right).all()


def test_label_transition_counts_separate_binary_flips_from_unknowns() -> None:
    legacy = pd.DataFrame(
        {
            "spd_binding_label": [0, 1, 0, pd.NA, 1],
        }
    )
    strict = pd.DataFrame(
        {
            "spd_binding_label": [1, 0, pd.NA, 1, 1],
        }
    )

    transitions = contract_cli._label_transition_counts(
        legacy,
        strict,
        pd.Series(True, index=legacy.index),
    )["spd_binding_label"]

    assert transitions == {
        "categorical_changes": 4,
        "binary_flips": 2,
        "to_unknown": 1,
        "from_unknown": 1,
    }
