from __future__ import annotations

import json
from pathlib import Path

from calibrator import chemdbl_calibrator


def test_calibrator_site_evidence_written(tmp_path: Path, monkeypatch) -> None:
    pdb_id = "TEST"
    out_root = tmp_path / "extracted_ligands"
    cache_dir = tmp_path / "cache"
    chain_map = {
        "A": {"uniprot": "P12345"},
        "B": {"uniprot": "Q99999"},
    }

    def fake_resolve(pdb_id_arg, deepcoy_root=chemdbl_calibrator.DEFAULT_DEEPCOY_ROOT):
        return (
            "P12345",
            [],
            {
                "source": "dbref",
                "chain_segments": chain_map,
                "selected_chain": "A",
                "selected_segment": {"uniprot": "P12345"},
            },
        )

    def fake_fetch(
        uniprot_id,
        pdb_id_arg,
        cache_dir_arg,
        timeout,
        retries,
        activity_types,
        chembl_max_phase,
        label_thresholds,
        return_records: bool = False,
        **_kwargs,
    ):
        labels = {"strong": ["CCC"], "weak": [], "non": ["CC"]}
        records = [
            {
                "molecule_chembl_id": "CHEMBL1",
                "canonical_smiles": "CCC",
                "inchi_key": "KEY1",
                "label": "strong",
                "supporting_activity": {
                    "pchembl_value": 7.2,
                    "standard_value": "12",
                    "standard_units": "nM",
                    "standard_type": "Ki",
                    "standard_relation": "=",
                },
                "assay_chembl_id": "A1",
                "activity_chembl_id": "ACT1",
                "document_chembl_id": "DOC1",
            },
            {
                "molecule_chembl_id": "CHEMBL2",
                "canonical_smiles": "CC",
                "inchi_key": None,
                "label": "non",
                "supporting_activity": {
                    "pchembl_value": 4.2,
                    "standard_value": "15000",
                    "standard_units": "nM",
                    "standard_type": "Ki",
                    "standard_relation": "=",
                },
                "assay_chembl_id": "A2",
                "activity_chembl_id": "ACT2",
                "document_chembl_id": "DOC2",
            },
        ]
        meta = {
            "counts": {"targets": 1, "assays": 1, "activities": 2, "molecules": 2},
            "target_selection": {"selected_ids": ["T1"]},
            "records": records if return_records else [],
        }
        return labels, meta

    def fake_mechanism_fetch(molecule_chembl_id, target_chembl_id):
        if molecule_chembl_id == "CHEMBL1":
            return [{"site_id": "SITE1"}]
        return []

    def fake_binding_site_fetch(site_id):
        return {
            "site_id": site_id,
            "site_name": "Active site",
            "comment": "Test site",
            "site_components": [
                {
                    "component_id": "COMP1",
                    "component_type": "protein",
                    "accession": "P12345",
                }
            ],
        }

    monkeypatch.setattr(chemdbl_calibrator, "resolve_target", fake_resolve)
    monkeypatch.setattr(
        chemdbl_calibrator,
        "_load_chain_uniprot_map",
        lambda _pdb_id, fallback=None: chain_map,
    )

    meta = chemdbl_calibrator.run_calibrator_for_pdb(
        pdb_id,
        out_root=out_root,
        cache_dir=cache_dir,
        timeout=1,
        retries=0,
        activity_types=["Ki"],
        chembl_max_phase=4,
        label_thresholds=chemdbl_calibrator.DEFAULT_LABEL_THRESHOLDS,
        log_dir=tmp_path / "logs",
        run_tag="pytest_run",
        fetch_fn=fake_fetch,
        mechanism_fetch_fn=fake_mechanism_fetch,
        binding_site_fetch_fn=fake_binding_site_fetch,
    )

    assert meta.get("status") == "ok"

    out_dir = out_root / f"{pdb_id}_calibrator"
    ligands_path = out_dir / "calibrator_ligands.jsonl"
    evidence_path = out_dir / "calibrator_site_evidence.jsonl"
    assert ligands_path.is_file()
    assert evidence_path.is_file()

    ligand_rows = [
        json.loads(line)
        for line in ligands_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(ligand_rows) == 2
    required_keys = {
        "pdb_id",
        "selected_chain",
        "selected_uniprot",
        "ligand_id",
        "molecule_chembl_id",
        "canonical_smiles",
        "inchi_key",
        "label",
        "supporting_activity",
        "assay_chembl_id",
        "activity_chembl_id",
        "document_chembl_id",
        "chembl_max_phase",
        "activity_types",
    }
    for row in ligand_rows:
        assert required_keys.issubset(row.keys())

    evidence_rows = [
        json.loads(line)
        for line in evidence_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(evidence_rows) == 2
    chembl1 = next(
        row for row in evidence_rows if row.get("molecule_chembl_id") == "CHEMBL1"
    )
    assert chembl1["evidence_chains"] == ["A"]
    assert chembl1["evidence_strength"] == "site_components_mapped"
    chembl2 = next(
        row for row in evidence_rows if row.get("molecule_chembl_id") == "CHEMBL2"
    )
    assert chembl2["evidence_strength"] == "none"
