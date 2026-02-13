from __future__ import annotations

import json
from pathlib import Path

import pocket_eval


def test_force_calibrator_invalidation(tmp_path, monkeypatch):
    pdb_id = "TEST"
    run_id = "force_run"
    docked_root = tmp_path / "docked"
    extracted_root = tmp_path / "extracted"
    cache_root = tmp_path / "cache"

    cfg = {
        "RUN_ID": run_id,
        "DOCKED_DIR": str(docked_root),
        "LIGAND_EXTRACTED_DIR": str(extracted_root),
        "CALIBRATOR_CACHE_DIR": cache_root,
        "FORCE_CALIBRATOR": True,
    }

    class DummyPaths:
        def __init__(self, root: Path):
            self._root = root

        def docked_pdb_root(self) -> Path:
            return self._root

    def fake_make_paths(cfg_obj, base_id, pdb_file):
        return DummyPaths(docked_root / run_id / base_id)

    monkeypatch.setattr(pocket_eval, "make_paths", fake_make_paths)

    dock_root = docked_root / run_id / pdb_id / "pocket_eval"
    dock_root.mkdir(parents=True, exist_ok=True)
    cache_path = dock_root / "calibrator_cache.json"
    cache_payload = {
        "pdb_id": pdb_id,
        "rows": [{"ligand_id": "old", "smiles": "OLD_SMILES", "label": "strong"}],
    }
    cache_path.write_text(json.dumps(cache_payload), encoding="utf-8")

    extracted_dir = extracted_root / f"{pdb_id}_calibrator"
    extracted_dir.mkdir(parents=True, exist_ok=True)
    (extracted_dir / "strong_binders.smi").write_text("OLD_SMILES old\n", encoding="utf-8")

    calls = {}

    def fake_run_calibrator(pdb_norm, **kwargs):
        calls["force_refresh"] = kwargs.get("force_refresh")
        calls["run_tag"] = kwargs.get("run_tag")
        return {}

    monkeypatch.setattr(pocket_eval, "run_calibrator_for_pdb", fake_run_calibrator)

    def fake_resolve_chain_uniprot_segments(pdb_norm, write_file=False):
        return {"A": {"uniprot": "PNEW", "unp_start": 1, "unp_end": 10}}

    def fake_select_primary_chain(chain_map):
        return "A", chain_map["A"]

    monkeypatch.setattr(
        pocket_eval,
        "resolve_chain_uniprot_segments",
        fake_resolve_chain_uniprot_segments,
    )
    monkeypatch.setattr(pocket_eval, "select_primary_chain", fake_select_primary_chain)

    def fake_fetch_chembl_labeled_smiles(
        uniprot_id,
        pdb_norm,
        cache_dir,
        timeout,
        retries,
        activity_types,
        chembl_max_phase,
        label_thresholds,
        **kwargs,
    ):
        return {"strong": ["NEW_SMILES"], "weak": [], "non": ["NON_BINDER"]}, {
            "cached": False,
            "counts": {},
        }

    monkeypatch.setattr(
        pocket_eval, "fetch_chembl_labeled_smiles", fake_fetch_chembl_labeled_smiles
    )

    rows = pocket_eval.get_calibration_set_for_pdb(pdb_id, cfg)
    smiles_set = {r["smiles"] for r in rows}

    assert "NEW_SMILES" in smiles_set
    assert "NON_BINDER" in smiles_set
    assert "OLD_SMILES" not in smiles_set
    assert calls.get("force_refresh") is True

    stored = json.loads(cache_path.read_text(encoding="utf-8"))
    stored_rows = stored.get("rows") or []
    assert all(row.get("smiles") != "OLD_SMILES" for row in stored_rows)
