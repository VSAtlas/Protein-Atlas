from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import docking
from fallback_recenter import RecenterParams


def test_control_redock_runs_after_ph_ensemble(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    call_order: list[str] = []

    monkeypatch.delenv("APO_HOLO_VARIANT", raising=False)
    monkeypatch.delenv("APO_HOLO_MODE", raising=False)

    root = tmp_path / "root"
    input_dir = root / "input_pdbs"
    output_dir = root / "processed_pdbs"
    docked_dir = root / "docked"
    ligands_dir = root / "prepped_ligands"
    for path in (input_dir, output_dir, docked_dir, ligands_dir):
        path.mkdir(parents=True, exist_ok=True)

    pdb_id = "TEST"
    pdb_path = input_dir / f"{pdb_id}.pdb"
    pdb_path.write_text(
        "ATOM      1  N   MET A   1      11.104  13.207  10.000  1.00 20.00           N\n",
        encoding="utf-8",
    )

    cfg = {
        "OVERALL_DIR": str(root),
        "INPUT_DIR": str(input_dir),
        "OUTPUT_DIR": str(output_dir),
        "DOCKED_DIR": str(docked_dir),
        "OUTPUT_LIGANDS_DIR": str(ligands_dir),
        "PH_ENSEMBLE": True,
        "RUN_ID": "",
    }

    def fake_prepare_receptor(cfg, paths, logger):
        cleaned = paths.receptor_cleaned_pdb(None)
        receptor = paths.receptor_pdbqt(None, None)
        cleaned.write_text(
            "ATOM      1  N   MET A   1      11.104  13.207  10.000  1.00 20.00           N\n",
            encoding="utf-8",
        )
        receptor.write_text("REMARK dummy receptor\n", encoding="utf-8")
        return str(cleaned), str(receptor)

    def fake_build_ph_ensemble(**_kwargs):
        call_order.append("ph_ensemble")
        return str(tmp_path / "ph_manifest.json")

    def fake_select_center(*_args, **_kwargs):
        call_order.append("control_redock")
        return (1.0, 2.0, 3.0), (24.0, 24.0, 24.0)

    monkeypatch.setattr(docking, "prepare_receptor", fake_prepare_receptor)
    monkeypatch.setattr(
        docking,
        "detect_active_site",
        lambda _pdb: ((0.0, 0.0, 0.0), (20.0, 20.0, 20.0), "test"),
    )
    monkeypatch.setattr(docking, "prep_ligands_from_pdb", lambda **_kwargs: None)
    monkeypatch.setattr(docking, "extract_ligands_to_nolig", lambda *_a, **_k: (0, []))
    monkeypatch.setattr(docking, "build_control_lookup", lambda *_a, **_k: {})
    monkeypatch.setattr(docking, "select_center_via_control_redock", fake_select_center)
    monkeypatch.setattr(docking, "_phase6_to8_ligands_and_docking", lambda *_a, **_k: None)
    monkeypatch.setattr(docking.protein_prep, "get_ion_probe_map", lambda *_a, **_k: {})
    monkeypatch.setattr(
        docking.protein_prep,
        "_holo_restore_from_input_if_needed",
        lambda *_a, **_k: (0, 0, False),
    )
    monkeypatch.setattr(docking.protein_prep, "run_metal_site_audit", lambda *_a, **_k: None)

    import context_ph
    import ph_ensemble

    monkeypatch.setattr(context_ph, "select_ph_values_for_protonation", lambda *_a, **_k: [7.0])
    monkeypatch.setattr(ph_ensemble, "build_ph_ensemble", fake_build_ph_ensemble)

    docking.process_one_protein(cfg, f"{pdb_id}.pdb", stages=[], params=RecenterParams())

    assert call_order == ["ph_ensemble", "control_redock"]
