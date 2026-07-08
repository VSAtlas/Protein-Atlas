from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import docking.docking as docking
import docking.docking_receptor as docking_receptor
import docking.docking_runtime_process as docking_runtime_process
import docking.docking_controls as docking_controls
from protein_prep import automate_protein_prep
import docking.active_site_detection as protein_functions
from docking.fallback_recenter import RecenterParams


def test_control_redock_runs_after_ph_ensemble(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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
        cleaned.parent.mkdir(parents=True, exist_ok=True)
        cleaned.write_text(
            "ATOM      1  N   MET A   1      11.104  13.207  10.000  1.00 20.00           N\n",
            encoding="utf-8",
        )
        receptor.parent.mkdir(parents=True, exist_ok=True)
        receptor.write_text("REMARK dummy receptor\n", encoding="utf-8")
        return str(cleaned), str(receptor)

    def fake_build_ph_ensemble(**_kwargs):
        call_order.append("ph_ensemble")
        return str(tmp_path / "ph_manifest.json")

    def fake_select_center(*_args, **_kwargs):
        call_order.append("control_redock")
        return (1.0, 2.0, 3.0), (24.0, 24.0, 24.0)

    monkeypatch.setattr(docking_receptor, "prepare_receptor", fake_prepare_receptor)
    monkeypatch.setattr(
        protein_functions,
        "detect_active_site",
        lambda _pdb: ((0.0, 0.0, 0.0), (20.0, 20.0, 20.0), "test"),
    )
    import prep_ligands.prep_ligands_crystal

    monkeypatch.setattr(
        prep_ligands.prep_ligands_crystal,
        "prep_ligands_from_pdb",
        lambda **_kwargs: None,
    )
    monkeypatch.setattr(
        docking_controls, "extract_ligands_to_nolig", lambda *_a, **_k: (0, [])
    )
    monkeypatch.setattr(docking_controls, "build_control_lookup", lambda *_a, **_k: {})
    monkeypatch.setattr(
        docking_controls, "select_center_via_control_redock", fake_select_center
    )
    monkeypatch.setattr(
        docking_runtime_process, "_phase6_to8_ligands_and_docking", lambda *_a, **_k: None
    )
    monkeypatch.setattr(
        automate_protein_prep, "get_ion_probe_map", lambda *_a, **_k: {}
    )
    monkeypatch.setattr(
        automate_protein_prep,
        "_holo_restore_from_input_if_needed",
        lambda *_a, **_k: (0, 0, False),
    )
    monkeypatch.setattr(
        automate_protein_prep, "run_metal_site_audit", lambda *_a, **_k: None
    )

    from path_router import context_ph

    def fake_phase2_to4(*_args, **_kwargs):
        # (cleaned_pdb, receptor_pdbqt, pdb_audit, clean_audit, center, box_size, center_source, control_stems, control_lookup)
        cleaned, receptor = fake_prepare_receptor(_args[0], _args[1], _args[2])
        return (
            cleaned,
            receptor,
            {},
            {},
            (0.0, 0.0, 0.0),
            (20.0, 20.0, 20.0),
            "mock",
            [],
            {},
        )

    def fake_phase5(*_args, **_kwargs):
        call_order.append("ph_ensemble")
        # receptor_pdbqt, active_ph_label
        return "dummy_ph.pdbqt", "pH7_0"

    def fake_phase5b(*_args, **_kwargs):
        call_order.append("control_redock")
        # center, box_size, center_by_ph, box_by_ph, center_source_by_ph, control_stems, control_lookup
        return (
            (0, 0, 0),
            (20, 20, 20),
            {"pH7_0": (0, 0, 0)},
            {"pH7_0": (20, 20, 20)},
            {"pH7_0": "mock"},
            [],
            {},
        )

    monkeypatch.setattr(docking_runtime_process, "_phase5_ph_ensemble_global", fake_phase5)
    monkeypatch.setattr(docking_runtime_process, "_phase5b_controls_and_control_redock", fake_phase5b)
    monkeypatch.setattr(docking_runtime_process, "_phase2_to4_receptor_and_center", fake_phase2_to4)
    monkeypatch.setattr(
        context_ph, "select_ph_values_for_protonation", lambda *_a, **_k: [7.0]
    )

    docking.process_one_protein(
        cfg, f"{pdb_id}.pdb", stages=[], params=RecenterParams()
    )

    assert call_order == ["ph_ensemble", "control_redock"]
