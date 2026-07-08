from __future__ import annotations

import logging
from pathlib import Path

import docking.docking as docking_mod
import docking.docking_runtime_process as runtime_process
from path_router import make_paths


def test_process_one_protein_reuses_combo_prep_geometry(monkeypatch) -> None:
    calls: dict[str, int] = {"phase5b": 0, "control_multi": 0, "phase6": 0}

    dummy_paths = type("DummyPaths", (), {"pdb_id": "BNJS"})()
    dummy_logger = logging.getLogger("test.reuse")

    monkeypatch.setattr(
        runtime_process,
        "_phase0_setup_paths_and_logger",
        lambda cfg, pdb_file: (dummy_paths, "BNJS", dummy_logger),
    )
    monkeypatch.setattr(
        runtime_process,
        "_phase1_variant_and_ion_context",
        lambda cfg, paths, logger: (
            "HOLO",
            "HOLO",
            "HOLO",
            {},
            {},
            False,
            Path("receptor.pdbqt"),
        ),
    )
    monkeypatch.setattr(
        runtime_process,
        "_phase2_to4_receptor_and_center",
        lambda *args, **kwargs: (
            "cleaned.pdb",
            "receptor.pdbqt",
            {},
            {},
            None,
            None,
            "control",
            [],
            {},
        ),
    )
    monkeypatch.setattr(
        runtime_process,
        "_phase5_ph_ensemble_global",
        lambda *args, **kwargs: ("receptor.pdbqt", "pH7_0"),
    )

    def _phase5b_should_not_run(*args, **kwargs):
        calls["phase5b"] += 1
        raise AssertionError("phase5b should be skipped when reuse payload is present")

    monkeypatch.setattr(
        runtime_process, "_phase5b_controls_and_control_redock", _phase5b_should_not_run
    )

    def _control_multi_should_not_run(*args, **kwargs):
        calls["control_multi"] += 1
        raise AssertionError("control multi should be skipped for reuse chunks")

    monkeypatch.setattr(
        runtime_process, "run_control_docking_multi_engine", _control_multi_should_not_run
    )
    monkeypatch.setattr(runtime_process, "_collect_control_pdbqts", lambda *a, **k: [])

    def _phase6(*args, **kwargs):
        calls["phase6"] += 1

    monkeypatch.setattr(runtime_process, "_phase6_to8_ligands_and_docking", _phase6)

    cfg = {
        "RUN_ID": "",
        "_SKIP_CONTROL_REDOCK": True,
        "_SKIP_CONTROL_DOCKING": True,
        "_COMBO_PREP_CENTER_BY_PH": {"pH7_0": [1.0, 2.0, 3.0]},
        "_COMBO_PREP_BOX_BY_PH": {"pH7_0": [24.0, 24.0, 24.0]},
        "_COMBO_PREP_SOURCE_BY_PH": {"pH7_0": "control"},
        "_COMBO_PREP_CONTROL_STEMS": ["KEU_A510", "KEU_B610"],
        "_COMBO_PREP_CONTROL_LOOKUP": {"KEU_A510": "foo.sdf"},
    }

    runtime_process.process_one_protein(cfg, "BNJS.pdb", stages=[], params={})  # type: ignore[arg-type]

    assert calls["phase6"] == 1
    assert calls["phase5b"] == 0
    assert calls["control_multi"] == 0


def test_distributed_prep_reuse_skips_receptor_and_ph_rebuild(
    tmp_path, monkeypatch
) -> None:
    root = tmp_path / "reuse_direct"
    cfg = {
        "OVERALL_DIR": str(root),
        "INPUT_DIR": str(root / "input_pdbs"),
        "OUTPUT_DIR": str(root / "processed_pdbs"),
        "DOCKED_DIR": str(root / "docked"),
        "PREPPED_LIGANDS_DIR": str(root / "prepped_ligands"),
        "CONFIGS_DIR": str(root / "configs"),
        "RUN_ID": "",
        "PH_ENSEMBLE": True,
        "_DISTRIBUTED_PREP_REUSE_ONLY": True,
        "_PH_TAG_OVERRIDE": "pH7_0",
        "_COMBO_PREP_CENTER_BY_PH": {"pH7_0": [1.0, 2.0, 3.0]},
        "_COMBO_PREP_BOX_BY_PH": {"pH7_0": [24.0, 24.0, 24.0]},
        "_COMBO_PREP_SOURCE_BY_PH": {"pH7_0": "control"},
        "_COMBO_PREP_CONTROL_STEMS": ["KEU_A510"],
        "_COMBO_PREP_CONTROL_LOOKUP": {"KEU_A510": "foo.sdf"},
    }
    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "PREPPED_LIGANDS_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)

    paths = make_paths(cfg, base_id="BNJS", pdb_file="BNJS.pdb")
    cleaned = paths.receptor_cleaned_pdb("HOLO")
    receptor = paths.receptor_pdbqt("HOLO", "pH7_0")
    cleaned.write_text("ATOM\n", encoding="utf-8")
    receptor.write_text("RECEPTOR\n", encoding="utf-8")

    monkeypatch.setenv("APO_HOLO_VARIANT", "HOLO")
    monkeypatch.setattr(
        runtime_process,
        "_phase2_to4_receptor_and_center",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("receptor prep should be skipped")
        ),
    )
    monkeypatch.setattr(
        runtime_process,
        "_phase5_ph_ensemble_global",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("pH ensemble rebuild should be skipped")
        ),
    )

    calls = {}

    def _phase6(
        cfg_arg,
        paths_arg,
        logger,
        pdb_id,
        variant_env,
        variant_token,
        variant_label,
        legacy_mode,
        cleaned_pdb,
        receptor_pdbqt,
        center,
        box_size,
        center_by_ph,
        box_by_ph,
        center_source_by_ph,
        stages,
        params,
        control_stems,
        control_lookup,
    ):
        calls["pdb_id"] = pdb_id
        calls["variant_token"] = variant_token
        calls["receptor_pdbqt"] = receptor_pdbqt
        calls["center"] = center
        calls["box_size"] = box_size
        calls["control_stems"] = control_stems

    monkeypatch.setattr(runtime_process, "_phase6_to8_ligands_and_docking", _phase6)

    runtime_process.process_one_protein(cfg, "BNJS.pdb", stages=[], params={})  # type: ignore[arg-type]

    assert calls["pdb_id"] == "BNJS"
    assert calls["variant_token"] == "HOLO"
    assert calls["receptor_pdbqt"] == str(receptor)
    assert calls["center"] == (1.0, 2.0, 3.0)
    assert calls["box_size"] == (24.0, 24.0, 24.0)
    assert calls["control_stems"] == ["KEU_A510"]
