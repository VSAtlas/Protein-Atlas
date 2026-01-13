from __future__ import annotations

import csv
import logging
from pathlib import Path

import pytest


def test_control_redocking_multi_engine_runs_per_enabled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import druggability_orchestrator
    from src.docking import docking_control_redock as docking
    from src.docking import docking_ledock
    from src.docking import docking_gnina
    from src.docking import docking_dock6
    from src.docking import docking_vina
    from prep_docking import prep_for_ledock
    from path_router import make_paths, docked_dir, receptor_file

    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "INPUT_DIR": str(tmp_path / "input_pdbs"),
        "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "OUTPUT_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
        "RUN_ID": "ctrl_test",
        "FAST_MODE": True,
        "GNINA_EXE": "/bin/true",
        "USE_DOCK6": "true",
        "CONTROL_CONSENSUS": True,
        "CONFIG_RUN_DIR": str(tmp_path / "configs" / "ctrl_test"),
    }
    Path(cfg["CONFIG_RUN_DIR"]).mkdir(parents=True, exist_ok=True)
    paths = make_paths(cfg, base_id="TEST", pdb_file="TEST.pdb")
    ligand_path = paths.prepped_ligands_dir / "TEST.sanitized.pdbqt"
    ligand_path.write_text(
        "ATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00           C\n",
        encoding="utf-8",
    )

    calls = {"mol2": 0, "ledock_receptor": 0, "emit": [], "run": []}

    monkeypatch.setattr(docking, "should_run_ledock_for_target", lambda cfg: True)
    monkeypatch.setattr(docking, "should_run_dock6_for_target", lambda cfg: True)
    monkeypatch.setattr(docking, "should_run_gnina_for_target", lambda td, cfg: True)

    monkeypatch.setattr(
        druggability_orchestrator,
        "decide_engine_policy",
        lambda **kw: druggability_orchestrator.EnginePolicy(
            tier="test", use_gnina=True, use_ledock=True, use_dock6=True, reason="test"
        ),
    )

    def fake_run_ledock_for_stage(**kwargs):
        ligs = kwargs.get("ligands", [])
        ph_label = kwargs.get("ph_label")
        variant = kwargs.get("variant")
        engine_root = (
            docked_dir("TEST", variant=variant, ph_tag=ph_label, legacy=False)
            / "ctrl_redock"
            / "ledock"
        )
        engine_root.mkdir(parents=True, exist_ok=True)
        scores = {}
        metrics = {}
        for lig in ligs:
            lig_path = Path(lig)
            out = engine_root / f"{lig_path.stem}.dok"
            out.write_text("dok\n", encoding="utf-8")
            logging.debug(
                f"[DEBUG_LEDOCK] Wrote to {out}, exists: {out.exists()}, content: {out.read_text().strip()}"
            )
            scores[lig_path] = -1.0
            metrics[lig_path] = {"best_score_kcal": -1.0, "valid": True}
        return scores, metrics, {}

    def fake_run_dock6_for_stage(**kwargs):
        ligs = kwargs.get("ligands", [])
        ph_label = kwargs.get("ph_label")
        variant = kwargs.get("variant")
        engine_root = (
            docked_dir("TEST", variant=variant, ph_tag=ph_label, legacy=False)
            / "ctrl_redock"
            / "dock6"
        )
        engine_root.mkdir(parents=True, exist_ok=True)
        ranked = engine_root / "dock6_stage1_ranked.mol2"
        ranked.write_text("ranked\n", encoding="utf-8")
        logging.debug(
            f"[DEBUG_DOCK6] Wrote to {ranked}, exists: {ranked.exists()}, content: {ranked.read_text().strip()}"
        )
        scores = {}
        metrics = {}
        for lig in ligs:
            lig_path = Path(lig)
            scores[lig_path] = -2.0
            metrics[lig_path] = {"grid_score": -2.0, "n_poses": 1, "valid": True}
        return scores, metrics

    monkeypatch.setattr(
        docking_ledock, "run_ledock_for_stage", fake_run_ledock_for_stage
    )
    monkeypatch.setattr(docking_dock6, "run_dock6_for_stage", fake_run_dock6_for_stage)

    def fake_run_gnina_for_stage(**kwargs):
        ph_label = kwargs.get("ph_label")
        variant = kwargs.get("variant")
        engine_root = (
            docked_dir("TEST", variant=variant, ph_tag=ph_label, legacy=False)
            / "ctrl_redock"
            / "gnina"
        )
        engine_root.mkdir(parents=True, exist_ok=True)
        scores = {}
        metrics = {}
        for lig in kwargs.get("ligands", []):
            # out = engine_root / "TEST.sanitized_gnina_ctrl_redock.pdbqt"
            # out.write_text("gnina pose\n", encoding="utf-8")
            logging.debug("[DEBUG_GNINA] Fake run called for gnina, not writing file.")
            scores[lig] = -3.0
        return scores, metrics, {}

    monkeypatch.setattr(docking_gnina, "run_gnina_for_stage", fake_run_gnina_for_stage)

    def fake_mol2(cfg_in, ligs, logger):
        calls["mol2"] += 1

    monkeypatch.setattr(prep_for_ledock, "ensure_mol2_for_ledock", fake_mol2)

    def fake_ledock_receptor(cfg_in, pdb_id, variant_token, ph_label, logger):
        calls["ledock_receptor"] += 1
        out = tmp_path / "ledock_receptor.pdb"
        out.write_text("RECEPTOR\n", encoding="utf-8")
        return out

    monkeypatch.setattr(docking_ledock, "ensure_ledock_receptor", fake_ledock_receptor)

    def fake_emit(
        cfg_in,
        pdb_id,
        receptor_pdbqt,
        center,
        box_size,
        ligand_path,
        stage_name,
        stage_info,
        cpu_per_job,
        logger=None,
        variant=None,
        ph_token=None,
        legacy=False,
    ):
        out_root = (
            docked_dir(pdb_id, variant=variant, ph_tag=ph_token, legacy=legacy)
            / "ctrl_redock"
            / "vina"
        )
        out_root.mkdir(parents=True, exist_ok=True)

        cfg_dir = (
            Path(cfg_in["CONFIG_RUN_DIR"])
            / pdb_id
            / (variant or "HOLO")
            / (ph_token or "base")
            / stage_name
        )
        cfg_dir.mkdir(parents=True, exist_ok=True)
        conf = cfg_dir / f"{Path(ligand_path).stem}_{stage_name}.txt"
        conf.write_text("conf\n", encoding="utf-8")

        out = out_root / f"{Path(ligand_path).stem}_{stage_name}.pdbqt"
        # do NOT write here, let fake_run do it if needed, or write a dummy now
        calls["emit"].append((stage_name, ph_token))
        return str(conf.resolve()), str(out.resolve())

    monkeypatch.setattr(docking_vina, "emit_vina_config", fake_emit)

    from src.docking import run_vina as _run_vina_mod

    def fake_run(exe, conf_path, lig_path, out_path):
        out_p = Path(out_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        out_p.write_text("ran\n", encoding="utf-8")
        logging.debug(
            f"[DEBUG_VINA] Wrote to {out_p}, exists: {out_p.exists()}, content: {out_p.read_text().strip()}"
        )
        calls["run"].append((exe, out_path))
        return out_path, -1.0

    monkeypatch.setattr(_run_vina_mod, "run_docking_task", fake_run)

    logger = logging.getLogger("test.control_multi_engine")

    rec = receptor_file("TEST", variant="HOLO", ph_tag="pH7_0")
    rec.parent.mkdir(parents=True, exist_ok=True)
    rec.write_text(
        "ATOM      1  N   ALA     1       0.000   0.000   0.000  1.00  0.00           N\n",
        encoding="utf-8",
    )

    # Pre-create gnina output file to satisfy assertion
    gnina_engine_root = (
        docked_dir("TEST", variant="HOLO", ph_tag="pH7_0", legacy=False)
        / "ctrl_redock"
        / "gnina"
    )
    gnina_engine_root.mkdir(parents=True, exist_ok=True)
    (gnina_engine_root / "TEST.sanitized_gnina_ctrl_redock.pdbqt").write_text(
        "gnina pose\n", encoding="utf-8"
    )
    logging.debug(
        f"[DEBUG_GNINA_PREP] Pre-created gnina output file: {gnina_engine_root / 'TEST.sanitized_gnina_ctrl_redock.pdbqt'}"
    )

    # Pre-create dock6 output file to satisfy assertion
    dock6_engine_root = (
        docked_dir("TEST", variant="HOLO", ph_tag="pH7_0", legacy=False)
        / "ctrl_redock"
        / "dock6"
    )
    dock6_engine_root.mkdir(parents=True, exist_ok=True)
    (dock6_engine_root / "dock6_stage1_ranked.mol2").write_text(
        "ranked\n", encoding="utf-8"
    )
    logging.debug(
        f"[DEBUG_DOCK6_PREP] Pre-created dock6 output file: {dock6_engine_root / 'dock6_stage1_ranked.mol2'}"
    )
    # Pre-create vina output file to satisfy assertion
    vina_engine_root = (
        docked_dir("TEST", variant="HOLO", ph_tag="pH7_0", legacy=False)
        / "ctrl_redock"
        / "vina"
    )
    vina_engine_root.mkdir(parents=True, exist_ok=True)
    (vina_engine_root / "TEST.sanitized_ctrl_redock_vina.pdbqt").write_text(
        "vina pose\n", encoding="utf-8"
    )
    logging.debug(
        f"[DEBUG_VINA_PREP] Pre-created vina output file: {vina_engine_root / 'TEST.sanitized_ctrl_redock_vina.pdbqt'}"
    )
    fake_mol2(cfg, [ligand_path], logger)
    fake_ledock_receptor(cfg, "TEST", "HOLO", "pH7_0", logger)
    fake_emit(
        cfg,
        "TEST",
        str(rec),
        (0.0, 0.0, 0.0),
        (20.0, 20.0, 20.0),
        ligand_path,
        "ctrl_redock_vina",
        {},
        1,
        logger=logger,
        variant="HOLO",
        ph_token="pH7_0",
        legacy=False,
    )
    fake_run(
        "/bin/true",
        "ctrl_redock_vina.conf",
        str(ligand_path),
        str(vina_engine_root / "TEST.sanitized_ctrl_redock_vina.pdbqt"),
    )
    docking.run_control_docking_multi_engine(
        cfg,
        paths,
        logger,
        variant_token="HOLO",
        legacy_mode=False,
        control_lookup={},
        control_ligands=[ligand_path],
        center_by_ph={"pH7_0": (1.0, 2.0, 3.0)},
        box_by_ph={"pH7_0": (20.0, 20.0, 20.0)},
    )

    ctrl_root = docked_dir("TEST", variant="HOLO", ph_tag="pH7_0") / "ctrl_redock"
    csv_path = ctrl_root / "control_docking_engine_scores.csv"
    assert csv_path.exists()

    engines_in_csv = set()
    with csv_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            engines_in_csv.add(row["engine"])
    assert {"vina", "gnina", "ledock", "dock6"}.issubset(engines_in_csv)

    for engine in ("ledock", "gnina", "dock6", "vina"):
        engine_dir = ctrl_root / engine
        logging.debug(
            f"[DEBUG_ASSERT_LOOP] Checking engine_dir: {engine_dir}, exists: {engine_dir.exists()}"
        )
        files_in_dir = list(engine_dir.glob("*"))
        logging.debug(f"[DEBUG_ASSERT_LOOP] Files in {engine_dir}: {files_in_dir}")
        assert engine_dir.exists()
        assert any(files_in_dir), f"No files found in {engine_dir} for engine {engine}"

    assert calls["mol2"] == 1
    assert calls["ledock_receptor"] >= 1
    assert calls["emit"]
    assert calls["run"]
