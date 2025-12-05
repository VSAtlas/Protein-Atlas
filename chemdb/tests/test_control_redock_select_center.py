import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

# Ensure project modules are importable when running from chemdb/tests
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import docking_controls  # noqa: E402
import run_vina  # noqa: E402


def _write_dummy_pdb(path: Path, coords: list[tuple[float, float, float]]):
    lines = []
    for idx, (x, y, z) in enumerate(coords, start=1):
        lines.append(
            f"ATOM  {idx:5d}  C   LIG A{idx:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
    path.write_text("\n".join(lines))


def _exploding_job(job):
    raise RuntimeError("boom")


@pytest.fixture
def fake_paths(tmp_path):
    lig_dir = tmp_path / "ligands_raw"
    lig_dir.mkdir()
    prepped_dir = tmp_path / "prepped_ligands"
    prepped_dir.mkdir()
    return SimpleNamespace(
        pdb_id="TEST",
        ligand_output_dir=lig_dir,
        prepped_ligands_dir=prepped_dir,
    )


def test_select_center_control_happy(monkeypatch, tmp_path, caplog, fake_paths):
    ctrl_ref = tmp_path / "CTRL1.sdf"
    ctrl_ref.write_text("fake sdf")
    receptor_pdbqt = tmp_path / "receptor.pdbqt"
    receptor_pdbqt.write_text("RECEPTOR")

    _write_dummy_pdb(fake_paths.ligand_output_dir / "CTRL1.pdb", [(10.0, 10.0, 10.0)])
    prepped_ctrl = fake_paths.prepped_ligands_dir / "CTRL1.sanitized.pdbqt"
    prepped_ctrl.write_text("ROOT\nENDROOT\n")

    def fake_emit_vina_config(cfg, pdb_id, receptor_pdbqt, center, box, lig_pdbqt, stage_name, stage_info, threads, logger=None, variant=None, ph_token=None, legacy=False):
        out_dir = tmp_path / "ctrl_redock"
        out_dir.mkdir(exist_ok=True)
        conf_path = out_dir / "vina.conf"
        out_path = out_dir / "docked.pdbqt"
        conf_path.write_text("# dummy")
        out_path.write_text("REMARK VINA RESULT: -7.5 0.0 0.0")
        return str(conf_path), str(out_path)

    def fake_best_model(out_path_str, obabel_cmd):
        best_pdb = tmp_path / "best_ctrl1.pdb"
        best_pdb.write_text("HETATM ...")
        return str(best_pdb), -7.5

    monkeypatch.setattr(docking_controls, "build_control_lookup", lambda _paths: {"CTRL1": ctrl_ref})
    monkeypatch.setattr(docking_controls, "emit_vina_config", fake_emit_vina_config)
    monkeypatch.setattr(docking_controls, "_ctrl_best_model_to_pdb", fake_best_model)
    monkeypatch.setattr(docking_controls, "compute_rmsd", lambda ref, dock: 1.2)
    monkeypatch.setattr(run_vina, "run_docking_task", lambda *args, **kwargs: ("ok", -7.5))

    cfg = {
        "RUN_ID": "TEST",
        "CONTROL_CENTER_POLICY": "best_redock",
        "CPU": 1,
        "MAX_PARALLEL_JOBS": 1,
    }
    logger = logging.getLogger("test-ctrl-redock")
    caplog.set_level(logging.INFO, logger="test-ctrl-redock")

    center, box = docking_controls.select_center_via_control_redock(
        cfg=cfg,
        paths=fake_paths,
        receptor_pdbqt=str(receptor_pdbqt),
        logger=logger,
        variant="HOLO",
        ph_token=None,
        legacy=False,
    )

    assert center is not None
    assert isinstance(center, tuple) and len(center) == 3
    assert box == (24.0, 24.0, 24.0)
    assert any("[control-redock] candidates=" in rec.getMessage() for rec in caplog.records)


def test_select_center_processpool_failure(monkeypatch, tmp_path, caplog, fake_paths):
    ctrl1_ref = tmp_path / "CTRL1.sdf"
    ctrl1_ref.write_text("fake sdf")
    ctrl2_ref = tmp_path / "CTRL2.sdf"
    ctrl2_ref.write_text("fake sdf")

    _write_dummy_pdb(fake_paths.ligand_output_dir / "CTRL1.pdb", [(0.0, 0.0, 0.0)])
    _write_dummy_pdb(fake_paths.ligand_output_dir / "CTRL2.pdb", [(1.0, 1.0, 1.0)])

    (fake_paths.prepped_ligands_dir / "CTRL1.sanitized.pdbqt").write_text("ROOT\nENDROOT\n")
    (fake_paths.prepped_ligands_dir / "CTRL2.sanitized.pdbqt").write_text("ROOT\nENDROOT\n")

    monkeypatch.setattr(docking_controls, "build_control_lookup", lambda _paths: {"CTRL1": ctrl1_ref, "CTRL2": ctrl2_ref})
    monkeypatch.setattr(docking_controls, "_ensure_ctrl_vina_manifest", lambda *args, **kwargs: tmp_path / "manifest.json")
    monkeypatch.setattr(docking_controls, "_finalize_ctrl_vina_manifest", lambda *args, **kwargs: None)
    monkeypatch.setattr(docking_controls, "_ctrl_redock_job", _exploding_job)

    cfg = {
        "RUN_ID": "TEST",
        "CONTROL_CENTER_POLICY": "best_redock",
        "CPU": 2,
        "MAX_PARALLEL_JOBS": 2,
    }
    logger = logging.getLogger("test-ctrl-redock")
    caplog.set_level(logging.ERROR, logger="test-ctrl-redock")

    center, box = docking_controls.select_center_via_control_redock(
        cfg=cfg,
        paths=fake_paths,
        receptor_pdbqt=str(tmp_path / "receptor.pdbqt"),
        logger=logger,
        variant="HOLO",
        ph_token=None,
        legacy=False,
    )

    assert center is None
    assert box is None
    assert any("[control-redock] ProcessPoolExecutor failed" in rec.getMessage() for rec in caplog.records)
