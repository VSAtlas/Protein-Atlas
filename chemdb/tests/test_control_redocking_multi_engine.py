from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src import docking
from path_router import make_paths, docked_dir, receptor_file


def test_control_redocking_multi_engine_runs_per_enabled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "INPUT_DIR": str(tmp_path / "input_pdbs"),
        "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "OUTPUT_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
        "RUN_ID": "ctrl_test",
        "FAST_MODE": True,
        "GNINA_EXE": "/bin/true",
    }
    paths = make_paths(cfg, base_id="TEST", pdb_file="TEST.pdb")
    ligand_path = paths.prepped_ligands_dir / "TEST.sanitized.pdbqt"
    ligand_path.write_text("TEST LIG\n", encoding="utf-8")

    calls = {"mol2": 0, "ledock_receptor": 0, "emit": [], "run": []}

    monkeypatch.setattr(docking, "should_run_ledock_for_target", lambda cfg: True)
    monkeypatch.setattr(docking, "should_run_dock6_for_target", lambda cfg: True)
    monkeypatch.setattr(docking, "should_run_gnina_for_target", lambda td, cfg: True)

    def fake_run_ledock_for_stage(**kwargs):
        ligs = kwargs.get("ligands", [])
        out_root = kwargs.get("output_root") or tmp_path / "ledock_out"
        Path(out_root).mkdir(parents=True, exist_ok=True)
        scores = {}
        metrics = {}
        for lig in ligs:
            lig_path = Path(lig)
            out = Path(out_root) / f"{lig_path.stem}.dok"
            out.write_text("dok\n", encoding="utf-8")
            scores[lig_path] = -1.0
            metrics[lig_path] = {"best_score_kcal": -1.0, "valid": True}
        return scores, metrics, {}

    def fake_run_dock6_for_stage(**kwargs):
        ligs = kwargs.get("ligands", [])
        out_root = kwargs.get("output_root") or tmp_path / "dock6_out"
        Path(out_root).mkdir(parents=True, exist_ok=True)
        ranked = Path(out_root) / "dock6_stage1_ranked.mol2"
        ranked.write_text("ranked\n", encoding="utf-8")
        scores = {}
        metrics = {}
        for lig in ligs:
            lig_path = Path(lig)
            scores[lig_path] = -2.0
            metrics[lig_path] = {"grid_score": -2.0, "n_poses": 1, "valid": True}
        return scores, metrics

    monkeypatch.setattr(docking, "run_ledock_for_stage", fake_run_ledock_for_stage)
    monkeypatch.setattr(docking, "run_dock6_for_stage", fake_run_dock6_for_stage)

    def fake_run_gnina_for_stage(**kwargs):
        paths_obj = kwargs["paths"]
        stage_name = kwargs["stage_name"]
        ph_label = kwargs.get("ph_label")
        variant = kwargs.get("variant")
        stage_dir = paths_obj.docked_stage_dir(variant, f"gnina_{stage_name}", ph_label)
        stage_dir.mkdir(parents=True, exist_ok=True)
        scores = {}
        metrics = {}
        for lig in kwargs.get("ligands", []):
            lig_path = Path(lig)
            out = stage_dir / f"{lig_path.stem}_gnina_{stage_name}.pdbqt"
            out.write_text("gnina pose\n", encoding="utf-8")
            scores[lig] = -3.0
            metrics[lig] = {"gnina_primary_score": -3.0}
        return scores, metrics, {}

    monkeypatch.setattr(docking, "run_gnina_for_stage", fake_run_gnina_for_stage)

    def fake_mol2(cfg_in, ligs, logger):
        calls["mol2"] += 1

    monkeypatch.setattr(docking, "ensure_mol2_for_ledock", fake_mol2)

    def fake_ledock_receptor(cfg_in, pdb_id, variant_token, ph_label, logger):
        calls["ledock_receptor"] += 1
        out = tmp_path / "ledock_receptor.pdb"
        out.write_text("RECEPTOR\n", encoding="utf-8")
        return out

    monkeypatch.setattr(docking, "ensure_ledock_receptor", fake_ledock_receptor)

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
        cfg_dir = Path(cfg_in["OUTPUT_DIR"]) / "configs" / stage_name
        cfg_dir.mkdir(parents=True, exist_ok=True)
        conf = cfg_dir / f"{Path(ligand_path).stem}.txt"
        conf.write_text("conf\n", encoding="utf-8")
        out_root = (
            docked_dir(pdb_id, variant=variant, ph_tag=ph_token, legacy=legacy)
            / stage_name
        )
        out_root.mkdir(parents=True, exist_ok=True)
        out = out_root / f"{Path(ligand_path).stem}.pdbqt"
        out.write_text("pose\n", encoding="utf-8")
        calls["emit"].append((stage_name, ph_token))
        return conf, out

    monkeypatch.setattr(docking, "emit_vina_config", fake_emit)

    def fake_run(exe, conf_path, lig_path, out_path):
        Path(out_path).write_text("ran\n", encoding="utf-8")
        calls["run"].append((exe, out_path))

    monkeypatch.setattr(docking, "run_docking_task", fake_run)

    logger = logging.getLogger("test.control_multi_engine")

    rec = receptor_file("TEST", variant="HOLO", ph_tag="pH7_0")
    rec.parent.mkdir(parents=True, exist_ok=True)
    rec.write_text("RECEPTOR\n", encoding="utf-8")

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
        assert engine_dir.exists()
        assert any(engine_dir.glob("*"))

    assert calls["mol2"] == 1
    assert calls["ledock_receptor"] >= 1
    assert calls["emit"]
    assert calls["run"]
