import builtins
import os
from pathlib import Path

import main
from post_docking.mmgbsa.run_mmgbsa import parse_mmpbsa_delta_total


def _write_mmpbsa_outputs(work_dir: Path, score: float, cfg: dict) -> dict:
    out_dat_name = cfg.get("MMGBSA_MMPBSA_OUT_DAT", "FINAL_RESULTS_MMPBSA.dat")
    out_csv_name = cfg.get("MMGBSA_MMPBSA_OUT_CSV", "FINAL_RESULTS_MMPBSA.csv")
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dat = work_dir / out_dat_name
    out_csv = work_dir / out_csv_name
    out_dat.write_text(f"FINAL RESULTS\nDELTA TOTAL       {score:.3f}\n", encoding="utf-8")
    out_csv.write_text(
        "DELTA Energy Terms\nVDWAALS,EELEC,EGB,ESURF,EPB,ECAVITY,DELTA TOTAL\n"
        f"0,0,0,0,0,0,{score:.3f}\n",
        encoding="utf-8",
    )
    return {"out_dat": str(out_dat), "out_csv": str(out_csv), "log_path": str(work_dir / "mmpbsa.log"), "input_path": str(work_dir / "mmpbsa.in"), "enabled": True}


def test_five_replicate_md_enabled(monkeypatch, tmp_path):
    cfg = {"MD_FIVE_REPLICATE": True}
    topo = {
        "complex_prmtop": str(tmp_path / "complex.prmtop"),
        "complex_inpcrd": str(tmp_path / "complex.inpcrd"),
        "receptor_prmtop": str(tmp_path / "receptor.prmtop"),
        "ligand_prmtop": str(tmp_path / "ligand.prmtop"),
    }
    for path in topo.values():
        Path(path).write_text("dummy", encoding="utf-8")

    def fake_md(complex_prmtop, complex_inpcrd, out_dir, cfg, replicate_index, seed, force, run):
        rep_dir = Path(out_dir) / f"rep{replicate_index}" / "md"
        rep_dir.mkdir(parents=True, exist_ok=True)
        traj_path = rep_dir / "prod.nc"
        traj_path.write_text(f"traj {replicate_index}", encoding="utf-8")
        return {"ok": True, "traj_path": str(traj_path)}

    monkeypatch.setattr(main, "run_implicit_md", fake_md)

    def fake_mmgbsa(complex_prmtop, receptor_prmtop, ligand_prmtop, trajectory_path, work_dir, cfg, force, run):
        rep_idx = int(Path(work_dir).parent.name.replace("rep", ""))
        score = 1.0 + rep_idx
        return _write_mmpbsa_outputs(Path(work_dir), score, cfg)

    monkeypatch.setattr(main, "run_mmgbsa", fake_mmgbsa)

    result = main._mmgbsa_five_replicate_runner(
        topo=topo,
        out_dir=tmp_path / "work",
        cfg=cfg,
        force=True,
        md_enabled=True,
        stage_dir="stage1",
        ligand_stem="ligA",
        pdb_id="EX4M",
        variant_dir="HOLO",
        ph_label="pH7_0",
        run_id="runX",
        logger=main.logging.getLogger("test"),
    )

    assert result["ok"] is True
    assert Path(result["results_csv"]).exists()
    mean_score = parse_mmpbsa_delta_total(Path(result["results_csv"]))
    assert mean_score == sum(1.0 + i for i in range(1, 6)) / 5
    for idx in range(1, 6):
        rep_csv = tmp_path / "work" / f"rep{idx}" / "mmpbsa" / "FINAL_RESULTS_MMPBSA.csv"
        assert rep_csv.exists()


def test_five_replicate_md_disabled(monkeypatch, tmp_path):
    cfg = {"MD_FIVE_REPLICATE": True}
    topo = {
        "complex_prmtop": str(tmp_path / "complex.prmtop"),
        "complex_inpcrd": str(tmp_path / "complex.inpcrd"),
        "receptor_prmtop": str(tmp_path / "receptor.prmtop"),
        "ligand_prmtop": str(tmp_path / "ligand.prmtop"),
    }
    for path in topo.values():
        Path(path).write_text("dummy", encoding="utf-8")

    def fake_traj(complex_prmtop, complex_inpcrd, out_dir, cfg, force, run_cpptraj):
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        traj_path = out_dir / "mdcrd"
        traj_path.write_text("traj", encoding="utf-8")
        return {"enabled": True, "trajout_path": str(traj_path)}

    monkeypatch.setattr(main, "make_mmgbsa_trajectory", fake_traj)

    def fake_mmgbsa(complex_prmtop, receptor_prmtop, ligand_prmtop, trajectory_path, work_dir, cfg, force, run):
        # single score reused for clones
        return _write_mmpbsa_outputs(Path(work_dir), 2.5, cfg)

    monkeypatch.setattr(main, "run_mmgbsa", fake_mmgbsa)

    result = main._mmgbsa_five_replicate_runner(
        topo=topo,
        out_dir=tmp_path / "work",
        cfg=cfg,
        force=True,
        md_enabled=False,
        stage_dir="stage1",
        ligand_stem="ligB",
        pdb_id="EX4M",
        variant_dir="HOLO",
        ph_label="pH7_0",
        run_id="runY",
        logger=main.logging.getLogger("test"),
    )

    assert result["ok"] is True
    agg_csv = Path(result["results_csv"])
    assert agg_csv.exists()
    assert parse_mmpbsa_delta_total(agg_csv) == 2.5
    for idx in range(1, 6):
        rep_csv = tmp_path / "work" / f"rep{idx}" / "mmpbsa" / "FINAL_RESULTS_MMPBSA.csv"
        assert rep_csv.exists()
