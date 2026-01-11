import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import post_docking.mmgbsa.prep_for_mmgbsa as prep


class DummyProc:
    def __init__(self, returncode: int) -> None:
        self.returncode = returncode


def _write_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _setup_paths(
    tmp_path: Path, ligand_name: str
) -> tuple[Path, Path, Path, Path, str]:
    base = tmp_path / "post_docked" / "RUNTEST" / "TESTP" / "HOLO" / "pH7_0"
    stage_dir = "stage1"
    sdf_path = base / stage_dir / f"{ligand_name}.sdf"
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    sdf_path.write_text("SDF", encoding="utf-8")

    mmgbsa_dir = base / "mmgbsa"
    mol2_path = mmgbsa_dir / "mol2" / stage_dir / f"{ligand_name}.mol2"
    frcmod_path = mmgbsa_dir / "frcmod" / stage_dir / f"{ligand_name}.frcmod"
    mol2_path.parent.mkdir(parents=True, exist_ok=True)
    frcmod_path.parent.mkdir(parents=True, exist_ok=True)

    return sdf_path, mol2_path, frcmod_path, mmgbsa_dir, stage_dir


def test_mmgbsa_bcc_sweep_then_success(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    sdf_path, mol2_path, frcmod_path, mmgbsa_dir, stage_dir = _setup_paths(
        tmp_path, "lig1"
    )

    monkeypatch.setattr(
        prep,
        "_select_ambertools",
        lambda amber_prefix, logger: ([], "antechamber", "parmchk2", "PATH"),
    )
    monkeypatch.setattr(prep, "_read_rdkit_info", lambda *args, **kwargs: (0, 2))

    def fake_run(cmd, cwd=None, stdout=None, stderr=None):
        if "antechamber" in cmd[0]:
            method = cmd[cmd.index("-c") + 1]
            charge = int(cmd[cmd.index("-nc") + 1])
            out_path = Path(cmd[cmd.index("-o") + 1])
            if method == "bcc" and charge == 1:
                _write_file(out_path, "mol2")
                return DummyProc(0)
            return DummyProc(1)
        if "parmchk2" in cmd[0]:
            out_path = Path(cmd[cmd.index("-o") + 1])
            _write_file(out_path, "frcmod")
            return DummyProc(0)
        return DummyProc(1)

    monkeypatch.setattr(prep.subprocess, "run", fake_run)

    cfg = {
        "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD": "bcc",
        "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD": "gas",
        "MMGBSA_LIGAND_NOMINAL_NET_CHARGE": 0,
        "MMGBSA_LIGAND_BCC_CHARGE_SWEEP": "-1,1",
        "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2": False,
        "MMGBSA_LIGAND_SQM_LEVEL": 2,
        "MMGBSA_RDKit_VALIDATE": True,
        "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD": 1,
    }

    result = prep.parameterize_ligand_with_fallback(
        sdf_path=str(sdf_path),
        out_mol2=str(mol2_path),
        out_frcmod=str(frcmod_path),
        cfg=cfg,
        stage_dir=stage_dir,
        ligand_stem="lig1",
        force=False,
    )

    assert result["ok"] is True
    metadata_path = mol2_path.with_suffix(".mmgbsa_prep.json")
    assert metadata_path.exists()
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert metadata["charge_method"] == "bcc"
    assert metadata["net_charge_used"] == 1
    assert metadata["low_confidence"] is True
    assert any(
        attempt["nc"] == 1 and attempt["ok"] for attempt in metadata["bcc_attempts"]
    )

    log_dir = mmgbsa_dir / "logs" / stage_dir / "lig1"
    assert (log_dir / "antechamber_attempt_bcc_nc0.log").exists()
    assert (log_dir / "antechamber_attempt_bcc_nc-1.log").exists()
    assert (log_dir / "antechamber_attempt_bcc_nc1.log").exists()
    assert (log_dir / "parmchk2.log").exists()


def test_mmgbsa_gas_fallback(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    sdf_path, mol2_path, frcmod_path, mmgbsa_dir, stage_dir = _setup_paths(
        tmp_path, "lig2"
    )

    monkeypatch.setattr(
        prep,
        "_select_ambertools",
        lambda amber_prefix, logger: ([], "antechamber", "parmchk2", "PATH"),
    )
    monkeypatch.setattr(prep, "_read_rdkit_info", lambda *args, **kwargs: (0, 0))

    def fake_run(cmd, cwd=None, stdout=None, stderr=None):
        if "antechamber" in cmd[0]:
            method = cmd[cmd.index("-c") + 1]
            charge = int(cmd[cmd.index("-nc") + 1])
            out_path = Path(cmd[cmd.index("-o") + 1])
            if method == "gas" and charge == 0:
                _write_file(out_path, "mol2")
                return DummyProc(0)
            return DummyProc(1)
        if "parmchk2" in cmd[0]:
            out_path = Path(cmd[cmd.index("-o") + 1])
            _write_file(out_path, "frcmod")
            return DummyProc(0)
        return DummyProc(1)

    monkeypatch.setattr(prep.subprocess, "run", fake_run)

    cfg = {
        "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD": "bcc",
        "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD": "gas",
        "MMGBSA_LIGAND_NOMINAL_NET_CHARGE": 0,
        "MMGBSA_LIGAND_BCC_CHARGE_SWEEP": "-1,1",
        "MMGBSA_LIGAND_SQM_LEVEL": 2,
        "MMGBSA_RDKit_VALIDATE": True,
    }

    result = prep.parameterize_ligand_with_fallback(
        sdf_path=str(sdf_path),
        out_mol2=str(mol2_path),
        out_frcmod=str(frcmod_path),
        cfg=cfg,
        stage_dir=stage_dir,
        ligand_stem="lig2",
        force=False,
    )

    assert result["ok"] is True
    metadata_path = mol2_path.with_suffix(".mmgbsa_prep.json")
    assert metadata_path.exists()
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert metadata["charge_method"] == "gas"
    assert metadata["gas_fallback_used"] is True
    assert metadata["net_charge_used"] == 0

    log_dir = mmgbsa_dir / "logs" / stage_dir / "lig2"
    assert (log_dir / "antechamber_attempt_gas_nc0.log").exists()
    assert (log_dir / "parmchk2.log").exists()
