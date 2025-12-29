import os
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import protein_prep_mmgbsa


def _pdb_line(record, serial, name, resname, chain, resseq, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain:1}"
        f"{resseq:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00          {element:>2}\n"
    )


def _tleap_available():
    env_prefix = Path(os.environ.get("AMBERTOOLS_PREFIX", "")).expanduser() if os.environ.get("AMBERTOOLS_PREFIX") else None
    if env_prefix:
        if "/micromamba/pkgs" in env_prefix.as_posix():
            return False
        if (env_prefix / "bin" / "tleap").is_file():
            return shutil.which("micromamba") is not None or Path("/home/michael/atlas/micromamba/bin/micromamba").is_file()
    return shutil.which("tleap") is not None


def test_mmgbsa_leap_generation_no_tleap(tmp_path):
    base_dir = tmp_path / "post_docked" / "RUNTEST" / "TESTP" / "HOLO" / "pH7_0" / "stage1"
    base_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdb = base_dir / "receptor.pdb"
    receptor_pdb.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
                _pdb_line("ATOM", 2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    ligand_sdf = base_dir / "lig1.sdf"
    ligand_sdf.write_text(
        "lig1\n  -OEChem-09222100002D\n\n  0  0  0  0  0  0  0  0  0  0  0  0 V2000\nM  END\n$$$$\n",
        encoding="utf-8",
    )

    mmgbsa_dir = base_dir.parent / "mmgbsa"
    mol2_path = mmgbsa_dir / "mol2" / "stage1" / "lig1.mol2"
    frcmod_path = mmgbsa_dir / "frcmod" / "stage1" / "lig1.frcmod"
    mol2_path.parent.mkdir(parents=True, exist_ok=True)
    frcmod_path.parent.mkdir(parents=True, exist_ok=True)
    mol2_path.write_text("@<TRIPOS>MOLECULE\nlig1\n", encoding="utf-8")
    frcmod_path.write_text("MASS\n", encoding="utf-8")

    result = protein_prep_mmgbsa.prep_mmgbsa_receptor_and_topologies(
        pdb_path=str(receptor_pdb),
        runid="RUNTEST",
        center=(0.0, 0.0, 0.0),
        radius=5.0,
        ligand_sdf_paths=[str(ligand_sdf)],
        force=False,
        run_tleap=False,
    )

    topo_root = Path(result["topology_root"])
    receptor_leap = topo_root / "receptor" / "build_receptor.leap"
    ligand_leap = topo_root / "stage1" / "lig1" / "build.leap"

    assert receptor_leap.exists()
    assert ligand_leap.exists()

    text = ligand_leap.read_text(encoding="utf-8")
    assert "source leaprc.protein.ff14SB" in text
    assert "source leaprc.gaff2" in text
    assert "lig1.frcmod" in text
    assert "lig1.mol2" in text
    assert "TESTP_mmgbsa_receptor.pdb" in text
    assert "saveamberparm COM complex.prmtop complex.inpcrd" in text


def test_mmgbsa_leap_generation_with_tleap(tmp_path):
    if not _tleap_available():
        pytest.skip("tleap not available")

    base_dir = tmp_path / "post_docked" / "RUNTEST" / "TESTP" / "HOLO" / "pH7_0" / "stage1"
    base_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdb = base_dir / "receptor.pdb"
    receptor_pdb.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
                _pdb_line("ATOM", 2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
                _pdb_line("ATOM", 3, "C", "ALA", "A", 1, 2.0, 1.5, 0.0, "C"),
                _pdb_line("ATOM", 4, "O", "ALA", "A", 1, 2.0, 2.5, 0.0, "O"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    result = protein_prep_mmgbsa.prep_mmgbsa_receptor_and_topologies(
        pdb_path=str(receptor_pdb),
        runid="RUNTEST",
        center=(0.0, 0.0, 0.0),
        radius=5.0,
        ligand_sdf_paths=[],
        force=False,
        run_tleap=True,
    )

    topo_root = Path(result["topology_root"])
    receptor_dir = topo_root / "receptor"
    assert (receptor_dir / "receptor.prmtop").exists()
    assert (receptor_dir / "receptor.inpcrd").exists()
