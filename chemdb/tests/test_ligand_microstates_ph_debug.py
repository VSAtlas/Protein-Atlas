import json
import shutil
import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from prep_ligands.prep_ligands_bulk import prep_ligands_with_mgltools, read_config
from prep_ligands.prep_ligands_microstates import (
    compute_microstate_id_from_pdbqt,
    enumerate_ligands_for_docking,
)


def _write_ph_debug_sdf(sdf_path: Path) -> None:
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    ligands = [
        ("ph_debug_acetic_acid", "CC(=O)O"),
        ("ph_debug_dimethylamine", "CN(C)C"),
        ("ph_debug_imidazole", "C1=NC=CN1"),
        ("ph_debug_aniline", "Nc1ccccc1"),
    ]
    writer = Chem.SDWriter(str(sdf_path))
    for name, smi in ligands:
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        mol.SetProp("_Name", name)
        writer.write(mol)
    writer.close()


def _ph_label(ph: float) -> str:
    return f"pH{ph:.1f}".replace(".", "_")


@pytest.mark.slow
def test_ph_debug_microstates_match_pdbqts_and_registry():
    root = Path(__file__).resolve().parents[2]
    lib_name = "ph_debug"
    sdf_path = root / "extracted_ligands" / lib_name / f"{lib_name}.sdf"
    lib_root = root / "prepped_ligands" / lib_name

    # Fresh inputs/outputs for deterministic assertions
    _write_ph_debug_sdf(sdf_path)
    for sub in ["microstates", "pH3_0", "pH7_4", "pH10_0"]:
        target = lib_root / sub
        if target.exists():
            shutil.rmtree(target, ignore_errors=True)
    for f in ["microstates.json", "ligand_prep_status.tsv"]:
        candidate = lib_root / f
        if candidate.exists():
            candidate.unlink()

    prep_ligands_with_mgltools(
        force=True,
        ph_values=[3.0, 7.4, 10.0],
        microstate_dedup=True,
        in_sdf_dir=sdf_path.parent,
        root_dir=lib_root,
    )

    ph_values = [3.0, 7.4, 10.0]
    all_ligand_id_maps = []
    for idx in range(1, 5):
        stem = f"{lib_name}_{idx:05d}"
        ph_to_id = {}
        for ph in ph_values:
            label = _ph_label(ph)
            pdbqt = lib_root / label / f"phpH{label[2:]}_debug_{idx:05d}.pdbqt"
            if pdbqt.is_file():
                ph_to_id[ph] = compute_microstate_id_from_pdbqt(pdbqt)
        all_ligand_id_maps.append(ph_to_id)

    assert any(len(set(m.values())) == 1 for m in all_ligand_id_maps if m)
    assert any(len(set(m.values())) > 1 for m in all_ligand_id_maps if m)

    reg = json.loads((lib_root / "microstates.json").read_text())
    entries = reg.get("microstates", []) or []
    micro_by_id = {e["microstate_id"]: e for e in entries}

    for idx, ph2id in enumerate(all_ligand_id_maps, start=1):
        stem = f"{lib_name}_{idx:05d}"
        for ph_value, mid in ph2id.items():
            assert mid in micro_by_id
            e = micro_by_id[mid]
            aliases = e.get("aliases") or []
            lig_aliases = [a for a in aliases if a.get("ligand_stem") == stem]
            assert any(abs(a.get("ph_value") - ph_value) < 1e-3 for a in lig_aliases)

        ids_for_ligand_from_registry = {
            e["microstate_id"]
            for e in entries
            for a in (e.get("aliases") or [])
            if a.get("ligand_stem") == stem
        }
        assert ids_for_ligand_from_registry.issubset(set(ph2id.values()))

    cfg = read_config()
    path_to_entry = {
        (lib_root / e["pdbqt_path"]).resolve(): e
        for e in entries
        if e.get("pdbqt_path")
    }
    for ph in ph_values:
        pdbqts = enumerate_ligands_for_docking(
            requested_ph_values=[ph],
            cfg=cfg,
            pdb_id="TEST",
            root_dir=lib_root,
            microstate_dedup=True,
            force=False,
        )
        for p in pdbqts:
            assert p.resolve() in path_to_entry
            e = path_to_entry[p.resolve()]
            aliases = e.get("aliases") or []
            assert any(abs(a.get("ph_value") - ph) < 1e-3 for a in aliases)

        stems_for_ph = {
            a.get("ligand_stem")
            for e in entries
            for a in (e.get("aliases") or [])
            if abs(a.get("ph_value") - ph) < 1e-3
        }
        assert len(pdbqts) == len(stems_for_ph)
