"""CCD-backed ligand chemistry helpers for protein-prep validation."""

from __future__ import annotations

import logging
import shlex
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from protein_prep.pdb_records import line_xyz


@dataclass(frozen=True)
class CcdAtom:
    atom_id: str
    element: str


@dataclass(frozen=True)
class CcdInstanceInput:
    token: str
    sdf_path: Path
    atoms: list[CcdAtom]
    coords: dict[str, tuple[float, float, float]]


def download_ccd_sdf(resname: str, work_dir: Path) -> Path | None:
    """Download an RCSB CCD ideal SDF for ``resname`` into ``work_dir``."""

    token = resname.strip().upper()
    if not token:
        return None
    return _download_ccd_file(
        token,
        work_dir / f"{token}_ideal.sdf",
        f"https://files.rcsb.org/ligands/download/{token}_ideal.sdf",
    )


def download_ccd_cif(resname: str, work_dir: Path) -> Path | None:
    """Download an RCSB CCD CIF for atom-name reconciliation."""

    token = resname.strip().upper()
    if not token:
        return None
    return _download_ccd_file(
        token,
        work_dir / f"{token}.cif",
        f"https://files.rcsb.org/ligands/download/{token}.cif",
    )


def write_ccd_instance_sdf(
    ligand_pdb: Path,
    resname: str,
    work_dir: Path,
    output_sdf: Path,
    *,
    selected_smiles: str | None = None,
) -> dict[str, object]:
    """Write an SDF with CCD bond orders and crystal heavy-atom coordinates.

    RCSB ideal SDF files carry chemically curated bond orders but not PDB atom
    names.  The companion CCD CIF carries atom names in the same component atom
    order.  We reconcile those names to the extracted crystal ligand PDB and
    then write an SDF suitable for Meeko ligand preparation.
    """

    loaded, error = _load_instance_inputs(ligand_pdb, resname, work_dir)
    if loaded is None:
        return error or _instance_status("missing_ccd_inputs")
    built, error = _build_crystal_heavy_mol(loaded)
    if built is None:
        return error or _instance_status("ccd_instance_build_failed")
    built, state_status = _apply_selected_smiles_template(built, selected_smiles)
    status = _write_instance_sdf(built, output_sdf)
    if status != "ok":
        return _instance_status(status)
    if not output_sdf.exists() or output_sdf.stat().st_size == 0:
        return _instance_status("write_failed")
    return _instance_status(
        "ok",
        ccd_resname=loaded.token,
        ccd_heavy_atoms=_mol_atom_count(built),
        output_sdf=str(output_sdf),
        selected_smiles_template_status=state_status,
    )


def _load_instance_inputs(
    ligand_pdb: Path,
    resname: str,
    work_dir: Path,
) -> tuple[CcdInstanceInput | None, dict[str, object] | None]:
    token = resname.strip().upper() or _first_ligand_resname(ligand_pdb)
    if not token:
        return None, _instance_status("missing_resname")
    sdf_path = download_ccd_sdf(token, work_dir)
    cif_path = download_ccd_cif(token, work_dir)
    if sdf_path is None:
        return None, _instance_status("missing_ccd_sdf")
    if cif_path is None:
        return None, _instance_status("missing_ccd_cif")
    atoms = parse_ccd_atom_order(cif_path)
    coords = _ligand_atom_coords(ligand_pdb)
    if not atoms:
        return None, _instance_status("missing_ccd_atom_order")
    if not coords:
        return None, _instance_status("missing_ligand_coords")
    return CcdInstanceInput(token, sdf_path, atoms, coords), None


def _build_crystal_heavy_mol(
    loaded: CcdInstanceInput,
) -> tuple[object | None, dict[str, object] | None]:
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except Exception as exc:
        return None, _instance_status(f"rdkit_unavailable:{str(exc)[:120]}")
    mol = _read_first_sdf_mol(loaded.sdf_path)
    if mol is None:
        return None, _instance_status("rdkit_ccd_sdf_failed")
    heavy_atoms = [atom for atom in loaded.atoms if atom.element.upper() != "H"]
    heavy_mol = Chem.RemoveHs(mol, sanitize=False)
    if heavy_mol.GetNumAtoms() != len(heavy_atoms):
        return None, _instance_status(
            "ccd_atom_count_mismatch",
            ccd_heavy_atoms=len(heavy_atoms),
            mol_heavy_atoms=heavy_mol.GetNumAtoms(),
        )
    atom_names = _ccd_to_crystal_atom_names(heavy_atoms, loaded.coords)
    missing = _missing_crystal_atom_names(heavy_atoms, atom_names, loaded.coords)
    if missing:
        return None, _instance_status(
            "missing_crystal_atom_names",
            missing_atoms=",".join(sorted(missing)[:12]),
            missing_atom_count=len(missing),
        )
    conformer = Chem.Conformer(heavy_mol.GetNumAtoms())
    for idx, crystal_name in enumerate(atom_names):
        xyz = loaded.coords[crystal_name.upper()]
        conformer.SetAtomPosition(idx, Point3D(*xyz))
        _set_rdkit_atom_name(heavy_mol.GetAtomWithIdx(idx), crystal_name)
    heavy_mol.RemoveAllConformers()
    heavy_mol.AddConformer(conformer, assignId=True)
    heavy_mol.SetProp("_Name", loaded.token)
    return heavy_mol, None


def _write_instance_sdf(mol: object, output_sdf: Path) -> str:
    from rdkit import Chem

    try:
        Chem.SanitizeMol(mol)
        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(output_sdf))
        writer.write(Chem.AddHs(mol, addCoords=True))
        writer.close()
    except Exception as exc:
        return f"write_failed:{str(exc)[:120]}"
    return "ok"


def _apply_selected_smiles_template(
    mol: Any,
    selected_smiles: str | None,
) -> tuple[Any, str]:
    if not selected_smiles:
        return mol, "not_requested"
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        template = Chem.MolFromSmiles(selected_smiles)
        if template is None:
            return mol, "selected_smiles_invalid"
        template = Chem.RemoveHs(template, sanitize=True)
        if template.GetNumAtoms() != mol.GetNumAtoms():
            return mol, "selected_smiles_atom_count_mismatch"
        assigned = AllChem.AssignBondOrdersFromTemplate(template, mol)
        Chem.SanitizeMol(assigned)
        return assigned, "applied"
    except Exception as exc:
        return mol, f"selected_smiles_template_failed:{str(exc)[:100]}"


def _mol_atom_count(mol: Any) -> int:
    return int(mol.GetNumAtoms())


def parse_ccd_atom_order(cif_path: Path) -> list[CcdAtom]:
    """Return CCD component atom rows in CIF order."""

    lines = cif_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    for fields, data in _iter_cif_loops(lines):
        try:
            atom_idx = fields.index("_chem_comp_atom.atom_id")
            elem_idx = fields.index("_chem_comp_atom.type_symbol")
        except ValueError:
            continue
        atoms: list[CcdAtom] = []
        for row in data:
            if len(row) <= max(atom_idx, elem_idx):
                continue
            atoms.append(CcdAtom(atom_id=_clean_cif_token(row[atom_idx]), element=_clean_cif_token(row[elem_idx]).upper()))
        return atoms
    return []


def _download_ccd_file(token: str, out_path: Path, url: str) -> Path | None:
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            out_path.write_bytes(response.read())
    except Exception as exc:
        logging.info("[ccd] download failed res=%s url=%s err=%s", token, url, exc)
        return None
    return out_path if out_path.exists() and out_path.stat().st_size > 0 else None


def _iter_cif_loops(lines: Iterable[str]) -> Iterable[tuple[list[str], list[list[str]]]]:
    line_list = list(lines)
    idx = 0
    while idx < len(line_list):
        if line_list[idx].strip() != "loop_":
            idx += 1
            continue
        idx += 1
        fields: list[str] = []
        while idx < len(line_list) and line_list[idx].strip().startswith("_"):
            fields.append(line_list[idx].strip().split()[0])
            idx += 1
        data: list[list[str]] = []
        while idx < len(line_list):
            text = line_list[idx].strip()
            if not text or text == "#":
                idx += 1
                break
            if text == "loop_" or text.startswith("_"):
                break
            data.append(shlex.split(text, posix=True))
            idx += 1
        yield fields, data


def _clean_cif_token(token: str) -> str:
    return token.strip().strip("'\"")


def _ligand_atom_coords(ligand_pdb: Path) -> dict[str, tuple[float, float, float]]:
    coords: dict[str, tuple[float, float, float]] = {}
    if not ligand_pdb.exists():
        return coords
    with ligand_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if _line_element(line) == "H":
                continue
            xyz = line_xyz(line)
            if xyz is None:
                continue
            coords.setdefault(line[12:16].strip().upper(), xyz)
    return coords


def _missing_crystal_atom_names(
    heavy_atoms: list[CcdAtom],
    atom_names: list[str],
    coords: dict[str, tuple[float, float, float]],
) -> list[str]:
    return [
        atom.atom_id
        for atom, crystal_name in zip(heavy_atoms, atom_names)
        if crystal_name.upper() not in coords
    ]


def _set_rdkit_atom_name(atom: Any, atom_name: str) -> None:
    atom.SetProp("atom_name", atom_name)
    atom.SetProp("_TriposAtomName", atom_name)


def _ccd_to_crystal_atom_names(
    heavy_atoms: list[CcdAtom],
    coords: dict[str, tuple[float, float, float]],
) -> list[str]:
    mapping: list[str] = []
    before_base_anchor = True
    for atom in heavy_atoms:
        atom_id = atom.atom_id.upper()
        mapping.append(
            _resolve_crystal_atom_name(
                atom_id,
                coords,
                before_base_anchor=before_base_anchor,
            )
        )
        if atom_id in {"N1", "N9"}:
            before_base_anchor = False
    return mapping


def _resolve_crystal_atom_name(
    atom_id: str,
    coords: dict[str, tuple[float, float, float]],
    *,
    before_base_anchor: bool,
) -> str:
    prime_name = f"{atom_id}'"
    nucleotide_sugar_names = {"C1", "C2", "C3", "C4", "C5", "O2", "O3", "O4", "O5"}
    if before_base_anchor and atom_id in nucleotide_sugar_names and prime_name in coords:
        return prime_name
    return atom_id


def _first_ligand_resname(ligand_pdb: Path) -> str:
    if not ligand_pdb.exists():
        return ""
    with ligand_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                return line[17:20].strip().upper()
    return ""


def _read_first_sdf_mol(sdf_path: Path) -> object | None:
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)
    return next((mol for mol in supplier if mol is not None), None)


def _line_element(line: str) -> str:
    return (line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[:1]).upper()


def _instance_status(status: str, **extra: object) -> dict[str, object]:
    return {"ccd_instance_status": status, **extra}
