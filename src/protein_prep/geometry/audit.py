"""Conservative heavy-atom geometry clash auditing and cleanup."""

from __future__ import annotations

import json
import math
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from protein_prep.pdb_records import line_xyz

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
COMMON_IONS = {
    "AG",
    "AL",
    "AU",
    "BA",
    "CL",
    "BR",
    "IOD",
    "NA",
    "K",
    "CS",
    "LI",
    "CA",
    "MG",
    "MN",
    "ZN",
    "FE",
    "CU",
    "CO",
    "NI",
    "PB",
    "PT",
    "SR",
}
COORDINATION_METALS = {
    "AG",
    "AL",
    "AU",
    "BA",
    "CA",
    "CD",
    "CO",
    "CU",
    "FE",
    "HG",
    "K",
    "MG",
    "MN",
    "NA",
    "NI",
    "PB",
    "PT",
    "SR",
    "ZN",
}
METAL_DONOR_ELEMENTS = {"N", "O", "S"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
CLASH_CUTOFF_A = 1.8
PEPTIDE_BOND_MAX_A = 1.7
METAL_COORDINATION_MAX_A = 3.2
METAL_COORDINATION_MIN_A = {
    "CA": 1.95,
    "CD": 1.9,
    "K": 2.2,
    "MG": 1.75,
    "NA": 1.8,
    "PT": 1.5,
}
DEFAULT_METAL_COORDINATION_MIN_A = 1.7


@dataclass(frozen=True)
class GeometryAtom:
    """Heavy atom parsed from ATOM/HETATM PDB records."""

    line_no: int
    serial: str
    resname: str
    chain: str
    resseq: str
    icode: str
    record: str
    atom_name: str
    element: str
    xyz: tuple[float, float, float]


def audit_geometry(
    path: Path,
    *,
    sidecar_path: Path | None = None,
    clash_cutoff: float = CLASH_CUTOFF_A,
) -> dict[str, object]:
    """Audit heavy-atom contacts and optionally write exact clash-pair sidecar JSON."""

    atoms = _collect_atoms(path)
    if len(atoms) < 2:
        result = _empty_geometry_result(path, sidecar_path)
        _write_sidecar(sidecar_path, path, result, [])
        return result

    contacts = _collect_close_contacts(atoms, clash_cutoff=clash_cutoff)
    actionable = [contact for contact in contacts if not contact["ignored"]]
    min_distance = _min_contact_distance(actionable)
    ignored_counts = _ignored_reason_counts(contacts)
    result = {
        "geometry_status": "ok" if not actionable else "clashes",
        "geometry_clash_count": len(actionable),
        "geometry_min_nonbonded_distance_a": min_distance,
        "geometry_total_close_contact_count": len(contacts),
        "geometry_ignored_close_contact_count": len(contacts) - len(actionable),
        "geometry_expected_metal_coordination_count": ignored_counts.get(
            "expected_metal_coordination",
            0,
        ),
        "geometry_peptide_covalent_contact_count": ignored_counts.get(
            "peptide_covalent_contact",
            0,
        ),
        "geometry_metal_ligand_adjacent_contact_count": ignored_counts.get(
            "metal_ligand_adjacent_contact",
            0,
        ),
        "geometry_clash_audit_path": str(sidecar_path) if sidecar_path else "",
    }
    _write_sidecar(sidecar_path, path, result, contacts)
    return result


def apply_conservative_geometry_fixes(
    input_pdb: Path,
    output_pdb: Path,
    *,
    reference_pdb: Path | None = None,
    sidecar_path: Path | None = None,
    ligand_center: tuple[float, float, float] | None = None,
    binding_site_radius: float = 8.0,
    clash_cutoff: float = CLASH_CUTOFF_A,
) -> dict[str, object]:
    """Drop only clashing waters and newly added clashing sidechain atoms."""

    atoms = _collect_atoms(input_pdb)
    contacts = _collect_close_contacts(atoms, clash_cutoff=clash_cutoff)
    actionable = [contact for contact in contacts if not contact["ignored"]]
    water_residues = _clashing_water_residues(actionable)
    added_sidechain_atoms = _clashing_added_sidechain_atoms(
        actionable,
        reference_pdb=reference_pdb,
    )
    binding_site_sidechain_atoms = _binding_site_sidechain_atoms(
        actionable,
        added_sidechain_atoms,
        ligand_center=ligand_center,
        binding_site_radius=binding_site_radius,
    )
    if not water_residues and not added_sidechain_atoms:
        summary = {
            "geometry_fix_status": "unchanged",
            "geometry_fix_removed_water_count": 0,
            "geometry_fix_dropped_sidechain_atom_count": 0,
            "geometry_fix_binding_site_sidechain_count": 0,
            "geometry_fix_requires_constructive_repair": False,
            "geometry_fix_audit_path": str(sidecar_path) if sidecar_path else "",
            "geometry_fix_output_pdb": "",
        }
        _write_fix_sidecar(
            sidecar_path,
            input_pdb=input_pdb,
            output_pdb=None,
            summary=summary,
            actionable=actionable,
            water_residues=water_residues,
            sidechain_atoms=added_sidechain_atoms,
            binding_site_sidechain_atoms=binding_site_sidechain_atoms,
            binding_site_radius=binding_site_radius,
        )
        return summary

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    _write_geometry_fixed_pdb(
        input_pdb,
        output_pdb,
        water_residues=water_residues,
        sidechain_atoms=added_sidechain_atoms,
    )
    summary = {
        "geometry_fix_status": _format_fix_status(water_residues, added_sidechain_atoms),
        "geometry_fix_removed_water_count": len(water_residues),
        "geometry_fix_dropped_sidechain_atom_count": len(added_sidechain_atoms),
        "geometry_fix_binding_site_sidechain_count": len(binding_site_sidechain_atoms),
        "geometry_fix_requires_constructive_repair": bool(binding_site_sidechain_atoms),
        "geometry_fix_audit_path": str(sidecar_path) if sidecar_path else "",
        "geometry_fix_output_pdb": str(output_pdb),
    }
    _write_fix_sidecar(
        sidecar_path,
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        summary=summary,
        actionable=actionable,
        water_residues=water_residues,
        sidechain_atoms=added_sidechain_atoms,
        binding_site_sidechain_atoms=binding_site_sidechain_atoms,
        binding_site_radius=binding_site_radius,
    )
    return summary


def _empty_geometry_result(path: Path, sidecar_path: Path | None) -> dict[str, object]:
    return {
        "geometry_status": "no_atoms",
        "geometry_clash_count": 0,
        "geometry_min_nonbonded_distance_a": "",
        "geometry_total_close_contact_count": 0,
        "geometry_ignored_close_contact_count": 0,
        "geometry_expected_metal_coordination_count": 0,
        "geometry_peptide_covalent_contact_count": 0,
        "geometry_clash_audit_path": str(sidecar_path) if sidecar_path else "",
        "geometry_input_pdb": str(path),
    }


def _collect_atoms(path: Path) -> list[GeometryAtom]:
    atoms: list[GeometryAtom] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            element = _line_element(line)
            if element == "H":
                continue
            xyz = line_xyz(line)
            if xyz is None:
                continue
            atoms.append(
                GeometryAtom(
                    line_no=line_no,
                    serial=line[6:11].strip(),
                    resname=line[17:20].strip().upper(),
                    chain=line[21:22].strip() or "-",
                    resseq=line[22:26].strip(),
                    icode=line[26:27].strip(),
                    record=line[:6].strip(),
                    atom_name=line[12:16].strip(),
                    element=element,
                    xyz=xyz,
                )
            )
    return atoms


def _line_element(line: str) -> str:
    return (line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[:1]).upper()


def _collect_close_contacts(
    atoms: list[GeometryAtom],
    *,
    clash_cutoff: float,
) -> list[dict[str, object]]:
    atoms_by_residue: dict[tuple[str, str, str, str], list[GeometryAtom]] = defaultdict(list)
    for atom in atoms:
        atoms_by_residue[(atom.chain, atom.resseq, atom.icode, atom.resname)].append(atom)
    cell_size = float(clash_cutoff)
    cutoff2 = clash_cutoff * clash_cutoff
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    contacts: list[dict[str, object]] = []
    for index, atom in enumerate(atoms):
        cell = _grid_cell(atom, cell_size)
        for neighbor_cell in _neighbor_cells(cell):
            for other_index in grid.get(neighbor_cell, []):
                other = atoms[other_index]
                if _same_residue(atom, other):
                    continue
                dist2 = _distance2(atom.xyz, other.xyz)
                if dist2 >= cutoff2:
                    continue
                contacts.append(
                    _contact_record(atom, other, math.sqrt(dist2), atoms_by_residue)
                )
        grid[cell].append(index)
    return contacts


def _contact_record(
    a: GeometryAtom,
    b: GeometryAtom,
    distance: float,
    atoms_by_residue: Mapping[tuple[str, str, str, str], list[GeometryAtom]],
) -> dict[str, object]:
    reason = _ignored_reason(a, b, distance, atoms_by_residue)
    categories = {_atom_category(a), _atom_category(b)}
    return {
        "distance_a": round(distance, 3),
        "ignored": bool(reason),
        "ignored_reason": reason,
        "atom_a": _atom_record(a),
        "atom_b": _atom_record(b),
        "categories": sorted(categories),
        "involves_water": "water" in categories,
        "involves_metal": "metal" in categories,
        "involves_ligand_or_cofactor": "ligand/cofactor" in categories,
        "involves_protein_backbone": "protein_backbone" in categories,
        "involves_protein_sidechain": "protein_sidechain" in categories,
    }


def _atom_record(atom: GeometryAtom) -> dict[str, object]:
    return {
        "line_no": atom.line_no,
        "serial": atom.serial,
        "record": atom.record,
        "resname": atom.resname,
        "chain": atom.chain,
        "resseq": atom.resseq,
        "icode": atom.icode,
        "atom_name": atom.atom_name,
        "element": atom.element,
        "category": _atom_category(atom),
        "xyz": [round(value, 3) for value in atom.xyz],
    }


def _atom_category(atom: GeometryAtom) -> str:
    if atom.record == "HETATM" and atom.resname in WATER_NAMES:
        return "water"
    if atom.record == "HETATM" and _is_metal(atom):
        return "metal"
    if atom.record == "ATOM" and atom.atom_name.upper() in BACKBONE_ATOMS:
        return "protein_backbone"
    if atom.record == "ATOM":
        return "protein_sidechain"
    if atom.resname in COMMON_IONS:
        return "buffer/salt"
    return "ligand/cofactor"


def _ignored_reason(
    a: GeometryAtom,
    b: GeometryAtom,
    distance: float,
    atoms_by_residue: Mapping[tuple[str, str, str, str], list[GeometryAtom]],
) -> str:
    if _is_peptide_covalent_contact(a, b, distance):
        return "peptide_covalent_contact"
    if _is_expected_metal_coordination(a, b, distance):
        return "expected_metal_coordination"
    if _is_metal_ligand_adjacent_contact(a, b, distance, atoms_by_residue):
        return "metal_ligand_adjacent_contact"
    return ""


def _is_peptide_covalent_contact(
    a: GeometryAtom,
    b: GeometryAtom,
    distance: float,
) -> bool:
    if a.chain != b.chain:
        return False
    if distance > PEPTIDE_BOND_MAX_A:
        return False
    if {a.atom_name.upper(), b.atom_name.upper()} != {"C", "N"}:
        return False
    if a.record == "ATOM" and b.record == "ATOM":
        return True
    try:
        return abs(int(a.resseq) - int(b.resseq)) == 1
    except ValueError:
        return False


def _is_expected_metal_coordination(
    a: GeometryAtom,
    b: GeometryAtom,
    distance: float,
) -> bool:
    metal, donor = _metal_and_donor(a, b)
    if metal is None or donor is None:
        return False
    min_distance = METAL_COORDINATION_MIN_A.get(
        metal.element or metal.resname,
        DEFAULT_METAL_COORDINATION_MIN_A,
    )
    return min_distance <= distance <= METAL_COORDINATION_MAX_A


def _is_metal_ligand_adjacent_contact(
    a: GeometryAtom,
    b: GeometryAtom,
    distance: float,
    atoms_by_residue: Mapping[tuple[str, str, str, str], list[GeometryAtom]],
) -> bool:
    metal, neighbor = _metal_and_non_donor(a, b)
    if metal is None or neighbor is None or distance > CLASH_CUTOFF_A:
        return False
    if neighbor.element != "C":
        return False
    residue_key = (neighbor.chain, neighbor.resseq, neighbor.icode, neighbor.resname)
    for candidate in atoms_by_residue.get(residue_key, []):
        if candidate is neighbor or candidate.element not in METAL_DONOR_ELEMENTS:
            continue
        donor_metal_distance = math.sqrt(_distance2(candidate.xyz, metal.xyz))
        if donor_metal_distance > METAL_COORDINATION_MAX_A:
            continue
        donor_neighbor_distance = math.sqrt(_distance2(candidate.xyz, neighbor.xyz))
        if 0.9 <= donor_neighbor_distance <= 2.1:
            return True
    return False


def _metal_and_donor(
    a: GeometryAtom,
    b: GeometryAtom,
) -> tuple[GeometryAtom | None, GeometryAtom | None]:
    if _is_metal(a) and b.element in METAL_DONOR_ELEMENTS:
        return a, b
    if _is_metal(b) and a.element in METAL_DONOR_ELEMENTS:
        return b, a
    return None, None


def _metal_and_non_donor(
    a: GeometryAtom,
    b: GeometryAtom,
) -> tuple[GeometryAtom | None, GeometryAtom | None]:
    if _is_metal(a) and b.element not in METAL_DONOR_ELEMENTS:
        return a, b
    if _is_metal(b) and a.element not in METAL_DONOR_ELEMENTS:
        return b, a
    return None, None


def _is_metal(atom: GeometryAtom) -> bool:
    return atom.resname in COORDINATION_METALS or atom.element in COORDINATION_METALS


def _same_residue(a: GeometryAtom, b: GeometryAtom) -> bool:
    return (a.chain, a.resseq, a.icode, a.resname) == (
        b.chain,
        b.resseq,
        b.icode,
        b.resname,
    )


def _clashing_water_residues(
    contacts: list[dict[str, object]],
) -> set[tuple[str, str, str, str]]:
    residues: set[tuple[str, str, str, str]] = set()
    for contact in contacts:
        for atom in (contact["atom_a"], contact["atom_b"]):
            if isinstance(atom, dict) and atom.get("category") == "water":
                residues.add(_atom_residue_key(atom))
    return residues


def _clashing_added_sidechain_atoms(
    contacts: list[dict[str, object]],
    *,
    reference_pdb: Path | None,
) -> set[tuple[str, str, str, str, str, str]]:
    if reference_pdb is None or not reference_pdb.exists():
        return set()
    reference_atoms = _reference_atom_keys(reference_pdb)
    selected: set[tuple[str, str, str, str, str, str]] = set()
    for contact in contacts:
        raw_categories = contact.get("categories", [])
        categories = set(raw_categories) if isinstance(raw_categories, list) else set()
        if not categories <= {"protein_backbone", "protein_sidechain"}:
            continue
        for atom in (contact["atom_a"], contact["atom_b"]):
            if not isinstance(atom, dict) or atom.get("category") != "protein_sidechain":
                continue
            key = _atom_key(atom)
            if key not in reference_atoms:
                selected.add(key)
    return selected


def _binding_site_sidechain_atoms(
    contacts: list[dict[str, object]],
    sidechain_atoms: set[tuple[str, str, str, str, str, str]],
    *,
    ligand_center: tuple[float, float, float] | None,
    binding_site_radius: float,
) -> set[tuple[str, str, str, str, str, str]]:
    if ligand_center is None or not sidechain_atoms:
        return set()
    selected: set[tuple[str, str, str, str, str, str]] = set()
    radius2 = float(binding_site_radius) * float(binding_site_radius)
    for contact in contacts:
        for atom in (contact["atom_a"], contact["atom_b"]):
            if not isinstance(atom, dict) or atom.get("category") != "protein_sidechain":
                continue
            key = _atom_key(atom)
            if key in sidechain_atoms and _atom_distance2(atom, ligand_center) <= radius2:
                selected.add(key)
    return selected


def _atom_distance2(
    atom: Mapping[str, object],
    point: tuple[float, float, float],
) -> float:
    raw = atom.get("xyz")
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        return float("inf")
    try:
        x, y, z = raw
        return (float(x) - point[0]) ** 2 + (float(y) - point[1]) ** 2 + (
            float(z) - point[2]
        ) ** 2
    except Exception:
        return float("inf")


def _reference_atom_keys(path: Path) -> set[tuple[str, str, str, str, str, str]]:
    return {
        _geometry_atom_key(atom)
        for atom in _collect_atoms(path)
        if atom.record == "ATOM" and _atom_category(atom) == "protein_sidechain"
    }


def _write_geometry_fixed_pdb(
    input_pdb: Path,
    output_pdb: Path,
    *,
    water_residues: set[tuple[str, str, str, str]],
    sidechain_atoms: set[tuple[str, str, str, str, str, str]],
) -> None:
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith(("ATOM  ", "HETATM")) and _drop_line(
                    line,
                    water_residues=water_residues,
                    sidechain_atoms=sidechain_atoms,
                ):
                    continue
                dst.write(line)
    if not output_pdb.exists():
        shutil.copy2(input_pdb, output_pdb)


def _drop_line(
    line: str,
    *,
    water_residues: set[tuple[str, str, str, str]],
    sidechain_atoms: set[tuple[str, str, str, str, str, str]],
) -> bool:
    resname = line[17:20].strip().upper()
    atom: dict[str, object] = {
        "record": line[:6].strip(),
        "resname": resname,
        "chain": line[21:22].strip() or "-",
        "resseq": line[22:26].strip(),
        "icode": line[26:27].strip(),
        "atom_name": line[12:16].strip(),
    }
    if resname in WATER_NAMES and _atom_residue_key(atom) in water_residues:
        return True
    return _atom_key(atom) in sidechain_atoms


def _atom_residue_key(atom: Mapping[str, object]) -> tuple[str, str, str, str]:
    return (
        str(atom.get("resname", "")),
        str(atom.get("chain", "")),
        str(atom.get("resseq", "")),
        str(atom.get("icode", "")),
    )


def _atom_key(atom: Mapping[str, object]) -> tuple[str, str, str, str, str, str]:
    return (
        str(atom.get("record", "")),
        str(atom.get("resname", "")),
        str(atom.get("chain", "")),
        str(atom.get("resseq", "")),
        str(atom.get("icode", "")),
        str(atom.get("atom_name", "")),
    )


def _geometry_atom_key(atom: GeometryAtom) -> tuple[str, str, str, str, str, str]:
    return (
        atom.record,
        atom.resname,
        atom.chain,
        atom.resseq,
        atom.icode,
        atom.atom_name,
    )


def _format_fix_status(
    water_residues: set[tuple[str, str, str, str]],
    sidechain_atoms: set[tuple[str, str, str, str, str, str]],
) -> str:
    parts: list[str] = []
    if water_residues:
        parts.append(f"removed_clashing_waters:{len(water_residues)}")
    if sidechain_atoms:
        parts.append(f"dropped_added_sidechain_atoms:{len(sidechain_atoms)}")
    return ";".join(parts)


def _min_contact_distance(contacts: list[dict[str, object]]) -> float | str:
    if not contacts:
        return ""
    return min(float(str(contact["distance_a"])) for contact in contacts)


def _ignored_reason_counts(contacts: list[dict[str, object]]) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for contact in contacts:
        reason = str(contact.get("ignored_reason", ""))
        if reason:
            counts[reason] += 1
    return dict(counts)


def _write_sidecar(
    sidecar_path: Path | None,
    pdb_path: Path,
    summary: dict[str, object],
    contacts: list[dict[str, object]],
) -> None:
    if sidecar_path is None:
        return
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "pdb_path": str(pdb_path),
        "summary": summary,
        "clashes": contacts,
    }
    sidecar_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _write_fix_sidecar(
    sidecar_path: Path | None,
    *,
    input_pdb: Path,
    output_pdb: Path | None,
    summary: dict[str, object],
    actionable: list[dict[str, object]],
    water_residues: set[tuple[str, str, str, str]],
    sidechain_atoms: set[tuple[str, str, str, str, str, str]],
    binding_site_sidechain_atoms: set[tuple[str, str, str, str, str, str]],
    binding_site_radius: float,
) -> None:
    if sidecar_path is None:
        return
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "input_pdb": str(input_pdb),
        "output_pdb": str(output_pdb) if output_pdb else "",
        "summary": summary,
        "binding_site_radius_a": float(binding_site_radius),
        "removed_water_residues": _sorted_tuple_lists(water_residues),
        "dropped_sidechain_atoms": _sorted_tuple_lists(sidechain_atoms),
        "binding_site_sidechain_atoms": _sorted_tuple_lists(
            binding_site_sidechain_atoms
        ),
        "actionable_clashes": actionable,
    }
    sidecar_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _sorted_tuple_lists(items: Iterable[tuple[str, ...]]) -> list[list[str]]:
    return [list(item) for item in sorted(items)]


def _grid_cell(atom: GeometryAtom, cell_size: float) -> tuple[int, int, int]:
    return (
        math.floor(atom.xyz[0] / cell_size),
        math.floor(atom.xyz[1] / cell_size),
        math.floor(atom.xyz[2] / cell_size),
    )


def _neighbor_cells(cell: tuple[int, int, int]) -> list[tuple[int, int, int]]:
    return [
        (cell[0] + dx, cell[1] + dy, cell[2] + dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
    ]


def _distance2(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
