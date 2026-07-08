"""Meeko receptor-input normalization and conservative parser repair helpers."""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from protein_prep.geometry.repair import (
    atom_drop_key,
    invalid_oxt_atom_key_reasons,
    repair_terminal_oxt_geometry_in_pdb,
)
from protein_prep.pdb_records import line_resseq_int as _line_resseq_int
from protein_prep.pdb_records import line_xyz as _line_xyz

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
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
PROTEIN_RECORD_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
PEPTIDE_BOND_MAX_A = 1.7
MEEKO_FALSE_BOND_CUTOFF_A = 1.9
MEEKO_SAME_RESIDUE_FALSE_BOND_CUTOFF_A = 2.1
_CLEAN_MEEKO_PRUNE_REASONS = frozenset({"internal_oxt_on_peptide_link"})
_COMMON_RESIDUE_BONDS = {
    frozenset(pair)
    for pair in (
        ("N", "CA"),
        ("CA", "C"),
        ("C", "O"),
        ("C", "OXT"),
        ("CA", "CB"),
    )
}
_SIDECHAIN_RESIDUE_BONDS = {
    "ARG": (("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2")),
    "ASN": (("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")),
    "ASP": (("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")),
    "CYS": (("CB", "SG"),),
    "GLN": (("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")),
    "GLU": (("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")),
    "HIS": (("CB", "CG"), ("CG", "ND1"), ("ND1", "CE1"), ("CE1", "NE2"), ("NE2", "CD2"), ("CD2", "CG")),
    "ILE": (("CB", "CG1"), ("CG1", "CD1"), ("CB", "CG2")),
    "LEU": (("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")),
    "LYS": (("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")),
    "MET": (("CB", "CG"), ("CG", "SD"), ("SD", "CE")),
    "PHE": (("CB", "CG"), ("CG", "CD1"), ("CD1", "CE1"), ("CE1", "CZ"), ("CZ", "CE2"), ("CE2", "CD2"), ("CD2", "CG")),
    "PRO": (("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")),
    "SER": (("CB", "OG"),),
    "THR": (("CB", "OG1"), ("CB", "CG2")),
    "TRP": (("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "NE1"), ("NE1", "CE2"), ("CE2", "CD2"), ("CD2", "CE3"), ("CE3", "CZ3"), ("CZ3", "CH2"), ("CH2", "CZ2"), ("CZ2", "CE2")),
    "TYR": (("CB", "CG"), ("CG", "CD1"), ("CD1", "CE1"), ("CE1", "CZ"), ("CZ", "CE2"), ("CE2", "CD2"), ("CD2", "CG"), ("CZ", "OH")),
    "VAL": (("CB", "CG1"), ("CB", "CG2")),
}
_RESIDUE_BONDS = {
    resname: _COMMON_RESIDUE_BONDS | {frozenset(pair) for pair in pairs}
    for resname, pairs in _SIDECHAIN_RESIDUE_BONDS.items()
}

Coord = tuple[float, float, float]


@dataclass(frozen=True)
class MeekoAtomPoint:
    resname: str
    chain: str
    resseq: str
    record: str
    atom_name: str
    element: str
    xyz: Coord


def write_meeko_ordered_copy(input_pdb: Path, output_pdb: Path) -> None:
    atom_records: list[tuple[tuple[str, str, int, str, int], str]] = []
    other_lines: list[str] = []
    serial = 0
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                other_lines.append(line)
                continue
            serial += 1
            chain = line[21:22].strip()
            resseq = _line_resseq_int(line)
            key = (line[0:6], chain, resseq, line[26:27].strip(), serial)
            atom_records.append((key, line))

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w", encoding="utf-8") as handle:
        for _, line in sorted(atom_records, key=lambda item: item[0]):
            handle.write(line)
        for line in other_lines:
            if line.startswith(("CONECT", "END")):
                handle.write(line)
        if not other_lines or not any(line.startswith("END") for line in other_lines):
            handle.write("END\n")


def write_meeko_normalized_copy(
    input_pdb: Path,
    output_pdb: Path,
    *,
    drop_waters: bool = True,
    drop_hydrogens: bool = True,
    repair_terminal_oxt: bool = True,
) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                normalized = _meeko_normalized_line(
                    line,
                    drop_waters=drop_waters,
                    drop_hydrogens=drop_hydrogens,
                )
                if normalized is not None:
                    dst.write(normalized)
            dst.write("END\n")
    if repair_terminal_oxt:
        repair_terminal_oxt_geometry_in_pdb(
            output_pdb,
            audit_path=output_pdb.with_suffix(".oxt_repair_audit.json"),
        )


def meeko_pruned_stage_label(
    pruned_pdb: Path,
    *,
    dropped: int,
    allow_bad: bool = False,
) -> str:
    suffix = "allow_bad_retry" if allow_bad else "retry"
    audit_path = pruned_pdb.with_suffix(".drop_audit.json")
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        reasons = {
            str(row.get("reason", ""))
            for row in audit.get("dropped_atoms", [])
            if isinstance(row, Mapping)
        }
    except Exception:
        reasons = set()
    if reasons and reasons.issubset(_CLEAN_MEEKO_PRUNE_REASONS):
        return f"internal_oxt_pruned_{suffix}:dropped_atoms={dropped}"
    return f"clash_pruned_{suffix}:dropped_atoms={dropped}"


def write_meeko_clash_pruned_copy(input_pdb: Path, output_pdb: Path) -> int:
    normalized = output_pdb.with_suffix(".normalized.pdb")
    write_meeko_normalized_copy(input_pdb, normalized)
    drop_reasons = _meeko_false_bond_atom_key_reasons(normalized)
    drop_keys = set(drop_reasons)
    dropped_records: list[dict[str, str]] = []
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    dropped = 0
    with normalized.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("END"):
                    continue
                if line.startswith(("ATOM  ", "HETATM")):
                    key = (
                        line[:6].strip(),
                        line[21:22].strip() or "-",
                        line[22:26].strip(),
                        _normalize_terminal_resname(line[17:21].strip()),
                        line[12:16].strip(),
                    )
                    if key in drop_keys:
                        dropped += 1
                        dropped_records.append(
                            {
                                "record": key[0],
                                "chain": key[1],
                                "resseq": key[2],
                                "resname": key[3],
                                "atom_name": key[4],
                                "reason": drop_reasons.get(
                                    key,
                                    "meeko_false_bond_or_invalid_oxt",
                                ),
                            }
                        )
                        continue
                dst.write(line)
            dst.write("END\n")
    audit_path = output_pdb.with_suffix(".drop_audit.json")
    audit_path.write_text(
        json.dumps(
            {
                "input_pdb": str(input_pdb),
                "normalized_pdb": str(normalized),
                "output_pdb": str(output_pdb),
                "dropped_atom_count": dropped,
                "dropped_atoms": dropped_records,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return dropped


def _normalize_terminal_resname(resname: str) -> str:
    token = resname.strip().upper()
    if len(token) == 4 and token[0] in {"N", "C"} and token[1:] in PROTEIN_RECORD_RESNAMES:
        return token[1:]
    return token


def _rewrite_resname(line: str, resname: str) -> str:
    return f"{line[:17]}{resname.rjust(3)}{line[20:]}"


def _meeko_normalized_line(
    line: str,
    *,
    drop_waters: bool,
    drop_hydrogens: bool,
) -> str | None:
    if line.startswith("CONECT") or line.startswith("USER  MOD"):
        return None
    if not line.startswith(("ATOM  ", "HETATM")):
        return line if line.startswith("TER") else None
    resname = _normalize_terminal_resname(line[17:21].strip())
    if drop_waters and line.startswith("HETATM") and resname in WATER_NAMES:
        return None
    if drop_hydrogens and _line_element(line) in {"H", "D"}:
        return None
    return _rewrite_resname(line, resname)


def _meeko_false_bond_atom_key_reasons(
    path: Path,
) -> dict[tuple[str, str, str, str, str], str]:
    atoms = _collect_atom_points(path)
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    cutoff2 = MEEKO_FALSE_BOND_CUTOFF_A * MEEKO_FALSE_BOND_CUTOFF_A
    drop = invalid_oxt_atom_key_reasons(atoms)
    for index, atom in enumerate(atoms):
        cell = _grid_cell(atom, MEEKO_FALSE_BOND_CUTOFF_A)
        for neighbor_cell in _neighbor_cells(cell):
            for other_index in grid.get(neighbor_cell, []):
                other = atoms[other_index]
                dist2 = _distance2(atom.xyz, other.xyz)
                selected = _false_bond_drop_atom(atom, other, dist2, cutoff2)
                if selected is not None:
                    drop.setdefault(atom_drop_key(selected), "meeko_false_bond_or_clash")
        grid[cell].append(index)
    return drop


def _false_bond_drop_atom(
    atom: MeekoAtomPoint,
    other: MeekoAtomPoint,
    dist2: float,
    cutoff2: float,
) -> MeekoAtomPoint | None:
    same_residue = (atom.chain, atom.resseq, atom.resname) == (
        other.chain,
        other.resseq,
        other.resname,
    )
    if same_residue:
        same_residue_cutoff2 = (
            MEEKO_SAME_RESIDUE_FALSE_BOND_CUTOFF_A
            * MEEKO_SAME_RESIDUE_FALSE_BOND_CUTOFF_A
        )
        if dist2 >= same_residue_cutoff2:
            return None
        if _is_expected_same_residue_bond(atom, other):
            return None
        return _choose_same_residue_false_bond_drop(atom, other)
    if dist2 >= cutoff2:
        return None
    if _is_valid_peptide_contact(atom, other, dist2):
        return None
    return _choose_meeko_clash_drop(atom, other)


def _is_valid_peptide_contact(
    a: MeekoAtomPoint,
    b: MeekoAtomPoint,
    dist2: float,
) -> bool:
    if a.record != "ATOM" or b.record != "ATOM" or a.chain != b.chain:
        return False
    try:
        adjacent = abs(int(a.resseq) - int(b.resseq)) == 1
    except Exception:
        return False
    if not adjacent or dist2 > PEPTIDE_BOND_MAX_A * PEPTIDE_BOND_MAX_A:
        return False
    names = {a.atom_name.upper(), b.atom_name.upper()}
    return names == {"C", "N"}


def _is_expected_same_residue_bond(a: MeekoAtomPoint, b: MeekoAtomPoint) -> bool:
    if a.record != "ATOM" or b.record != "ATOM":
        return False
    if (a.chain, a.resseq, a.resname) != (b.chain, b.resseq, b.resname):
        return False
    names = frozenset((a.atom_name.upper(), b.atom_name.upper()))
    if "OXT" in names:
        return names == frozenset(("C", "OXT"))
    return names in _RESIDUE_BONDS.get(a.resname, _COMMON_RESIDUE_BONDS)


def _choose_same_residue_false_bond_drop(
    a: MeekoAtomPoint,
    b: MeekoAtomPoint,
) -> MeekoAtomPoint | None:
    if a.atom_name.upper() == "OXT":
        return a
    if b.atom_name.upper() == "OXT":
        return b
    return _choose_meeko_clash_drop(a, b)


def _choose_meeko_clash_drop(a: MeekoAtomPoint, b: MeekoAtomPoint) -> MeekoAtomPoint | None:
    a_priority = _atom_drop_priority(a)
    b_priority = _atom_drop_priority(b)
    if max(a_priority, b_priority) <= 0:
        return None
    if a_priority == b_priority:
        return b if b.record == "ATOM" else a
    return a if a_priority > b_priority else b


def _atom_drop_priority(atom: MeekoAtomPoint) -> int:
    if atom.record == "HETATM" and atom.resname in WATER_NAMES:
        return 50
    if atom.record == "HETATM" and atom.resname not in COORDINATION_METALS:
        return 40
    if atom.record == "ATOM" and atom.atom_name.upper() not in BACKBONE_ATOMS:
        return 30
    if atom.record == "HETATM":
        return 20
    return 0


def _collect_atom_points(path: Path) -> list[MeekoAtomPoint]:
    atoms: list[MeekoAtomPoint] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or not _is_heavy_atom(line):
                continue
            xyz = _line_xyz(line)
            if xyz is None:
                continue
            atoms.append(
                MeekoAtomPoint(
                    resname=_normalize_terminal_resname(line[17:21].strip()),
                    chain=line[21:22].strip() or "-",
                    resseq=line[22:26].strip(),
                    record=line[:6].strip(),
                    atom_name=line[12:16].strip(),
                    element=_line_element(line),
                    xyz=xyz,
                )
            )
    return atoms


def _is_heavy_atom(line: str) -> bool:
    return _line_element(line) not in {"H", "D"}


def _line_element(line: str) -> str:
    return (line[76:78].strip() or line[12:16].strip()[:1]).upper()


def _grid_cell(atom: MeekoAtomPoint, cell_size: float) -> tuple[int, int, int]:
    return (
        int(atom.xyz[0] // cell_size),
        int(atom.xyz[1] // cell_size),
        int(atom.xyz[2] // cell_size),
    )


def _neighbor_cells(cell: tuple[int, int, int]) -> list[tuple[int, int, int]]:
    x, y, z = cell
    return [
        (x + dx, y + dy, z + dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
    ]


def _distance2(a: Coord, b: Coord) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


__all__ = [
    "meeko_pruned_stage_label",
    "write_meeko_clash_pruned_copy",
    "write_meeko_normalized_copy",
    "write_meeko_ordered_copy",
]
