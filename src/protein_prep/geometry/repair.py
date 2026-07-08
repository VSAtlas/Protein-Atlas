"""Small geometry repairs for clean Meeko receptor parsing."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Protocol, Sequence

Coord = tuple[float, float, float]
ResidueKey = tuple[str, str, str, str]
ResidueAtom = tuple[int, str, Coord]
ResidueAtoms = dict[str, ResidueAtom]


class AtomPointLike(Protocol):
    @property
    def record(self) -> str: ...

    @property
    def chain(self) -> str: ...

    @property
    def resseq(self) -> str: ...

    @property
    def resname(self) -> str: ...

    @property
    def atom_name(self) -> str: ...

    @property
    def xyz(self) -> Coord: ...


def atom_drop_key(atom: AtomPointLike) -> tuple[str, str, str, str, str]:
    return (atom.record, atom.chain, atom.resseq, atom.resname, atom.atom_name)


def _distance2(
    a: Coord,
    b: Coord,
) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _distance(
    a: Coord,
    b: Coord,
) -> float:
    return math.sqrt(_distance2(a, b))


def _sub(
    a: Coord,
    b: Coord,
) -> Coord:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(
    a: Coord,
    b: Coord,
) -> Coord:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(
    a: Coord,
    factor: float,
) -> Coord:
    return (a[0] * factor, a[1] * factor, a[2] * factor)


def _dot(
    a: Coord,
    b: Coord,
) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(
    a: Coord,
    b: Coord,
) -> Coord:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Coord) -> float:
    return math.sqrt(_dot(a, a))


def _unit(
    a: Coord,
) -> Coord | None:
    norm = _norm(a)
    if norm < 1.0e-8:
        return None
    return _scale(a, 1.0 / norm)


def _line_xyz(line: str) -> Coord | None:
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except Exception:
        return None


def _replace_line_xyz(
    line: str,
    xyz: Coord,
) -> str:
    return f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"


def _residue_key(line: str) -> ResidueKey:
    return (
        line[21:22].strip() or "-",
        line[22:26].strip(),
        line[26:27].strip(),
        line[17:21].strip().upper(),
    )


def _candidate_terminal_oxt_positions(
    *,
    c_xyz: Coord,
    ca_xyz: Coord,
    o_xyz: Coord,
) -> list[Coord]:
    c_o = _unit(_sub(o_xyz, c_xyz))
    c_ca = _unit(_sub(ca_xyz, c_xyz))
    if c_o is None or c_ca is None:
        return []
    away_from_ca = _scale(c_ca, -1.0)
    in_plane = _sub(away_from_ca, _scale(c_o, _dot(away_from_ca, c_o)))
    in_plane_u = _unit(in_plane)
    if in_plane_u is None:
        fallback = _cross(c_o, (1.0, 0.0, 0.0))
        in_plane_u = _unit(fallback) or _unit(_cross(c_o, (0.0, 1.0, 0.0)))
    if in_plane_u is None:
        return []

    target_len = 1.25
    cos_angle = math.cos(math.radians(120.0))
    sin_angle = math.sin(math.radians(120.0))
    candidates: list[Coord] = []
    for sign in (1.0, -1.0):
        direction = _add(_scale(c_o, cos_angle), _scale(in_plane_u, sign * sin_angle))
        direction_u = _unit(direction)
        if direction_u is not None:
            candidates.append(_add(c_xyz, _scale(direction_u, target_len)))
    return candidates


def _oxt_geometry_bad(
    residue_atoms: ResidueAtoms,
) -> tuple[bool, str]:
    required = {"C", "CA", "O", "OXT"}
    if not required.issubset(residue_atoms):
        return False, "missing_terminal_reference_atoms"
    oxt_xyz = residue_atoms["OXT"][2]
    c_xyz = residue_atoms["C"][2]
    c_oxt = _distance(c_xyz, oxt_xyz)
    if not 1.05 <= c_oxt <= 1.45:
        return True, f"bad_c_oxt_distance:{c_oxt:.3f}"
    close_contacts = [
        atom_name
        for atom_name, (_, _, xyz) in residue_atoms.items()
        if atom_name not in {"C", "O", "OXT"} and _distance(oxt_xyz, xyz) < 1.70
    ]
    if close_contacts:
        return True, "same_residue_oxt_close_contact:" + ",".join(sorted(close_contacts))
    return False, "ok"


def _score_oxt_candidate(
    candidate: Coord,
    residue_atoms: ResidueAtoms,
) -> float:
    nonbonded_distances = [
        _distance(candidate, xyz)
        for atom_name, (_, _, xyz) in residue_atoms.items()
        if atom_name not in {"C", "O", "OXT"}
    ]
    min_nonbonded = min(nonbonded_distances, default=9.0)
    oxygen_distance = _distance(candidate, residue_atoms["O"][2])
    return min_nonbonded - 0.15 * abs(oxygen_distance - 2.15)


def _collect_oxt_reference_atoms(lines: Sequence[str]) -> dict[ResidueKey, ResidueAtoms]:
    residues: dict[ResidueKey, ResidueAtoms] = defaultdict(dict)
    for index, line in enumerate(lines):
        if not line.startswith("ATOM  "):
            continue
        atom_name = line[12:16].strip().upper()
        if atom_name not in {"C", "CA", "O", "OXT"}:
            continue
        xyz = _line_xyz(line)
        if xyz is not None:
            residues[_residue_key(line)][atom_name] = (index, line, xyz)
    return residues


def _oxt_audit_record(key: ResidueKey, reason: str) -> dict[str, object]:
    return {
        "chain": key[0],
        "resseq": key[1],
        "icode": key[2],
        "resname": key[3],
        "reason": reason,
    }


def _repair_oxt_residue(
    key: ResidueKey,
    residue_atoms: ResidueAtoms,
    lines: list[str],
) -> tuple[dict[str, object] | None, dict[str, object] | None]:
    if "OXT" not in residue_atoms:
        return None, None
    bad, reason = _oxt_geometry_bad(residue_atoms)
    if not bad:
        return None, None
    if not {"C", "CA", "O", "OXT"}.issubset(residue_atoms):
        return None, _oxt_audit_record(key, reason)

    c_xyz = residue_atoms["C"][2]
    old_xyz = residue_atoms["OXT"][2]
    candidates = _candidate_terminal_oxt_positions(
        c_xyz=c_xyz,
        ca_xyz=residue_atoms["CA"][2],
        o_xyz=residue_atoms["O"][2],
    )
    if not candidates:
        return None, _oxt_audit_record(key, "no_repair_candidate")

    new_xyz = max(candidates, key=lambda item: _score_oxt_candidate(item, residue_atoms))
    oxt_index, oxt_line, _ = residue_atoms["OXT"]
    lines[oxt_index] = _replace_line_xyz(oxt_line, new_xyz)
    record = _oxt_audit_record(key, reason)
    record.update(
        {
            "old_xyz": [round(value, 3) for value in old_xyz],
            "new_xyz": [round(value, 3) for value in new_xyz],
            "old_c_oxt_a": round(_distance(c_xyz, old_xyz), 3),
            "new_c_oxt_a": round(_distance(c_xyz, new_xyz), 3),
        }
    )
    return record, None


def repair_terminal_oxt_geometry_in_pdb(
    pdb_path: Path,
    *,
    audit_path: Path | None = None,
) -> dict[str, object]:
    """Repair terminal OXT positions that would otherwise be dropped for Meeko."""

    lines = pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines(
        keepends=True
    )
    repairs: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []
    for key, residue_atoms in _collect_oxt_reference_atoms(lines).items():
        repair, skip = _repair_oxt_residue(key, residue_atoms, lines)
        if repair is not None:
            repairs.append(repair)
        if skip is not None:
            skipped.append(skip)

    if repairs:
        pdb_path.write_text("".join(lines), encoding="utf-8")

    summary: dict[str, object] = {
        "pdb_path": str(pdb_path),
        "repaired_oxt_count": len(repairs),
        "skipped_oxt_count": len(skipped),
        "repairs": repairs,
        "skipped": skipped,
    }
    if audit_path is not None:
        audit_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    return summary


def invalid_oxt_atom_keys(
    atoms: Sequence[AtomPointLike],
) -> set[tuple[str, str, str, str, str]]:
    """Return terminal OXT atoms whose geometry would create false RDKit bonds."""

    return set(invalid_oxt_atom_key_reasons(atoms))


def invalid_oxt_atom_key_reasons(
    atoms: Sequence[AtomPointLike],
) -> dict[tuple[str, str, str, str, str], str]:
    """Return invalid OXT atom keys with chemistry-preserving prune reasons."""

    all_atoms = list(atoms)
    residues: dict[tuple[str, str, str], list[AtomPointLike]] = defaultdict(list)
    for atom in all_atoms:
        residues[(atom.chain, atom.resseq, atom.resname)].append(atom)
    drop: dict[tuple[str, str, str, str, str], str] = {}
    for residue_atoms in residues.values():
        drop.update(_invalid_oxt_atom_key_reasons_for_residue(residue_atoms, all_atoms))
    return drop


def _invalid_oxt_atom_key_reasons_for_residue(
    residue_atoms: Sequence[AtomPointLike],
    all_atoms: Sequence[AtomPointLike],
) -> dict[tuple[str, str, str, str, str], str]:
    oxt_atoms = [atom for atom in residue_atoms if atom.atom_name.upper() == "OXT"]
    c_atoms = [atom for atom in residue_atoms if atom.atom_name.upper() == "C"]
    drop: dict[tuple[str, str, str, str, str], str] = {}
    for oxt in oxt_atoms:
        internal_reason = _internal_oxt_reason(oxt, c_atoms, all_atoms)
        if internal_reason:
            drop[atom_drop_key(oxt)] = internal_reason
        elif _invalid_oxt_geometry(oxt, c_atoms, residue_atoms):
            drop[atom_drop_key(oxt)] = "invalid_oxt_geometry"
    return drop


def _invalid_oxt_geometry(
    oxt: AtomPointLike,
    c_atoms: Sequence[AtomPointLike],
    residue_atoms: Sequence[AtomPointLike],
) -> bool:
    c_oxt = min((_distance2(oxt.xyz, atom.xyz) for atom in c_atoms), default=999.0)
    if not (1.05 * 1.05 <= c_oxt <= 1.45 * 1.45):
        return True
    return any(
        atom.atom_name.upper() not in {"C", "OXT"}
        and _distance2(oxt.xyz, atom.xyz) < 1.7 * 1.7
        for atom in residue_atoms
    )


def _internal_oxt_reason(
    oxt: AtomPointLike,
    c_atoms: Sequence[AtomPointLike],
    all_atoms: Sequence[AtomPointLike],
) -> str:
    residue_id = (oxt.chain, oxt.resseq, oxt.resname)
    peptide_min2 = 1.05 * 1.05
    peptide_max2 = 1.70 * 1.70
    for atom in all_atoms:
        if (atom.chain, atom.resseq, atom.resname) == residue_id:
            continue
        if atom.record != "ATOM" or atom.atom_name.upper() != "N":
            continue
        if any(
            peptide_min2 <= _distance2(c_atom.xyz, atom.xyz) <= peptide_max2
            for c_atom in c_atoms
        ):
            return "internal_oxt_on_peptide_link"
    return ""
