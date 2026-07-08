"""Binding-site water policy audit helpers for receptor preparation."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from protein_prep.pdb_records import line_xyz, residue_key

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
METAL_NAMES = {
    "AG",
    "AL",
    "BA",
    "CA",
    "CD",
    "CO",
    "CS",
    "CU",
    "FE",
    "HG",
    "K",
    "LI",
    "MG",
    "MN",
    "NA",
    "NI",
    "RB",
    "SR",
    "ZN",
}
POLAR_ELEMENTS = {"N", "O", "S"}


@dataclass(frozen=True)
class WaterRecord:
    key: tuple[str, str, str, str]
    atom_name: str
    xyz: tuple[float, float, float]
    occupancy: float | None
    bfactor: float | None


@dataclass(frozen=True)
class ContactAtom:
    resname: str
    chain: str
    resseq: str
    atom_name: str
    element: str
    xyz: tuple[float, float, float]


def audit_water_policy(
    input_pdb: Path,
    prepared_pdb: Path,
    ligand_pdb: Path,
    *,
    sidecar_path: Path | None = None,
    ligand_radius: float = 3.6,
    protein_radius: float = 3.6,
    metal_radius: float = 3.2,
    low_occupancy_cutoff: float = 0.5,
    high_bfactor_cutoff: float = 60.0,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Audit whether retained/removed waters have obvious structural support."""

    input_waters = _collect_waters(input_pdb)
    retained_keys = {water.key for water in _collect_waters(prepared_pdb)}
    ligand_atoms = _collect_contact_atoms(ligand_pdb, include_all_hets=True)
    protein_atoms = _collect_contact_atoms(input_pdb, include_protein_polar=True)
    metal_atoms = _collect_contact_atoms(input_pdb, include_metals=True)

    rows = [
        _water_row(
            water,
            retained=water.key in retained_keys,
            ligand_atoms=ligand_atoms,
            protein_atoms=protein_atoms,
            metal_atoms=metal_atoms,
            ligand_radius=ligand_radius,
            protein_radius=protein_radius,
            metal_radius=metal_radius,
            low_occupancy_cutoff=low_occupancy_cutoff,
            high_bfactor_cutoff=high_bfactor_cutoff,
        )
        for water in input_waters
    ]
    input_water_keys = {water.key for water in input_waters}
    summary = _summarize_water_rows(
        rows,
        input_water_count=len(input_waters),
        retained_water_count=len(input_water_keys & retained_keys),
        sidecar_path=sidecar_path,
    )
    _write_sidecar(sidecar_path, input_pdb, prepared_pdb, ligand_pdb, summary, rows)
    return summary, rows


def write_supported_water_receptor(
    input_pdb: Path,
    output_pdb: Path,
    rows: Iterable[Mapping[str, object]],
    *,
    source_water_pdb: Path | None = None,
) -> dict[str, object]:
    """Write a conservative receptor variant with only evidence-supported waters."""

    row_list = [row if isinstance(row, dict) else dict(row) for row in rows]
    supported = {
        _row_residue_key(row)
        for row in row_list
        if _supported_water_row(row)
    }
    _mark_selected_rows(row_list, supported)
    input_count = len(row_list)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    selected_count = _copy_with_water_selection(
        input_pdb,
        output_pdb,
        supported,
        source_water_pdb=source_water_pdb,
    )
    status = "no_input_waters" if input_count <= 0 else "selected"
    return {
        "water_policy_selection_status": status,
        "water_policy_selected_water_count": selected_count,
        "water_policy_removed_water_count": max(0, input_count - selected_count),
        "water_policy_selected_receptor": str(output_pdb),
    }


def _collect_waters(path: Path) -> list[WaterRecord]:
    waters: list[WaterRecord] = []
    if not path.exists():
        return waters
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper()
            if resname not in WATER_NAMES or _line_element(line) == "H":
                continue
            xyz = line_xyz(line)
            if xyz is None:
                continue
            waters.append(
                WaterRecord(
                    key=residue_key(line),
                    atom_name=line[12:16].strip(),
                    xyz=xyz,
                    occupancy=_safe_float(line[54:60]),
                    bfactor=_safe_float(line[60:66]),
                )
            )
    return waters


def _collect_contact_atoms(
    path: Path,
    *,
    include_all_hets: bool = False,
    include_protein_polar: bool = False,
    include_metals: bool = False,
) -> list[ContactAtom]:
    atoms: list[ContactAtom] = []
    if not path.exists():
        return atoms
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            element = _line_element(line)
            if element == "H":
                continue
            resname = line[17:20].strip().upper()
            record = line[:6].strip()
            if not _include_contact_atom(
                record=record,
                resname=resname,
                element=element,
                include_all_hets=include_all_hets,
                include_protein_polar=include_protein_polar,
                include_metals=include_metals,
            ):
                continue
            xyz = line_xyz(line)
            if xyz is None:
                continue
            atoms.append(
                ContactAtom(
                    resname=resname,
                    chain=line[21:22].strip() or "-",
                    resseq=line[22:26].strip(),
                    atom_name=line[12:16].strip(),
                    element=element,
                    xyz=xyz,
                )
            )
    return atoms


def _include_contact_atom(
    *,
    record: str,
    resname: str,
    element: str,
    include_all_hets: bool,
    include_protein_polar: bool,
    include_metals: bool,
) -> bool:
    if record == "HETATM" and resname in WATER_NAMES:
        return False
    if include_all_hets:
        return record == "HETATM"
    if include_protein_polar:
        return record == "ATOM" and element in POLAR_ELEMENTS
    if include_metals:
        return record == "HETATM" and resname in METAL_NAMES
    return False


def _water_row(
    water: WaterRecord,
    *,
    retained: bool,
    ligand_atoms: Iterable[ContactAtom],
    protein_atoms: Iterable[ContactAtom],
    metal_atoms: Iterable[ContactAtom],
    ligand_radius: float,
    protein_radius: float,
    metal_radius: float,
    low_occupancy_cutoff: float,
    high_bfactor_cutoff: float,
) -> dict[str, object]:
    nearest_ligand = _nearest_distance(water.xyz, ligand_atoms)
    nearest_protein = _nearest_distance(water.xyz, protein_atoms)
    nearest_metal = _nearest_distance(water.xyz, metal_atoms)
    near_ligand = _within(nearest_ligand, ligand_radius)
    near_protein = _within(nearest_protein, protein_radius)
    near_metal = _within(nearest_metal, metal_radius)
    bridging = (near_ligand and near_protein) or (near_ligand and near_metal) or (
        near_protein and near_metal
    )
    flags = _quality_flags(
        water,
        retained=retained,
        bridging=bridging,
        low_occupancy_cutoff=low_occupancy_cutoff,
        high_bfactor_cutoff=high_bfactor_cutoff,
    )
    return {
        "residue_key": list(water.key),
        "residue": _format_water_key(water.key),
        "retained": retained,
        "occupancy": _rounded_or_blank(water.occupancy),
        "bfactor": _rounded_or_blank(water.bfactor),
        "nearest_ligand_a": _rounded_or_blank(nearest_ligand),
        "nearest_protein_polar_a": _rounded_or_blank(nearest_protein),
        "nearest_metal_a": _rounded_or_blank(nearest_metal),
        "bridging_candidate": bridging,
        "near_ligand": near_ligand,
        "near_protein_polar": near_protein,
        "near_metal": near_metal,
        "quality_flags": flags,
    }


def _summarize_water_rows(
    rows: list[dict[str, object]],
    *,
    input_water_count: int,
    retained_water_count: int,
    sidecar_path: Path | None,
) -> dict[str, object]:
    removed_bridging = _count(rows, retained=False, bridging=True)
    retained_bridging = _count(rows, retained=True, bridging=True)
    retained_unjustified = sum(
        1
        for row in rows
        if bool(row.get("retained"))
        and not bool(row.get("bridging_candidate"))
        and not _flagged(row, "missing_quality_metadata")
    )
    retained_low_occupancy = _flag_count(rows, "retained_low_occupancy")
    retained_high_bfactor = _flag_count(rows, "retained_high_bfactor")
    review_count = (
        removed_bridging
        + retained_unjustified
        + retained_low_occupancy
        + retained_high_bfactor
    )
    status = _water_policy_status(input_water_count, review_count)
    return {
        "water_policy_status": status,
        "water_policy_input_water_count": input_water_count,
        "water_policy_retained_input_water_count": retained_water_count,
        "water_policy_retained_bridging_count": retained_bridging,
        "water_policy_removed_bridging_count": removed_bridging,
        "water_policy_retained_unjustified_count": retained_unjustified,
        "water_policy_low_occupancy_retained_count": retained_low_occupancy,
        "water_policy_high_bfactor_retained_count": retained_high_bfactor,
        "water_policy_review_count": review_count,
        "water_policy_audit_path": str(sidecar_path) if sidecar_path else "",
    }


def _water_policy_status(input_water_count: int, review_count: int) -> str:
    if input_water_count <= 0:
        return "no_input_waters"
    if review_count:
        return "review"
    return "ok"


def _quality_flags(
    water: WaterRecord,
    *,
    retained: bool,
    bridging: bool,
    low_occupancy_cutoff: float,
    high_bfactor_cutoff: float,
) -> list[str]:
    flags: list[str] = []
    _add_missing_quality_flag(flags, water)
    _add_occupancy_flags(
        flags,
        water,
        retained=retained,
        low_occupancy_cutoff=low_occupancy_cutoff,
    )
    _add_bfactor_flags(
        flags,
        water,
        retained=retained,
        high_bfactor_cutoff=high_bfactor_cutoff,
    )
    if bridging and not retained:
        flags.append("removed_bridging_candidate")
    return flags


def _add_missing_quality_flag(flags: list[str], water: WaterRecord) -> None:
    if water.occupancy is None or water.bfactor is None:
        flags.append("missing_quality_metadata")


def _add_occupancy_flags(
    flags: list[str],
    water: WaterRecord,
    *,
    retained: bool,
    low_occupancy_cutoff: float,
) -> None:
    if water.occupancy is not None and water.occupancy < low_occupancy_cutoff:
        flags.append("low_occupancy")
        if retained:
            flags.append("retained_low_occupancy")


def _add_bfactor_flags(
    flags: list[str],
    water: WaterRecord,
    *,
    retained: bool,
    high_bfactor_cutoff: float,
) -> None:
    if water.bfactor is not None and water.bfactor > high_bfactor_cutoff:
        flags.append("high_bfactor")
        if retained:
            flags.append("retained_high_bfactor")


def _count(rows: Iterable[Mapping[str, object]], *, retained: bool, bridging: bool) -> int:
    return sum(
        1
        for row in rows
        if bool(row.get("retained")) is retained
        and bool(row.get("bridging_candidate")) is bridging
    )


def _flag_count(rows: Iterable[Mapping[str, object]], flag: str) -> int:
    return sum(1 for row in rows if _flagged(row, flag))


def _flagged(row: Mapping[str, object], flag: str) -> bool:
    flags = row.get("quality_flags", [])
    return isinstance(flags, list) and flag in flags


def _supported_water_row(row: Mapping[str, object]) -> bool:
    if not _active_site_bridge(row):
        return False
    flags = row.get("quality_flags", [])
    if not isinstance(flags, list):
        return False
    blocked = {"missing_quality_metadata", "low_occupancy", "high_bfactor"}
    return not any(flag in flags for flag in blocked)


def _active_site_bridge(row: Mapping[str, object]) -> bool:
    return (
        bool(row.get("bridging_candidate"))
        and bool(row.get("near_ligand"))
        and (bool(row.get("near_protein_polar")) or bool(row.get("near_metal")))
    )


def _mark_selected_rows(
    rows: list[dict[str, object]],
    supported: set[tuple[str, ...]],
) -> None:
    for row in rows:
        selected = _row_residue_key(row) in supported
        row["water_selected_for_receptor"] = selected
        row["water_selection_reason"] = _water_selection_reason(row, selected)


def _water_selection_reason(row: Mapping[str, object], selected: bool) -> str:
    if selected:
        return "retained_active_site_bridge_quality_supported"
    if not bool(row.get("retained")):
        return "not_present_after_preparation"
    if not bool(row.get("near_ligand")):
        return "outside_ligand_site"
    if not bool(row.get("bridging_candidate")):
        return "not_bridging"
    flags = row.get("quality_flags", [])
    if isinstance(flags, list) and flags:
        return "quality_flags:" + ",".join(str(flag) for flag in flags)
    return "unsupported"


def _row_residue_key(row: Mapping[str, object]) -> tuple[str, ...]:
    raw = row.get("residue_key", [])
    if not isinstance(raw, list | tuple):
        return ()
    return tuple(str(part) for part in raw)


def _copy_with_water_selection(
    input_pdb: Path,
    output_pdb: Path,
    supported: set[tuple[str, ...]],
    *,
    source_water_pdb: Path | None = None,
) -> int:
    if source_water_pdb is not None and source_water_pdb != input_pdb:
        return _copy_receptor_with_source_waters(
            input_pdb,
            output_pdb,
            supported,
            source_water_pdb,
        )
    selected: set[tuple[str, ...]] = set()
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                if _is_water_line(line):
                    key = tuple(residue_key(line))
                    if key not in supported:
                        continue
                    selected.add(key)
                dst.write(line)
    return len(selected)


def _copy_receptor_with_source_waters(
    receptor_pdb: Path,
    output_pdb: Path,
    supported: set[tuple[str, ...]],
    source_water_pdb: Path,
) -> int:
    water_lines, selected = _selected_water_lines(source_water_pdb, supported)
    end_lines: list[str] = []
    with receptor_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                if _is_water_line(line):
                    continue
                if line.startswith("END"):
                    end_lines.append(line)
                    continue
                dst.write(line)
            dst.writelines(water_lines)
            dst.writelines(end_lines)
    return len(selected)


def _selected_water_lines(
    source_water_pdb: Path,
    supported: set[tuple[str, ...]],
) -> tuple[list[str], set[tuple[str, ...]]]:
    selected: set[tuple[str, ...]] = set()
    lines: list[str] = []
    if not supported or not source_water_pdb.exists():
        return lines, selected
    with source_water_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        for line in src:
            if not _is_water_line(line):
                continue
            key = tuple(residue_key(line))
            if key not in supported:
                continue
            lines.append(line)
            selected.add(key)
    return lines, selected


def _is_water_line(line: str) -> bool:
    return line.startswith("HETATM") and line[17:20].strip().upper() in WATER_NAMES


def _nearest_distance(
    xyz: tuple[float, float, float],
    atoms: Iterable[ContactAtom],
) -> float | None:
    best: float | None = None
    for atom in atoms:
        dist = math.sqrt(_distance2(xyz, atom.xyz))
        if best is None or dist < best:
            best = dist
    return best


def _within(distance: float | None, radius: float) -> bool:
    return distance is not None and distance <= radius


def _distance2(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _line_element(line: str) -> str:
    return (line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[:1]).upper()


def _safe_float(text: str) -> float | None:
    try:
        return float(text)
    except Exception:
        return None


def _rounded_or_blank(value: float | None) -> float | str:
    if value is None or not math.isfinite(value):
        return ""
    return round(value, 3)


def _format_water_key(key: tuple[str, str, str, str]) -> str:
    resname, chain, resseq, icode = key
    suffix = icode if icode else ""
    return f"{resname}:{chain}:{resseq}{suffix}"


def _write_sidecar(
    sidecar_path: Path | None,
    input_pdb: Path,
    prepared_pdb: Path,
    ligand_pdb: Path,
    summary: Mapping[str, object],
    rows: list[dict[str, object]],
) -> None:
    if sidecar_path is None:
        return
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "input_pdb": str(input_pdb),
        "prepared_pdb": str(prepared_pdb),
        "ligand_pdb": str(ligand_pdb),
        "summary": dict(summary),
        "waters": rows,
    }
    sidecar_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
