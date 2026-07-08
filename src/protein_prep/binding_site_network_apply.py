"""Apply conservative binding-site network recommendations to PDB receptors."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from protein_prep.pdb_records import line_xyz, residue_key

Coord = tuple[float, float, float]
ResidueKey = tuple[str, str, str, str]

POLAR_ELEMENTS = {"N", "O", "S"}
WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
APPLY_RADIUS_A = 4.5
WATER_H_BOND_A = 0.957
WATER_H_ANGLE_DEG = 104.5


def apply_binding_site_network_recommendations(
    input_pdb: Path,
    output_pdb: Path,
    *,
    selection_sidecar: Path,
    ligand_pdb: Path | None = None,
    water_policy_rows: Sequence[Mapping[str, object]] = (),
    apply_protein_modes: bool = True,
    apply_water_hydrogens: bool = False,
    apply_low_confidence: bool = False,
    audit_path: Path | None = None,
) -> dict[str, object]:
    """Apply only low-risk local network recommendations.

    Protein-side recommendations are intentionally conservative.  We do not
    override explicit PDB2PQR/Reduce hydrogens because that can silently create
    inconsistent residue names and hydrogen coordinates.  Water hydrogens are
    added only to waters selected by the water-policy audit.
    """

    if not input_pdb.exists():
        return _write_summary(
            audit_path,
            "missing_input",
            input_pdb,
            output_pdb,
            [],
            low_confidence_skipped=0,
        )
    records = _collect_records(input_pdb)
    if not selection_sidecar.exists():
        _write_lines(output_pdb, [record.line for record in records])
        return _write_summary(
            audit_path,
            "missing_selection_sidecar",
            input_pdb,
            output_pdb,
            [],
            low_confidence_skipped=0,
        )

    payload = _load_json(selection_sidecar)
    recommendations, low_confidence_skipped = _collect_recommendations(
        payload,
        apply_low_confidence=apply_low_confidence,
    )
    ligand_polar = _collect_polar_points(ligand_pdb) if ligand_pdb else []
    residue_records = _records_by_residue(records)
    actions: list[dict[str, object]] = []
    line_rewrites: dict[int, str] = {}

    if apply_protein_modes:
        actions.extend(
            _plan_protein_mode_rewrites(
                residue_records,
                recommendations,
                ligand_polar,
                line_rewrites,
            )
        )

    selected_waters = _selected_water_keys(water_policy_rows)
    if apply_water_hydrogens and selected_waters:
        actions.extend(
            _plan_selected_water_hydrogens(
                records,
                residue_records,
                selected_waters,
                ligand_polar,
                line_rewrites,
            )
        )

    output_lines = [line_rewrites.get(record.index, record.line) for record in records]
    _write_lines(output_pdb, output_lines)
    status = "applied" if _applied_count(actions) else "no_changes"
    return _write_summary(
        audit_path,
        status,
        input_pdb,
        output_pdb,
        actions,
        low_confidence_skipped=low_confidence_skipped,
    )


class _Record:
    def __init__(
        self,
        *,
        index: int,
        line: str,
        key: ResidueKey | None = None,
        atom_name: str = "",
        element: str = "",
        xyz: Coord | None = None,
    ) -> None:
        self.index = index
        self.line = line
        self.key = key
        self.atom_name = atom_name
        self.element = element
        self.xyz = xyz


def _collect_records(path: Path) -> list[_Record]:
    records: list[_Record] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for index, line in enumerate(handle):
            if not line.startswith(("ATOM  ", "HETATM")):
                records.append(_Record(index=index, line=line))
                continue
            atom_name = line[12:16].strip().upper()
            records.append(
                _Record(
                    index=index,
                    line=line,
                    key=residue_key(line),
                    atom_name=atom_name,
                    element=_line_element(line),
                    xyz=line_xyz(line),
                )
            )
    return records


def _records_by_residue(
    records: Sequence[_Record],
) -> dict[ResidueKey, list[_Record]]:
    grouped: dict[ResidueKey, list[_Record]] = defaultdict(list)
    for record in records:
        if record.key is not None:
            grouped[record.key].append(record)
    return grouped


def _collect_recommendations(
    payload: Mapping[str, object],
    *,
    apply_low_confidence: bool,
) -> tuple[dict[ResidueKey, dict[str, object]], int]:
    recommendations: dict[ResidueKey, dict[str, object]] = {}
    low_confidence_skipped = 0
    selections = payload.get("selections", [])
    if not isinstance(selections, list):
        return recommendations, 0
    for selection in selections:
        rows, skipped = _recommendations_for_selection(
            selection,
            apply_low_confidence=apply_low_confidence,
        )
        low_confidence_skipped += skipped
        for key, mode in rows:
            recommendations[key] = mode
    return recommendations, low_confidence_skipped


def _recommendations_for_selection(
    selection: object,
    *,
    apply_low_confidence: bool,
) -> tuple[list[tuple[ResidueKey, dict[str, object]]], int]:
    if not isinstance(selection, Mapping):
        return [], 0
    if str(selection.get("status", "")) != "selected":
        return [], 0
    if str(selection.get("confidence", "")).lower() == "low" and not apply_low_confidence:
        return [], 1
    modes = _recommended_site_modes(selection)
    rows: list[tuple[ResidueKey, dict[str, object]]] = []
    for mode in modes:
        key = _parse_residue_key(str(mode.get("residue", "")))
        if key is not None:
            rows.append((key, dict(mode)))
    return rows, 0


def _recommended_site_modes(selection: Mapping[str, object]) -> list[Mapping[str, object]]:
    network = selection.get("network_optimizer", {})
    if not isinstance(network, Mapping):
        return []
    modes = network.get("recommended_site_modes", [])
    if not isinstance(modes, list):
        return []
    return [
        mode
        for mode in modes
        if isinstance(mode, Mapping) and _eligible_site_mode(mode)
    ]


def _eligible_site_mode(mode: Mapping[str, object]) -> bool:
    try:
        distance = float(str(mode.get("distance_to_het_a", "999")))
    except ValueError:
        distance = 999.0
    try:
        contact_count = int(str(mode.get("network_contact_count", "0")))
    except ValueError:
        contact_count = 0
    return distance <= APPLY_RADIUS_A or contact_count > 0


def _parse_residue_key(text: str) -> ResidueKey | None:
    parts = text.split(":")
    if len(parts) != 3:
        return None
    resname, chain, resseq = parts
    return (resname.upper(), "" if chain == "-" else chain, resseq, "")


def _plan_protein_mode_rewrites(
    residue_records: Mapping[ResidueKey, Sequence[_Record]],
    recommendations: Mapping[ResidueKey, Mapping[str, object]],
    ligand_polar: Sequence[Coord],
    line_rewrites: dict[int, str],
) -> list[dict[str, object]]:
    actions: list[dict[str, object]] = []
    for key, recommendation in recommendations.items():
        records = residue_records.get(key, ())
        if not records:
            continue
        site_type = str(recommendation.get("site_type", ""))
        if site_type == "histidine_tautomer_protomer":
            actions.append(
                _apply_histidine_mode(records, ligand_polar, recommendation, line_rewrites)
            )
        elif site_type == "sidechain_flip":
            actions.append(
                _apply_sidechain_flip(records, ligand_polar, recommendation, line_rewrites)
            )
        elif site_type == "borderline_carboxylate_protonation":
            actions.append(
                _skipped_action(
                    records,
                    recommendation,
                    "carboxylate_protonation_left_to_pdb2pqr_propka",
                )
            )
    return actions


def _apply_histidine_mode(
    records: Sequence[_Record],
    ligand_polar: Sequence[Coord],
    recommendation: Mapping[str, object],
    line_rewrites: dict[int, str],
) -> dict[str, object]:
    current = records[0].key[0] if records[0].key else ""
    recommended = _choose_histidine_state(records, ligand_polar)
    if not ligand_polar:
        return _skipped_action(
            records,
            recommendation,
            "missing_ligand_polar_context",
            recommended=recommended,
        )
    if current in {"HID", "HIE", "HIP"}:
        return _skipped_action(
            records,
            recommendation,
            f"already_typed_by_pdb2pqr_reduce:{current}",
            recommended=recommended,
        )
    if _has_any_sidechain_hydrogen(records):
        return _skipped_action(
            records,
            recommendation,
            "explicit_histidine_hydrogens_present",
            recommended=recommended,
        )
    for record in records:
        line_rewrites[record.index] = _replace_resname(record.line, recommended)
    return _applied_action(records, recommendation, "histidine_renamed", recommended)


def _choose_histidine_state(records: Sequence[_Record], ligand_polar: Sequence[Coord]) -> str:
    nd1 = _atom_xyz(records, "ND1")
    ne2 = _atom_xyz(records, "NE2")
    if nd1 is None and ne2 is None:
        return "HIE"
    nd1_dist = _min_distance(nd1, ligand_polar) if nd1 else math.inf
    ne2_dist = _min_distance(ne2, ligand_polar) if ne2 else math.inf
    if nd1_dist + 0.25 < ne2_dist:
        return "HID"
    if ne2_dist + 0.25 < nd1_dist:
        return "HIE"
    return "HIE"


def _apply_sidechain_flip(
    records: Sequence[_Record],
    ligand_polar: Sequence[Coord],
    recommendation: Mapping[str, object],
    line_rewrites: dict[int, str],
) -> dict[str, object]:
    if _has_any_sidechain_hydrogen(records):
        return _skipped_action(
            records,
            recommendation,
            "explicit_sidechain_hydrogens_present_reduce_authoritative",
        )
    resname = records[0].key[0] if records[0].key else ""
    pair = ("OD1", "ND2") if resname == "ASN" else ("OE1", "NE2") if resname == "GLN" else ()
    if not pair or not _amide_flip_recommended(records, ligand_polar, pair):
        return _skipped_action(records, recommendation, "current_amide_orientation_preferred")
    for record in records:
        if record.atom_name == pair[0]:
            line_rewrites[record.index] = _replace_atom_name_and_element(record.line, pair[1], "N")
        elif record.atom_name == pair[1]:
            line_rewrites[record.index] = _replace_atom_name_and_element(record.line, pair[0], "O")
    return _applied_action(records, recommendation, "amide_atom_names_flipped", ",".join(pair))


def _amide_flip_recommended(
    records: Sequence[_Record],
    ligand_polar: Sequence[Coord],
    pair: tuple[str, str],
) -> bool:
    if not ligand_polar:
        return False
    oxygen = _atom_xyz(records, pair[0])
    nitrogen = _atom_xyz(records, pair[1])
    if oxygen is None or nitrogen is None:
        return False
    oxygen_dist = _min_distance(oxygen, ligand_polar)
    nitrogen_dist = _min_distance(nitrogen, ligand_polar)
    return oxygen_dist + 0.40 < nitrogen_dist and oxygen_dist <= 3.2


def _plan_selected_water_hydrogens(
    records: Sequence[_Record],
    residue_records: Mapping[ResidueKey, Sequence[_Record]],
    selected_waters: set[ResidueKey],
    ligand_polar: Sequence[Coord],
    line_rewrites: dict[int, str],
) -> list[dict[str, object]]:
    actions: list[dict[str, object]] = []
    serial = _max_serial(records)
    protein_polar = _protein_polar_points(records)
    for key in sorted(selected_waters):
        water_records = residue_records.get(key, ())
        oxygen = _water_oxygen_record(water_records)
        if oxygen is None:
            continue
        if any(record.element == "H" for record in water_records):
            actions.append(_water_action(key, "skipped", "water_hydrogens_already_present"))
            continue
        oxygen_xyz = oxygen.xyz or (0.0, 0.0, 0.0)
        h1, h2 = _water_hydrogen_positions(
            oxygen_xyz,
            ligand_polar=ligand_polar,
            protein_polar=protein_polar,
        )
        serial += 1
        h1_line = _format_water_hydrogen_line(oxygen.line, serial, "H1", h1)
        serial += 1
        h2_line = _format_water_hydrogen_line(oxygen.line, serial, "H2", h2)
        line_rewrites[oxygen.index] = oxygen.line + h1_line + h2_line
        actions.append(_water_action(key, "applied", "water_hydrogens_added"))
    return actions


def _protein_polar_points(records: Sequence[_Record]) -> list[Coord]:
    return [
        record.xyz
        for record in records
        if record.xyz is not None
        and record.key is not None
        and record.key[0] not in WATER_NAMES
        and record.element in POLAR_ELEMENTS
        and record.line.startswith("ATOM")
    ]


def _water_oxygen_record(water_records: Sequence[_Record]) -> _Record | None:
    for record in water_records:
        if record.element == "O" and record.xyz is not None:
            return record
    return None


def _selected_water_keys(rows: Sequence[Mapping[str, object]]) -> set[ResidueKey]:
    keys: set[ResidueKey] = set()
    for row in rows:
        if not bool(row.get("water_selected_for_receptor")):
            continue
        raw = row.get("residue_key", [])
        if isinstance(raw, list | tuple) and len(raw) == 4:
            keys.add(tuple(str(part) for part in raw))  # type: ignore[arg-type]
    return keys


def _water_hydrogen_positions(
    oxygen: Coord,
    *,
    ligand_polar: Sequence[Coord],
    protein_polar: Sequence[Coord],
) -> tuple[Coord, Coord]:
    targets = _nearest_points(oxygen, ligand_polar, limit=1) + _nearest_points(
        oxygen,
        protein_polar,
        limit=1,
    )
    vectors = [_unit(_sub(point, oxygen)) for point in targets if _norm(_sub(point, oxygen)) > 0.001]
    if not vectors:
        vectors = [(1.0, 0.0, 0.0)]
    first = vectors[0]
    if len(vectors) > 1 and _dot(first, vectors[1]) < 0.75:
        second = vectors[1]
    else:
        second = _rotate_around_z(first, WATER_H_ANGLE_DEG)
    return (_add_scaled(oxygen, first, WATER_H_BOND_A), _add_scaled(oxygen, second, WATER_H_BOND_A))


def _nearest_points(origin: Coord, points: Sequence[Coord], *, limit: int) -> list[Coord]:
    return [
        point
        for _, point in sorted((_distance2(origin, point), point) for point in points)[:limit]
    ]


def _collect_polar_points(path: Path | None) -> list[Coord]:
    points: list[Coord] = []
    if path is None or not path.exists():
        return points
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if _line_element(line) not in POLAR_ELEMENTS:
                continue
            xyz = line_xyz(line)
            if xyz is not None:
                points.append(xyz)
    return points


def _has_any_sidechain_hydrogen(records: Sequence[_Record]) -> bool:
    return any(record.element == "H" and record.atom_name not in {"H", "H1", "H2", "H3"} for record in records)


def _atom_xyz(records: Sequence[_Record], atom_name: str) -> Coord | None:
    for record in records:
        if record.atom_name == atom_name and record.xyz is not None:
            return record.xyz
    return None


def _min_distance(point: Coord, targets: Sequence[Coord]) -> float:
    if not targets:
        return math.inf
    return math.sqrt(min(_distance2(point, target) for target in targets))


def _max_serial(records: Sequence[_Record]) -> int:
    serial = 0
    for record in records:
        if record.line.startswith(("ATOM  ", "HETATM")):
            try:
                serial = max(serial, int(record.line[6:11]))
            except ValueError:
                continue
    return serial


def _format_water_hydrogen_line(template: str, serial: int, atom_name: str, xyz: Coord) -> str:
    return (
        f"HETATM{serial:5d} {atom_name:>4} {template[17:20]}"
        f" {template[21:22]}{template[22:26]}{template[26:27]}   "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"{template[54:60] if len(template) >= 60 else '  1.00'}"
        f"{template[60:66] if len(template) >= 66 else ' 20.00'}"
        "           H  \n"
    )


def _replace_resname(line: str, resname: str) -> str:
    return line[:17] + f"{resname:>3}" + line[20:]


def _replace_atom_name_and_element(line: str, atom_name: str, element: str) -> str:
    padded = f"{atom_name:>4}"
    if len(line) >= 78:
        return line[:12] + padded + line[16:76] + f"{element:>2}" + line[78:]
    return line[:12] + padded + line[16:].rstrip("\n").ljust(76) + f"{element:>2}\n"


def _line_element(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    atom = line[12:16].strip()
    return (element or atom[:1]).upper()


def _applied_action(
    records: Sequence[_Record],
    recommendation: Mapping[str, object],
    action: str,
    mode: str,
) -> dict[str, object]:
    key = records[0].key if records and records[0].key else ("", "", "", "")
    return {
        "residue": _format_key(key),
        "site_type": str(recommendation.get("site_type", "")),
        "status": "applied",
        "action": action,
        "mode": mode,
    }


def _skipped_action(
    records: Sequence[_Record],
    recommendation: Mapping[str, object],
    reason: str,
    *,
    recommended: str = "",
) -> dict[str, object]:
    key = records[0].key if records and records[0].key else ("", "", "", "")
    return {
        "residue": _format_key(key),
        "site_type": str(recommendation.get("site_type", "")),
        "status": "skipped",
        "reason": reason,
        "recommended": recommended,
    }


def _water_action(key: ResidueKey, status: str, action: str) -> dict[str, object]:
    return {
        "residue": _format_key(key),
        "site_type": "selected_water_orientation",
        "status": status,
        "action": action,
    }


def _applied_count(actions: Sequence[Mapping[str, object]]) -> int:
    return sum(1 for action in actions if action.get("status") == "applied")


def _write_summary(
    audit_path: Path | None,
    status: str,
    input_pdb: Path,
    output_pdb: Path,
    actions: Sequence[Mapping[str, object]],
    *,
    low_confidence_skipped: int,
) -> dict[str, object]:
    applied = [action for action in actions if action.get("status") == "applied"]
    summary = {
        "network_application_status": status,
        "network_application_input_pdb": str(input_pdb),
        "network_application_output_pdb": str(output_pdb),
        "network_application_action_count": len(actions),
        "network_application_applied_count": len(applied),
        "network_application_histidine_renamed_count": _count_action(applied, "histidine_renamed"),
        "network_application_sidechain_flip_count": _count_action(applied, "amide_atom_names_flipped"),
        "network_application_water_hydrogen_added_count": _count_action(applied, "water_hydrogens_added"),
        "network_application_low_confidence_skipped_count": low_confidence_skipped,
        "network_application_audit_path": str(audit_path) if audit_path else "",
    }
    if audit_path is not None:
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        audit_path.write_text(
            json.dumps({"summary": summary, "actions": list(actions)}, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    return summary


def _count_action(actions: Sequence[Mapping[str, object]], action_name: str) -> int:
    return sum(1 for action in actions if action.get("action") == action_name)


def _write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.writelines(lines)


def _load_json(path: Path) -> Mapping[str, object]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, Mapping) else {}


def _format_key(key: ResidueKey) -> str:
    resname, chain, resseq, icode = key
    suffix = icode if icode else ""
    return f"{resname}:{chain or '-'}:{resseq}{suffix}"


def _sub(a: Coord, b: Coord) -> Coord:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _norm(a: Coord) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Coord) -> Coord:
    n = _norm(a)
    if n <= 0.0:
        return (1.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def _dot(a: Coord, b: Coord) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _add_scaled(origin: Coord, vector: Coord, scale: float) -> Coord:
    return (
        origin[0] + vector[0] * scale,
        origin[1] + vector[1] * scale,
        origin[2] + vector[2] * scale,
    )


def _rotate_around_z(vector: Coord, angle_deg: float) -> Coord:
    angle = math.radians(angle_deg)
    x, y, z = vector
    rotated = (
        x * math.cos(angle) - y * math.sin(angle),
        x * math.sin(angle) + y * math.cos(angle),
        z,
    )
    return _unit(rotated)


def _distance2(a: Coord, b: Coord) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
