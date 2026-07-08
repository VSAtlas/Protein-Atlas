"""Density, map, and conservation evidence for binding-site waters."""

from __future__ import annotations

import json
import math
import os
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from protein_prep.pdb_records import (
    atom_name,
    chain_id,
    line_element,
    line_xyz,
    residue_key,
    residue_name,
    resseq,
)

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
POLAR_ELEMENTS = {"N", "O", "S"}
METAL_NAMES = {
    "CA",
    "CD",
    "CO",
    "CU",
    "FE",
    "MG",
    "MN",
    "NI",
    "ZN",
}
DEFAULT_CONSERVATION_REFS: dict[str, tuple[str, ...]] = {
    "3EML": ("3PWH", "4EIY", "5IU4"),
    "1UYG": ("1UY6", "1YER", "1YC4"),
    "1XL2": ("1HXW", "1HSG", "2Q3K"),
}


@dataclass(frozen=True)
class AtomPoint:
    resname: str
    chain: str
    resseq: str
    element: str
    record: str
    xyz: tuple[float, float, float]
    bfactor: float | None = None


@dataclass(frozen=True)
class WaterPoint:
    key: tuple[str, str, str, str]
    xyz: tuple[float, float, float]
    occupancy: float | None
    bfactor: float | None


def augment_water_policy_evidence(
    *,
    raw_pdb: Path,
    prepared_pdb: Path,
    ligand_pdb: Path,
    rows: Sequence[Mapping[str, object]],
    target_dir: Path,
    pdb_id: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Add publication-style water evidence fields to water-policy rows."""

    sidecar_path = target_dir / "water_evidence_audit.json"
    raw_waters = {water.key: water for water in _collect_waters(raw_pdb)}
    protein_atoms = _collect_atoms(raw_pdb, protein_only=True)
    resolution = _parse_resolution(raw_pdb)
    refs = _conservation_refs_for_pdb(pdb_id)
    conservation = _conservation_matches(raw_pdb, raw_waters.values(), refs, target_dir)
    maps = _map_paths_for_pdb(pdb_id, target_dir)
    map_values = _sample_maps(raw_waters, maps)
    enriched = [
        _enrich_row(
            dict(row),
            water=raw_waters.get(_row_key(row)),
            protein_atoms=protein_atoms,
            resolution=resolution,
            conservation=conservation,
            map_values=map_values,
        )
        for row in rows
    ]
    summary = _summarize(enriched, refs, resolution, sidecar_path, maps)
    _write_sidecar(
        sidecar_path,
        raw_pdb=raw_pdb,
        prepared_pdb=prepared_pdb,
        ligand_pdb=ligand_pdb,
        summary=summary,
        rows=enriched,
    )
    return summary, enriched


def _enrich_row(
    row: dict[str, object],
    *,
    water: WaterPoint | None,
    protein_atoms: Sequence[AtomPoint],
    resolution: float | None,
    conservation: Mapping[tuple[str, str, str, str], Mapping[str, object]],
    map_values: Mapping[tuple[str, str, str, str], Mapping[str, object]],
) -> dict[str, object]:
    local_b = _local_protein_bfactor(water, protein_atoms)
    density = _density_status(row, water, resolution, local_b)
    conserved = conservation.get(_row_key(row), {})
    maps = map_values.get(_row_key(row), {})
    row.update(
        {
            "resolution_a": _rounded_or_blank(resolution),
            "local_protein_bfactor": _rounded_or_blank(local_b),
            "density_bfactor_delta": _rounded_or_blank(_delta(water.bfactor if water else None, local_b)),
            "density_quality_status": density,
            "density_quality_supported": density == "supported",
            "conserved_match_count": _int_value(conserved.get("match_count", 0)),
            "conserved_role_match_count": _int_value(conserved.get("role_match_count", 0)),
            "conserved_reference_pdbs": ";".join(_string_list(conserved.get("reference_pdbs", []))),
            "conservation_supported": bool(conserved.get("supported", False)),
            "map_support_status": str(maps.get("status", "not_configured")),
            "map_2fofc_value": _rounded_or_blank(_as_float_or_none(maps.get("two_fofc"))),
            "map_fofc_value": _rounded_or_blank(_as_float_or_none(maps.get("fofc"))),
            "map_supported": bool(maps.get("supported", False)),
            "water_publication_evidence_supported": _row_evidence_supported(row, density, conserved, maps),
        }
    )
    return row


def _density_status(
    row: Mapping[str, object],
    water: WaterPoint | None,
    resolution: float | None,
    local_b: float | None,
) -> str:
    review = _density_review_reason(row, water, resolution, local_b)
    return review or "supported"


def _density_review_reason(
    row: Mapping[str, object],
    water: WaterPoint | None,
    resolution: float | None,
    local_b: float | None,
) -> str:
    if water is None:
        return "missing_water_coordinates"
    metadata_review = _density_metadata_review(water, resolution, local_b)
    if metadata_review:
        return metadata_review
    if _flagged(row, "low_occupancy") or _flagged(row, "high_bfactor"):
        return "quality_flag_review"
    if not (bool(row.get("near_protein_polar")) or bool(row.get("near_metal"))):
        return "geometry_review"
    return ""


def _density_metadata_review(
    water: WaterPoint,
    resolution: float | None,
    local_b: float | None,
) -> str:
    if resolution is None:
        return "missing_resolution"
    if resolution > 2.5:
        return "resolution_review"
    if water.occupancy is None or water.occupancy < 0.8:
        return "occupancy_review"
    if water.bfactor is None:
        return "missing_water_bfactor"
    if local_b is None:
        return "missing_local_bfactor"
    if water.bfactor > max(local_b + 15.0, local_b * 1.5, 45.0):
        return "bfactor_review"
    return ""


def _row_evidence_supported(
    row: Mapping[str, object],
    density: str,
    conserved: Mapping[str, object],
    maps: Mapping[str, object],
) -> bool:
    if not _policy_selected(row):
        return False
    return density == "supported" or bool(conserved.get("supported")) or bool(maps.get("supported"))


def _summarize(
    rows: Sequence[Mapping[str, object]],
    refs: Sequence[str],
    resolution: float | None,
    sidecar_path: Path,
    maps: Mapping[str, Path],
) -> dict[str, object]:
    selected = [row for row in rows if _policy_selected(row)]
    counts = _selected_evidence_counts(selected)
    status = _evidence_status(len(selected), counts["any"])
    return {
        "water_evidence_status": status,
        "water_evidence_supported": status in {"no_selected_waters", "supported"},
        "water_density_resolution_a": _rounded_or_blank(resolution),
        "water_density_supported_selected_count": counts["density"],
        "water_conservation_reference_count": len(refs),
        "water_conservation_supported_selected_count": counts["conservation"],
        "water_map_input_count": len(maps),
        "water_map_supported_selected_count": counts["map"],
        "water_selected_evidence_supported_count": counts["any"],
        "water_selected_evidence_required_count": len(selected),
        "water_evidence_audit_path": str(sidecar_path),
    }


def _selected_evidence_counts(rows: Sequence[Mapping[str, object]]) -> dict[str, int]:
    return {
        "any": _count_bool(rows, "water_publication_evidence_supported"),
        "density": _count_bool(rows, "density_quality_supported"),
        "conservation": _count_bool(rows, "conservation_supported"),
        "map": _count_bool(rows, "map_supported"),
    }


def _count_bool(rows: Sequence[Mapping[str, object]], key: str) -> int:
    return sum(1 for row in rows if bool(row.get(key)))


def _evidence_status(
    selected_count: int,
    supported_count: int,
) -> str:
    if selected_count <= 0:
        return "no_selected_waters"
    if supported_count == selected_count:
        return "supported"
    if supported_count:
        return "partial_review"
    return "review"


def _conservation_refs_for_pdb(pdb_id: str) -> tuple[str, ...]:
    raw = os.environ.get("ATLAS_WATER_CONSERVATION_REFS", "").strip()
    parsed = _parse_conservation_refs(raw, pdb_id)
    if parsed:
        return parsed
    return DEFAULT_CONSERVATION_REFS.get(pdb_id.strip().upper(), ())


def _parse_conservation_refs(raw: str, pdb_id: str) -> tuple[str, ...]:
    if not raw:
        return ()
    key = pdb_id.strip().upper()
    parsed = _parse_conservation_refs_json(raw, key)
    if parsed:
        return parsed
    return _parse_conservation_refs_blocks(raw, key)


def _parse_conservation_refs_json(raw: str, pdb_id: str) -> tuple[str, ...]:
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        return ()
    if isinstance(data, Mapping):
        value = data.get(pdb_id) or data.get(pdb_id.lower()) or []
        if isinstance(value, str):
            return _split_refs(value)
        if isinstance(value, list | tuple):
            return tuple(str(item).strip().upper() for item in value if str(item).strip())
    return ()


def _parse_conservation_refs_blocks(raw: str, pdb_id: str) -> tuple[str, ...]:
    for block in raw.split(";"):
        if ":" not in block:
            continue
        block_key, block_refs = block.split(":", 1)
        if block_key.strip().upper() == pdb_id:
            return _split_refs(block_refs)
    return ()


def _split_refs(raw: str) -> tuple[str, ...]:
    return tuple(part.strip().upper() for part in re.split(r"[, ]+", raw) if part.strip())


def _conservation_matches(
    target_pdb: Path,
    waters: Iterable[WaterPoint],
    refs: Sequence[str],
    target_dir: Path,
) -> dict[tuple[str, str, str, str], dict[str, object]]:
    if not refs:
        return {}
    target_ca = _ca_by_chain(target_pdb)
    rows = {
        water.key: {
            "_water": water,
            "match_count": 0,
            "role_match_count": 0,
            "reference_pdbs": [],
        }
        for water in waters
    }
    for ref in refs:
        ref_path = _reference_pdb_path(ref, target_dir)
        if not _ensure_reference_pdb(ref, ref_path):
            continue
        transform = _best_ca_transform(target_ca, _ca_by_chain(ref_path))
        if transform is None:
            continue
        ref_waters = _transformed_ref_waters(ref_path, transform)
        ref_roles = _reference_water_roles(ref_path)
        _update_conservation_rows(rows, ref, ref_waters, ref_roles)
    for row in rows.values():
        row["supported"] = (
            _int_value(row.get("match_count", 0)) >= 1
            and _int_value(row.get("role_match_count", 0)) >= 1
        )
        row.pop("_water", None)
    return rows


def _update_conservation_rows(
    rows: dict[tuple[str, str, str, str], dict[str, object]],
    ref: str,
    ref_waters: Sequence[tuple[WaterPoint, tuple[float, float, float]]],
    ref_roles: Mapping[tuple[str, str, str, str], set[str]],
) -> None:
    for key, row in rows.items():
        target_water = _target_water_from_key(key, rows)
        if target_water is None:
            continue
        match = _nearest_ref_water(target_water.xyz, ref_waters)
        if match is None:
            continue
        water, _ = match
        row["match_count"] = _int_value(row.get("match_count", 0)) + 1
        reference_pdbs = row.get("reference_pdbs")
        if not isinstance(reference_pdbs, list):
            reference_pdbs = []
            row["reference_pdbs"] = reference_pdbs
        reference_pdbs.append(ref.upper())
        if ref_roles.get(water.key):
            row["role_match_count"] = _int_value(row.get("role_match_count", 0)) + 1


def _target_water_from_key(
    key: tuple[str, str, str, str],
    rows: Mapping[tuple[str, str, str, str], Mapping[str, object]],
) -> WaterPoint | None:
    raw = rows.get(key, {}).get("_water")
    return raw if isinstance(raw, WaterPoint) else None


def _nearest_ref_water(
    target_xyz: tuple[float, float, float],
    waters: Sequence[tuple[WaterPoint, tuple[float, float, float]]],
) -> tuple[WaterPoint, float] | None:
    best: tuple[WaterPoint, float] | None = None
    for water, xyz in waters:
        dist = _distance(target_xyz, xyz)
        if dist <= 1.5 and (best is None or dist < best[1]):
            best = (water, dist)
    return best


def _best_ca_transform(
    target: Mapping[str, Mapping[str, tuple[float, float, float]]],
    ref: Mapping[str, Mapping[str, tuple[float, float, float]]],
) -> tuple[object, object, object] | None:
    best_pair: tuple[list[tuple[float, float, float]], list[tuple[float, float, float]]] | None = None
    for target_atoms in target.values():
        for ref_atoms in ref.values():
            common = sorted(set(target_atoms) & set(ref_atoms))
            if len(common) < 12:
                continue
            pair = ([ref_atoms[pos] for pos in common], [target_atoms[pos] for pos in common])
            if best_pair is None or len(pair[0]) > len(best_pair[0]):
                best_pair = pair
    if best_pair is None:
        return None
    return _kabsch_transform(*best_pair)


def _kabsch_transform(
    ref_coords: Sequence[tuple[float, float, float]],
    target_coords: Sequence[tuple[float, float, float]],
) -> tuple[object, object, object] | None:
    try:
        import numpy as np
    except Exception:
        return None
    ref_arr = np.asarray(ref_coords, dtype=float)
    target_arr = np.asarray(target_coords, dtype=float)
    ref_centroid = ref_arr.mean(axis=0)
    target_centroid = target_arr.mean(axis=0)
    ref_centered = ref_arr - ref_centroid
    target_centered = target_arr - target_centroid
    u, _, vt = np.linalg.svd(ref_centered.T @ target_centered)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    return ref_centroid, target_centroid, rot


def _apply_transform(
    xyz: tuple[float, float, float],
    transform: tuple[object, object, object],
) -> tuple[float, float, float]:
    import numpy as np

    ref_centroid, target_centroid, rot = transform
    value = (np.asarray(xyz, dtype=float) - ref_centroid) @ rot + target_centroid
    return float(value[0]), float(value[1]), float(value[2])


def _transformed_ref_waters(
    ref_path: Path,
    transform: tuple[object, object, object],
) -> list[tuple[WaterPoint, tuple[float, float, float]]]:
    return [(water, _apply_transform(water.xyz, transform)) for water in _collect_waters(ref_path)]


def _reference_water_roles(ref_path: Path) -> dict[tuple[str, str, str, str], set[str]]:
    waters = _collect_waters(ref_path)
    atoms = _collect_atoms(ref_path, protein_only=False)
    roles: dict[tuple[str, str, str, str], set[str]] = {}
    for water in waters:
        roles[water.key] = _roles_for_water(water.xyz, atoms)
    return roles


def _roles_for_water(
    xyz: tuple[float, float, float],
    atoms: Sequence[AtomPoint],
) -> set[str]:
    roles: set[str] = set()
    for atom in atoms:
        role = _water_contact_role(atom, _distance(xyz, atom.xyz))
        if role:
            roles.add(role)
    return roles


def _water_contact_role(atom: AtomPoint, dist: float) -> str:
    if atom.record == "ATOM" and atom.element in POLAR_ELEMENTS and dist <= 3.6:
        return "protein_polar"
    if atom.record == "HETATM" and atom.resname in METAL_NAMES and dist <= 3.2:
        return "metal"
    if atom.record == "HETATM" and atom.resname not in WATER_NAMES and dist <= 3.8:
        return "het"
    return ""


def _reference_pdb_path(ref: str, target_dir: Path) -> Path:
    return target_dir / "water_conservation_refs" / f"{ref.upper()}.pdb"


def _ensure_reference_pdb(ref: str, path: Path) -> bool:
    if path.exists():
        return True
    path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{ref.upper()}.pdb"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            text = response.read().decode("utf-8", errors="ignore")
        path.write_text(text, encoding="utf-8")
        return path.exists() and path.stat().st_size > 1000
    except Exception:
        return False


def _ca_by_chain(path: Path) -> dict[str, dict[str, tuple[float, float, float]]]:
    chains: dict[str, dict[str, tuple[float, float, float]]] = {}
    for line in _atom_lines(path):
        if line.startswith("ATOM") and atom_name(line) == "CA":
            xyz = line_xyz(line)
            if xyz is not None:
                chains.setdefault(chain_id(line), {})[resseq(line)] = xyz
    return chains


def _collect_waters(path: Path) -> list[WaterPoint]:
    waters: list[WaterPoint] = []
    for line in _atom_lines(path):
        if not line.startswith("HETATM") or residue_name(line) not in WATER_NAMES:
            continue
        if line_element(line) == "H":
            continue
        xyz = line_xyz(line)
        if xyz is None:
            continue
        waters.append(
            WaterPoint(
                key=residue_key(line),
                xyz=xyz,
                occupancy=_safe_float(line[54:60]),
                bfactor=_safe_float(line[60:66]),
            )
        )
    return waters


def _collect_atoms(path: Path, *, protein_only: bool) -> list[AtomPoint]:
    atoms: list[AtomPoint] = []
    for line in _atom_lines(path):
        if not line.startswith(("ATOM", "HETATM")) or line_element(line) == "H":
            continue
        if protein_only and not line.startswith("ATOM"):
            continue
        xyz = line_xyz(line)
        if xyz is None:
            continue
        atoms.append(
            AtomPoint(
                resname=residue_name(line),
                chain=chain_id(line),
                resseq=resseq(line),
                element=line_element(line),
                record=line[:6].strip(),
                xyz=xyz,
                bfactor=_safe_float(line[60:66]),
            )
        )
    return atoms


def _atom_lines(path: Path) -> Iterable[str]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return [line for line in handle if line.startswith(("ATOM", "HETATM", "REMARK"))]


def _parse_resolution(path: Path) -> float | None:
    pattern = re.compile(r"RESOLUTION\.\s+([0-9]+(?:\.[0-9]+)?)\s+ANGSTROMS", re.I)
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM", "MODEL")):
                break
            if not line.startswith("REMARK   2"):
                continue
            match = pattern.search(line)
            if match:
                return _safe_float(match.group(1))
    return None


def _local_protein_bfactor(
    water: WaterPoint | None,
    protein_atoms: Sequence[AtomPoint],
    radius: float = 4.0,
) -> float | None:
    if water is None:
        return None
    values = [
        atom.bfactor
        for atom in protein_atoms
        if atom.bfactor is not None and _distance(water.xyz, atom.xyz) <= radius
    ]
    if not values:
        return None
    return sum(values) / len(values)


def _map_paths_for_pdb(pdb_id: str, target_dir: Path) -> dict[str, Path]:
    key = pdb_id.strip().upper()
    paths: dict[str, Path] = {}
    for label, env_name in (("two_fofc", "ATLAS_WATER_2FOFC_MAP"), ("fofc", "ATLAS_WATER_FOFC_MAP")):
        raw = os.environ.get(env_name, "").strip()
        candidate = _path_from_mapping(raw, key) if raw else None
        if candidate is None:
            suffix = "2fofc" if label == "two_fofc" else "fofc"
            candidate = target_dir / f"{key}_{suffix}.ccp4"
        if candidate.exists():
            paths[label] = candidate
    return paths


def _path_from_mapping(raw: str, pdb_id: str) -> Path | None:
    if not raw:
        return None
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        data = None
    if isinstance(data, Mapping):
        value = data.get(pdb_id) or data.get(pdb_id.lower())
        return Path(str(value)) if value else None
    return Path(raw)


def _sample_maps(
    waters: Mapping[tuple[str, str, str, str], WaterPoint],
    paths: Mapping[str, Path],
) -> dict[tuple[str, str, str, str], dict[str, object]]:
    if not paths:
        return {}
    grids = _read_ccp4_grids(paths)
    out: dict[tuple[str, str, str, str], dict[str, object]] = {}
    for key, water in waters.items():
        two_fofc = _sample_grid(grids.get("two_fofc"), water.xyz)
        fofc = _sample_grid(grids.get("fofc"), water.xyz)
        out[key] = {
            "status": "checked" if grids else "map_read_failed",
            "two_fofc": two_fofc,
            "fofc": fofc,
            "supported": two_fofc is not None and two_fofc > 1.0 and (fofc is None or fofc > -3.0),
        }
    return out


def _read_ccp4_grids(paths: Mapping[str, Path]) -> dict[str, object]:
    try:
        import gemmi
    except Exception:
        return {}
    grids: dict[str, object] = {}
    for label, path in paths.items():
        try:
            ccp4 = gemmi.read_ccp4_map(str(path))
            ccp4.setup(0.0)
            grids[label] = ccp4.grid
        except Exception:
            continue
    return grids


def _sample_grid(grid: object | None, xyz: tuple[float, float, float]) -> float | None:
    if grid is None:
        return None
    try:
        import gemmi

        interpolate = getattr(grid, "interpolate_value", None)
        if callable(interpolate):
            return float(interpolate(gemmi.Position(*xyz)))
        return None
    except Exception:
        return None


def _policy_selected(row: Mapping[str, object]) -> bool:
    if "water_selected_for_receptor" in row:
        return bool(row.get("water_selected_for_receptor"))
    if not bool(row.get("retained")):
        return False
    if not bool(row.get("bridging_candidate")):
        return False
    return not _flagged(row, "low_occupancy") and not _flagged(row, "high_bfactor")


def _flagged(row: Mapping[str, object], flag: str) -> bool:
    flags = row.get("quality_flags", [])
    return isinstance(flags, list) and flag in flags


def _row_key(row: Mapping[str, object]) -> tuple[str, str, str, str]:
    raw = row.get("residue_key", [])
    if not isinstance(raw, list | tuple):
        return ("", "", "", "")
    parts = [str(part) for part in raw]
    while len(parts) < 4:
        parts.append("")
    return tuple(parts[:4])  # type: ignore[return-value]


def _safe_float(text: object) -> float | None:
    try:
        value = float(str(text).strip())
    except Exception:
        return None
    return value if math.isfinite(value) else None


def _int_value(value: object) -> int:
    try:
        return int(str(value))
    except Exception:
        return 0


def _string_list(value: object) -> list[str]:
    if not isinstance(value, list | tuple):
        return []
    return [str(item) for item in value if str(item)]


def _as_float_or_none(value: object) -> float | None:
    return _safe_float(value)


def _rounded_or_blank(value: float | None) -> float | str:
    if value is None or not math.isfinite(value):
        return ""
    return round(float(value), 3)


def _delta(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return a - b


def _distance(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _write_sidecar(
    sidecar_path: Path,
    *,
    raw_pdb: Path,
    prepared_pdb: Path,
    ligand_pdb: Path,
    summary: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
) -> None:
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.write_text(
        json.dumps(
            {
                "raw_pdb": str(raw_pdb),
                "prepared_pdb": str(prepared_pdb),
                "ligand_pdb": str(ligand_pdb),
                "summary": dict(summary),
                "waters": list(rows),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
