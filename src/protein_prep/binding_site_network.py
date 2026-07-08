"""Open, bounded binding-site H-bond network optimizer.

This is a Protoss-inspired heuristic: it optimizes only the local component
around retained HET groups instead of attempting whole-protein enumeration.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from protein_prep.pdb_records import line_xyz, residue_key

Coord = tuple[float, float, float]

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
METAL_ELEMENTS = {"CA", "CD", "CO", "CU", "FE", "MG", "MN", "NI", "ZN"}
POLAR_ELEMENTS = {"N", "O", "S"}
FLIPPABLE_RESIDUES = {"ASN", "GLN", "HIS", "HID", "HIE", "HIP"}
HBOND_CUTOFF_A = 3.5
DEFAULT_RADIUS_A = 8.0


@dataclass(frozen=True)
class NetworkAtom:
    key: tuple[str, str, str, str]
    atom_name: str
    resname: str
    record: str
    element: str
    xyz: Coord


@dataclass(frozen=True)
class NetworkSite:
    key: tuple[str, str, str, str]
    site_type: str
    resname: str
    distance_to_het_a: float
    mode_count: int


def optimize_binding_site_network(
    *,
    pdb_path: Path,
    het_resname: str,
    het_instances: Sequence[Any],
    candidates: Sequence[Mapping[str, object]],
    target_ph: float,
    radius_a: float = DEFAULT_RADIUS_A,
    max_exact_modes: int = 4096,
) -> dict[str, object]:
    """Score HET states against a local variable H-bond/contact network."""

    het_points = _het_points(het_instances)
    if not het_points or not candidates:
        return _empty_network("no_het_points_or_candidates", radius_a, target_ph)

    atoms = _collect_network_atoms(pdb_path)
    shell = _binding_site_shell(atoms, het_points, radius_a, het_resname)
    sites = _variable_sites(shell, het_points)
    contacts = _site_contacts(shell, sites)
    components = _connected_components(sites, contacts)
    site_modes = _estimate_mode_count(sites, candidate_count=len(candidates))
    candidate_scores = {
        str(row.get("smiles", "")): _score_candidate_against_sites(row, sites, contacts)
        for row in candidates
        if str(row.get("smiles", ""))
    }
    return {
        "status": "optimized" if sites else "no_variable_sites",
        "method": "open_protoss_inspired_bounded_binding_site_network",
        "not_protoss_equivalent": True,
        "target_ph": round(float(target_ph), 2),
        "radius_a": radius_a,
        "hbond_cutoff_a": HBOND_CUTOFF_A,
        "het_resname": het_resname,
        "variable_site_count": len(sites),
        "site_type_counts": _format_counter(Counter(site.site_type for site in sites)),
        "component_count": len(components),
        "largest_component_site_count": max((len(item) for item in components), default=0),
        "mode_count_estimate": site_modes,
        "mode_search": "factorized_exact" if site_modes <= max_exact_modes else "pruned_factorized",
        "candidate_scores": candidate_scores,
        "recommended_site_modes": _recommended_site_modes(sites, contacts),
    }


def _empty_network(reason: str, radius_a: float, target_ph: float) -> dict[str, object]:
    return {
        "status": reason,
        "method": "open_protoss_inspired_bounded_binding_site_network",
        "not_protoss_equivalent": True,
        "target_ph": round(float(target_ph), 2),
        "radius_a": radius_a,
        "hbond_cutoff_a": HBOND_CUTOFF_A,
        "variable_site_count": 0,
        "site_type_counts": "",
        "component_count": 0,
        "largest_component_site_count": 0,
        "mode_count_estimate": 0,
        "mode_search": "not_applicable",
        "candidate_scores": {},
        "recommended_site_modes": [],
    }


def _het_points(het_instances: Sequence[Any]) -> list[Coord]:
    points: list[Coord] = []
    for instance in het_instances:
        for _atom_name, xyz in getattr(instance, "atoms", ()):
            points.append(xyz)
    return points


def _collect_network_atoms(path: Path) -> list[NetworkAtom]:
    atoms: list[NetworkAtom] = []
    if not path.exists():
        return atoms
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            element = _line_element(line)
            if element == "H":
                continue
            xyz = line_xyz(line)
            if xyz is None:
                continue
            key = residue_key(line)
            atoms.append(
                NetworkAtom(
                    key=key,
                    atom_name=line[12:16].strip().upper(),
                    resname=key[0],
                    record=line[:6].strip(),
                    element=element,
                    xyz=xyz,
                )
            )
    return atoms


def _binding_site_shell(
    atoms: Sequence[NetworkAtom],
    het_points: Sequence[Coord],
    radius_a: float,
    het_resname: str,
) -> list[NetworkAtom]:
    radius2 = radius_a * radius_a
    shell: list[NetworkAtom] = []
    for atom in atoms:
        if atom.record == "HETATM" and atom.resname == het_resname:
            continue
        if any(_distance2(atom.xyz, point) <= radius2 for point in het_points):
            shell.append(atom)
    return shell


def _variable_sites(
    shell: Sequence[NetworkAtom],
    het_points: Sequence[Coord],
) -> list[NetworkSite]:
    by_residue: dict[tuple[str, str, str, str], list[NetworkAtom]] = defaultdict(list)
    for atom in shell:
        by_residue[atom.key].append(atom)
    sites: list[NetworkSite] = []
    for key, atoms in sorted(by_residue.items()):
        site_type, mode_count = _site_type_and_modes(atoms, het_points)
        if not site_type:
            continue
        sites.append(
            NetworkSite(
                key=key,
                site_type=site_type,
                resname=atoms[0].resname,
                distance_to_het_a=_min_distance_to_points(atoms, het_points),
                mode_count=mode_count,
            )
        )
    return sites


def _site_type_and_modes(
    atoms: Sequence[NetworkAtom],
    het_points: Sequence[Coord],
) -> tuple[str, int]:
    resname = atoms[0].resname
    record = atoms[0].record
    element = atoms[0].element
    if record == "HETATM":
        return _het_site_type_and_modes(resname, element, atoms, het_points)
    if record == "ATOM":
        return _protein_site_type_and_modes(resname, atoms, het_points)
    return "", 0


def _het_site_type_and_modes(
    resname: str,
    element: str,
    atoms: Sequence[NetworkAtom],
    het_points: Sequence[Coord],
) -> tuple[str, int]:
    if element in METAL_ELEMENTS:
        return "metal_coordination", 1
    if resname in WATER_NAMES and _water_bridges(atoms, het_points):
        return "bridging_water_orientation", 3
    return "", 0


def _protein_site_type_and_modes(
    resname: str,
    atoms: Sequence[NetworkAtom],
    het_points: Sequence[Coord],
) -> tuple[str, int]:
    if resname in {"ASN", "GLN"}:
        return "sidechain_flip", 2
    if resname in {"HIS", "HID", "HIE", "HIP"}:
        return "histidine_tautomer_protomer", 3
    if resname in {"ASP", "GLU"} and _has_polar_contact_to_het(atoms, het_points):
        return "borderline_carboxylate_protonation", 2
    return "", 0


def _water_bridges(atoms: Sequence[NetworkAtom], het_points: Sequence[Coord]) -> bool:
    return _has_polar_contact_to_het(atoms, het_points)


def _has_polar_contact_to_het(
    atoms: Sequence[NetworkAtom],
    het_points: Sequence[Coord],
) -> bool:
    cutoff2 = HBOND_CUTOFF_A * HBOND_CUTOFF_A
    return any(
        atom.element in POLAR_ELEMENTS
        and any(_distance2(atom.xyz, point) <= cutoff2 for point in het_points)
        for atom in atoms
    )


def _site_contacts(
    shell: Sequence[NetworkAtom],
    sites: Sequence[NetworkSite],
) -> list[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]]:
    site_keys = {site.key for site in sites}
    atoms_by_site: dict[tuple[str, str, str, str], list[NetworkAtom]] = defaultdict(list)
    for atom in shell:
        if atom.key in site_keys and atom.element in POLAR_ELEMENTS | METAL_ELEMENTS:
            atoms_by_site[atom.key].append(atom)
    contacts: set[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]] = set()
    keys = sorted(atoms_by_site)
    cutoff2 = HBOND_CUTOFF_A * HBOND_CUTOFF_A
    for left_index, left_key in enumerate(keys):
        for right_key in keys[left_index + 1 :]:
            if _sites_contact(atoms_by_site[left_key], atoms_by_site[right_key], cutoff2):
                contacts.add((left_key, right_key))
    return sorted(contacts)


def _sites_contact(
    left: Sequence[NetworkAtom],
    right: Sequence[NetworkAtom],
    cutoff2: float,
) -> bool:
    return any(_distance2(a.xyz, b.xyz) <= cutoff2 for a in left for b in right)


def _connected_components(
    sites: Sequence[NetworkSite],
    contacts: Sequence[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]],
) -> list[set[tuple[str, str, str, str]]]:
    graph: dict[tuple[str, str, str, str], set[tuple[str, str, str, str]]] = {
        site.key: set() for site in sites
    }
    for left, right in contacts:
        graph.setdefault(left, set()).add(right)
        graph.setdefault(right, set()).add(left)
    remaining = set(graph)
    components: list[set[tuple[str, str, str, str]]] = []
    while remaining:
        start = remaining.pop()
        component = _walk_component(start, graph)
        remaining.difference_update(component)
        components.append(component)
    return components


def _walk_component(
    start: tuple[str, str, str, str],
    graph: Mapping[tuple[str, str, str, str], set[tuple[str, str, str, str]]],
) -> set[tuple[str, str, str, str]]:
    seen = {start}
    queue: deque[tuple[str, str, str, str]] = deque([start])
    while queue:
        current = queue.popleft()
        for neighbor in graph.get(current, set()):
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    return seen


def _score_candidate_against_sites(
    candidate: Mapping[str, object],
    sites: Sequence[NetworkSite],
    contacts: Sequence[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]],
) -> dict[str, object]:
    charge = _safe_int(candidate.get("formal_charge", 0))
    hbd = _safe_int(candidate.get("hbd", 0))
    hba = _safe_int(candidate.get("hba", 0))
    counts = Counter(site.site_type for site in sites)
    score = 0.0
    score += min(hbd, counts["borderline_carboxylate_protonation"] + counts["sidechain_flip"]) * 0.45
    score += min(hba, counts["histidine_tautomer_protomer"] + counts["metal_coordination"]) * 0.50
    score += min(hbd + hba, counts["bridging_water_orientation"]) * 0.35
    score += min(abs(charge), counts["borderline_carboxylate_protonation"]) * 0.20
    score -= max(0, abs(charge) - counts["metal_coordination"] - counts["borderline_carboxylate_protonation"]) * 0.30
    score += min(len(contacts), hbd + hba + abs(charge)) * 0.03
    return {
        "network_score": round(score, 3),
        "best_local_mode_score": round(score, 3),
        "mode_explanation": _mode_explanation(charge, hbd, hba, counts),
    }


def _mode_explanation(
    charge: int,
    hbd: int,
    hba: int,
    counts: Counter[str],
) -> str:
    terms: list[str] = []
    if counts["metal_coordination"]:
        terms.append("metal-adjacent")
    if counts["bridging_water_orientation"]:
        terms.append("water-bridged")
    if counts["sidechain_flip"] or counts["histidine_tautomer_protomer"]:
        terms.append("flip-aware")
    if counts["borderline_carboxylate_protonation"]:
        terms.append("titratable-neighbor")
    terms.append(f"charge={charge}")
    terms.append(f"hbd={hbd}")
    terms.append(f"hba={hba}")
    return ",".join(terms)


def _recommended_site_modes(
    sites: Sequence[NetworkSite],
    contacts: Sequence[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]],
) -> list[dict[str, object]]:
    contact_counts = Counter(key for pair in contacts for key in pair)
    return [
        {
            "residue": _format_key(site.key),
            "site_type": site.site_type,
            "recommended_mode": _default_mode_for_site(site),
            "mode_count": site.mode_count,
            "distance_to_het_a": round(site.distance_to_het_a, 3),
            "network_contact_count": contact_counts[site.key],
        }
        for site in sites
    ]


def _default_mode_for_site(site: NetworkSite) -> str:
    if site.site_type == "histidine_tautomer_protomer":
        return "score_HID_HIE_HIP_against_local_acceptors_donors"
    if site.site_type == "sidechain_flip":
        return "score_current_vs_180deg_flip"
    if site.site_type == "borderline_carboxylate_protonation":
        return "prefer_deprotonated_unless_unsatisfied_local_donor_network"
    if site.site_type == "bridging_water_orientation":
        return "orient_to_maximize_two_center_bridge"
    return "retain_geometry"


def _estimate_mode_count(sites: Sequence[NetworkSite], *, candidate_count: int) -> int:
    total = max(1, candidate_count)
    for site in sites:
        total *= max(1, site.mode_count)
        if total > 1_000_000_000:
            return total
    return total


def _min_distance_to_points(
    atoms: Sequence[NetworkAtom],
    points: Sequence[Coord],
) -> float:
    return math.sqrt(min(_distance2(atom.xyz, point) for atom in atoms for point in points))


def _distance2(a: Coord, b: Coord) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _line_element(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    atom_name = line[12:16].strip()
    return (element or atom_name[:2] or atom_name[:1]).upper()


def _safe_int(value: object) -> int:
    try:
        return int(str(value))
    except Exception:
        return 0


def _format_counter(counter: Counter[str]) -> str:
    return ";".join(f"{key}:{value}" for key, value in sorted(counter.items()) if value)


def _format_key(key: tuple[str, str, str, str]) -> str:
    resname, chain, resseq, icode = key
    suffix = icode if icode else ""
    return f"{resname}:{chain or '-'}:{resseq}{suffix}"
