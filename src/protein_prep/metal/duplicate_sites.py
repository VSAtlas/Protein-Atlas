"""Conservative duplicate/overlapping metal-site resolution.

This module only handles physically impossible metal-metal overlaps. Ambiguous
mixed-metal sites are left as review failures instead of being guessed.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from protein_prep.pdb_records import line_xyz

DEFAULT_OVERLAP_CUTOFF_A = 0.30
OCCUPANCY_DECISION_MARGIN = 0.05
B_FACTOR_DECISION_MARGIN = 15.0
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


@dataclass(frozen=True)
class MetalRecord:
    """One metal atom record from a PDB file."""

    line_no: int
    serial: str
    element: str
    resname: str
    chain: str
    resseq: str
    icode: str
    atom_name: str
    altloc: str
    occupancy: float | None
    b_factor: float | None
    xyz: tuple[float, float, float]


def resolve_duplicate_metal_sites(
    input_pdb: Path,
    output_pdb: Path,
    *,
    sidecar_path: Path,
    overlap_cutoff_a: float | None = None,
) -> dict[str, object]:
    """Resolve obvious duplicate metals or request review for ambiguous sites.

    The resolver is intentionally narrow. It drops duplicate records only when
    there is direct record-level support such as identical atom identity or a
    clear occupancy winner. Mixed-metal overlaps without a clear occupancy
    decision are flagged as poor Atlas candidates for review.
    """

    cutoff = (
        float(overlap_cutoff_a)
        if overlap_cutoff_a is not None
        else _env_float("ATLAS_DUPLICATE_METAL_CUTOFF_A", DEFAULT_OVERLAP_CUTOFF_A)
    )
    lines = input_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(True)
    metals = _collect_metals(lines)
    clusters = _overlap_clusters(metals, cutoff)
    decisions = [_cluster_decision(cluster, cutoff) for cluster in clusters]
    review_decisions = [
        decision for decision in decisions if bool(decision.get("review_required"))
    ]
    drop_line_numbers = _drop_line_numbers(decisions)
    summary = {
        "duplicate_metal_site_status": _status(decisions, review_decisions),
        "duplicate_metal_site_review_required": bool(review_decisions),
        "duplicate_metal_site_fail_for_review": bool(review_decisions),
        "duplicate_metal_site_overlap_cutoff_a": round(cutoff, 3),
        "duplicate_metal_site_cluster_count": len(decisions),
        "duplicate_metal_site_resolved_cluster_count": len(decisions)
        - len(review_decisions),
        "duplicate_metal_site_review_cluster_count": len(review_decisions),
        "duplicate_metal_site_dropped_atom_count": len(drop_line_numbers),
        "duplicate_metal_site_output_pdb": str(output_pdb) if drop_line_numbers else "",
        "duplicate_metal_site_audit_path": str(sidecar_path),
        "duplicate_metal_site_decisions": decisions,
    }
    _write_sidecar(sidecar_path, input_pdb, output_pdb, summary)
    if review_decisions:
        return summary
    if drop_line_numbers:
        _write_without_lines(lines, output_pdb, drop_line_numbers)
    return summary


def _status(
    decisions: list[dict[str, object]],
    review_decisions: list[dict[str, object]],
) -> str:
    if not decisions:
        return "no_duplicate_metal_sites"
    if review_decisions:
        return "review_required"
    return "resolved"


def _drop_line_numbers(decisions: list[dict[str, object]]) -> set[int]:
    selected: set[int] = set()
    for decision in decisions:
        raw = decision.get("drop_line_numbers", [])
        if not isinstance(raw, list):
            continue
        for line_no in raw:
            if isinstance(line_no, int):
                selected.add(line_no)
    return selected


def _collect_metals(lines: Iterable[str]) -> list[MetalRecord]:
    metals: list[MetalRecord] = []
    for line_no, line in enumerate(lines, start=1):
        if not line.startswith("HETATM"):
            continue
        element = _line_element(line)
        resname = line[17:20].strip().upper()
        if element not in COORDINATION_METALS and resname not in COORDINATION_METALS:
            continue
        xyz = line_xyz(line)
        if xyz is None:
            continue
        metals.append(
            MetalRecord(
                line_no=line_no,
                serial=line[6:11].strip(),
                element=element or resname,
                resname=resname,
                chain=line[21:22].strip() or "-",
                resseq=line[22:26].strip(),
                icode=line[26:27].strip(),
                atom_name=line[12:16].strip(),
                altloc=line[16:17].strip(),
                occupancy=_line_float(line, 54, 60),
                b_factor=_line_float(line, 60, 66),
                xyz=xyz,
            )
        )
    return metals


def _overlap_clusters(metals: list[MetalRecord], cutoff: float) -> list[list[MetalRecord]]:
    if len(metals) < 2:
        return []
    parent = list(range(len(metals)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    cutoff2 = cutoff * cutoff
    for left_index, left in enumerate(metals):
        for right_index in range(left_index + 1, len(metals)):
            right = metals[right_index]
            if _distance2(left.xyz, right.xyz) <= cutoff2:
                union(left_index, right_index)
    grouped: dict[int, list[MetalRecord]] = {}
    for index, metal in enumerate(metals):
        grouped.setdefault(find(index), []).append(metal)
    return [cluster for cluster in grouped.values() if len(cluster) > 1]


def _cluster_decision(
    cluster: list[MetalRecord],
    cutoff: float,
) -> dict[str, object]:
    sorted_cluster = sorted(cluster, key=lambda metal: metal.line_no)
    keep, reason, confidence = _choose_keep_record(sorted_cluster)
    review_required = keep is None
    drop = [] if keep is None else [metal for metal in sorted_cluster if metal != keep]
    return {
        "status": "review_required" if review_required else "resolved",
        "review_required": review_required,
        "reason": reason,
        "confidence": confidence,
        "max_pair_distance_a": round(_max_pair_distance(sorted_cluster), 3),
        "overlap_cutoff_a": round(cutoff, 3),
        "elements": sorted({metal.element for metal in sorted_cluster}),
        "keep": _metal_payload(keep) if keep else None,
        "drop": [_metal_payload(metal) for metal in drop],
        "drop_line_numbers": [metal.line_no for metal in drop],
        "members": [_metal_payload(metal) for metal in sorted_cluster],
    }


def _choose_keep_record(
    cluster: list[MetalRecord],
) -> tuple[MetalRecord | None, str, str]:
    same_atom_identity = {
        (
            metal.element,
            metal.resname,
            metal.chain,
            metal.resseq,
            metal.icode,
            metal.atom_name,
        )
        for metal in cluster
    }
    if len(same_atom_identity) == 1:
        return cluster[0], "identical_metal_atom_duplicate", "high"

    occupancy_choice = _highest_occupancy_choice(cluster)
    if occupancy_choice is not None:
        return occupancy_choice, "highest_occupancy_metal_site", "high"

    same_element = len({metal.element for metal in cluster}) == 1
    if same_element:
        b_factor_choice = _lowest_b_factor_choice(cluster)
        if b_factor_choice is not None:
            return b_factor_choice, "same_element_lower_b_factor_site", "medium"

    return (
        None,
        "ambiguous_overlapping_metal_identity_or_altloc",
        "low",
    )


def _highest_occupancy_choice(cluster: list[MetalRecord]) -> MetalRecord | None:
    with_occ = [metal for metal in cluster if metal.occupancy is not None]
    if len(with_occ) != len(cluster):
        return None
    ranked = sorted(
        with_occ,
        key=lambda metal: (float(metal.occupancy or 0.0), -metal.line_no),
        reverse=True,
    )
    if len(ranked) < 2:
        return ranked[0] if ranked else None
    top = float(ranked[0].occupancy or 0.0)
    second = float(ranked[1].occupancy or 0.0)
    if top >= second + OCCUPANCY_DECISION_MARGIN:
        return ranked[0]
    return None


def _lowest_b_factor_choice(cluster: list[MetalRecord]) -> MetalRecord | None:
    with_b = [metal for metal in cluster if metal.b_factor is not None]
    if len(with_b) != len(cluster):
        return None
    ranked = sorted(with_b, key=lambda metal: (float(metal.b_factor or 0.0), metal.line_no))
    if len(ranked) < 2:
        return ranked[0] if ranked else None
    best = float(ranked[0].b_factor or 0.0)
    second = float(ranked[1].b_factor or 0.0)
    if second >= best + B_FACTOR_DECISION_MARGIN:
        return ranked[0]
    return None


def _metal_payload(metal: MetalRecord | None) -> dict[str, object] | None:
    if metal is None:
        return None
    return {
        "line_no": metal.line_no,
        "serial": metal.serial,
        "element": metal.element,
        "resname": metal.resname,
        "chain": metal.chain,
        "resseq": metal.resseq,
        "icode": metal.icode,
        "atom_name": metal.atom_name,
        "altloc": metal.altloc,
        "occupancy": metal.occupancy if metal.occupancy is not None else "",
        "b_factor": metal.b_factor if metal.b_factor is not None else "",
        "xyz": [round(value, 3) for value in metal.xyz],
    }


def _write_without_lines(
    lines: list[str],
    output_pdb: Path,
    drop_line_numbers: set[int],
) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w", encoding="utf-8") as handle:
        for line_no, line in enumerate(lines, start=1):
            if line_no in drop_line_numbers:
                continue
            handle.write(line)


def _write_sidecar(
    sidecar_path: Path,
    input_pdb: Path,
    output_pdb: Path,
    summary: dict[str, object],
) -> None:
    payload = {
        "input_pdb": str(input_pdb),
        "output_pdb": str(output_pdb),
        **summary,
    }
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _line_element(line: str) -> str:
    token = line[76:78].strip().upper() if len(line) >= 78 else ""
    if token:
        return token
    atom_name = line[12:16].strip().upper()
    if len(atom_name) >= 2 and atom_name[:2] in COORDINATION_METALS:
        return atom_name[:2]
    return atom_name[:1]


def _line_float(line: str, start: int, end: int) -> float | None:
    try:
        token = line[start:end].strip()
        return float(token) if token else None
    except ValueError:
        return None


def _distance2(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> float:
    return sum((left[index] - right[index]) ** 2 for index in range(3))


def _max_pair_distance(cluster: list[MetalRecord]) -> float:
    distances: list[float] = []
    for left_index, left in enumerate(cluster):
        for right in cluster[left_index + 1 :]:
            distances.append(math.sqrt(_distance2(left.xyz, right.xyz)))
    return max(distances) if distances else 0.0


def _env_float(name: str, default: float) -> float:
    try:
        return float(str(os.environ.get(name, default)).strip())
    except (TypeError, ValueError):
        return default
