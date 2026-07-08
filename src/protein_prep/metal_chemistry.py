"""Heuristic metal-site chemistry typing for audit and routing.

These helpers make receptor-prep reports more explicit about metal charge and
geometry assumptions. They do not parameterize metal centers for docking or MD.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Mapping, Sequence

Coord3 = tuple[float, float, float]

COMMON_FORMAL_CHARGES: dict[str, tuple[int, ...]] = {
    "AG": (1,),
    "AL": (3,),
    "BA": (2,),
    "CA": (2,),
    "CD": (2,),
    "CO": (2, 3),
    "CS": (1,),
    "CU": (1, 2),
    "FE": (2, 3),
    "HG": (2,),
    "K": (1,),
    "LI": (1,),
    "MG": (2,),
    "MN": (2, 3),
    "NA": (1,),
    "NI": (2,),
    "RB": (1,),
    "SR": (2,),
    "ZN": (2,),
}

METALS_WITH_DIRECTIONAL_AUTODOCK_MODEL = {"ZN"}


def infer_metal_formal_charge(element: str) -> dict[str, object]:
    """Return common oxidation-state guesses for an elemental metal token."""

    token = element.strip().upper()
    states = COMMON_FORMAL_CHARGES.get(token, ())
    if not states:
        return {
            "formal_charge_status": "unknown",
            "formal_charge_states": "",
            "formal_charge_primary": "",
            "formal_charge_note": "no_common_charge_rule",
        }
    primary = states[0]
    return {
        "formal_charge_status": "common_single_state"
        if len(states) == 1
        else "ambiguous_common_states",
        "formal_charge_states": ",".join(_signed_charge(state) for state in states),
        "formal_charge_primary": _signed_charge(primary),
        "formal_charge_note": "heuristic_common_oxidation_state",
    }


def classify_coordination_geometry(
    metal_xyz: Coord3 | None,
    donor_xyzs: Sequence[Coord3],
) -> dict[str, object]:
    """Classify first-shell coordination geometry using donor count and angles."""

    donor_count = len(donor_xyzs)
    if metal_xyz is None or donor_count < 2:
        return _geometry_result(donor_count, "undercoordinated_or_unknown", "", "")
    angles = _coordination_angles(metal_xyz, donor_xyzs)
    if donor_count == 2:
        return _linear_geometry(donor_count, angles)
    if donor_count == 3:
        return _geometry_result(donor_count, "trigonal_or_pyramidal", angles, "low")
    if donor_count == 4:
        return _four_coordinate_geometry(donor_count, angles)
    if donor_count == 5:
        return _geometry_result(
            donor_count,
            "square_pyramidal_or_trigonal_bipyramidal",
            angles,
            "low",
        )
    if donor_count == 6:
        return _geometry_result(donor_count, "octahedral_like", angles, "medium")
    return _geometry_result(donor_count, "high_coordination_or_multisite", angles, "low")


def metal_treatment_recommendation(
    element: str,
    geometry: str,
    donor_count: int,
    donor_categories: Sequence[str],
) -> dict[str, object]:
    """Choose an auditable treatment route without pretending PDBQT solves metals."""

    token = element.strip().upper()
    categories = {category.strip().lower() for category in donor_categories}
    ligand_bound = "ligand" in categories
    if token in METALS_WITH_DIRECTIONAL_AUTODOCK_MODEL and donor_count in {3, 4}:
        return {
            "metal_treatment": "autodock4zn_candidate",
            "metal_treatment_note": (
                "Zn has an established AutoDock4Zn/TZ pseudoatom route; use only "
                "with AD4Zn maps/parameters, not as generic Vina PDBQT chemistry."
            ),
            "requires_parameterization": True,
            "supports_pdbqt_only": False,
            "ligand_bound": ligand_bound,
        }
    if token == "ZN":
        return {
            "metal_treatment": "preserve_audit_zn_nonideal",
            "metal_treatment_note": "Zn site not in the simple AutoDock4Zn candidate geometry/count window.",
            "requires_parameterization": True,
            "supports_pdbqt_only": False,
            "ligand_bound": ligand_bound,
        }
    if token in COMMON_FORMAL_CHARGES:
        return {
            "metal_treatment": "preserve_audit_and_parameterize_if_scored",
            "metal_treatment_note": (
                f"{token} has no Atlas directional docking model; preserve the ion, "
                "audit geometry, and use AmberTools/MMGBSA-style parameterization "
                "or curated parameters for minimization/scoring."
            ),
            "requires_parameterization": True,
            "supports_pdbqt_only": False,
            "ligand_bound": ligand_bound,
        }
    return {
        "metal_treatment": "unknown_metal_review_required",
        "metal_treatment_note": "Unknown metal token; manual chemistry review required.",
        "requires_parameterization": True,
        "supports_pdbqt_only": False,
        "ligand_bound": ligand_bound,
    }


def summarize_metal_chemistry(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    """Summarize per-metal chemistry audit rows for benchmark CSV output."""

    if not rows:
        return {
            "metal_chemistry_status": "no_metals",
            "metal_formal_charge_states": "",
            "metal_coordination_geometry_counts": "",
            "metal_treatment_counts": "",
            "metal_ad4zn_candidate_count": 0,
            "metal_parameterization_required_count": 0,
            "metal_bound_protonation_status": "no_metals",
            "metal_bound_protonation_review_count": 0,
            "metal_bound_protonation_unresolved_count": 0,
            "metal_bound_protonation_evidence_count": 0,
        }
    geometry_counts = _row_counter(rows, "coordination_geometry")
    treatment_counts = _row_counter(rows, "metal_treatment")
    return {
        "metal_chemistry_status": "heuristic_audit_only",
        "metal_formal_charge_states": ",".join(_formal_charge_states(rows)),
        "metal_coordination_geometry_counts": _format_counter(geometry_counts),
        "metal_treatment_counts": _format_counter(treatment_counts),
        "metal_ad4zn_candidate_count": _count_value(
            rows,
            "metal_treatment",
            "autodock4zn_candidate",
        ),
        "metal_parameterization_required_count": _count_truthy(
            rows,
            "requires_parameterization",
        ),
        "metal_bound_protonation_status": _summarize_protonation_status(rows),
        "metal_bound_protonation_review_count": _sum_int(
            rows,
            "metal_bound_protonation_unresolved_count",
        ),
        "metal_bound_protonation_unresolved_count": _sum_int(
            rows,
            "metal_bound_protonation_unresolved_count",
        ),
        "metal_bound_protonation_evidence_count": _sum_int(
            rows,
            "metal_bound_protonation_evidence_count",
        ),
    }


def audit_metal_chemistry(
    *,
    element: str,
    metal_xyz: Coord3 | None,
    donors: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Build one heuristic metal chemistry audit row."""

    donor_coords: list[Coord3] = []
    for donor in donors:
        coords = _coords_from_donor(donor)
        if coords is not None:
            donor_coords.append(coords)
    donor_categories = [str(donor.get("category", "")) for donor in donors]
    charge = infer_metal_formal_charge(element)
    geometry = classify_coordination_geometry(metal_xyz, donor_coords)
    treatment = metal_treatment_recommendation(
        element,
        str(geometry.get("coordination_geometry", "")),
        int(str(geometry.get("coordination_number", 0) or 0)),
        donor_categories,
    )
    protonation = metal_bound_protonation_recommendations(element, donors)
    return {
        "element": element.strip().upper(),
        **charge,
        **geometry,
        **treatment,
        **protonation,
    }


def metal_bound_protonation_recommendations(
    element: str,
    donors: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Flag metal-bound residues whose protonation state must be curated."""

    recommendations: list[dict[str, object]] = []
    for donor in donors:
        if str(donor.get("category", "")).lower() != "protein":
            continue
        rec = _protein_donor_protonation_recommendation(element, donor)
        if rec:
            recommendations.append(rec)
    recommendations = _dedupe_protonation_recommendations(recommendations)
    unresolved = [
        row for row in recommendations if bool(row.get("review_required", True))
    ]
    evidence = [
        row for row in recommendations if not bool(row.get("review_required", True))
    ]
    if not recommendations:
        return {
            "metal_bound_protonation_status": "no_metal_bound_titratable_residues",
            "metal_bound_protonation_review_count": 0,
            "metal_bound_protonation_unresolved_count": 0,
            "metal_bound_protonation_evidence_count": 0,
            "metal_bound_protonation_recommendations": [],
        }
    if not unresolved:
        return {
            "metal_bound_protonation_status": "resolved_by_standard_residue_names",
            "metal_bound_protonation_review_count": 0,
            "metal_bound_protonation_unresolved_count": 0,
            "metal_bound_protonation_evidence_count": len(evidence),
            "metal_bound_protonation_recommendations": recommendations,
        }
    return {
        "metal_bound_protonation_status": "review_required",
        "metal_bound_protonation_review_count": len(unresolved),
        "metal_bound_protonation_unresolved_count": len(unresolved),
        "metal_bound_protonation_evidence_count": len(evidence),
        "metal_bound_protonation_recommendations": recommendations,
    }


def _coords_from_donor(donor: Mapping[str, object]) -> Coord3 | None:
    raw = donor.get("coords")
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        return None
    try:
        x, y, z = raw
        return float(x), float(y), float(z)
    except Exception:
        return None


def _protein_donor_protonation_recommendation(
    element: str,
    donor: Mapping[str, object],
) -> dict[str, object] | None:
    resname = str(donor.get("resname", "")).upper()
    atom_name = str(donor.get("atom_name", "")).upper()
    donor_element = str(donor.get("element", "")).upper()
    residue_id = _donor_residue_id(resname, atom_name, donor)
    if resname in {"HIS", "HID", "HIE", "HIP"} or atom_name in {"ND1", "NE2"}:
        return _histidine_protonation_row(residue_id, atom_name)
    if resname in {"CYS", "CYM"} and donor_element == "S":
        return _cysteine_protonation_row(residue_id, element)
    if resname in {"ASP", "ASH", "GLU", "GLH"} and donor_element == "O":
        return _carboxylate_protonation_row(residue_id, resname)
    if resname in {"TYR", "LYS"}:
        return _manual_titration_row(residue_id, resname)
    return None


def _donor_residue_id(
    resname: str,
    atom_name: str,
    donor: Mapping[str, object],
) -> str:
    return ":".join(
        part
        for part in (
            resname,
            str(donor.get("chain", "")),
            str(donor.get("resseq", "")),
            atom_name,
        )
        if part
    )


def _histidine_protonation_row(
    residue_id: str,
    atom_name: str,
) -> dict[str, object]:
    suggested = "HIE" if atom_name == "ND1" else "HID" if atom_name == "NE2" else ""
    return {
        "residue": residue_id,
        "recommendation": "select_histidine_tautomer_with_coordinating_nitrogen_neutral",
        "suggested_amber_resname": suggested,
        "review_required": True,
        "reason": (
            "Metal-bound histidine should generally donate through the coordinating "
            "nitrogen without forcing HIP unless the local electrostatics justify it."
        ),
    }


def _cysteine_protonation_row(
    residue_id: str,
    element: str,
) -> dict[str, object]:
    return {
        "residue": residue_id,
        "recommendation": "review_deprotonated_cysteinate_and_bonded_metal_model",
        "suggested_amber_resname": "CYM",
        "review_required": True,
        "reason": (
            f"{element.strip().upper()}-bound cysteine sulfur often requires CYM-like "
            "naming and metal-site parameters rather than ordinary neutral CYS."
        ),
    }


def _carboxylate_protonation_row(
    residue_id: str,
    resname: str,
) -> dict[str, object]:
    suggested = "ASP" if resname in {"ASP", "ASH"} else "GLU"
    if resname in {"ASP", "GLU"}:
        return {
            "residue": residue_id,
            "recommendation": "standard_deprotonated_carboxylate_retained",
            "suggested_amber_resname": suggested,
            "review_required": False,
            "reason": (
                "Metal-bound Asp/Glu donor is already named as a charged "
                "carboxylate in the prepared receptor."
            ),
        }
    return {
        "residue": residue_id,
        "recommendation": "prefer_deprotonated_carboxylate_unless_validated_otherwise",
        "suggested_amber_resname": suggested,
        "review_required": True,
        "reason": "Metal-bound Asp/Glu oxygen donors usually need charged carboxylate review.",
    }


def _manual_titration_row(residue_id: str, resname: str) -> dict[str, object]:
    return {
        "residue": residue_id,
        "recommendation": "manual_metal_bound_titration_review",
        "review_required": True,
        "reason": f"{resname} metal donation is context-sensitive and not solved by PROPKA alone.",
    }


def _summarize_protonation_status(rows: Sequence[Mapping[str, object]]) -> str:
    if any(_row_int(row, "metal_bound_protonation_unresolved_count") > 0 for row in rows):
        return "review_required"
    if any(_row_int(row, "metal_bound_protonation_evidence_count") > 0 for row in rows):
        return "resolved_by_standard_residue_names"
    return "no_metal_bound_titratable_residues"


def _row_int(row: Mapping[str, object], key: str) -> int:
    return _safe_int(row.get(key, 0))


def _dedupe_protonation_recommendations(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    deduped: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        residue = str(row.get("residue", ""))
        residue_key = ":".join(residue.split(":")[:3])
        recommendation = str(row.get("recommendation", ""))
        key = (residue_key, recommendation)
        deduped.setdefault(key, dict(row) | {"residue": residue_key})
    return list(deduped.values())


def _row_counter(rows: Sequence[Mapping[str, object]], key: str) -> Counter[str]:
    return Counter(str(row.get(key, "")) for row in rows)


def _formal_charge_states(rows: Sequence[Mapping[str, object]]) -> list[str]:
    return sorted(
        {
            str(row.get("formal_charge_primary", ""))
            for row in rows
            if str(row.get("formal_charge_primary", ""))
        }
    )


def _count_value(
    rows: Sequence[Mapping[str, object]],
    key: str,
    value: str,
) -> int:
    return sum(1 for row in rows if row.get(key) == value)


def _count_truthy(rows: Sequence[Mapping[str, object]], key: str) -> int:
    return sum(1 for row in rows if bool(row.get(key)))


def _sum_int(rows: Sequence[Mapping[str, object]], key: str) -> int:
    total = 0
    for row in rows:
        total += _safe_int(row.get(key, 0))
    return total


def _safe_int(value: object) -> int:
    try:
        return int(str(value or 0))
    except Exception:
        return 0


def _coordination_angles(metal_xyz: Coord3, donor_xyzs: Sequence[Coord3]) -> list[float]:
    vectors: list[Coord3] = []
    for donor_xyz in donor_xyzs:
        vector = _unit_vector(_subtract(donor_xyz, metal_xyz))
        if vector is not None:
            vectors.append(vector)
    angles: list[float] = []
    for idx, left in enumerate(vectors):
        for right in vectors[idx + 1 :]:
            dot = max(-1.0, min(1.0, _dot(left, right)))
            angles.append(round(math.degrees(math.acos(dot)), 1))
    return angles


def _linear_geometry(donor_count: int, angles: Sequence[float]) -> dict[str, object]:
    max_angle = max(angles) if angles else 0.0
    geometry = "linear_like" if max_angle >= 150.0 else "bent"
    return _geometry_result(donor_count, geometry, angles, "medium")


def _four_coordinate_geometry(
    donor_count: int,
    angles: Sequence[float],
) -> dict[str, object]:
    max_angle = max(angles) if angles else 0.0
    tetra_dev = _mean_abs_deviation(angles, 109.5)
    if max_angle >= 150.0:
        return _geometry_result(donor_count, "square_planar_or_distorted", angles, "medium")
    if tetra_dev <= 25.0:
        return _geometry_result(donor_count, "tetrahedral_like", angles, "medium")
    return _geometry_result(donor_count, "four_coordinate_distorted", angles, "low")


def _geometry_result(
    donor_count: int,
    geometry: str,
    angles: Sequence[float] | str,
    confidence: str,
) -> dict[str, object]:
    angle_text = angles if isinstance(angles, str) else ",".join(f"{angle:.1f}" for angle in angles)
    return {
        "coordination_number": donor_count,
        "coordination_geometry": geometry,
        "coordination_geometry_confidence": confidence,
        "coordination_angles_deg": angle_text,
    }


def _mean_abs_deviation(values: Sequence[float], target: float) -> float:
    if not values:
        return float("inf")
    return sum(abs(value - target) for value in values) / len(values)


def _subtract(left: Coord3, right: Coord3) -> Coord3:
    return left[0] - right[0], left[1] - right[1], left[2] - right[2]


def _unit_vector(vector: Coord3) -> Coord3 | None:
    norm = math.sqrt(_dot(vector, vector))
    if norm <= 1e-8:
        return None
    return vector[0] / norm, vector[1] / norm, vector[2] / norm


def _dot(left: Coord3, right: Coord3) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _signed_charge(charge: int) -> str:
    return f"+{charge}" if charge > 0 else str(charge)


def _format_counter(counter: Counter[str]) -> str:
    return ";".join(f"{key}:{counter[key]}" for key in sorted(counter) if key)
