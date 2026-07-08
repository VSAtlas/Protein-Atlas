"""Publication-readiness decisions for retained metal sites."""

from __future__ import annotations

from typing import Mapping


def classify_metal_publication_policy(row: Mapping[str, object]) -> dict[str, object]:
    """Separate default docking readiness from metal-aware physics requirements."""

    if _row_int(row, "metals_baseline") <= 0:
        return _result("no_metals", False, "", "")
    if _row_int(row, "metals_baseline") < _row_int(row, "metals_input"):
        return _result("metal_loss", True, "metal_count_decreased", "")
    if not _truthy(row.get("metal_coord_retained")):
        return _result("coordination_not_retained", True, "metal_coordination_not_retained", "")
    unresolved_protonation = _row_int(
        row,
        "metal_bound_protonation_unresolved_count",
        fallback_key="metal_bound_protonation_review_count",
    )
    if unresolved_protonation > 0:
        return _result(
            "unresolved_protonation",
            True,
            "metal_bound_protonation_review",
            "",
        )
    if _row_int(row, "metal_bound_competitive_ligand_retained_count") > 0:
        return _result(
            "competitive_ligand_retained",
            True,
            "metal_bound_competitive_ligand_retained_in_receptor",
            "",
        )
    if _row_int(row, "metal_bound_het_expected_removed_donor_count") > 0:
        return _result(
            "ligand_metal_coordination_removed",
            False,
            "",
            "metal_bound_competitive_ligand_removed_for_redocking",
        )
    if _row_int(row, "metal_bound_het_unexpected_removed_count") > 0:
        return _result(
            "unexpected_structural_metal_contact_loss",
            True,
            "unexpected_metal_bound_het_removed",
            "",
        )
    if _row_int(row, "metal_parameterization_required_count") > 0:
        return _result(
            "retained_structural_metal_audit",
            False,
            "",
            "metal_retained_with_charge_geometry_coordination_audit_not_forcefield",
        )
    return _result("metal_ready", False, "", "metal_retained_no_extra_review")


def _result(
    route: str,
    review_required: bool,
    failure: str,
    disclosure: str,
) -> dict[str, object]:
    return {
        "metal_publication_route": route,
        "metal_publication_review_required": review_required,
        "metal_publication_failure": failure,
        "metal_publication_disclosure": disclosure,
    }


def _row_int(
    row: Mapping[str, object],
    key: str,
    *,
    fallback_key: str | None = None,
) -> int:
    raw = row.get(key, row.get(fallback_key, 0) if fallback_key else 0)
    try:
        return int(float(str(raw or 0)))
    except Exception:
        return 0


def _truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "ok"}
