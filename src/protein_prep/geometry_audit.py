"""Compatibility wrapper for protein_prep.geometry.audit."""

from protein_prep.geometry.audit import (
    apply_conservative_geometry_fixes,
    audit_geometry,
)

__all__ = ["apply_conservative_geometry_fixes", "audit_geometry"]
