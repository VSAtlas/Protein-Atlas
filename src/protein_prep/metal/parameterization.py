"""Metal parameterization staging API.

This remains review-only for default docking prep; it is exposed here so MD
and MMGBSA workflows can depend on a dedicated metal namespace.
"""

from protein_prep.metal_parameterization import (
    build_mcpb_parameterization_plan,
    summarize_parameterization_plan,
    write_parameterization_plan,
)

__all__ = [
    "build_mcpb_parameterization_plan",
    "summarize_parameterization_plan",
    "write_parameterization_plan",
]
