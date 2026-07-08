"""Narrow exception groups for MMGBSA and post-docking I/O.

Broad ``except Exception`` hides programming errors. These tuples capture the
expected failure modes for file reads, config/CSV parsing, subprocess edges,
and RDKit-wrapped Boost errors (often ``RuntimeError`` / ``ValueError``).
"""

from __future__ import annotations

import json
import statistics
import subprocess
from typing import Final

# Plain-text and encoded file reads (config, CSV, MMPBSA outputs)
TEXT_FILE_READ_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    UnicodeDecodeError,
)

JSON_DECODE_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    UnicodeDecodeError,
    json.JSONDecodeError,
)

FLOAT_COERCE_ERRORS: Final[tuple[type[BaseException], ...]] = (ValueError, TypeError)

INT_COERCE_ERRORS: Final[tuple[type[BaseException], ...]] = (ValueError, TypeError)

# statistics.mean/median/pstdev on edge-case inputs
STAT_AGG_ERRORS: Final[tuple[type[BaseException], ...]] = (
    statistics.StatisticsError,
    ValueError,
    TypeError,
)

# Optional third-party imports (RDKit, extension modules)
OPTIONAL_IMPORT_ERRORS: Final[tuple[type[BaseException], ...]] = (
    ImportError,
    ModuleNotFoundError,
)

PATH_RESOLVE_ERRORS: Final[tuple[type[BaseException], ...]] = (OSError, RuntimeError)

# SDMolSupplier / Mol readers: filesystem + RDKit internal failures
RDKIT_SUPPLIER_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    RuntimeError,
    ValueError,
)

# Atom/getter property paths (conformer positions, formal charge, etc.)
RDKIT_PROPERTY_ERRORS: Final[tuple[type[BaseException], ...]] = (
    RuntimeError,
    ValueError,
    TypeError,
    AttributeError,
)

SUBPROCESS_RUN_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    subprocess.SubprocessError,
    PermissionError,
)

# CLI entrypoints: user-facing failures we translate to exit code 1
MMGBSA_CLI_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    ValueError,
    RuntimeError,
    subprocess.SubprocessError,
    PermissionError,
    FileNotFoundError,
    MemoryError,
)

# End-to-end pose SDF export: I/O + RDKit + JSON sidecars
MMGBSA_POSE_EXPORT_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    UnicodeDecodeError,
    RuntimeError,
    ValueError,
    TypeError,
    AttributeError,
    ImportError,
    json.JSONDecodeError,
)

# Best-effort QC sidecar writes (must not abort MMGBSA aggregation)
OPTIONAL_QC_STEP_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    ValueError,
    RuntimeError,
    TypeError,
    PermissionError,
    KeyError,
)
