"""Narrow exception groups for optional protein-prep pipeline steps."""

from __future__ import annotations

import json
from typing import Final

TEXT_FILE_READ_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    UnicodeDecodeError,
)

JSON_LOAD_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    UnicodeDecodeError,
    json.JSONDecodeError,
)

# Chain prune, water centroid prefilter, element post-checks: best-effort
OPTIONAL_PIPELINE_STEP_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,
    ValueError,
    RuntimeError,
    TypeError,
    KeyError,
    AttributeError,
    ZeroDivisionError,
)
