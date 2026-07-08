"""Small scalar helpers shared by reporting modules.

Boolean coercion contract (reporting stack):
- ``parse_bool``: Config / ``config.txt`` values and other strings that may be
  quote-wrapped. Uses ``strip_quotes``, tri-state with ``default``, and explicit
  token sets (includes ``1.0`` / ``0.0``).
- ``as_report_row_bool``: Heterogeneous master-row / CSV cells: native ``bool``;
  ``int`` / ``float`` via ``bool(x)`` (same as Python, so ``NaN`` is truthy);
  strings only match ``("1", "true", "yes", "on")`` after strip (no quote layer).
- ``is_truthy_matrix_cell``: Heatmap / parquet matrix flags — ``bool``; numeric
  non-zero with ``NaN`` → false; strings in ``MATRIX_BOOL_TRUE_STRINGS``.
- ``truthy_string_mask``: Lowercased stripped string membership for vectorized
  publication columns (same tokens as ``parse_bool`` true side, including
  ``1.0``).
"""

from __future__ import annotations

import math
from typing import Any, Iterable, Iterator, cast

# True/false after lowercasing; used by parse_bool (with strip_quotes on input).
CONFIG_BOOL_TRUE_STRINGS = frozenset({"1", "1.0", "true", "yes", "y", "on"})
CONFIG_BOOL_FALSE_STRINGS = frozenset({"0", "0.0", "false", "no", "n", "off"})

# Heatmap matrix and publication vectors reuse the same true-token set as ``parse_bool``.
MATRIX_BOOL_TRUE_STRINGS = CONFIG_BOOL_TRUE_STRINGS
TRUTHY_PUBLICATION_STRINGS = CONFIG_BOOL_TRUE_STRINGS

# Legacy master-row string tokens (no "y", "1.0", etc.) for as_report_row_bool.
_REPORT_ROW_BOOL_TRUE_STRINGS = frozenset({"1", "true", "yes", "on"})


def normalize_text(value: Any) -> str:
    return str(value or "").strip()


def normalize_side_effect_label(value: Any) -> str:
    """Collapse underscores and repeated whitespace (side-effect string cleanup)."""
    return " ".join(str(value or "").replace("_", " ").split()).strip()


def _iter_normalized_side_effect_labels(values: object) -> Iterator[str]:
    # str/bytes/dict: exclude so we do not iterate strings character-by-character.
    if values is None:
        return
    if isinstance(values, (str, bytes, dict)):
        return
    if not hasattr(values, "__iter__"):
        return
    for value in cast(Iterable[object], values):
        clean = normalize_side_effect_label(value)
        if clean:
            yield clean


def normalized_side_effect_list(values: object) -> list[str]:
    """Non-empty `normalize_side_effect_label` outputs, preserving source order."""
    return list(_iter_normalized_side_effect_labels(values))


def normalized_side_effect_set(values: object) -> set[str]:
    """Deduplicated non-empty normalized side-effect-ish strings from a JSON-ish list."""
    return set(_iter_normalized_side_effect_labels(values))


def normalize_whitespace_runs(value: Any) -> str:
    """Collapse internal whitespace runs while trimming ends."""
    return " ".join(str(value or "").split()).strip()


def strip_quotes(value: Any) -> str:
    stripped = normalize_text(value)
    if len(stripped) >= 2 and stripped[0] == stripped[-1] and stripped[0] in ("'", '"'):
        return stripped[1:-1].strip()
    return stripped


def clean_report_text(value: Any) -> str:
    """Trim; treat empty and case-insensitive ``nan`` as empty (reporting cells)."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def as_float(value: Any) -> float | None:
    try:
        if value is None:
            return None
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def parse_bool(value: Any, default: bool | None = None) -> bool | None:
    if isinstance(value, bool):
        return value
    text = strip_quotes(value).lower()
    if not text:
        return default
    if text in CONFIG_BOOL_TRUE_STRINGS:
        return True
    if text in CONFIG_BOOL_FALSE_STRINGS:
        return False
    return default


def parse_bool_or(value: Any, *, default: bool) -> bool:
    """Like ``parse_bool`` but never returns None (unknown → ``default``)."""
    parsed = parse_bool(value, default=None)
    return default if parsed is None else parsed


def as_report_row_bool(x: Any) -> bool:
    """Infer boolean for heterogeneous master CSV / row dict cells (legacy semantics)."""
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    s = str(x).strip().lower()
    return s in _REPORT_ROW_BOOL_TRUE_STRINGS


def is_truthy_matrix_cell(value: Any) -> bool:
    """Truthy for heatmap matrix / promiscuity cells (numeric + string rules)."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            numeric = 0.0
        if math.isnan(numeric):
            return False
        return numeric != 0.0
    return normalize_text(value).lower() in MATRIX_BOOL_TRUE_STRINGS


def median_floats(values: Iterable[float]) -> float | None:
    """Median of finite floats; empty input → ``None`` (even-count average of middles)."""
    materialized = [float(v) for v in values]
    if not materialized:
        return None
    sorted_values = sorted(materialized)
    mid = len(sorted_values) // 2
    if len(sorted_values) % 2 == 1:
        return sorted_values[mid]
    return (sorted_values[mid - 1] + sorted_values[mid]) / 2.0
