from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation, localcontext
from typing import Any

import pandas as pd


SPD_ACTIVITY_PARSER_VERSION = "spd_activity_interval_v3"
SPD_LABEL_POSITIVE = 1
SPD_LABEL_NEGATIVE = 0
SPD_LABEL_UNKNOWN: None = None
SPD_LABEL_EXCLUDED_CONFLICT = -1

_SUPPORTED_RELATIONS = frozenset({"=", ">", ">=", "<", "<="})
_UNIT_ALIASES: dict[str, tuple[str, Decimal]] = {
    "M": ("M", Decimal("1000000")),
    "molar": ("M", Decimal("1000000")),
    "mm": ("mM", Decimal("1000")),
    "millimolar": ("mM", Decimal("1000")),
    "um": ("uM", Decimal("1")),
    "µm": ("uM", Decimal("1")),
    "μm": ("uM", Decimal("1")),
    "micromolar": ("uM", Decimal("1")),
    "nm": ("nM", Decimal("0.001")),
    "nanomolar": ("nM", Decimal("0.001")),
    "pm": ("pM", Decimal("0.000001")),
    "picomolar": ("pM", Decimal("0.000001")),
}


@dataclass(frozen=True)
class ActivityInterval:
    activity_value_uM: float | None
    activity_lower_bound_uM: float | None
    activity_upper_bound_uM: float | None
    lower_inclusive: bool | None
    upper_inclusive: bool | None
    normalized_relation: str | None
    normalized_unit: str | None
    activity_parse_status: str
    activity_parse_reason: str
    _activity_value_uM_decimal: Decimal | None = field(
        default=None, repr=False, compare=False
    )
    _activity_lower_bound_uM_decimal: Decimal | None = field(
        default=None, repr=False, compare=False
    )
    _activity_upper_bound_uM_decimal: Decimal | None = field(
        default=None, repr=False, compare=False
    )

    @property
    def is_valid(self) -> bool:
        return self.activity_parse_status == "valid"


@dataclass(frozen=True)
class BindingLabelResult:
    numeric_label: int | None
    binding_status: str
    binding_reason: str


@dataclass(frozen=True)
class ExposureLabelResult:
    exposure_margin_lower: float | None
    exposure_margin_upper: float | None
    margin_lower_inclusive: bool | None
    margin_upper_inclusive: bool | None
    numeric_label: int | None
    exposure_status: str
    exposure_reason: str


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    try:
        return bool(missing)
    except (TypeError, ValueError):
        return False


def _parse_finite_decimal(value: Any) -> tuple[Decimal | None, str | None]:
    if _is_missing(value):
        return None, "missing"
    if isinstance(value, bool):
        return None, "malformed"
    if isinstance(value, str):
        value = value.strip()
        if not value:
            return None, "missing"
    try:
        parsed = Decimal(str(value))
    except (InvalidOperation, TypeError, ValueError):
        return None, "malformed"
    if not parsed.is_finite():
        return None, "nonfinite"
    return parsed, None


def _exact_decimal_product(left: Decimal, right: Decimal) -> Decimal:
    """Multiply finite source decimals without context-induced rounding."""

    precision = max(
        28,
        len(left.as_tuple().digits) + len(right.as_tuple().digits) + 4,
    )
    with localcontext() as context:
        context.prec = precision
        return left * right


def _reporting_decimal_quotient(
    numerator: Decimal,
    denominator: Decimal,
) -> Decimal:
    """Return a high-precision quotient used only for public float reporting."""

    precision = max(
        28,
        len(numerator.as_tuple().digits)
        + len(denominator.as_tuple().digits)
        + 10,
    )
    with localcontext() as context:
        context.prec = precision
        return numerator / denominator


def _exact_interval_bound(
    exact: Decimal | None,
    public: float | None,
) -> Decimal | None:
    if exact is not None:
        return exact
    if public is None:
        return None
    return Decimal(str(public))


def _normalize_relation(relation: Any) -> tuple[str | None, str | None]:
    if _is_missing(relation):
        return None, "missing_relation"
    normalized = str(relation).strip()
    if not normalized:
        return None, "blank_relation"
    if normalized not in _SUPPORTED_RELATIONS:
        return None, "unsupported_relation"
    return normalized, None


def _normalize_unit(
    unit: Any,
) -> tuple[str | None, Decimal | None, str | None]:
    if _is_missing(unit):
        return None, None, "missing_activity_unit"
    token = str(unit).strip()
    if not token:
        return None, None, "missing_activity_unit"
    lookup = token if token == "M" else token.lower()
    normalized = _UNIT_ALIASES.get(lookup)
    if normalized is None:
        return None, None, "unsupported_activity_unit"
    return normalized[0], normalized[1], None


def _invalid_interval(
    *,
    value_uM: float | None,
    value_uM_decimal: Decimal | None,
    relation: str | None,
    unit: str | None,
    reason: str,
) -> ActivityInterval:
    return ActivityInterval(
        activity_value_uM=value_uM,
        activity_lower_bound_uM=None,
        activity_upper_bound_uM=None,
        lower_inclusive=None,
        upper_inclusive=None,
        normalized_relation=relation,
        normalized_unit=unit,
        activity_parse_status="invalid",
        activity_parse_reason=reason,
        _activity_value_uM_decimal=value_uM_decimal,
        _activity_lower_bound_uM_decimal=None,
        _activity_upper_bound_uM_decimal=None,
    )


def parse_spd_activity_interval(
    activity_value: Any,
    activity_relation: Any,
    activity_unit: Any,
) -> ActivityInterval:
    """Parse one SPD activity observation into a closed/open feasible interval.

    Only the five explicit relation symbols are recognized. Blank, missing, or
    approximate aliases fail closed instead of being interpreted as exact.
    """

    relation, relation_error = _normalize_relation(activity_relation)
    numeric_value, numeric_error = _parse_finite_decimal(activity_value)
    normalized_unit, unit_factor, unit_error = _normalize_unit(activity_unit)
    value_uM_decimal = (
        _exact_decimal_product(numeric_value, unit_factor)
        if numeric_value is not None and unit_factor is not None
        else None
    )
    value_uM = (
        float(value_uM_decimal) if value_uM_decimal is not None else None
    )

    if relation_error is not None:
        return _invalid_interval(
            value_uM=value_uM,
            value_uM_decimal=value_uM_decimal,
            relation=relation,
            unit=normalized_unit,
            reason=relation_error,
        )
    if numeric_error is not None:
        return _invalid_interval(
            value_uM=None,
            value_uM_decimal=None,
            relation=relation,
            unit=normalized_unit,
            reason=f"{numeric_error}_activity_value",
        )
    if unit_error is not None:
        return _invalid_interval(
            value_uM=None,
            value_uM_decimal=None,
            relation=relation,
            unit=normalized_unit,
            reason=unit_error,
        )
    if value_uM is None or value_uM_decimal is None:
        raise AssertionError("validated activity value unexpectedly missing")
    if value_uM_decimal < 0:
        return _invalid_interval(
            value_uM=value_uM,
            value_uM_decimal=value_uM_decimal,
            relation=relation,
            unit=normalized_unit,
            reason="negative_activity_value",
        )

    if relation == "=":
        lower_decimal, upper_decimal = value_uM_decimal, value_uM_decimal
        lower_inclusive, upper_inclusive = True, True
        reason = "exact_activity"
    elif relation == ">":
        lower_decimal, upper_decimal = value_uM_decimal, Decimal("Infinity")
        lower_inclusive, upper_inclusive = False, False
        reason = "open_lower_bound"
    elif relation == ">=":
        lower_decimal, upper_decimal = value_uM_decimal, Decimal("Infinity")
        lower_inclusive, upper_inclusive = True, False
        reason = "closed_lower_bound"
    elif relation == "<":
        lower_decimal, upper_decimal = Decimal("0"), value_uM_decimal
        lower_inclusive, upper_inclusive = True, False
        reason = "open_upper_bound"
    else:
        lower_decimal, upper_decimal = Decimal("0"), value_uM_decimal
        lower_inclusive, upper_inclusive = True, True
        reason = "closed_upper_bound"

    if lower_decimal == upper_decimal and not (
        lower_inclusive and upper_inclusive
    ):
        return _invalid_interval(
            value_uM=value_uM,
            value_uM_decimal=value_uM_decimal,
            relation=relation,
            unit=normalized_unit,
            reason="empty_activity_interval",
        )
    lower = float(lower_decimal)
    upper = float(upper_decimal)
    return ActivityInterval(
        activity_value_uM=value_uM,
        activity_lower_bound_uM=lower,
        activity_upper_bound_uM=upper,
        lower_inclusive=lower_inclusive,
        upper_inclusive=upper_inclusive,
        normalized_relation=relation,
        normalized_unit=normalized_unit,
        activity_parse_status="valid",
        activity_parse_reason=reason,
        _activity_value_uM_decimal=value_uM_decimal,
        _activity_lower_bound_uM_decimal=lower_decimal,
        _activity_upper_bound_uM_decimal=upper_decimal,
    )


def label_spd_binding_observation(
    interval: ActivityInterval,
    *,
    active_um: float = 1.0,
    inactive_um: float = 10.0,
) -> BindingLabelResult:
    """Assign the canonical SPD binding label from the complete feasible interval."""

    active_decimal, active_error = _parse_finite_decimal(active_um)
    inactive_decimal, inactive_error = _parse_finite_decimal(inactive_um)
    if (
        active_error is not None
        or inactive_error is not None
        or active_decimal is None
        or inactive_decimal is None
        or active_decimal < 0
        or inactive_decimal <= active_decimal
    ):
        raise ValueError("binding thresholds require 0 <= active_um < inactive_um")
    if not interval.is_valid:
        return BindingLabelResult(
            numeric_label=None,
            binding_status="unknown_invalid_activity",
            binding_reason=interval.activity_parse_reason,
        )
    lower = _exact_interval_bound(
        interval._activity_lower_bound_uM_decimal,
        interval.activity_lower_bound_uM,
    )
    upper = _exact_interval_bound(
        interval._activity_upper_bound_uM_decimal,
        interval.activity_upper_bound_uM,
    )
    if lower is None or upper is None:
        raise AssertionError("valid activity interval is missing bounds")
    if upper <= active_decimal:
        return BindingLabelResult(
            SPD_LABEL_POSITIVE,
            "positive",
            "interval_definitely_at_or_below_active_threshold",
        )
    if lower >= inactive_decimal:
        return BindingLabelResult(
            SPD_LABEL_NEGATIVE,
            "negative",
            "interval_definitely_at_or_above_inactive_threshold",
        )
    return BindingLabelResult(
        SPD_LABEL_UNKNOWN,
        "unknown",
        "activity_interval_crosses_or_lies_between_thresholds",
    )


def label_spd_exposure_observation(
    interval: ActivityInterval,
    free_cmax_um: Any,
    *,
    threshold: float = 10.0,
) -> ExposureLabelResult:
    """Assign the canonical primary exposure label from the feasible margin interval."""

    threshold_decimal, threshold_error = _parse_finite_decimal(threshold)
    if (
        threshold_error is not None
        or threshold_decimal is None
        or threshold_decimal <= 0
    ):
        raise ValueError("exposure threshold must be positive")
    if not interval.is_valid:
        return ExposureLabelResult(
            None,
            None,
            None,
            None,
            None,
            "unknown_invalid_activity",
            interval.activity_parse_reason,
        )
    free_cmax, free_cmax_error = _parse_finite_decimal(free_cmax_um)
    if free_cmax_error is not None:
        return ExposureLabelResult(
            None,
            None,
            None,
            None,
            None,
            "unknown_invalid_free_cmax",
            f"{free_cmax_error}_free_cmax",
        )
    if free_cmax is None:
        raise AssertionError("validated free Cmax unexpectedly missing")
    if free_cmax <= 0:
        return ExposureLabelResult(
            None,
            None,
            None,
            None,
            None,
            "unknown_invalid_free_cmax",
            "nonpositive_free_cmax",
        )

    lower = _exact_interval_bound(
        interval._activity_lower_bound_uM_decimal,
        interval.activity_lower_bound_uM,
    )
    upper = _exact_interval_bound(
        interval._activity_upper_bound_uM_decimal,
        interval.activity_upper_bound_uM,
    )
    if lower is None or upper is None:
        raise AssertionError("valid activity interval is missing bounds")
    margin_lower_decimal = _reporting_decimal_quotient(lower, free_cmax)
    margin_upper_decimal = _reporting_decimal_quotient(upper, free_cmax)
    margin_lower = float(margin_lower_decimal)
    margin_upper = float(margin_upper_decimal)
    lower_inclusive = interval.lower_inclusive
    upper_inclusive = interval.upper_inclusive

    label: int | None
    status: str
    reason: str
    scaled_threshold = _exact_decimal_product(threshold_decimal, free_cmax)
    if upper <= scaled_threshold:
        label, status, reason = (
            SPD_LABEL_POSITIVE,
            "positive",
            "margin_interval_definitely_at_or_below_threshold",
        )
    elif lower > scaled_threshold or (
        lower == scaled_threshold and lower_inclusive is False
    ):
        label, status, reason = (
            SPD_LABEL_NEGATIVE,
            "negative",
            "margin_interval_definitely_above_threshold",
        )
    else:
        label, status, reason = (
            SPD_LABEL_UNKNOWN,
            "unknown",
            "margin_interval_crosses_or_can_equal_threshold",
        )
    return ExposureLabelResult(
        exposure_margin_lower=margin_lower,
        exposure_margin_upper=margin_upper,
        margin_lower_inclusive=lower_inclusive,
        margin_upper_inclusive=upper_inclusive,
        numeric_label=label,
        exposure_status=status,
        exposure_reason=reason,
    )


def interval_record(interval: ActivityInterval) -> dict[str, Any]:
    """Return the stable field names used by audit and DataFrame adapters."""

    return {
        "activity_value_uM": interval.activity_value_uM,
        "activity_lower_bound_uM": interval.activity_lower_bound_uM,
        "activity_upper_bound_uM": interval.activity_upper_bound_uM,
        "lower_inclusive": interval.lower_inclusive,
        "upper_inclusive": interval.upper_inclusive,
        "normalized_relation": interval.normalized_relation,
        "normalized_unit": interval.normalized_unit,
        "activity_parse_status": interval.activity_parse_status,
        "activity_parse_reason": interval.activity_parse_reason,
    }
