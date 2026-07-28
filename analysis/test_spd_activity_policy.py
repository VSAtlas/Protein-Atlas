from __future__ import annotations

import math

import pandas as pd
import pytest

from analysis.spd_activity_policy import (
    SPD_ACTIVITY_PARSER_VERSION,
    label_spd_binding_observation,
    label_spd_exposure_observation,
    parse_spd_activity_interval,
)


@pytest.mark.parametrize(
    ("relation", "lower", "upper", "lower_inclusive", "upper_inclusive"),
    [
        ("=", 2.0, 2.0, True, True),
        (">", 2.0, math.inf, False, False),
        (">=", 2.0, math.inf, True, False),
        ("<", 0.0, 2.0, True, False),
        ("<=", 0.0, 2.0, True, True),
    ],
)
def test_parse_activity_relations(
    relation: str,
    lower: float,
    upper: float,
    lower_inclusive: bool,
    upper_inclusive: bool,
) -> None:
    parsed = parse_spd_activity_interval(2.0, relation, "uM")

    assert parsed.activity_parse_status == "valid"
    assert parsed.normalized_relation == relation
    assert parsed.activity_value_uM == 2.0
    assert parsed.activity_lower_bound_uM == lower
    assert parsed.activity_upper_bound_uM == upper
    assert parsed.lower_inclusive is lower_inclusive
    assert parsed.upper_inclusive is upper_inclusive


@pytest.mark.parametrize(
    ("value", "unit", "expected_um"),
    [
        (1.0, "M", 1_000_000.0),
        (1.0, "mM", 1_000.0),
        (1.0, "uM", 1.0),
        (1.0, "µM", 1.0),
        (1.0, "μM", 1.0),
        (1.0, "nM", 0.001),
        (1.0, "pM", 0.000001),
        (" 2.5 ", "micromolar", 2.5),
    ],
)
def test_parse_activity_units(value: object, unit: str, expected_um: float) -> None:
    parsed = parse_spd_activity_interval(value, "=", unit)

    assert parsed.activity_parse_status == "valid"
    assert parsed.activity_value_uM == expected_um


@pytest.mark.parametrize(
    ("value", "relation", "unit", "reason"),
    [
        (1.0, "", "uM", "blank_relation"),
        (1.0, "   ", "uM", "blank_relation"),
        (1.0, None, "uM", "missing_relation"),
        (1.0, pd.NA, "uM", "missing_relation"),
        (1.0, "≤", "uM", "unsupported_relation"),
        (1.0, "eq", "uM", "unsupported_relation"),
        (1.0, "~", "uM", "unsupported_relation"),
        ("bad", "=", "uM", "malformed_activity_value"),
        ("", "=", "uM", "missing_activity_value"),
        (None, "=", "uM", "missing_activity_value"),
        (math.inf, "=", "uM", "nonfinite_activity_value"),
        (1.0, "=", "mg/L", "unsupported_activity_unit"),
        (1.0, "=", None, "missing_activity_unit"),
        (-1.0, "=", "uM", "negative_activity_value"),
        (0.0, "<", "uM", "empty_activity_interval"),
    ],
)
def test_parse_activity_fails_closed(
    value: object,
    relation: object,
    unit: object,
    reason: str,
) -> None:
    parsed = parse_spd_activity_interval(value, relation, unit)

    assert parsed.activity_parse_status == "invalid"
    assert parsed.activity_parse_reason == reason
    assert parsed.activity_lower_bound_uM is None
    assert parsed.activity_upper_bound_uM is None


@pytest.mark.parametrize(
    ("value", "relation", "expected_label", "expected_status"),
    [
        (0.5, "=", 1, "positive"),
        (1.0, "=", 1, "positive"),
        (1.0, "<", 1, "positive"),
        (1.0, "<=", 1, "positive"),
        (0.5, ">", None, "unknown"),
        (5.0, "=", None, "unknown"),
        (20.0, "<", None, "unknown"),
        (1.0, ">=", None, "unknown"),
        (10.0, "=", 0, "negative"),
        (10.0, ">", 0, "negative"),
        (10.0, ">=", 0, "negative"),
    ],
)
def test_binding_contract(
    value: float,
    relation: str,
    expected_label: int | None,
    expected_status: str,
) -> None:
    interval = parse_spd_activity_interval(value, relation, "uM")
    result = label_spd_binding_observation(interval)

    assert result.numeric_label == expected_label
    assert result.binding_status == expected_status
    assert result.binding_reason


def test_binding_invalid_activity_is_unknown_with_reason() -> None:
    interval = parse_spd_activity_interval(2.0, "", "uM")
    result = label_spd_binding_observation(interval)

    assert result.numeric_label is None
    assert result.binding_status == "unknown_invalid_activity"
    assert result.binding_reason == "blank_relation"


@pytest.mark.parametrize(
    ("value", "relation", "free_cmax", "expected_label", "expected_status"),
    [
        (10.0, "=", 1.0, 1, "positive"),
        (10.0, "<=", 1.0, 1, "positive"),
        (10.0, "<", 1.0, 1, "positive"),
        (10.0, ">", 1.0, 0, "negative"),
        (20.0, "=", 1.0, 0, "negative"),
        (10.0, ">=", 1.0, None, "unknown"),
        (20.0, "<", 1.0, None, "unknown"),
        (5.0, ">", 1.0, None, "unknown"),
        (5.0, "=", None, None, "unknown_invalid_free_cmax"),
        (5.0, "=", 0.0, None, "unknown_invalid_free_cmax"),
        (5.0, "=", -1.0, None, "unknown_invalid_free_cmax"),
    ],
)
def test_exposure_contract(
    value: float,
    relation: str,
    free_cmax: object,
    expected_label: int | None,
    expected_status: str,
) -> None:
    interval = parse_spd_activity_interval(value, relation, "uM")
    result = label_spd_exposure_observation(interval, free_cmax)

    assert result.numeric_label == expected_label
    assert result.exposure_status == expected_status
    assert result.exposure_reason


def test_exposure_preserves_margin_interval_and_inclusivity() -> None:
    interval = parse_spd_activity_interval(10.0, ">=", "uM")
    result = label_spd_exposure_observation(interval, 2.0)

    assert result.exposure_margin_lower == 5.0
    assert result.exposure_margin_upper == math.inf
    assert result.margin_lower_inclusive is True
    assert result.margin_upper_inclusive is False


def test_exposure_invalid_activity_is_unknown_with_reason() -> None:
    interval = parse_spd_activity_interval(2.0, "unexpected", "uM")
    result = label_spd_exposure_observation(interval, 1.0)

    assert result.numeric_label is None
    assert result.exposure_status == "unknown_invalid_activity"
    assert result.exposure_reason == "unsupported_relation"


@pytest.mark.parametrize(
    (
        "activity_value",
        "free_cmax",
        "relation",
        "expected_label",
        "expected_status",
    ),
    [
        ("0.7", "0.07", ">", 0, "negative"),
        ("0.7", "0.07", ">=", None, "unknown"),
        ("0.7", "0.07", "=", 1, "positive"),
        ("0.3", "0.03", ">", 0, "negative"),
        ("0.3", "0.03", ">=", None, "unknown"),
        ("0.3", "0.03", "=", 1, "positive"),
    ],
)
def test_exposure_uses_exact_decimal_threshold_adjudication(
    activity_value: str,
    free_cmax: str,
    relation: str,
    expected_label: int | None,
    expected_status: str,
) -> None:
    interval = parse_spd_activity_interval(activity_value, relation, "uM")

    result = label_spd_exposure_observation(interval, free_cmax)

    assert result.numeric_label == expected_label
    assert result.exposure_status == expected_status
    assert result.exposure_margin_lower == 10.0


@pytest.mark.parametrize(
    ("relation", "expected_label", "expected_status"),
    [
        (">", 0, "negative"),
        (">=", None, "unknown"),
        ("=", 1, "positive"),
        ("<", 1, "positive"),
        ("<=", 1, "positive"),
    ],
)
def test_exposure_exact_threshold_survives_nm_to_um_conversion(
    relation: str,
    expected_label: int | None,
    expected_status: str,
) -> None:
    interval = parse_spd_activity_interval("700", relation, "nM")

    result = label_spd_exposure_observation(interval, "0.07")

    assert interval.activity_value_uM == 0.7
    assert result.numeric_label == expected_label
    assert result.exposure_status == expected_status


@pytest.mark.parametrize(
    ("activity_value", "unit", "relation", "expected_label"),
    [
        ("1000", "nM", "=", 1),
        ("1000", "nM", ">", None),
        ("10000", "nM", "=", 0),
        ("10000", "nM", ">=", 0),
        ("1000000", "pM", "=", 1),
        ("10000000", "pM", "=", 0),
    ],
)
def test_binding_exact_decimal_unit_conversion_boundaries(
    activity_value: str,
    unit: str,
    relation: str,
    expected_label: int | None,
) -> None:
    interval = parse_spd_activity_interval(activity_value, relation, unit)

    result = label_spd_binding_observation(interval)

    assert result.numeric_label == expected_label


@pytest.mark.parametrize("threshold", [math.nan, math.inf, -math.inf])
def test_nonfinite_thresholds_are_rejected(threshold: float) -> None:
    interval = parse_spd_activity_interval("1", "=", "uM")

    with pytest.raises(ValueError):
        label_spd_exposure_observation(interval, "0.1", threshold=threshold)


def test_parser_implementation_version_tracks_decimal_adjudication() -> None:
    assert SPD_ACTIVITY_PARSER_VERSION == "spd_activity_interval_v3"


def test_threshold_validation() -> None:
    interval = parse_spd_activity_interval(1.0, "=", "uM")

    with pytest.raises(ValueError):
        label_spd_binding_observation(interval, active_um=10.0, inactive_um=10.0)
    with pytest.raises(ValueError):
        label_spd_exposure_observation(interval, 1.0, threshold=0.0)
