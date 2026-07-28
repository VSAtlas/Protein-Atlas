from __future__ import annotations

from analysis.cli.audit_spd_decimal_boundary import (
    BOUNDARY_COLUMNS,
    COMPARISON_COLUMNS,
    EXPECTED_COUNTS,
    _boundary_cases,
    _is_pass_status,
)


def test_upstream_success_statuses_are_documented_variants() -> None:
    assert _is_pass_status("pass")
    assert _is_pass_status("passed")
    assert not _is_pass_status("failed")


def test_decimal_boundary_cases_are_independently_expected_and_pass() -> None:
    cases = _boundary_cases()

    assert cases.columns.tolist() == BOUNDARY_COLUMNS
    assert bool(cases["passed"].all())
    assert cases.set_index("case_id").loc[
        "exposure_07_open", "actual_exposure_label"
    ] == "negative"
    assert cases.set_index("case_id").loc[
        "exposure_07_closed", "actual_exposure_label"
    ] == "unknown"
    assert cases.set_index("case_id").loc[
        "exposure_07_exact", "actual_exposure_label"
    ] == "positive"


def test_decimal_audit_schemas_and_frozen_counts_are_explicit() -> None:
    assert COMPARISON_COLUMNS[0] == "comparison_level"
    assert "binding_prompt2c_changed" in COMPARISON_COLUMNS
    assert "exposure_prompt2c_changed" in COMPARISON_COLUMNS
    assert EXPECTED_COUNTS["raw_exposure"] == {
        "positive": 1_680,
        "negative": 50_095,
        "unknown": 69_322,
    }
    assert EXPECTED_COUNTS["pair_exposure"] == {
        "positive": 1_397,
        "negative": 37_684,
        "unknown": 56_431,
    }
