import csv
import os
from pathlib import Path

import pytest

from chemdb.tools import fill_fda_mapping_names as tool


def _fail_on_low_quality_enabled() -> bool:
    text = str(os.getenv("FAIL_ON_LOW_QUALITY", "1")).strip().lower()
    return text in {"1", "true", "yes", "y", "on"}


def test_fda_mapping_has_no_low_quality_display_names():
    if not _fail_on_low_quality_enabled():
        pytest.skip("FAIL_ON_LOW_QUALITY is disabled.")

    csv_path = os.getenv("FDA_MAPPING_CSV", "").strip()
    if not csv_path:
        pytest.fail("FDA_MAPPING_CSV is required for this gate.")

    path = Path(csv_path)
    if not path.exists():
        pytest.fail(f"FDA_MAPPING_CSV does not exist: {path}")

    bad_rows = 0
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            display_name = (row.get("display_name") or "").strip()
            if not display_name:
                bad_rows += 1
                continue
            if tool.is_cas_only_name(display_name):
                bad_rows += 1
                continue
            if tool.is_formula_like_name(display_name):
                bad_rows += 1
                continue
            if tool.is_unii_like_name(display_name):
                bad_rows += 1
                continue
            if tool.is_spaced_iupacish_name(display_name):
                bad_rows += 1
                continue

    assert bad_rows == 0, f"Found {bad_rows} low-quality display names."
