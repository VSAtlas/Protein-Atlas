import csv
import os
import re
from pathlib import Path

import pytest

FDA_PLACEHOLDER_RE = re.compile(r"^fda_\d+$", re.IGNORECASE)


def _fail_on_unk_enabled() -> bool:
    text = str(os.getenv("FAIL_ON_UNK", "1")).strip().lower()
    return text in {"1", "true", "yes", "y", "on"}


def test_fda_mapping_has_no_placeholders():
    if not _fail_on_unk_enabled():
        pytest.skip("FAIL_ON_UNK is disabled.")

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
            if display_name.upper().startswith("UNK_"):
                bad_rows += 1
                continue
            if FDA_PLACEHOLDER_RE.match(display_name):
                bad_rows += 1
                continue

    assert bad_rows == 0, f"Found {bad_rows} unresolved placeholder display names."
