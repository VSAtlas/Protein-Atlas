"""Shared ``TEST_MODE_ENABLE`` parsing for master export and run reporting."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import List

from analysis.reporting.report_config import read_config_key, report_config_candidates

TEST_MODE_OFF_VALUES = frozenset({"", "0", "false", "no", "off", "none", "null"})
TEST_MODE_ON_VALUES = frozenset({"true", "yes", "on", "1"})
TEST_MODE_BOTH_VALUES = frozenset(
    {
        "both",
        "fda_dud",
        "dud_fda",
        "fda+dud",
        "dud+fda",
        "fda-dud",
        "dud-fda",
    }
)


def read_test_mode_enable_from_file(path: Path) -> str | None:
    return read_config_key(path, "TEST_MODE_ENABLE")


def parse_test_mode_value(raw: object) -> List[str]:
    if isinstance(raw, bool):
        return ["dud", "fda"] if raw else ["fda"]
    s = str(raw).strip()
    if not s:
        return ["fda"]
    lowered = s.lower()
    if lowered in TEST_MODE_OFF_VALUES:
        return ["fda"]
    if lowered in TEST_MODE_ON_VALUES:
        return ["dud", "fda"]
    if lowered == "default":
        return ["fda"]
    if lowered in TEST_MODE_BOTH_VALUES:
        return ["dud", "fda"]

    tokens = [tok for tok in re.split(r"[+,\s]+", lowered) if tok]
    normalized: List[str] = []
    for tok in tokens:
        if tok in {"and", "off", "none", "null"}:
            continue
        if tok == "default":
            tok = "fda"
        normalized.append(tok)
    if not normalized:
        return ["fda"]
    deduped: List[str] = []
    seen: set[str] = set()
    for tok in normalized:
        if tok not in seen:
            seen.add(tok)
            deduped.append(tok)
    return deduped


def resolve_test_mode_tokens(repo_root: Path, run_id: str) -> List[str]:
    env_raw = os.environ.get("TEST_MODE_ENABLE")
    if env_raw is not None:
        return parse_test_mode_value(env_raw)
    for path in report_config_candidates(repo_root, run_id):
        raw = read_test_mode_enable_from_file(path)
        if raw is not None:
            return parse_test_mode_value(raw)
    return ["fda"]
