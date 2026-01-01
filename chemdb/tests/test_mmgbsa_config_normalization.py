import logging
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import main


def test_md_enabled_true_sets_mode(caplog: pytest.LogCaptureFixture) -> None:
    main._MMGBSA_DEPRECATED_MD_WARNED = False
    cfg = {"MMGBSA_MD_ENABLED": "true"}
    with caplog.at_level(logging.WARNING):
        result = main._mmgbsa_effective_md_config(cfg, logging.getLogger("mmgbsa.test"))
    assert result["md_enabled"] is True
    assert result["effective_mode"] == "IMPLICIT_MD"
    assert "Deprecated MMGBSA keys detected" not in caplog.text


def test_legacy_md_run_warns_and_enables(caplog: pytest.LogCaptureFixture) -> None:
    main._MMGBSA_DEPRECATED_MD_WARNED = False
    cfg = {"MMGBSA_MD_RUN": "true"}
    with caplog.at_level(logging.WARNING):
        result = main._mmgbsa_effective_md_config(cfg, logging.getLogger("mmgbsa.test"))
    assert result["md_enabled"] is True
    assert result["effective_mode"] == "IMPLICIT_MD"
    assert "Deprecated MMGBSA keys detected" in caplog.text


def test_md_disabled_defaults_oneframe() -> None:
    main._MMGBSA_DEPRECATED_MD_WARNED = False
    cfg = {}
    result = main._mmgbsa_effective_md_config(cfg, logging.getLogger("mmgbsa.test"))
    assert result["md_enabled"] is False
    assert result["effective_mode"] == "ONEFRAME"
