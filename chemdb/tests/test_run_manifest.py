"""Tests for run_manifest._refresh_summary variant/pH dedup logic."""

from __future__ import annotations

from run_manifest import _refresh_summary


def test_refresh_summary_dedupes_base_when_ph_present():
    manifest = {
        "proteins": {
            "TEST|APO|base": {"status": "completed"},
            "TEST|APO|pH7_0": {"status": "completed"},
        },
        "summary": {},
    }

    _refresh_summary(manifest)

    proteins = manifest["proteins"]
    assert "TEST|APO|base" not in proteins
    assert "TEST|APO|pH7_0" in proteins
    assert manifest["summary"]["total_proteins_scheduled"] == 1
    assert manifest["summary"]["total_proteins_completed"] == 1
    assert manifest["summary"]["total_proteins_failed"] == 0
