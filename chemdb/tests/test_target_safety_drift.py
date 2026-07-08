from analysis.reporting.target_safety_drift import build_target_safety_drift_summary


def test_build_target_safety_drift_summary_detects_bucket_changes() -> None:
    previous = {
        "generated_at": "2026-03-01T00:00:00+00:00",
        "entries": [
            {
                "pdb_id": "1ABC",
                "target_name": "Alpha",
                "primary_display_safety": "Endocrine",
                "safety_confidence": "high",
                "safety_sources": ["Open Targets safety"],
            }
        ],
    }
    current = {
        "generated_at": "2026-03-02T00:00:00+00:00",
        "entries": [
            {
                "pdb_id": "1ABC",
                "target_name": "Alpha",
                "primary_display_safety": "Hepatotoxicity",
                "safety_confidence": "medium",
                "safety_sources": ["Open Targets drug adverse events"],
            }
        ],
    }

    summary = build_target_safety_drift_summary(previous, current)

    assert summary["changed_count"] == 1
    assert summary["changes"][0]["pdb_id"] == "1ABC"
    assert summary["changes"][0]["old_primary_display_safety"] == "Endocrine"
    assert summary["changes"][0]["new_primary_display_safety"] == "Hepatotoxicity"
