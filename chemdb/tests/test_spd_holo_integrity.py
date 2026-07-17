from __future__ import annotations

import json
import os
from pathlib import Path

from analysis.reporting.spd_holo_integrity import _CohortRow, _audit_target


def _write_layout(layout: Path, *, retention_name: str) -> tuple[Path, Path]:
    ligands = layout / "ligands_raw"
    receptor_dir = layout / "receptor"
    work = layout / "work"
    ligands.mkdir(parents=True)
    receptor_dir.mkdir()
    work.mkdir()

    ligand = ligands / "ATP.pdb"
    ligand.write_text(
        "HETATM    1  C1  ATP A   1      11.104  13.207   8.678  1.00 20.00           C\n",
        encoding="utf-8",
    )
    receptor = receptor_dir / "prepared.pdbqt"
    receptor.write_text("REMARK prepared receptor\n", encoding="utf-8")
    meeko = receptor_dir / "prepared.meeko_input_stage.json"
    meeko.write_text(
        json.dumps({"ok": True, "output_pdbqt": str(receptor)}),
        encoding="utf-8",
    )

    source_pdb = work / "retention_source.pdb"
    target_pdb = work / "retention_target.pdb"
    source_pdb.write_text("ATOM\n", encoding="utf-8")
    target_pdb.write_text("ATOM\n", encoding="utf-8")
    retention = work / retention_name
    retention.write_text(
        json.dumps(
            {
                "source_pdb": str(source_pdb),
                "target_pdb": str(target_pdb),
                "after": {"cofactors": 1, "metals": 0},
                "missing_after_count": 0,
                "candidates": [{"token": "ATP"}],
            }
        ),
        encoding="utf-8",
    )
    return receptor, retention


def _manifest_entry(receptor: Path) -> dict[str, object]:
    return {
        "pdb_id": "1ABC",
        "stages": {
            "prep": {
                "status": "completed",
                "details": {"receptor_pdbqt": str(receptor)},
            }
        },
    }


def _cohort() -> _CohortRow:
    return _CohortRow(
        gene="GENE1",
        pdb_id="1ABC",
        expected_input_class="strict_holo",
        expected_components=("ATP",),
    )


def test_manifest_lineage_ignores_unrelated_newest_retention_sidecar(
    tmp_path: Path,
) -> None:
    processed = tmp_path / "processed"
    target = processed / "1ABC"
    selected = target / "variant_a"
    unrelated = target / "variant_b"
    receptor, selected_retention = _write_layout(
        selected, retention_name="selected.retained_hets_audit.json"
    )
    _other_receptor, unrelated_retention = _write_layout(
        unrelated, retention_name="unrelated.retained_hets_audit.json"
    )
    newer = selected_retention.stat().st_mtime_ns + 10_000_000
    os.utime(unrelated_retention, ns=(newer, newer))

    row = _audit_target(
        _cohort(),
        repo_root=tmp_path,
        run_id="run_a",
        processed_root=processed,
        manifest_entries=(("1ABC|A", _manifest_entry(receptor)),),
    )

    assert row["cohort_requirement_certified"] is True
    assert row["artifact_layouts"] == [str(selected)]
    assert row["retention_sidecar_path"] == str(selected_retention)
    assert str(unrelated_retention) not in row["evidence_sha256s"]
    assert row["reason_codes"] == []


def test_conflicting_manifest_layouts_are_rejected_without_mixing_artifacts(
    tmp_path: Path,
) -> None:
    processed = tmp_path / "processed"
    target = processed / "1ABC"
    first = target / "variant_a"
    second = target / "variant_b"
    first_receptor, _first_retention = _write_layout(
        first, retention_name="first.retained_hets_audit.json"
    )
    second_receptor, _second_retention = _write_layout(
        second, retention_name="second.retained_hets_audit.json"
    )

    row = _audit_target(
        _cohort(),
        repo_root=tmp_path,
        run_id="run_b",
        processed_root=processed,
        manifest_entries=(
            ("1ABC|A", _manifest_entry(first_receptor)),
            ("1ABC|B", _manifest_entry(second_receptor)),
        ),
    )

    assert row["cohort_requirement_certified"] is False
    assert "artifact_layout_ambiguous" in row["reason_codes"]
    assert row["ligand_component_paths"] == []
    assert row["receptor_pdbqt_paths"] == []
    assert row["meeko_sidecar_paths"] == []
    assert row["retention_status"] == "not_reported"
