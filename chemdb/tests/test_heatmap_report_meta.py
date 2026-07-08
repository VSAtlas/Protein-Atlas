import csv
import json
from pathlib import Path

from analysis.reporting.heatmap_report_meta import build_heatmap_report_meta


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def test_build_heatmap_report_meta_includes_coverage_and_validation(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    run_id = "META"
    data_dir = repo_root / "data" / run_id
    cache_dir = repo_root / "pathways" / "cache"
    _write_csv(
        data_dir / "master_rows.csv",
        [
            {"ligand_display": "LigA", "ligand_base": "LigA", "target_id": "T1", "pdb_id": "1ABC"},
            {"ligand_display": "LigB", "ligand_base": "LigB", "target_id": "T2", "pdb_id": "2XYZ"},
            {"ligand_display": "LigC", "ligand_base": "LigC", "target_id": "T3", "pdb_id": "3HIJ"},
        ],
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "target_safety_aggregated.json").write_text(
        json.dumps(
            {
                "generated_at": "2026-03-01T12:00:00+00:00",
                "entries": [{"pdb_id": "1ABC"}],
            }
        ),
        encoding="utf-8",
    )
    (cache_dir / "target_safety_drift_summary.json").write_text(
        json.dumps({"generated_at": "2026-03-02T12:00:00+00:00", "changed_count": 2}),
        encoding="utf-8",
    )
    meta = build_heatmap_report_meta(
        repo_root=repo_root,
        run_id=run_id,
        rows=[
            {"ligand_display": "LigA", "pdb_id": "1ABC"},
            {"ligand_display": "LigB", "pdb_id": "2XYZ"},
        ],
        full_row_labels=["LigA", "LigB"],
        full_col_labels=["1ABC\nAlpha", "2XYZ\nBeta"],
        display_top_k=1,
        target_group_meta={
            "1ABC\nAlpha": {"pdb": "1ABC", "target_name": "Alpha", "uniprot_accessions": ["P11111"]},
            "2XYZ\nBeta": {"pdb": "2XYZ", "target_name": "Beta", "uniprot_accessions": []},
        },
    )

    assert meta["report_summary"]["run_id"] == run_id
    assert meta["report_summary"]["generated_at"].endswith("Z")
    assert meta["coverage_summary"]["default_visible_ligands"] == 1
    assert meta["coverage_summary"]["interactive_ligands"] == 2
    assert meta["coverage_summary"]["full_run_ligands"] == 3
    assert meta["coverage_summary"]["full_run_targets"] == 3
    assert meta["coverage_summary"]["matrix_is_subset"] is True
    assert meta["validation_summary"]["cache"]["status"] in {"fresh", "stale"}
    assert meta["validation_summary"]["missing_target_mappings"]["count"] == 1
    assert meta["validation_summary"]["missing_target_mappings"]["targets"][0]["pdb"] == "2XYZ"
    assert meta["validation_summary"]["safety_bucket_drift"]["changed_count"] == 2
