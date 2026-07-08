import csv
from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa
import pyarrow.parquet as pq

from analysis.cli import convert_interactions_to_parquet as converter
from analysis.reporting import heatmap_html
from analysis.reporting import run_report_core as run_report


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _part_counts(out_dir: Path, run_id: str, targets: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for target in targets:
        part_dir = out_dir / f"run_id={run_id}" / f"target_id={target}"
        counts[target] = len(list(part_dir.glob("part-*.parquet")))
    return counts


def test_convert_interactions_to_parquet_writes_partitions(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    run_id = "RUN123"
    master_csv = repo_root / "data" / run_id / "master_rows.csv"
    rows = [
        {
            "pdb_id": "PDB1",
            "variant": "APO",
            "ph_label": "7",
            "ligand_base": "LIG1",
            "ligand_display": "Lig One",
            "z_selected": "2.5",
            "is_decoy": "0",
            "is_control": "0",
            "pose_valid_any": "true",
            "rank": "1",
            "pct_rank": "0.10",
        },
        {
            "pdb_id": "PDB2",
            "variant": "HOLO",
            "ph_label": "8",
            "ligand_base": "LIG2",
            "ligand_display": "",
            "z_selected": "1.25",
            "is_decoy": "1",
            "is_control": "0",
            "pose_valid_any": "false",
            "rank": "2",
            "pct_rank": "0.20",
        },
        {
            "pdb_id": "PDB2",
            "variant": "HOLO",
            "ph_label": "8",
            "ligand_base": "LIG3",
            "ligand_display": "Lig Three",
            "z_selected": "not-a-number",
            "is_decoy": "0",
            "is_control": "0",
            "pose_valid_any": "true",
            "rank": "3",
            "pct_rank": "0.30",
        },
    ]
    _write_csv(master_csv, rows)

    exit_code = converter.main(["--run-id", run_id, "--repo-root", str(repo_root)])
    assert exit_code == 0

    out_dir = repo_root / "data" / run_id / "dataset" / "interactions"
    target_a = out_dir / "run_id=RUN123" / "target_id=PDB1|APO|7"
    target_b = out_dir / "run_id=RUN123" / "target_id=PDB2|HOLO|8"
    assert list(target_a.glob("part-*.parquet"))
    assert list(target_b.glob("part-*.parquet"))

    part_path = next(target_a.glob("part-*.parquet"))
    table = pq.ParquetFile(part_path).read()
    schema = table.schema
    for col in (
        "run_id",
        "target_id",
        "pdb_id",
        "variant",
        "ph_label",
        "ligand_base",
        "ligand_display",
        "z_selected",
    ):
        assert col in schema.names
    assert schema.field("z_selected").type == pa.float64()
    assert schema.field("run_id").type == pa.string()
    assert schema.field("is_decoy").type == pa.bool_()
    assert schema.field("rank").type == pa.int32()
    assert schema.field("pct_rank").type == pa.float32()


def test_heatmap_reads_parquet_same_as_csv(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    run_id = "RUN_HEAT"
    heatmap_csv = repo_root / "data" / run_id / "heatmap_input.csv"
    rows = [
        {
            "target_id": "T1|APO|7",
            "ligand_display": "LIG1",
            "ligand_base": "LIG1",
            "z_selected": "1.0",
            "rank": "2",
            "pct_rank": "0.2",
            "pose_valid_any": "1",
            "library": "FDA",
            "is_decoy": "0",
        },
        {
            "target_id": "T1|APO|7",
            "ligand_display": "LIG1",
            "ligand_base": "LIG1",
            "z_selected": "2.0",
            "rank": "1",
            "pct_rank": "0.1",
            "pose_valid_any": "1",
            "library": "FDA",
            "is_decoy": "0",
        },
        {
            "target_id": "T2|HOLO|8",
            "ligand_display": "LIG2",
            "ligand_base": "LIG2",
            "z_selected": "0.5",
            "rank": "3",
            "pct_rank": "0.3",
            "pose_valid_any": "1",
            "library": "FDA",
            "is_decoy": "0",
        },
        {
            "target_id": "T2|HOLO|8",
            "ligand_display": "DECOY1",
            "ligand_base": "DECOY1",
            "z_selected": "3.5",
            "rank": "1",
            "pct_rank": "0.1",
            "pose_valid_any": "1",
            "library": "DECOY",
            "is_decoy": "1",
        },
    ]
    _write_csv(heatmap_csv, rows)

    exit_code = converter.main(
        [
            "--run-id",
            run_id,
            "--repo-root",
            str(repo_root),
            "--input",
            str(heatmap_csv),
        ]
    )
    assert exit_code == 0

    parquet_dir = repo_root / "data" / run_id / "dataset" / "interactions"
    csv_rows = heatmap_html._load_heatmap_rows(heatmap_csv, include_decoys=False)
    parquet_rows = heatmap_html._load_heatmap_rows(parquet_dir, include_decoys=False)

    csv_agg, _ = heatmap_html._aggregate_rows(csv_rows)
    parquet_agg, _ = heatmap_html._aggregate_rows(parquet_rows)

    csv_scores = {key: value["z_selected"] for key, value in csv_agg.items()}
    parquet_scores = {key: value["z_selected"] for key, value in parquet_agg.items()}
    assert set(csv_scores.keys()) == set(parquet_scores.keys())
    assert csv_scores == parquet_scores


def test_partial_rebuild_targets(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    run_id = "RUN_PARTIAL"
    master_csv = repo_root / "data" / run_id / "master_rows.csv"
    rows = [
        {
            "target_id": "PDB1|APO|7",
            "ligand_display": "L1",
            "ligand_base": "L1",
            "z_selected": "1.0",
            "is_decoy": "0",
            "is_control": "0",
        },
        {
            "target_id": "PDB2|HOLO|8",
            "ligand_display": "L2",
            "ligand_base": "L2",
            "z_selected": "2.0",
            "is_decoy": "0",
            "is_control": "0",
        },
    ]
    _write_csv(master_csv, rows)

    first_exit = converter.main(["--run-id", run_id, "--repo-root", str(repo_root)])
    assert first_exit == 0

    out_dir = repo_root / "data" / run_id / "dataset" / "interactions"
    targets = ["PDB1|APO|7", "PDB2|HOLO|8"]
    before_counts = _part_counts(out_dir, run_id, targets)

    second_exit = converter.main(
        [
            "--run-id",
            run_id,
            "--repo-root",
            str(repo_root),
            "--targets",
            "PDB1|APO|7",
        ]
    )
    assert second_exit == 0

    after_counts = _part_counts(out_dir, run_id, targets)
    assert after_counts["PDB1|APO|7"] > before_counts["PDB1|APO|7"]
    assert after_counts["PDB2|HOLO|8"] == before_counts["PDB2|HOLO|8"]


def test_run_report_prefers_parquet_dataset_for_heatmap(
    tmp_path: Path, monkeypatch
) -> None:
    repo_root = tmp_path / "repo"
    run_id = "RUN_REPORT_PARQUET"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True, exist_ok=True)

    parquet_dir = data_dir / "dataset" / "interactions"
    parquet_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "heatmap_input.csv").write_text("placeholder\n", encoding="utf-8")
    (data_dir / "master_rows.csv").write_text("placeholder\n", encoding="utf-8")

    seen: dict[str, Path] = {}

    def _fake_render(_repo_root, _run_id, source_path, *_args, **_kwargs):
        seen["path"] = Path(source_path)
        return "<div>ok</div>"

    monkeypatch.setattr(run_report, "render_interactive_heatmap_html", _fake_render)

    out_path = data_dir / "report.html"
    run_report._write_html_report(
        {
            "run_id": run_id,
            "generated_at": "2026-01-01T00:00:00Z",
            "targets": {},
            "multi_target_hits": [],
        },
        out_path,
        [],
        repo_root,
    )

    assert seen["path"] == parquet_dir
