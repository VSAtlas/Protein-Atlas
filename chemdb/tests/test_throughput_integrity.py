from __future__ import annotations

import csv
import json
from pathlib import Path
import datetime as dt
import os

import yaml  # type: ignore[import-untyped]

from analysis.reporting.throughput_integrity import main as throughput_integrity_main


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _build_fixture(repo_root: Path, *, include_all_post_rows: bool) -> str:
    run_id = "ti_smoke"
    combo = ("ABCD", "HOLO", "pH7_0")
    pdb_id, variant, ph = combo

    manifests_dir = repo_root / "outputs" / "manifests" / run_id
    configs_dir = repo_root / "outputs" / "configs" / run_id / pdb_id / variant / ph / "stage1"
    docked_combo_dir = repo_root / "outputs" / "docked" / run_id / pdb_id / variant / ph
    data_dir = repo_root / "outputs" / "data" / run_id
    prepped_lib_dir = repo_root / "prepped_ligands" / "libA"
    receptor_root = repo_root / "outputs" / "processed_pdbs" / run_id / pdb_id / variant / "receptor"
    ph_root = receptor_root / "ph_ensemble"

    manifests_dir.mkdir(parents=True, exist_ok=True)
    configs_dir.mkdir(parents=True, exist_ok=True)
    docked_combo_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    (prepped_lib_dir / "actives").mkdir(parents=True, exist_ok=True)
    (prepped_lib_dir / "decoys").mkdir(parents=True, exist_ok=True)
    ph_root.mkdir(parents=True, exist_ok=True)

    entries = {
        "actives_final_00001": "actives/actives_final_00001.pdbqt",
        "decoys_final_00001": "decoys/decoys_final_00001.pdbqt",
        "decoys_final_00002": "decoys/decoys_final_00002.pdbqt",
    }
    for rel in entries.values():
        target = prepped_lib_dir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("MODEL\nENDMDL\n", encoding="utf-8")
    (prepped_lib_dir / "_manifest.json").write_text(
        json.dumps({"entries": entries, "filenames": {}}, indent=2), encoding="utf-8"
    )

    for lig in entries.keys():
        (configs_dir / f"{lig}_stage1.txt").write_text("cpu = 1\n", encoding="utf-8")
    (receptor_root / f"{pdb_id}_cleaned.pdb").write_text("ATOM\n", encoding="utf-8")
    (ph_root / f"{pdb_id}_{ph}.withH.pdb").write_text("ATOM\n", encoding="utf-8")
    (ph_root / f"{pdb_id}_{ph}.pdbqt").write_text("ATOM\n", encoding="utf-8")

    _write_csv(
        docked_combo_dir / "docking_score_summary.csv",
        ["run_id", "variant", "Ligand", "stage1", "stage2", "stage3"],
        [
            {
                "run_id": run_id,
                "variant": variant,
                "Ligand": f"{lig}.pdbqt",
                "stage1": "-7.0",
                "stage2": "",
                "stage3": "",
            }
            for lig in entries.keys()
        ],
    )

    post_ligs = list(entries.keys())
    if not include_all_post_rows:
        post_ligs = post_ligs[:-1]
    _write_csv(
        data_dir / "master_rows.csv",
        [
            "run_id",
            "pdb_id",
            "variant",
            "ph_label",
            "is_control",
            "ligand_base",
            "ligand_file",
            "final_score",
        ],
        [
            {
                "run_id": run_id,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_label": ph,
                "is_control": "false",
                "ligand_base": lig,
                "ligand_file": f"{lig}.pdbqt",
                "final_score": "1.0",
            }
            for lig in post_ligs
        ],
    )

    log_dir = repo_root / "outputs" / "docked" / run_id / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "pipeline.log").write_text(
        "2026-02-16 08:00:00,000 - INFO - start\n"
        "2026-02-16 08:10:00,000 - INFO - done\n",
        encoding="utf-8",
    )

    manifest_payload = {
        "run_id": run_id,
        "status": "completed",
        "paths": {
            "run_dir": str(repo_root / "outputs" / "configs" / run_id),
            "docked_dir": str(repo_root / "outputs" / "docked"),
            "prepped_ligands_dir": str(repo_root / "prepped_ligands"),
        },
        "proteins": {
            f"{pdb_id}|{variant}|{ph}": {
                "pdb_id": pdb_id,
                "variant": variant,
                "ph": ph,
                "library": "libA",
                "status": "completed",
            }
        },
    }
    (manifests_dir / "run_manifest.yaml").write_text(
        yaml.safe_dump(manifest_payload, sort_keys=False), encoding="utf-8"
    )
    return run_id


def _write_distributed_chunk_artifacts(
    repo_root: Path,
    *,
    run_id: str,
    chunk_status_by_id: dict[str, str],
) -> None:
    dist_dir = repo_root / "outputs" / "manifests" / run_id / "distributed"
    result_dir = dist_dir / "chunk_results"
    result_dir.mkdir(parents=True, exist_ok=True)
    plan = {
        "chunks": [
            {
                "chunk_id": "chunk_000",
                "pdb_id": "ABCD",
                "variant_label": "HOLO",
                "ph_tag": "pH7_0",
                "library_name": "libA",
                "ligand_bases": ["actives_final_00001"],
            },
            {
                "chunk_id": "chunk_001",
                "pdb_id": "ABCD",
                "variant_label": "HOLO",
                "ph_tag": "pH7_0",
                "library_name": "libA",
                "ligand_bases": ["decoys_final_00001", "decoys_final_00002"],
            },
        ]
    }
    (dist_dir / "combo_chunks_holo.json").write_text(
        json.dumps(plan, indent=2), encoding="utf-8"
    )

    for chunk_id, status in chunk_status_by_id.items():
        payload = {"chunk_id": chunk_id, "status": status, "attempt": 1}
        (result_dir / f"{chunk_id}.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )


def test_throughput_integrity_pass(tmp_path: Path) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 0

    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    out_csv = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.csv"
    assert out_json.exists()
    assert out_csv.exists()

    payload = json.loads(out_json.read_text(encoding="utf-8"))
    assert payload["integrity_pass"] is True
    assert payload["totals"]["expected_ligands"] == 3
    assert payload["totals"]["post_scored_ligands"] == 3


def test_throughput_integrity_honors_relocated_manifest_and_data_roots(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    relocated_manifest_dir = tmp_path / "scratch" / "outputs" / "manifests" / run_id
    relocated_data_dir = tmp_path / "scratch" / "outputs" / "data" / run_id
    relocated_manifest_dir.mkdir(parents=True)
    relocated_data_dir.mkdir(parents=True)
    (relocated_manifest_dir / "run_manifest.yaml").write_text(
        (tmp_path / "outputs" / "manifests" / run_id / "run_manifest.yaml").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )
    (relocated_data_dir / "master_rows.csv").write_text(
        (tmp_path / "outputs" / "data" / run_id / "master_rows.csv").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("MANIFESTS_DIR", str(relocated_manifest_dir.parent))
    monkeypatch.setenv("DATA_DIR", str(relocated_data_dir.parent))

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )

    out_json = relocated_data_dir / "throughput_integrity.json"
    assert rc == 0
    assert out_json.exists()
    assert not (tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json").exists()


def test_throughput_integrity_strict_fails_on_missing_post(tmp_path: Path) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=False)
    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 2

    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    assert payload["integrity_pass"] is False
    combo = payload["combos"][0]
    assert combo["missing_post_unexplained_n"] == 1
    assert "missing_post_unexplained" in combo["failure_reasons"]


def test_throughput_integrity_fails_when_expected_unresolved(tmp_path: Path) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    manifest_path = tmp_path / "outputs" / "manifests" / run_id / "run_manifest.yaml"
    payload = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    key = next(iter(payload["proteins"].keys()))
    payload["proteins"][key]["library"] = "missing_library"
    manifest_path.write_text(
        yaml.safe_dump(payload, sort_keys=False), encoding="utf-8"
    )
    stage1_dir = tmp_path / "outputs" / "configs" / run_id / "ABCD" / "HOLO" / "pH7_0" / "stage1"
    for f in stage1_dir.glob("*_stage1.txt"):
        f.unlink()
    (
        tmp_path
        / "outputs"
        / "docked"
        / run_id
        / "ABCD"
        / "HOLO"
        / "pH7_0"
        / "docking_score_summary.csv"
    ).unlink()
    (tmp_path / "outputs" / "data" / run_id / "master_rows.csv").unlink()

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 2
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    result = json.loads(out_json.read_text(encoding="utf-8"))
    combo = result["combos"][0]
    assert combo["expected_n"] == 0
    assert combo["expected_resolved"] is False
    assert "expected_set_unresolved" in combo["failure_reasons"]


def test_throughput_integrity_detects_partial_scorch_coverage_then_passes_when_full(
    tmp_path: Path,
) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)

    pdb2, variant2, ph2 = ("EFGH", "HOLO", "pH7_2")
    ligands2 = ["actives_final_00011", "decoys_final_00011"]
    cfg_stage_dir = tmp_path / "outputs" / "configs" / run_id / pdb2 / variant2 / ph2 / "stage1"
    docked_combo_dir = tmp_path / "outputs" / "docked" / run_id / pdb2 / variant2 / ph2
    prepped_lib_dir = tmp_path / "prepped_ligands" / "libB"
    receptor_root = tmp_path / "outputs" / "processed_pdbs" / run_id / pdb2 / variant2 / "receptor"
    ph_root = receptor_root / "ph_ensemble"

    cfg_stage_dir.mkdir(parents=True, exist_ok=True)
    docked_combo_dir.mkdir(parents=True, exist_ok=True)
    (prepped_lib_dir / "actives").mkdir(parents=True, exist_ok=True)
    (prepped_lib_dir / "decoys").mkdir(parents=True, exist_ok=True)
    ph_root.mkdir(parents=True, exist_ok=True)

    entries2 = {
        "actives_final_00011": "actives/actives_final_00011.pdbqt",
        "decoys_final_00011": "decoys/decoys_final_00011.pdbqt",
    }
    for rel in entries2.values():
        target = prepped_lib_dir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("MODEL\nENDMDL\n", encoding="utf-8")
    (prepped_lib_dir / "_manifest.json").write_text(
        json.dumps({"entries": entries2, "filenames": {}}, indent=2), encoding="utf-8"
    )
    for lig in ligands2:
        (cfg_stage_dir / f"{lig}_stage1.txt").write_text("cpu = 1\n", encoding="utf-8")
    (receptor_root / f"{pdb2}_cleaned.pdb").write_text("ATOM\n", encoding="utf-8")
    (ph_root / f"{pdb2}_{ph2}.withH.pdb").write_text("ATOM\n", encoding="utf-8")
    (ph_root / f"{pdb2}_{ph2}.pdbqt").write_text("ATOM\n", encoding="utf-8")
    _write_csv(
        docked_combo_dir / "docking_score_summary.csv",
        ["run_id", "variant", "Ligand", "stage1", "stage2", "stage3"],
        [
            {
                "run_id": run_id,
                "variant": variant2,
                "Ligand": f"{lig}.pdbqt",
                "stage1": "-7.0",
                "stage2": "",
                "stage3": "",
            }
            for lig in ligands2
        ],
    )

    manifest_path = tmp_path / "outputs" / "manifests" / run_id / "run_manifest.yaml"
    manifest_payload = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest_payload["proteins"][f"{pdb2}|{variant2}|{ph2}"] = {
        "pdb_id": pdb2,
        "variant": variant2,
        "ph": ph2,
        "library": "libB",
        "status": "completed",
    }
    manifest_path.write_text(
        yaml.safe_dump(manifest_payload, sort_keys=False), encoding="utf-8"
    )

    rc_partial = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc_partial == 2
    partial_payload = json.loads(
        (tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json").read_text(
            encoding="utf-8"
        )
    )
    combo2_partial = next(
        row
        for row in partial_payload["combos"]
        if row["pdb_id"] == pdb2 and row["ph"] == ph2
    )
    assert combo2_partial["missing_post_unexplained_n"] == len(ligands2)

    master_rows_path = tmp_path / "outputs" / "data" / run_id / "master_rows.csv"
    with master_rows_path.open("r", encoding="utf-8", newline="") as handle:
        existing_rows = list(csv.DictReader(handle))
    existing_rows.extend(
        [
            {
                "run_id": run_id,
                "pdb_id": pdb2,
                "variant": variant2,
                "ph_label": ph2,
                "is_control": "false",
                "ligand_base": lig,
                "ligand_file": f"{lig}.pdbqt",
                "final_score": "1.0",
            }
            for lig in ligands2
        ]
    )
    _write_csv(
        master_rows_path,
        [
            "run_id",
            "pdb_id",
            "variant",
            "ph_label",
            "is_control",
            "ligand_base",
            "ligand_file",
            "final_score",
        ],
        existing_rows,
    )

    rc_full = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc_full == 0
    full_payload = json.loads(
        (tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json").read_text(
            encoding="utf-8"
        )
    )
    assert full_payload["integrity_pass"] is True


def test_throughput_integrity_fails_when_prep_artifacts_missing(tmp_path: Path) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    ph_path = (
        tmp_path
        / "outputs"
        / "processed_pdbs"
        / run_id
        / "ABCD"
        / "HOLO"
        / "receptor"
        / "ph_ensemble"
        / "ABCD_pH7_0.pdbqt"
    )
    ph_path.unlink()

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 2
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    result = json.loads(out_json.read_text(encoding="utf-8"))
    combo = result["combos"][0]
    assert combo["prep_artifacts_found_n"] == 2
    assert combo["prep_artifacts_required_n"] == 3
    assert "prep_artifacts_missing" in combo["failure_reasons"]


def test_stale_failure_marker_does_not_explain_missing_docking(tmp_path: Path) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)

    summary_path = (
        tmp_path / "outputs" / "docked" / run_id / "ABCD" / "HOLO" / "pH7_0" / "docking_score_summary.csv"
    )
    rows: list[dict[str, str]] = []
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            if row["Ligand"].startswith("decoys_final_00002"):
                row["stage1"] = ""
                row["stage2"] = ""
                row["stage3"] = ""
            rows.append(row)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    marker_dir = tmp_path / "outputs" / "docked" / run_id / "ABCD" / "HOLO" / "pH7_0" / "stage1"
    marker_dir.mkdir(parents=True, exist_ok=True)
    marker = marker_dir / "decoys_final_00002_stage1.pdbqt.failed.txt"
    marker.write_text("failed\n", encoding="utf-8")
    stale_ts = dt.datetime(2020, 1, 1, 0, 0, 0).timestamp()
    marker.touch()
    # Ensure this marker is outside run window and ignored.
    os.utime(marker, (stale_ts, stale_ts))

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 2
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    result = json.loads(out_json.read_text(encoding="utf-8"))
    combo = result["combos"][0]
    assert combo["known_failed_stage1_n"] == 0
    assert combo["missing_docking_unexplained_n"] == 1
    assert "missing_docking_unexplained" in combo["failure_reasons"]


def test_distributed_chunk_plan_is_used_for_expected_and_coverage(tmp_path: Path) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    _write_distributed_chunk_artifacts(
        tmp_path,
        run_id=run_id,
        chunk_status_by_id={"chunk_000": "completed", "chunk_001": "completed"},
    )

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 0
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    combo = payload["combos"][0]
    assert combo["expected_source"] == "distributed_chunk_plan"
    assert combo["chunk_plan_total_n"] == 2
    assert combo["chunk_plan_completed_n"] == 2
    assert combo["chunk_plan_failed_n"] == 0
    assert combo["chunk_plan_missing_result_n"] == 0
    assert payload["totals"]["chunk_coverage_pass"] is True


def test_distributed_chunk_plan_incomplete_fails_integrity(tmp_path: Path) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    _write_distributed_chunk_artifacts(
        tmp_path,
        run_id=run_id,
        chunk_status_by_id={"chunk_000": "completed"},
    )

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 2
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    combo = payload["combos"][0]
    assert combo["expected_source"] == "distributed_chunk_plan"
    assert combo["chunk_plan_total_n"] == 2
    assert combo["chunk_plan_completed_n"] == 1
    assert combo["chunk_plan_missing_result_n"] == 1
    assert "chunk_plan_incomplete" in combo["failure_reasons"]
    assert "chunk_results_missing" in combo["failure_reasons"]


def test_throughput_integrity_ignores_chunk_plan_when_combo_chunks_disabled_env(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    _write_distributed_chunk_artifacts(
        tmp_path,
        run_id=run_id,
        chunk_status_by_id={"chunk_000": "completed"},
    )
    monkeypatch.setenv("ATLAS_DISABLE_COMBO_CHUNKS", "1")

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 0
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    combo = payload["combos"][0]
    assert combo["expected_source"] != "distributed_chunk_plan"
    assert combo["chunk_plan_total_n"] == 0
    assert combo["chunk_plan_completed_n"] == 0
    assert combo["chunk_plan_missing_result_n"] == 0
    assert "chunk_plan_incomplete" not in combo["failure_reasons"]


def test_throughput_integrity_prefers_master_rows_over_stale_snapshot(
    tmp_path: Path,
) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    coverage_dir = tmp_path / "outputs" / "manifests" / run_id / "coverage"
    coverage_dir.mkdir(parents=True, exist_ok=True)
    # Snapshot intentionally under-counts post coverage (drops active ligand),
    # simulating refresh from a mode-specific post CSV.
    (coverage_dir / "ABCD__HOLO__pH7_0.json").write_text(
        json.dumps(
            {
                "run_id": run_id,
                "pdb_id": "ABCD",
                "variant": "HOLO",
                "ph": "pH7_0",
                "library": "libA",
                "expected_ligands": [
                    "actives_final_00001",
                    "decoys_final_00001",
                    "decoys_final_00002",
                ],
                "docking_scored_ligands": [
                    "actives_final_00001",
                    "decoys_final_00001",
                    "decoys_final_00002",
                ],
                "post_scored_ligands": [
                    "decoys_final_00001",
                    "decoys_final_00002",
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 0
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    combo = payload["combos"][0]
    assert combo["expected_source"] == "library_manifest"
    assert combo["post_scored_n"] == 3
    assert combo["missing_post_unexplained_n"] == 0


def test_throughput_integrity_prefers_docking_csv_over_stale_snapshot(
    tmp_path: Path,
) -> None:
    run_id = _build_fixture(tmp_path, include_all_post_rows=True)
    coverage_dir = tmp_path / "outputs" / "manifests" / run_id / "coverage"
    coverage_dir.mkdir(parents=True, exist_ok=True)
    (coverage_dir / "ABCD__HOLO__pH7_0.json").write_text(
        json.dumps(
            {
                "run_id": run_id,
                "pdb_id": "ABCD",
                "variant": "HOLO",
                "ph": "pH7_0",
                "library": "libA",
                "expected_ligands": [
                    "actives_final_00001",
                    "decoys_final_00001",
                    "decoys_final_00002",
                ],
                "docking_scored_ligands": [
                    "actives_final_00001",
                    "decoys_final_00001",
                ],
                "post_scored_ligands": [
                    "actives_final_00001",
                    "decoys_final_00001",
                    "decoys_final_00002",
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    rc = throughput_integrity_main(
        ["--run-id", run_id, "--repo-root", str(tmp_path), "--strict"]
    )
    assert rc == 0
    out_json = tmp_path / "outputs" / "data" / run_id / "throughput_integrity.json"
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    combo = payload["combos"][0]
    assert combo["docking_scored_n"] == 3
    assert combo["missing_docking_unexplained_n"] == 0
