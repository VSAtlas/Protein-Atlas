from __future__ import annotations

import csv
import json
from pathlib import Path

from cli.qol import manifest_repair as manifest_repair_mod
from cli.qol.manifest_repair import repair_run_manifest
from cli.resume_manifest_repair import auto_repair_resume_manifest
from cli.qol.slurm import _cmd_slurm
from cli.run_manifest_support import distributed_enabled, distributed_mode

try:
    import yaml  # type: ignore[import-untyped]
except Exception:  # pragma: no cover
    yaml = None  # type: ignore[assignment]


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_yaml(path: Path, payload: object) -> None:
    assert yaml is not None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def _read_yaml(path: Path) -> dict:
    assert yaml is not None
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _read_csv_values(path: Path, column: str) -> list[str]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [str(row.get(column) or "") for row in csv.DictReader(handle)]


def _write_complete_fda_scorch(root: Path, run_id: str, pdb_id: str, *, nested: bool = False) -> None:
    base = root / "outputs" / "post_docked" / run_id / pdb_id
    if nested:
        base = base / "LEGACY" / "base"
    base.mkdir(parents=True, exist_ok=True)
    (base / "scorch_scores_all.csv").write_text(
        "pdb_id,variant,ph,source,run_mode,Ligand_ID,SCORCH_score,SCORCH_certainty\n"
        f"{pdb_id},,,vina,fda,a,1,1\n",
        encoding="utf-8",
    )
    (base / "consensus_reranked_scorch.csv").write_text(
        "ligand,run_mode,rescored_flag,SCORCH_score_used\n"
        "a,fda,1,1\n",
        encoding="utf-8",
    )
    _write_json(
        base / "scorch" / "scorch_coverage_summary.json",
        {
            "all_complete": True,
            "streams": {
                "fda": {
                    "complete": True,
                    "expected_rows": 1,
                    "raw_rows": 1,
                    "final_rescored_rows": 1,
                    "invalid_rows": 0,
                },
                "dud": {
                    "complete": True,
                    "expected_rows": 0,
                    "raw_rows": 0,
                    "final_rescored_rows": 0,
                    "invalid_rows": 0,
                },
            },
        },
    )


def _build_repair_fixture(root: Path) -> Path:
    run_id = "repair_rid"
    manifest_run = root / "outputs" / "manifests" / run_id
    dist = manifest_run / "distributed"
    _write_yaml(
        manifest_run / "run_manifest.yaml",
        {
            "run_id": run_id,
            "status": "completed",
            "summary": {
                "total_proteins_scheduled": 99,
                "total_protein_list": ["STALE"],
                "total_proteins_completed": 99,
                "total_proteins_failed": 0,
            },
            "proteins": {
                "PART|LEGACY|base": {
                    "pdb_id": "PART",
                    "variant": "LEGACY",
                    "ph": "base",
                    "status": "completed",
                }
            },
        },
    )
    _write_json(
        dist / "combo_chunks_0.json",
        {
            "chunks": [
                {
                    "chunk_id": "full_0",
                    "pdb_id": "FULL",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "library_name": "lib",
                    "ligand_bases": ["a"],
                },
                {
                    "chunk_id": "full_1",
                    "pdb_id": "FULL",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "library_name": "lib",
                    "ligand_bases": ["b"],
                },
                {
                    "chunk_id": "part_0",
                    "pdb_id": "PART",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "library_name": "lib",
                    "ligand_bases": ["a"],
                },
                {
                    "chunk_id": "part_1",
                    "pdb_id": "PART",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "library_name": "lib",
                    "ligand_bases": ["b"],
                },
                {
                    "chunk_id": "fail_0",
                    "pdb_id": "FAIL",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "library_name": "lib",
                    "ligand_bases": ["a"],
                },
            ]
        },
    )
    for chunk_id, status in {
        "full_0": "completed",
        "full_1": "completed",
        "part_0": "completed",
        "fail_0": "failed",
    }.items():
        _write_json(
            dist / "chunk_results" / f"{chunk_id}.json",
            {"chunk_id": chunk_id, "status": status},
        )

    for pdb_id in ("FULL", "PART"):
        path = root / "outputs" / "docked" / run_id / pdb_id / "consensus_docking_scores.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("ligand,consensus_score\na,1\n", encoding="utf-8")
    reranked = (
        root
        / "outputs"
        / "post_docked"
        / run_id
        / "FULL"
        / "LEGACY"
        / "base"
        / "consensus_reranked_scorch.csv"
    )
    reranked.parent.mkdir(parents=True, exist_ok=True)
    _write_complete_fda_scorch(root, run_id, "FULL", nested=True)
    return manifest_run / "run_manifest.yaml"


def test_repair_manifest_rebuilds_from_chunk_coverage(tmp_path: Path) -> None:
    assert yaml is not None
    manifest_path = _build_repair_fixture(tmp_path)

    summary = repair_run_manifest("repair_rid", root=tmp_path)

    repaired = _read_yaml(manifest_path)
    proteins = repaired["proteins"]
    assert proteins["FULL|LEGACY|base"]["status"] == "completed"
    assert proteins["PART|LEGACY|base"]["status"] == "running"
    assert proteins["FAIL|LEGACY|base"]["status"] == "failed"
    assert repaired["summary"]["total_proteins_scheduled"] == 3
    assert repaired["summary"]["total_proteins_completed"] == 1
    assert repaired["summary"]["total_proteins_failed"] == 1
    assert summary.planned_chunks == 5
    assert summary.completed_chunks == 3
    assert summary.failed_chunks == 1
    assert summary.missing_chunks == 1
    assert summary.consensus_files == 2
    assert summary.scorch_reranked_files == 1
    assert (
        tmp_path
        / "outputs"
        / "post_docked"
        / "repair_rid"
        / "FULL"
        / "scorch"
        / "_DONE"
    ).exists()
    assert (
        tmp_path
        / "outputs"
        / "post_docked"
        / "repair_rid"
        / "FULL"
        / "LEGACY"
        / "base"
        / "scorch"
        / "_DONE"
    ).exists()
    assert summary.backup_path is not None
    assert summary.backup_path.exists()
    backup = _read_yaml(summary.backup_path)
    assert backup["summary"]["total_proteins_scheduled"] == 99


def test_repair_manifest_dry_run_does_not_write_or_backup(tmp_path: Path) -> None:
    assert yaml is not None
    manifest_path = _build_repair_fixture(tmp_path)
    before = manifest_path.read_text(encoding="utf-8")

    summary = repair_run_manifest("repair_rid", root=tmp_path, dry_run=True)

    assert summary.dry_run is True
    assert summary.backup_path is None
    assert manifest_path.read_text(encoding="utf-8") == before
    assert list(manifest_path.parent.glob("run_manifest.yaml.bak.*")) == []


def test_repair_manifest_auto_requires_scorch_when_partial_scorch_exists(
    tmp_path: Path,
) -> None:
    assert yaml is not None
    run_id = "scorch_auto_rid"
    manifest_run = tmp_path / "outputs" / "manifests" / run_id
    dist = manifest_run / "distributed"
    _write_yaml(
        manifest_run / "run_manifest.yaml",
        {"run_id": run_id, "status": "completed", "proteins": {}},
    )
    _write_json(
        dist / "combo_chunks_0.json",
        {
            "chunks": [
                {
                    "chunk_id": "done_0",
                    "pdb_id": "DONE",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "ligand_bases": ["a"],
                },
                {
                    "chunk_id": "miss_0",
                    "pdb_id": "MISS",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "ligand_bases": ["a"],
                },
            ]
        },
    )
    for chunk_id in ("done_0", "miss_0"):
        _write_json(
            dist / "chunk_results" / f"{chunk_id}.json",
            {"chunk_id": chunk_id, "status": "completed"},
        )
    for pdb_id in ("DONE", "MISS"):
        score = (
            tmp_path
            / "outputs"
            / "docked"
            / run_id
            / pdb_id
            / "consensus_docking_scores.csv"
        )
        score.parent.mkdir(parents=True, exist_ok=True)
        score.write_text("ligand,consensus_score\na,1\n", encoding="utf-8")
    _write_complete_fda_scorch(tmp_path, run_id, "DONE")

    summary = repair_run_manifest(run_id, root=tmp_path, require_scorch=None)

    repaired = _read_yaml(manifest_run / "run_manifest.yaml")
    assert summary.require_scorch is True
    assert summary.targets_completed == 1
    assert summary.targets_running == 1
    assert repaired["proteins"]["DONE|LEGACY|base"]["status"] == "completed"
    assert repaired["proteins"]["MISS|LEGACY|base"]["status"] == "running"
    assert (
        tmp_path
        / "outputs"
        / "post_docked"
        / run_id
        / "DONE"
        / "scorch"
        / "_DONE"
    ).exists()


def test_repair_manifest_dry_run_does_not_clear_complete_scorch_outputs(
    tmp_path: Path,
) -> None:
    assert yaml is not None
    run_id = "dry_complete_scorch"
    manifest_run = tmp_path / "outputs" / "manifests" / run_id
    dist = manifest_run / "distributed"
    _write_yaml(
        manifest_run / "run_manifest.yaml",
        {"run_id": run_id, "status": "running", "proteins": {}},
    )
    _write_json(
        dist / "combo_chunks_0.json",
        {
            "chunks": [
                {
                    "chunk_id": "done_0",
                    "pdb_id": "DONE",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "ligand_bases": ["a"],
                },
            ],
        },
    )
    _write_json(
        dist / "chunk_results" / "done_0.json",
        {"chunk_id": "done_0", "status": "completed"},
    )
    consensus = (
        tmp_path
        / "outputs"
        / "docked"
        / run_id
        / "DONE"
        / "consensus_docking_scores.csv"
    )
    consensus.parent.mkdir(parents=True, exist_ok=True)
    consensus.write_text("ligand,consensus_score\na,1\n", encoding="utf-8")
    _write_complete_fda_scorch(tmp_path, run_id, "DONE")

    summary = repair_run_manifest(
        run_id,
        root=tmp_path,
        dry_run=True,
        require_scorch=True,
        clear_incomplete_scorch_outputs=True,
    )

    assert summary.targets_completed == 1
    assert summary.scorch_output_files_cleared == 0
    assert summary.scorch_done_sentinels_cleared == 0
    assert summary.scorch_coverage_summaries_cleared == 0


def test_repair_manifest_clears_stale_scorch_coordination_state(
    tmp_path: Path,
) -> None:
    assert yaml is not None
    manifest_path = _build_repair_fixture(tmp_path)
    run_id = "repair_rid"
    dist = tmp_path / "outputs" / "manifests" / run_id / "distributed"
    post_root = tmp_path / "outputs" / "post_docked" / run_id

    stale_state = dist / "scorch_scope_state" / "stale_scope.json"
    stale_claim = dist / "scorch_scope_claims" / "stale_scope.claim.json"
    done_state = dist / "scorch_scope_state" / "done_scope.json"
    done_claim = dist / "scorch_scope_claims" / "done_scope.claim.json"
    _write_json(
        stale_state,
        {
            "status": "running",
            "scope_key": "stale_scope",
            "pdb_id": "PART",
            "variant": "LEGACY",
            "ph": "base",
        },
    )
    _write_json(
        stale_claim,
        {
            "scope_key": "stale_scope",
            "pdb_id": "PART",
            "variant": "LEGACY",
            "ph": "base",
        },
    )
    _write_json(
        done_state,
        {
            "status": "completed",
            "scope_key": "done_scope",
            "pdb_id": "FULL",
            "variant": "LEGACY",
            "ph": "base",
        },
    )
    _write_json(
        done_claim,
        {
            "scope_key": "done_scope",
            "pdb_id": "FULL",
            "variant": "LEGACY",
            "ph": "base",
        },
    )

    missing_result = dist / "scorch_shard_results" / "missing_shard.json"
    stale_shard_claim = dist / "scorch_shard_claims" / "missing_shard.claim.json"
    running_result = dist / "scorch_shard_results" / "running_shard.json"
    running_shard_claim = dist / "scorch_shard_claims" / "running_shard.claim.json"
    valid_output = post_root / "FULL" / "LEGACY" / "base" / "scorch_scores_vina_best.part000.csv"
    valid_output.parent.mkdir(parents=True, exist_ok=True)
    valid_output.write_text("Ligand_ID,SCORCH_score\na,1\n", encoding="utf-8")
    partial_output = post_root / "PART" / "LEGACY" / "base" / "scorch_scores_vina_best.part000.csv"
    partial_output.parent.mkdir(parents=True, exist_ok=True)
    partial_output.write_text("Ligand_ID,SCORCH_score\na,1\n", encoding="utf-8")
    partial_reranked = post_root / "PART" / "LEGACY" / "base" / "consensus_reranked_scorch.csv"
    partial_reranked.write_text(
        "ligand,run_mode,rescored_flag,SCORCH_score_used\na,fda,1,1\n",
        encoding="utf-8",
    )
    partial_done = post_root / "PART" / "LEGACY" / "base" / "scorch" / "_DONE"
    partial_done.parent.mkdir(parents=True, exist_ok=True)
    partial_done.write_text("ok\n", encoding="utf-8")
    partial_summary = (
        post_root / "PART" / "LEGACY" / "base" / "scorch" / "scorch_coverage_summary.json"
    )
    partial_summary.write_text("{}\n", encoding="utf-8")
    partial_input = post_root / "PART" / "LEGACY" / "base" / ".scorch_inputs" / "a.pdbqt"
    partial_input.parent.mkdir(parents=True, exist_ok=True)
    partial_input.write_text("MODEL\n", encoding="utf-8")
    valid_result = dist / "scorch_shard_results" / "valid_shard.json"
    _write_json(
        missing_result,
        {
            "status": "failed",
            "shard_id": "missing_shard",
            "output_csv": str(post_root / "PART" / "missing.csv"),
        },
    )
    _write_json(
        stale_shard_claim,
        {
            "shard_id": "missing_shard",
        },
    )
    _write_json(
        running_result,
        {
            "status": "running",
            "shard_id": "running_shard",
            "output_csv": str(post_root / "PART" / "running.csv"),
        },
    )
    _write_json(
        running_shard_claim,
        {
            "shard_id": "running_shard",
        },
    )
    _write_json(
        valid_result,
        {
            "status": "completed",
            "shard_id": "valid_shard",
            "output_csv": str(valid_output),
        },
    )

    summary = repair_run_manifest(run_id, root=tmp_path, require_scorch=True)
    repaired = _read_yaml(manifest_path)

    assert summary.targets_completed == 1
    assert repaired["proteins"]["PART|LEGACY|base"]["status"] == "running"
    assert not stale_state.exists()
    assert not stale_claim.exists()
    assert done_state.exists()
    assert not done_claim.exists()
    assert missing_result.exists()
    assert not stale_shard_claim.exists()
    assert not running_result.exists()
    assert not running_shard_claim.exists()
    assert valid_result.exists()
    assert valid_output.exists()
    assert partial_output.exists()
    assert partial_reranked.exists()
    assert partial_done.exists()
    assert partial_summary.exists()
    assert partial_input.exists()
    assert summary.scorch_output_files_cleared == 0
    assert summary.scorch_done_sentinels_cleared == 0
    assert summary.scorch_coverage_summaries_cleared == 0
    assert summary.scorch_input_dirs_cleared == 0


def test_repair_manifest_repairs_missing_scorch_shard_result_from_part_csv(
    tmp_path: Path,
) -> None:
    assert yaml is not None
    manifest_path = _build_repair_fixture(tmp_path)
    run_id = "repair_rid"
    dist = tmp_path / "outputs" / "manifests" / run_id / "distributed"
    post_root = tmp_path / "outputs" / "post_docked" / run_id
    part_csv = post_root / "PART" / "LEGACY" / "base" / "scorch_scores_vina_best.part001.csv"
    part_csv.parent.mkdir(parents=True, exist_ok=True)
    part_csv.write_text("Ligand_ID,SCORCH_score\nb,1\n", encoding="utf-8")
    plan = dist / "scorch_shard_plans" / "PART.json"
    _write_json(
        plan,
        {
            "combo": ["PART", "LEGACY", "base"],
            "shards": [
                {
                    "shard_id": "part_repair",
                    "combo": ["PART", "LEGACY", "base"],
                    "run_mode": "fda",
                    "chunk_tag": "part001",
                    "output_csv": str(part_csv),
                    "allowed_bases": ["b"],
                    "control_bases": [],
                    "expected_output_rows": 1,
                    "spec": {"source": "vina_best", "stage_dir": "vina_best"},
                }
            ],
        },
    )

    summary = repair_run_manifest(run_id, root=tmp_path, require_scorch=True)
    repaired_result = _read_yaml(manifest_path)["repair"]
    result_path = dist / "scorch_shard_results" / "part_repair.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))

    assert summary.targets_running >= 1
    assert plan.exists()
    assert part_csv.exists()
    assert result["status"] == "completed"
    assert result["repaired_from_output"] is True
    assert repaired_result["scorch_shard_results_repaired"] == 1


def test_repair_manifest_reaggregates_complete_scorch_shards(
    tmp_path: Path,
) -> None:
    assert yaml is not None
    run_id = "scorch_reaggregate_rid"
    manifest_run = tmp_path / "outputs" / "manifests" / run_id
    dist = manifest_run / "distributed"
    _write_yaml(
        manifest_run / "run_manifest.yaml",
        {
            "run_id": run_id,
            "status": "running",
            "proteins": {
                "REAG|LEGACY|base": {
                    "pdb_id": "REAG",
                    "variant": "LEGACY",
                    "ph": "base",
                    "status": "running",
                }
            },
        },
    )
    _write_json(
        dist / "combo_chunks_0.json",
        {
            "chunks": [
                {
                    "chunk_id": "reag_0",
                    "pdb_id": "REAG",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "library_name": "fda+dud",
                    "ligand_bases": ["lig_a", "decoy_a"],
                }
            ]
        },
    )
    _write_json(dist / "chunk_results" / "reag_0.json", {"chunk_id": "reag_0", "status": "completed"})
    dock_root = tmp_path / "outputs" / "docked" / run_id / "REAG"
    dock_root.mkdir(parents=True, exist_ok=True)
    (dock_root / "consensus_docking_scores.csv").write_text(
        "pdb_id,variant,ph_label,ligand,consensus_score\n"
        "REAG,,,lig_a,1\n",
        encoding="utf-8",
    )
    (dock_root / "decoys_consensus_docking_scores.csv").write_text(
        "pdb_id,variant,ph_label,ligand,consensus_score\n"
        "REAG,,,decoy_a,1\n",
        encoding="utf-8",
    )
    post_root = tmp_path / "outputs" / "post_docked" / run_id / "REAG"
    fda_part = post_root / "scorch_scores_vina_best.part000.csv"
    dud_part = post_root / "decoys_scorch_scores_vina_best.part000.csv"
    fda_part.parent.mkdir(parents=True, exist_ok=True)
    fda_part.write_text(
        "Ligand_ID,source,run_mode,SCORCH_score,SCORCH_certainty\n"
        "lig_a,vina,fda,1,1\n",
        encoding="utf-8",
    )
    dud_part.write_text(
        "Ligand_ID,source,run_mode,SCORCH_score,SCORCH_certainty\n"
        "decoy_a,vina,dud,1,1\n",
        encoding="utf-8",
    )
    _write_json(
        dist / "scorch_shard_plans" / "REAG__base__base__decoys.json",
        {
            "combo": ["REAG", "", ""],
            "decoy_prefix": "decoys",
            "shards": [
                {
                    "shard_id": "reag_fda",
                    "combo": ["REAG", "", ""],
                    "run_mode": "fda",
                    "chunk_tag": "part000",
                    "output_csv": str(fda_part),
                    "allowed_bases": ["lig_a"],
                    "control_bases": [],
                    "expected_output_rows": 1,
                    "spec": {
                        "source": "vina",
                        "stage_dir": "vina_best",
                        "output_name": "scorch_scores_vina_best.csv",
                    },
                },
                {
                    "shard_id": "reag_dud",
                    "combo": ["REAG", "", ""],
                    "run_mode": "dud",
                    "chunk_tag": "part000",
                    "output_csv": str(dud_part),
                    "allowed_bases": ["decoy_a"],
                    "control_bases": [],
                    "expected_output_rows": 1,
                    "spec": {
                        "source": "vina",
                        "stage_dir": "vina_best",
                        "output_name": "scorch_scores_vina_best.csv",
                    },
                },
            ],
        },
    )

    summary = repair_run_manifest(run_id, root=tmp_path, require_scorch=True)
    repaired = _read_yaml(manifest_run / "run_manifest.yaml")
    coverage = json.loads(
        (post_root / "scorch" / "scorch_coverage_summary.json").read_text(
            encoding="utf-8"
        )
    )

    assert summary.targets_completed == 1
    assert repaired["proteins"]["REAG|LEGACY|base"]["status"] == "completed"
    assert (post_root / "scorch_scores_all.csv").exists()
    assert (post_root / "decoys_scorch_scores_all.csv").exists()
    assert (post_root / "consensus_reranked_scorch.csv").exists()
    assert (post_root / "decoys_consensus_reranked_scorch.csv").exists()
    assert coverage["all_complete"] is True
    assert (post_root / "scorch" / "_DONE").exists()


def test_repair_manifest_repairs_scorch_selection_semantics(tmp_path: Path) -> None:
    assert yaml is not None
    manifest_path = _build_repair_fixture(tmp_path)
    run_id = "repair_rid"
    dist = tmp_path / "outputs" / "manifests" / run_id / "distributed"
    docked_root = tmp_path / "outputs" / "docked" / run_id
    post_root = tmp_path / "outputs" / "post_docked" / run_id

    docking_artifact = docked_root / "FULL" / "vina_best" / "keep.pdbqt"
    docking_artifact.parent.mkdir(parents=True, exist_ok=True)
    docking_artifact.write_text("MODEL\n", encoding="utf-8")
    fda_consensus = docked_root / "FULL" / "consensus_docking_scores.csv"
    fda_consensus.write_text(
        "ligand,consensus_score\na,1\nb,3\nc,2\n",
        encoding="utf-8",
    )
    dud_consensus = docked_root / "FULL" / "dud_consensus_docking_scores.csv"
    dud_consensus.write_text(
        "ligand,consensus_score\nd1,-5\nd2,-1\n",
        encoding="utf-8",
    )

    combo_post = post_root / "FULL" / "LEGACY" / "base"
    (combo_post / "scorch").mkdir(parents=True, exist_ok=True)
    (combo_post / "scorch" / "_DONE").write_text("ok\n", encoding="utf-8")
    (combo_post / "scorch" / "scorch_coverage_summary.json").write_text(
        "{}\n",
        encoding="utf-8",
    )
    (combo_post / "scorch_scores_vina_best.part000.csv").write_text(
        "Ligand_ID,SCORCH_score\na,1\n",
        encoding="utf-8",
    )
    (combo_post / "dud_scorch_scores_all.csv").write_text(
        "Ligand_ID,SCORCH_score\nd1,1\n",
        encoding="utf-8",
    )
    (combo_post / "dud_consensus_reranked_scorch.csv").write_text(
        "ligand,run_mode,rescored_flag,SCORCH_score_used\nd1,dud,1,1\n",
        encoding="utf-8",
    )
    scorch_input = combo_post / ".scorch_inputs" / "fda" / "vina_best" / "chunk0" / "a.pdbqt"
    scorch_input.parent.mkdir(parents=True, exist_ok=True)
    scorch_input.write_text("MODEL\n", encoding="utf-8")
    _write_json(dist / "scorch_scope_state" / "FULL.json", {"status": "completed"})
    _write_json(dist / "scorch_scope_claims" / "FULL.claim.json", {"status": "completed"})
    _write_json(
        dist / "scorch_shard_plans" / "FULL.json",
        {"combo": ["FULL", "LEGACY", "base"], "shards": [{"shard_id": "shard0"}]},
    )
    _write_json(
        dist / "scorch_shard_results" / "shard0.json",
        {"status": "completed", "shard_id": "shard0"},
    )
    _write_json(dist / "scorch_shard_claims" / "shard0.claim.json", {"shard_id": "shard0"})
    _write_json(dist / "scorch_shard_summary.json", {"all_complete": True})

    summary = repair_run_manifest(
        run_id,
        root=tmp_path,
        require_scorch=True,
        repair_scorch_selection=True,
    )
    repaired = _read_yaml(manifest_path)

    assert summary.repair_scorch_selection is True
    assert summary.require_scorch is True
    assert summary.consensus_files_checked >= 3
    assert summary.consensus_files_sorted >= 2
    assert _read_csv_values(fda_consensus, "ligand") == ["b", "c", "a"]
    assert _read_csv_values(dud_consensus, "ligand") == ["d2", "d1"]
    assert docking_artifact.exists()
    assert repaired["proteins"]["FULL|LEGACY|base"]["status"] == "running"
    assert repaired["proteins"]["FULL|LEGACY|base"]["stages"]["docking"]["status"] == "completed"
    assert not list(post_root.rglob("*scorch_scores*.csv"))
    assert not list(post_root.rglob("*reranked_scorch*.csv"))
    assert not list(post_root.rglob("scorch_coverage_summary.json"))
    assert not list(post_root.rglob(".scorch_inputs"))
    assert not list(post_root.rglob("_DONE"))
    assert not list((dist / "scorch_scope_state").glob("*.json"))
    assert not list((dist / "scorch_scope_claims").glob("*.claim.json"))
    assert not list((dist / "scorch_shard_plans").glob("*.json"))
    assert not list((dist / "scorch_shard_results").glob("*.json"))
    assert not list((dist / "scorch_shard_claims").glob("*.claim.json"))
    assert not (dist / "scorch_shard_summary.json").exists()


def test_auto_repair_resume_manifest_uses_atlas_all_dirs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    assert yaml is not None
    run_id = "auto_resume_rid"
    manifest_run = tmp_path / "outputs" / "manifests" / run_id
    dist = manifest_run / "distributed"
    _write_yaml(
        manifest_run / "run_manifest.yaml",
        {"run_id": run_id, "status": "completed", "proteins": {}},
    )
    _write_json(
        dist / "combo_chunks_0.json",
        {
            "chunks": [
                {
                    "chunk_id": "full_0",
                    "pdb_id": "FULL",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "ligand_bases": ["a"],
                },
                {
                    "chunk_id": "need_0",
                    "pdb_id": "NEED",
                    "variant_label": "LEGACY",
                    "ph_tag": "base",
                    "ligand_bases": ["a"],
                },
            ]
        },
    )
    for chunk_id in ("full_0", "need_0"):
        _write_json(
            dist / "chunk_results" / f"{chunk_id}.json",
            {"chunk_id": chunk_id, "status": "completed"},
        )
    for pdb_id in ("FULL", "NEED"):
        score = (
            tmp_path
            / "outputs"
            / "docked"
            / run_id
            / pdb_id
            / "consensus_docking_scores.csv"
        )
        score.parent.mkdir(parents=True, exist_ok=True)
        score.write_text("ligand,consensus_score\na,1\n", encoding="utf-8")
    _write_complete_fda_scorch(tmp_path, run_id, "FULL")
    partial_need = (
        tmp_path
        / "outputs"
        / "post_docked"
        / run_id
        / "NEED"
        / "LEGACY"
        / "base"
        / "scorch_scores_vina_best.part000.csv"
    )
    partial_need.parent.mkdir(parents=True, exist_ok=True)
    partial_need.write_text("Ligand_ID,SCORCH_score\na,1\n", encoding="utf-8")
    monkeypatch.setenv("ATLAS_ALL_DIRS", str(tmp_path))
    reaggregated: list[str] = []
    real_reaggregate = manifest_repair_mod._reaggregate_scorch_from_shards

    def tracking_reaggregate(*args, **kwargs):
        scope = args[1]
        reaggregated.append(scope.pdb_id)
        return real_reaggregate(*args, **kwargs)

    monkeypatch.setattr(
        manifest_repair_mod,
        "_reaggregate_scorch_from_shards",
        tracking_reaggregate,
    )

    summary = auto_repair_resume_manifest({"USE_SCORCH": True}, run_id)

    repaired = _read_yaml(manifest_run / "run_manifest.yaml")
    assert summary is not None
    assert summary.targets_completed == 1
    assert summary.targets_running == 1
    assert summary.scorch_output_files_cleared == 1
    assert repaired["status"] == "running"
    assert not partial_need.exists()
    assert (
        tmp_path
        / "outputs"
        / "post_docked"
        / run_id
        / "FULL"
        / "scorch"
        / "_DONE"
    ).exists()
    assert "FULL" not in reaggregated
    assert "NEED" in reaggregated


def test_auto_repair_resume_manifest_can_repair_scorch_selection(
    tmp_path: Path,
    monkeypatch,
) -> None:
    assert yaml is not None
    run_id = "repair_rid"
    manifest_path = _build_repair_fixture(tmp_path)
    full_consensus = (
        tmp_path
        / "outputs"
        / "docked"
        / run_id
        / "FULL"
        / "consensus_docking_scores.csv"
    )
    full_consensus.write_text(
        "ligand,consensus_score\na,1\nb,2\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("ATLAS_ALL_DIRS", str(tmp_path))
    monkeypatch.setenv("ATLAS_REPAIR_SCORCH_SELECTION", "1")

    summary = auto_repair_resume_manifest({"USE_SCORCH": True}, run_id)
    repaired = _read_yaml(manifest_path)

    assert summary is not None
    assert summary.repair_scorch_selection is True
    assert _read_csv_values(full_consensus, "ligand") == ["b", "a"]
    assert repaired["proteins"]["FULL|LEGACY|base"]["status"] == "running"


def test_slurm_repair_manifest_cli_accepts_all_dirs(tmp_path: Path, capsys) -> None:
    assert yaml is not None
    _build_repair_fixture(tmp_path)

    rc = _cmd_slurm(
        [
            "repair-manifest",
            "repair_rid",
            "--all-dirs",
            str(tmp_path),
            "--dry-run",
            "--repair-scorch-selection",
        ]
    )

    captured = capsys.readouterr()
    assert rc == 0
    assert "repair_scorch_selection=True" in captured.out
    assert "planned=5" in captured.out
    assert "missing=1" in captured.out


def test_manifest_distributed_mode_infers_slurm_array(monkeypatch) -> None:
    monkeypatch.delenv("ATLAS_DISTRIBUTED_MODE", raising=False)
    monkeypatch.delenv("ATLAS_DIST_TASK_COUNT", raising=False)
    monkeypatch.delenv("ATLAS_DIST_TASK_ID", raising=False)
    monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "8")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "3")

    assert distributed_mode({}) == "slurm_array"
    assert distributed_enabled({}) is True

    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "off")
    assert distributed_mode({}) == "off"
    assert distributed_enabled({}) is False
