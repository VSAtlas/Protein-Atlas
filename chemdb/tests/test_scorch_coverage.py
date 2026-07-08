from __future__ import annotations

import csv
import logging
from pathlib import Path

from cli.qol.manifest_repair import repair_run_manifest
from post_docking.rescoring.rescore_reranker import rerank_consensus_with_scorch
from post_docking.rescoring.scorch_coverage import (
    evaluate_scorch_coverage,
    write_coverage_summary,
)

try:
    import yaml  # type: ignore[import-untyped]
except Exception:  # pragma: no cover
    yaml = None  # type: ignore[assignment]


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _write_yaml(path: Path, payload: object) -> None:
    assert yaml is not None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def _write_manifest_plan(root: Path, run_id: str, pdb_id: str) -> None:
    manifest_dir = root / "outputs" / "manifests" / run_id
    _write_yaml(
        manifest_dir / "run_manifest.yaml",
        {"run_id": run_id, "status": "running", "proteins": {}},
    )
    plan = {
        "chunks": [
            {
                "chunk_id": f"{pdb_id.lower()}_0",
                "pdb_id": pdb_id,
                "variant_label": "LEGACY",
                "ph_tag": "base",
                "ligand_bases": ["a", "b"],
            }
        ]
    }
    (manifest_dir / "distributed").mkdir(parents=True, exist_ok=True)
    (manifest_dir / "distributed" / "combo_chunks_0.json").write_text(
        __import__("json").dumps(plan),
        encoding="utf-8",
    )
    (manifest_dir / "distributed" / "chunk_results").mkdir(parents=True, exist_ok=True)
    (manifest_dir / "distributed" / "chunk_results" / f"{pdb_id.lower()}_0.json").write_text(
        __import__("json").dumps({"chunk_id": f"{pdb_id.lower()}_0", "status": "completed"}),
        encoding="utf-8",
    )


def test_coverage_requires_fda_and_dud_raw_and_final_rows(tmp_path: Path) -> None:
    run_id = "cov_rid"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "consensus_docking_scores.csv",
        [{"ligand": f"a{i}.pdbqt", "consensus_score": "1"} for i in range(10)],
    )
    _write_csv(
        docked / "dud_consensus_docking_scores.csv",
        [{"ligand": f"d{i}.pdbqt", "consensus_score": "1"} for i in range(10)],
    )
    _write_csv(
        post / "scorch_scores_all.csv",
        [{"run_mode": "fda", "Ligand_ID": "a0", "SCORCH_score": "1"}],
    )
    _write_csv(
        post / "consensus_reranked_scorch.csv",
        [{"ligand": "a0.pdbqt", "run_mode": "fda", "rescored_flag": "1"}],
    )
    (post / "scorch" / "_DONE").parent.mkdir(parents=True, exist_ok=True)
    (post / "scorch" / "_DONE").write_text("ok\n", encoding="utf-8")

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=0.10,
    )
    write_coverage_summary(tmp_path / "outputs" / "post_docked" / run_id, summary)

    assert summary.fda.complete is True
    assert summary.dud.complete is False
    assert summary.all_complete is False
    assert "dud_raw_under_coverage:0/1" in summary.reasons
    assert (post / "scorch" / "scorch_coverage_summary.json").exists()


def test_coverage_rejects_quarantined_placeholder_rows(tmp_path: Path) -> None:
    run_id = "cov_quarantine_placeholder"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "consensus_docking_scores.csv",
        [
            {"ligand": "a0.pdbqt", "run_mode": "fda", "consensus_score": "10"},
            {"ligand": "a1.pdbqt", "run_mode": "fda", "consensus_score": "9"},
        ],
    )
    _write_csv(
        post / "scorch_scores_all.csv",
        [
            {
                "Ligand_ID": "a0",
                "run_mode": "fda",
                "SCORCH_score": "0",
                "scorch_status": "quarantined",
                "scorch_failure_reason": "timeout",
                "scorch_returncode": "124",
                "scorch_unscorable_flag": "1",
            }
        ],
    )
    _write_csv(
        post / "consensus_reranked_scorch.csv",
        [
            {
                "ligand": "a0.pdbqt",
                "run_mode": "fda",
                "rescored_flag": "1",
                "SCORCH_score_used": "0",
                "scorch_status": "quarantined",
                "scorch_failure_reason": "timeout",
                "scorch_returncode": "124",
                "scorch_unscorable_flag": "1",
            }
        ],
    )

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=0.50,
    )

    assert summary.fda.expected_rows == 1
    assert summary.fda.raw_rows == 0
    assert summary.fda.final_rescored_rows == 0
    assert summary.fda.invalid_rows == 2
    assert summary.fda.complete is False
    assert summary.all_complete is False
    assert "fda_raw_under_coverage:0/1" in summary.reasons
    assert "fda_final_under_coverage:0/1" in summary.reasons
    assert "fda_invalid_scorch_rows:2" in summary.reasons


def test_coverage_accepts_explicit_quarantine_when_enabled(tmp_path: Path) -> None:
    run_id = "cov_accept_quarantine"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "consensus_docking_scores.csv",
        [
            {"ligand": "a0.pdbqt", "run_mode": "fda", "consensus_score": "10"},
            {"ligand": "a1.pdbqt", "run_mode": "fda", "consensus_score": "9"},
        ],
    )
    quarantine_row = {
        "Ligand_ID": "a0",
        "ligand": "a0.pdbqt",
        "run_mode": "fda",
        "SCORCH_score": "0",
        "rescored_flag": "1",
        "scorch_status": "quarantined",
        "scorch_failure_reason": "explicit_ligand_quarantine",
        "scorch_returncode": "",
        "scorch_unscorable_flag": "1",
    }
    _write_csv(post / "scorch_scores_all.csv", [quarantine_row])
    _write_csv(post / "consensus_reranked_scorch.csv", [quarantine_row])

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=0.50,
        accept_quarantined=True,
    )

    assert summary.fda.expected_rows == 1
    assert summary.fda.raw_rows == 1
    assert summary.fda.final_rescored_rows == 1
    assert summary.fda.invalid_rows == 0
    assert summary.fda.quarantined_rows == 2
    assert summary.fda.complete is True
    assert summary.all_complete is True
    assert summary.reasons == ()


def test_coverage_accepts_global_timeout_quarantine_by_default(
    tmp_path: Path,
) -> None:
    run_id = "cov_accept_global_timeout_quarantine"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "dud_consensus_docking_scores.csv",
        [
            {"ligand": "d0.pdbqt", "run_mode": "dud", "consensus_score": "10"},
            {"ligand": "d1.pdbqt", "run_mode": "dud", "consensus_score": "9"},
        ],
    )
    quarantine_row = {
        "Ligand_ID": "d0",
        "ligand": "d0.pdbqt",
        "run_mode": "dud",
        "SCORCH_score": "0",
        "rescored_flag": "1",
        "scorch_status": "quarantined",
        "scorch_failure_reason": "global_ligand_timeout",
        "scorch_returncode": "124",
        "scorch_unscorable_flag": "1",
    }
    _write_csv(post / "dud_scorch_scores_all.csv", [quarantine_row])
    _write_csv(post / "consensus_reranked_scorch.csv", [quarantine_row])

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=0.50,
    )

    assert summary.dud.expected_rows == 1
    assert summary.dud.raw_rows == 1
    assert summary.dud.final_rescored_rows == 1
    assert summary.dud.invalid_rows == 0
    assert summary.dud.quarantined_rows == 2
    assert summary.dud.complete is True
    assert summary.all_complete is True
    assert summary.reasons == ()


def test_coverage_accepts_controls_only_consensus_with_zero_expected(
    tmp_path: Path,
) -> None:
    run_id = "cov_controls_only"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "consensus_docking_scores.csv",
        [{"ligand": "ctrl_fda.pdbqt", "consensus_score": "1"}],
    )
    _write_csv(
        docked / "dud_consensus_docking_scores.csv",
        [{"ligand": "ctrl_dud.pdbqt", "consensus_score": "1"}],
    )
    _write_csv(
        post / "scorch_scores_all.csv",
        [{"run_mode": "fda", "Ligand_ID": "ctrl_fda", "SCORCH_score": "1"}],
    )
    _write_csv(
        post / "dud_scorch_scores_all.csv",
        [{"run_mode": "dud", "Ligand_ID": "ctrl_dud", "SCORCH_score": "1"}],
    )
    _write_csv(
        post / "consensus_reranked_scorch.csv",
        [{"ligand": "ctrl_fda.pdbqt", "run_mode": "fda", "rescored_flag": "1"}],
    )
    _write_csv(
        post / "dud_consensus_reranked_scorch.csv",
        [{"ligand": "ctrl_dud.pdbqt", "run_mode": "dud", "rescored_flag": "1"}],
    )

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=1.0,
        control_bases={"ctrl_fda", "ctrl_dud"},
    )

    assert summary.fda.expected_rows == 0
    assert summary.dud.expected_rows == 0
    assert summary.all_complete is True
    assert "missing_docking_consensus" not in summary.reasons


def test_coverage_uses_expected_ligand_set_when_raw_mode_is_wrong(tmp_path: Path) -> None:
    run_id = "cov_dud_raw_mode"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "consensus_docking_scores.csv",
        [
            {
                "ligand": "actives_final_00001.pdbqt",
                "run_mode": "dud",
                "library": "dud",
                "is_decoy": "0",
                "consensus_score": "1",
            },
            {
                "ligand": "decoys_final_00001.pdbqt",
                "run_mode": "dud",
                "library": "dud",
                "is_decoy": "1",
                "consensus_score": "1",
            },
        ],
    )
    _write_csv(
        post / "scorch_scores_all.csv",
        [
            {"run_mode": "fda", "Ligand_ID": "actives_final_00001_stage1", "SCORCH_score": "1"},
            {"run_mode": "fda", "Ligand_ID": "decoys_final_00001_stage1", "SCORCH_score": "1"},
        ],
    )
    _write_csv(
        post / "consensus_reranked_scorch.csv",
        [
            {
                "ligand": "actives_final_00001.pdbqt",
                "run_mode": "dud",
                "library": "dud",
                "is_decoy": "0",
                "rescored_flag": "1",
            },
            {
                "ligand": "decoys_final_00001.pdbqt",
                "run_mode": "dud",
                "library": "dud",
                "is_decoy": "1",
                "rescored_flag": "1",
            },
        ],
    )

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=1.0,
    )

    assert summary.fda.expected_rows == 0
    assert summary.dud.expected_rows == 2
    assert summary.dud.raw_rows == 2
    assert summary.dud.final_rescored_rows == 2
    assert summary.all_complete is True


def test_coverage_normalizes_custom_library_stage_suffixes(tmp_path: Path) -> None:
    run_id = "cov_custom_prefix"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "test_library_20_consensus_docking_scores.csv",
        [
            {
                "ligand": "19A_A360.sanitized.pdbqt",
                "run_mode": "test_library_20",
                "library": "test_library_20",
                "is_decoy": "0",
                "consensus_score": "1",
            }
        ],
    )
    _write_csv(
        post / "test_library_20_scorch_scores_all.csv",
        [
            {
                "run_mode": "dud",
                "Ligand_ID": "19A_A360.sanitized_test_library_20_stage3",
                "SCORCH_score": "1",
            }
        ],
    )
    _write_csv(
        post / "test_library_20_consensus_reranked_scorch.csv",
        [
            {
                "ligand": "19A_A360.sanitized_test_library_20_stage3",
                "run_mode": "test_library_20",
                "rescored_flag": "1",
            }
        ],
    )

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=1.0,
        decoy_prefix="test_library_20",
    )

    assert summary.dud.expected_rows == 1
    assert summary.dud.raw_rows == 1
    assert summary.dud.final_rescored_rows == 1
    assert summary.all_complete is True


def test_coverage_prefers_shard_plan_selected_ligands(tmp_path: Path) -> None:
    run_id = "cov_shard_plan"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    manifest = tmp_path / "outputs" / "manifests" / run_id
    _write_csv(
        docked / "dud_consensus_docking_scores.csv",
        [
            {"ligand": "d0.pdbqt", "run_mode": "dud", "consensus_score": "1"},
            {"ligand": "d1.pdbqt", "run_mode": "dud", "consensus_score": "1"},
            {"ligand": "d2.pdbqt", "run_mode": "dud", "consensus_score": "1"},
            {"ligand": "d3.pdbqt", "run_mode": "dud", "consensus_score": "1"},
        ],
    )
    _write_csv(
        post / "dud_scorch_scores_all.csv",
        [
            {"run_mode": "dud", "Ligand_ID": "d2_dud_stage1", "SCORCH_score": "1"},
            {"run_mode": "dud", "Ligand_ID": "d3_dud_stage2", "SCORCH_score": "1"},
            {"run_mode": "dud", "Ligand_ID": "extra_dud_stage1", "SCORCH_score": "1"},
        ],
    )
    _write_csv(
        post / "dud_consensus_reranked_scorch.csv",
        [
            {"ligand": "d2_dud_stage1", "run_mode": "dud", "rescored_flag": "1"},
            {"ligand": "d3_dud_stage2", "run_mode": "dud", "rescored_flag": "1"},
            {"ligand": "extra_dud_stage1", "run_mode": "dud", "rescored_flag": "1"},
        ],
    )
    plan_dir = manifest / "distributed" / "scorch_shard_plans"
    plan_dir.mkdir(parents=True, exist_ok=True)
    (plan_dir / "P1__base__base__dud.json").write_text(
        __import__("json").dumps(
            {
                "schema": "atlas.scorch_shard_plan.v1",
                "run_id": run_id,
                "combo": ["P1", "", ""],
                "decoy_prefix": "dud",
                "shards": [
                    {
                        "run_mode": "dud",
                        "allowed_bases": ["d2", "d3"],
                        "control_bases": [],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=0.50,
        manifest_run_root=manifest,
    )

    assert summary.dud.expected_rows == 2
    assert summary.dud.raw_rows == 2
    assert summary.dud.final_rescored_rows == 2
    assert summary.all_complete is True


def test_coverage_uses_same_combo_shard_plan_when_prefix_alias_differs(
    tmp_path: Path,
) -> None:
    run_id = "cov_shard_plan_prefix_alias"
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    manifest = tmp_path / "outputs" / "manifests" / run_id
    _write_csv(
        docked / "consensus_docking_scores.csv",
        [
            {"ligand": "FEX_A1.pdbqt", "run_mode": "fda", "consensus_score": "99"},
            {"ligand": "f0.pdbqt", "run_mode": "fda", "consensus_score": "1"},
        ],
    )
    _write_csv(
        docked / "dud_consensus_docking_scores.csv",
        [{"ligand": "d0.pdbqt", "run_mode": "dud", "consensus_score": "1"}],
    )
    _write_csv(
        post / "scorch_scores_all.csv",
        [{"run_mode": "fda", "Ligand_ID": "f0", "SCORCH_score": "1"}],
    )
    _write_csv(
        post / "dud_scorch_scores_all.csv",
        [{"run_mode": "dud", "Ligand_ID": "d0", "SCORCH_score": "1"}],
    )
    _write_csv(
        post / "consensus_reranked_scorch.csv",
        [
            {"ligand": "f0.pdbqt", "run_mode": "fda", "rescored_flag": "1"},
            {"ligand": "d0.pdbqt", "run_mode": "dud", "rescored_flag": "1"},
        ],
    )
    plan_dir = manifest / "distributed" / "scorch_shard_plans"
    plan_dir.mkdir(parents=True, exist_ok=True)
    (plan_dir / "P1__base__base__dud.json").write_text(
        __import__("json").dumps(
            {
                "schema": "atlas.scorch_shard_plan.v1",
                "run_id": run_id,
                "combo": ["P1", "", ""],
                "decoy_prefix": "dud",
                "shards": [
                    {
                        "run_mode": "fda",
                        "allowed_bases": ["FEX_A1", "f0"],
                        "control_bases": ["FEX_A1"],
                    },
                    {
                        "run_mode": "dud",
                        "allowed_bases": ["d0"],
                        "control_bases": [],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=0.50,
        decoy_prefix="dud_fda",
        manifest_run_root=manifest,
    )

    assert summary.fda.expected_rows == 1
    assert summary.fda.raw_rows == 1
    assert summary.fda.final_rescored_rows == 1
    assert summary.dud.expected_rows == 1
    assert summary.dud.raw_rows == 1
    assert summary.dud.final_rescored_rows == 1
    assert summary.all_complete is True
    assert summary.reasons == ()


def test_coverage_requires_docking_consensus_rows(tmp_path: Path) -> None:
    run_id = "cov_missing_consensus"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        post / "consensus_reranked_scorch.csv",
        [{"ligand": "a.pdbqt", "run_mode": "fda", "rescored_flag": "1"}],
    )

    summary = evaluate_scorch_coverage(
        docked_run_root=tmp_path / "outputs" / "docked" / run_id,
        post_run_root=tmp_path / "outputs" / "post_docked" / run_id,
        pdb_id="P1",
        top_fraction=1.0,
    )

    assert summary.all_complete is False
    assert "missing_docking_consensus" in summary.reasons


def test_manifest_repair_requires_coverage_and_clears_stale_done(tmp_path: Path) -> None:
    assert yaml is not None
    run_id = "repair_cov"
    _write_manifest_plan(tmp_path, run_id, "P1")
    docked = tmp_path / "outputs" / "docked" / run_id / "P1"
    post = tmp_path / "outputs" / "post_docked" / run_id / "P1"
    _write_csv(
        docked / "consensus_docking_scores.csv",
        [{"ligand": "a.pdbqt", "consensus_score": "1"}],
    )
    _write_csv(
        post / "consensus_reranked_scorch.csv",
        [{"ligand": "a.pdbqt", "run_mode": "fda", "rescored_flag": "0"}],
    )
    (post / "scorch" / "_DONE").parent.mkdir(parents=True, exist_ok=True)
    (post / "scorch" / "_DONE").write_text("ok\n", encoding="utf-8")
    dist = tmp_path / "outputs" / "manifests" / run_id / "distributed"
    stale_plan = dist / "scorch_shard_plans" / "P1__base__base__dud.json"
    stale_result = dist / "scorch_shard_results" / "scorch_stale.json"
    stale_claim = dist / "scorch_shard_claims" / "scorch_stale.claim.json"
    stale_plan.parent.mkdir(parents=True, exist_ok=True)
    stale_result.parent.mkdir(parents=True, exist_ok=True)
    stale_claim.parent.mkdir(parents=True, exist_ok=True)
    stale_plan.write_text(
        __import__("json").dumps(
            {
                "combo": ["P1", "", ""],
                "shards": [{"shard_id": "scorch_stale"}],
            }
        ),
        encoding="utf-8",
    )
    stale_result.write_text(
        __import__("json").dumps(
            {
                "shard_id": "scorch_stale",
                "status": "completed",
                "output_csv": str(post / "scorch_scores_all.csv"),
            }
        ),
        encoding="utf-8",
    )
    stale_claim.write_text("{}", encoding="utf-8")

    summary = repair_run_manifest(run_id, root=tmp_path, require_scorch=True)

    assert summary.targets_completed == 0
    assert summary.targets_running == 1
    assert not (post / "scorch" / "_DONE").exists()
    assert stale_plan.exists()
    assert stale_result.exists()
    assert not stale_claim.exists()


def test_reranker_matches_legacy_empty_variant_and_ph(tmp_path: Path) -> None:
    consensus = tmp_path / "consensus_docking_scores.csv"
    scorch = tmp_path / "scorch_scores_all.csv"
    out = tmp_path / "consensus_reranked_scorch.csv"
    _write_csv(
        consensus,
        [
            {
                "run_id": "r1",
                "pdb_id": "1B41",
                "variant": "legacy",
                "ph_label": "",
                "ligand": "rdk_0000018.pdbqt",
                "consensus_score": "1.0",
                "best_engine": "vina",
                "run_mode": "fda",
            }
        ],
    )
    _write_csv(
        scorch,
        [
            {
                "pdb_id": "1B41",
                "variant": "",
                "ph": "",
                "source": "vina",
                "run_mode": "fda",
                "Ligand_ID": "rdk_0000018_stage2",
                "SCORCH_score": "0.25",
                "SCORCH_certainty": "0.9",
            }
        ],
    )

    ok = rerank_consensus_with_scorch(
        consensus,
        scorch,
        out,
        logging.getLogger("test_reranker_matches_legacy"),
        overwrite=True,
    )

    assert ok
    with out.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["rescored_flag"] == "1"
    assert row["SCORCH_score_used"] == "0.25"


def test_reranker_infers_scope_for_metadata_less_scorch_rows(tmp_path: Path) -> None:
    consensus = tmp_path / "consensus_docking_scores.csv"
    scorch = tmp_path / "scorch_scores_all.csv"
    out = tmp_path / "consensus_reranked_scorch.csv"
    _write_csv(
        consensus,
        [
            {
                "run_id": "r1",
                "pdb_id": "5WIU",
                "variant": "legacy",
                "ph_label": "",
                "ligand": "rdk_0000443.pdbqt",
                "consensus_score": "1.0",
                "best_engine": "vina",
                "run_mode": "fda",
            }
        ],
    )
    _write_csv(
        scorch,
        [
            {
                "Ligand_ID": "rdk_0000443_stage1",
                "SCORCH_score": "0.04432",
                "SCORCH_certainty": "0.96197",
            }
        ],
    )

    ok = rerank_consensus_with_scorch(
        consensus,
        scorch,
        out,
        logging.getLogger("test_reranker_infers_scope"),
        overwrite=True,
    )

    assert ok
    with out.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["rescored_flag"] == "1"
    assert row["SCORCH_score_used"] == "0.04432"
