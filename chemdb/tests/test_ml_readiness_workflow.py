from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from analysis.cli import run_ml_readiness as readiness_cli
from analysis.cli.run_ml_training_pass import main as training_pass_main
from analysis.ml.leaderboard import write_ml_leaderboard
from analysis.ml.labels import binary_label_series
from analysis.ml.split_manifest import load_locked_split_manifest
from analysis.ml.preflight import build_ml_doctor_report
from analysis.ml.split_locking import lock_run_splits, validate_run_splits
from cli.qol import ml as ml_cli


def _binding_rows(n: int = 96) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for idx in range(n):
        label = int(idx % 3 == 0)
        rows.append(
            {
                "canonical_pair_key": f"D{idx % 24}|T{idx % 12}",
                "drug_id": f"D{idx % 24}",
                "target_id": f"T{idx % 12}",
                "pdb_id": f"P{idx % 8}",
                "spd_binding_label": label,
                "consensus_score": float(label) + (idx % 7) / 20.0,
                "banana_score": float(label) + (idx % 5) / 20.0,
                "banana_binding_probability": 0.2 + 0.6 * label,
                "banana_score_normalized": 0.2 + 0.5 * label,
                "banana_atlas_blend_score": 0.3 + 0.5 * label,
                "binding_expert_score": 0.3 + 0.6 * label,
                "structure_quality": 0.8,
                "protein_class": f"PCLASS{idx % 3}",
                "target_family": f"TFAM{idx % 4}",
                "ligand_chemotype": f"CHEM{idx % 6}",
                "scaffold_key": f"SCAF{idx % 9}",
                "chemical_cluster": f"CLUST{idx % 8}",
                "label_source": ["spd", "chembl", "toxcast"][idx % 3],
                "source_family": ["spd", "chembl", "toxcast"][idx % 3],
                "upstream_source": "synthetic",
                "assay_type": "binding",
                "endpoint_type": "activity",
                "activity_type": "IC50",
                "activity_publication_year": 2012 + (idx % 6),
                "database_release_year": 2020 + (idx % 3),
                "smiles": ["CCO", "CCN", "c1ccccc1"][idx % 3],
                "rdkit_mol_wt": 100 + idx % 20,
                "rdkit_mol_logp": (idx % 8) / 3.0,
                "rdkit_tpsa": 20 + idx % 30,
                "rdkit_hbd": idx % 3,
                "rdkit_hba": idx % 5,
                "rdkit_rotatable_bonds": idx % 6,
                "rdkit_formal_charge": 0,
                "rdkit_aromatic_rings": idx % 2,
                "rdkit_fraction_csp3": (idx % 5) / 5.0,
                "rdkit_qed": 0.4 + (idx % 4) / 10.0,
            }
        )
    return rows


def _write_run(tmp_path: Path) -> Path:
    run_dir = tmp_path / "outputs" / "data" / "R1"
    dataset = run_dir / "spd_four_experts" / "model_ready" / "ml_spd_binding_model_ready.csv"
    dataset.parent.mkdir(parents=True)
    pd.DataFrame(_binding_rows()).to_csv(dataset, index=False)
    return run_dir


def test_doctor_split_lock_train_and_leaderboard_contract(tmp_path):
    run_dir = _write_run(tmp_path)
    ml_root = run_dir / "ml"

    doctor = build_ml_doctor_report(run_id="R1", run_dir=run_dir)
    assert doctor["status"] == "ready_with_warnings"
    assert not any("binding" in blocker for blocker in doctor["blockers"])

    lock = lock_run_splits(run_id="R1", run_dir=run_dir, out_dir=ml_root / "splits", experts=["binding"], split_modes=["drug_holdout"], seed=3)
    assert lock["status"] == "written"
    split_manifest = ml_root / "splits" / "binding" / "drug_holdout"
    assert (split_manifest / "split_manifest.json").exists()

    validation = validate_run_splits(run_dir=run_dir, split_root=ml_root / "splits")
    assert validation["status"] == "valid"

    rc = training_pass_main(
        [
            "--repo-root", str(tmp_path),
            "--run-id", "R1",
            "--run-dir", str(run_dir),
            "--out-dir", str(ml_root / "locked_train"),
            "--experts", "binding",
            "--split-manifest", str(split_manifest),
            "--bootstraps", "0",
            "--skip-audit",
            "--no-mlflow",
        ]
    )
    assert rc == 0
    trained_manifest = json.loads((ml_root / "locked_train" / "ml_training_pass_manifest.json").read_text())
    assert trained_manifest["experts"][0]["split_manifest_input"] == str(split_manifest)

    leaderboard = write_ml_leaderboard(ml_root)
    assert leaderboard["rows"] >= 1
    board = pd.read_csv(ml_root / "leaderboard" / "model_leaderboard.csv")
    assert "model" in set(board["kind"])


def test_atlas_ml_new_wrapper_commands_forward(monkeypatch, tmp_path):
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(ml_cli._bindings, "repo_root", lambda: tmp_path)
    monkeypatch.setattr(ml_cli, "_run_module_main", lambda module, forwarded: calls.append((module, forwarded)) or 0)

    assert ml_cli._cmd_ml(["prepare-run", "--run-id", "R1", "--sample-rows", "123"]) == 0
    assert calls[-1][0] == "analysis.cli.run_ml_readiness"
    assert "--sample-rows" in calls[-1][1]

    assert ml_cli._cmd_ml(["splits", "lock", "--run-id", "R1", "--splits", "drug_holdout"]) == 0
    assert calls[-1][0] == "analysis.cli.run_ml_split_lock"
    assert calls[-1][1][0] == "lock"

    assert ml_cli._cmd_ml(["leaderboard", "--run-id", "R1"]) == 0
    assert calls[-1][0] == "analysis.cli.build_ml_leaderboard"

    assert ml_cli._cmd_ml(["external-eval", "--run-id", "R1", "--model-dir", "m", "--dataset", "d.csv", "--label", "y", "--name", "bench"]) == 0
    assert calls[-1][0] == "analysis.cli.run_ml_external_eval"


def test_prepare_run_blocks_when_no_default_datasets(tmp_path):
    run_dir = tmp_path / "outputs" / "data" / "EMPTY"
    run_dir.mkdir(parents=True)
    ml_root = run_dir / "ml"

    rc = readiness_cli.main(
        [
            "--repo-root", str(tmp_path),
            "--run-id", "EMPTY",
            "--run-dir", str(run_dir),
            "--ml-root", str(ml_root),
        ]
    )

    assert rc == 1
    manifest = json.loads((ml_root / "readiness" / "ml_readiness_manifest.json").read_text())
    assert manifest["status"] == "blocked"
    doctor = json.loads((ml_root / "readiness" / "ml_doctor_report.json").read_text())
    assert any("no default ML datasets" in blocker for blocker in doctor["blockers"])


def test_prepare_run_blocks_when_sample_training_manifest_failed(monkeypatch, tmp_path):
    run_dir = _write_run(tmp_path)
    ml_root = run_dir / "ml"

    def fake_training_main(forwarded: list[str]) -> int:
        sample_out = Path(forwarded[forwarded.index("--out-dir") + 1])
        sample_out.mkdir(parents=True)
        (sample_out / "ml_training_pass_manifest.json").write_text(
            json.dumps({"experts": [{"expert": "binding", "status": "failed", "error": "boom"}]}),
            encoding="utf-8",
        )
        return 0

    monkeypatch.setattr(readiness_cli, "training_main", fake_training_main)
    rc = readiness_cli.main(
        [
            "--repo-root", str(tmp_path),
            "--run-id", "R1",
            "--run-dir", str(run_dir),
            "--ml-root", str(ml_root),
            "--skip-splits",
        ]
    )

    assert rc == 1
    manifest = json.loads((ml_root / "readiness" / "ml_readiness_manifest.json").read_text())
    assert manifest["status"] == "blocked"
    sample_step = next(step for step in manifest["steps"] if step["name"] == "sample_train")
    assert sample_step["status"] == "failed"
    assert sample_step["failed_experts"] == ["binding"]


def test_prepare_run_wires_default_locked_split_to_sample_training(monkeypatch, tmp_path):
    run_dir = _write_run(tmp_path)
    ml_root = run_dir / "ml"
    calls: list[list[str]] = []

    def fake_training_main(forwarded: list[str]) -> int:
        calls.append(forwarded)
        sample_out = Path(forwarded[forwarded.index("--out-dir") + 1])
        sample_out.mkdir(parents=True)
        (sample_out / "ml_training_pass_manifest.json").write_text(
            json.dumps({"experts": [{"expert": "binding", "status": "trained"}]}),
            encoding="utf-8",
        )
        return 0

    monkeypatch.setattr(readiness_cli, "training_main", fake_training_main)
    rc = readiness_cli.main(
        [
            "--repo-root", str(tmp_path),
            "--run-id", "R1",
            "--run-dir", str(run_dir),
            "--ml-root", str(ml_root),
            "--sample-rows", "24",
        ]
    )

    assert rc == 0
    assert calls
    forwarded = calls[-1]
    assert "--expert-split-manifest" in forwarded
    split_arg = forwarded[forwarded.index("--expert-split-manifest") + 1]
    assert split_arg.startswith("binding=")
    assert split_arg.endswith("/splits/binding/drug_holdout")
    manifest = json.loads((ml_root / "readiness" / "ml_readiness_manifest.json").read_text())
    assert any(step["name"] == "locked_split_training_overrides" for step in manifest["steps"])


def test_locked_split_manifest_rejects_wrong_split_mode(tmp_path):
    run_dir = _write_run(tmp_path)
    split_root = tmp_path / "splits"
    lock_run_splits(run_id="R1", run_dir=run_dir, out_dir=split_root, experts=["binding"], split_modes=["drug_holdout"], seed=3)
    metadata = json.loads((split_root / "binding" / "drug_holdout" / "split_manifest.json").read_text())
    locked_dataset = Path(metadata["split_summary"]["locked_dataset"])
    frame = pd.read_csv(locked_dataset)
    labels = binary_label_series(frame["spd_binding_label"])
    frame["_atlas_observed_label"] = labels
    frame = frame.loc[labels.notna()].copy()
    frame["spd_binding_label"] = labels.loc[labels.notna()].astype(int)

    with pytest.raises(ValueError, match="mode mismatch"):
        load_locked_split_manifest(
            frame,
            split_root / "binding" / "drug_holdout",
            label_col="spd_binding_label",
            split_mode="target_holdout",
        )


def test_default_split_lock_includes_source_holdout_when_feasible(tmp_path):
    run_dir = _write_run(tmp_path)
    lock = lock_run_splits(run_id="R1", run_dir=run_dir, out_dir=tmp_path / "splits", experts=["binding"], seed=5)
    rows = [row for row in lock["results"] if row.get("expert") == "binding" and row.get("split_mode") == "source_holdout"]
    assert rows
    assert rows[0]["status"] == "locked"


def test_training_pass_derives_sampled_locked_manifest(tmp_path):
    run_dir = _write_run(tmp_path)
    ml_root = run_dir / "ml"
    lock_run_splits(run_id="R1", run_dir=run_dir, out_dir=ml_root / "splits", experts=["binding"], split_modes=["drug_holdout"], seed=3)
    split_manifest = ml_root / "splits" / "binding" / "drug_holdout"

    rc = training_pass_main(
        [
            "--repo-root", str(tmp_path),
            "--run-id", "R1",
            "--run-dir", str(run_dir),
            "--out-dir", str(ml_root / "sample_locked_train"),
            "--experts", "binding",
            "--split-manifest", str(split_manifest),
            "--sample-rows", "48",
            "--bootstraps", "0",
            "--skip-audit",
            "--no-mlflow",
        ]
    )

    assert rc == 0
    trained_manifest = json.loads((ml_root / "sample_locked_train" / "ml_training_pass_manifest.json").read_text())
    expert = trained_manifest["experts"][0]
    assert expert["status"] == "trained"
    sample_manifest = expert["sample_manifest"]
    assert sample_manifest["split_locked"] is True
    assert sample_manifest["source_split_manifest"] == str(split_manifest)
    assert expert["split_manifest_input"] == sample_manifest["sampled_split_manifest"]

