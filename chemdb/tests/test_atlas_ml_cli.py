from __future__ import annotations

from pathlib import Path

from cli.qol import ml as ml_cli


def test_atlas_ml_train_forwards_outputs_data_defaults(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(ml_cli._bindings, "repo_root", lambda: tmp_path)
    monkeypatch.setattr(ml_cli, "_run_module_main", lambda module, forwarded: calls.append((module, forwarded)) or 0)

    assert ml_cli._cmd_ml(["train", "--run-id", "R1", "--bootstraps", "0"]) == 0

    module, forwarded = calls[0]
    assert module == "analysis.cli.run_ml_training_pass"
    assert "--repo-root" in forwarded
    assert str(tmp_path / "outputs" / "data" / "R1") in forwarded
    assert str(tmp_path / "outputs" / "data" / "R1" / "ml" / "training_pass") in forwarded
    assert "--bootstraps" in forwarded


def test_atlas_ml_audit_and_source_pu_forward(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(ml_cli._bindings, "repo_root", lambda: tmp_path)
    monkeypatch.setattr(ml_cli, "_run_module_main", lambda module, forwarded: calls.append((module, forwarded)) or 0)

    assert ml_cli._cmd_ml(["audit", "--dataset", "d.csv", "--label", "y", "--splits", "random"]) == 0
    assert calls[-1][0] == "analysis.cli.run_ml_audit_suite"
    assert str(tmp_path / "outputs" / "data" / "pilotstudy" / "ml" / "audit_suite") in calls[-1][1]
    assert "--splits" in calls[-1][1]

    assert ml_cli._cmd_ml(["source-pu", "--run-id", "R2", "--skip-pu-suite"]) == 0
    assert calls[-1][0] == "analysis.cli.run_ml_source_pu_suite"
    assert str(tmp_path / "outputs" / "data" / "R2" / "ml" / "source_pu_suite") in calls[-1][1]
    assert "--skip-pu-suite" in calls[-1][1]


def test_atlas_ml_status_json_reads_manifests(tmp_path, capsys):
    ml_root = tmp_path / "ml"
    manifest = ml_root / "training_pass" / "ml_training_pass_manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text('{"experts":[{"expert":"binding","status":"trained"}]}')

    rc = ml_cli._cmd_ml_status(type("Args", (), {"run_id": "R1", "out_dir": str(ml_root), "json_out": True})(), Path("."))
    assert rc == 0
    assert "binding" in capsys.readouterr().out
