from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from cli import qol_cli
from cli.qol import ligands as ligand_cli
from prep_ligands import fda_identity_audit


def _audit_result(output_dir: Path, *, blockers: int = 0) -> SimpleNamespace:
    return SimpleNamespace(
        row_audit_csv=output_dir / "fda_identity_audit.csv",
        evidence_jsonl=output_dir / "fda_identity_evidence.jsonl",
        summary_json=output_dir / "fda_identity_summary.json",
        review_csv=output_dir / "fda_identity_review.csv",
        quarantine_manifest_csv=output_dir / "fda_identity_quarantine.csv",
        eligible_manifest_csv=output_dir / "fda_docking_eligible.csv",
        source_catalog_csv=output_dir / "fda_source_substances.csv",
        total_rows=3,
        blocker_count=blockers,
        review_count=1,
        has_blockers=blockers > 0,
    )


def test_fda_audit_uses_repo_defaults_and_reports_artifacts(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    captured: list[fda_identity_audit.FDAIdentityAuditConfig] = []
    expected_out = tmp_path / "outputs" / "data" / "fda_identity_audit"

    def _run(config: fda_identity_audit.FDAIdentityAuditConfig) -> SimpleNamespace:
        captured.append(config)
        return _audit_result(config.output_dir)

    monkeypatch.setattr(ligand_cli._bindings, "repo_root", lambda: tmp_path)
    monkeypatch.setattr(fda_identity_audit, "run_fda_identity_audit", _run)

    rc = qol_cli.dispatch(["ligands", "audit", "fda"])

    assert rc == 0
    assert len(captured) == 1
    config = captured[0]
    assert config.mapping_csv == (
        tmp_path / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv"
    )
    assert config.source_sdf is None
    assert config.approval_manifest is None
    assert config.library_dir is None
    assert config.output_dir == expected_out
    assert config.classified_nonmed_csv == (
        tmp_path / "docs" / "fda_non_medication_substances.csv"
    )
    output = capsys.readouterr().out
    assert "ligand_audit: fda" in output
    assert "blocker_count: 0" in output
    assert "output_dir: outputs/data/fda_identity_audit" in output
    assert (
        "summary_json: outputs/data/fda_identity_audit/fda_identity_summary.json"
        in output
    )


@pytest.mark.parametrize("failure_flag", ["--fail-on-blockers", "--fail-on-blocking"])
def test_fda_audit_forwards_explicit_paths_and_can_fail_on_blockers(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    failure_flag: str,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    approval_manifest = tmp_path / "fda_filter_manifest.json"
    library = tmp_path / "library"
    output = tmp_path / "audit"
    classified = tmp_path / "classified.csv"
    captured: list[fda_identity_audit.FDAIdentityAuditConfig] = []

    def _run(config: fda_identity_audit.FDAIdentityAuditConfig) -> SimpleNamespace:
        captured.append(config)
        return _audit_result(config.output_dir, blockers=2)

    monkeypatch.setattr(ligand_cli._bindings, "repo_root", lambda: tmp_path)
    monkeypatch.setattr(fda_identity_audit, "run_fda_identity_audit", _run)

    rc = ligand_cli._cmd_ligands(
        [
            "audit",
            "fda",
            "--mapping-csv",
            str(mapping),
            "--source-sdf",
            str(source),
            "--approval-manifest",
            str(approval_manifest),
            "--library-dir",
            str(library),
            "--out-dir",
            str(output),
            "--classified-non-medication-csv",
            str(classified),
            failure_flag,
        ]
    )

    assert rc == 1
    config = captured[0]
    assert config.mapping_csv == mapping
    assert config.source_sdf == source
    assert config.approval_manifest == approval_manifest
    assert config.library_dir == library
    assert config.output_dir == output
    assert config.classified_nonmed_csv == classified


def test_fda_audit_input_error_returns_two(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(ligand_cli._bindings, "repo_root", lambda: tmp_path)

    def _raise(_config: fda_identity_audit.FDAIdentityAuditConfig) -> None:
        raise ValueError("mapping has no ligand identifier column")

    monkeypatch.setattr(fda_identity_audit, "run_fda_identity_audit", _raise)

    rc = ligand_cli._cmd_ligands(["audit", "fda"])

    assert rc == 2
    assert "mapping has no ligand identifier column" in capsys.readouterr().err


def test_fda_audit_help_describes_offline_manual_inputs(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as exc_info:
        ligand_cli._cmd_ligands(["audit", "fda", "--help"])

    assert exc_info.value.code == 0
    output = capsys.readouterr().out
    assert "offline and never installs a library" in output
    for option in (
        "--mapping-csv",
        "--source-sdf",
        "--approval-manifest",
        "--library-dir",
        "--out-dir",
        "--classified-non-medication-csv",
        "--fail-on-blockers",
    ):
        assert option in output
    assert "--fail-on-blocking" not in output
