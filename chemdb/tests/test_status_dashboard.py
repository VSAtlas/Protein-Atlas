from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from analysis.cli.summarize_run_errors import build_error_summary, classify_root_cause
from cli import qol_cli
from cli.status_dashboard import (
    ETA_HISTORY_NAME,
    build_runs_index,
    build_status_dashboard,
    default_status_html_path,
    explain_status,
    format_runs_table,
    format_status_dashboard,
    parse_sacct_output,
    parse_squeue_output,
    render_status_html,
    write_status_html,
)


FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "slurm"


def _write_manifest(root: Path, run_id: str, payload: dict) -> Path:
    manifest_dir = root / "manifests" / run_id
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_dir / "run_manifest.yaml"
    manifest_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return manifest_path


def _manifest(
    *,
    status: str,
    scheduled: int,
    completed: int,
    failed: int = 0,
    proteins: dict | None = None,
    wall_time_sec: float | None = None,
    slurm_job_id: str | None = None,
) -> dict:
    payload: dict[str, Any] = {
        "status": status,
        "summary": {
            "total_proteins_scheduled": scheduled,
            "total_proteins_completed": completed,
            "total_proteins_failed": failed,
        },
        "timing": {
            "started_at": "2026-05-02T00:00:00Z",
        },
        "proteins": proteins or {},
    }
    if wall_time_sec is not None:
        payload["timing"]["wall_time_sec"] = wall_time_sec
    if slurm_job_id:
        payload["resources"] = {"slurm": {"job_id": slurm_job_id}}
    return payload


def _install_fake_slurm_bin(
    tmp_path: Path,
    monkeypatch,
    *,
    squeue_fixture: str | None = None,
    sacct_fixture: str | None = None,
    squeue_empty: bool = False,
) -> Path:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    if squeue_fixture:
        squeue_path = FIXTURE_ROOT / squeue_fixture
        fake_squeue = fake_bin / "squeue"
        fake_squeue.write_text(
            "#!/bin/sh\n" f"cat {squeue_path}\n",
            encoding="utf-8",
        )
        fake_squeue.chmod(0o755)
    elif squeue_empty:
        fake_squeue = fake_bin / "squeue"
        fake_squeue.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        fake_squeue.chmod(0o755)
    if sacct_fixture:
        sacct_path = FIXTURE_ROOT / sacct_fixture
        fake_sacct = fake_bin / "sacct"
        fake_sacct.write_text(
            "#!/bin/sh\n" f"cat {sacct_path}\n",
            encoding="utf-8",
        )
        fake_sacct.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}")
    return fake_bin


def _eta_history_cache_path(root: Path) -> Path:
    return root / "manifests" / ETA_HISTORY_NAME


def _read_eta_history_cache(root: Path) -> dict[str, Any]:
    return json.loads(_eta_history_cache_path(root).read_text(encoding="utf-8"))


def test_slurm_fixture_parsers_cover_squeue_and_sacct_samples() -> None:
    squeue_jobs = parse_squeue_output(
        (FIXTURE_ROOT / "squeue_pipe_sample.txt").read_text(encoding="utf-8")
    )
    assert [job["state"] for job in squeue_jobs] == ["RUNNING", "PENDING", "PENDING"]
    assert squeue_jobs[1]["reason_or_nodelist"] == "(Resources)"

    sacct_jobs = parse_sacct_output(
        (FIXTURE_ROOT / "sacct_pipe_sample.txt").read_text(encoding="utf-8")
    )
    assert [job["state"] for job in sacct_jobs] == [
        "COMPLETED",
        "COMPLETED",
        "FAILED",
        "OUT_OF_MEMORY",
    ]
    assert sacct_jobs[-1]["elapsed_sec"] == 222.0
    assert sacct_jobs[-1]["exit_code"] == "0:125"


def test_status_dashboard_uses_fake_squeue_and_explains_pending_reason(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = "run_slurm"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(
            status="running",
            scheduled=2,
            completed=0,
            proteins={
                "TEST": {
                    "status": "running",
                    "stages": {"prep": {"status": "running"}},
                }
            },
            slurm_job_id="12345",
        ),
    )
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_squeue = fake_bin / "squeue"
    fake_squeue.write_text(
        "#!/bin/sh\n"
        f"cat {str(FIXTURE_ROOT / 'squeue_pipe_sample.txt')}\n",
        encoding="utf-8",
    )
    fake_squeue.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}")

    status = build_status_dashboard(tmp_path, run_id, include_slurm=True)
    assert status["slurm"]["source"] == "squeue"
    assert status["slurm"]["state_counts"] == {"PENDING": 2, "RUNNING": 1}
    assert any("Resources" in item for item in status["explanations"])

    rendered = format_status_dashboard(status, explain=True)
    assert "slurm: PENDING=2, RUNNING=1" in rendered
    assert "explain:" in rendered


def test_eta_history_calibrates_running_run_from_completed_manifest(tmp_path: Path) -> None:
    _write_manifest(
        tmp_path,
        "completed_run",
        _manifest(
            status="completed",
            scheduled=2,
            completed=2,
            wall_time_sec=100.0,
            proteins={
                "A": {
                    "status": "completed",
                    "stages": {"prep": {"status": "completed", "timing": {"wall_time_sec": 5}}},
                },
                "B": {
                    "status": "completed",
                    "stages": {"prep": {"status": "completed", "timing": {"wall_time_sec": 7}}},
                },
            },
        ),
    )
    _write_manifest(
        tmp_path,
        "running_run",
        _manifest(
            status="running",
            scheduled=4,
            completed=1,
            proteins={
                "A": {"status": "completed"},
                "B": {"status": "running"},
                "C": {"status": "pending"},
                "D": {"status": "pending"},
            },
        ),
    )

    status = build_status_dashboard(tmp_path, "running_run", include_slurm=False)
    assert status["eta"]["basis"] == "history_per_protein"
    assert status["eta"]["eta_sec"] == 150.0
    assert status["eta_history"]["per_protein_wall_sec_median"] == 50.0
    assert _eta_history_cache_path(tmp_path).exists()


def test_eta_history_cache_skips_running_and_incomplete_completed_runs(tmp_path: Path) -> None:
    _write_manifest(
        tmp_path,
        "completed_ok",
        _manifest(status="completed", scheduled=2, completed=2, wall_time_sec=80.0),
    )
    _write_manifest(
        tmp_path,
        "running_now",
        _manifest(status="running", scheduled=3, completed=1),
    )
    _write_manifest(
        tmp_path,
        "completed_no_wall",
        _manifest(status="completed", scheduled=1, completed=1),
    )

    build_status_dashboard(tmp_path, "running_now", include_slurm=False)

    cache = _read_eta_history_cache(tmp_path)
    assert cache["schema_version"] == 1
    assert cache["run_count"] == 1
    assert {row["run_id"] for row in cache["runs"]} == {"completed_ok"}
    assert cache["per_protein_wall_sec"]["median"] == 40.0
    assert cache["per_protein_wall_sec"]["sample_count"] == 1


def test_eta_history_cache_refreshes_when_new_completed_run_arrives(tmp_path: Path) -> None:
    _write_manifest(
        tmp_path,
        "hist_one",
        _manifest(status="completed", scheduled=2, completed=2, wall_time_sec=40.0),
    )
    build_status_dashboard(tmp_path, "hist_one", include_slurm=False)
    first = _read_eta_history_cache(tmp_path)
    assert first["run_count"] == 1

    _write_manifest(
        tmp_path,
        "hist_two",
        _manifest(status="completed", scheduled=4, completed=4, wall_time_sec=120.0),
    )
    build_status_dashboard(tmp_path, "hist_two", include_slurm=False)
    second = _read_eta_history_cache(tmp_path)

    assert second["run_count"] == 2
    assert {row["run_id"] for row in second["runs"]} == {"hist_one", "hist_two"}
    assert second["per_protein_wall_sec"]["median"] == 25.0
    assert second["updated_at"] >= first["updated_at"]


def test_root_cause_rules_yaml_can_be_overridden_per_repo(tmp_path: Path) -> None:
    rules_path = tmp_path / "analysis" / "config" / "root_cause_rules.yaml"
    rules_path.parent.mkdir(parents=True)
    rules_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "root_cause": "custom_fixture_bucket",
                        "patterns": ["very specific sentinel"],
                    }
                ],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    run_id = "run_errors"
    docked_dir = tmp_path / "docked" / run_id
    docked_dir.mkdir(parents=True)
    (docked_dir / "main.log").write_text(
        "2026-05-02 00:00:00 ERROR very specific sentinel broke receptor prep\n",
        encoding="utf-8",
    )

    summary = build_error_summary(root=tmp_path, run_id=run_id, recent_hours=0)
    assert summary["root_cause_counts"] == {"custom_fixture_bucket": 1}
    assert classify_root_cause("very specific sentinel", (("x", ("sentinel",)),)) == "x"


def test_runs_index_and_schema_docs_are_available(tmp_path: Path) -> None:
    _write_manifest(
        tmp_path,
        "run_a",
        _manifest(status="completed", scheduled=1, completed=1, wall_time_sec=10.0),
    )
    _write_manifest(
        tmp_path,
        "run_b",
        _manifest(status="failed", scheduled=2, completed=1, failed=1),
    )
    rows = build_runs_index(tmp_path, limit=10)
    assert {row["run_id"] for row in rows} == {"run_a", "run_b"}
    table = format_runs_table(rows)
    assert "run_id" in table
    assert "failed" in table

    schema = json.loads(
        (Path("docs") / "schemas" / "status_dashboard.schema.json").read_text(
            encoding="utf-8"
        )
    )
    assert "explanations" in schema["required"]
    assert "eta_history" in schema["properties"]


def test_explain_status_flags_missing_report(tmp_path: Path) -> None:
    run_id = "needs_report"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(status="completed", scheduled=1, completed=1, wall_time_sec=5.0),
    )
    status = build_status_dashboard(tmp_path, run_id, include_slurm=False)
    explanations = explain_status(status)
    assert any("atlas report needs_report --vina" in item for item in explanations)


def test_status_dashboard_surfaces_scorch_postprocessing_failures(
    tmp_path: Path,
) -> None:
    run_id = "scorch_failed"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(
            status="completed",
            scheduled=1,
            completed=1,
            proteins={
                "BNJS|HOLO|pH7_0": {
                    "status": "completed",
                    "stages": {
                        "postprocessing": {
                            "status": "failed",
                            "error": "returncode=1",
                            "timing": {"wall_time_sec": 2.7},
                            "details": {
                                "scorch_returncode": 1,
                                "scorch_mode": "subprocess",
                                "scorch_phase": "final",
                                "scorch_device_effective": "cpu",
                            },
                        }
                    },
                }
            },
        ),
    )

    status = build_status_dashboard(tmp_path, run_id, include_slurm=False)

    assert status["scorch_failures"]["failed"] == 1
    assert status["scorch_failures"]["returncodes"] == {"1": 1}
    rendered = format_status_dashboard(status, explain=True)
    assert "scorch_failures: failed=1" in rendered
    assert "SCORCH postprocessing failed for 1 combo" in rendered
    assert "atlas debug scorch_failed --deep" in rendered


def test_status_dashboard_html_renderer_uses_status_json_contract(tmp_path: Path) -> None:
    run_id = "html_run"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(
            status="completed",
            scheduled=1,
            completed=1,
            wall_time_sec=5.0,
            proteins={
                "TEST": {
                    "status": "completed",
                    "stages": {
                        "prep": {"status": "completed", "timing": {"wall_time_sec": 5}},
                        "docking": {"status": "completed", "timing": {"wall_time_sec": 12}},
                    },
                }
            },
        ),
    )
    status = build_status_dashboard(tmp_path, run_id, include_slurm=False, include_errors=True)
    html_text = render_status_html(status)
    assert "Atlas Status html_run" in html_text
    assert "Status JSON" in html_text
    assert "&quot;run_id&quot;: &quot;html_run&quot;" in html_text
    assert "Stages" in html_text

    out_path = write_status_html(status, tmp_path / "data" / run_id / "status.html")
    assert out_path.exists()
    written = out_path.read_text(encoding="utf-8")
    assert "Stages" in written
    assert "Explain" in written


def test_cmd_status_explain_prints_operator_actions(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    run_id = "cli_explain_run"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(status="completed", scheduled=1, completed=1, wall_time_sec=5.0),
    )
    monkeypatch.setattr(qol_cli, "_repo_root", lambda: tmp_path)

    rc = qol_cli._cmd_status([run_id, "--no-slurm"])
    plain = capsys.readouterr().out
    assert rc == 0
    assert "explain:" not in plain

    rc = qol_cli._cmd_status([run_id, "--explain", "--no-slurm"])
    explained = capsys.readouterr().out
    assert rc == 0
    assert "explain:" in explained
    assert "atlas report cli_explain_run --vina" in explained


def test_cmd_status_explain_uses_fake_squeue_fixture(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    run_id = "cli_slurm_explain"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(
            status="running",
            scheduled=2,
            completed=0,
            slurm_job_id="12345",
            proteins={
                "TEST": {
                    "status": "running",
                    "stages": {"prep": {"status": "running"}},
                }
            },
        ),
    )
    _install_fake_slurm_bin(tmp_path, monkeypatch, squeue_fixture="squeue_pipe_sample.txt")
    monkeypatch.setattr(qol_cli, "_repo_root", lambda: tmp_path)

    rc = qol_cli._cmd_status([run_id, "--explain"])
    output = capsys.readouterr().out

    assert rc == 0
    assert "slurm: PENDING=2, RUNNING=1" in output
    assert "explain:" in output
    assert "Resources" in output


def test_cmd_status_explain_uses_sacct_fixture_when_squeue_empty(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    run_id = "cli_sacct_explain"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(
            status="failed",
            scheduled=3,
            completed=2,
            failed=1,
            slurm_job_id="12345",
        ),
    )
    _install_fake_slurm_bin(
        tmp_path,
        monkeypatch,
        squeue_empty=True,
        sacct_fixture="sacct_pipe_sample.txt",
    )
    monkeypatch.setattr(qol_cli, "_repo_root", lambda: tmp_path)

    rc = qol_cli._cmd_status([run_id, "--explain"])
    output = capsys.readouterr().out

    assert rc == 0
    assert "slurm: COMPLETED=2, FAILED=1, OUT_OF_MEMORY=1" in output
    assert "OUT_OF_MEMORY" in output
    assert "0:125" in output


def test_cmd_status_html_writes_default_and_custom_paths(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    run_id = "cli_html_run"
    _write_manifest(
        tmp_path,
        run_id,
        _manifest(status="completed", scheduled=1, completed=1, wall_time_sec=5.0),
    )
    monkeypatch.setattr(qol_cli, "_repo_root", lambda: tmp_path)

    rc = qol_cli._cmd_status([run_id, "--html", "--no-slurm"])
    output = capsys.readouterr().out
    default_path = default_status_html_path(tmp_path, run_id)

    assert rc == 0
    assert default_path.exists()
    assert f"html: {default_path}" in output
    assert "Atlas Status cli_html_run" in default_path.read_text(encoding="utf-8")

    custom_path = tmp_path / "reports" / f"{run_id}_status.html"
    rc = qol_cli._cmd_status([run_id, "--html", str(custom_path), "--no-slurm"])
    output = capsys.readouterr().out

    assert rc == 0
    assert custom_path.exists()
    assert f"html: {custom_path}" in output


def test_cmd_status_json_includes_eta_history_summary(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    _write_manifest(
        tmp_path,
        "eta_hist_done",
        _manifest(status="completed", scheduled=2, completed=2, wall_time_sec=60.0),
    )
    _write_manifest(
        tmp_path,
        "eta_hist_running",
        _manifest(status="running", scheduled=4, completed=1),
    )
    monkeypatch.setattr(qol_cli, "_repo_root", lambda: tmp_path)

    rc = qol_cli._cmd_status(["eta_hist_running", "--json", "--no-slurm"])
    output = capsys.readouterr().out

    assert rc == 0
    payload = json.loads(output)
    assert payload["eta"]["basis"] == "history_per_protein"
    assert payload["eta_history"]["run_count"] == 1
    assert payload["eta_history"]["per_protein_wall_sec_median"] == 30.0
    assert _eta_history_cache_path(tmp_path).exists()
