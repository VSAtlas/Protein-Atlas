from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import cli.postrun_hooks_runtime as postrun_hooks

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_module(module: str, args: list[str]) -> None:
    subprocess.run(
        [sys.executable, "-m", module, *args],
        cwd=str(REPO_ROOT),
        check=True,
    )


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_dud_eval_cli_relocated_roots_writes_summary(tmp_path: Path) -> None:
    relocated_root = tmp_path / "relocated_repo"
    run_id = "RELOC_DUD"

    _write_csv(
        relocated_root / "docked" / run_id / "TEST" / "docking_score_long.csv",
        ["ligand_file", "score"],
        [
            {"ligand_file": "actives_final_00001.pdbqt", "score": "-9.0"},
            {"ligand_file": "decoys_final_00001.pdbqt", "score": "-4.0"},
        ],
    )
    _write_csv(
        relocated_root
        / "post_docked"
        / run_id
        / "TEST"
        / "HOLO"
        / "pH7_0"
        / "consensus_reranked_scorch.csv",
        ["ligand", "final_score", "run_id", "pdb_id", "variant", "ph_label"],
        [
            {
                "ligand": "actives_final_00001.pdbqt",
                "final_score": "0.9",
                "run_id": run_id,
                "pdb_id": "TEST",
                "variant": "HOLO",
                "ph_label": "pH7_0",
            },
            {
                "ligand": "decoys_final_00001.pdbqt",
                "final_score": "0.1",
                "run_id": run_id,
                "pdb_id": "TEST",
                "variant": "HOLO",
                "ph_label": "pH7_0",
            },
        ],
    )

    out_root = relocated_root / "analysis" / "dud_eval"
    _run_module(
        "analysis.cli.dud_eval",
        [
            "--run-id",
            run_id,
            "--pdb-id",
            "TEST",
            "--docked-root",
            str(relocated_root / "docked"),
            "--post-docked-root",
            str(relocated_root / "post_docked"),
            "--out-dir",
            str(out_root),
            "--log-level",
            "INFO",
        ],
    )

    run_out = out_root / run_id
    summary_candidates = list(run_out.glob("summary*.tsv"))
    assert summary_candidates
    assert (run_out / "TEST" / "metrics.tsv").exists()


def test_run_report_cli_relocated_repo_root_writes_report_html(tmp_path: Path) -> None:
    relocated_root = tmp_path / "relocated_repo"
    run_id = "RELOC_REPORT"

    _write_csv(
        relocated_root / "data" / run_id / "master_rows.csv",
        [
            "pdb_id",
            "variant",
            "ph_label",
            "ligand_base",
            "z_selected",
            "z_selected_source",
            "is_decoy",
            "is_control",
            "ef1",
            "consensus_mu_decoy",
            "consensus_sigma_decoy",
            "consensus_n_decoys",
            "blend_mu_decoy",
            "blend_sigma_decoy",
            "blend_n_decoys",
        ],
        [
            {
                "pdb_id": "TEST",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "ligand_base": "LIG1",
                "z_selected": "2.0",
                "z_selected_source": "stage2",
                "is_decoy": "0",
                "is_control": "0",
                "ef1": "5.0",
                "consensus_mu_decoy": "0.5",
                "consensus_sigma_decoy": "0.2",
                "consensus_n_decoys": "10",
                "blend_mu_decoy": "0.4",
                "blend_sigma_decoy": "0.2",
                "blend_n_decoys": "10",
            },
            {
                "pdb_id": "TEST",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "ligand_base": "DECOY1",
                "z_selected": "-1.0",
                "z_selected_source": "stage1",
                "is_decoy": "1",
                "is_control": "0",
                "ef1": "5.0",
                "consensus_mu_decoy": "0.5",
                "consensus_sigma_decoy": "0.2",
                "consensus_n_decoys": "10",
                "blend_mu_decoy": "0.4",
                "blend_sigma_decoy": "0.2",
                "blend_n_decoys": "10",
            },
        ],
    )

    _run_module(
        "analysis.cli.run_report",
        [
            "--run-id",
            run_id,
            "--repo-root",
            str(relocated_root),
            "--overwrite",
            "--emit-html",
        ],
    )

    assert (relocated_root / "data" / run_id / "report.html").exists()


def test_master_schema_export_cli_relocated_repo_root_writes_master_rows(
    tmp_path: Path,
) -> None:
    relocated_root = tmp_path / "relocated_repo"
    run_id = "RELOC_MASTER"

    _write_csv(
        relocated_root
        / "post_docked"
        / run_id
        / "TEST"
        / "HOLO"
        / "pH7_0"
        / "consensus_reranked_scorch.csv",
        [
            "pdb_id",
            "variant",
            "ph_label",
            "ligand",
            "ligand_file",
            "library",
            "z_vs_decoys_consensus",
            "z_vs_decoys_blend",
            "z_vs_decoys_consensus_pre",
            "run_mode",
        ],
        [
            {
                "pdb_id": "TEST",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "ligand": "LIG_A.sanitized.pdbqt",
                "ligand_file": "LIG_A.sanitized.pdbqt",
                "library": "LIB1",
                "z_vs_decoys_consensus": "1.0",
                "z_vs_decoys_blend": "1.2",
                "z_vs_decoys_consensus_pre": "0.8",
                "run_mode": "fda",
            },
            {
                "pdb_id": "TEST",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "ligand": "decoys_LIG_A.pdbqt",
                "ligand_file": "decoys_LIG_A.pdbqt",
                "library": "LIB1",
                "z_vs_decoys_consensus": "-1.0",
                "z_vs_decoys_blend": "-1.2",
                "z_vs_decoys_consensus_pre": "-0.8",
                "run_mode": "dud",
            },
        ],
    )
    (relocated_root / "processed_pdbs" / "TEST" / "ligands_raw").mkdir(
        parents=True, exist_ok=True
    )

    _run_module(
        "analysis.cli.master_schema_export",
        [
            "--run-id",
            run_id,
            "--repo-root",
            str(relocated_root),
            "--overwrite",
        ],
    )

    output_csv = relocated_root / "outputs" / "data" / run_id / "master_rows.csv"
    assert output_csv.exists()
    with output_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows


def test_postrun_analysis_hooks_use_relocated_roots(monkeypatch, tmp_path: Path) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False, **kwargs):
        del cwd, check, kwargs
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)

    relocated_root = tmp_path / "relocated_repo"
    cfg = {
        "OVERALL_DIR": str(relocated_root),
        "DOCKED_DIR": str(relocated_root / "docked"),
        "POST_DOCKED_DIR": str(relocated_root / "post_docked"),
        "TEST_MODE_ENABLE": "dud",
    }
    run_id = "RELOC_HOOK"

    postrun_hooks._maybe_run_dud_eval(cfg, run_id, ["TEST.pdb"])
    postrun_hooks._maybe_run_report_generation(cfg, run_id)
    postrun_hooks._maybe_run_master_schema_export(cfg, run_id)

    assert len(captured) == 3
    dud_cmd, report_cmd, master_cmd = captured

    assert dud_cmd[dud_cmd.index("--docked-root") + 1] == str(relocated_root / "docked")
    assert dud_cmd[dud_cmd.index("--post-docked-root") + 1] == str(
        relocated_root / "post_docked"
    )
    assert report_cmd[report_cmd.index("--repo-root") + 1] == str(relocated_root)
    assert master_cmd[master_cmd.index("--repo-root") + 1] == str(relocated_root)


def test_postrun_analysis_hooks_resolve_relative_cfg_roots(
    monkeypatch, tmp_path: Path
) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False, **kwargs):
        del cwd, check, kwargs
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)

    relocated_root = tmp_path / "relocated_repo"
    cfg = {
        "OVERALL_DIR": str(relocated_root),
        "DOCKED_DIR": "docked",
        "POST_DOCKED_DIR": "post_docked",
        "TEST_MODE_ENABLE": "dud",
    }
    run_id = "RELOC_HOOK_REL"

    postrun_hooks._maybe_run_dud_eval(cfg, run_id, ["TEST.pdb"])
    postrun_hooks._maybe_run_report_generation(cfg, run_id)
    postrun_hooks._maybe_run_master_schema_export(cfg, run_id)

    assert len(captured) == 3
    dud_cmd, report_cmd, master_cmd = captured

    assert dud_cmd[dud_cmd.index("--docked-root") + 1] == str(relocated_root / "docked")
    assert dud_cmd[dud_cmd.index("--post-docked-root") + 1] == str(
        relocated_root / "post_docked"
    )
    assert report_cmd[report_cmd.index("--repo-root") + 1] == str(relocated_root)
    assert master_cmd[master_cmd.index("--repo-root") + 1] == str(relocated_root)


def test_postrun_analysis_hooks_infer_relocated_outputs_root(
    monkeypatch, tmp_path: Path
) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False, **kwargs):
        del cwd, check, kwargs
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)

    work_root = tmp_path / "work_repo"
    scratch_root = tmp_path / "scratch_repo"
    cfg = {
        "OVERALL_DIR": str(work_root),
        "DOCKED_DIR": str(scratch_root / "outputs" / "docked"),
        "POST_DOCKED_DIR": str(scratch_root / "outputs" / "post_docked"),
        "TEST_MODE_ENABLE": "dud",
    }
    run_id = "RELOC_HOOK_OUTPUTS"

    postrun_hooks._maybe_run_dud_eval(cfg, run_id, ["TEST.pdb"])
    postrun_hooks._maybe_run_report_generation(cfg, run_id)
    postrun_hooks._maybe_run_master_schema_export(cfg, run_id)
    assert postrun_hooks._maybe_run_throughput_integrity(cfg, run_id, strict=False)

    assert len(captured) == 4
    dud_cmd, report_cmd, master_cmd, throughput_cmd = captured

    assert dud_cmd[dud_cmd.index("--docked-root") + 1] == str(
        scratch_root / "outputs" / "docked"
    )
    assert dud_cmd[dud_cmd.index("--post-docked-root") + 1] == str(
        scratch_root / "outputs" / "post_docked"
    )
    assert report_cmd[report_cmd.index("--repo-root") + 1] == str(scratch_root)
    assert master_cmd[master_cmd.index("--repo-root") + 1] == str(scratch_root)
    assert throughput_cmd[throughput_cmd.index("--repo-root") + 1] == str(scratch_root)


def test_postrun_analysis_hooks_use_overall_dir_env_defaults(
    monkeypatch, tmp_path: Path
) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False, **kwargs):
        del cwd, check, kwargs
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)

    relocated_root = tmp_path / "relocated_repo"
    monkeypatch.setenv("OVERALL_DIR", str(relocated_root))
    monkeypatch.delenv("DOCKED_DIR", raising=False)
    monkeypatch.delenv("POST_DOCKED_DIR", raising=False)

    cfg = {
        "TEST_MODE_ENABLE": "dud",
    }
    run_id = "RELOC_HOOK_ENV"

    postrun_hooks._maybe_run_dud_eval(cfg, run_id, ["TEST.pdb"])
    postrun_hooks._maybe_run_report_generation(cfg, run_id)
    postrun_hooks._maybe_run_master_schema_export(cfg, run_id)

    assert len(captured) == 3
    dud_cmd, report_cmd, master_cmd = captured

    assert dud_cmd[dud_cmd.index("--docked-root") + 1] == str(relocated_root / "docked")
    assert dud_cmd[dud_cmd.index("--post-docked-root") + 1] == str(
        relocated_root / "post_docked"
    )
    assert report_cmd[report_cmd.index("--repo-root") + 1] == str(relocated_root)
    assert master_cmd[master_cmd.index("--repo-root") + 1] == str(relocated_root)
