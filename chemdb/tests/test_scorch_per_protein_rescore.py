from __future__ import annotations

import csv
import logging
import os
import subprocess
from pathlib import Path
from types import SimpleNamespace

import main
import cli.postrun_hooks_runtime as postrun_hooks
import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def _touch(path: Path, content: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def test_scorch_runs_per_protein_after_docking(monkeypatch, tmp_path: Path) -> None:
    run_id = "scorch_per_protein"
    input_dir = tmp_path / "input_pdbs"
    docked_root = tmp_path / "docked"
    post_docked_root = tmp_path / "post_docked"
    output_dir = tmp_path / "processed_pdbs"
    logs_dir = tmp_path / "run_logs"
    config_run_dir = tmp_path / "configs" / run_id
    overall_dir = tmp_path

    _touch(input_dir / "test.pdb", "HEADER TEST\n")
    _touch(input_dir / "t3mp.pdb", "HEADER T3MP\n")

    cfg = {
        "RUN_ID": run_id,
        "OVERALL_DIR": str(overall_dir),
        "INPUT_DIR": str(input_dir),
        "DOCKED_DIR": str(docked_root),
        "OUTPUT_DIR": str(output_dir),
        "OUTPUT_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
        "CONFIGS_DIR": str(tmp_path / "configs"),
        "DOCKING_MODE": "discovery",
        "SPECIFIED_PROTEINS": "test,t3mp",
        "APO_HOLO_MODE": "APO",
        "PH_ENSEMBLE": False,
        "TEST_MODE_ENABLE": "off",
        "USE_DOCK6": True,
        "USE_GNINA": True,
        "USE_LEDOCK": True,
        "USE_SCORCH": True,
        "CPU": 1,
        "MAX_PARALLEL_JOBS": 1,
        "SCORCH_THREADS": 1,
        "MMGBSA_ENABLED": False,
    }

    event_log: list[tuple[str, str]] = []
    scorch_calls: list[tuple[str, str, str, str]] = []

    def _noop(*_args, **_kwargs) -> None:
        return None

    def fake_init_config_run_dir(cfg_obj, run_id, reset, logger) -> None:
        del run_id, reset, logger
        config_run_dir.mkdir(parents=True, exist_ok=True)
        cfg_obj["CONFIG_RUN_DIR"] = str(config_run_dir)

    monkeypatch.setattr(main, "load_inputs", lambda: cfg.copy())
    monkeypatch.setattr(main, "validate_config", lambda _cfg: None)
    monkeypatch.setattr(
        main,
        "get_atom_rules",
        lambda: SimpleNamespace(
            alias_sets=SimpleNamespace(waters=set()),
            waters=set(),
            cofactors=set(),
            elem_tokens_canonical=set(),
            policy_mode="LEGACY",
            ),
        raising=False,
    )
    monkeypatch.setattr(main, "init_config_run_dir", fake_init_config_run_dir)
    monkeypatch.setattr(main, "run_logs_dir", lambda _cfg: logs_dir, raising=False)
    monkeypatch.setattr(
        main,
        "_prepare_run_logfile",
        lambda _run_id: str(tmp_path / "logs" / f"main_{_run_id}.log"),
        raising=False,
    )
    monkeypatch.setattr(main, "_tee_stdio_to", _noop, raising=False)
    monkeypatch.setattr(main, "bootstrap_root_logging", _noop, raising=False)
    monkeypatch.setattr(main, "ensure_file_handler", _noop, raising=False)
    monkeypatch.setattr(main, "init_run_manifest", _noop, raising=False)
    monkeypatch.setattr(main, "update_manifest_for_config_hash", _noop, raising=False)
    monkeypatch.setattr(
        main, "update_manifest_for_scheduled_proteins", _noop, raising=False
    )
    monkeypatch.setattr(main, "update_manifest_for_run_config", _noop, raising=False)
    monkeypatch.setattr(main, "update_manifest_for_protein_start", _noop, raising=False)
    monkeypatch.setattr(main, "update_manifest_for_protein_success", _noop, raising=False)
    monkeypatch.setattr(main, "update_manifest_for_protein_failure", _noop, raising=False)
    monkeypatch.setattr(main, "apply_pocket_detection_events", _noop, raising=False)
    monkeypatch.setattr(main, "finalize_run_manifest", _noop, raising=False)
    monkeypatch.setattr(
        main, "define_docking_stages", lambda _mode: ["stage1"], raising=False
    )
    monkeypatch.setattr(main, "resolve_apo_holo_mode", lambda _cfg: ("apo", ["APO"]))
    monkeypatch.setattr(main, "_maybe_run_mmgbsa_for_pdb", _noop, raising=False)
    monkeypatch.setattr(main, "_maybe_run_dud_eval", _noop, raising=False)
    monkeypatch.setattr(main, "_maybe_run_master_schema_export", _noop, raising=False)
    monkeypatch.setattr(main, "_maybe_run_report_generation", _noop, raising=False)
    monkeypatch.setattr(main, "_log_rescore_verification", _noop, raising=False)

    def fake_process_one_protein(cfg_for_pdb, pdb_file, stages, params) -> None:
        del stages, params
        pdb_id = Path(pdb_file).stem.upper()
        variant = os.environ.get("APO_HOLO_VARIANT", "APO") or "APO"
        ph = "ph_7_0"
        combo_root = Path(cfg_for_pdb["DOCKED_DIR"]) / run_id / pdb_id / variant / ph

        _touch(
            combo_root / "ledock_stage1" / f"{pdb_id.lower()}_ledock_pose.dok",
            "REMARK ledock\n",
        )
        _touch(
            combo_root / "dock6_stage1" / f"{pdb_id.lower()}_dock6_pose.mol2",
            "@<TRIPOS>MOLECULE\nPOSE\n",
        )
        _touch(
            combo_root / "gnina_stage1" / f"{pdb_id.lower()}_gnina_pose.sdf",
            "$$$$\n",
        )
        _write_csv(
            combo_root / "consensus_docking_scores.csv",
            ["ligand", "consensus_score"],
            [[f"{pdb_id.lower()}_pose_stage1.pdbqt", "-7.0"]],
        )
        _touch(Path(cfg_for_pdb["DOCKED_DIR"]) / run_id / pdb_id / "DOCKING_DONE", "ok\n")
        event_log.append(("dock_done", pdb_id))

    monkeypatch.setattr(main, "process_one_protein", fake_process_one_protein)

    def fake_postrun_subprocess(cmd, cwd=None, check=False, **kwargs):
        del cwd, check, kwargs
        cmd_list = [str(part) for part in cmd]
        if not any(part.endswith("rescoring_scorch.py") for part in cmd_list):
            return subprocess.CompletedProcess(cmd_list, 0)
        assert "--pdb-id" in cmd_list, f"Expected per-protein scope, got: {cmd_list}"
        call_run_id = cmd_list[cmd_list.index("--run-id") + 1]
        call_pdb = cmd_list[cmd_list.index("--pdb-id") + 1].upper()
        call_variant = (
            cmd_list[cmd_list.index("--variant") + 1] if "--variant" in cmd_list else "APO"
        )
        call_ph = cmd_list[cmd_list.index("--ph") + 1] if "--ph" in cmd_list else "ph_7_0"
        scorch_calls.append((call_run_id, call_pdb, call_variant, call_ph))

        done_marker = docked_root / call_run_id / call_pdb / "DOCKING_DONE"
        assert done_marker.exists(), f"SCORCH started before docking for {call_pdb}"
        event_log.append(("scorch_start", call_pdb))

        post_combo = post_docked_root / call_run_id / call_pdb / call_variant / call_ph
        _touch(
            post_combo / "ledock_pdbqt" / f"{call_pdb.lower()}__ledock_stage1.pdbqt",
            "REMARK ledock pdbqt\n",
        )
        _touch(
            post_combo / "dock6_pdbqt" / f"{call_pdb.lower()}__dock6_stage1.pdbqt",
            "REMARK dock6 pdbqt\n",
        )
        _touch(
            post_combo / "gnina_best" / f"{call_pdb.lower()}_gnina_stage1.pdbqt",
            "REMARK gnina pdbqt\n",
        )
        _write_csv(
            post_combo / "scorch_scores_ledock.csv",
            ["ligand", "SCORCH_score"],
            [[f"{call_pdb.lower()}__ledock_stage1.pdbqt", "-5.0"]],
        )
        _write_csv(
            post_combo / "scorch_scores_dock6.csv",
            ["ligand", "SCORCH_score"],
            [[f"{call_pdb.lower()}__dock6_stage1.pdbqt", "-5.1"]],
        )
        _write_csv(
            post_combo / "scorch_scores_gnina_best.csv",
            ["ligand", "SCORCH_score"],
            [[f"{call_pdb.lower()}_gnina_stage1.pdbqt", "-4.9"]],
        )
        _write_csv(
            post_combo / "scorch_scores_all.csv",
            ["ligand", "SCORCH_score"],
            [[f"{call_pdb.lower()}__ledock_stage1.pdbqt", "-5.0"]],
        )
        _touch(post_combo / "scorch" / "_DONE", "ok\n")

        event_log.append(("scorch_done", call_pdb))
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", fake_postrun_subprocess)

    argv = ["main.py", "--run-id", run_id, "--pdb", "test", "--pdb", "t3mp", "-fast"]
    monkeypatch.setattr(main.sys, "argv", argv)

    main.main()

    assert len(scorch_calls) == 2
    assert {entry[1] for entry in scorch_calls} == {"TEST", "T3MP"}

    for pdb_id in ("TEST", "T3MP"):
        combo_root = post_docked_root / run_id / pdb_id / "APO" / "ph_7_0"
        assert (combo_root / "scorch" / "_DONE").exists()
        assert (combo_root / "ledock_pdbqt" / f"{pdb_id.lower()}__ledock_stage1.pdbqt").exists()
        assert (combo_root / "dock6_pdbqt" / f"{pdb_id.lower()}__dock6_stage1.pdbqt").exists()
        assert (combo_root / "gnina_best" / f"{pdb_id.lower()}_gnina_stage1.pdbqt").exists()
        assert (combo_root / "scorch_scores_ledock.csv").exists()
        assert (combo_root / "scorch_scores_dock6.csv").exists()
        assert (combo_root / "scorch_scores_gnina_best.csv").exists()
        assert (combo_root / "scorch_scores_all.csv").exists()

        dock_idx = event_log.index(("dock_done", pdb_id))
        scorch_start_idx = event_log.index(("scorch_start", pdb_id))
        scorch_done_idx = event_log.index(("scorch_done", pdb_id))
        assert dock_idx < scorch_start_idx < scorch_done_idx


def test_scorch_done_sentinel_scopes_all_ph_when_unfiltered(tmp_path: Path) -> None:
    post_run_root = tmp_path / "post_docked" / "run_x"
    combos = {
        ("TEST", "APO", "ph_6_5"),
        ("TEST", "APO", "ph_7_0"),
    }

    _touch(post_run_root / "TEST" / "APO" / "ph_6_5" / "scorch" / "_DONE", "ok\n")
    pending = rescoring_scorch._filter_done_combos(
        combos,
        post_run_root,
        overwrite=False,
        logger=logging.getLogger("scorch_done_filter"),
    )
    assert pending == {("TEST", "APO", "ph_7_0")}

    _touch(post_run_root / "TEST" / "APO" / "ph_7_0" / "scorch" / "_DONE", "ok\n")
    pending_after_all_done = rescoring_scorch._filter_done_combos(
        combos,
        post_run_root,
        overwrite=False,
        logger=logging.getLogger("scorch_done_filter"),
    )
    assert pending_after_all_done == set()
