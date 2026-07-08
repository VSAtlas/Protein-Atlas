from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path

import docking.docking_stage_runner_support as stage_runner
import docking.fallback_recenter as fallback_recenter


def test_validate_first_valid_pose_uses_unique_temp_files(
    tmp_path: Path, monkeypatch
) -> None:
    receptor = tmp_path / "rec.pdbqt"
    ligand = tmp_path / "lig_stage3.pdbqt"
    receptor.write_text("REMARK receptor\n", encoding="utf-8")
    ligand.write_text(
        "MODEL        1\n"
        "ATOM      1  C   LIG A   1       1.000   2.000   3.000  0.00  0.00      C\n"
        "ENDMDL\n",
        encoding="utf-8",
    )

    seen_tmp_paths: list[str] = []

    def _fake_validate_pose_pdbqt(**kwargs):
        seen_tmp_paths.append(str(kwargs["ligand_pdbqt"]))
        return {"valid": False, "reason": "mock_invalid"}

    monkeypatch.setattr(
        fallback_recenter,
        "validate_pose_pdbqt",
        _fake_validate_pose_pdbqt,
    )

    def _one_call(_idx: int) -> None:
        fallback_recenter.validate_first_valid_pose(
            receptor_pdbqt=str(receptor),
            ligand_pdbqt=str(ligand),
            pocket_center=(0.0, 0.0, 0.0),
            surface_coords=[],
            max_models=1,
        )

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(_one_call, range(32)))

    assert len(seen_tmp_paths) == 32
    assert len(set(seen_tmp_paths)) == 32
    assert all(path.endswith(".tmp.pdbqt") for path in seen_tmp_paths)
    assert list(tmp_path.glob("*.tmp.pdbqt")) == []


def test_run_one_stage_marks_missing_pose_file_non_fatal(
    tmp_path: Path, monkeypatch
) -> None:
    receptor = tmp_path / "rec.pdbqt"
    receptor.write_text("REMARK receptor\n", encoding="utf-8")
    ligand = tmp_path / "lig_a.sdf"
    ligand.write_text("ligand\n", encoding="utf-8")

    cfg = {
        "CPU": 1,
        "THREADS_PER_VINA": 1,
        "RUN_ID": "",
        "VINA_EXE": "vina",
        "OVERALL_DIR": str(tmp_path),
        "INPUT_DIR": str(tmp_path / "input_pdbs"),
        "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "CONFIG_RUN_DIR": str(tmp_path / "configs" / "run"),
    }
    stage = {"name": "stage3", "exhaustiveness": 1, "num_modes": 1}
    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "CONFIG_RUN_DIR"):
        Path(str(cfg[key])).mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(stage_runner, "extract_surface_atoms", lambda **_kwargs: [])
    monkeypatch.setattr(
        stage_runner,
        "emit_vina_config",
        lambda *_args, **_kwargs: (
            str(Path(str(cfg["CONFIG_RUN_DIR"])) / "cfg.conf"),
            str(tmp_path / "missing_out.pdbqt"),
        ),
    )
    monkeypatch.setattr(
        stage_runner,
        "run_docking_task",
        lambda *_args, **_kwargs: ("lig_a", -7.2),
    )
    monkeypatch.setattr(
        stage_runner,
        "filter_and_rewrite_poses_by_rmsd",
        lambda *_args, **_kwargs: (1, 0),
    )

    def _raise_missing_pose(**_kwargs):
        raise FileNotFoundError("simulated transient missing tmp pose")

    monkeypatch.setattr(stage_runner, "validate_first_valid_pose", _raise_missing_pose)

    @contextmanager
    def _fake_acquire_global_cores(*_args, **_kwargs):
        yield 1

    monkeypatch.setattr(stage_runner, "acquire_global_cores", _fake_acquire_global_cores)

    scores, validated, distances, raw_docked, invalids = stage_runner.run_one_stage(
        cfg=cfg,
        pdb_id="TEST",
        receptor_pdbqt=str(receptor),
        center=(0.0, 0.0, 0.0),
        box_size=(20.0, 20.0, 20.0),
        stage=stage,
        ligands=[str(ligand)],
        logger=logging.getLogger("test.pose.transient"),
        retry_mgr=stage_runner.RetryManager(max_retries=1),
        control_lookup={},
        ph_label=None,
    )

    assert scores == {}
    assert validated == []
    assert distances == []
    assert len(raw_docked) == 1
    assert str(ligand) in raw_docked
    assert len(invalids) == 1
    (_, reason), = invalids.values()
    assert reason == "pose_file_missing_transient"
