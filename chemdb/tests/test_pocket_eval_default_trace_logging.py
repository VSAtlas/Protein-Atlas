from __future__ import annotations

import json
import logging
from pathlib import Path

import pocket_eval


def test_pocket_eval_default_trace_logging(tmp_path: Path, monkeypatch, caplog) -> None:
    pdb_id = "TEST"
    run_id = "trace_run"
    docked_root = tmp_path / "docked"
    prepped_root = tmp_path / "prepped_ligands"

    cfg = {
        "RUN_ID": run_id,
        "DOCKED_DIR": str(docked_root),
        "POCKET_EVAL_FOLDS": 2,
        "POCKET_EVAL_SEED": 0,
    }

    class DummyPaths:
        def __init__(self, docked: Path, prepped: Path) -> None:
            self._docked = docked
            self.prepped_ligands_dir = prepped

        def docked_pdb_root(self) -> Path:
            return self._docked

    def fake_make_paths(cfg_obj, base_id, pdb_file):
        return DummyPaths(
            docked_root / run_id / base_id,
            prepped_root / base_id,
        )

    monkeypatch.setattr(pocket_eval, "make_paths", fake_make_paths)

    dock_root = docked_root / run_id / pdb_id / "pocket_eval"
    dock_root.mkdir(parents=True, exist_ok=True)
    cache_payload = {
        "pdb_id": pdb_id,
        "rows": [
            {"ligand_id": "lig_strong", "smiles": "C", "label": "strong"},
            {"ligand_id": "lig_non", "smiles": "CC", "label": "non"},
        ],
    }
    (dock_root / "calibrator_cache.json").write_text(
        json.dumps(cache_payload), encoding="utf-8"
    )

    pockets_json_path = tmp_path / "pockets.json"
    pockets_payload = {
        "pockets": [
            {
                "pocket_id": "pocket_1",
                "center": [0.0, 0.0, 0.0],
                "bounds": {"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
            }
        ]
    }
    pockets_json_path.write_text(json.dumps(pockets_payload), encoding="utf-8")

    receptor_pdbqt = tmp_path / "receptor.pdbqt"
    receptor_pdbqt.write_text("RECEPTOR\n", encoding="utf-8")

    def fake_prepare_calibrator_ligands(
        cal_rows, out_dir, cfg_obj, logger=None, inputs_dir_override=None
    ):
        return [
            {
                "ligand_id": "lig_strong",
                "label": "strong",
                "pdbqt_path": str(tmp_path / "lig_strong.pdbqt"),
            },
            {
                "ligand_id": "lig_non",
                "label": "non",
                "pdbqt_path": str(tmp_path / "lig_non.pdbqt"),
            },
        ]

    monkeypatch.setattr(
        pocket_eval, "prepare_calibrator_ligands", fake_prepare_calibrator_ligands
    )

    def fake_run_docking_task(vina_exe, config_path, ligand_id, out_path):
        score = -7.5 if "strong" in ligand_id else -5.0
        return out_path, score

    monkeypatch.setattr(pocket_eval, "run_docking_task", fake_run_docking_task)

    logger = logging.getLogger("pocket-eval-test")
    caplog.set_level(logging.INFO, logger=logger.name)

    pocket_eval.select_pocket_with_eval(
        pdb_id=pdb_id,
        cfg=cfg,
        logger=logger,
        pockets_json_path=pockets_json_path,
        receptor_pdbqt_path=receptor_pdbqt,
    )

    assert "[pocket-eval] calibrator_cache loaded" in caplog.text
    assert "[pocket-eval] pockets_json loaded" in caplog.text
    assert "[pocket-eval] dock.start" in caplog.text
    assert "[pocket-eval] dock.done" in caplog.text
