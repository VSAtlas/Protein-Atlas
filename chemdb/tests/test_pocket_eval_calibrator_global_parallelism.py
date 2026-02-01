from __future__ import annotations

import logging
import threading
from pathlib import Path

import pocket_eval


def test_pocket_eval_calibrator_global_parallelism(tmp_path: Path, monkeypatch) -> None:
    cfg = {"CPU": 2, "THREADS_PER_VINA": 1}
    prepared_scoring = [
        {
            "ligand_id": "lig_1",
            "label": "strong",
            "pdbqt_path": str(tmp_path / "lig_1.pdbqt"),
        }
    ]
    receptor_pdbqt = tmp_path / "receptor.pdbqt"
    receptor_pdbqt.write_text("RECEPTOR\n", encoding="utf-8")

    barrier = threading.Barrier(2)
    lock = threading.Lock()
    inflight = 0
    max_inflight = 0

    def fake_run_docking_task(vina_exe, config_path, ligand_id, out_path):
        nonlocal inflight, max_inflight
        with lock:
            inflight += 1
            max_inflight = max(max_inflight, inflight)
        try:
            barrier.wait(timeout=2)
        except threading.BrokenBarrierError as exc:
            raise RuntimeError("Barrier timeout; docking not parallel") from exc
        finally:
            with lock:
                inflight -= 1
        score = -7.0 if "P1" in str(config_path) else -6.0
        return out_path, score

    def fake_write_vina_config(*_args, **_kwargs) -> None:
        return None

    monkeypatch.setattr(pocket_eval, "run_docking_task", fake_run_docking_task)
    monkeypatch.setattr(pocket_eval, "_write_vina_config", fake_write_vina_config)

    pocket1_dir = tmp_path / "P1"
    pocket2_dir = tmp_path / "P2"
    pockets_plan = [
        {
            "pocket_id": "P1",
            "dock_dir": pocket1_dir,
            "scores_path": pocket1_dir / "calibration_scores.csv",
            "center": (0.0, 0.0, 0.0),
            "box_size": (10.0, 10.0, 10.0),
        },
        {
            "pocket_id": "P2",
            "dock_dir": pocket2_dir,
            "scores_path": pocket2_dir / "calibration_scores.csv",
            "center": (0.0, 0.0, 0.0),
            "box_size": (10.0, 10.0, 10.0),
        },
    ]

    (
        scores_by_pocket,
        no_score_by_pocket,
        dock_events_by_pocket,
    ) = pocket_eval._dock_calibrators_globally(
        pockets_plan=pockets_plan,
        prepared_scoring=prepared_scoring,
        receptor_pdbqt_path=receptor_pdbqt,
        vina_exe="vina",
        exhaustiveness=1,
        num_modes=1,
        verbosity=0,
        seed=0,
        cfg=cfg,
        logger=logging.getLogger("test"),
    )

    assert no_score_by_pocket.get("P1", 0) == 0
    assert no_score_by_pocket.get("P2", 0) == 0
    assert max_inflight >= 2

    p1_rows = scores_by_pocket.get("P1", [])
    p2_rows = scores_by_pocket.get("P2", [])
    assert len(p1_rows) == 1
    assert len(p2_rows) == 1
    assert p1_rows[0]["score"] == -7.0
    assert p2_rows[0]["score"] == -6.0

    p1_scores = pocket_eval._read_scores_csv(pocket1_dir / "calibration_scores.csv")
    p2_scores = pocket_eval._read_scores_csv(pocket2_dir / "calibration_scores.csv")
    assert p1_scores is not None
    assert p2_scores is not None
    assert len(p1_scores) == 1
    assert len(p2_scores) == 1
    assert dock_events_by_pocket["P1"][0]["status"] == "scored"
    assert dock_events_by_pocket["P2"][0]["status"] == "scored"
