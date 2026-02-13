from __future__ import annotations

import logging
import threading
from pathlib import Path

import pocket_eval


def test_pocket_eval_calibrator_parallelism(tmp_path: Path, monkeypatch) -> None:
    cfg = {"CPU": 2, "THREADS_PER_VINA": 1}
    prepared_scoring = [
        {
            "ligand_id": "lig_1",
            "label": "strong",
            "pdbqt_path": str(tmp_path / "lig_1.pdbqt"),
        },
        {
            "ligand_id": "lig_2",
            "label": "non",
            "pdbqt_path": str(tmp_path / "lig_2.pdbqt"),
        },
    ]
    dock_dir = tmp_path / "dock"
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
        score = -7.0 if ligand_id.endswith("1") else -6.0
        return out_path, score

    def fake_write_vina_config(*_args, **_kwargs) -> None:
        return None

    monkeypatch.setattr(pocket_eval, "run_docking_task", fake_run_docking_task)
    monkeypatch.setattr(pocket_eval, "_write_vina_config", fake_write_vina_config)

    max_workers, threads_per_vina = pocket_eval._calc_vina_parallel_workers(
        cfg, len(prepared_scoring), logging.getLogger("test")
    )
    scores_rows, no_score = pocket_eval._dock_calibrators_for_pocket(
        pocket_id="P1",
        prepared_scoring=prepared_scoring,
        receptor_pdbqt_path=receptor_pdbqt,
        center=(0.0, 0.0, 0.0),
        box_size=(10.0, 10.0, 10.0),
        dock_dir=dock_dir,
        vina_exe="vina",
        max_workers=max_workers,
        threads_per_vina=threads_per_vina,
        exhaustiveness=1,
        num_modes=1,
        verbosity=0,
        seed=0,
        logger=logging.getLogger("test"),
    )

    assert max_workers >= 2
    assert len(scores_rows) == 2
    assert no_score == 0
    assert max_inflight >= 2
