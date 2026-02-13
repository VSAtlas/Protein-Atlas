from __future__ import annotations

import logging
from pathlib import Path

import pocket_eval


def test_calibration_emits_log_file(tmp_path: Path, monkeypatch, caplog) -> None:
    pdb_id = "TEST"
    run_id = "cal_run"
    docked_root = tmp_path / "docked"
    extracted_root = tmp_path / "extracted_ligands"

    cfg = {
        "RUN_ID": run_id,
        "DOCKED_DIR": str(docked_root),
        "LIGAND_EXTRACTED_DIR": str(extracted_root),
    }

    class DummyPaths:
        def __init__(self, docked: Path) -> None:
            self._docked = docked

        def docked_pdb_root(self) -> Path:
            return self._docked

    def fake_make_paths(cfg_obj, base_id, pdb_file):
        return DummyPaths(docked_root / run_id / base_id)

    monkeypatch.setattr(pocket_eval, "make_paths", fake_make_paths)

    def fake_run_calibrator_for_pdb(
        pdb_id: str,
        *,
        out_root: Path,
        log_dir: Path | None = None,
        run_tag: str | None = None,
        **_kwargs,
    ):
        pdb_norm = pdb_id.upper()
        out_dir = Path(out_root) / f"{pdb_norm}_calibrator"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "strong_binders.smi").write_text("C\n", encoding="utf-8")
        (out_dir / "non_binders.smi").write_text("CC\n", encoding="utf-8")
        (out_dir / "weak_binders.smi").write_text("CCC\n", encoding="utf-8")
        log_file = None
        if log_dir and run_tag:
            run_log_dir = Path(log_dir) / str(run_tag)
            run_log_dir.mkdir(parents=True, exist_ok=True)
            log_file = run_log_dir / f"{pdb_norm}.log"
            log_file.write_text("calibration log\n", encoding="utf-8")
        return {"status": "ok", "log_file": str(log_file) if log_file else None}

    monkeypatch.setattr(pocket_eval, "run_calibrator_for_pdb", fake_run_calibrator_for_pdb)

    logger = logging.getLogger("pocket-eval-calibration-test")
    caplog.set_level(logging.INFO, logger=logger.name)

    rows = pocket_eval.get_calibration_set_for_pdb(pdb_id, cfg, logger)
    assert rows

    expected_log = docked_root / run_id / "logs" / run_id / f"{pdb_id}.log"
    assert expected_log.exists()
    assert "[calibration.start]" in caplog.text
    assert "[calibration.end]" in caplog.text
