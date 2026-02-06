from __future__ import annotations

import json
import uuid
from pathlib import Path

import pandas as pd

from ml.train_and_eval import run_pipeline


def _write_split(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def _row(smiles: str, active: int, ex_rec_pdb: str, pocket: str) -> dict[str, object]:
    return {
        "lig_smiles": smiles,
        "active": active,
        "ex_rec_pdb": ex_rec_pdb,
        "pocket": pocket,
    }


def test_train_and_eval_outputs_include_calibration_artifacts(tmp_path):
    bigbind_root = tmp_path / "BigBindV1.5"
    bigbind_root.mkdir(parents=True, exist_ok=True)

    _write_split(
        bigbind_root / "activities_train.csv",
        [
            _row("CCO", 1, "TRN1", "PKT1"),
            _row("CCN", 0, "TRN1", "PKT1"),
            _row("c1ccccc1", 1, "TRN2", "PKT2"),
            _row("CCCC", 0, "TRN2", "PKT2"),
            _row("CC(C)O", 1, "TRN3", "PKT3"),
            _row("CC(C)C", 0, "TRN3", "PKT3"),
        ],
    )
    _write_split(
        bigbind_root / "activities_val.csv",
        [
            _row("c1ccncc1", 1, "HOLD", "PKTH"),
            _row("CCC", 0, "HOLD", "PKTH"),
        ],
    )
    _write_split(
        bigbind_root / "activities_test.csv",
        [
            _row("c1ncccc1", 1, "HOLD", "PKTH"),
            _row("CCCO", 0, "HOLD", "PKTH"),
        ],
    )

    run_id = f"ml_train_eval_cal_{uuid.uuid4().hex[:10]}"
    cfg = {
        "bigbind_dir": str(bigbind_root),
        "train_pdb": "HOLD",
        "splits": ["train", "val", "test"],
        "dataset_split_mode": "legacy",
        "random_seed": 7,
        "max_rows": None,
        "run_id": run_id,
        "inner_scaffold_folds": 2,
        "model_family": "logreg",
        "model_params": {
            "C": 1.0,
            "class_weight": "balanced",
            "max_iter": 500,
            "solver": "liblinear",
            "penalty": "l2",
        },
        "calibration": {
            "enabled": True,
            "method": "sigmoid",
            "cv_folds": 2,
            "seed": 17,
        },
        "metrics": {
            "report_brier": True,
            "report_ece": True,
            "ece_bins": 10,
        },
        "registry": {"enabled": True},
        "features": {
            "ligand_descriptors": True,
            "ligand_extra_descriptors": False,
            "ligand_morgan_fp_bits": 64,
            "pocket_fpocket": False,
            "vina_score": False,
        },
    }
    cfg_path = tmp_path / "ml_train_eval_config.json"
    cfg_path.write_text(json.dumps(cfg), encoding="utf-8")

    output_dir, _metrics = run_pipeline(str(cfg_path))

    assert (output_dir / "model.joblib").exists()
    assert (output_dir / "metrics_report.json").exists()
    assert (output_dir / "metrics_report.csv").exists()
    assert (output_dir / "holdout_predictions.csv").exists()
    assert (output_dir / "config_snapshot.txt").exists()
    assert (output_dir / "config_snapshot.json").exists()

    assert (output_dir / "calibrated_holdout_predictions.csv").exists()
    assert (output_dir / "calibration_oof_train.csv").exists()
    assert (output_dir / "metrics_calibration.json").exists()
    assert (output_dir / "reliability_plot.png").exists()
    assert (output_dir / "registry.json").exists()

    pred_df = pd.read_csv(output_dir / "calibrated_holdout_predictions.csv")
    assert {
        "row_id",
        "lig_smiles",
        "active",
        "ex_rec_pdb",
        "pocket",
        "scaffold",
        "raw_score",
        "prob_uncalibrated",
        "prob_calibrated",
    }.issubset(set(pred_df.columns))
    assert len(pred_df) > 0
