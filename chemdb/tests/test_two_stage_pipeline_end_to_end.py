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


def _make_dataset(root: Path) -> None:
    _write_split(
        root / "activities_train.csv",
        [
            _row("c1ccccc1O", 1, "TRN1", "PK1"),
            _row("c1ccncc1", 1, "TRN1", "PK1"),
            _row("c1ccccc1N", 1, "TRN2", "PK2"),
            _row("CCO", 0, "TRN1", "PK1"),
            _row("CCN", 0, "TRN1", "PK1"),
            _row("CCCC", 0, "TRN2", "PK2"),
            _row("CCCO", 0, "TRN2", "PK2"),
            _row("CC(C)C", 0, "TRN2", "PK2"),
        ],
    )
    _write_split(
        root / "activities_val.csv",
        [
            _row("c1ncccc1", 1, "HOLD", "PKH"),
            _row("c1ccccc1Cl", 1, "HOLD", "PKH"),
            _row("CCC", 0, "HOLD", "PKH"),
            _row("CCCOC", 0, "HOLD", "PKH"),
        ],
    )
    _write_split(root / "activities_test.csv", [])


def _config(base: dict, *, run_id: str, pipeline_enabled: bool) -> dict:
    cfg = dict(base)
    cfg["run_id"] = run_id
    cfg["pipeline"] = {
        "enabled": pipeline_enabled,
        "stage1": {
            "model_family": "logreg",
            "model_params": {"C": 1.0, "solver": "liblinear", "penalty": "l2", "max_iter": 500},
            "keep_top_pct": 60.0,
            "keep_within_group": True,
            "group_key": "target_pocket",
        },
        "stage2": {
            "model_family": "logreg",
            "model_params": {"C": 1.0, "solver": "liblinear", "penalty": "l2", "max_iter": 500},
            "group_key": "target_pocket",
        },
    }
    return cfg


def test_two_stage_pipeline_end_to_end(tmp_path):
    bigbind_root = tmp_path / "BigBindV1.5"
    bigbind_root.mkdir(parents=True, exist_ok=True)
    _make_dataset(bigbind_root)

    base = {
        "bigbind_dir": str(bigbind_root),
        "train_pdb": "HOLD",
        "dataset_split_mode": "legacy",
        "splits": ["train", "val"],
        "random_seed": 7,
        "features": {
            "ligand_descriptors": True,
            "ligand_extra_descriptors": False,
            "ligand_morgan_fp_bits": 64,
            "pocket_fpocket": False,
            "pocket_features": False,
            "vina_score": False,
        },
        "model_family": "logreg",
        "model_params": {
            "C": 1.0,
            "solver": "liblinear",
            "penalty": "l2",
            "max_iter": 500,
            "class_weight": "balanced",
        },
        "calibration": {"enabled": False},
        "registry": {"enabled": False},
        "diagnostics": {
            "adversarial_validation": False,
            "drop_feature_tests": False,
            "permutation_importance": False,
        },
    }

    cfg_baseline = _config(base, run_id=f"two_stage_base_{uuid.uuid4().hex[:8]}", pipeline_enabled=False)
    cfg_pipeline = _config(base, run_id=f"two_stage_pipe_{uuid.uuid4().hex[:8]}", pipeline_enabled=True)

    cfg_base_path = tmp_path / "cfg_base.json"
    cfg_pipe_path = tmp_path / "cfg_pipe.json"
    cfg_base_path.write_text(json.dumps(cfg_baseline), encoding="utf-8")
    cfg_pipe_path.write_text(json.dumps(cfg_pipeline), encoding="utf-8")

    base_dir, _ = run_pipeline(str(cfg_base_path))
    pipe_dir, _ = run_pipeline(str(cfg_pipe_path))

    base_metrics = json.loads((base_dir / "metrics_report.json").read_text(encoding="utf-8"))
    pipe_metrics = json.loads((pipe_dir / "metrics_report.json").read_text(encoding="utf-8"))

    assert float(pipe_metrics["EF@1%"]) >= float(base_metrics["EF@1%"]) - 1e-9
    assert (pipe_dir / "stage1_predictions_holdout.csv").exists()
    assert (pipe_dir / "stage2_rank_scores_holdout.csv").exists()
    assert (pipe_dir / "pipeline_summary.json").exists()
