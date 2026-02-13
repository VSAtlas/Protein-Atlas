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


def test_adversarial_validation_outputs(tmp_path):
    bigbind_root = tmp_path / "BigBindV1.5"
    bigbind_root.mkdir(parents=True, exist_ok=True)

    train_rows = [
        _row("CCO", 1, "TRN1", "PK1"),
        _row("CCN", 0, "TRN1", "PK1"),
        _row("CCC", 0, "TRN2", "PK2"),
        _row("CCCO", 1, "TRN2", "PK2"),
        _row("CC(C)O", 0, "TRN3", "PK3"),
        _row("CC(C)C", 1, "TRN3", "PK3"),
    ]
    # Holdout shifted toward aromatic molecules.
    holdout_rows = [
        _row("c1ccccc1", 1, "HOLD", "PKH"),
        _row("c1ccncc1", 0, "HOLD", "PKH"),
        _row("c1ccc(cc1)Cl", 1, "HOLD", "PKH"),
        _row("c1ncccc1", 0, "HOLD", "PKH"),
    ]
    _write_split(bigbind_root / "activities_train.csv", train_rows)
    _write_split(bigbind_root / "activities_val.csv", holdout_rows)
    _write_split(bigbind_root / "activities_test.csv", [])

    cfg = {
        "bigbind_dir": str(bigbind_root),
        "train_pdb": "HOLD",
        "dataset_split_mode": "legacy",
        "splits": ["train", "val"],
        "random_seed": 7,
        "run_id": f"adv_val_{uuid.uuid4().hex[:8]}",
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
        "registry": {"enabled": False},
        "diagnostics": {
            "adversarial_validation": True,
            "drop_feature_tests": False,
            "permutation_importance": False,
        },
    }
    cfg_path = tmp_path / "adv_cfg.json"
    cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
    out_dir, _ = run_pipeline(str(cfg_path))

    adv_path = out_dir / "adversarial_validation.json"
    assert adv_path.exists()
    payload = json.loads(adv_path.read_text(encoding="utf-8"))
    assert float(payload["adversarial_auc"]) > 0.6
