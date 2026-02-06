from __future__ import annotations

import json
import uuid
from pathlib import Path

import pandas as pd

from ml.automl import run_automl


def _build_row(
    *,
    smiles: str,
    active: int,
    ex_rec_pdb: str,
    pocket: str,
) -> dict[str, object]:
    return {
        "lig_smiles": smiles,
        "active": active,
        "pocket": pocket,
        "ex_rec_pdb": ex_rec_pdb,
        "pocket_center_x": 1.0,
        "pocket_center_y": 2.0,
        "pocket_center_z": 3.0,
    }


def _write_split(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def test_ml_automl_grid_smoke(monkeypatch, tmp_path):
    def _mock_fpocket_loader(cfg, pdb_id, variant, ph_label, center, logger):
        if center is None:
            return None
        is_holdout = str(pdb_id).strip().upper() == "HOLD"
        return {
            "fpocket_druggability": 0.6 if is_holdout else 0.4,
            "fpocket_volume": 220.0,
            "fpocket_openness": 0.33,
            "fpocket_polar_fraction": 0.25,
            "fpocket_has_metal": 0.0,
            "fpocket_tier": 2.0,
        }

    monkeypatch.setattr(
        "ml.featurize.load_fpocket_metrics_for_ml",
        _mock_fpocket_loader,
    )

    bigbind_dir = tmp_path / "BigBindV1.5"
    bigbind_dir.mkdir(parents=True, exist_ok=True)

    train_rows = [
        _build_row(smiles="c1ccccc1O", active=1, ex_rec_pdb="TRN1", pocket="PKT1"),
        _build_row(smiles="c1ccccc1Cl", active=0, ex_rec_pdb="TRN1", pocket="PKT1"),
        _build_row(smiles="c1ccncc1", active=1, ex_rec_pdb="TRN2", pocket="PKT2"),
        _build_row(smiles="c1ccncc1C", active=0, ex_rec_pdb="TRN2", pocket="PKT2"),
        _build_row(smiles="N1CCCCC1", active=1, ex_rec_pdb="TRN3", pocket="PKT3"),
        _build_row(smiles="C1CCCCC1", active=0, ex_rec_pdb="TRN3", pocket="PKT3"),
    ]
    val_rows = [
        _build_row(smiles="c1ccccc1N", active=1, ex_rec_pdb="HOLD", pocket="PKTH"),
        _build_row(smiles="CC(C)C", active=0, ex_rec_pdb="HOLD", pocket="PKTH"),
    ]
    test_rows = [
        _build_row(smiles="c1ncccc1", active=1, ex_rec_pdb="TRN4", pocket="PKT4"),
        _build_row(smiles="CCC", active=0, ex_rec_pdb="HOLD", pocket="PKTH"),
    ]
    _write_split(bigbind_dir / "activities_train.csv", train_rows)
    _write_split(bigbind_dir / "activities_val.csv", val_rows)
    _write_split(bigbind_dir / "activities_test.csv", test_rows)

    run_id = f"ml_automl_smoke_{uuid.uuid4().hex[:10]}"
    cfg = {
        "bigbind_dir": str(bigbind_dir),
        "train_pdb": "HOLD",
        "splits": ["train", "val", "test"],
        "random_seed": 7,
        "max_rows": None,
        "run_id": run_id,
        "inner_scaffold_folds": 2,
        "atlas_cfg_path": None,
        "fpocket_center_columns": [
            "pocket_center_x",
            "pocket_center_y",
            "pocket_center_z",
        ],
        "fpocket_variant_column": "variant",
        "fpocket_ph_column": "pH",
        "model": {
            "type": "logreg",
            "C": 1.0,
            "class_weight": "balanced",
            "max_iter": 300,
        },
        "features": {
            "ligand_descriptors": True,
            "ligand_extra_descriptors": False,
            "ligand_morgan_fp_bits": 64,
            "pocket_fpocket": True,
            "vina_score": False,
        },
        "automl": {
            "feature_variants": ["baseline", "no_fp"],
            "model_families": ["logreg"],
            "grid": {
                "logreg": {
                    "C": [0.1, 1.0],
                    "penalty": ["l2"],
                    "solver": ["lbfgs"],
                    "class_weight": ["balanced"],
                }
            },
        },
    }
    cfg_path = tmp_path / "ml_automl_config.json"
    cfg_path.write_text(json.dumps(cfg), encoding="utf-8")

    run_dir = run_automl(str(cfg_path), mode="grid", top_k_report=2)

    assert (run_dir / "automl_trials.csv").exists()
    assert (run_dir / "best_model.joblib").exists()
    assert (run_dir / "best_trial.json").exists()
    assert (run_dir / "best_holdout_report.csv").exists()

    trials = pd.read_csv(run_dir / "automl_trials.csv")
    assert len(trials) == 4
