from __future__ import annotations

import json
import uuid
from pathlib import Path

import pandas as pd

from ml.experiment_runner import run_feature_ablation_experiments


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
        "ex_rec_pocket_file": f"{pocket}/{ex_rec_pdb}_pocket.pdb",
        "num_pocket_residues": 20,
        "pocket_size_x": 10.0,
        "pocket_size_y": 10.0,
        "pocket_size_z": 10.0,
        "pocket_center_x": 1.0,
        "pocket_center_y": 2.0,
        "pocket_center_z": 3.0,
    }


def _write_split(path: Path, rows: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(rows)
    frame.to_csv(path, index=False)


def test_experiment_runner_smoke(monkeypatch, tmp_path):
    def _mock_fpocket_loader(cfg, pdb_id, variant, ph_label, center, logger):
        if center is None:
            return None
        base = 0.62 if str(pdb_id).strip().upper() == "HOLD" else 0.41
        return {
            "fpocket_druggability": base,
            "fpocket_volume": 210.0,
            "fpocket_openness": 0.34,
            "fpocket_polar_fraction": 0.29,
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
        _build_row(smiles="c1ccccc1O", active=1, ex_rec_pdb="TRN1", pocket="PKT_TRN1"),
        _build_row(smiles="CCO", active=0, ex_rec_pdb="TRN1", pocket="PKT_TRN1"),
        _build_row(smiles="CCN(CC)CC", active=1, ex_rec_pdb="TRN2", pocket="PKT_TRN2"),
        _build_row(smiles="CCCC", active=0, ex_rec_pdb="TRN2", pocket="PKT_TRN2"),
        _build_row(smiles="c1ncccc1", active=1, ex_rec_pdb="TRN3", pocket="PKT_TRN3"),
        _build_row(smiles="CC(C)O", active=0, ex_rec_pdb="TRN3", pocket="PKT_TRN3"),
        _build_row(smiles="c1ccccc1N", active=1, ex_rec_pdb="HOLD", pocket="PKT_HOLD"),
        _build_row(smiles="CC(C)C", active=0, ex_rec_pdb="HOLD", pocket="PKT_HOLD"),
    ]
    val_rows = [
        _build_row(smiles="c1ccncc1", active=1, ex_rec_pdb="TRN4", pocket="PKT_TRN4"),
        _build_row(smiles="CCN", active=0, ex_rec_pdb="HOLD", pocket="PKT_HOLD"),
    ]
    test_rows = [
        _build_row(smiles="c1ccccc1", active=1, ex_rec_pdb="TRN5", pocket="PKT_TRN5"),
        _build_row(smiles="CCC", active=0, ex_rec_pdb="HOLD", pocket="PKT_HOLD"),
    ]
    _write_split(bigbind_dir / "activities_train.csv", train_rows)
    _write_split(bigbind_dir / "activities_val.csv", val_rows)
    _write_split(bigbind_dir / "activities_test.csv", test_rows)

    run_id = f"ml_experiment_smoke_{uuid.uuid4().hex[:10]}"
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
            "max_iter": 500,
        },
        "features": {
            "ligand_descriptors": True,
            "ligand_extra_descriptors": False,
            "ligand_morgan_fp_bits": 64,
            "pocket_fpocket": True,
            "vina_score": False,
        },
    }
    config_path = tmp_path / "ml_config.json"
    config_path.write_text(json.dumps(cfg), encoding="utf-8")

    output_dir, results = run_feature_ablation_experiments(str(config_path))
    results_path = output_dir / "experiment_results.csv"

    assert results_path.exists()
    assert not results.empty
    assert (output_dir / "features_train_baseline.csv").exists()
    assert (output_dir / "features_holdout_baseline.csv").exists()
    assert (output_dir / "morgan_onbits_train_baseline.csv").exists()
    assert (output_dir / "morgan_onbits_holdout_baseline.csv").exists()
    assert (output_dir / "features_all_variants.csv").exists()
    assert (output_dir / "morgan_onbits_all_variants.csv").exists()
    all_features = pd.read_csv(output_dir / "features_all_variants.csv")
    all_onbits = pd.read_csv(output_dir / "morgan_onbits_all_variants.csv")
    assert {"variant", "split", "row_number"}.issubset(set(all_features.columns))
    assert {"variant", "split", "row_number", "bit"}.issubset(set(all_onbits.columns))
    assert set(all_features["split"].unique().tolist()) == {"train", "holdout"}
    assert set(all_features["variant"].unique().tolist()) == {
        "baseline",
        "ligand_only",
        "pocket_only",
        "no_fp",
        "add_ligand_extras",
    }

    produced = set(results["variant"].tolist())
    assert produced == {
        "baseline",
        "ligand_only",
        "pocket_only",
        "no_fp",
        "add_ligand_extras",
    }
    baseline = results.loc[results["variant"] == "baseline"].iloc[0]
    assert int(baseline["fpocket_loaded_count"]) > 0
