from __future__ import annotations

import pandas as pd

from ml.data.bigbind import load_train_holdout_from_bigbind


def _write_split(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False)


def test_bigbind_standard_split_mode_uses_ordered_splits(tmp_path):
    root = tmp_path / "BigBindV1.5"
    root.mkdir(parents=True, exist_ok=True)
    base = {
        "lig_smiles": "CCO",
        "active": 1,
        "pocket": "TGT_A/1abc_A_rec_pocket.pdb",
        "ex_rec_pdb": "NOT_HOLD",
    }
    _write_split(root / "activities_train.csv", [base | {"lig_smiles": "CCO"}])
    _write_split(root / "activities_val.csv", [base | {"lig_smiles": "CCC"}])
    _write_split(root / "activities_test.csv", [base | {"lig_smiles": "CCN"}])

    ds = load_train_holdout_from_bigbind(
        bigbind_dir=root,
        train_pdb="HOLD",
        splits=("train", "val", "test"),
        max_rows=None,
        random_seed=7,
        dataset_split_mode="standard",
    )
    assert ds.dataset_split_mode == "standard"
    assert len(ds.train_df) == 1
    assert len(ds.val_df) == 1
    assert len(ds.test_df) == 1
    assert len(ds.holdout_df) == 1
    assert ds.holdout_match_column == "split:val"
    assert ds.train_df.iloc[0]["split"] == "train"
    assert ds.holdout_df.iloc[0]["split"] == "val"


def test_bigbind_legacy_split_mode_uses_train_pdb_holdout(tmp_path):
    root = tmp_path / "BigBindV1.5"
    root.mkdir(parents=True, exist_ok=True)
    train_rows = [
        {
            "lig_smiles": "CCO",
            "active": 1,
            "pocket": "TGT_A/1abc_A_rec_pocket.pdb",
            "ex_rec_pdb": "TRN1",
        },
        {
            "lig_smiles": "CCC",
            "active": 0,
            "pocket": "TGT_H/2abc_A_rec_pocket.pdb",
            "ex_rec_pdb": "HOLD",
        },
    ]
    _write_split(root / "activities_train.csv", train_rows)
    _write_split(root / "activities_val.csv", [])
    _write_split(root / "activities_test.csv", [])

    ds = load_train_holdout_from_bigbind(
        bigbind_dir=root,
        train_pdb="HOLD",
        splits=("train",),
        max_rows=None,
        random_seed=7,
        dataset_split_mode="legacy",
    )
    assert ds.dataset_split_mode == "legacy"
    assert len(ds.train_df) == 1
    assert len(ds.holdout_df) == 1
    assert ds.holdout_match_column == "ex_rec_pdb"
