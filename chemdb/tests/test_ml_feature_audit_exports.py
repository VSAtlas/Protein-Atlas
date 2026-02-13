from __future__ import annotations

import csv
import json

import pandas as pd

from ml.config import FeaturesConfig
from ml.featurize import BigBindFeaturizer


def test_ml_feature_audit_exports(tmp_path):
    features = FeaturesConfig(
        ligand_descriptors=True,
        ligand_extra_descriptors=False,
        ligand_morgan_fp_bits=64,
        pocket_fpocket=False,
        vina_score=False,
    )
    featurizer = BigBindFeaturizer(
        bigbind_root=tmp_path,
        features=features,
        atlas_cfg={},
    )

    frame = pd.DataFrame(
        [
            {"lig_smiles": "CCO", "active": 1, "ex_rec_pdb": "EX1", "pocket": "PK1"},
            {"lig_smiles": "c1ccccc1", "active": 0, "ex_rec_pdb": "EX2", "pocket": "PK2"},
            {"lig_smiles": "CCN", "active": 1, "ex_rec_pdb": "EX3", "pocket": "PK3"},
        ]
    )

    rows = featurizer.transform(frame, audit_dir=tmp_path, audit_tag="train")

    features_csv = tmp_path / "features_train.csv"
    onbits_csv = tmp_path / "morgan_onbits_train.csv"
    summary_json = tmp_path / "feature_audit_train_summary.json"
    assert features_csv.exists()
    assert onbits_csv.exists()
    assert summary_json.exists()

    with features_csv.open("r", encoding="utf-8", newline="") as handle:
        feature_rows = list(csv.reader(handle))
    feature_header = feature_rows[0]
    feature_data_rows = feature_rows[1:]
    assert feature_header[:7] == [
        "row_number",
        "source_index",
        "active",
        "lig_smiles",
        "ex_rec_pdb",
        "pocket",
        "murcko_scaffold",
    ]
    assert feature_header[7:] == featurizer.metadata()["continuous_feature_names"]
    assert len(feature_data_rows) == rows.X.shape[0]

    with onbits_csv.open("r", encoding="utf-8", newline="") as handle:
        onbit_rows = list(csv.reader(handle))
    assert onbit_rows[0] == ["row_number", "bit"]
    bit_values = [int(entry[1]) for entry in onbit_rows[1:]]
    assert all(0 <= bit < featurizer.fingerprint_bits for bit in bit_values)

    cont_cols = len(featurizer.metadata()["continuous_feature_names"])
    assert len(bit_values) == int(rows.X[:, cont_cols:].nnz)

    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    assert summary["dropped_invalid_smiles"] == 0
    assert summary["number_of_rows_written"] == rows.X.shape[0]
    assert summary["fingerprint_bits"] == featurizer.fingerprint_bits
