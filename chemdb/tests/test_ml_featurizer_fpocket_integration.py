from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from ml.config import FeaturesConfig
from ml.featurize import BigBindFeaturizer


def test_featurizer_includes_fpocket_columns_and_values(monkeypatch, tmp_path):
    def _mock_fpocket_loader(cfg, pdb_id, variant, ph_label, center, logger):
        assert isinstance(cfg, dict)
        assert pdb_id == "EX4M"
        assert center == (1.0, 2.0, 3.0)
        return {
            "fpocket_druggability": 0.45,
            "fpocket_volume": 230.0,
            "fpocket_openness": 0.32,
            "fpocket_polar_fraction": 0.28,
            "fpocket_has_metal": 1.0,
            "fpocket_tier": 2.0,
        }

    monkeypatch.setattr(
        "ml.featurize.load_fpocket_metrics_for_ml",
        _mock_fpocket_loader,
    )

    features = FeaturesConfig(
        ligand_descriptors=False,
        ligand_extra_descriptors=False,
        ligand_morgan_fp_bits=0,
        pocket_fpocket=True,
        vina_score=False,
    )
    featurizer = BigBindFeaturizer(
        bigbind_root=tmp_path,
        features=features,
        atlas_cfg={"OVERALL_DIR": str(tmp_path)},
        logger=logging.getLogger("test.ml.featurizer.fpocket"),
    )

    frame = pd.DataFrame(
        [
            {
                "lig_smiles": "CCO",
                "active": 1,
                "ex_rec_pdb": "EX4M",
                "pocket_center_x": 1.0,
                "pocket_center_y": 2.0,
                "pocket_center_z": 3.0,
            }
        ]
    )

    rows = featurizer.transform(frame)
    assert featurizer.feature_names == [
        "fpocket_druggability",
        "fpocket_volume",
        "fpocket_openness",
        "fpocket_polar_fraction",
        "fpocket_has_metal",
        "fpocket_tier",
    ]
    np.testing.assert_allclose(
        rows.X.toarray()[0],
        np.asarray([0.45, 230.0, 0.32, 0.28, 1.0, 2.0], dtype=np.float32),
        rtol=0,
        atol=1e-6,
    )
    assert featurizer.loaded_fpocket_count == 1
    assert featurizer.missing_center_count == 0
    assert featurizer.missing_info_count == 0
