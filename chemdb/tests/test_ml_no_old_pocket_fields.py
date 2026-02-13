from __future__ import annotations

from ml.config import FeaturesConfig
from ml.featurize import BigBindFeaturizer


def test_featurizer_feature_names_do_not_include_legacy_pocket_fields(tmp_path):
    featurizer = BigBindFeaturizer(
        bigbind_root=tmp_path,
        features=FeaturesConfig(
            ligand_descriptors=False,
            ligand_extra_descriptors=False,
            ligand_morgan_fp_bits=0,
            pocket_fpocket=True,
            vina_score=False,
        ),
        atlas_cfg={},
    )
    feature_names = set(featurizer.feature_names)

    assert "pocket_size_x" not in feature_names
    assert "pocket_size_y" not in feature_names
    assert "pocket_size_z" not in feature_names
    assert "pocket_volume" not in feature_names
    assert "num_pocket_residues" not in feature_names
    assert not any(name.startswith("pocket_frac_") for name in feature_names)

    assert "fpocket_druggability" in feature_names
    assert "fpocket_volume" in feature_names
    assert "fpocket_openness" in feature_names
