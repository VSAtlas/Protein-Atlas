from __future__ import annotations

import numpy as np

from ml.pocket_features import (
    compute_pocket_features,
    load_pdb_atoms,
    metal_context_features,
    pocket_geometry_features,
    pocket_residue_composition,
)


def test_pocket_features_from_tiny_pdb(tmp_path):
    pdb_path = tmp_path / "tiny_pocket.pdb"
    pdb_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   ASP A  10      10.000  10.000  10.000  1.00 20.00           N",
                "ATOM      2  CA  ASP A  10      11.000  10.200  10.300  1.00 20.00           C",
                "ATOM      3  N   LYS A  11      12.100  11.000  10.500  1.00 20.00           N",
                "ATOM      4  CA  PHE A  12      13.200  11.800  10.800  1.00 20.00           C",
                "HETATM    5 ZN   ZN  A 500      14.000  12.000  11.000  1.00 20.00          ZN",
                "END",
            ]
        ),
        encoding="utf-8",
    )

    atoms = load_pdb_atoms(pdb_path)
    comp = pocket_residue_composition(atoms)
    geom = pocket_geometry_features(atoms)
    metals = metal_context_features(atoms)
    combined = compute_pocket_features(pdb_path, cache_root=tmp_path / "cache")

    frac_sum = (
        comp["pocket_frac_polar"]
        + comp["pocket_frac_hydrophobic"]
        + comp["pocket_frac_aromatic"]
        + comp["pocket_frac_charged_pos"]
        + comp["pocket_frac_charged_neg"]
    )
    assert np.isclose(frac_sum, 1.0, atol=1e-6)
    assert metals["pocket_metal_zn"] >= 1.0
    assert np.isfinite(geom["pocket_rgyr"])
    assert np.isfinite(geom["pocket_bbox_vol"])
    assert np.isfinite(geom["pocket_pca_ratio1"])
    assert np.isfinite(geom["pocket_pca_ratio2"])
    assert combined["pocket_metal_zn"] >= 1.0
