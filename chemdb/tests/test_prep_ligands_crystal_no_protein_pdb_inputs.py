from __future__ import annotations

from pathlib import Path

from prep_ligands.prep_ligands_crystal import prep_ligands_from_pdb

PROTEIN_PDB = """\
HEADER    TEST PDB WITH PROTEIN + LIGAND
ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00 20.00           N
HETATM    2  C1  LIG A   3      10.000  10.000   1.000  1.00 20.00           C
END
"""

LIGAND_PDB = """\
HETATM    1  C1  LIG A   3      10.000  10.000   1.000  1.00 20.00           C
HETATM    2  O1  LIG A   3      10.500  10.500   1.500  1.00 20.00           O
END
"""


def test_prep_ligands_from_pdb_filters_protein_pdb_inputs(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()

    (root / "ABCD.pdb").write_text(PROTEIN_PDB)
    (root / "LIG_A3.pdb").write_text(LIGAND_PDB)

    out_dir = tmp_path / "out"
    mol2_dir = tmp_path / "mol2"

    result = prep_ligands_from_pdb(
        ligand_output_dir=root,
        ligands_mol2_dir=mol2_dir,
        prepped_ligands_dir=out_dir,
        dry_run=True,
    )

    assert result is not None
    selected_inputs = result["selected_inputs"]
    planned_outputs = result["planned_outputs"]

    assert "LIG_A3.pdb" in selected_inputs
    assert "ABCD.pdb" not in selected_inputs
    assert "LIG_A3.sanitized.pdbqt" in planned_outputs
    assert "ABCD.sanitized.pdbqt" not in planned_outputs
