# ChemDB Test Fixtures

This directory contains small deterministic files used by acceptance and
regression tests. They are source fixtures, not runtime outputs.

- `prepped_ligands/fda_test_library_10/`: ten lightweight PDBQT ligands used by
  fast smoke and full-run tests.
- `prepped_ligands/test_library_10/`: manifest-only legacy DUD fixture record for
  tests that exercise library naming without depending on runtime roots.
- `extracted_ligands/`: paired source SDF fixtures for ligand-prep and synthetic
  benchmark helpers.
- `input_pdbs/TEST.pdb`: tiny receptor fixture used by `--test` smoke paths.

Generated run products still belong in ignored runtime directories such as
`input_pdbs/`, `prepped_ligands/`, `extracted_ligands/`, `docked/`, and `logs/`.
