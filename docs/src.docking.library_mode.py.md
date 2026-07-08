# `src/docking/library_mode.py` compatibility note

## 2026-02-15 — ligand root fallback compatibility

- Restored fallback behavior for test/runtime configs that still provide `OUTPUT_LIGANDS_DIR` instead of `PREPPED_LIGANDS_DIR`.
- Root selection now resolves as:
  1. `PREPPED_LIGANDS_DIR`
  2. `OUTPUT_LIGANDS_DIR`
  3. `prepped_ligands` (default)
- This keeps mapped-library root resolution compatible with existing tests and legacy config shapes.
