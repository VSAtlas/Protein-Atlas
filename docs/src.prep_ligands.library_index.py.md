# `src/prep_ligands/library_index.py` bugfix notes

## 2026-02-15 — manifest freshness hardening

- Issue: Strict manifest validation only compared root `mtime`, which can miss stale manifest content in some nested-library scenarios.
- Change:
  - Added `child_dir_mtimes` freshness metadata to manifests.
  - On load with strict checks, require both root `mtime` and `child_dir_mtimes` to match before trusting cached entries.
  - Persist the new freshness metadata from both manifest build paths.
- Effect: Library manifest cache is less likely to trust stale or partial data when ligand directories changed under the library root.
