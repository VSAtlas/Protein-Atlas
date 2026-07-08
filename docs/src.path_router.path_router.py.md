# src/path_router/path_router.py

## 2026-02-16: Run-Scoped Docked Path Guard for Relocated Distributed Runs

- **Issue:** In relocated distributed execution, if router roots were initialized before `ATLAS_RUN_ID` was set, later `docked_dir(...)` calls could keep using a flat `DOCKED_DIR/<PDB>` path.
- **Fix:** `docked_dir(...)` now enforces run scoping at call time by appending `ATLAS_RUN_ID` when present and missing from the resolved docked root.
- **Impact:** Prevents mixed layouts where some worker writes land in `docked/<PDB>` instead of `docked/<RUN_ID>/<PDB>`.
