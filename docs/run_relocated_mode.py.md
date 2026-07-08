# tools/run_relocated_mode.py

## 2026-02-16

- Added an exclusive extraction lock (`.prepped_extract.lock`) around archive extraction so concurrent Slurm array workers do not race while materializing `PREPPED_LIGANDS_DIR`.
- Kept extraction idempotent: workers re-check payload presence after lock acquisition and skip extraction when already populated.
- Updated temporary extraction directory placement to use `prepped_root.parent` instead of default system temp roots.
