# analysis/reporting/master_schema_export.py

## Pathing fix (2026-02-14)

- Updated DUD summary discovery to prefer canonical repo-root outputs:
  - `<repo>/analysis/<decoy_prefix>_eval/<run_id>`
- Added compatibility fallbacks for historical layouts:
  - `<repo>/<decoy_prefix>_eval/<run_id>`
  - `<repo>/analysis/dud_eval/<run_id>`
  - `<repo>/docked/<decoy_prefix>_eval/<run_id>`
