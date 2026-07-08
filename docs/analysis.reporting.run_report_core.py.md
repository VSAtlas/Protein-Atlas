# analysis/reporting/run_report_core.py source-policy note (2026-02-14)

## Problem
- Heatmap source resolution preferred CSV over parquet dataset, conflicting with canonical dataset-first behavior.

## Change
- Updated `_resolve_heatmap_source` precedence to:
  1. `data/<run_id>/dataset/interactions`
  2. `data/<run_id>/heatmap_input.csv`
  3. `data/<run_id>/master_rows.csv`

## Why
- Standardizes source resolution for report generation.
- Ensures reports prefer normalized parquet datasets when present.

## Vina fallback mini-heatmap (2026-02-27)

- Added `-vina/--vina` to `analysis.cli.run_report` (wired in `run_report_core.py`).
- New fallback behavior when `data/<run_id>/master_rows.csv` is missing:
  - Instead of hard-failing, the command can build a temporary Vina-derived master table from stage outputs under `docked/<run_id>/<pdb>/<variant>/<ph>/stage3|stage2|stage1`.
  - Scores are parsed from pose `.pdbqt` files (`REMARK VINA RESULT`) and converted to heatmap-compatible polarity with `z_selected = -score`.
- Fallback targets are limited to unprefixed `stage3`-ready directories so the pilot heatmap reflects FDA-progressed targets rather than DUD-only control artifacts.
- Emits:
  - `data/<run_id>/heatmap_input.csv` (or user-provided `--heatmap-csv-path`)
  - `data/<run_id>/master_rows_vina.csv` as the intermediate source for heatmap CSV generation
  - `data/<run_id>/report.html` (or user-provided `--html-path`) using existing interactive heatmap renderer.
