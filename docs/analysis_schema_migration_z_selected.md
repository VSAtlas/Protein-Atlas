# Analysis Score Schema Migration: `z_selected` (2026-02-14)

## Scope
This migration standardizes selected-score naming in analysis/reporting paths.

## Breaking changes
- Canonical selected-score fields are now:
  - `z_selected`
  - `z_selected_source`
- Legacy `t_selected`/`t_selected_source` names are no longer accepted in analysis readers/writers.

## Updated modules
- `analysis/cli/convert_interactions_to_parquet.py`
- `analysis/reporting/heatmap_html.py`
- `analysis/reporting/run_report_core.py`
- `analysis/reporting/master_schema_export.py`
- `analysis/r/HeatMap.R`
- `analysis/HeatMap.R`

## Source policy
Run report heatmap source preference is now parquet-first:
1. `data/<run_id>/dataset/interactions`
2. `data/<run_id>/heatmap_input.csv`
3. `data/<run_id>/master_rows.csv`

## Caller migration
Use canonical imports/modules:
- `analysis.cli.*`
- `analysis.dud_eval_core.*`
- `analysis.reporting.*`

Legacy flat wrappers remain temporarily with deprecation warnings and are scheduled for removal after 2 release cycles.
