# Analysis Architecture and Import Rules

## Canonical module boundaries

- `analysis/cli/`: CLI entrypoints only.
- `analysis/dud_eval_core/`: DUD evaluation core logic.
- `analysis/reporting/`: report/export/heatmap logic.
- `analysis/r/`: R assets.

## Dependency direction

- `analysis/cli` may import from `analysis/dud_eval_core` and `analysis/reporting`.
- `analysis/dud_eval_core` must not import from `analysis/cli`.
- `analysis/reporting` must not import from `analysis/cli`.
- New code must not depend on legacy flat wrappers in `analysis/*.py`.

## Legacy wrapper policy

Flat wrapper modules in `analysis/*.py` are compatibility shims only.

- No new feature logic in wrappers.
- No new imports should target wrappers.
- Wrappers emit `DeprecationWarning` and are scheduled for removal after 2 release cycles.

## Score schema policy

`z_selected` is the canonical selected-score column.

- New readers/writers must require or emit `z_selected`.
- Do not introduce new `t_selected` columns or aliases.

## Source selection policy

Run report heatmap source preference is:

1. `data/<run_id>/dataset/interactions` (parquet)
2. `data/<run_id>/heatmap_input.csv`
3. `data/<run_id>/master_rows.csv`

## Drift prevention

Run `tools/check_analysis_imports.sh` to enforce canonical imports and block legacy wrapper imports in non-wrapper code.
