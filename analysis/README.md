# Analysis Package Layout

This directory is now split into canonical subpackages with compatibility wrappers.

## Canonical modules

- `analysis/cli/`: command-line entrypoints.
- `analysis/dud_eval_core/`: DUD/DUD-E evaluation internals.
- `analysis/reporting/`: report/export/heatmap support code.
- `analysis/r/`: R assets used by analysis tooling.
- `analysis/ARCHITECTURE.md`: architecture + import/dependency rules.

## Compatibility policy (Phase 1)

Legacy flat modules (for example `analysis/run_report.py`, `analysis/dud_eval_metrics.py`) are thin wrappers that re-export from canonical modules.

- Wrappers emit `DeprecationWarning`.
- Planned removal: after 2 release cycles.
- New code must import canonical modules only.

Use `tools/check_analysis_imports.sh` to enforce canonical imports and prevent drift.

## Output policy

Do not write generated artifacts into `analysis/`.

Use run-scoped output roots, typically under:

- `docked/<RUN_ID>/...`
- `post_docked/<RUN_ID>/...`

`analysis/` should remain source-only.
