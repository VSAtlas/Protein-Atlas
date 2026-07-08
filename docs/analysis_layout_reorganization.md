# analysis/ layout reorganization note (2026-02-14)

## Goal
- Reduce maintenance cost in `analysis/` by separating CLI entrypoints, reusable evaluation/reporting logic, and non-Python assets.
- Preserve backward compatibility for existing imports and module invocations.

## Structural change
- Added canonical packages:
  - `analysis/dud_eval_core/`
  - `analysis/reporting/`
  - `analysis/cli/`
  - `analysis/r/`
- Copied previous implementation modules into those canonical paths.
- Converted legacy flat modules in `analysis/*.py` into compatibility wrappers that alias to canonical modules.

## Compatibility guarantees kept
- Existing imports like `analysis.run_report`, `analysis.heatmap_html`, and `analysis.dud_eval_metrics` remain valid.
- Existing script/module entrypoints remain valid:
  - `python analysis/dud_eval.py`
  - `python -m analysis.run_report`
  - `python -m analysis.master_schema_export`
- New canonical module entrypoints are available:
  - `python -m analysis.cli.dud_eval`
  - `python -m analysis.cli.run_report`
  - `python -m analysis.cli.master_schema_export`
- New user-facing workflows should prefer Atlas wrappers:
  - `atlas analysis report <RUN_ID> --status-html`
  - `atlas analysis export-master <RUN_ID>`
  - `atlas analysis dud-eval <RUN_ID>`
  - `atlas analysis throughput integrity <RUN_ID>`
  - `atlas analysis interactions export <RUN_ID>`

## Additional behavior updates
- Post-run hooks now call canonical CLI modules (`analysis.cli.*`).
- DUD post-run output root is now repo-root `dud_eval` in post-run hooks and CLI defaults (legacy `docked/dud_eval` and `analysis/dud_eval` are still read for compatibility).
- Added a compatibility sync in `analysis/dud_eval_core/orchestrate.py` so explicit `--out-dir` + `--run-id` still emits legacy root-level summary/target artifacts.
- Heatmap source selection now prefers parquet dataset first:
  - `data/<run_id>/dataset/interactions`
  - then `data/<run_id>/heatmap_input.csv`
  - then `data/<run_id>/master_rows.csv`
- Selected-score schema is now `z_selected` / `z_selected_source` (hard migration from `t_selected` naming).

## Deprecation timeline
- Legacy flat wrappers in `analysis/*.py` now emit `DeprecationWarning`.
- Target removal: after 2 release cycles.
- Canonical imports are enforced by `tools/check_analysis_imports.sh`.

## Notes
- `analysis/dud_eval.py` remains as a file entrypoint for compatibility, so the canonical package name is `analysis.dud_eval_core` instead of `analysis.dud_eval`.
- `analysis/HeatMap.R` is retained; canonical copy is in `analysis/r/HeatMap.R`.
