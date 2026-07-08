# analysis/dud_eval_core/orchestrate.py

## Pathing + Aggregation update (2026-02-14)

- Updated CLI default `--out-dir` to canonical `analysis/dud_eval`.
- Retired legacy flat mirroring; run-level artifacts stay under `analysis/dud_eval/<run_id>/`.
- Added class-level DUD-E aggregation map (Protease, Nuclear receptor, GPCR, Ion channel, Cytochrome P450, Kinase, Miscellaneous, Other enzyme).
- Added pooled class-level plotting under:
  - `analysis/dud_eval/<run_id>/class_aggregates/docking/<class>/`
  - `analysis/dud_eval/<run_id>/class_aggregates/consensus/<class>/`
  - `analysis/dud_eval/<run_id>/class_aggregates/post_docked_scorch/<class>/`
- Extended `summary_macro.tsv` with class aggregate rows (`class_<class_slug>`) after `macro_avg`.
- Added placeholder per-target `metrics.tsv` emission and pre-creation of `consensus/` and `post_docked/` run folders so run layout remains consistent even when source CSVs are missing or degenerate.

## Class figure contract hardening (2026-02-16)

- Class aggregate generation is now deterministic per run and mode:
  - Always creates `class_aggregates/docking|consensus|post_docked_scorch/<class_slug>/`.
- Sparse/degenerate class cases no longer silently skip:
  - Writes placeholder `metrics.tsv` with explicit `status_reason` and target-count context.
  - Optional placeholder figure emission controlled by `--class-placeholder-images`.
- Added run-level class generation ledger:
  - `analysis/dud_eval/<run_id>/class_aggregates/class_generation_report.tsv`
  - Includes mode, class, status, present/evaluable target counts, pooled row counts.
- Added CLI controls:
  - `--class-min-targets` to require minimum class coverage before computing class metrics.
  - `--class-placeholder-images/--no-class-placeholder-images` for sparse-class PNG placeholders.
- `summary_macro.tsv` now always appends all class rows (`class_<slug>`) and includes `status_reason` to distinguish real metrics from placeholders.

## Missing input_pdb fallback discovery fix (2026-02-16)

- Non-manifest discovery no longer hard-skips a target directory when `input_pdbs/<PDB>.pdb` is absent.
- New behavior:
  - Logs a warning (`discover.fallback`) for missing input PDB.
  - Continues scanning run outputs for `docking_score_long.csv` / `dud_docking_score_long.csv`.
- This prevents interrupted sharded runs (empty/incomplete manifest) from failing discovery when output CSVs exist but local `input_pdbs` is incomplete.
- Post-docked reranked discovery remains enabled through existing target specs and `_resolve_reranked_scorch_path(...)`; no control-log source changes were introduced in this patch.

## CLI root precedence for CSV discovery (2026-02-16)

- `_discover_docking_csvs(...)` now prioritizes the explicit scan root passed from CLI args (`--docked-root` and run-id-derived scan roots).
- Path-router/config-derived `docked_pdb_root()` is retained as a secondary fallback only.
- This avoids false "no docking_score_long.csv" failures when config path roots point to a different filesystem (for example `/work2/...`) than the run artifacts being evaluated (for example `/scratch/...`).
