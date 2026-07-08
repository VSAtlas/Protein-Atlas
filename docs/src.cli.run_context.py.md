# src/cli/run_context.py

## 2026-02-15

- Updated `_prepare_run_logfile()` to include Slurm array identifiers in the log filename when `SLURM_ARRAY_TASK_ID` is present, preventing multi-task log collisions.
