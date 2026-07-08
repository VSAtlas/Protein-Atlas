# tools/simulate_slurm_array.py

## 2026-02-15

- Added a local multi-task runner that simulates Slurm array environment variables (`SLURM_ARRAY_TASK_*`, `SLURM_CPUS_PER_TASK`) for multi-node workflow validation on a single machine.
- Supports dry-run, serial, and parallel local execution modes.

## 2026-03-02

- Refactored duplicate command construction/log emission into shared helpers.
- Behavior and CLI flags are unchanged; this was a maintainability-only cleanup.

## 2026-05-15

- Added `--main-args-file` for reproducible command profiles and `--summary-json` for compact run summaries.
- Added `--array` parsing for Slurm-style task specs, including concurrency throttles such as `0-31%8`.
- Dry-run and executed simulations now share the same command/env summary schema.
