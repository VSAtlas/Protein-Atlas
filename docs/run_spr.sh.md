# tools/slurm/run_spr.sh

## 2026-02-15

- Added a Slurm array submission template named `tools/slurm/run_spr.sh` that launches `tools/run_relocated_mode.py` with shared `RUN_ID`.
- Exports distributed runtime env vars (`ATLAS_DISTRIBUTED_MODE`, task metadata, barrier settings) so the pipeline can shard work across array tasks and finalize on the leader task.

## 2026-02-16

- Updated `tools/slurm/run_spr.sh` to keep worker tasks free of in-array barrier/finalization behavior.
- Simplified to a 2-file flow:
  - `tools/slurm/run_spr.sh` now supports submit mode (`bash tools/slurm/run_spr.sh`) and automatically submits:
    1. the worker array job
    2. the dependency finalizer job (`afterany`) pointing to `tools/slurm/run_spr_finalize.sh`.
  - `tools/slurm/run_spr_finalize.sh` runs distributed reduction/finalization hooks once per run.
- Updated practical defaults for simple DUD-E scale runs:
  - worker CPU request default is `--cpus-per-task=128`
  - no default memory request line is set in the script
  - default worker args are now `--run-id <RUN_ID> -dude` (no `--fast`)
- Runtime CPU still follows scheduler allocation in distributed mode, so if granted CPUs differ from requested, Atlas uses granted CPUs automatically.
- Submit mode now passes `-p <partition>` to `sbatch` for both array and finalizer jobs.
  - Default partition is `normal`
  - Override with `PARTITION=<queue>` when launching.

## 2026-05-09

- Prefer `atlas slurm submit` over calling this template directly for new runs:
  - `atlas slurm submit --run-id <RUN_ID> --array 0-31%4 --cpus-per-task 8 --dry-run`
  - `atlas slurm progress <RUN_ID>`
  - `atlas slurm finalize <RUN_ID> --reconcile-only`
- `atlas slurm submit` passes `RUN_ID`, `ATLAS_MAIN_ARGS`, and `ATLAS_ALL_DIRS` into this script. Choose CPU requests and array concurrency for the target cluster allocation and queue policy.
- Worker execution now honors `ATLAS_MAIN_ARGS` and `ATLAS_ALL_DIRS`; the previous hardcoded resume run is no longer the default worker payload.
