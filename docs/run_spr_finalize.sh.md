# tools/slurm/run_spr_finalize.sh

## 2026-02-16

- Added a Slurm finalizer job script that calls `tools/finalize_distributed_run.py --run-id <RUN_ID>`.
- Intended to run once per distributed run as a dependency job after worker-array completion.
- In the simplified 2-file flow, this script is submitted automatically by `tools/slurm/run_spr.sh` submit mode.
