# Notes: `src/docking/docking_stage_runner.py`

- Added early split for very large per-ligand core requests (`THREADS_PER_VINA >= 5`) to reduce fragmentation pressure.
- Switched retry docking admission to core-grant context (`acquire_global_cores`) and pass granted cores to Vina runtime.
- Main stage docking now forwards scheduler-granted cores to `run_docking_task(cpu_override=...)` so runtime CPU usage matches admission.
