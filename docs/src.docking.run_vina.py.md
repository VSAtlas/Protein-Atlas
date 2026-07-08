# Notes: `src/docking/run_vina.py`

- Added optional `cpu_override` parameter to `run_docking_task(...)`.
- When provided, Vina now uses a single effective CPU setting:
  - If config already contains `cpu = ...`, a temporary rewritten config is used and CLI `--cpu` is not appended.
  - If config does not contain `cpu = ...`, CLI `--cpu <n>` is appended once.
- Temporary rewritten config filenames are unique per task invocation to avoid collisions under parallel docking.
- This fixes the Vina parse failure `option '--cpu' cannot be specified more than once` and prevents thread-count drift against granted cores.
