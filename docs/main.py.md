# main.py change note

## 2026-02-13

- Updated `--verify-tools` to include `SCORCH` as an optional check target.
- Verification still fails only on required tools; missing SCORCH now reports as `optional missing` and does not return a failure code.
- Added optional SCORCH env probe (`micromamba run -n/-p ... python -V`) during `--verify-tools`.
- Updated SCORCH env probe to use the shared runner resolver (`MICROMAMBA_EXE`/`MAMBA_EXE`, PATH, `~/micromamba/bin/micromamba`, then conda), avoiding PATH-only false negatives.
- Updated help text to clarify required vs optional verification behavior.

## 2026-02-14

- Refactored runtime side effects into `_bootstrap_runtime()` to avoid import-time mutation.
- Moved warning filters, stream reconfigure calls, and `install_debug_makedirs()` into bootstrap.
- Removed file-level Ruff suppression and made imports top-level with no `E402` violations.
- Called `_bootstrap_runtime()` from `main()` before pipeline logic and from script wrapper for parity.

## 2026-02-15

- Added distributed Slurm-array execution support (env-driven): deterministic per-PDB sharding using `ATLAS_DISTRIBUTED_MODE=slurm_array` and `SLURM_ARRAY_TASK_*`.
- Preserved current per-node behavior: each task still uses existing adaptive intra-node parallelism, now capped to scheduler CPU allocation.
- Added leader-only run finalization/report hooks with worker completion markers under `manifests/<run_id>/distributed/markers`.
- Added distributed barrier/reduce flow so leader waits for workers, merges per-protein state, and finalizes a single `run_manifest.yaml`.

## 2026-02-16

- Switched distributed array behavior to worker-only completion inside `main.py`: each task writes marker state and exits without in-array barrier waiting.
- Run-level reduction/finalization/report hooks are now intended for a dedicated dependency job (`afterok` on the array), not for a leader task inside the array.
- Hardened DUD-only selection (`TEST_MODE_ENABLE=dud` / `-dude`): `TEST_LIBRARY_MAP` must be non-empty and cover all selected PDB IDs, otherwise the run exits with status 2.
- Mixed modes that include DUD (for example `fda+dud`, `hmdb+dud`) no longer hard-fail on unmapped proteins; unmapped PDBs continue with non-DUD libraries and emit a warning.
- Updated CPU policy for distributed mode: when `ATLAS_DISTRIBUTED_MODE=slurm_array` is enabled, runtime `CPU` and `GLOBAL_SCHEDULER_CPUS` are set to detected allocated CPUs (scheduler env), so per-run config edits are not required to use full assigned cores.
- Non-distributed mode keeps prior behavior: runtime CPU remains capped by `min(config CPU, allocated CPU)`.
- Added persistent mixed SCORCH queue startup and drain points in main flow:
  - queue starts once per run when SCORCH is enabled
  - proteins submit SCORCH work to shared queue after docking
  - queue is drained before final reporting/integrity enforcement
- Switched global scheduling work units from per-PDB to per-(PDB,variant,pH) combo:
  - main builds explicit combo items per variant using discovered pH tags
  - each combo runs as an independent schedulable task with `ph_tag`-scoped manifest updates
  - SCORCH submit path now forwards `--ph` to avoid cross-combo duplicate rescoring
- Routed pocket-detection event consolidation through global admission (`acquire_global_slot`) as a low-core housekeeping task.

## 2026-02-17

- Added `-bench-small/--bench-small` profile:
  - keeps `bench2` unchanged
  - uses `test_library_20` for the standard benchmark PDB set
  - runs Vina + SCORCH only with holo + pH ensemble.
- Updated distributed combo-chunk mode behavior in `main.py`:
  - combo-chunk assignment is now the default for multi-task distributed runs (`slurm_array` task_count > 1)
  - legacy per-PDB sharding is available as an explicit opt-out via `ATLAS_DISABLE_COMBO_CHUNKS=1`
  - worker assignment is now chunk-plan based (LPT weighted initial assignment)
  - workers claim chunks via atomic marker files and can steal unclaimed work
  - chunk plans are persisted under `manifests/<run_id>/distributed/`.
- Enabled combo chunk splitting for large distributed runs (not only small-profile capped runs).
- Kept run-history-informed chunk weights with deterministic fallback and deterministic small-profile ligand cap.
