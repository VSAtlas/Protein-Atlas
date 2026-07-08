# run_manifest.py change note

## 2026-02-15

- Added distributed per-protein state writing under `manifests/<run_id>/distributed/proteins/*.json` for multi-node array workers.
- Added merge helper `reduce_distributed_protein_states(...)` to consolidate distributed state files into the canonical `run_manifest.yaml`.
- Updated manifest write paths for protein/stage/pocket/druggability updates to use per-protein distributed state when `ATLAS_DISTRIBUTED_MODE=slurm_array` is active.
- Updated finalize path to merge distributed state before marking run completion.

## 2026-02-16

- Added `reconcile_distributed_manifest_state(...)`:
  - detects distributed state presence/staleness
  - reduces shard state into `run_manifest.yaml` when stale or forced
  - reports reconciliation diagnostics (`action`, `state_files`, `merged_entries`, stale flags).
- `load_run_manifest(...)` now performs best-effort reconcile-on-read, so interrupted distributed runs can recover canonical manifest state without manual edits.
- `finalize_run_manifest(...)` now reconciles distributed shard state via the same helper before writing terminal run status.
- Improved reducer metadata stability:
  - avoids downgrading `command.DISTRIBUTED_TASK_COUNT` to `1` when reducer runs outside array env
  - preserves existing count or infers from observed shard task IDs when env hints are absent.
- Reducer writes are now no-op aware (`_write_manifest` only when merged manifest content actually changed).

## 2026-02-16 (distributed pH base-key suppression)

- Added distributed write gating so distributed manifest mutations do not emit `...|base` shard entries when:
  - distributed mode is active (`slurm_array`), and
  - `PH_ENSEMBLE=true`, and
  - incoming manifest update has `ph_tag=None`.
- pH-tagged entries (for example `...|pH7_0`) are still written normally.
- Non-pH distributed runs keep existing base-key behavior for backward compatibility.
