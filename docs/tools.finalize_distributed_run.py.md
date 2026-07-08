# tools/finalize_distributed_run.py

## 2026-02-16

- Added a dedicated distributed finalizer entrypoint for dependency-job execution after array workers complete.
- The tool:
  - loads run config snapshot for `RUN_ID`
  - merges distributed manifest shards
  - aggregates worker failure markers
  - finalizes `run_manifest.yaml`
  - runs run-level post hooks (unless `--no-hooks`).
- Added explicit repo/src path bootstrap at script startup so `python tools/finalize_distributed_run.py ...` works reliably when launched from `tools/` (for example via Slurm scripts).
- Hook execution is best-effort (exceptions are logged and skipped), matching main-run behavior so finalization is durable.

## 2026-02-16 (reconciliation hardening)

- Finalizer now logs distributed input diagnostics before finalize:
  - marker file/payload counts
  - distributed per-protein shard-state file count
  - present task IDs, expected task count (from marker payloads), and marker gaps.
- Switched from direct reducer call to `reconcile_distributed_manifest_state(..., force=True)` to guarantee manifest reconciliation even when canonical YAML is stale/missing.
- Added `--reconcile-only` mode for interrupted runs:
  - merges distributed shard state into `run_manifest.yaml`
  - skips run-status finalization and all post hooks.
