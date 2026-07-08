# src/cli/distributed_context.py

## 2026-02-15

- Added distributed runtime context helpers for Slurm-array mode:
  - scheduler allocation CPU detection
  - deterministic per-PDB sharding
  - per-task completion marker write/read
  - leader barrier wait and failure aggregation helpers

## 2026-02-17

- Added combo-chunk filesystem primitives for distributed execution:
  - chunk plan path helpers
  - atomic JSON writer helper
  - chunk claim/result path helpers
  - `try_claim_chunk(...)` with stale-lease reclaim support
  - `write_chunk_result(...)` and result-exists helper.
- Kept existing per-PDB sharding behavior intact for runs not using combo-chunk mode.
