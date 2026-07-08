# Notes: `chemdb/bench/util_bench/run_bench.py`

- Backfill admission now performs size-aware best-fit + short-duration tie-breaking instead of first-fit only.
- Synthetic benchmark tasks now emit explicit `want_cores`, `min_cores`, and `granted_cores` in task lifecycle events.
- Added elastic handling for large synthetic tasks (`cores >= 5`) in the runner:
  - admissible with reduced grant (min-core floor),
  - command execution now uses a stable granted-worker environment envelope (`UTIL_BENCH_GRANTED_WORKERS`) instead of brittle regex rewriting.
