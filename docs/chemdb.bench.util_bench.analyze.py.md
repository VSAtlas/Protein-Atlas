# Notes: `chemdb/bench/util_bench/analyze.py`

- Analyzer now models pending tasks with `want_cores`, `min_cores`, and `granted_cores` (`TaskShape`) instead of only `cores`.
- Fragmentation attribution now prefers `min_cores` from events, with fallback to legacy `cores`.
- Added event-field coverage counters to `summary.json`/`summary.csv`:
  - `task_events_total`
  - `task_events_with_min_cores`
  - `task_starts_total`
  - `task_starts_with_granted_cores`
- Added hedge diagnostics in summaries:
  - `hedge_started`
  - `hedge_won`
  - `hedge_canceled`
  - `hedge_extra_core_seconds_est`
- Purpose: ensure future regressions can distinguish scheduler packing behavior from missing instrumentation.
- Updated idle attribution to be elastic-aware: fragmentation now considers `want_cores`/`granted_cores` pressure instead of relying only on `free_cores < min_runnable`.
