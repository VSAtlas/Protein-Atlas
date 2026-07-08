# throughput_acceptance.py

## Summary

Added `analysis.reporting.throughput_acceptance`. For user-facing runs, prefer
the Atlas wrapper:

```bash
atlas analysis throughput acceptance <BASELINE_RUN_ID> <CANDIDATE_RUN_ID>
```

## What it enforces

1. Hard integrity gate (no skipped completed-combo ligands):
   - expected resolved
   - expected == docking_scored == post_scored
   - no known-failed masking and no unexplained missing counts
2. Wall-time improvement target (default 10%).
3. Optional util-bench targets (idle-core, idle-frac, fragmentation+dispatch improvements).

## Typical use

Run after baseline and candidate runs to generate `data/<candidate_run_id>/throughput_acceptance.json` and fail fast on acceptance regressions.
