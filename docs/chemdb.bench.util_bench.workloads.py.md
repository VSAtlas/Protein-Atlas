# `chemdb/bench/util_bench/workloads.py` change note

## Context
`test_stage_collapse_workload_shape` failed because the generated stage-3 count was lower than expected.

## Root cause
`generate_stage_collapse_workload()` computed `stage3` as a fraction of `stage2`, which compounds rounding and can collapse too aggressively.

## Fix
`stage3_keep_frac` is now applied against `stage1` directly:
- before: `s3 = round(s2 * stage3_keep_frac)`
- after: `s3 = round(s1 * stage3_keep_frac)`

This keeps stage-shape behavior stable and aligned with the regression test expectation.

## Additional fix: synthetic command validity + core realism
The synthetic task command previously failed with `SyntaxError` and each task used only one core regardless of `task.cores`, which made utilization numbers misleading.

Updated `_busy_loop_command(duration, cores)` to use `multiprocessing` and start `cores` worker processes running a timed busy loop script. This makes task core requests match actual CPU pressure and produces meaningful utilization deltas.

## Elastic worker envelope compatibility
Synthetic commands now read `UTIL_BENCH_GRANTED_WORKERS` (when set) to select runtime worker count. This lets the bench runner enforce granted-core elasticity without command-text rewriting.
