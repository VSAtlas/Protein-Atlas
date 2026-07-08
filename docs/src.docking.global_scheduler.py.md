# Notes: `src/docking/global_scheduler.py`

- Added tail-mode detection to shift admission toward tighter fit and shorter estimated jobs near run tail.
- Strengthened large-request elastic capping (`want_cores >= 5`) under pressure/tail so oversized requests shrink earlier.
- Preserved starvation protection via aging while changing tie-break order for tail packing.
- Tightened tail-mode activation criteria so tail heuristics only engage during true endgame conditions (small queue + meaningful stranded capacity + poor immediate fit mix).
- Added a critical-tail mode for straggler control:
  - triggers only in narrow late-run conditions (small queue, high free-core fraction, large chunks still pending)
  - prioritizes aged large fitting requests before filler tasks
  - reduces long-tail stranding without changing normal mixed-load ordering.
