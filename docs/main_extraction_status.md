## Main Extraction Status (Phase 6)

Behavior-preserving seam extraction policy:

- one seam at a time;
- regression tests with each seam;
- no output path/schema drift without explicit migration.

## Current status

`main.py` already delegates substantial orchestration to owner modules in `src/cli/`, including:

- run bootstrap/context (`run_bootstrap.py`, `run_context.py`)
- run lifecycle/finalization (`run_lifecycle.py`)
- distributed chunk planning/runtime (`distributed_*`)
- per-protein process wrapper wiring (`run_process_one.py`, `run_execution.py`)

## Remaining rule

Keep reducing root `main.py` only via tested seam moves.
No new behavior should be introduced as part of seam movement.
