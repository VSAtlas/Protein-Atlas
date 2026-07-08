# Notes: `src/docking/docking_control_redock.py`

- Control docking loop now uses core-grant admission (`acquire_global_cores`) and propagates granted cores into Vina execution.
- This prevents oversized control-redock runs from consuming more CPU threads than the scheduler granted.
- Added pH override support via `_PH_TAG_OVERRIDE`:
  - per-combo runs can restrict control-center/control-redock work to a single pH label
  - avoids redundant cross-pH control docking when combos are scheduled independently.
