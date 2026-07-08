# Notes: `src/docking/docking_receptor_phases.py`

- Removed eager active-site detection from `_phase2_to4_receptor_and_center(...)`.
- Active-site/P2Rank is now explicitly deferred to later fallback handling after controls-first centering attempts.
- This avoids running pocket fallback logic for proteins that can be centered from controls.
