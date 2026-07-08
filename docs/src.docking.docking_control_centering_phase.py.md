# Notes: `src/docking/docking_control_centering_phase.py`

- `_phase5b_controls_and_control_redock(...)` now performs a controls-only centering pass first (no fallback center provided).
- Active-site fallback (`detect_active_site`) is invoked only when controls fail to produce any per-pH center.
- After fallback center resolution, centering is retried with fallback enabled to preserve existing downstream behavior.
