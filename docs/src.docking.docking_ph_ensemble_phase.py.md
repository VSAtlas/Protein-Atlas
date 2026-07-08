# Notes: `src/docking/docking_ph_ensemble_phase.py`

- Removed eager `detect_active_site(...)` call from `_phase5_ph_ensemble_global(...)`.
- When no precomputed center is available, pH ensemble now uses global-default center/radius directly.
- This prevents pre-control active-site/P2Rank work from being triggered during pH ensemble setup.
- Added `_PH_TAG_OVERRIDE` handling:
  - combo-level scheduling can request a specific pH label
  - pH ensemble generation uses override-derived pH values for that combo instead of broad context enumeration.
