# Notes: `src/docking/docking_ph_subruns.py`

- Added `_PH_TAG_OVERRIDE` support in `resolve_ph_tags_and_root(...)`.
- When override is set, subrun pH iteration is constrained to that single tag.
- This enables per-(PDB,variant,pH) combo scheduling while preserving existing run-mode and library-root resolution.
