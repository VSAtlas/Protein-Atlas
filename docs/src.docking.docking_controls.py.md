# Notes: `src/docking/docking_controls.py`

- Control redock serial path now acquires explicit granted cores via `acquire_global_cores(...)` instead of slot-only admission.
- Vina control redock invocations now pass `cpu_override=<granted>` to align runtime thread count with global admission.
- Restored controls-first reliability by relaxing control reference lookup fallback:
  - Keep non-empty `.pdb` refs as last-resort when RDKit readability checks fail for `.sdf/.mol2`.
  - Candidate control PDBQTs are now admitted when either `control_lookup` **or** extracted `ctrl_pdb_map` has a reference for that base.
- Redock center and RMSD reference resolution now use a unified fallback (`control_lookup` first, then extracted crystal PDB), reducing false `fallback_no_controls` outcomes.
- Removed eager active-site fallback from `_resolve_ph_scope(...)` for pocket-scoped pH selection:
  - If no control centroid is available yet, scope now falls back to global directly instead of calling `detect_active_site` early.
