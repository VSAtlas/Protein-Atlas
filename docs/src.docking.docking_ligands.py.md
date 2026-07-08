# `src/docking/docking_ligands.py` bugfix notes

## 2026-02-15 — mapped library undercount / random ligand subset

- Issue: In manifest-backed enumeration, non-control candidates were sourced from `manifest.entries`, which is lookup-oriented and can under-represent full inventory in stale/partial states.
- Change:
  - Prefer `manifest.filenames` as the canonical full ligand inventory.
  - Keep `entries` only as compatibility fallback when `filenames` is unavailable.
  - Add per-root undercount guard: if manifest count is lower than filesystem-discovered `.pdbqt` count, augment from filesystem and log `[lib-index.manifest.undercount]`.
  - Sort deduped candidates deterministically by normalized path.
  - Add `[ligands.coverage]` log line with roots/candidate/filtered/control counts.
- Effect: Stage-1 ligand attempt pool is deterministic and no longer vulnerable to partial-manifest under-enumeration in mapped DUD/bench runs.

## 2026-02-17 — distributed chunk ligand selection hook

- Added `_CHUNK_LIGAND_KEYS` filtering for non-control ligands in `prepare_and_filter_ligands(...)`.
- When chunk keys are provided, only matching non-control ligands are admitted; controls remain untouched.
- Added small-profile deterministic cap path for non-control ligands (used by bench-small combo-chunk planning fallback).
- Purpose: allow distributed combo-chunk workers to run disjoint ligand subsets without large pipeline refactors.
