# analysis/dud_eval_discovery.py import-fix note (2026-02-13)

## Problem
- `analysis.run_report` failed at runtime with `ModuleNotFoundError: No module named 'path_router'`.
- Root cause was reliance on top-level imports (`path_router.*`, `prep_ligands.*`) that require `src/` to already be on `PYTHONPATH`.

## Change
- Updated imports in `analysis/dud_eval_discovery.py` to explicit package-qualified forms:
  - `src.path_router.path_router`
  - `src.prep_ligands.library_index`

## Why
- `analysis` tools can be launched as standalone modules/subprocesses during post-run hooks.
- Explicit `src.*` imports avoid environment-dependent import failures and remove dependence on path injection side effects.
