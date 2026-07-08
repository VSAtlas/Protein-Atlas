# analysis/dud_eval_orchestrate.py import-fix note (2026-02-13)

## Problem
- Orchestration layer imported `expand_variants` / `make_paths` via `path_router.path_router`.
- This can fail when `src/` is not pre-injected into Python search paths.

## Change
- Switched to `from src.path_router.path_router import expand_variants, make_paths`.

## Why
- Keeps import behavior consistent with package layout and prevents module-resolution regressions in report/evaluation subprocesses.
