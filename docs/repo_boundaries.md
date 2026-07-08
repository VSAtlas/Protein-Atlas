# Repo Boundaries

## Maintained First-Party Code
- `src/`
- `analysis/`
- `chemdb/`
- thin root entrypoints such as `main.py`, `run_manifest.py`, and `postrun_hooks.py`

New implementation code should live under an owning package, primarily in `src/`.
Root-level modules should stay as wrappers or compatibility entrypoints unless there
is a strong reason to keep library code at the repo root.
For maintained command-line workflows, prefer packaging entry points via
`[project.scripts]` that target owner modules under `src/cli`.

## Quarantined Legacy Or Vendor Code
- `molprobity_coot.py`
- `propka_wire.py`
- `benchmark_mode.py`
- `benchmark_auto_analysis.py`
- `identify_ligand.py`
- other legacy one-off root scripts that are not part of the maintained runtime

These paths are kept for compatibility, benchmarking, or integration purposes, but
they are not the default place for new work.

Current policy for quarantined root scripts:
- they may remain import-compatible or runnable for legacy workflows
- they do not block the maintained-code quality bar
- if a maintained module needs functionality from them, that functionality should be
  moved behind a maintained owner module instead of importing the script directly

## Generated Or Runtime Data
- `analysis/_tmp/`
- `outputs/docked/`
- `outputs/logs/`
- `node_modules/`
- `outputs/post_docked/`
- `prepped_ligands/`
- `outputs/processed_pdbs/`
- `input_pdbs/`

These directories are runtime outputs or inputs, not source code. They should stay
out of routine linting, typing, and review diffs.
Do not track host-specific pointer files for these directories (for example a
root `node_modules` path pointer). Keep dependency ownership in lock/config
files under an owning tool directory, such as
`tools/report_assets/npm/package.json`, instead.

## Engineering Rules
- No new maintained module should grow beyond `1000` lines. Once a file crosses that
  threshold, new behavior should go into a new owner module.
- Maintained runtime code should import from owner modules, not compatibility facades,
  when the owner module already exists.
- Refactors should finish ownership moves in the same patch. Avoid half-migrated states
  where the new owner exists but callers still treat the old facade as canonical.
- Quality gates should be credible:
  - `ruff check --fix .` should pass for maintained code
  - `mypy` should pass for changed maintained modules
  - vendor, generated, and protected-test paths should be explicitly excluded rather
    than silently tolerated
