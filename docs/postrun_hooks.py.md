# postrun_hooks.py

## 2026-02-16: Relocated Analysis Root Resolution for Automatic Hooks

- **Issue:** Post-run analysis hooks used the code repository root as `--repo-root` and DUD eval defaulted to local `docked/`, which breaks automatic analysis when runs are relocated (for example via `tools/run_relocated_mode.py --all-dirs ...`).
- **Fix:** Added centralized hook root resolution that derives `analysis_root`, `docked_root`, and `post_docked_root` from config/env (`OVERALL_DIR`, `DOCKED_DIR`, `POST_DOCKED_DIR`) with safe fallbacks.
- **Behavior change:** Automatic hooks now pass relocated-aware roots to:
  - `analysis.cli.dud_eval` (`--docked-root`, `--post-docked-root`, relocated `--out-dir`)
  - `analysis.cli.run_report` (`--repo-root <relocated_root>`)
  - `analysis.cli.master_schema_export` (`--repo-root <relocated_root>`)
- **Result:** Relocated runs no longer require manual report/schema invocations to point at relocated directories.

## 2026-02-14: Pathing Fix

- Fixed DUD post-run output location from legacy locations to canonical `analysis/dud_eval`.
- Updated `_maybe_run_dud_eval(...)` command args:
  - `--out-dir analysis/dud_eval`
- Updated completion log message to report `out_root=analysis/dud_eval`.

## Compatibility

- Legacy data under `dud_eval` and `docked/dud_eval` remains discoverable by reporting fallback search paths.
