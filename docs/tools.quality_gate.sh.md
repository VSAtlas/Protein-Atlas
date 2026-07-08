# tools/quality_gate.sh

## 2026-03-02

- Added a single-command quality gate for maintained pipeline modules.
- Default checks:
  - `ruff check` on core orchestration/runtime/tooling paths.
  - `mypy --follow-imports skip` on distributed runtime/planner and benchmark tools.
- Optional flags:
  - `--fix` to apply Ruff autofixes on the scoped targets.
  - `--smoke` to append `python main.py --test -fast` after lint/type checks.
- Purpose: keep a reproducible, fast professionalism gate that avoids noisy legacy/vendor areas while preserving existing runtime behavior.
