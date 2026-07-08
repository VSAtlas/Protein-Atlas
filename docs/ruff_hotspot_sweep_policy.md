# Ruff Hotspot Sweep Policy

This repository should not receive full-repo Ruff rewrites in a single patch.

## Review policy
- Prefer small, subsystem-scoped lint patches.
- Avoid behavior-changing refactors in lint-only changes.
- If import-order (`E402`) is intentionally tied to bootstrap/runtime setup, use a file-level explicit Ruff waiver rather than risky import movement.

## Execution policy
- Run Ruff on target files first.
- Keep non-target Ruff debt for later bounded sweeps.
- Pair each sweep with focused smoke tests for touched modules.
