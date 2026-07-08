# artifact_retention.py Bug Note (2026-02-15)

## Issue

Running:

```bash
python3 -m src.post_docking.artifact_retention ...
```

failed in environments that do not have the full `config` module wiring on `PYTHONPATH`:

- `ModuleNotFoundError: No module named 'config'`

The failure happened at module import time via `input_and_export_functions`.

## Fix

- Switched artifact retention config loading to a lazy/best-effort helper.
- If importing `load_inputs()` fails, retention now falls back to `{}` and continues.
- This preserves CLI-driven operation (`--docked-root`, `--post-docked-root`, `--run-id`) without hard-failing.

## Why this is safe

- Root resolution already supports explicit CLI overrides.
- Missing optional config now degrades to defaults instead of aborting.
- When config is available, behavior is unchanged.
