## Atlas2 Baseline Snapshot (Phase 0)

Date: 2026-04-28
Goal: Capture pre-cleanup behavior before structural cleanup.

### Required command results

- `tools/quality_gate.sh`: **failed**
  - Ruff reported missing files in target list:
    - `postrun_hooks.py`
    - `run_manifest.py`
- `tools/quality_gate.sh --smoke`: **failed** (same missing-file errors before smoke run)
- `pytest chemdb/tests/test_ligand_run_modes.py -q`: **passed** (`2 passed`)
- `pytest chemdb/tests/test_path_router.py -q`: **passed** (`27 passed`)
- `pytest chemdb/tests/test_record_data_csv.py -q`: **passed** (`10 passed`)
- `pytest chemdb/tests/test_main_full_run.py -q`: **in progress beyond 10 min timeout**
  - Command continued in background after timeout; no final pytest summary captured yet.

### Current public CLI surface (observed)

- Packaging entrypoint (`pyproject.toml`) maps:
  - `atlas` -> `cli.atlas_main_cli:main`
  - `atlas-main` -> `cli.atlas_main_cli:main`
- `src/cli/atlas_main_cli.py` currently delegates to `main.main()`.
- Root script `python main.py` remains operational and appears to hold primary orchestration.

### Known baseline issues to preserve/document (not silently ignored)

- Quality gate target list references missing root files (`postrun_hooks.py`, `run_manifest.py`).
- Main full-run acceptance test runtime currently exceeds the 10-minute command timeout used in this baseline capture.

### Supported vs unsupported (baseline statement)

- Supported user-facing entrypoint in README: `atlas [OPTIONS] [PDB_IDS...]`.
- `python main.py` remains documented for developer compatibility.
- Legacy one-off scripts are present and not yet fully classified in a dedicated inventory doc.
