## Contributing

Thanks for contributing to Atlas2.

## Ground rules

- Preserve behavior by default; avoid docking/scoring logic changes unless backed by regression tests.
- Keep new maintained implementation in owner modules under `src/`, `analysis/`, or `chemdb/`.
- Treat root scripts as compatibility wrappers unless explicitly promoted.

## Local checks before opening a PR

```bash
tools/quality_gate.sh
pytest chemdb/tests/test_main_full_run.py -q
pytest chemdb/tests/test_ligand_run_modes.py -q
pytest chemdb/tests/test_path_router.py -q
pytest chemdb/tests/test_record_data_csv.py -q
```

## Documentation updates

Please update docs when changing:

- public CLI/config/output surface (`docs/public_api.md`)
- ownership boundaries (`docs/module_ownership.md`)
- schemas (`docs/schemas/*.md`)
- reproducibility assumptions (`docs/reproducibility.md`)

## Pull request guidance

- keep changes narrow and reviewable;
- include rationale and verification notes;
- avoid mixing large refactors with behavior changes.
