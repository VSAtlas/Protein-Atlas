# Testing

Atlas2 uses pytest markers to make common verification slices explicit. The
default suite should not silently skip passing tests unless they require
unavailable external services or tools.

## Marker Reference

- `fast`: lightweight unit or documentation checks.
- `integration`: end-to-end or multi-module pipeline checks.
- `external_tool`: wrappers or behavior around external binaries.
- `report`: generated report, heatmap, provenance, and reporting metadata tests.
- `slow`: full-run or long-running integration checks.
- `network`: tests that require external network access.

## Common Suites

Smoke suite:

```bash
pytest chemdb/tests/test_path_router.py -q
pytest chemdb/tests/test_ligand_run_modes.py -q
python main.py --test -fast
```

Report suite:

```bash
pytest -m report chemdb/tests -q
```

Pathing suite:

```bash
pytest chemdb/tests/test_path_router.py -q
pytest chemdb/tests/test_config_inline_comments_parsing.py -q
```

Full local quality gate:

```bash
atlas dev verify --full --fix --smoke
```

External-tool-dependent suite:

```bash
pytest -m external_tool chemdb/tests -q
```

Use `atlas --doctor` before long runs when environment activation is unclear.
