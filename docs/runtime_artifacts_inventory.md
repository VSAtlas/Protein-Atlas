## Runtime Artifacts Inventory (Phase 4)

Goal: keep the repository reviewable by separating generated runtime products from intentional fixtures.

## Classification

| Path | Classification | Policy |
| --- | --- | --- |
| `outputs/docked/` | runtime output | ignore from VCS; regenerate per run |
| `outputs/post_docked/` | runtime output | ignore from VCS; regenerate per run |
| `prepped_ligands/` | runtime output | ignore from VCS; bundled fixtures live under `chemdb/tests/fixtures/prepped_ligands/` |
| `extracted_ligands/` | runtime output | ignore from VCS; bundled fixtures live under `chemdb/tests/fixtures/extracted_ligands/` |
| `input_pdbs/` | user input working directory | ignore from VCS; bundled `TEST.pdb` fixture lives under `chemdb/tests/fixtures/input_pdbs/` |
| `input_ligands/` | user input working directory | ignore from VCS; do not review local FDA/library dumps |
| `outputs/processed_pdbs/<run_id>/` | runtime output | ignore from VCS |
| `outputs/logs/` | runtime output | ignore from VCS |
| `outputs/configs/<run_id>/` | runtime output | ignore from VCS |
| `outputs/manifests/<run_id>/` | runtime output | ignore from VCS |
| `outputs/data/<run_id>/` | runtime output | ignore from VCS |
| `chemdb/bench/runs/` | mixed (historical benchmark corpus + runtime) | quarantine; do not add new large run outputs unless explicitly documented |
| `*.stdout.json` | runtime output | ignore from VCS |

## Current Audit

This cleanup pass audited tracked generated paths with:

```bash
git ls-files outputs docked post_docked prepped_ligands processed_pdbs logs chemdb/bench/runs data configs '*.stdout.json'
```

Findings:

| Path family | Tracked status | Decision |
| --- | --- | --- |
| `outputs/`, legacy `docked/`, legacy `post_docked/`, legacy `processed_pdbs/`, legacy `logs/`, legacy `configs/`, legacy `data/` | no intentional source files in this audit scope | keep ignored; do not review runtime products in normal PRs |
| `chemdb/tests/fixtures/prepped_ligands/fda_test_library_10/` | small smoke fixture | canonical checked-in prepped ligand fixture for fast acceptance tests |
| `chemdb/tests/fixtures/prepped_ligands/test_library_10/` | manifest-only legacy fixture record | retained for tests that still exercise DUD library naming without relying on runtime roots |
| `chemdb/tests/fixtures/extracted_ligands/` | small source-ligand fixtures | canonical checked-in extracted ligand fixture root |
| `chemdb/tests/fixtures/input_pdbs/TEST.pdb` | small receptor smoke fixture | canonical checked-in input PDB fixture for `--test` acceptance paths |
| `chemdb/bench/runs/` | generated benchmark outputs were tracked | untrack generated run products; regenerate with `chemdb/bench/util_bench/run_bench.py` when needed |

The preferred final fixture locations are `chemdb/tests/fixtures/`, `analysis/tests/fixtures/`, and `docs/examples/`.

## Fixture policy

Small deterministic fixtures belong under:

- `chemdb/tests/fixtures/`
- `analysis/tests/fixtures/`
- `docs/examples/`

Rules:

- keep fixture files small;
- document fixture origin and intent in nearby README or test docstring;
- avoid storing machine-specific absolute paths in fixture metadata.
