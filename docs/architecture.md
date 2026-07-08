## Atlas2 Architecture

Maintained owners:

- CLI and runtime orchestration: `src/cli/`
- Path routing: `src/path_router/`
- Docking engines/orchestration: `src/docking/`
- Protein preparation: `src/protein_prep/`
- Ligand preparation: `src/prep_ligands/`
- Post-docking and rescoring: `src/post_docking/`
- Reporting and analysis: `analysis/reporting/`, `analysis/dud_eval_core/`
- Chemical annotations and mappings: `chemdb/`

Root-level policy:

- `atlas` is the public CLI.
- `main.py` is compatibility-only and should not accumulate new feature logic.
- Root wrappers delegate to owner modules.

For deeper analysis implementation details, see `analysis/ARCHITECTURE.md`.
