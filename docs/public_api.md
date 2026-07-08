## Public API Surface

This document defines the publication-facing API for Atlas2.

Behavior-preserving rule: this file describes existing supported behavior; it does not authorize new behavior.

## Supported CLI entrypoints

Official user-facing command:

- `atlas [OPTIONS] [PDB_IDS...]`

Compatibility/developer commands:

- `python main.py [OPTIONS] [PDB_IDS...]`
- `atlas-main [OPTIONS] [PDB_IDS...]`

Specialized developer utilities installed from `pyproject.toml`:

- `atlas-cap-driver`: cap/driver utility entrypoint.
- `atlas-phenix-clean`: Phenix cleaning utility entrypoint.
- `atlas-generate-skeletons`: documentation/skeleton generation utility.

These utilities are supported only for their documented purpose. The general docking workflow remains `atlas`.

## Supported CLI workflows

The stable public workflows are:

- Print help: `atlas --help`
- Guided local setup/config detection: `atlas init`
- Generate no-docking demo report: `atlas demo`
- Run the publication-safe smoke: `atlas smoke public`
- Plan or launch the first target/library run: `atlas first-run --uniprot <ACCESSION> --ligands <library>`
- Summarize a run: `atlas status [RUN_ID]`
- Generate report artifacts with master export: `atlas report <RUN_ID>`
- Generate analysis/report/status artifacts together: `atlas analysis report <RUN_ID> --status-html`
- Write the reproducibility artifact index under the run data directory: `atlas reproduce bundle <RUN_ID>`
- Debug run failures through structured summaries: `atlas debug <RUN_ID> --deep`
- Submit/track/finalize Slurm array runs: `atlas slurm submit|progress|finalize ...`
- Render pose screenshots: `atlas screenshot <RUN_ID> --pdb <ID> --top 20`; default is full-protein cartoon context, with `--view-context pocket` and `--protein-style surface|mesh|lines|cartoon|none` for alternate receptor context.
- Render a pose review gallery: `atlas screenshot <RUN_ID> --pdb <ID> --representatives --gallery-html`; writes representative best/median/worst PNG cards to a self-contained `gallery.html`.
- Query target structures from genes: `atlas targets from-genes EGFR ABL1 --out <CSV>`
- Query target structures by topic text: `atlas targets search "BRCA DNA repair" --out <CSV>`
- List built-in target panels: `atlas targets panels`
- Select and install target PDBs: `atlas targets install --panel kinases`
- One-command panel plus ligand workflow: `atlas run-panel kinases --ligands chembl --fast`
- Preflight a target/ligand combo: `atlas doctor --pdb 1ABC --ligands fda`
- Verify external tool resolution: `atlas --verify-tools`
- Print resolved configuration: `atlas --print-effective-config`
- Run selected PDB IDs: `atlas --pdb 1BN1 --pdb 2OJ9 --fast`
- Run bare PDB ID arguments: `atlas 1BN1 2OJ9 --fast`
- Run all configured/input PDBs: `atlas`
- Run the internal smoke mode: `atlas --test -fast`
- Run the public smoke fixture: `bash tools/public_smoke_check.sh`
- Run the small FDA library mode: `atlas -test-fda --fast`
- Resume a run: `atlas --resume --run-id <RUN_ID>`
- Rebuild ligand manifests: `atlas -rebuild`

Benchmark, DUD, pH, retention, and distributed modes are maintained developer/publication workflows when invoked through `atlas` and covered by current tests. They are not a license to call legacy helper scripts directly.

## Supported analysis/reporting commands

Maintained reporting and analysis code lives under `analysis/` and is normally invoked through documented CLI wrappers. Current supported commands include:

- `atlas analysis report <RUN_ID> --status-html` for master export, report generation, and optional status HTML.
- `atlas analysis export-master <RUN_ID>` for canonical table export when report generation is not needed.
- `atlas analysis dud-eval <RUN_ID>` for DUD/benchmark evaluation.
- `atlas analysis interactions export <RUN_ID>` for interaction export conversion.
- `atlas analysis throughput integrity <RUN_ID>` and `atlas analysis throughput acceptance <BASELINE_RUN_ID> <CANDIDATE_RUN_ID>` for throughput checks.
- `atlas debug <RUN_ID> [--deep]` for summarized run debugging.

For the manuscript-oriented side-effect/ADR workflow, see `docs/atlas_analysis_pipeline.md`.

## Supported config templates and keys

Primary supported configuration is `config.txt`, usually copied from `config.example.txt`:

```bash
cp config.example.txt config.txt
```

`config.example.txt` is the portable template, and `config.full.example.txt` is
the expanded optional-key reference. They stay at the repository root so the
documented `cp config.example.txt config.txt` workflow is obvious. Local
`config.txt` files and config backups are machine-specific and ignored.

Supported key families:

- Repository/run roots: `OVERALL_DIR` and configured input/output roots.
- External tool roots and prefixes: `ADFRSUITE_BIN`, `PHENIX_DIR`, `PHENIX_LIB_PATH`, `P2RANK_PATH`, `MMGBSA_AMBERTOOLS_PREFIX`, `AMBERTOOLS_PREFIX`.
- Tool verification: `TOOL_VERIFY_ON_START`.
- SCORCH environment/device selection: `SCORCH_DEVICE`, `SCORCH_GPU_IDS`,
  `SCORCH_ENV_PREFIX`, `SCORCH_ENV`.
- Optional PDB query defaults: `PDB_QUERY_GENES_FILE`.
- Library and test-mode selection: `LIBRARY_SUBDIR_DEFAULT`, `HMDB_LIBRARY_SUBDIR`, `TEST_MODE_ENABLE`, `TEST_LIBRARY_MAP`.
- Downloaded ligand-source cache: `LIGAND_SOURCE_CACHE_DIR`.
- Runtime controls used by the CLI: run ID, resume mode, fast mode, DUD library override, distributed mode, pH mode, retention mode, and scheduler-related knobs.

Publication policy: newly supported config keys must be added to `config.example.txt` and documented here or in README-linked docs.

## Expected input layout

Minimal expected layout for supported runs:

- `config.txt` at the repository root, copied from `config.example.txt`.
- Input proteins under the configured input root; the default project layout uses `input_pdbs/`.
- Ligand libraries under the configured ligand roots.
- Optional benchmark or DUD libraries under their configured library roots.
- Bundled smoke fixtures live under `chemdb/tests/fixtures/`, not under root runtime directories.

Ligand library helpers:

- List supported downloadable libraries: `atlas ligands sources`
- Download and prepare a library: `atlas ligands install chembl`
- Plan library install without fetching or preparing: `atlas ligands install chembl --dry-run`
- Other built-in libraries: `chebi`, `coconut`, `fda`, `hmdb`
- Prepare a manually downloaded HMDB archive: `atlas ligands install hmdb --source-sdf <path>`
- One-command install-if-needed plus docking: `atlas chembl --pdb <ID> --fast`

Target install helpers:

- Guided target setup: `atlas targets guide`
- Write a reusable gene CSV: `atlas targets genes --genes EGFR ABL1 --gene-type KINASE`
- List built-in target panels: `atlas targets panels`
- List organism aliases: `atlas targets species`
- Select/download panel targets: `atlas targets install --panel kinases`
- Select/download custom gene targets: `atlas targets install EGFR ABL1`
- Select/download non-human targets: `atlas targets install Trp53 --species mouse`
- Select/download direct UniProt targets: `atlas targets install --uniprot P04637`
- Script-compatible custom gene targets: `atlas targets install --genes EGFR ABL1`
- Filter custom gene searches: `atlas targets install --genes EGFR --species "Homo sapiens" --taxonomy-id 9606 --methods X-RAY --resolution-max 2.8 --gene-type KINASE`
- Require co-crystal ligands and query quality controls: `atlas targets install EGFR --ligand ATP --dedupe sequence-identity 90 --quality publication`
- Disable default sequence de-duplication when needed: `atlas targets install EGFR --no-dedupe`
- Include apo structures when needed: `atlas targets from-genes EGFR --allow-apo`
- Dry-run without downloads: `atlas targets install --panel gpcrs --dry-run`
- Run a complete panel/library workflow: `atlas run-panel kinases --ligands chembl --fast`
- Validate prepared receptor/ligand inputs before docking: `atlas doctor --pdb <ID> --ligands <library>`

Do not treat local runtime data as source. Inputs such as `input_pdbs/` may be required to run Atlas locally, but they are protected working data and are not routine code-review targets.

## Supported output directories and files

Supported run outputs are run-scoped artifacts under `outputs/` by default, including:

- `outputs/manifests/<RUN_ID>/run_manifest.yaml`: canonical run audit manifest.
- `outputs/configs/<RUN_ID>/`: per-run configuration snapshots.
- `outputs/logs/main_<RUN_ID>.log` or distributed variant names: top-level run logs.
- `outputs/processed_pdbs/<RUN_ID>/<PDB>/`: cleaned and prepared receptor structures.
- `prepped_ligands/`: prepared ligand libraries.
- `outputs/docked/<RUN_ID>/<PDB>/`: docking scores and poses.
- `outputs/post_docked/<RUN_ID>/<PDB>/`: post-docking analysis, reranking, and reports.
- `outputs/data/<RUN_ID>/`: reports, heatmap inputs, master rows, status HTML, and run summaries.

Schema contracts promoted for publication live under `docs/schemas/`. Changing those schemas should be reviewed as a scientific-method change.

## Unsupported or quarantined paths

Legacy utility scripts outside the documented CLI are not part of the public API unless explicitly listed in README and tested in acceptance workflows. Current quarantined examples are tracked in `docs/legacy_inventory.md`.

The following directories are runtime/protected data rather than source API:

- `input_pdbs/`
- `outputs/docked/`
- `outputs/post_docked/`
- `prepped_ligands/`
- `outputs/processed_pdbs/`
- `outputs/logs/`

## Compatibility wrapper policy

Root-level and bridge entrypoints may remain for import or script compatibility. For those files:

- Compatibility wrapper only.
- No new logic here.
- Use the canonical owner module for new code.
- New implementation code goes to owner modules under `src/`, `analysis/`, or `chemdb/`.
