# Atlas Protein Automation

Atlas is a Python docking pipeline for preparing protein structures, selecting
active-site, docking ligand libraries, rescoring/reranking poses, with MMGBSA and analysis-ready outputs.
The maintained workflow is the `atlas` CLI and its owner modules under `src/`.

The public entrypoint is:
```bash
atlas [OPTIONS] [PDB_IDS...]
``
## Requirements

- Linux HPC environment
- Python 3.10
- Conda/micromamba environment from `environment.core.yml` (extras in `environment.optional.yml`)
- Core external tools configured in `config.txt` or available on `PATH`
- Input PDB files under `input_pdbs/`

One-command install:

```bash
bash tools/installers/install_atlas_publication_stack.sh
```
This installs the environment, installs Atlas, creates
and auto-populates local `config.txt`, and runs the public no-docking demo.
It also installs a small `atlas` shim into a writable directory already
on `PATH` when possible. The shim runs
`atlas ...` commands inside `docking-env`

After install, run:
```bash
atlas setup-report
atlas doctor
```

`atlas setup-report` explains which tools were detected
`atlas doctor --pdb <PDB_ID> --ligands <library>` also checks a prepared
target/ligand combo before docking and reports fixes for common PDBQT issues.

Download shortcut:

```bash
bash tools/install_atlas_dev.sh
```


Or run the steps manually:
```bash
micromamba create -f environment.core.yml
micromamba activate docking-env
python -m pip install -e .
```

Optional extras:
```bash
micromamba env update -n docking-env -f environment.optional.yml
```


Start from the portable template:
```bash
cp config.example.txt config.txt
```

Check tool resolution before a real run:
```bash
atlas --verify-tools
```

Verification success criteria:
- exit code is `0`;
- output includes `Tools all successfully verified`;
- optional tool warnings may appear as `<N> optional missing` and do not fail verification.


Inspect the resolved runtime config:
```bash
atlas --print-effective-config
```

## Quick Start

Generate a tiny interactive report:

```bash
atlas demo --status-html
atlas status atlas_demo
```

Monitor a real run:
```bash
atlas status <RUN_ID>
atlas status <RUN_ID> --errors --explain
atlas status <RUN_ID> --watch 30
```

`atlas status` reports per-stage progress, chunk progress, run ETA, available Slurm state,
top root causes, key output paths, and next actions

Build a first real command before launching docking. If you do not already have
a PDB under `input_pdbs/`, let Atlas select one target and download it first:

```bash
atlas first-run --uniprot P04637 --ligands fda --small-library --fast --dry-run
```

Run every PDB in `input_pdbs/`:
```bash
atlas
```

Run two PDB IDs in fast mode:
```bash
atlas --pdb 1BN1 --pdb 2OJ9 --fast
atlas --pdb 1BN1,2OJ9 --fast
```

Use the small FDA test library:
```bash
atlas -test-fda --fast
```

Resume an existing run:
```bash
atlas --resume --run-id 20260425_120000
```


Select structures from gene symbols before installing input PDBs:

```bash
atlas targets panels
atlas targets species
atlas targets guide
atlas targets genes --genes EGFR ABL1 --gene-type KINASE --out analysis/gene_list/kinase_genes.csv
atlas targets install --panel kinases --max-return 3
atlas targets from-genes EGFR ABL1 --out analysis/gene_list/kinase_targets.csv
atlas targets search "BRCA DNA repair" --out analysis/gene_list/brca_topic_targets.csv
atlas targets install BRCA1 BRCA2
atlas targets install Trp53 --species mouse
atlas targets install --uniprot P04637
atlas targets install EGFR --ligand ATP --dedupe sequence-identity 90 --quality publication
atlas targets install EGFR --no-dedupe
atlas targets install --genes EGFR ABL1 --species "Homo sapiens" --taxonomy-id 9606 --methods X-RAY --resolution-max 2.8 --gene-type KINASE
```


Run a built-in target panel against an installable ligand library:
```bash
atlas run-panel kinases --ligands chembl --limit-ligands 1000 --fast
```

## Main Options

| Option | Purpose |
| --- | --- |
| `--pdb ID` | Run one PDB ID. Repeat for multiple proteins. |
| `--pdbs "ID1,ID2"` | Run a comma- or space-separated PDB list. |
| `--run-id RUN_ID` | Set an explicit run identifier. |
| `--resume` | Resume a previous run; use with `--run-id`. |
| `--fast` / `-fast` | Use low-exhaustiveness docking for tests. |
| `-test-fda` | Use the small FDA test library. |
| `--single NAME` | Dock one ligand matching `NAME`. |
| `--no-docking` | Prepare receptors without library docking. |
| `-dude` | Force DUD test-library mode |
| `--dud-library LIBRARY_SUBDIR` | Use one DUD/decoy library subdirectory for all proteins in the run. |
| `-bench`, `-bench2` | Run small reproducible benchmark profiles. |
| `-bench-micro`, `-bench-small` | Run reduced benchmark profiles for scheduler iteration. |
| `--verify-tools` | Validate required external tools and exit. |
| `--print-effective-config` | Print resolved config JSON and exit. |
| `-rebuild` | Rebuild ligand-library manifests and exit. |
| `--retain`, `--noretain`, `--retainmode MODE` | Override artifact retention for this run. |

Task-oriented aliases:

| Command | Purpose |
| --- | --- |
| `atlas init` | Create/update local `config.txt` and auto-detect tool paths. |
| `atlas demo` | Generate `outputs/data/atlas_demo/report.html` without docking. |
| `atlas smoke public` | Run the smoke fixture. |
| `atlas first-run --uniprot P04637 --ligands fda --dry-run` | Select a target, plan ligand prep, and print the first launch command. |
| `atlas status [RUN_ID] --explain` | Show operational progress, scheduler state, ETA, and next actions. |
| `atlas runs` | List recent run IDs with progress, failures, age, and report availability. |
| `atlas report RUN_ID` | Export `master_rows.csv` if needed and generate `outputs/data/<RUN_ID>/report.html`. - runs do this automatically|
| `atlas analysis report RUN_ID --status-html` | One command for master export, report HTML/YAML, and status HTML. |
| `atlas debug RUN_ID --deep` | Summarize structured run errors and next actions without reading raw logs. |
| `atlas artifacts measure --run-id RUN_ID` | Measure run-scoped output size under existing run output directories. |
| `atlas slurm submit --run-id RUN_ID --array 0-31%4 --cpus-per-task 8 --dry-run` | Plan a Slurm array submission with finalizer wiring. |
| `atlas slurm progress RUN_ID` | Show run progress with live Slurm probes when available. |
| `atlas slurm finalize RUN_ID --reconcile-only` | Reconcile distributed worker state into the run manifest. |
| `atlas screenshot RUN_ID --pdb 1ABC --top 20` | Render PNGs for top docked poses; default is full-protein cartoon context, with `--view-context pocket` and `--protein-style surface/none` available. |
| `atlas screenshot RUN_ID --pdb 1ABC --representatives --gallery-html` | Render best/median/worst scored examples and write a self-contained HTML gallery with embedded PNGs. |
| `atlas targets guide` | Prompt through gene/panel, organism, resolution, and query/install choices. |
| `atlas targets species` | List common organism aliases and inferred NCBI taxonomy IDs. |
| `atlas targets genes ...` | Write a reusable `gene,category` CSV for custom or built-in panels. |
| `atlas targets from-genes ...` | Query candidate PDB structures from gene symbols. |
| `atlas targets search "BRCA DNA repair"` | Full-text RCSB topic search for target discovery. |
| `atlas targets install BRCA1 BRCA2` | Select/download target PDBs from gene symbols and write target-install provenance. |
| `atlas targets install --panel kinases` | Select/download target PDBs from a built-in panel. |
| `atlas run-panel kinases --ligands chembl` | Install targets, install/prep ligands, and run the panel. |

Target gene searches default to human protein structures (`Homo sapiens`, taxonomy
`9606`, entity type `Protein`), all experimental methods, up to 10 PDBs per
gene, ligand-bound targets, and the resolution cutoff for entries where RCSB
reports resolution.
They can be narrowed with `--species`, `--taxonomy-id`,
`--entity-type`, `--methods`, `--resolution-max`, `--gene-type`, direct
`--uniprot`, required co-crystal ligands such as `--ligand ATP`, and
`--quality publication`; pass `--allow-apo` when apo structures are acceptable.
Search results are de-duplicated by 90% sequence identity by default; use
`--dedupe sequence-identity 70` or `--no-dedupe` to change that behavior.
Common organism aliases such as `human`, `mouse`, `rat`, `zebrafish`, `fly`,
`worm`, `yeast`, and `ecoli` infer the matching taxonomy ID; use
`atlas targets species` to list aliases.

For the complete built-in help:
```bash
atlas --help
```

### Library Selection

Atlas resolves the main library from `LIBRARY_SUBDIR_DEFAULT`
and DUD/decoy libraries from `TEST_LIBRARY_MAP` when `TEST_MODE_ENABLE`
includes `dud`. For a one-off run that should use the same decoy library for
every selected protein, use:

```bash
atlas --pdb 1BN1 --dud-library aa2ar --fast
```

This updates the in-memory `TEST_LIBRARY_MAP` for the selected proteins and
ensures `dud` is active for the run. Keep using `TEST_LIBRARY_MAP` in
`config.txt` when each target needs its own DUD-E library.

## Outputs

Outputs are run-scoped and organized in this manner :
Common output roots:

- `outputs/logs/`: top-level run logs
- `outputs/manifests/<RUN_ID>/`: run manifest
- `outputs/processed_pdbs/<RUN_ID>/<PDB>/`: cleaned and prepared receptor structures
- `prepped_ligands/`: prepared ligand libraries
- `outputs/docked/<RUN_ID>/<PDB>/`: docking scores and poses
- `outputs/post_docked/<RUN_ID>/<PDB>/`: post-docking analysis, reranking, and reports
- `outputs/configs/<RUN_ID>/`: per-run config snapshots
- `outputs/data/<RUN_ID>/`: reports, heatmap inputs, master rows, and run summaries

Audit file:
```text
outputs/manifests/<RUN_ID>/run_manifest.yaml
```

## Supported Scope
Maintained:
- `main.py`
- `src/`
- `analysis/`
- `chemdb/`
- `tools/` for developer/operations utilities
- config templates and docs

Tracked fixtures live under `chemdb/tests/fixtures/`; root input/output
