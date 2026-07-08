# Atlas New User Quickstart

This path is for a fresh clone on Linux or an HPC login node.

## One Command Install

Choose the path that matches your role:

| Persona | Command | Use when |
| --- | --- | --- |
| Public researcher | `bash tools/installers/install_atlas_publication_stack.sh --public-smoke` | You want the supported fresh-clone install and publication-safe smoke check. |
| Developer | `bash tools/install_atlas_dev.sh` | You will edit Atlas code and run the stricter developer gate. |
| Minimal reproducer | `micromamba create -f environment.core.yml && python -m pip install -e .` | You only need the core Python package and will provide external tools yourself. |

```bash
bash tools/installers/install_atlas_publication_stack.sh --public-smoke
```

The default install creates the core `docking-env` environment, installs Atlas
editable, creates `config.txt`, auto-detects legal/open tools available on the
machine, installs a small `atlas` shim into a writable directory already on
`PATH` when possible, and runs the public no-docking smoke check. The shim runs
each `atlas ...` command inside `docking-env`, so users usually do not need to
activate the environment first. If no writable PATH directory is available, the
installer falls back to `${HOME}/.local/bin`.

On clusters with separate HOME, work, and scratch filesystems, you can force a
site-appropriate shim directory or environment location with two optional
settings:

```bash
ATLAS_SHIM_DIR=/path/on/your/PATH bash tools/installers/install_atlas_publication_stack.sh
ATLAS_ENV_PREFIX=/scratch/$USER/micromamba/envs/docking-env atlas smoke public
```

Atlas does not bundle proprietary, restricted, or bring-your-own-license tools.
The installer reports those as optional/BYOL instead of treating them as a
failed public install.

## Check Readiness

```bash
atlas setup-report
atlas doctor
```

If `atlas` is not found after installation, add the fallback shim directory to
`PATH` or set `ATLAS_SHIM_DIR` to another writable PATH directory before
rerunning the installer:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Use `atlas setup-report` for a wetlab-facing summary: what works now, which
tool paths were detected, and which optional tools still need registration.
Use `atlas doctor` when you need the lower-level environment table with a
`fix` column.

Before a real docking run, you can point doctor at a target and prepared ligand
library:

```bash
atlas doctor --pdb 1ABC --ligands fda
```

This checks that the raw input PDB, prepared receptor PDBQT, ligand PDBQT files,
and optional Vina config paths line up before Vina reports harder-to-read parse
errors.

## Try Atlas Without Docking

```bash
atlas demo
atlas status atlas_demo --html
```

This creates a tiny local report under `outputs/data/atlas_demo/` without requiring
external docking tools or private data.

## Prepare First Inputs

```bash
atlas init --write-example-inputs
```

This writes `docs/examples/new_user_inputs.md` with input naming notes. If you
already have receptor structures, Atlas expects them in the configured input
root, usually `input_pdbs/<PDB_ID>.pdb`.

## Build A First Real Run

```bash
atlas first-run --uniprot P04637 --ligands fda --small-library --fast --dry-run
```

This selects one TP53 structure, plans the FDA small-library run, and prints the
command it would run. Use a different UniProt accession, `--genes`, `--panel`,
or `--pdb <ID>` for your own target.

When the dry run is ready:

```bash
atlas first-run --uniprot P04637 --ligands fda --small-library --fast --yes --no-dry-run
```

## Run Artifacts And Slurm

For reporting, use the bundled analysis wrapper:

```bash
atlas analysis report <RUN_ID> --status-html
atlas reproduce bundle <RUN_ID>
```

`atlas reproduce bundle` reuses `outputs/data/<RUN_ID>/` and writes an index of
the existing manifest, config snapshot, master rows, report, status dashboard,
and throughput files.

For array runs, use the Slurm wrapper instead of calling the site template
directly:

```bash
atlas slurm submit --run-id <RUN_ID> --array 0-31%4 --cpus-per-task 8 --dry-run
atlas slurm progress <RUN_ID>
atlas slurm finalize <RUN_ID> --reconcile-only
```

The wrapper passes run-specific arguments into the maintained `tools/slurm/`
scripts. Choose array concurrency and `--cpus-per-task` for the actual cluster
allocation and queue policy.

## Optional Tool Registration

Only configure these when your workflow needs them:

```bash
bash tools/installers/install_atlas_publication_stack.sh --adfrsuite --adfrsuite-root /path/to/adfrsuite
bash tools/installers/install_atlas_publication_stack.sh --p2rank --p2rank-root /path/to/p2rank
bash tools/installers/install_atlas_publication_stack.sh --scorch --scorch-source /path/to/scorch
```

After any registration step:

```bash
atlas setup-report
```

## Running Without Activating A Shell

The installer-created `atlas` shim is the supported "activate for this call
only" path. It dynamically discovers `micromamba` or `conda` at runtime and
dispatches into `docking-env`. This makes the shim more robust across HPC
systems where HOME, work, and scratch paths differ between machines.

If the cluster exposes the environment manager through modules, load that module
before using the shim. If the environment itself is under scratch/work, set
`ATLAS_ENV_PREFIX` to the environment prefix.

Keep using the installer for one-command setup. `pip install` still should not
create or mutate a micromamba/conda environment as a side effect; pip is best
kept to installing the Atlas package inside an already selected environment.
