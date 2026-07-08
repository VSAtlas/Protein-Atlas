# Atlas Publication Installers

These scripts provide publication-oriented setup tiers while keeping restricted
or licensed tools optional and user supplied.

## One-Command Install

```bash
bash tools/installers/install_atlas_publication_stack.sh
```

This is the recommended fresh-clone command. It creates or updates
`docking-env` from `environment.core.yml`, installs Atlas in editable mode,
creates local `config.txt`, auto-detects legal/open tools from the active
environment and `PATH`, installs the pinned Clustergrammer browser assets under
the external Atlas tools root, and runs the no-proprietary public demo smoke
that writes `data/atlas_demo/report.html`.

The core environment installs the legal/open conda-forge pieces needed for
standard Vina-based runs, including Python dependencies, RDKit, OpenMM/PDBFixer,
Open Babel, AutoDock Vina, Meeko, propka, pdb2pqr, Reduce/Probe, and the
open Dimorphite-DL HET protonation-state enumeration helper, PyMOL, and the
development quality tools.

The report asset installer uses `tools/report_assets/npm/package-lock.json` and
creates a repo-local `report_assets` pointer to
`$ATLAS_TOOLS_DIR/report_assets` (default `../../tools/report_assets`). Run it
directly when you only need to repair report assets:

```bash
bash tools/installers/install_clustergrammer_assets.sh
```

## Core Atlas Environment Only

```bash
bash tools/installers/install_atlas_core.sh
```

Use this when you want only the environment/editable install/config detection
and do not need the public smoke run.

## Optional External Tools

Atlas does not redistribute restricted tools or license files. Install or
register optional tools with explicit local sources:

```bash
bash tools/installers/install_adfrsuite.sh --existing-root /path/to/adfrsuite --configure
bash tools/installers/install_p2rank.sh --archive /path/to/p2rank-release.zip --configure
bash tools/installers/install_scorch.sh --source-dir /path/to/SCORCH --configure
bash tools/installers/install_ambertools.sh --configure
bash tools/installers/install_banana.sh
bash tools/installers/install_mlops_tools.sh --mlops
bash tools/installers/install_graph_ml_backends.sh
bash tools/installers/install_tabular_ml_backends.sh
```

`--configure` updates only local `config.txt`; it never modifies
`config.example.txt`.

By default, open-source external tools are installed outside this source tree
under `../../tools` relative to the repository, which is
`/stor/work/VDS_Beckham/atlas/tools` on the shared Atlas layout. Override this
with `ATLAS_TOOLS_DIR=/path/to/tools` when needed.

Optional tiers are opt-in:

```bash
bash tools/installers/install_atlas_publication_stack.sh \
  --core \
  --report-assets \
  --adfrsuite --adfrsuite-root /path/to/adfrsuite \
  --p2rank --p2rank-archive /path/to/p2rank-release.zip \
  --scorch --scorch-source /path/to/SCORCH \
  --ambertools \
  --banana \
  --configure
```

BANANA defaults to the upstream
`https://github.com/molecularmodelinglab/banana.git` checkout pinned at the
revision used by the current Atlas integration. Override with
`--banana-git-url` and `--banana-git-ref` to reproduce another vetted revision.

AmberTools defaults to a conda-forge `ambertools=24.8` environment under
`$ATLAS_TOOLS_DIR/envs/ambertools`. Override with `--ambertools-package` or
`--ambertools-prefix` if your platform needs a different build.

ML/MLOps tooling is intentionally optional and isolated from `docking-env`. Use
`bash tools/installers/install_mlops_tools.sh --mlops` for MLflow, Optuna, DVC,
DVCLive, W&B, and ClearML; add `--molecular`, `--rl`, `--generative`, or
`--all` to install Chemprop/TDC/DeepChem, Gymnasium/SB3/Ray, and molecule
generation tooling under `$ATLAS_TOOLS_DIR`. The generative tier installs
GuacaMol plus editable REINVENT4 and MOSES checkouts, defaults REINVENT to the
CPU PyTorch wheel index, and reapplies `numpy<2 pandas<2` compatibility pins for
RDKit/MOSES unless `ATLAS_GENERATIVE_COMPAT_PACKAGES` is overridden. The
installer also writes
per-env freeze/spec files and generative git SHAs under
`$ATLAS_TOOLS_DIR/ml_tool_env_state/`; its pip, conda, XDG, Matplotlib, and
installer-local home caches default to the tools root.

Graph/KGE backends are optional candidates for sparse mechanism experiments:

```bash
bash tools/installers/install_graph_ml_backends.sh
bash tools/installers/install_tabular_ml_backends.sh
```

The graph installer adds CPU PyTorch, PyKEEN, PyTorch Geometric, DGL, and
pykan into `docking-env`, then writes `outputs/data/graph_model_backends.csv`.
It pins `torch==2.2.1+cpu` and `torchdata==0.7.1` because the current DGL
wheel bundles GraphBolt libraries for the PyTorch 2.0-2.2 line, not newer 2.12
wheels.

The tabular installer also installs the optional read-only EDA stack used by
`atlas analysis dataset-eda`: ydata-profiling, Phi_K/phik, dython,
scikit-learn mutual information, NetworkX, and PyVis. Core column summaries and
numeric network fallback run without all optional packages; missing extras are
recorded in `dataset_eda_manifest.json`.

After configuring external tools, verify the local installation:

```bash
atlas --verify-tools
```
