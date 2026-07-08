# External Tools

Atlas is a Python workflow that shells out to several established docking,
preparation, and analysis tools. For a fresh clone, the recommended one-command
install is:

```bash
bash tools/installers/install_atlas_publication_stack.sh
```

That command installs the core legal/open environment, installs Atlas editable,
creates and auto-populates `config.txt`, and runs the no-proprietary public demo
smoke.

Publication-tier installer scripts are documented in
[`tools/installers/INSTALLERS.md`](../tools/installers/INSTALLERS.md).

Atlas does not redistribute proprietary/commercial binaries or license files.
For licensed tools, users must install them independently and provide local
paths in `config.txt` (BYOL model).

The helper creates or activates the `docking-env` environment with micromamba
or conda, installs `atlas` in editable mode, creates `config.txt` if needed,
auto-detects local tool paths, installs an outer `atlas` shim, and checks that
the CLI starts. The shim is written to a writable PATH directory when possible
and falls back to `${HOME}/.local/bin`; set `ATLAS_SHIM_DIR` to force a site
specific location. The shim dynamically discovers micromamba/conda and runs
`atlas ...` commands inside `docking-env` for that call only, so users do not
need to activate the environment first when the shim directory is on `PATH`.
For HPC layouts, use `ATLAS_SHIM_DIR` when the shim belongs in a site-specific
PATH directory and `ATLAS_ENV_PREFIX` when the environment lives under
scratch/work.
Detection uses `PATH`, the active conda/micromamba environment, and common
local environment prefixes; it does not recursively scan the filesystem.

Disable auto-detection when you want a fully manual config:

```bash
ATLAS_AUTO_CONFIG=0 bash tools/install_atlas_dev.sh
```

Run or rerun the same detection later:

```bash
atlas init
```

The older developer shortcut remains available:

```bash
bash tools/install_atlas_dev.sh
```

To also run external tool validation:

```bash
ATLAS_VERIFY_TOOLS=1 bash tools/install_atlas_dev.sh
```

If you are already inside a prepared Python 3.10 environment:

```bash
ATLAS_SKIP_ENV_CREATE=1 bash tools/install_atlas_dev.sh
```

Manual core install is:

```bash
micromamba create -f environment.core.yml
micromamba activate docking-env
python -m pip install -e .
```

Optional extras can be layered on after core install:

```bash
micromamba env update -n docking-env -f environment.optional.yml
```

Legacy compatibility note:

- `environment.yml` remains in-repo for historical/internal workflows.
- Publication-facing setup should treat `environment.core.yml` as the default.

Then copy `config.example.txt` to `config.txt`, set only the paths that apply to
your machine, and check resolution:

```bash
cp config.example.txt config.txt
atlas --verify-tools
```

`atlas --verify-tools` success criteria:

- process exits with code `0`;
- stdout includes `Tools all successfully verified`;
- stdout includes `[verify-tools]` runtime probe details for Meeko ligand prep;
- optional tool warnings (`<N> optional missing`) are non-blocking and do not
  fail the command.

Publication-safe public smoke fixture (no proprietary tool installs required):

```bash
bash tools/public_smoke_check.sh
```

The public smoke now imports the CLI and writes a small generated report under
`outputs/data/atlas_demo/report.html` on fresh clones.

## Required For Standard Docking

| Tool | Used For | Config Key |
| --- | --- | --- |
| AutoDock Vina | Primary docking engine | `VINA_EXE` or `VINA_PATH` |
| Meeko | SDF/restored-SDF ligand PDBQT and receptor PDBQT preparation | installed in environment |
| P2Rank | Pocket detection fallback/planning | `P2RANK_PATH` |

Vina-GPU probe code is retained only as an isolated experimental module. The
standard Atlas docking pipeline is CPU Vina only and does not expose Vina-GPU
configuration keys.

## Built-In Ligand Library Sources

Atlas can cache source SDFs and prepare PDBQT libraries with Meeko:

```bash
atlas ligands sources
atlas ligands install chembl --limit 1000
atlas ligands install chebi --limit 1000
atlas ligands install coconut --limit 1000
atlas ligands install hmdb --source-sdf ~/Downloads/structures.zip --limit 100
atlas chembl --pdb 1ABC --fast
```

Downloaded source files are stored under `LIGAND_SOURCE_CACHE_DIR` and prepared
PDBQT libraries are written under the configured `PREPPED_LIGANDS_DIR`. Built-in
sources include ChEMBL bioactive drug-like molecules, ChEBI high-confidence
3-star chemical structures, COCONUT natural products, DrugCentral/FDA drugs, and
HMDB metabolite structures. Source manifests include source URLs,
license/citation notes, raw SDF checksums, and prepared ligand counts.

HMDB bulk downloads may be protected by an interactive browser/form or
Cloudflare challenge. When the direct HMDB URL is blocked, download the HMDB
structures SDF or ZIP from the HMDB downloads page in a browser, then pass it to
Atlas with `--source-sdf`; Meeko preparation and manifest writing are unchanged.
ChEMBL and ChEBI are discovered from EMBL-EBI HTTPS indexes. COCONUT is
discovered from the current COCONUT download page and then fetched from its S3
archive URL.

## Optional Engines And Analysis

| Tool | Used For | Config Key |
| --- | --- | --- |
| Open Babel | Optional format conversion/debug utility; not used by the active ligand-prep writer | `OPENBABEL_PATH` |
| GNINA | CNN/minimized-affinity docking and rescoring | `GNINA_EXE`, `USE_GNINA` |
| LeDock | Additional docking consensus signal | `USE_LEDOCK` |
| DOCK6 | Additional docking consensus signal | `DOCK6_EXE`, `USE_DOCK6` |
| PHENIX | Structure cleanup workflows | `PHENIX_DIR`, `PHENIX_LIB_PATH` |
| AmberTools | MMGBSA and trajectory/topology utilities | `AMBERTOOLS_PREFIX`, `MMGBSA_AMBERTOOLS_PREFIX` |
| SCORCH | Post-docking rescoring | `SCORCH`, `SCORCH_DEVICE`, `SCORCH_GPU_IDS`, `SCORCH_ENV`, `SCORCH_ENV_PREFIX` |

SCORCH device selection is conservative. `SCORCH_DEVICE=cpu` hides GPUs from
SCORCH child processes. `SCORCH_DEVICE=auto` uses GPUs only when allocated or
visible GPUs are detected and the selected SCORCH TensorFlow environment reports
at least one usable GPU. In `auto`, Atlas keeps small SCORCH chunks on CPU to
avoid TensorFlow/ROCm startup overhead, while GPU runs use larger backend
default scoring batches. `SCORCH_DEVICE=amd` or `SCORCH_DEVICE=nvidia` forces
backend selection while still probing before use. `SCORCH_GPU_IDS` can restrict
Atlas to a subset of the visible GPUs. Backend env prefixes are auto-detected
from common names such as `scorch-rocm`, `scorch-amd`, `scorch-cuda`, and
`scorch-nvidia` beside the configured SCORCH env; advanced overrides remain
available as `SCORCH_ENV_PREFIX_AMD` and `SCORCH_ENV_PREFIX_NVIDIA`.

MM/GBSA requires AmberTools-compatible `antechamber`, `parmchk2`, `tleap`,
`cpptraj`, and `MMPBSA.py`/`MMPBSA.py.MPI`. The local install path is
`atlas/tools/envs/ambertools`; Atlas auto-discovers this prefix from
`OVERALL_DIR` when `MMGBSA_AMBERTOOLS_PREFIX` and `AMBERTOOLS_PREFIX` are blank.
Use `bash tools/installers/install_ambertools.sh --configure` to create this
environment reproducibly under the shared Atlas tools root.

BANANA inference is an optional binding-prior workflow. Use
`bash tools/installers/install_banana.sh` to install the pinned open-source
checkout under `ATLAS_TOOLS_DIR/banana` with its isolated Python environment at
`ATLAS_TOOLS_DIR/envs/banana`. The Atlas BANANA wrapper auto-discovers those
paths unless `BANANA_ROOT` or `BANANA_PYTHON` is set.

## Optional ML/MLOps And Molecular Design Tools

Use the ML/MLOps installer to keep research tooling out of the core docking
environment while making it reproducible on the shared Atlas layout:

```bash
bash tools/installers/install_mlops_tools.sh --mlops
bash tools/installers/install_mlops_tools.sh --all
```

The default root is `$ATLAS_TOOLS_DIR` or `../../tools` relative to this repo
(`/stor/work/VDS_Beckham/atlas/tools` on the shared filesystem). Installed
environments are written under `$ATLAS_TOOLS_DIR/envs/atlas-ml-*`. Pip,
conda/micromamba package caches, Matplotlib/XDG caches, and the installer-local
home directory are also routed under `$ATLAS_TOOLS_DIR` by default so live
cluster installs do not write quota-heavy metadata into `$HOME`.

| Tier | Installs | Intended use |
| --- | --- | --- |
| `--mlops` | MLflow, Optuna, DVC/DVCLive, W&B, ClearML | experiment tracking, model registry handoff, sweeps, dataset manifests |
| `--molecular` | RDKit, Chemprop, PyTDC, DeepChem | molecular baselines and external benchmark loaders |
| `--rl` | Gymnasium, Stable-Baselines3, Ray Tune/RLlib | general RL interfaces when molecule-specific loops are insufficient |
| `--generative` | REINVENT4 and MOSES editable installs, GuacaMol, source checkouts, and an RDKit env | generative/RL molecule design and benchmark reuse |

The installer writes `$ATLAS_TOOLS_DIR/ml_tools_manifest.txt` plus per-env
`pip freeze`, conda explicit specs, and generative source SHAs under
`$ATLAS_TOOLS_DIR/ml_tool_env_state/`. Generative source checkouts use
`GIT_LFS_SKIP_SMUDGE=1` by default because MOSES includes large LFS benchmark
data that may be unavailable on GitHub; override with
`ATLAS_ML_GIT_LFS_SKIP_SMUDGE=0` only when you explicitly need those data blobs.
The generative tier installs REINVENT4 and MOSES from their source checkouts into
`$ATLAS_TOOLS_DIR/envs/atlas-ml-generative`, then applies
`numpy<2 pandas<2` by default so RDKit and the legacy MOSES metrics import
cleanly in the shared CPU environment. Override
`ATLAS_GENERATIVE_COMPAT_PACKAGES` only after testing REINVENT and MOSES imports;
set `ATLAS_REINVENT_TORCH_INDEX` if you need a CUDA-specific PyTorch wheel index.

The generative adapter CLIs are capture-first wrappers. They do not import
REINVENT4, GuacaMol, or MOSES inside Atlas. Capture an existing REINVENT4
sample file with:

```bash
python -m analysis.cli.run_reinvent_generation \
  --generated-smiles /path/to/reinvent_samples.csv \
  --out-dir outputs/data/<RUN_ID>/ml/generative/reinvent4
```

External execution is opt-in. Add `--execute --command ...` only when launching
an installed tool; the command is passed without a shell and stdout/stderr are
written beside the manifest:

```bash
python -m analysis.cli.run_reinvent_generation \
  --generated-smiles /path/to/reinvent_samples.csv \
  --out-dir outputs/data/<RUN_ID>/ml/generative/reinvent4 \
  --execute --command -- "$ATLAS_TOOLS_DIR/envs/atlas-ml-generative/bin/reinvent" /path/to/config.toml -d cpu
```

Prepare GuacaMol or MOSES benchmark inputs and capture benchmark result files
without importing either package:

```bash
python -m analysis.cli.benchmark_generated_smiles \
  --suite guacamol \
  --generated-smiles outputs/data/<RUN_ID>/ml/generative/reinvent4/generated_smiles.csv \
  --benchmark-results /path/to/guacamol_results.json \
  --out-dir outputs/data/<RUN_ID>/ml/generative/guacamol
```

Atlas training logs local `experiment_runs.jsonl`, `registry.json`, and
`model_registry_index.csv` by default. External trackers are opt-in:

```bash
ATLAS_ENABLE_MLFLOW=1 ATLAS_MLFLOW_TRACKING_URI=file:/path/to/mlruns atlas ml train --run-id <RUN_ID>
ATLAS_WANDB_PROJECT=atlas-ml WANDB_MODE=offline atlas ml train --run-id <RUN_ID>
ATLAS_CLEARML_PROJECT=AtlasML atlas ml train --run-id <RUN_ID>
```

## Publication-Friendly Install Notes

Keep one tested recipe per supported platform or HPC image that includes:

- exact upstream download links or package-manager commands for each external
  tool;
- expected executable names on `PATH`;
- the corresponding `config.txt` keys;
- a copy-paste `atlas --verify-tools` example that satisfies the success
  criteria above;
- execution of `bash tools/public_smoke_check.sh`.
