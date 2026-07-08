## Reproducibility Guide

This guide defines the minimum reproducible Atlas2 workflow for publication review.

## 1) Record the code state

Before a publication run, record:

- Git commit hash: `git rev-parse HEAD`
- Git tag or release label, if present
- Atlas package version from `pyproject.toml`
- Whether the worktree had local changes

Publication/preprint builds should use the tagged manuscript state, for example `atlas2-preprint-v0.1` once that tag exists.

## 2) Environment creation

```bash
bash tools/install_atlas_dev.sh
```

Equivalent manual setup:

```bash
micromamba create -f environment.core.yml
micromamba activate docking-env
python -m pip install -e .
```

Optional extras:

```bash
micromamba env update -n docking-env -f environment.optional.yml
```

## 3) Config setup

Start from the portable template:

```bash
cp config.example.txt config.txt
```

Set only machine-local tool roots and run roots in `config.txt`. Do not commit host-specific defaults back to `config.example.txt`.

Minimum keys to review:

- `OVERALL_DIR`
- `ADFRSUITE_BIN`
- `PHENIX_DIR`
- `PHENIX_LIB_PATH`
- `P2RANK_PATH`
- `MMGBSA_AMBERTOOLS_PREFIX`
- `AMBERTOOLS_PREFIX`
- `SCORCH_DEVICE`, optional `SCORCH_GPU_IDS`, and `SCORCH_ENV_PREFIX` /
  `SCORCH_ENV` if SCORCH rescoring is used

Inspect the resolved runtime config:

```bash
atlas --print-effective-config
```

## 4) Required external tools

Atlas can run with different docking backends, but publication runs must report exact tool stack:

- AutoDock Vina or the selected docking backend
- Meeko (ligand and receptor PDBQT preparation)
- Phenix / Reduce (for specific protein-prep paths)
- P2Rank (pocket detection path)
- AmberTools MM/GBSA components (if MM/GBSA rescoring used)
- RDKit-dependent chemistry tooling used by ligand prep
- SCORCH/GNINA/CNN rescoring environments if those rescoring paths are enabled

Verify local resolution:

```bash
atlas --verify-tools
```

Save the command output with the run notes because external tool versions and paths affect reproducibility.

## 5) Baseline acceptance checks

Before structural cleanup or manuscript claims, run:

```bash
atlas dev verify --full --fix --smoke
pytest chemdb/tests/test_main_full_run.py -q
pytest chemdb/tests/test_ligand_run_modes.py -q
pytest chemdb/tests/test_path_router.py -q
pytest chemdb/tests/test_record_data_csv.py -q
```

Known failures should be recorded in the handoff or release checklist with the exact failing command and likely cause. Do not silently skip failures.

## 6) Minimal smoke run

```bash
bash tools/public_smoke_check.sh
```

Expected outputs:

- smoke status logs under `outputs/logs/`
- run manifest under `outputs/manifests/<RUN_ID>/run_manifest.yaml`
- prepared receptor files under `outputs/processed_pdbs/<RUN_ID>/<PDB>/`
- any generated public-smoke artifacts documented by the smoke script output

Internal smoke mode:

```bash
atlas --test -fast
```

## 7) Small FDA test-library run

```bash
atlas -test-fda --fast
```

Expected output families:

- `outputs/docked/<RUN_ID>/...` score/pose outputs
- `outputs/post_docked/<RUN_ID>/...` post-docking artifacts (when enabled)
- `outputs/processed_pdbs/<RUN_ID>/...` receptor preparation artifacts
- `outputs/manifests/<RUN_ID>/run_manifest.yaml`

Capture the run ID, config snapshot, and a concise summary of generated score/report files.

## 8) Regenerate report / heatmap

Use the maintained Atlas wrapper for run reports. It exports master rows before
generating the report so a completed run has one normal reporting command:

```bash
atlas analysis report <RUN_ID> --status-html
atlas reproduce bundle <RUN_ID>
```

`atlas reproduce bundle` writes `reproducibility_bundle.json` and
`status_dashboard.json` under `outputs/data/<RUN_ID>/` by default. It indexes
the existing manifest, config snapshot, master rows, report, status dashboard,
and throughput files instead of creating a separate top-level `outputs/repro`
tree.

For the side-effect/ADR analysis workflow:

```bash
python -m analysis.pipeline \
  --screening-table <screening.csv> \
  --metadata-table <metadata.csv> \
  --outdir outputs/analysis
```

See `docs/atlas_analysis_pipeline.md` and `docs/schemas/` for the promoted table contracts.

## 9) Expected output files

Minimum reproducibility package for a run:

- `outputs/manifests/<RUN_ID>/run_manifest.yaml`
- `outputs/configs/<RUN_ID>/` config snapshot
- prepared receptor structures under `outputs/processed_pdbs/<RUN_ID>/`
- docking score summaries under `outputs/docked/<RUN_ID>/`
- post-docking/reranking/report artifacts under `outputs/post_docked/<RUN_ID>/`, when enabled
- exported reporting tables or heatmaps generated for the manuscript
- tool verification output
- exact commands used

## 10) Hardware assumptions

- Linux or Linux-like HPC environment
- Sufficient CPU for docking workloads
- Do not exceed 32 cores for docking runs in this repository policy
- GPU availability only if selected optional ML/GNINA/SCORCH workflows require it

Record CPU count, GPU model if used, scheduler mode, and whether the run was local or distributed.

## 11) Known nondeterminism

Potential sources:

- stochastic docking/search heuristics across engines
- filesystem ordering and distributed task timing
- optional external tool version differences
- pH/protomer/tautomer enumeration differences across tool versions
- floating-point and thread scheduling differences across CPU/GPU hardware
- missing-loop, water, ion, protonation, and binding-site choices

Reproducibility reports should record:

- git revision
- config snapshot
- run ID
- tool versions from `--verify-tools` output

## 12) External tool citation

Manuscript/preprint should cite each external docking/prep/scoring dependency used in the run stack, including the docking backend, ligand preparation tools, receptor preparation tools, pocket detection tools, MM/GBSA tools, and ML rescoring tools when enabled.

## 13) Manuscript/preprint reproduction checklist

For manuscript/preprint, include:

- exact commit hash
- tagged release identifier
- config snapshot used
- command lines used for smoke + main benchmark/demo run
- external tool versions and citations
- hardware/scheduler description
- input protein IDs and ligand-library provenance
- expected output files and schema versions
- known limitations relevant to the claims being made
