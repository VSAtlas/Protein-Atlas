# Atlas Molecular Baseline Adapters

These adapters keep Chemprop and PyTDC optional. Importing Atlas ML code does
not import either package; only the adapter command that needs an external tool
touches it.

## Chemprop Baseline Handoff

Stage an Atlas model-ready CSV into the minimal Chemprop shape:

```bash
python -m analysis.cli.run_chemprop_training \
  --dataset outputs/data/<RUN_ID>/ml/training_pass/_spd_expert_tables/model_ready/ml_spd_binding_model_ready.csv \
  --label spd_binding_label \
  --out-dir outputs/data/<RUN_ID>/ml/chemprop/chemprop_binding_baseline
```

Outputs:

- `chemprop_input.csv`: `smiles,<label>` input for Chemprop.
- `chemprop_metadata.csv`: Atlas row IDs and metadata for joining predictions
  back to `drug_id`, `target_id`, label source, split metadata, and source
  lineage.
- `chemprop_training_manifest.json`: row counts, inferred task type, column
  mapping, and policy notes.

To launch Chemprop after staging:

```bash
python -m analysis.cli.run_chemprop_training \
  --dataset outputs/data/<RUN_ID>/ml/training_pass/_spd_expert_tables/model_ready/ml_spd_binding_model_ready.csv \
  --label spd_binding_label \
  --out-dir outputs/data/<RUN_ID>/ml/chemprop/chemprop_binding_baseline \
  --run
```

The default run shape is Chemprop v2-style `chemprop train`. Use
`--command-style v1` for `chemprop_train`, `--chemprop-bin` for a custom
executable, and repeat `--chemprop-extra-arg` for version-specific options.

## TDC Dataset Staging

Stage a PyTDC dataset directly when PyTDC is installed:

```bash
python -m analysis.cli.stage_tdc_dataset \
  --tdc-module single_pred \
  --tdc-class Tox \
  --name hERG \
  --out outputs/data/<RUN_ID>/ml/tdc/herg.tsv
```

Stage an already downloaded CSV without importing PyTDC:

```bash
python -m analysis.cli.stage_tdc_dataset \
  --input-csv ~/datasets/tdc_herg.csv \
  --dataset-name hERG \
  --out outputs/data/<RUN_ID>/ml/tdc/herg.tsv
```

Outputs are Atlas-style benchmark source rows plus
`<output>.manifest.json`. Rows are marked `benchmark_only=True` by default; use
`--training-allowed` only after source lineage, licensing, overlap, leakage,
and split-policy review.

Install optional molecular tooling outside the repo with:

```bash
bash tools/installers/install_mlops_tools.sh --molecular
```
