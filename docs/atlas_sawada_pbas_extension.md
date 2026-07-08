# Atlas Sawada/PBAS Extension

Sawada et al. introduced PBAS-style proteome-wide docking profiles:
drug-level vectors of predicted binding affinities across human proteins,
used for therapeutic indication and side-effect prediction.

Atlas does not attempt to replicate Sawada's proteome-wide scale by default.
Instead, Atlas uses PBAS/Sawada-style profiles as a baseline and extends the
framework with pair-level statistical significance, Vina/SCORCH consensus,
MM/GBSA refinement, exposure plausibility, tissue expression, target-ADR
evidence, pathway evidence, and mechanism-graph interpretation.

Preferred interpretation:

Atlas provides complementary pair-level mechanistic prioritization on top of
PBAS-style docking-profile side-effect prediction.

Avoid claiming that Atlas proves Sawada wrong, validates all predicted side
effects, or predicts clinical ADRs without experimental or prospective support.

## Commands

Build PBAS pair scores:

```bash
python -m analysis.cli.build_pbas_baseline \
  --pbas-matrix data/external/pbas/pbas_matrix.parquet \
  --pbas-pocket-info data/external/pbas/pocket_info.csv \
  --mapping data/mappings/pbas_atlas_mapping.csv \
  --out outputs/analysis/pbas/pbas_pair_scores.csv
```

Build Atlas drug x target profiles:

```bash
python -m analysis.cli.build_binding_profiles \
  --pair-table outputs/analysis/pair_feature_table.csv \
  --value atlas_score \
  --out outputs/analysis/profiles/atlas_profile.parquet
```

Run the Sawada-style drug-level baseline:

```bash
python -m analysis.cli.run_sawada_baseline \
  --profile outputs/analysis/profiles/atlas_profile.parquet \
  --labels data/external/sider/drug_side_effect_labels.csv \
  --model logistic_regression \
  --split drug_holdout \
  --out-dir outputs/analysis/sawada_baseline
```

Compare PBAS and Atlas pair scores:

```bash
python -m analysis.cli.compare_pbas_atlas \
  --pbas-pairs outputs/analysis/pbas/pbas_pair_scores.csv \
  --atlas-pairs outputs/analysis/pair_feature_table.csv \
  --out-dir outputs/analysis/pbas_vs_atlas
```
