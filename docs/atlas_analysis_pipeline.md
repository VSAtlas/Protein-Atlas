# Atlas Analysis Pipeline

The `analysis` package turns existing Atlas screening outputs into manuscript-oriented pair tables, heatmaps, empirical significance tables, and enrichment summaries.

## Single Command

```bash
python -m analysis.pipeline \
  --screening-results outputs/docking/ranked_pairs.csv \
  --annotations data/annotations \
  --out-dir outputs \
  --permutations 10000
```

Outputs:

- `outputs/analysis/pair_feature_table.csv`
- `outputs/analysis/pair_feature_table.json`
- `outputs/analysis/heatmap_matrix.csv`
- `outputs/analysis/pair_significance.csv`
- `outputs/analysis/enrichment_summary.csv`
- `outputs/analysis/permutation_results.csv`
- `outputs/figures/atlas_heatmap.png`
- `outputs/figures/atlas_heatmap.svg`
- `outputs/figures/significance_heatmap.png`
- `outputs/figures/enrichment_curve.png`

## Stepwise Commands

```bash
python -m analysis.build_pair_table \
  --screening-results outputs/docking/ranked_pairs.csv \
  --annotations data/annotations \
  --out outputs/analysis/pair_feature_table.csv

python -m analysis.run_pair_significance \
  --pair-table outputs/analysis/pair_feature_table.csv \
  --score atlas_score \
  --out outputs/analysis/pair_significance.csv

python -m analysis.generate_heatmap \
  --pair-table outputs/analysis/pair_feature_table.csv \
  --value atlas_score \
  --row-order protein_class \
  --col-order ligand_chemotype \
  --significance-table outputs/analysis/pair_significance.csv \
  --matrix-out outputs/analysis/heatmap_matrix.csv \
  --out outputs/figures/atlas_heatmap.png

python -m analysis.run_enrichment \
  --pair-table outputs/analysis/pair_feature_table.csv \
  --score atlas_score \
  --label literature_supported_label \
  --permutations 10000 \
  --out outputs/analysis/enrichment_summary.csv \
  --permutation-out outputs/analysis/permutation_results.csv \
  --curve-out outputs/figures/enrichment_curve.png
```

`atlas_score` is a higher-is-better alias for the existing Atlas z-score export (`z_selected` when present). Raw lower-is-better docking scores are negated only when no z-score field is available. Heatmap SVG cells include hover titles with the z-score value; missing MM/GBSA values remain masked rather than imputed.

`literature_supported_label` is used only as an evaluation label for enrichment analysis. It is not included in `priority_score`.
