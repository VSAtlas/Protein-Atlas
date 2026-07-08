# input_and_export_functions.py

## Migration note (2026-02-14)

- Canonical FDA decoy normalization output is now `z_vs_decoys`.
- Legacy compatibility alias `t_vs_decoys` is still written with identical values.
- New canonical entrypoint: `annotate_fda_long_csv_with_z_scores_vs_decoys(...)`.
- Legacy entrypoint `annotate_fda_long_csv_with_t_scores_vs_decoys(...)` is retained as a wrapper.
