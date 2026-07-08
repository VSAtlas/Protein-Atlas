# analysis/reporting/heatmap_html.py schema-hardening note (2026-02-14)

## Problem
- Heatmap loaders and aggregators still had transitional logic around legacy selected-score naming.
- During migration, this created inconsistent semantics and brittle test behavior.

## Change
- Enforced canonical selected-score key: `z_selected`.
- Removed legacy fallback behavior for `t_selected` in CSV/parquet loader paths.
- Kept aggregation and rendering paths strictly on `z_selected`.

## Why
- Makes score semantics unambiguous.
- Prevents silent mixed-schema behavior.
- Aligns with canonical reporting and parquet schema rules.
