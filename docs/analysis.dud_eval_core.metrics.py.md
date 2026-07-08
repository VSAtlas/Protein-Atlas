# analysis/dud_eval_core/metrics.py bug-fix note (2026-02-14)

## Summary
- Fixed two metric-definition bugs for DUD eval reporting:
  - BEDROC now uses **1-based active ranks** in the exponential term.
  - EF@x% now uses `k = ceil(f * N)` for top-fraction cutoff.

## What Changed
- File: `analysis/dud_eval_core/metrics.py`
- Function `bedroc(...)`:
  - Changed active rank indexing from 0-based to 1-based before computing
    `sum(exp(-alpha * rank / N))`.
- Function `ef_at_fractions(...)`:
  - Changed cutoff from `round(f * N)` to `ceil(f * N)`.
  - Clarified semantics: EF@x% is computed on the top `x%` of the **total ranked list**
    (not ROC/FPR-based EF).

## Why
- 1-based BEDROC ranking is consistent with standard BEDROC/RIE formulations and
  common toolkit implementations.
- `ceil(f * N)` avoids dropping boundary compounds at small sample sizes and aligns
  with common enrichment reporting practice.

## Impact
- BEDROC values may decrease compared with prior outputs that used 0-based ranks.
- EF values can change for cases where `round(f * N) != ceil(f * N)`.
- Metric column names are unchanged (`BEDROC_alpha_20`, `EF@1%`, `EF@2%`, `EF@5%`, `EF@10%`).
