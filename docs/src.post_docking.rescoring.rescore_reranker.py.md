# src/post_docking/rescoring/rescore_reranker.py

## Migration note (2026-02-14)

- Canonical decoy-normalized columns are now:
  - `z_vs_decoys_consensus`
  - `z_vs_decoys_consensus_pre`
  - `z_vs_decoys_blend`
- Legacy aliases are still emitted and synchronized:
  - `t_vs_decoys_consensus`
  - `t_vs_decoys_consensus_pre`
  - `t_vs_decoys_blend`

## Edge-case fix

- Input rows that only had legacy `t_*` columns are now normalized on ingest and copied into canonical `z_*` columns.
- If both `z_*` and `t_*` are present but differ, canonical `z_*` wins and aliases are overwritten to match.
- `final_score` now prefers canonical `z_*` fields with legacy fallback.

## Consensus discovery hardening (2026-02-14)

- `find_consensus_csv(...)` now falls back to token-prefixed consensus filenames, including `*_consensus_docking_scores.csv`.
- Selection skips known non-primary files (`dud_consensus`, `rerank`, `pretty`) and ranks deterministic best candidates.
- This prevents false `missing_consensus` skips when docking emits library-token-prefixed consensus files.
