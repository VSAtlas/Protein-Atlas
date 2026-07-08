# Notes: `src/post_docking/rescoring/rescoring_cnn.py`

- Switched consensus input resolution to reuse `find_consensus_csv(...)` from reranker when available, with canonical fallback.
- Fixed root resolution call site to use `resolve_roots(args, cfg)` and include `processed_root`.
- Hardened `resolve_roots(...)` to tolerate missing optional CLI fields (`docked_root`, `post_docked_root`) via `getattr`, avoiding `AttributeError` in default CLI usage.
