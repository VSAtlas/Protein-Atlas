# tools/progress_snapshot.sh

## 2026-03-02

- Hardened status counting logic to avoid hanging when no combo chunk-plan files exist.
- Switched recursive post-docked file counting to an `rg --files` fast path (with `find` fallback when `rg` is unavailable).
- Preserves existing output contract:
  - `run=<id> chunks=<done>/<planned> pct=<...> claims=<...> phase=<...> scorch=<...> consensus=<...> combos=<...>/<...> failed=<...>`
