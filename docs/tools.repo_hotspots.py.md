# tools/repo_hotspots.py

## 2026-03-02

- Added a lightweight hotspot reporter for Python file size/line-count ranking.
- Default output focuses on large files (`--min-bytes 50000`) to guide safe simplification work.
- Excludes runtime artifact and cache roots by default (for example `docked/`, `post_docked/`, `.git/`, `node_modules/`).
- Example:
  - `python tools/repo_hotspots.py --top 12 --min-bytes 60000`
