# main.py change note

## 2026-02-13

- Updated `--verify-tools` to include `SCORCH` as an optional check target.
- Verification still fails only on required tools; missing SCORCH now reports as `optional missing` and does not return a failure code.
- Added optional SCORCH env probe (`micromamba run -n/-p ... python -V`) during `--verify-tools`.
- Updated help text to clarify required vs optional verification behavior.
