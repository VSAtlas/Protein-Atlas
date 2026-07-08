# tool_resolver.py change note

## 2026-02-13

- Updated `resolve_micromamba()` to also honor `MAMBA_EXE` (in addition to `MICROMAMBA_EXE`), matching micromamba shell-hook environments where activation works without `micromamba` on `PATH`.
