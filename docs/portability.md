# Portability Guidance

## Goal
Run the same pipeline across machines with minimal machine-specific config.

## Templates
- `config.example.txt`: minimal portable template (recommended starting point)
- `config.full.example.txt`: full optional template with advanced knobs

Keep these templates at the repository root. Local `config.txt` files and
`config.txt.backup*` files are host-specific working state and should not be
tracked.

## Preferred config shape
Use a minimal `config.txt` with root/prefix keys:

- `OVERALL_DIR`
- `ADFRSUITE_BIN`
- `PHENIX_DIR`
- `PHENIX_LIB_PATH`
- `P2RANK_PATH`
- `MMGBSA_AMBERTOOLS_PREFIX` (or `AMBERTOOLS_PREFIX`)

All other executable paths are resolved through:

1. config-derived paths from prefix keys
2. environment overrides
3. `PATH` fallback

## Precedence
Runtime config precedence is:

1. defaults
2. `config.txt` values
3. environment overrides
4. CLI flags

## Legacy compatibility
Legacy MGLTools keys are still accepted by config normalization for old local
configs, but receptor PDBQT preparation now uses Meeko in the open stack.

## Effective config inspection
Use:

```bash
atlas --print-effective-config
```

This prints resolved runtime config as pretty JSON and exits.

## Verification
Use CLI verification only:

```bash
atlas --verify-tools
```

Expected output:

- Success: `Tools all successfully verified`
- Failure:
  - `<N> missing`
  - one line per missing tool key

Optional toolchains (`MMGBSA`, `GNINA`, `DOCK6`, `LeDock`) are configurable but do not fail the core `--verify-tools` check. `SCORCH` is also checked as optional/non-blocking: missing SCORCH script and/or unusable SCORCH env are reported but do not fail verification. SCORCH env probing resolves the runner via `MICROMAMBA_EXE`/`MAMBA_EXE`, PATH, `~/micromamba/bin/micromamba`, then conda fallback.
