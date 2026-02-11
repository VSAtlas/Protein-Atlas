# Portability Guidance

## Goal
Run the same pipeline across machines with minimal machine-specific config.

## Templates
- `config.example.txt`: minimal portable template (recommended starting point)
- `config.full.example.txt`: full optional template with advanced knobs

## Preferred config shape
Use a minimal `config.txt` with root/prefix keys:

- `OVERALL_DIR`
- `MGLTOOLS_PATH`
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
Legacy per-tool keys are still accepted for backward compatibility. Recommended migration:

- `MGLTOOLS_PYTHON` -> derived from `MGLTOOLS_PATH`
- `PREPARE_RECEPTOR_SCRIPT` -> derived from `MGLTOOLS_PATH`
- `PREPARE_LIGAND_SCRIPT` -> derived from `MGLTOOLS_PATH` (or PATH fallback)

## Effective config inspection
Use:

```bash
python main.py --print-effective-config
```

This prints resolved runtime config as pretty JSON and exits.

## Verification
Use CLI verification only:

```bash
python main.py --verify-tools
```

Expected output:

- Success: `Tools all successfully verified`
- Failure:
  - `<N> missing`
  - one line per missing tool key

Optional toolchains (`MMGBSA`, `GNINA`, `DOCK6`, `LeDock`) are configurable but do not fail the core `--verify-tools` check.
