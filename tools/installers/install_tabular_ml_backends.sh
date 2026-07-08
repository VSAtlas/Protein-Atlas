#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ATLAS_TABULAR_ML_ENV:-docking-env}"
MAMBA_BIN="${MAMBA_BIN:-micromamba}"

if ! command -v "$MAMBA_BIN" >/dev/null 2>&1; then
  echo "micromamba not found. Set MAMBA_BIN=/path/to/micromamba." >&2
  exit 1
fi

PYTHON="$($MAMBA_BIN run -n "$ENV_NAME" python -c 'import sys; print(sys.executable)')"

"$PYTHON" -m pip install numpy==1.26.4 lightgbm catboost interpret-core scikit-learn ydata-profiling phik dython networkx pyvis
"$PYTHON" -m pip install numpy==1.26.4
"$PYTHON" - <<'PY'
import importlib
import numpy

print(f"numpy {numpy.__version__}")
for module in [
    "lightgbm",
    "catboost",
    "interpret.glassbox",
    "sklearn",
    "ydata_profiling",
    "phik",
    "dython.nominal",
    "networkx",
    "pyvis.network",
]:
    importlib.import_module(module)
    print(f"{module} ok")
PY

echo "Tabular ML backend install complete for environment: $ENV_NAME"
