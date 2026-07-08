#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ATLAS_GRAPH_ML_ENV:-docking-env}"
MAMBA_BIN="${MAMBA_BIN:-micromamba}"

if ! command -v "$MAMBA_BIN" >/dev/null 2>&1; then
  echo "micromamba not found. Set MAMBA_BIN=/path/to/micromamba." >&2
  exit 1
fi

PYTHON="$($MAMBA_BIN run -n "$ENV_NAME" python -c 'import sys; print(sys.executable)')"

"$PYTHON" -m pip install torch==2.2.1+cpu --index-url https://download.pytorch.org/whl/cpu
"$PYTHON" -m pip install pykeen torch-geometric dgl pykan torchdata==0.7.1
"$PYTHON" -m analysis.cli.audit_graph_model_backends --out "outputs/data/graph_model_backends.csv"

echo "Graph/KGE/KAN backend install complete for environment: $ENV_NAME"
