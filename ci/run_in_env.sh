#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ATLAS_ENV_NAME:-atlas-ci}"

if command -v micromamba >/dev/null 2>&1; then
  echo "[ci] using micromamba env: $ENV_NAME"
  micromamba create -y -n "$ENV_NAME" -f environment.yml || true
  exec micromamba run -n "$ENV_NAME" pytest -q tests/test_ions_acceptance.py -k ions_acceptance --maxfail=1 -vv
elif command -v conda >/dev/null 2>&1; then
  echo "[ci] using conda env: $ENV_NAME"
  conda env update -n "$ENV_NAME" -f environment.yml || conda env create -n "$ENV_NAME" -f environment.yml
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
  exec pytest -q tests/test_ions_acceptance.py -k ions_acceptance --maxfail=1 -vv
else
  echo "Neither micromamba nor conda found. Install micromamba: https://mamba.readthedocs.io/"
  exit 1
fi
