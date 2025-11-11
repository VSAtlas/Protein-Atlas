#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ATLAS_ENV_NAME:-atlas-ci}"
if [ "$#" -gt 0 ]; then
  CMD="$*"
else
  CMD="pytest -q tests/test_ions_acceptance.py -k ions_acceptance --maxfail=1 -vv"
fi

if command -v micromamba >/dev/null 2>&1; then
  echo "[ci] using micromamba env: $ENV_NAME"
  micromamba create -y -n "$ENV_NAME" -f environment.yml || true
  exec micromamba run -n "$ENV_NAME" bash -lc "$CMD"
elif command -v conda >/dev/null 2>&1; then
  echo "[ci] using conda env: $ENV_NAME"
  conda env update -n "$ENV_NAME" -f environment.yml || conda env create -n "$ENV_NAME" -f environment.yml
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
  exec bash -lc "$CMD"
else
  echo "[ci] neither micromamba nor conda found; running in current environment"
  exec bash -lc "$CMD"
fi
