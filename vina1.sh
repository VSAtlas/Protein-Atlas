#!/usr/bin/env bash
set -euo pipefail

CPU="${VINA_CPU:-1}"        # override with: export VINA_CPU=8
ENV_NAME="${VINA_ENV:-docking-env}"

# If caller explicitly provided a path to Vina, honor it
if [ -n "${VINA_EXE:-}" ] && [ -x "${VINA_EXE}" ]; then
  exec "$VINA_EXE" --cpu "$CPU" "$@"
fi

# If Vina is already on PATH (env activated), use it
if command -v vina >/dev/null 2>&1; then
  exec vina --cpu "$CPU" "$@"
fi

# Prefer micromamba-run, which does not require activation
if command -v micromamba >/dev/null 2>&1; then
  exec micromamba run -n "$ENV_NAME" --no-capture-output vina --cpu "$CPU" "$@"
fi

# Accept mamba/conda as fallback
if command -v mamba >/dev/null 2>&1; then
  exec mamba run -n "$ENV_NAME" vina --cpu "$CPU" "$@"
fi
if command -v conda >/dev/null 2>&1; then
  exec conda run -n "$ENV_NAME" vina --cpu "$CPU" "$@"
fi

# Last resort: try MAMBA_ROOT_PREFIX layout
if [ -n "${MAMBA_ROOT_PREFIX:-}" ] && [ -x "$MAMBA_ROOT_PREFIX/envs/$ENV_NAME/bin/vina" ]; then
  exec "$MAMBA_ROOT_PREFIX/envs/$ENV_NAME/bin/vina" --cpu "$CPU" "$@"
fi

echo "vina1.sh: could not locate 'vina' or a conda/mamba runner. Install micromamba and create env '$ENV_NAME'." >&2
exit 127
