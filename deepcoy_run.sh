#!/usr/bin/env bash
set -euo pipefail

# Optional: first arg can be "-j N" or "--jobs N" to run N instances in parallel
JOBS=1
if [[ $# -ge 2 && ( "${1:-}" == "-j" || "${1:-}" == "--jobs" ) ]]; then
    JOBS="$2"
    shift 2
fi

# Directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Atlas root:
# - Default: two levels up from this script (…/atlas2 or …/atlas)
# - Override: export ATLAS_ROOT=/path/to/atlas_root
ATLAS_ROOT="${ATLAS_ROOT:-"$(cd "$SCRIPT_DIR/../.." && pwd)"}"

# DeepCoy install location:
# - Default: $ATLAS_ROOT/tools/DeepCoy
# - Override: export DEEPCOY_DIR=/some/other/DeepCoy
DEEPCOY_DIR="${DEEPCOY_DIR:-"${ATLAS_ROOT}/tools/DeepCoy"}"

# Environment name:
# - Default: deepcoy-env
# - Override: export DEEPCOY_ENV=my_deepcoy_env
DEEPCOY_ENV="${DEEPCOY_ENV:-deepcoy-env}"

# Micromamba command:
# - Default: "micromamba" (whatever your PATH provides)
# - Override: export MICROMAMBA_BIN=/full/path/to/micromamba
MM_BIN="${MICROMAMBA_BIN:-micromamba}"

if [[ ! -d "$DEEPCOY_DIR" ]]; then
  echo "Error: DeepCoy directory not found at '$DEEPCOY_DIR'." >&2
  echo "Set ATLAS_ROOT or DEEPCOY_DIR to the correct path before running." >&2
  exit 1
fi

cd "$DEEPCOY_DIR"

if [[ "$JOBS" -le 1 ]]; then
  # Single DeepCoy run (original behavior)
  exec "$MM_BIN" run -n "$DEEPCOY_ENV" python DeepCoy.py "$@"
else
  echo "Launching ${JOBS} DeepCoy jobs in parallel..."
  for i in $(seq 1 "$JOBS"); do
    echo "  -> DeepCoy job $i"
    "$MM_BIN" run -n "$DEEPCOY_ENV" python DeepCoy.py "$@" &
  done
  wait
  echo "All DeepCoy jobs finished."
fi
