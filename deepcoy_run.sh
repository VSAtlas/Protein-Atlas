#!/usr/bin/env bash
set -euo pipefail

# Directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Atlas root:
# - Default: two levels up from this script (…/atlas2)
# - Override: export ATLAS_ROOT=/path/to/atlas2
ATLAS_ROOT="${ATLAS_ROOT:-"$(cd "$SCRIPT_DIR/../.." && pwd)"}"

# DeepCoy install location:
# - Default: $ATLAS_ROOT/tools/DeepCoy
# - Override: export DEEPCOY_DIR=/some/other/DeepCoy
DEEPCOY_DIR="${DEEPCOY_DIR:-"${ATLAS_ROOT}/tools/DeepCoy"}"

# Environment name:
# - Default: deepcoy-env
# - Override: export DEEPCOY_ENV=my_deepcoy_env
DEEPCOY_ENV="${DEEPCOY_ENV:-deepcoy-env}"

if [[ ! -d "$DEEPCOY_DIR" ]]; then
  echo "Error: DeepCoy directory not found at '$DEEPCOY_DIR'." >&2
  echo "Set ATLAS_ROOT or DEEPCOY_DIR to the correct path before running." >&2
  exit 1
fi

cd "$DEEPCOY_DIR"

# Forward all arguments to DeepCoy.py inside the DeepCoy env
exec micromamba run -n "$DEEPCOY_ENV" python DeepCoy.py "$@"
