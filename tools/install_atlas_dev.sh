#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
env_name="${ATLAS_ENV_NAME:-docking-env}"
skip_env="${ATLAS_SKIP_ENV_CREATE:-0}"
verify_tools="${ATLAS_VERIFY_TOOLS:-0}"
env_file="${ATLAS_ENV_FILE:-environment.core.yml}"
install_optional="${ATLAS_INSTALL_OPTIONAL:-0}"
auto_config="${ATLAS_AUTO_CONFIG:-1}"
configure_args="${ATLAS_CONFIGURE_ARGS:-}"

cd "$repo_root"

if [[ "$skip_env" != "1" ]]; then
  if [[ ! -f "$env_file" ]]; then
    echo "Environment file not found: $env_file" >&2
    echo "Set ATLAS_ENV_FILE to a valid file (for example environment.core.yml)." >&2
    exit 2
  fi
  if command -v micromamba >/dev/null 2>&1; then
    if ! micromamba env list | awk '{print $1}' | grep -qx "$env_name"; then
      micromamba create -y -n "$env_name" -f "$env_file"
    fi
    if [[ "$install_optional" == "1" ]]; then
      if [[ -f "environment.optional.yml" ]]; then
        micromamba env update -n "$env_name" -f environment.optional.yml
      else
        echo "ATLAS_INSTALL_OPTIONAL=1 requested but environment.optional.yml is missing." >&2
        exit 2
      fi
    fi
    # shellcheck disable=SC1091
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "$env_name"
  elif command -v conda >/dev/null 2>&1; then
    if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
      conda env create -n "$env_name" -f "$env_file"
    fi
    if [[ "$install_optional" == "1" ]]; then
      if [[ -f "environment.optional.yml" ]]; then
        conda env update -n "$env_name" -f environment.optional.yml
      else
        echo "ATLAS_INSTALL_OPTIONAL=1 requested but environment.optional.yml is missing." >&2
        exit 2
      fi
    fi
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$env_name"
  else
    echo "Neither micromamba nor conda was found on PATH." >&2
    echo "Install micromamba or rerun with ATLAS_SKIP_ENV_CREATE=1 inside an active Python 3.10 environment." >&2
    exit 2
  fi
fi

python -m pip install -e .

if [[ ! -f config.txt ]]; then
  cp config.example.txt config.txt
  echo "Created config.txt from config.example.txt"
fi

if [[ "$auto_config" == "1" ]]; then
  # Best-effort local path detection. Keep manual overrides in config.txt intact.
  # shellcheck disable=SC2086
  python tools/installers/configure_local_tools.py $configure_args
else
  echo "Config auto-detection skipped. Run atlas init later, or set ATLAS_AUTO_CONFIG=1."
fi

atlas --help >/dev/null
echo "Atlas CLI installed: $(command -v atlas)"

if [[ "$verify_tools" == "1" ]]; then
  atlas --verify-tools
else
  echo "External tool check skipped. Run ATLAS_VERIFY_TOOLS=1 bash tools/install_atlas_dev.sh or atlas --verify-tools when tool paths are configured."
fi
