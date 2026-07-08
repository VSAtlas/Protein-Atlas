#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_scorch.sh [--source-dir DIR | --git-url URL] [--prefix DIR] [--env-name NAME] [--env-prefix DIR] [--configure]

Install or register SCORCH and create the separate scorch-env environment used
by Atlas rescoring. Network access is used only when --git-url is supplied.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
atlas_root="$(cd "$repo_root/../.." && pwd)"
tools_root="${ATLAS_TOOLS_DIR:-$atlas_root/tools}"
source_dir="${SCORCH_SOURCE_DIR:-}"
git_url="${SCORCH_GIT_URL:-}"
prefix="${SCORCH_PREFIX:-$tools_root/SCORCH}"
env_name="${SCORCH_ENV:-scorch-env}"
env_prefix="${SCORCH_ENV_PREFIX:-$tools_root/envs/scorch-env}"
configure="${ATLAS_CONFIGURE_TOOLS:-0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-dir)
      source_dir="${2:?--source-dir requires a value}"
      shift 2
      ;;
    --git-url)
      git_url="${2:?--git-url requires a value}"
      shift 2
      ;;
    --prefix)
      prefix="${2:?--prefix requires a value}"
      shift 2
      ;;
    --env-name)
      env_name="${2:?--env-name requires a value}"
      shift 2
      ;;
    --env-prefix)
      env_prefix="${2:?--env-prefix requires a value}"
      shift 2
      ;;
    --configure)
      configure="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

cd "$repo_root"
mkdir -p "$(dirname "$prefix")" "$(dirname "$env_prefix")"

if [[ -n "$source_dir" ]]; then
  mkdir -p "$prefix"
  cp -a "$source_dir"/. "$prefix"/
elif [[ -n "$git_url" ]]; then
  if [[ -e "$prefix" ]]; then
    echo "SCORCH prefix already exists: $prefix" >&2
    echo "Remove it or choose --prefix." >&2
    exit 2
  fi
  git clone "$git_url" "$prefix"
elif [[ -d "$prefix" ]]; then
  :
else
  echo "Provide --source-dir DIR or --git-url URL, or set SCORCH_PREFIX to an existing checkout." >&2
  exit 2
fi

scorch_script="$(find "$prefix" -maxdepth 4 -type f -name scorch.py 2>/dev/null | head -n 1 || true)"
if [[ -z "$scorch_script" ]]; then
  echo "Could not locate scorch.py under $prefix" >&2
  exit 2
fi

if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  if [[ -n "$env_prefix" ]]; then
    if [[ ! -d "$env_prefix" ]]; then
      micromamba create -y -p "$env_prefix" python=3.10 pip
    fi
    micromamba activate "$env_prefix"
  else
    if ! micromamba env list | awk '{print $1}' | grep -qx "$env_name"; then
      micromamba create -y -n "$env_name" python=3.10 pip
    fi
    micromamba activate "$env_name"
  fi
elif command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ -n "$env_prefix" ]]; then
    if [[ ! -d "$env_prefix" ]]; then
      conda create -y -p "$env_prefix" python=3.10 pip
    fi
    conda activate "$env_prefix"
  else
    if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
      conda create -y -n "$env_name" python=3.10 pip
    fi
    conda activate "$env_name"
  fi
else
  echo "Neither micromamba nor conda was found on PATH." >&2
  exit 2
fi

if [[ -f "$prefix/requirements.txt" ]]; then
  python -m pip install -r "$prefix/requirements.txt"
fi
if [[ -f "$prefix/pyproject.toml" || -f "$prefix/setup.py" ]]; then
  python -m pip install -e "$prefix"
fi

echo "SCORCH script: $scorch_script"
echo "SCORCH env: ${env_prefix:-$env_name}"

if [[ "$configure" == "1" ]]; then
  args=(--set "SCORCH=$scorch_script" --set "SCORCH_ENV=$env_name")
  if [[ -n "$env_prefix" ]]; then
    args+=(--set "SCORCH_ENV_PREFIX=$env_prefix")
  fi
  tools/installers/configure_local_tools.py "${args[@]}"
fi
