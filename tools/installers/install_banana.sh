#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_banana.sh [--git-url URL | --source-dir DIR | --existing-root DIR] [--git-ref REF] [--prefix DIR] [--env-prefix DIR]

Install or register the BANANA inference checkout and its isolated Python
environment. Defaults install under $ATLAS_TOOLS_DIR/banana and
$ATLAS_TOOLS_DIR/envs/banana.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
atlas_root="$(cd "$repo_root/../.." && pwd)"
tools_root="${ATLAS_TOOLS_DIR:-$atlas_root/tools}"
prefix="${BANANA_ROOT:-$tools_root/banana}"
env_prefix="${BANANA_ENV_PREFIX:-$tools_root/envs/banana}"
source_dir="${BANANA_SOURCE_DIR:-}"
existing_root="${BANANA_EXISTING_ROOT:-}"
git_url="${BANANA_GIT_URL:-https://github.com/molecularmodelinglab/banana.git}"
git_ref="${BANANA_GIT_REF:-06b3ee49894ec7c73e3f835cca31d2d5bae5763f}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-dir)
      source_dir="${2:?--source-dir requires a value}"
      shift 2
      ;;
    --existing-root)
      existing_root="${2:?--existing-root requires a value}"
      shift 2
      ;;
    --git-url)
      git_url="${2:?--git-url requires a value}"
      shift 2
      ;;
    --git-ref)
      git_ref="${2:?--git-ref requires a value}"
      shift 2
      ;;
    --prefix)
      prefix="${2:?--prefix requires a value}"
      shift 2
      ;;
    --env-prefix)
      env_prefix="${2:?--env-prefix requires a value}"
      shift 2
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

if [[ -n "$existing_root" ]]; then
  prefix="$existing_root"
elif [[ -n "$source_dir" ]]; then
  mkdir -p "$prefix"
  cp -a "$source_dir"/. "$prefix"/
elif [[ -d "$prefix" ]]; then
  :
else
  git clone "$git_url" "$prefix"
fi

if [[ ! -f "$prefix/inference.py" ]]; then
  echo "BANANA inference.py not found under $prefix" >&2
  exit 2
fi

if [[ -d "$prefix/.git" && -n "$git_ref" ]]; then
  git -C "$prefix" fetch --tags origin || true
  git -C "$prefix" checkout "$git_ref"
fi

if command -v micromamba >/dev/null 2>&1; then
  if [[ -d "$env_prefix" ]]; then
    micromamba install -y -p "$env_prefix" -c conda-forge python=3.10 pip rdkit
  else
    micromamba create -y -p "$env_prefix" -c conda-forge python=3.10 pip rdkit
  fi
elif command -v conda >/dev/null 2>&1; then
  if [[ -d "$env_prefix" ]]; then
    conda install -y -p "$env_prefix" -c conda-forge python=3.10 pip rdkit
  else
    conda create -y -p "$env_prefix" -c conda-forge python=3.10 pip rdkit
  fi
else
  echo "Neither micromamba nor conda was found on PATH." >&2
  exit 2
fi

python_bin="$env_prefix/bin/python"
if [[ ! -x "$python_bin" ]]; then
  echo "BANANA Python not found: $python_bin" >&2
  exit 3
fi

"$python_bin" -m pip install --upgrade pip
"$python_bin" -m pip install torch dgl
if [[ -f "$prefix/requirements.txt" ]]; then
  "$python_bin" -m pip install -r "$prefix/requirements.txt"
fi

"$python_bin" -c "import torch, dgl, rdkit, dgllife, omegaconf" >/dev/null

echo "BANANA root: $prefix"
echo "BANANA Python: $python_bin"
echo "Use with: BANANA_ROOT=$prefix BANANA_PYTHON=$python_bin"
