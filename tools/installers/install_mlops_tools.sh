#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_mlops_tools.sh [--all] [--mlops] [--molecular] [--rl] [--generative] [--dry-run]
                              [--tools-root DIR] [--env-prefix-dir DIR] [--pip-cache-dir DIR]
                              [--conda-pkgs-dir DIR] [--cache-home DIR] [--home-dir DIR] [--skip-verify]

Install optional Atlas ML/MLOps helper environments outside the repository.
Defaults:
  tools root: $ATLAS_TOOLS_DIR or ../../tools relative to this repository
  envs:       $ATLAS_TOOLS_DIR/envs/atlas-ml-*

Tiers:
  --mlops       MLflow, Optuna, DVC, DVCLive, W&B, ClearML
  --molecular   Chemprop, PyTDC, DeepChem, RDKit
  --rl          Gymnasium, Stable-Baselines3, Ray Tune/RLlib
  --generative  REINVENT4, GuacaMol, MOSES source checkouts plus RDKit env

No tier selected means --mlops.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
atlas_root="$(cd "$repo_root/../.." && pwd)"
original_args="$*"
tools_root="${ATLAS_TOOLS_DIR:-$atlas_root/tools}"
env_prefix_dir="${ATLAS_ML_ENVS_DIR:-$tools_root/envs}"
pip_cache_dir="${ATLAS_ML_PIP_CACHE_DIR:-$tools_root/pip_cache}"
conda_pkgs_dir="${ATLAS_ML_CONDA_PKGS_DIR:-$tools_root/conda_pkgs}"
cache_home="${ATLAS_ML_CACHE_HOME:-$tools_root/cache}"
home_dir="${ATLAS_ML_HOME:-$cache_home/home}"
mamba_root_prefix="${ATLAS_ML_MAMBA_ROOT_PREFIX:-$tools_root/micromamba_root}"
mpl_config_dir="${ATLAS_ML_MPLCONFIGDIR:-$cache_home/matplotlib}"
dry_run="0"
verify="1"
install_mlops="0"
install_molecular="0"
install_rl="0"
install_generative="0"

mlops_package_spec="${ATLAS_MLOPS_PACKAGES:-mlflow optuna dvc[s3,ssh] dvclive wandb clearml}"
molecular_package_spec="${ATLAS_MOLECULAR_ML_PACKAGES:-chemprop PyTDC deepchem}"
rl_package_spec="${ATLAS_RL_PACKAGES:-gymnasium stable-baselines3 ray[tune,rllib]}"
read -r -a mlops_packages <<< "$mlops_package_spec"
read -r -a molecular_packages <<< "$molecular_package_spec"
read -r -a rl_packages <<< "$rl_package_spec"
generative_reinvent_url="${ATLAS_REINVENT4_GIT_URL:-https://github.com/MolecularAI/REINVENT4.git}"
generative_guacamol_url="${ATLAS_GUACAMOL_GIT_URL:-https://github.com/BenevolentAI/guacamol.git}"
generative_moses_url="${ATLAS_MOSES_GIT_URL:-https://github.com/molecularsets/moses.git}"
generative_reinvent_ref="${ATLAS_REINVENT4_GIT_REF:-}"
generative_guacamol_ref="${ATLAS_GUACAMOL_GIT_REF:-}"
generative_moses_ref="${ATLAS_MOSES_GIT_REF:-}"
git_lfs_skip_smudge="${ATLAS_ML_GIT_LFS_SKIP_SMUDGE:-1}"
generative_reinvent_torch_index="${ATLAS_REINVENT_TORCH_INDEX:-https://download.pytorch.org/whl/cpu}"
generative_install_reinvent="${ATLAS_GENERATIVE_INSTALL_REINVENT:-1}"
generative_install_moses="${ATLAS_GENERATIVE_INSTALL_MOSES:-1}"
generative_compat_package_spec="${ATLAS_GENERATIVE_COMPAT_PACKAGES:-numpy<2 pandas<2}"
read -r -a generative_compat_packages <<< "$generative_compat_package_spec"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all)
      install_mlops="1"
      install_molecular="1"
      install_rl="1"
      install_generative="1"
      shift
      ;;
    --mlops)
      install_mlops="1"
      shift
      ;;
    --molecular)
      install_molecular="1"
      shift
      ;;
    --rl)
      install_rl="1"
      shift
      ;;
    --generative)
      install_generative="1"
      shift
      ;;
    --dry-run)
      dry_run="1"
      shift
      ;;
    --skip-verify)
      verify="0"
      shift
      ;;
    --tools-root)
      tools_root="${2:?--tools-root requires a value}"
      shift 2
      ;;
    --env-prefix-dir)
      env_prefix_dir="${2:?--env-prefix-dir requires a value}"
      shift 2
      ;;
    --pip-cache-dir)
      pip_cache_dir="${2:?--pip-cache-dir requires a value}"
      shift 2
      ;;
    --conda-pkgs-dir)
      conda_pkgs_dir="${2:?--conda-pkgs-dir requires a value}"
      shift 2
      ;;
    --cache-home)
      cache_home="${2:?--cache-home requires a value}"
      if [[ -z "${ATLAS_ML_HOME:-}" ]]; then
        home_dir="$cache_home/home"
      fi
      mpl_config_dir="$cache_home/matplotlib"
      shift 2
      ;;
    --home-dir)
      home_dir="${2:?--home-dir requires a value}"
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

if [[ "$install_mlops$install_molecular$install_rl$install_generative" == "0000" ]]; then
  install_mlops="1"
fi

run_cmd() {
  echo "+ $*"
  if [[ "$dry_run" == "0" ]]; then
    "$@"
  fi
}

managed_cache_cmd() {
  if [[ "$dry_run" == "1" ]]; then
    echo "+ HOME=$home_dir MAMBA_ROOT_PREFIX=$mamba_root_prefix MAMBA_PKGS_DIRS=$conda_pkgs_dir CONDA_PKGS_DIRS=$conda_pkgs_dir XDG_CACHE_HOME=$cache_home MPLCONFIGDIR=$mpl_config_dir $*"
    return 0
  fi
  HOME="$home_dir" \
    MAMBA_ROOT_PREFIX="$mamba_root_prefix" \
    MAMBA_PKGS_DIRS="$conda_pkgs_dir" \
    CONDA_PKGS_DIRS="$conda_pkgs_dir" \
    XDG_CACHE_HOME="$cache_home" \
    MPLCONFIGDIR="$mpl_config_dir" \
    "$@"
}

git_cmd() {
  if [[ "$dry_run" == "1" ]]; then
    echo "+ GIT_LFS_SKIP_SMUDGE=$git_lfs_skip_smudge git $*"
    return 0
  fi
  GIT_LFS_SKIP_SMUDGE="$git_lfs_skip_smudge" git "$@"
}

manager() {
  if [[ "$dry_run" == "1" ]]; then
    echo "${ATLAS_ENV_MANAGER:-micromamba}"
    return 0
  fi
  if command -v micromamba >/dev/null 2>&1; then
    echo micromamba
  elif command -v conda >/dev/null 2>&1; then
    echo conda
  else
    echo "Neither micromamba nor conda was found on PATH." >&2
    exit 2
  fi
}

create_env() {
  local prefix="$1"
  shift
  local mgr
  mgr="$(manager)"
  run_cmd mkdir -p "$(dirname "$prefix")" "$conda_pkgs_dir" "$cache_home" "$home_dir" "$mamba_root_prefix" "$mpl_config_dir"
  if [[ "$mgr" == "micromamba" ]]; then
    if [[ -d "$prefix" ]]; then
      managed_cache_cmd micromamba install -y -p "$prefix" -c conda-forge "$@"
    else
      managed_cache_cmd micromamba create -y -p "$prefix" -c conda-forge "$@"
    fi
  elif [[ -d "$prefix" ]]; then
    managed_cache_cmd conda install -y -p "$prefix" -c conda-forge "$@"
  else
    managed_cache_cmd conda create -y -p "$prefix" -c conda-forge "$@"
  fi
}

pip_install() {
  local prefix="$1"
  shift
  run_cmd mkdir -p "$pip_cache_dir" "$cache_home" "$home_dir" "$mpl_config_dir"
  if [[ "$dry_run" == "1" ]]; then
    echo "+ HOME=$home_dir PIP_CACHE_DIR=$pip_cache_dir XDG_CACHE_HOME=$cache_home MPLCONFIGDIR=$mpl_config_dir $prefix/bin/python -m pip install --upgrade pip"
    echo "+ HOME=$home_dir PIP_CACHE_DIR=$pip_cache_dir XDG_CACHE_HOME=$cache_home MPLCONFIGDIR=$mpl_config_dir $prefix/bin/python -m pip install $*"
    return 0
  fi
  HOME="$home_dir" \
    PIP_CACHE_DIR="$pip_cache_dir" \
    XDG_CACHE_HOME="$cache_home" \
    MPLCONFIGDIR="$mpl_config_dir" \
    "$prefix/bin/python" -m pip install --upgrade pip
  HOME="$home_dir" \
    PIP_CACHE_DIR="$pip_cache_dir" \
    XDG_CACHE_HOME="$cache_home" \
    MPLCONFIGDIR="$mpl_config_dir" \
    "$prefix/bin/python" -m pip install "$@"
}

pip_install_requirements() {
  local prefix="$1"
  local requirements="$2"
  run_cmd mkdir -p "$pip_cache_dir" "$cache_home" "$home_dir" "$mpl_config_dir"
  if [[ "$dry_run" == "1" ]]; then
    echo "+ HOME=$home_dir PIP_CACHE_DIR=$pip_cache_dir XDG_CACHE_HOME=$cache_home MPLCONFIGDIR=$mpl_config_dir $prefix/bin/python -m pip install -r $requirements"
    return 0
  fi
  HOME="$home_dir" \
    PIP_CACHE_DIR="$pip_cache_dir" \
    XDG_CACHE_HOME="$cache_home" \
    MPLCONFIGDIR="$mpl_config_dir" \
    "$prefix/bin/python" -m pip install -r "$requirements"
}

pip_install_project_editable() {
  local prefix="$1"
  local project_dir="$2"
  shift 2
  run_cmd mkdir -p "$pip_cache_dir" "$cache_home" "$home_dir" "$mpl_config_dir"
  if [[ "$dry_run" == "1" ]]; then
    echo "+ HOME=$home_dir PIP_CACHE_DIR=$pip_cache_dir XDG_CACHE_HOME=$cache_home MPLCONFIGDIR=$mpl_config_dir $prefix/bin/python -m pip install -e $project_dir $*"
    return 0
  fi
  HOME="$home_dir" \
    PIP_CACHE_DIR="$pip_cache_dir" \
    XDG_CACHE_HOME="$cache_home" \
    MPLCONFIGDIR="$mpl_config_dir" \
    "$prefix/bin/python" -m pip install -e "$project_dir" "$@"
}

verify_imports() {
  local prefix="$1"
  local imports="$2"
  if [[ "$verify" != "1" ]]; then
    return 0
  fi
  managed_cache_cmd "$prefix/bin/python" -c "$imports"
}

clone_or_update() {
  local url="$1"
  local dest="$2"
  local ref="$3"
  run_cmd mkdir -p "$(dirname "$dest")"
  if [[ -d "$dest/.git" ]]; then
    git_cmd -C "$dest" fetch --tags origin
  elif [[ ! -e "$dest" ]]; then
    git_cmd clone "$url" "$dest"
  else
    echo "Existing non-git path left untouched: $dest" >&2
  fi
  if [[ -n "$ref" && -d "$dest/.git" ]]; then
    git_cmd -C "$dest" checkout "$ref"
  fi
}

capture_env_state() {
  local name="$1"
  local prefix="$2"
  local state_dir="$tools_root/ml_tool_env_state"
  if [[ "$dry_run" == "1" ]]; then
    echo "+ write env state for $name under $state_dir"
    return 0
  fi
  mkdir -p "$state_dir"
  HOME="$home_dir" \
    XDG_CACHE_HOME="$cache_home" \
    MPLCONFIGDIR="$mpl_config_dir" \
    "$prefix/bin/python" -m pip freeze --all > "$state_dir/$name-pip-freeze.txt" 2>/dev/null || true
  managed_cache_cmd "$(manager)" list -p "$prefix" --explicit > "$state_dir/$name-conda-explicit.txt" 2>/dev/null || true
}

capture_git_state() {
  local name="$1"
  local path="$2"
  local state_dir="$tools_root/ml_tool_env_state"
  if [[ "$dry_run" == "1" ]]; then
    echo "+ write git state for $name under $state_dir"
    return 0
  fi
  if [[ -d "$path/.git" ]]; then
    mkdir -p "$state_dir"
    {
      echo "path=$path"
      echo "remote=$(git -C "$path" remote get-url origin 2>/dev/null || true)"
      echo "head=$(git -C "$path" rev-parse HEAD 2>/dev/null || true)"
      echo "describe=$(git -C "$path" describe --tags --always --dirty 2>/dev/null || true)"
    } > "$state_dir/$name-git.txt"
  fi
}

write_manifest() {
  local manifest="$tools_root/ml_tools_manifest.txt"
  local history="$tools_root/ml_tools_install_history.tsv"
  run_cmd mkdir -p "$tools_root"
  if [[ "$dry_run" == "1" ]]; then
    echo "+ write $manifest"
    return 0
  fi
  {
    echo "atlas_ml_tools_installed_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "tools_root=$tools_root"
    echo "env_prefix_dir=$env_prefix_dir"
    echo "pip_cache_dir=$pip_cache_dir"
    echo "conda_pkgs_dir=$conda_pkgs_dir"
    echo "cache_home=$cache_home"
    echo "home_dir=$home_dir"
    echo "mamba_root_prefix=$mamba_root_prefix"
    echo "mpl_config_dir=$mpl_config_dir"
    echo "state_dir=$tools_root/ml_tool_env_state"
    echo "git_lfs_skip_smudge=$git_lfs_skip_smudge"
    echo "generative_reinvent_torch_index=$generative_reinvent_torch_index"
    echo "generative_install_reinvent=$generative_install_reinvent"
    echo "generative_install_moses=$generative_install_moses"
    echo "generative_compat_packages=$generative_compat_package_spec"
    echo "invocation=$original_args"
    echo "mlops=$install_mlops"
    echo "molecular=$install_molecular"
    echo "rl=$install_rl"
    echo "generative=$install_generative"
  } > "$manifest"
  printf '%s\tmlops=%s\tmolecular=%s\trl=%s\tgenerative=%s\targs=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "$install_mlops" \
    "$install_molecular" \
    "$install_rl" \
    "$install_generative" \
    "$original_args" >> "$history"
  echo "ML tools manifest: $manifest"
}

cd "$repo_root"

if [[ "$install_mlops" == "1" ]]; then
  mlops_prefix="$env_prefix_dir/atlas-ml-mlops"
  create_env "$mlops_prefix" python=3.11 pip
  pip_install "$mlops_prefix" "${mlops_packages[@]}"
  verify_imports "$mlops_prefix" "import mlflow, optuna, dvc, dvclive, wandb, clearml"
  capture_env_state "atlas-ml-mlops" "$mlops_prefix"
  echo "Atlas MLOps Python: $mlops_prefix/bin/python"
fi

if [[ "$install_molecular" == "1" ]]; then
  molecular_prefix="$env_prefix_dir/atlas-ml-molecular"
  create_env "$molecular_prefix" python=3.11 pip rdkit
  pip_install "$molecular_prefix" "${molecular_packages[@]}"
  verify_imports "$molecular_prefix" "import chemprop, deepchem; import tdc"
  capture_env_state "atlas-ml-molecular" "$molecular_prefix"
  echo "Atlas molecular ML Python: $molecular_prefix/bin/python"
fi

if [[ "$install_rl" == "1" ]]; then
  rl_prefix="$env_prefix_dir/atlas-ml-rl"
  create_env "$rl_prefix" python=3.11 pip
  pip_install "$rl_prefix" "${rl_packages[@]}"
  verify_imports "$rl_prefix" "import gymnasium, stable_baselines3, ray"
  capture_env_state "atlas-ml-rl" "$rl_prefix"
  echo "Atlas RL Python: $rl_prefix/bin/python"
fi

if [[ "$install_generative" == "1" ]]; then
  generative_prefix="$env_prefix_dir/atlas-ml-generative"
  generative_src="$tools_root/ml/generative"
  create_env "$generative_prefix" python=3.11 pip rdkit
  clone_or_update "$generative_reinvent_url" "$generative_src/REINVENT4" "$generative_reinvent_ref"
  clone_or_update "$generative_guacamol_url" "$generative_src/guacamol" "$generative_guacamol_ref"
  clone_or_update "$generative_moses_url" "$generative_src/moses" "$generative_moses_ref"
  pip_install "$generative_prefix" guacamol
  if [[ -f "$generative_src/REINVENT4/requirements.txt" ]]; then
    pip_install_requirements "$generative_prefix" "$generative_src/REINVENT4/requirements.txt"
  fi
  if [[ "$generative_install_reinvent" == "1" && -f "$generative_src/REINVENT4/pyproject.toml" ]]; then
    pip_install_project_editable "$generative_prefix" "$generative_src/REINVENT4" --extra-index-url "$generative_reinvent_torch_index"
  fi
  if [[ -f "$generative_src/moses/requirements.txt" ]]; then
    pip_install_requirements "$generative_prefix" "$generative_src/moses/requirements.txt"
  fi
  if [[ "$generative_install_moses" == "1" && -f "$generative_src/moses/setup.py" ]]; then
    pip_install_project_editable "$generative_prefix" "$generative_src/moses"
  fi
  if [[ ${#generative_compat_packages[@]} -gt 0 ]]; then
    pip_install "$generative_prefix" "${generative_compat_packages[@]}"
  fi
  generative_verify_imports="import guacamol"
  if [[ "$generative_install_reinvent" == "1" ]]; then
    generative_verify_imports="$generative_verify_imports, reinvent"
  fi
  if [[ "$generative_install_moses" == "1" ]]; then
    generative_verify_imports="$generative_verify_imports, moses"
  fi
  verify_imports "$generative_prefix" "$generative_verify_imports"
  capture_env_state "atlas-ml-generative" "$generative_prefix"
  capture_git_state "REINVENT4" "$generative_src/REINVENT4"
  capture_git_state "guacamol" "$generative_src/guacamol"
  capture_git_state "moses" "$generative_src/moses"
  echo "Atlas generative ML Python: $generative_prefix/bin/python"
  echo "Generative source checkouts: $generative_src"
fi

write_manifest
