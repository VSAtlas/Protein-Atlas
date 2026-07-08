#!/usr/bin/env bash
#SBATCH -J atlas_relocated_spr
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH -t 10:00:00
#SBATCH -p spr
#SBATCH -A MCB24041
#SBATCH -o logs/slurm/%x_%A_%a.out
#SBATCH -e logs/slurm/%x_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${ATLAS_REPO_ROOT:-}" ]]; then
  REPO_ROOT="$(cd "${ATLAS_REPO_ROOT}" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/src/cli" ]]; then
  REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
  REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi

eval "$(micromamba shell hook --shell bash)"
RUN_ROOT="${ATLAS_REPO_ROOT:-${REPO_ROOT}}"
cd "${RUN_ROOT}"
micromamba activate docking-env
mkdir -p "${RUN_ROOT}/logs/slurm"

: "${ATLAS_ALL_DIRS:?Set ATLAS_ALL_DIRS to the run-specific relocated output root.}"

RUN_ID="${RUN_ID:-atlas_relocated_${SLURM_ARRAY_JOB_ID:-local}}"
if [[ -z "${ATLAS_USE_NODE_CPUS+x}" && -z "${SLURM_CPUS_PER_TASK:-}" ]]; then
  export ATLAS_USE_NODE_CPUS=1
fi

has_prepped_payload() {
  local root="$1"
  [[ -n "${root}" && -d "${root}" ]] || return 1
  find "${root}" -type f \( -name '_manifest.json' -o -name '*.pdbqt' \) -print -quit 2>/dev/null | grep -q .
}

resolve_prepped_root() {
  local candidate
  local -a candidates=()
  if [[ -n "${ATLAS_PREPPED_ROOT:-}" ]]; then
    candidates+=("${ATLAS_PREPPED_ROOT}")
  fi
  candidates+=(
    "${ATLAS_ALL_DIRS}/prepped_ligands"
    "${REPO_ROOT}/prepped_ligands"
  )
  if [[ -n "${SCRATCH:-}" ]]; then
    candidates+=("${SCRATCH}/atlas/code/protein_automation/prepped_ligands")
  fi
  if [[ "${REPO_ROOT}" =~ ^/work[0-9]*/([^/]+)/([^/]+)/(.*)$ ]]; then
    candidates+=("/scratch/${BASH_REMATCH[1]}/${BASH_REMATCH[2]}/${BASH_REMATCH[3]}/prepped_ligands")
  fi
  for candidate in "${candidates[@]}"; do
    if has_prepped_payload "${candidate}"; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

if [[ -z "${ATLAS_PREPPED_ROOT:-}" ]]; then
  AUTO_PREPPED_ROOT="$(resolve_prepped_root || true)"
  if [[ -n "${AUTO_PREPPED_ROOT}" ]]; then
    export ATLAS_PREPPED_ROOT="${AUTO_PREPPED_ROOT}"
  fi
fi

PREPPED_ARGS=()
if [[ -n "${ATLAS_PREPPED_ROOT:-}" ]]; then
  PREPPED_ARGS+=(--prepped-root "${ATLAS_PREPPED_ROOT}")
fi

EXTRACT_ARGS=()
if [[ "${ATLAS_SKIP_EXTRACT:-1}" == "1" || "${ATLAS_SKIP_EXTRACT:-1}" == "true" ]]; then
  EXTRACT_ARGS+=(--skip-extract)
fi

echo "[atlas-relocated-spr] host=$(hostname) run_id=${RUN_ID} task=${SLURM_ARRAY_TASK_ID:-local}/${SLURM_ARRAY_TASK_COUNT:-1} cpus_per_task=${SLURM_CPUS_PER_TASK:-auto} cpus_on_node=${SLURM_CPUS_ON_NODE:-} all_dirs=${ATLAS_ALL_DIRS} prepped_root=${ATLAS_PREPPED_ROOT:-<all-dirs>/prepped_ligands}"

python tools/run_relocated_mode.py \
  --all-dirs "${ATLAS_ALL_DIRS}" \
  "${PREPPED_ARGS[@]}" \
  "${EXTRACT_ARGS[@]}" \
  --main-args "${ATLAS_MAIN_ARGS:---resume --run-id ${RUN_ID}}"
