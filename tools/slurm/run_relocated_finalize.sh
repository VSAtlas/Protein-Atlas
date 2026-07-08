#!/usr/bin/env bash
#SBATCH -J atlas_relocated_finalize
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 04:00:00
#SBATCH -p spr
#SBATCH -A MCB24041
#SBATCH -o logs/slurm/%x_%j.out
#SBATCH -e logs/slurm/%x_%j.err

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

RUN_ID="${RUN_ID:-atlas_relocated_${SLURM_JOB_ID:-local}}"

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

export PREPPED_LIGANDS_DIR="${ATLAS_PREPPED_ROOT:-${ATLAS_ALL_DIRS}/prepped_ligands}"
export EXTRACTED_LIGANDS_DIR="${ATLAS_ALL_DIRS}/extracted_ligands"
export DOCKED_DIR="${ATLAS_ALL_DIRS}/outputs/docked"
export POST_DOCKED_DIR="${ATLAS_ALL_DIRS}/outputs/post_docked"
export OUTPUT_DIR="${ATLAS_ALL_DIRS}/outputs/processed_pdbs"
export CONFIGS_DIR="${ATLAS_ALL_DIRS}/outputs/configs"
export MANIFESTS_DIR="${ATLAS_ALL_DIRS}/outputs/manifests"
export LOGS_DIR="${ATLAS_ALL_DIRS}/outputs/logs"
export DATA_DIR="${ATLAS_ALL_DIRS}/outputs/data"

echo "[atlas-relocated-finalize] run_id=${RUN_ID} all_dirs=${ATLAS_ALL_DIRS} prepped_root=${PREPPED_LIGANDS_DIR}"

python tools/finalize_distributed_run.py --run-id "${RUN_ID}"
