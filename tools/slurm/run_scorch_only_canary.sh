#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: tools/slurm/run_scorch_only_canary.sh RUN_ID --all-dirs PATH [--pdb-id PDB]

Runs SCORCH rescoring only inside an existing Slurm/idev allocation. This reuses
existing docking outputs and does not launch Atlas docking.

Overrides:
  ATLAS_IDEV_NODES, ATLAS_IDEV_TASKS, ATLAS_IDEV_TASKS_PER_NODE
  ATLAS_IDEV_CPUS_PER_TASK, SCORCH_TOP_FRACTION, ATLAS_SCORCH_CHUNK_SIZE
  SCORCH_PYTHON, SCORCH_ENV_PREFIX, SCORCH_ENV
EOF
}

if [[ $# -lt 1 ]]; then
  usage
  exit 2
fi

RUN_ID="$1"
shift

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ALL_DIRS_ROOT="${ATLAS_ALL_DIRS:-}"
PDB_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all-dirs)
      ALL_DIRS_ROOT="${2:-}"
      shift 2
      ;;
    --pdb-id)
      PDB_ARGS+=(--pdb-id "${2:-}")
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[atlas-scorch-only] unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "${RUN_ID}" || -z "${ALL_DIRS_ROOT}" ]]; then
  usage
  exit 2
fi

resolve_scorch_python() {
  local env_name="${SCORCH_ENV:-scorch-env}"
  local candidate
  local candidates=()
  if [[ -n "${SCORCH_PYTHON:-}" && -x "${SCORCH_PYTHON}" ]]; then
    echo "${SCORCH_PYTHON}"
    return 0
  fi
  if [[ -n "${SCORCH_ENV_PREFIX:-}" ]]; then
    candidates+=(
      "${SCORCH_ENV_PREFIX}/bin/python"
      "${SCORCH_ENV_PREFIX}/bin/python3.7"
    )
  fi
  if [[ -n "${MAMBA_ROOT_PREFIX:-}" ]]; then
    candidates+=(
      "${MAMBA_ROOT_PREFIX}/envs/${env_name}/bin/python"
      "${MAMBA_ROOT_PREFIX}/envs/${env_name}/bin/python3.7"
    )
  fi
  candidates+=(
    "${REPO_ROOT}/../../../.mamba/envs/${env_name}/bin/python"
    "${REPO_ROOT}/../../../.mamba/envs/${env_name}/bin/python3.7"
    "${HOME}/.mamba/envs/${env_name}/bin/python"
    "${HOME}/.mamba/envs/${env_name}/bin/python3.7"
  )
  for candidate in "${candidates[@]}"; do
    if [[ -n "${candidate}" && -x "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

SCORCH_PY="$(resolve_scorch_python || true)"
if [[ -z "${SCORCH_PY}" ]]; then
  echo "[atlas-scorch-only] could not resolve scorch python; set SCORCH_PYTHON or SCORCH_ENV_PREFIX" >&2
  exit 2
fi

NODES="${ATLAS_IDEV_NODES:-${SLURM_NNODES:-1}}"
TASKS="${ATLAS_IDEV_TASKS:-${SLURM_NTASKS:-48}}"
if [[ -n "${ATLAS_IDEV_TASKS_PER_NODE:-}" ]]; then
  TASKS_PER_NODE="${ATLAS_IDEV_TASKS_PER_NODE}"
elif [[ "${NODES}" =~ ^[0-9]+$ && "${TASKS}" =~ ^[0-9]+$ && "${NODES}" -gt 0 ]]; then
  TASKS_PER_NODE="$((TASKS / NODES))"
else
  TASKS_PER_NODE="${TASKS}"
fi
CPUS_PER_TASK="${ATLAS_IDEV_CPUS_PER_TASK:-1}"
SLURM_LOG_ROOT="${ATLAS_SLURM_LOG_ROOT:-${ALL_DIRS_ROOT}/outputs/logs/slurm}"
TAG="${ATLAS_SCORCH_ONLY_TAG:-scorch_only_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "${SLURM_LOG_ROOT}"

export ATLAS_ALL_DIRS="${ALL_DIRS_ROOT}"
export PREPPED_LIGANDS_DIR="${ALL_DIRS_ROOT}/prepped_ligands"
export EXTRACTED_LIGANDS_DIR="${ALL_DIRS_ROOT}/extracted_ligands"
export DOCKED_DIR="${ALL_DIRS_ROOT}/outputs/docked"
export POST_DOCKED_DIR="${ALL_DIRS_ROOT}/outputs/post_docked"
export OUTPUT_DIR="${ALL_DIRS_ROOT}/outputs/processed_pdbs"
export CONFIGS_DIR="${ALL_DIRS_ROOT}/outputs/configs"
export MANIFESTS_DIR="${ALL_DIRS_ROOT}/outputs/manifests"
export LOGS_DIR="${ALL_DIRS_ROOT}/outputs/logs"
export DATA_DIR="${ALL_DIRS_ROOT}/outputs/data"
export USE_SCORCH=true
export SCORCH_TOP_FRACTION="${SCORCH_TOP_FRACTION:-1.0}"
export ATLAS_SCORCH_SHARDS_ENABLE="${ATLAS_SCORCH_SHARDS_ENABLE:-1}"
export ATLAS_SCORCH_SHARD_LOCAL_WORKERS="${ATLAS_SCORCH_SHARD_LOCAL_WORKERS:-1}"
export ATLAS_SCORCH_CHUNK_SIZE="${ATLAS_SCORCH_CHUNK_SIZE:-8}"
export ATLAS_SCORCH_CHUNK_MIN="${ATLAS_SCORCH_CHUNK_MIN:-1}"
export ATLAS_SCORCH_PARALLEL_PROFILE="${ATLAS_SCORCH_PARALLEL_PROFILE:-conservative}"
export CPU="${CPU:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

echo "[atlas-scorch-only] run_id=${RUN_ID} all_dirs=${ALL_DIRS_ROOT} job=${SLURM_JOB_ID:-none} nodes=${NODES} tasks=${TASKS} tasks_per_node=${TASKS_PER_NODE} cpus_per_task=${CPUS_PER_TASK}"
echo "[atlas-scorch-only] scorch_python=${SCORCH_PY}"
echo "[atlas-scorch-only] out=${SLURM_LOG_ROOT}/${RUN_ID}_${TAG}.out err=${SLURM_LOG_ROOT}/${RUN_ID}_${TAG}.err"

srun \
  --exclusive \
  --nodes="${NODES}" \
  --ntasks="${TASKS}" \
  --ntasks-per-node="${TASKS_PER_NODE}" \
  --cpus-per-task="${CPUS_PER_TASK}" \
  bash -lc '
    set -euo pipefail
    export PYTHONPATH="'"${REPO_ROOT}"'/src:${PYTHONPATH:-}"
    export CPU="${CPU:-1}"
    export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
    export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
    export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
    unset LD_PRELOAD
    unset LD_AUDIT
    SCORCH_LIB="$(cd "$(dirname "'"${SCORCH_PY}"'")/.." && pwd)/lib"
    export LD_LIBRARY_PATH="${SCORCH_LIB}:${LD_LIBRARY_PATH:-}"
    export ATLAS_DISTRIBUTED_MODE=slurm_array
    export SLURM_ARRAY_TASK_ID="${SLURM_PROCID:-0}"
    export SLURM_ARRAY_TASK_COUNT="${SLURM_NTASKS:-1}"
    export SLURM_ARRAY_TASK_MIN=0
    export SLURM_ARRAY_TASK_MAX=$((SLURM_ARRAY_TASK_COUNT - 1))
    echo "[atlas-scorch-only-task] host=$(hostname) task=${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} cpu=${CPU}"
    "'"${SCORCH_PY}"'" -m post_docking.rescoring.rescoring_scorch \
      --run-id "'"${RUN_ID}"'" \
      --repo-root "'"${REPO_ROOT}"'" \
      --docked-root "${DOCKED_DIR}" \
      --post-docked-root "${POST_DOCKED_DIR}" \
      --threads 1 \
      --jobs 1 \
      --decoy-prefix decoys \
      --verbose \
      '"${PDB_ARGS[*]}"'
  ' >"${SLURM_LOG_ROOT}/${RUN_ID}_${TAG}.out" \
    2>"${SLURM_LOG_ROOT}/${RUN_ID}_${TAG}.err"
