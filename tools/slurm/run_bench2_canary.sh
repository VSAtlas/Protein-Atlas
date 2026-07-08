#!/usr/bin/env bash
#SBATCH -J atlas_bench2_canary
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH -t 00:15:00
#SBATCH --array=0-2%3
#SBATCH -o logs/slurm/bench2_canary_%A_%a.out
#SBATCH -e logs/slurm/bench2_canary_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${ATLAS_REPO_ROOT:-}" ]]; then
  REPO_ROOT="$(cd "${ATLAS_REPO_ROOT}" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/src/cli" ]]; then
  REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
  REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  cd "${REPO_ROOT}"
  export PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
  python -m cli.atlas_main_cli slurm submit --bench2-canary "$@"
  exit $?
fi

eval "$(micromamba shell hook --shell bash)"
RUN_ROOT="${ATLAS_REPO_ROOT:-${REPO_ROOT}}"
cd "${RUN_ROOT}"
micromamba activate docking-env
cd "${RUN_ROOT}"
mkdir -p "${RUN_ROOT}/logs/slurm"

RUN_ID="${RUN_ID:-bench2_canary_${SLURM_ARRAY_JOB_ID:-local}}"
if [[ -z "${ATLAS_USE_NODE_CPUS+x}" && -z "${SLURM_CPUS_PER_TASK:-}" ]]; then
  export ATLAS_USE_NODE_CPUS=1
fi

echo "[atlas-bench2-canary] host=$(hostname) pwd=$(pwd) run_id=${RUN_ID} task=${SLURM_ARRAY_TASK_ID:-local}/${SLURM_ARRAY_TASK_COUNT:-1} cpus_per_task=${SLURM_CPUS_PER_TASK:-} cpus_on_node=${SLURM_CPUS_ON_NODE:-} job_cpus_per_node=${SLURM_JOB_CPUS_PER_NODE:-}"

python tools/run_relocated_mode.py \
  --skip-extract \
  --all-dirs "${ATLAS_ALL_DIRS:-${REPO_ROOT}}" \
  --main-args "${ATLAS_MAIN_ARGS:--bench2 -fast --run-id ${RUN_ID}}"
