#!/usr/bin/env bash
#SBATCH -J atlas_spr
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH -t 12:00:00
#SBATCH -p spr
#SBATCH -A MCB24041
#SBATCH --array=0-12
#SBATCH -o logs/slurm/dud_e_%a.out
#SBATCH -e logs/slurm/dud_e_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${ATLAS_REPO_ROOT:-}" ]]; then
  REPO_ROOT="$(cd "${ATLAS_REPO_ROOT}" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/src/cli" ]]; then
  REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
  REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
RUN_TAG="${RUN_TAG:-sprDUD_24}"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  RUN_ID="${RUN_ID:-dud24_$(date +%Y%m%d_%H%M%S)}"

  ARRAY_JOB_RAW="$(
      sbatch \
      --parsable \
      --chdir="${REPO_ROOT}" \
      --export=ALL,RUN_ID="${RUN_ID}",ATLAS_REPO_ROOT="${REPO_ROOT}" \
      "${SCRIPT_DIR}/run_spr.sh"
  )"
  ARRAY_JOB_ID="${ARRAY_JOB_RAW%%;*}"

  FINAL_JOB_RAW="$(
    sbatch \
      --parsable \
      --chdir="${REPO_ROOT}" \
      --dependency="afterany:${ARRAY_JOB_ID}" \
      --export=ALL,RUN_ID="${RUN_ID}",ATLAS_REPO_ROOT="${REPO_ROOT}" \
      "${SCRIPT_DIR}/run_spr_finalize.sh"
  )"
  FINAL_JOB_ID="${FINAL_JOB_RAW%%;*}"

  echo "Submitted array job: ${ARRAY_JOB_ID}"
  echo "Submitted finalize job: ${FINAL_JOB_ID} (afterany:${ARRAY_JOB_ID})"
  echo "RUN_ID=${RUN_ID}"
  exit 0
fi

eval "$(micromamba shell hook --shell bash)"
cd "${ATLAS_REPO_ROOT:-${REPO_ROOT}}"
micromamba activate docking-env
mkdir -p logs/slurm

RUN_ID="${RUN_ID:-${RUN_TAG}_${SLURM_ARRAY_JOB_ID:-local}}"

relocated_args=(
  --all-dirs "${ATLAS_ALL_DIRS:-${REPO_ROOT}}"
  --main-args "${ATLAS_MAIN_ARGS:---resume --run-id ${RUN_ID}}"
)
if [[ -n "${ATLAS_PREPPED_ROOT:-}" ]]; then
  relocated_args+=(--prepped-root "${ATLAS_PREPPED_ROOT}" --skip-extract)
fi

python tools/run_relocated_mode.py "${relocated_args[@]}"
