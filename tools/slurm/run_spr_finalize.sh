#!/usr/bin/env bash
#SBATCH -J atlas_spr_finalize
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

cd "${ATLAS_REPO_ROOT:-${REPO_ROOT}}"
mkdir -p logs/slurm

RUN_TAG="${RUN_TAG:-spr}"
RUN_ID="${RUN_ID:-${RUN_TAG}_${SLURM_JOB_ID:-local}}"

python tools/finalize_distributed_run.py --run-id "${RUN_ID}"
