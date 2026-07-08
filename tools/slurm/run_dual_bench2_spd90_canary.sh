#!/usr/bin/env bash
#SBATCH -J atlas_dual_canary
#SBATCH -N 6
#SBATCH --exclusive
#SBATCH -t 00:30:00
#SBATCH -p skx-dev
#SBATCH -o logs/slurm/dual_canary_%j.out
#SBATCH -e logs/slurm/dual_canary_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${ATLAS_REPO_ROOT:-}" ]]; then
  REPO_ROOT="$(cd "${ATLAS_REPO_ROOT}" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/src/cli" ]]; then
  REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
  REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
RUN_ROOT="${ATLAS_REPO_ROOT:-${REPO_ROOT}}"

eval "$(micromamba shell hook --shell bash)"
cd "${RUN_ROOT}"
micromamba activate docking-env
cd "${RUN_ROOT}"
mkdir -p "${RUN_ROOT}/logs/slurm"

BENCH_RUN_ID="${BENCH_RUN_ID:-bench2_canary_3n_30m_${SLURM_JOB_ID:-local}}"
SPD_RUN_ID="${SPD_RUN_ID:-spd90_dry_bench64_3n_30m_${SLURM_JOB_ID:-local}}"
SPD_PDB_FILE="${SPD_PDB_FILE:-analysis/gene_list/spd_targets_dry90_pdbs.txt}"
SPD_LIBRARY="${SPD_LIBRARY:-dry_bench_64}"

mapfile -t ALLOC_NODES < <(scontrol show hostnames "${SLURM_NODELIST}")
if (( ${#ALLOC_NODES[@]} < 6 )); then
  echo "[dual-canary] need at least 6 allocated nodes, got ${#ALLOC_NODES[@]}" >&2
  exit 2
fi

join_by_comma() {
  local IFS=,
  echo "$*"
}

BENCH_NODES="$(join_by_comma "${ALLOC_NODES[@]:0:3}")"
SPD_NODES="$(join_by_comma "${ALLOC_NODES[@]:3:3}")"
SPD_PDBS="$(tr '\n' ' ' < "${SPD_PDB_FILE}" | sed 's/[[:space:]]*$//')"

run_atlas_srun() {
  local node_list="$1"
  local run_id="$2"
  local main_args="$3"
  local log_prefix="$4"
  shift 4

  (
    export ATLAS_ALL_DIRS="${RUN_ROOT}"
    export ATLAS_USE_NODE_CPUS=1
    export ATLAS_MAIN_ARGS="${main_args}"
    "$@" srun \
      --exclusive \
      --nodes=3 \
      --ntasks=3 \
      --ntasks-per-node=1 \
      --nodelist="${node_list}" \
      bash -lc '
        set -euo pipefail
        export ATLAS_DISTRIBUTED_MODE=slurm_array
        export SLURM_ARRAY_TASK_ID="${SLURM_PROCID:-0}"
        export SLURM_ARRAY_TASK_COUNT="${SLURM_NTASKS:-3}"
        export SLURM_ARRAY_TASK_MIN=0
        export SLURM_ARRAY_TASK_MAX=$((SLURM_ARRAY_TASK_COUNT - 1))
        export RUN_ID="'"${run_id}"'"
        python tools/run_relocated_mode.py \
          --skip-extract \
          --all-dirs "${ATLAS_ALL_DIRS}" \
          --main-args "${ATLAS_MAIN_ARGS}"
      '
  ) >"logs/slurm/${log_prefix}.out" 2>"logs/slurm/${log_prefix}.err"
}

echo "[dual-canary] job=${SLURM_JOB_ID:-local} bench_run=${BENCH_RUN_ID} bench_nodes=${BENCH_NODES}"
echo "[dual-canary] job=${SLURM_JOB_ID:-local} spd_run=${SPD_RUN_ID} spd_nodes=${SPD_NODES} pdb_file=${SPD_PDB_FILE} library=${SPD_LIBRARY}"

run_atlas_srun \
  "${BENCH_NODES}" \
  "${BENCH_RUN_ID}" \
  "-bench2 -fast --run-id ${BENCH_RUN_ID}" \
  "${BENCH_RUN_ID}_srun" &
BENCH_PID=$!

run_atlas_srun \
  "${SPD_NODES}" \
  "${SPD_RUN_ID}" \
  "${SPD_PDBS} --dud-library ${SPD_LIBRARY} --run-id ${SPD_RUN_ID}" \
  "${SPD_RUN_ID}_srun" \
  env USE_SCORCH=true USE_GNINA=false USE_LEDOCK=false USE_DOCK6=false SCORCH_TOP_FRACTION=1.0 ATLAS_SCORCH_PARALLEL_PROFILE=aggressive &
SPD_PID=$!

BENCH_RC=0
SPD_RC=0
wait "${BENCH_PID}" || BENCH_RC=$?
wait "${SPD_PID}" || SPD_RC=$?

echo "[dual-canary] bench_rc=${BENCH_RC} spd_rc=${SPD_RC}"

python tools/finalize_distributed_run.py --run-id "${BENCH_RUN_ID}" >"logs/slurm/${BENCH_RUN_ID}_finalize.out" 2>"logs/slurm/${BENCH_RUN_ID}_finalize.err" || true
python tools/finalize_distributed_run.py --run-id "${SPD_RUN_ID}" >"logs/slurm/${SPD_RUN_ID}_finalize.out" 2>"logs/slurm/${SPD_RUN_ID}_finalize.err" || true

if (( BENCH_RC != 0 )); then
  exit "${BENCH_RC}"
fi
exit "${SPD_RC}"
