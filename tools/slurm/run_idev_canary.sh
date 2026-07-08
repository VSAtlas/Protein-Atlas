#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: tools/slurm/run_idev_canary.sh bench2|spd90-dry|spd10-dry-fast|resume

Runs one Atlas canary inside an existing idev/Slurm allocation. The launcher uses
whole allocated nodes by default: 3 nodes, 3 tasks, 1 task per node, no explicit
cpus-per-task cap. Override with ATLAS_IDEV_NODES, ATLAS_IDEV_TASKS, and
ATLAS_IDEV_TASKS_PER_NODE.
EOF
}

if [[ $# -lt 1 ]]; then
  usage
  exit 2
fi

MODE="$1"
shift || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RUN_ROOT="${ATLAS_REPO_ROOT:-${REPO_ROOT}}"
ALL_DIRS_ROOT="${ATLAS_ALL_DIRS:-${RUN_ROOT}}"
PREPPED_ROOT="${ATLAS_PREPPED_ROOT:-}"
NODES="${ATLAS_IDEV_NODES:-3}"
TASKS="${ATLAS_IDEV_TASKS:-3}"
TASKS_PER_NODE="${ATLAS_IDEV_TASKS_PER_NODE:-1}"
SPD_PDB_FILE="${SPD_PDB_FILE:-}"
SPD_LIBRARY="${SPD_LIBRARY:-dry_bench_64}"
SLURM_LOG_ROOT="${ATLAS_SLURM_LOG_ROOT:-${ALL_DIRS_ROOT}/outputs/logs/slurm}"

detect_cpus_per_task() {
  local raw
  for raw in "${SLURM_JOB_CPUS_PER_NODE:-}" "${SLURM_CPUS_ON_NODE:-}"; do
    if [[ "${raw}" =~ ([0-9]+) ]]; then
      echo "${BASH_REMATCH[1]}"
      return 0
    fi
  done
  return 1
}

CPUS_PER_TASK="${ATLAS_IDEV_CPUS_PER_TASK:-$(detect_cpus_per_task || true)}"
SRUN_CPU_ARGS=()
if [[ -n "${CPUS_PER_TASK}" ]]; then
  SRUN_CPU_ARGS+=(--cpus-per-task="${CPUS_PER_TASK}")
fi
USE_NODE_CPUS_DEFAULT=1
if [[ -n "${ATLAS_IDEV_CPUS_PER_TASK:-}" ]]; then
  USE_NODE_CPUS_DEFAULT=0
fi

cd "${RUN_ROOT}"
mkdir -p "${SLURM_LOG_ROOT}"

timestamp="$(date +%Y%m%d_%H%M%S)"
case "${MODE}" in
  bench2)
    RUN_ID="${RUN_ID:-bench2_idev${NODES}_${timestamp}}"
    MAIN_ARGS="-bench2 -fast --run-id ${RUN_ID}"
    ;;
  spd90-dry)
    RUN_ID="${RUN_ID:-spd90_dry64_idev${NODES}_${timestamp}}"
    SPD_PDB_FILE="${SPD_PDB_FILE:-analysis/gene_list/spd_targets_dry90_pdbs.txt}"
    if [[ ! -f "${SPD_PDB_FILE}" ]]; then
      echo "[atlas-idev-canary] missing SPD_PDB_FILE=${SPD_PDB_FILE}" >&2
      exit 2
    fi
    SPD_PDBS="$(tr '\n' ' ' < "${SPD_PDB_FILE}" | sed 's/[[:space:]]*$//')"
    MAIN_ARGS="${SPD_PDBS} --dud-library ${SPD_LIBRARY} -test-fda --run-id ${RUN_ID}"
    ;;
  spd10-dry-fast)
    RUN_ID="${RUN_ID:-spd10_dry64_fast_idev${NODES}_${timestamp}}"
    SPD_PDB_FILE="${SPD_PDB_FILE:-analysis/gene_list/spd_targets_dry10_pdbs.txt}"
    if [[ ! -f "${SPD_PDB_FILE}" ]]; then
      echo "[atlas-idev-canary] missing SPD_PDB_FILE=${SPD_PDB_FILE}" >&2
      exit 2
    fi
    SPD_PDBS="$(tr '\n' ' ' < "${SPD_PDB_FILE}" | sed 's/[[:space:]]*$//')"
    MAIN_ARGS="${SPD_PDBS} --dud-library ${SPD_LIBRARY} -test-fda -fast --run-id ${RUN_ID}"
    ;;
  resume)
    if [[ -z "${RUN_ID:-}" ]]; then
      echo "[atlas-idev-canary] RUN_ID is required for resume mode" >&2
      exit 2
    fi
    MAIN_ARGS="--resume --run-id ${RUN_ID}"
    ;;
  -h|--help|help)
    usage
    exit 0
    ;;
  *)
    usage
    exit 2
    ;;
esac

echo "[atlas-idev-canary] mode=${MODE} run_id=${RUN_ID} job=${SLURM_JOB_ID:-none} nodes=${NODES} tasks=${TASKS} tasks_per_node=${TASKS_PER_NODE} cpus_per_task=${CPUS_PER_TASK:-auto} nodelist=${SLURM_NODELIST:-none} all_dirs=${ALL_DIRS_ROOT}"
if [[ -n "${PREPPED_ROOT}" ]]; then
  echo "[atlas-idev-canary] prepped_root=${PREPPED_ROOT}"
fi
echo "[atlas-idev-canary] start=$(date --iso-8601=seconds)"
echo "[atlas-idev-canary] out=${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.out err=${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.err"

run_srun() {
  export ATLAS_ALL_DIRS="${ALL_DIRS_ROOT}"
  export ATLAS_PREPPED_ROOT="${PREPPED_ROOT}"
  export ATLAS_USE_NODE_CPUS="${ATLAS_USE_NODE_CPUS:-${USE_NODE_CPUS_DEFAULT}}"
  export ATLAS_MAIN_ARGS="${MAIN_ARGS}"

  srun \
    --exclusive \
    --nodes="${NODES}" \
    --ntasks="${TASKS}" \
    --ntasks-per-node="${TASKS_PER_NODE}" \
    "${SRUN_CPU_ARGS[@]}" \
    bash -lc '
      set -euo pipefail
      PYTHON_EXE="$(command -v python || true)"
      if [[ -n "${PYTHON_EXE}" ]]; then
        PYTHON_ENV_PREFIX="$(cd "$(dirname "${PYTHON_EXE}")/.." && pwd)"
        if [[ -d "${PYTHON_ENV_PREFIX}/lib" ]]; then
          export LD_LIBRARY_PATH="${PYTHON_ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
        fi
      fi
      unset LD_PRELOAD
      unset LD_AUDIT
      export ATLAS_DISTRIBUTED_MODE=slurm_array
      export SLURM_ARRAY_TASK_ID="${SLURM_PROCID:-0}"
      export SLURM_ARRAY_TASK_COUNT="${SLURM_NTASKS:-1}"
      export SLURM_ARRAY_TASK_MIN=0
      export SLURM_ARRAY_TASK_MAX=$((SLURM_ARRAY_TASK_COUNT - 1))
      echo "[atlas-idev-canary-task] host=$(hostname) task=${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} cpus_on_node=${SLURM_CPUS_ON_NODE:-} job_cpus_per_node=${SLURM_JOB_CPUS_PER_NODE:-}"
      PREPPED_ARGS=()
      if [[ -n "${ATLAS_PREPPED_ROOT:-}" ]]; then
        PREPPED_ARGS+=(--prepped-root "${ATLAS_PREPPED_ROOT}")
      fi
      python tools/run_relocated_mode.py \
        --skip-extract \
        --all-dirs "${ATLAS_ALL_DIRS}" \
        "${PREPPED_ARGS[@]}" \
        --main-args "${ATLAS_MAIN_ARGS}"
    '
}

RC=0
if [[ "${MODE}" == "spd90-dry" ]]; then
  export TEST_MODE_ENABLE="${TEST_MODE_ENABLE:-dud}"
  export TEST_FDA_LIBRARY_SUBDIR="${TEST_FDA_LIBRARY_SUBDIR:-${SPD_LIBRARY}}"
  export USE_SCORCH=true
  export USE_GNINA=false
  export USE_LEDOCK=false
  export USE_DOCK6=false
  export SCORCH_TOP_FRACTION="${SCORCH_TOP_FRACTION:-1.0}"
  export ATLAS_EXECUTION_MODE="${ATLAS_EXECUTION_MODE:-distributed_combo_hybrid}"
  export ATLAS_SCORCH_PARALLEL_PROFILE="${ATLAS_SCORCH_PARALLEL_PROFILE:-conservative}"
  export ATLAS_SCORCH_QUEUE_MAX_WORKERS="${ATLAS_SCORCH_QUEUE_MAX_WORKERS:-1}"
  run_srun \
    >"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.out" \
    2>"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.err" || RC=$?
elif [[ "${MODE}" == "spd10-dry-fast" ]]; then
  export TEST_MODE_ENABLE="${TEST_MODE_ENABLE:-dud}"
  export TEST_FDA_LIBRARY_SUBDIR="${TEST_FDA_LIBRARY_SUBDIR:-${SPD_LIBRARY}}"
  export USE_SCORCH="${USE_SCORCH:-false}"
  export USE_GNINA=false
  export USE_LEDOCK=false
  export USE_DOCK6=false
  export SCORCH_TOP_FRACTION="${SCORCH_TOP_FRACTION:-0.10}"
  export ATLAS_EXECUTION_MODE="${ATLAS_EXECUTION_MODE:-distributed_combo_hybrid}"
  run_srun \
    >"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.out" \
    2>"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.err" || RC=$?
elif [[ "${MODE}" == "resume" ]]; then
  export USE_SCORCH="${USE_SCORCH:-true}"
  export SCORCH_TOP_FRACTION="${SCORCH_TOP_FRACTION:-0.10}"
  export ATLAS_EXECUTION_MODE="${ATLAS_EXECUTION_MODE:-distributed_combo_hybrid}"
  run_srun \
    >"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.out" \
    2>"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.err" || RC=$?
else
  run_srun \
    >"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.out" \
    2>"${SLURM_LOG_ROOT}/${RUN_ID}_idev_srun.err" || RC=$?
fi

if [[ "${ATLAS_IDEV_SKIP_FINALIZE:-0}" == "1" ]]; then
  echo "[atlas-idev-canary] finish=$(date --iso-8601=seconds) rc=${RC} finalize=skipped"
  exit "${RC}"
fi

FINALIZE_RC=0
export PREPPED_LIGANDS_DIR="${ALL_DIRS_ROOT}/prepped_ligands"
if [[ -n "${PREPPED_ROOT}" ]]; then
  export PREPPED_LIGANDS_DIR="${PREPPED_ROOT}"
fi
export EXTRACTED_LIGANDS_DIR="${ALL_DIRS_ROOT}/extracted_ligands"
export DOCKED_DIR="${ALL_DIRS_ROOT}/outputs/docked"
export POST_DOCKED_DIR="${ALL_DIRS_ROOT}/outputs/post_docked"
export OUTPUT_DIR="${ALL_DIRS_ROOT}/outputs/processed_pdbs"
export CONFIGS_DIR="${ALL_DIRS_ROOT}/outputs/configs"
export MANIFESTS_DIR="${ALL_DIRS_ROOT}/outputs/manifests"
export LOGS_DIR="${ALL_DIRS_ROOT}/outputs/logs"
export DATA_DIR="${ALL_DIRS_ROOT}/outputs/data"
python tools/finalize_distributed_run.py \
  --run-id "${RUN_ID}" \
  >"${SLURM_LOG_ROOT}/${RUN_ID}_finalize.out" \
  2>"${SLURM_LOG_ROOT}/${RUN_ID}_finalize.err" || FINALIZE_RC=$?

MANIFEST_PATH="${MANIFESTS_DIR}/${RUN_ID}/run_manifest.yaml"
MANIFEST_STATUS="missing"
if [[ -s "${MANIFEST_PATH}" ]]; then
  MANIFEST_STATUS="$(python - "${MANIFEST_PATH}" <<'PY'
from pathlib import Path
import sys
try:
    import yaml
    payload = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8")) or {}
    print(str(payload.get("status") or "unknown").strip().lower())
except Exception:
    print("unreadable")
PY
)"
fi
if [[ "${MANIFEST_STATUS}" != "completed" && "${MANIFEST_STATUS}" != "success" ]]; then
  echo "[atlas-idev-canary] manifest_status=${MANIFEST_STATUS} path=${MANIFEST_PATH}" >&2
  if (( FINALIZE_RC == 0 )); then
    FINALIZE_RC=2
  fi
fi

echo "[atlas-idev-canary] finish=$(date --iso-8601=seconds) rc=${RC} finalize_rc=${FINALIZE_RC}"
if (( RC != 0 )); then
  exit "${RC}"
fi
exit "${FINALIZE_RC}"
