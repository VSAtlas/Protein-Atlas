#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  tools/bench2_refactor_verify.sh \
    --baseline-run-id <RUN_ID> \
    --candidate-run-id <RUN_ID> \
    [--reference <CSV_OR_TSV>] \
    [--simulate-distributed] \
    [--tasks <N>] \
    [--cpus-per-task <N>] \
    [--skip-runs]

Description:
  Runs a baseline and candidate bench2-fast benchmark, captures both with
  benchmark_bench2_compare.py, and executes compare.

Defaults:
  tasks=3
  cpus-per-task=4
  max distributed CPU budget (tasks * cpus-per-task) <= 32

Notes:
  - Without --simulate-distributed, runs:
      python main.py -bench2 -fast --run-id <RUN_ID>
  - With --simulate-distributed, runs:
      python tools/simulate_slurm_array.py --tasks ... --cpus-per-task ... \
        --run-id <RUN_ID> --main-args "-bench2 -fast"
  - Use --skip-runs to only run capture + compare for existing run IDs.
EOF
}

BASELINE_RUN_ID=""
CANDIDATE_RUN_ID=""
REFERENCE_PATH=""
SIMULATE_DISTRIBUTED=0
TASKS=3
CPUS_PER_TASK=4
SKIP_RUNS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseline-run-id)
      BASELINE_RUN_ID="${2:-}"
      shift 2
      ;;
    --candidate-run-id)
      CANDIDATE_RUN_ID="${2:-}"
      shift 2
      ;;
    --reference)
      REFERENCE_PATH="${2:-}"
      shift 2
      ;;
    --simulate-distributed)
      SIMULATE_DISTRIBUTED=1
      shift
      ;;
    --tasks)
      TASKS="${2:-}"
      shift 2
      ;;
    --cpus-per-task)
      CPUS_PER_TASK="${2:-}"
      shift 2
      ;;
    --skip-runs)
      SKIP_RUNS=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "$BASELINE_RUN_ID" || -z "$CANDIDATE_RUN_ID" ]]; then
  echo "ERROR: --baseline-run-id and --candidate-run-id are required." >&2
  usage
  exit 2
fi

if ! [[ "$TASKS" =~ ^[0-9]+$ ]] || ! [[ "$CPUS_PER_TASK" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --tasks and --cpus-per-task must be positive integers." >&2
  exit 2
fi
if (( TASKS <= 0 || CPUS_PER_TASK <= 0 )); then
  echo "ERROR: --tasks and --cpus-per-task must be > 0." >&2
  exit 2
fi

if (( SIMULATE_DISTRIBUTED == 1 )); then
  TOTAL_CPU=$(( TASKS * CPUS_PER_TASK ))
  if (( TOTAL_CPU > 32 )); then
    echo "ERROR: distributed simulation CPU budget exceeds 32 (${TOTAL_CPU})." >&2
    exit 2
  fi
fi

run_one() {
  local run_id="$1"
  if (( SIMULATE_DISTRIBUTED == 1 )); then
    echo "[bench2.verify] run_id=${run_id} mode=simulate tasks=${TASKS} cpus_per_task=${CPUS_PER_TASK}"
    python tools/simulate_slurm_array.py \
      --tasks "${TASKS}" \
      --cpus-per-task "${CPUS_PER_TASK}" \
      --run-id "${run_id}" \
      --main-args "-bench2 -fast"
  else
    echo "[bench2.verify] run_id=${run_id} mode=direct"
    python main.py -bench2 -fast --run-id "${run_id}"
  fi
}

capture_one() {
  local run_id="$1"
  if [[ -n "$REFERENCE_PATH" ]]; then
    python tools/benchmark_bench2_compare.py capture --run-id "${run_id}" --reference "${REFERENCE_PATH}"
  else
    python tools/benchmark_bench2_compare.py capture --run-id "${run_id}"
  fi
}

compare_runs() {
  if [[ -n "$REFERENCE_PATH" ]]; then
    python tools/benchmark_bench2_compare.py compare \
      --baseline-run-id "${BASELINE_RUN_ID}" \
      --optimized-run-id "${CANDIDATE_RUN_ID}" \
      --reference "${REFERENCE_PATH}"
  else
    python tools/benchmark_bench2_compare.py compare \
      --baseline-run-id "${BASELINE_RUN_ID}" \
      --optimized-run-id "${CANDIDATE_RUN_ID}"
  fi
}

if (( SKIP_RUNS == 0 )); then
  run_one "${BASELINE_RUN_ID}"
  run_one "${CANDIDATE_RUN_ID}"
else
  echo "[bench2.verify] --skip-runs enabled; using existing run artifacts."
fi

echo "[bench2.verify] capture baseline=${BASELINE_RUN_ID}"
capture_one "${BASELINE_RUN_ID}"
echo "[bench2.verify] capture candidate=${CANDIDATE_RUN_ID}"
capture_one "${CANDIDATE_RUN_ID}"
echo "[bench2.verify] compare baseline=${BASELINE_RUN_ID} candidate=${CANDIDATE_RUN_ID}"
compare_runs
echo "[bench2.verify] complete"
