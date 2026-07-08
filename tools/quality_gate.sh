#!/usr/bin/env bash
set -euo pipefail

RUN_SMOKE=0
APPLY_FIX=0
RUN_FULL=0
EXIT_CODE=0

run_step() {
  "$@" || EXIT_CODE=$?
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --smoke)
      RUN_SMOKE=1
      shift
      ;;
    --fix)
      APPLY_FIX=1
      shift
      ;;
    --full)
      RUN_FULL=1
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 [--fix] [--smoke] [--full]" >&2
      exit 2
      ;;
  esac
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}" || exit 1

RUFF_ARGS=(check)
if [[ "${APPLY_FIX}" -eq 1 ]]; then
  RUFF_ARGS+=(--fix)
fi

RUFF_TARGETS=(
  src/cli
  src/config
  src/path_router
  analysis/reporting
  analysis/dud_eval_core
  src/post_docking/rescoring
  tools/simulate_slurm_array.py
  tools/benchmark_bench2_compare.py
  tools/finalize_distributed_run.py
  main.py
)

MYPY_TARGETS=(
  src/path_router
  src/config
  src/cli/distributed_chunk_planner.py
  src/cli/distributed_chunk_runtime.py
  src/cli/distributed_mode_policy.py
  src/cli/hybrid_chunk_scheduler.py
  src/cli/main_script_entrypoint.py
  tools/simulate_slurm_array.py
  tools/benchmark_bench2_compare.py
  tools/finalize_distributed_run.py
)

# Gate scope policy:
# - Focus on maintained owner modules first.
# - Keep quarantined/legacy scripts excluded until intentionally adopted.
# - Removed stale root targets (postrun_hooks.py, run_manifest.py) after relocation to src/cli/*.
# - Ruff covers stable reporting/evaluation owners before broader type adoption;
#   MyPy remains package-level only where import contracts are stable.

if [[ "${RUN_FULL}" -eq 1 ]]; then
  RUFF_TARGETS=(".")
fi

echo "[quality] ruff ${RUFF_ARGS[*]} ${RUFF_TARGETS[*]}"
run_step ruff "${RUFF_ARGS[@]}" "${RUFF_TARGETS[@]}"

if [[ "${RUN_FULL}" -eq 1 ]]; then
  echo "[quality] vulture . --min-confidence 80"
  run_step vulture . --min-confidence 80 --exclude ".git,.mypy_cache,.pytest_cache,.ruff_cache,__pycache__,chemdb/tests,docked,input_pdbs,logs,post_docked,prepped_ligands,processed_pdbs,tools/autodock4zn,tools/Vina-GPU-2.1"

  echo "[quality] tach check"
  run_step tach check

  echo "[quality] python tools/architecture_gate.py"
  run_step python tools/architecture_gate.py

  echo "[quality] mypy ."
  run_step mypy .
else
  echo "[quality] mypy --follow-imports skip ${MYPY_TARGETS[*]}"
  run_step mypy --follow-imports skip "${MYPY_TARGETS[@]}"
fi

if [[ "${RUN_SMOKE}" -eq 1 ]]; then
  echo "[quality] python main.py --test -fast"
  run_step python main.py --test -fast
fi

echo "[quality] complete"
exit "${EXIT_CODE}"
