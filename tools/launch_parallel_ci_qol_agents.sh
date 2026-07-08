#!/usr/bin/env bash
# Launch headless Cursor agents in parallel worktrees for green CI + QoL CLI.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_ROOT="${ATLAS_AGENT_LOG_ROOT:-/stor/home/mpg2352/atlas_parallel_agents}"
RUN_ID="${ATLAS_AGENT_RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${LOG_ROOT}/${RUN_ID}"
MODEL="${ATLAS_AGENT_MODEL:-composer-2.5-fast}"

AGENT_COMMON=(
  agent
  --print
  --force
  --trust
  --yolo
  --model
  "${MODEL}"
  --workspace
  "${REPO_ROOT}"
)

COMMON_RULES='Repo: protein_automation. Never edit tests/, ci/, input_pdbs/, docked/. New tests only under chemdb/tests/. Use src.path_router. Append one line to docs/agent_handoffs.md when done. Do not git commit.'

run_stream() {
  local stream_id="$1"
  local prompt="$2"
  local log="${RUN_DIR}/${stream_id}.log"
  local pidfile="${RUN_DIR}/${stream_id}.pid"

  if [[ -f "${pidfile}" ]] && kill -0 "$(cat "${pidfile}")" 2>/dev/null; then
    echo "[skip] ${stream_id} already running (pid $(cat "${pidfile}"))"
    return 0
  fi

  echo "[launch] ${stream_id} -> ${log}"
  (
    cd "${REPO_ROOT}"
    "${AGENT_COMMON[@]}" \
      --worktree "${stream_id}" \
      "${prompt}" \
      2>&1 | tee "${log}"
    echo "AGENT_EXIT:${stream_id}:$?" >>"${log}"
  ) &
  echo $! >"${pidfile}"
  # Avoid concurrent writes to ~/.cursor/cli-config.json when many agents start at once.
  sleep "${ATLAS_AGENT_LAUNCH_STAGGER_SEC:-4}"
}

prompt_ci_ruff() {
  cat <<EOF
You are stream ci-ruff. ${COMMON_RULES}
Task: Make 'ruff check .' pass. Run 'ruff check . --fix' first, then fix remaining issues manually.
Verify: ruff check .
On success print AGENT_DONE:ci-ruff; else AGENT_FAILED:ci-ruff with the failing command output summary.
EOF
}

prompt_ci_mypy_gate() {
  cat <<EOF
You are stream ci-mypy-gate. ${COMMON_RULES}
Task: Fix mypy errors for the scoped gate targets in tools/quality_gate.sh (MYPY_TARGETS), especially tools/finalize_distributed_run.py:334 assignment error.
Verify: cd ${REPO_ROOT} && mypy --follow-imports skip src/path_router src/config src/cli/distributed_chunk_planner.py src/cli/distributed_chunk_runtime.py src/cli/distributed_mode_policy.py src/cli/hybrid_chunk_scheduler.py src/cli/main_script_entrypoint.py tools/simulate_slurm_array.py tools/benchmark_bench2_compare.py tools/finalize_distributed_run.py
On success print AGENT_DONE:ci-mypy-gate; else AGENT_FAILED:ci-mypy-gate.
EOF
}

prompt_ci_mypy_analysis() {
  cat <<EOF
You are stream ci-mypy-analysis. ${COMMON_RULES}
Task: Reduce mypy errors under analysis/ that block 'mypy .'. Prefer minimal typing fixes; do not refactor report HTML templates unless required.
Verify: mypy analysis --follow-imports skip (fix until clean or document remaining count < 5 blockers with file:line list).
On success print AGENT_DONE:ci-mypy-analysis; else AGENT_FAILED:ci-mypy-analysis.
EOF
}

prompt_ci_vulture() {
  cat <<EOF
You are stream ci-vulture. ${COMMON_RULES}
Task: Address vulture findings in tools/ with --min-confidence 80 first (delete dead code or fix logic). For 60% confidence public API aliases in tools/pdb_query.py, add a short comment and vulture whitelist entry in pyproject.toml or setup.cfg if the project has one; otherwise remove only if truly unused.
Do not touch src/ except tools/.
Verify: vulture tools --min-confidence 80
On success print AGENT_DONE:ci-vulture; else AGENT_FAILED:ci-vulture.
EOF
}

prompt_ci_xenon() {
  cat <<EOF
You are stream ci-xenon. ${COMMON_RULES}
Task: Make 'python tools/architecture_gate.py' exit 0. Prefer small refactors for src/cli/doctor.py and src/cli/main_pipeline.py if feasible; otherwise add the violating file paths to tools/xenon_baseline.txt (one path per line, sorted) per docs/tools.architecture_gate.py.md policy.
Verify: python tools/architecture_gate.py
On success print AGENT_DONE:ci-xenon; else AGENT_FAILED:ci-xenon.
EOF
}

prompt_qol_test_init() {
  cat <<EOF
You are stream qol-test-init. ${COMMON_RULES}
Task: Add chemdb/tests/test_cli_init_new_run.py with hermetic tests for atlas init, atlas init --write-example-inputs, and atlas new-run --dry-run (subprocess or API hooks; follow test_cli_help.py patterns). Touch src/cli/qol_cli.py only if a test hook is required.
Verify: pytest chemdb/tests/test_cli_init_new_run.py -q
On success print AGENT_DONE:qol-test-init; else AGENT_FAILED:qol-test-init.
EOF
}

prompt_qol_test_status() {
  cat <<EOF
You are stream qol-test-status. ${COMMON_RULES}
Task: Expand chemdb/tests/test_status_dashboard.py for atlas status --explain, --html, and ETA/history cache behavior using fixtures (see existing slurm fixtures). Do not require live Slurm.
Verify: pytest chemdb/tests/test_status_dashboard.py -q
On success print AGENT_DONE:qol-test-status; else AGENT_FAILED:qol-test-status.
EOF
}

prompt_qol_test_doctor() {
  cat <<EOF
You are stream qol-test-doctor. ${COMMON_RULES}
Task: Expand chemdb/tests/test_cli_doctor.py for atlas doctor --pdb/--ligands paths (mock prep outputs). Improve error messages in src/cli/doctor.py only if tests expose gaps.
Verify: pytest chemdb/tests/test_cli_doctor.py -q
On success print AGENT_DONE:qol-test-doctor; else AGENT_FAILED:qol-test-doctor.
EOF
}

prompt_qol_test_help_dev() {
  cat <<EOF
You are stream qol-test-help-dev. ${COMMON_RULES}
Task: Expand chemdb/tests/test_cli_help.py to assert atlas dev verify help text and that tools/quality_gate.sh --full is invoked by the dev verify code path (mock subprocess). No live smoke in unit tests.
Verify: pytest chemdb/tests/test_cli_help.py -q
On success print AGENT_DONE:qol-test-help-dev; else AGENT_FAILED:qol-test-help-dev.
EOF
}

launch_all() {
  mkdir -p "${RUN_DIR}"
  echo "${RUN_ID}" >"${RUN_DIR}/RUN_ID"
  echo "[orchestrator] RUN_ID=${RUN_ID} LOG_DIR=${RUN_DIR} MODEL=${MODEL}"

  run_stream "ci-ruff" "$(prompt_ci_ruff)"
  run_stream "ci-mypy-gate" "$(prompt_ci_mypy_gate)"
  run_stream "ci-mypy-analysis" "$(prompt_ci_mypy_analysis)"
  run_stream "ci-vulture" "$(prompt_ci_vulture)"
  run_stream "ci-xenon" "$(prompt_ci_xenon)"
  run_stream "qol-test-init" "$(prompt_qol_test_init)"
  run_stream "qol-test-status" "$(prompt_qol_test_status)"
  run_stream "qol-test-doctor" "$(prompt_qol_test_doctor)"
  run_stream "qol-test-help-dev" "$(prompt_qol_test_help_dev)"

  echo "[orchestrator] launched 9 streams; logs in ${RUN_DIR}"
}

launch_one() {
  mkdir -p "${RUN_DIR}"
  local name="$1"
  case "${name}" in
    ci-ruff) run_stream "ci-ruff" "$(prompt_ci_ruff)" ;;
    ci-mypy-gate) run_stream "ci-mypy-gate" "$(prompt_ci_mypy_gate)" ;;
    ci-mypy-analysis) run_stream "ci-mypy-analysis" "$(prompt_ci_mypy_analysis)" ;;
    ci-vulture) run_stream "ci-vulture" "$(prompt_ci_vulture)" ;;
    ci-xenon) run_stream "ci-xenon" "$(prompt_ci_xenon)" ;;
    qol-test-init) run_stream "qol-test-init" "$(prompt_qol_test_init)" ;;
    qol-test-status) run_stream "qol-test-status" "$(prompt_qol_test_status)" ;;
    qol-test-doctor) run_stream "qol-test-doctor" "$(prompt_qol_test_doctor)" ;;
    qol-test-help-dev) run_stream "qol-test-help-dev" "$(prompt_qol_test_help_dev)" ;;
    *) echo "Unknown stream: ${name}" >&2; exit 2 ;;
  esac
}

status() {
  local dir="${1:-}"
  if [[ -z "${dir}" ]]; then
    dir="$(ls -1dt "${LOG_ROOT}"/*/RUN_ID 2>/dev/null | head -1 | xargs dirname 2>/dev/null || true)"
  fi
  if [[ -z "${dir}" || ! -d "${dir}" ]]; then
    echo "No run dir found under ${LOG_ROOT}" >&2
    exit 1
  fi
  echo "=== status ${dir} ==="
  for f in "${dir}"/*.pid; do
    [[ -f "${f}" ]] || continue
    local base
    base="$(basename "${f}" .pid)"
    local pid
    pid="$(cat "${f}")"
    local state="stopped"
    kill -0 "${pid}" 2>/dev/null && state="running"
    local marker=""
    grep -q "AGENT_DONE:${base}" "${dir}/${base}.log" 2>/dev/null && marker=" DONE"
    grep -q "AGENT_FAILED:${base}" "${dir}/${base}.log" 2>/dev/null && marker=" FAILED"
    echo "  ${base}: pid=${pid} ${state}${marker}"
  done
  echo "--- markers ---"
  grep -h 'AGENT_DONE:\|AGENT_FAILED:\|AGENT_EXIT:' "${dir}"/*.log 2>/dev/null || true
}

usage() {
  cat <<EOF
Usage:
  ATLAS_AGENT_RUN_ID=<id> $0 launch          # start all 9 streams
  ATLAS_AGENT_RUN_ID=<id> $0 launch <stream> # start one stream
  $0 status [LOG_DIR]

Env:
  ATLAS_AGENT_RUN_ID   run folder name (default: timestamp)
  ATLAS_AGENT_LOG_ROOT base log dir (default: /stor/home/mpg2352/atlas_parallel_agents)
  ATLAS_AGENT_MODEL    agent model (default: composer-2.5-fast)
EOF
}

main() {
  chmod +x "${BASH_SOURCE[0]}" 2>/dev/null || true
  case "${1:-}" in
    launch)
      if [[ -n "${2:-}" ]]; then
        launch_one "${2}"
      else
        launch_all
      fi
      ;;
    status) status "${2:-}" ;;
    *) usage; exit 2 ;;
  esac
}

main "$@"
