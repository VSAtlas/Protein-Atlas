#!/usr/bin/env bash
# Launch simplify-round headless agents (optional nested subagents per stream).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_ROOT="${ATLAS_AGENT_LOG_ROOT:-/stor/home/mpg2352/atlas_parallel_agents}"
RUN_ID="${ATLAS_AGENT_RUN_ID:-simplify_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${LOG_ROOT}/${RUN_ID}"
MODEL="${ATLAS_AGENT_MODEL:-composer-2.5-fast}"
MAX_SUBAGENTS="${ATLAS_SIMPLIFY_MAX_SUBAGENTS:-3}"

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

COMMON_RULES='Repo: protein_automation. Never edit tests/, ci/, input_pdbs/, docked/. No behavior changes. Split files over 100KB or 800+ lines. Canonical imports only (cli, path_router, config, post_docking). Subagents: spawn at most '"${MAX_SUBAGENTS}"' via nested agent --worktree <stream>-subN. Merge sub-worktrees into your stream before AGENT_DONE. Append docs/agent_handoffs.md. Do not git commit unless all verify steps pass.'

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
    "${AGENT_COMMON[@]}" --worktree "${stream_id}" "${prompt}" 2>&1 | tee "${log}"
    echo "AGENT_EXIT:${stream_id}:$?" >>"${log}"
  ) &
  echo $! >"${pidfile}"
  sleep "${ATLAS_AGENT_LAUNCH_STAGGER_SEC:-5}"
}

subagent_block() {
  local stream="$1"
  cat <<EOF

You MAY spawn up to ${MAX_SUBAGENTS} subagents. For each sub-task:
  agent --print --force --trust --yolo --model ${MODEL} \\
    --workspace ${REPO_ROOT} --worktree ${stream}-sub<N> \\
    "SUBAGENT for ${stream}: <task>. ${COMMON_RULES} Print AGENT_DONE:${stream}-sub<N> when done."
Merge subagent branches into worktree ${stream} before finishing.
Read docs/parallel_agent_plan_simplify_rounds.md subagent table for ${stream}.
EOF
}

prompt_simp_01() {
  cat <<EOF
You are stream simp-01-report. ${COMMON_RULES}
Task: Split analysis/reporting/run_report_core.py into analysis/reporting/run_report/ package. Keep public API stable via re-exports. Target: no file >800 lines.
$(subagent_block simp-01-report)
Verify: pytest chemdb/tests/test_run_report.py chemdb/tests/test_run_report_highlights_multi.py -q && ruff check analysis/reporting/run_report_core.py analysis/reporting/run_report/
On success print AGENT_DONE:simp-01-report; else AGENT_FAILED:simp-01-report.
EOF
}

prompt_simp_02() {
  cat <<EOF
You are stream simp-02-heatmap. ${COMMON_RULES}
Task: Split heatmap_html.py and heatmap_html_clustergrammer_template.py — templates vs payload/build logic.
$(subagent_block simp-02-heatmap)
Verify: pytest chemdb/tests/test_heatmap_html_clustergrammer.py chemdb/tests/test_heatmap_outputs.py chemdb/tests/test_report_html_portable.py -q
On success print AGENT_DONE:simp-02-heatmap; else AGENT_FAILED:simp-02-heatmap.
EOF
}

prompt_simp_03() {
  cat <<EOF
You are stream simp-03-dudeval. ${COMMON_RULES}
Task: Split analysis/dud_eval_core/orchestrate.py into orchestrate.py (thin) + orchestrate_*.py helpers (discovery, metrics, cli).
$(subagent_block simp-03-dudeval)
Verify: pytest chemdb/tests/test_dud_eval_class_aggregates.py -q && ruff check analysis/dud_eval_core/
On success print AGENT_DONE:simp-03-dudeval; else AGENT_FAILED:simp-03-dudeval.
EOF
}

prompt_simp_04() {
  cat <<EOF
You are stream simp-04-qol. ${COMMON_RULES}
Task: Split src/cli/qol_cli.py into src/cli/qol/ package; keep cli.qol_cli.dispatch and atlas entrypoints working.
$(subagent_block simp-04-qol)
Verify: pytest chemdb/tests/test_cli_init_new_run.py chemdb/tests/test_cli_help.py chemdb/tests/test_cli_doctor.py -q
On success print AGENT_DONE:simp-04-qol; else AGENT_FAILED:simp-04-qol.
EOF
}

prompt_simp_05() {
  cat <<EOF
You are stream simp-05-manifest. ${COMMON_RULES}
Task: Extract helpers from run_manifest_runtime.py and run_manifest_support.py; reduce main file below 1200 lines.
$(subagent_block simp-05-manifest)
Verify: pytest chemdb/tests/test_run_manifest.py -q
On success print AGENT_DONE:simp-05-manifest; else AGENT_FAILED:simp-05-manifest.
EOF
}

prompt_simp_06() {
  cat <<EOF
You are stream simp-06-postrun. ${COMMON_RULES}
Task: Split postrun_hooks_runtime.py by hook phase (scorch, report, retention, integrity).
$(subagent_block simp-06-postrun)
Verify: pytest chemdb/tests/test_run_lifecycle_scorch_barrier.py chemdb/tests/test_scorch_per_protein_rescore.py -q
On success print AGENT_DONE:simp-06-postrun; else AGENT_FAILED:simp-06-postrun.
EOF
}

prompt_simp_07() {
  cat <<EOF
You are stream simp-07-dist-plan. ${COMMON_RULES}
Task: Simplify distributed_chunk_planner.py, planner_chunking.py, planner_manifest.py via extraction; no planner behavior change.
$(subagent_block simp-07-dist-plan)
Verify: pytest chemdb/tests/test_distributed_chunk_planning.py -q
On success print AGENT_DONE:simp-07-dist-plan; else AGENT_FAILED:simp-07-dist-plan.
EOF
}

prompt_simp_08() {
  cat <<EOF
You are stream simp-08-dist-run. ${COMMON_RULES}
Task: Split distributed_runtime_claim_loop.py and slim distributed_chunk_runtime.py claim/completion paths.
$(subagent_block simp-08-dist-run)
Verify: pytest chemdb/tests/test_distributed_chunk_regressions.py chemdb/tests/test_hybrid_chunk_scheduler.py -q
On success print AGENT_DONE:simp-08-dist-run; else AGENT_FAILED:simp-08-dist-run.
EOF
}

prompt_simp_09() {
  cat <<EOF
You are stream simp-09-status. ${COMMON_RULES}
Task: Extract status_dashboard.py panels and doctor.py display helpers; reduce complexity for xenon.
$(subagent_block simp-09-status)
Verify: pytest chemdb/tests/test_status_dashboard.py chemdb/tests/test_cli_doctor.py -q
On success print AGENT_DONE:simp-09-status; else AGENT_FAILED:simp-09-status.
EOF
}

prompt_simp_10() {
  cat <<EOF
You are stream simp-10-schema. ${COMMON_RULES}
Task: Split analysis/reporting/master_schema_export.py into analysis/reporting/master_schema/ package; keep public API via re-exports.
$(subagent_block simp-10-schema)
Verify: pytest chemdb/tests/test_master_schema_export.py chemdb/tests/test_master_schema_export_decoy_prefix.py -q
On success print AGENT_DONE:simp-10-schema; else AGENT_FAILED:simp-10-schema.
EOF
}

prompt_simp_11() {
  cat <<EOF
You are stream simp-11-ml. ${COMMON_RULES}
Task: Split analysis/ml/train_classifier.py and slim analysis/ml/audit_suite.py via extraction; no training behavior change.
$(subagent_block simp-11-ml)
Verify: pytest chemdb/tests/test_ml_eval_metrics_wrapper.py chemdb/tests/test_ml_metrics_extra_ece_brier.py -q
On success print AGENT_DONE:simp-11-ml; else AGENT_FAILED:simp-11-ml.
EOF
}

prompt_simp_12() {
  cat <<EOF
You are stream simp-12-scorch. ${COMMON_RULES}
Task: Further split src/post_docking/rescoring/rescoring_scorch.py and scorch_orchestration.py; preserve SCORCH hook/test compatibility.
$(subagent_block simp-12-scorch)
Verify: pytest chemdb/tests/test_scorch_parallelism_allocator.py chemdb/tests/test_scorch_materialize_inputs.py chemdb/tests/test_scorch_no_combos_logging.py chemdb/tests/test_scorch_rescore_decoy_prefix.py -q
On success print AGENT_DONE:simp-12-scorch; else AGENT_FAILED:simp-12-scorch.
EOF
}

WAVE1_STREAMS=(
  simp-01-report
  simp-02-heatmap
  simp-03-dudeval
  simp-04-qol
  simp-05-manifest
  simp-06-postrun
  simp-07-dist-plan
  simp-08-dist-run
  simp-09-status
)

launch_wave1() {
  mkdir -p "${RUN_DIR}"
  echo "${RUN_ID}" >"${RUN_DIR}/RUN_ID"
  echo "[orchestrator] wave1 RUN_ID=${RUN_ID} MODEL=${MODEL} MAX_SUBAGENTS=${MAX_SUBAGENTS}"

  run_stream simp-01-report "$(prompt_simp_01)"
  run_stream simp-02-heatmap "$(prompt_simp_02)"
  run_stream simp-03-dudeval "$(prompt_simp_03)"
  run_stream simp-04-qol "$(prompt_simp_04)"
  run_stream simp-05-manifest "$(prompt_simp_05)"
  run_stream simp-06-postrun "$(prompt_simp_06)"
  run_stream simp-07-dist-plan "$(prompt_simp_07)"
  run_stream simp-08-dist-run "$(prompt_simp_08)"
  run_stream simp-09-status "$(prompt_simp_09)"

  echo "[orchestrator] launched ${#WAVE1_STREAMS[@]} simplify streams -> ${RUN_DIR}"
}

launch_wave2() {
  mkdir -p "${RUN_DIR}"
  echo "${RUN_ID}" >"${RUN_DIR}/RUN_ID"
  echo "[orchestrator] wave2 RUN_ID=${RUN_ID} MODEL=${MODEL} MAX_SUBAGENTS=${MAX_SUBAGENTS}"

  run_stream simp-10-schema "$(prompt_simp_10)"
  run_stream simp-11-ml "$(prompt_simp_11)"
  run_stream simp-12-scorch "$(prompt_simp_12)"

  echo "[orchestrator] launched 3 wave2 streams -> ${RUN_DIR}"
}

launch_one() {
  mkdir -p "${RUN_DIR}"
  local name="$1"
  case "${name}" in
    simp-01-report) run_stream "${name}" "$(prompt_simp_01)" ;;
    simp-02-heatmap) run_stream "${name}" "$(prompt_simp_02)" ;;
    simp-03-dudeval) run_stream "${name}" "$(prompt_simp_03)" ;;
    simp-04-qol) run_stream "${name}" "$(prompt_simp_04)" ;;
    simp-05-manifest) run_stream "${name}" "$(prompt_simp_05)" ;;
    simp-06-postrun) run_stream "${name}" "$(prompt_simp_06)" ;;
    simp-07-dist-plan) run_stream "${name}" "$(prompt_simp_07)" ;;
    simp-08-dist-run) run_stream "${name}" "$(prompt_simp_08)" ;;
    simp-09-status) run_stream "${name}" "$(prompt_simp_09)" ;;
    simp-10-schema) run_stream "${name}" "$(prompt_simp_10)" ;;
    simp-11-ml) run_stream "${name}" "$(prompt_simp_11)" ;;
    simp-12-scorch) run_stream "${name}" "$(prompt_simp_12)" ;;
    *) echo "Unknown stream: ${name}" >&2; exit 2 ;;
  esac
}

status() {
  local dir="${1:-}"
  if [[ -z "${dir}" ]]; then
    dir="$(ls -1dt "${LOG_ROOT}"/*/RUN_ID 2>/dev/null | head -1 | xargs dirname 2>/dev/null || true)"
  fi
  [[ -d "${dir}" ]] || { echo "No run dir under ${LOG_ROOT}" >&2; exit 1; }
  echo "=== ${dir} ==="
  for f in "${dir}"/*.pid; do
    [[ -f "${f}" ]] || continue
    base="$(basename "${f}" .pid)"
    pid="$(cat "${f}")"
    state="stopped"
    kill -0 "${pid}" 2>/dev/null && state="running"
    marker=""
    grep -q "AGENT_DONE:${base}" "${dir}/${base}.log" 2>/dev/null && marker=" DONE"
    grep -q "AGENT_FAILED:${base}" "${dir}/${base}.log" 2>/dev/null && marker=" FAILED"
    echo "  ${base}: pid=${pid} ${state}${marker}"
  done
}

usage() {
  cat <<EOF
Usage:
  ATLAS_AGENT_RUN_ID=<id> $0 launch [wave1|wave2]  # wave1: 9 streams, wave2: 3 streams
  ATLAS_AGENT_RUN_ID=<id> $0 launch <stream>   # one stream
  $0 status [LOG_DIR]

Docs: docs/parallel_agent_plan_simplify_rounds.md
Merge: tools/merge_parallel_simplify_worktrees.sh wave1|wave2
EOF
}

main() {
  chmod +x "${BASH_SOURCE[0]}" 2>/dev/null || true
  case "${1:-}" in
    launch)
      if [[ "${2:-wave1}" == "wave1" || -z "${2:-}" ]]; then
        launch_wave1
      elif [[ "${2:-}" == "wave2" ]]; then
        launch_wave2
      else
        launch_one "${2}"
      fi
      ;;
    status) status "${2:-}" ;;
    *) usage; exit 2 ;;
  esac
}

main "$@"
