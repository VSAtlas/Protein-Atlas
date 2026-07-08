#!/usr/bin/env bash
set -u

RUN_ID="${1:-${RUN_ID:-}}"
if [[ -z "${RUN_ID}" ]]; then
  echo "Usage: $0 <RUN_ID>   (or RUN_ID=<id> $0)"
  exit 0
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}" || exit 0

D="manifests/${RUN_ID}/distributed"
P="post_docked/${RUN_ID}"
M="manifests/${RUN_ID}/run_manifest.yaml"

shopt -s nullglob

count_glob() {
  local -a _matches=($1)
  echo "${#_matches[@]}"
}

count_named_files() {
  local base_dir="$1"
  local filename="$2"
  if command -v rg >/dev/null 2>&1; then
    rg --files "${base_dir}" -g "${filename}" 2>/dev/null | wc -l | tr -d ' '
  else
    find "${base_dir}" -name "${filename}" 2>/dev/null | wc -l | tr -d ' '
  fi
}

if [[ ! -d "${D}" ]]; then
  echo "run=${RUN_ID} status=missing_distributed_dir path=${D} cwd=${PWD}"
  exit 0
fi

DONE=$(count_glob "${D}/chunk_results/*.json")
CLAIMS=$(count_glob "${D}/chunk_claims/*.json")
chunk_plan_files=("${D}"/combo_chunks_*.json)
if [[ ${#chunk_plan_files[@]} -gt 0 ]]; then
  PLANNED=$(grep -rho '"chunk_id"' "${chunk_plan_files[@]}" 2>/dev/null | wc -l | tr -d ' ' || true)
else
  PLANNED="0"
fi
PHASE=$(ls "${D}/phase_markers" 2>/dev/null | paste -sd, - || true)
[[ -z "${PHASE}" ]] && PHASE="-"

SCORCH=$(count_named_files "${P}" 'scorch_scores_all.csv')
CONS=$(count_named_files "${P}" 'consensus_reranked_scorch.csv')
PCT=$(awk -v d="${DONE:-0}" -v p="${PLANNED:-0}" 'BEGIN{if(p>0) printf "%.1f",100*d/p; else printf "0.0"}')

SCHED=$(grep -m1 -E 'total_combos_scheduled:' "${M}" 2>/dev/null | awk -F: '{gsub(/ /,"",$2); print $2}' || true)
COMP=$(grep -m1 -E 'total_combos_completed:' "${M}" 2>/dev/null | awk -F: '{gsub(/ /,"",$2); print $2}' || true)
FAIL=$(grep -m1 -E 'total_combos_failed:' "${M}" 2>/dev/null | awk -F: '{gsub(/ /,"",$2); print $2}' || true)
[[ -z "${SCHED}" ]] && SCHED="?"
[[ -z "${COMP}" ]] && COMP="?"
[[ -z "${FAIL}" ]] && FAIL="?"

echo "run=${RUN_ID} chunks=${DONE}/${PLANNED} pct=${PCT}% claims=${CLAIMS} phase=${PHASE} scorch=${SCORCH} consensus=${CONS} combos=${COMP}/${SCHED} failed=${FAIL}"
