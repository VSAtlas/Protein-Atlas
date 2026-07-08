#!/usr/bin/env bash
# Merge simplify-round worktree branches into main.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

WAVE="${1:-wave1}"

if [[ "${WAVE}" == "wave1" ]]; then
  BRANCHES=(
    simp-01-report
    simp-02-heatmap
    simp-03-dudeval
    simp-09-status
    simp-04-qol
    simp-05-manifest
    simp-06-postrun
    simp-07-dist-plan
    simp-08-dist-run
  )
elif [[ "${WAVE}" == "wave2" ]]; then
  BRANCHES=(
    simp-10-schema
    simp-11-ml
    simp-12-scorch
  )
else
  echo "Unknown wave: ${WAVE}" >&2
  exit 2
fi

echo "[merge-simplify] repo=${REPO_ROOT} wave=${WAVE}"
git checkout main

for branch in "${BRANCHES[@]}"; do
  if ! git show-ref --verify --quiet "refs/heads/${branch}"; then
    echo "[skip] no branch ${branch}"
    continue
  fi
  echo "[merge] ${branch} -> main"
  git merge --no-edit "${branch}" || {
    echo "[FAIL] conflict on ${branch}; resolve and re-run from this branch"
    exit 1
  }
done

echo "[merge-simplify] done; run: atlas dev verify --full --fix --smoke"
