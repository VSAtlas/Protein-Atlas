#!/usr/bin/env bash
# Merge Cursor agent worktree branches back into main (integrator helper).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

BRANCHES=(
  ci-ruff
  ci-mypy-gate
  ci-mypy-analysis
  ci-vulture
  ci-xenon
  qol-test-init
  qol-test-status
  qol-test-doctor
  qol-test-help-dev
)

echo "[merge] repo=${REPO_ROOT} target=main"
git checkout main

for branch in "${BRANCHES[@]}"; do
  if ! git show-ref --verify --quiet "refs/heads/${branch}"; then
    echo "[skip] no branch ${branch}"
    continue
  fi
  echo "[merge] ${branch} -> main"
  if ! git merge --no-edit "${branch}"; then
    echo "[FAIL] merge conflict on ${branch}; resolve manually then re-run from this branch"
    exit 1
  fi
done

echo "[merge] complete; run: atlas dev verify --full --fix --smoke"
