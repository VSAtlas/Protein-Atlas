#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
export PYTHONPATH="${repo_root}/src:${repo_root}:${PYTHONPATH:-}"
python_bin="${ATLAS_PYTHON:-python3}"

echo "[public-smoke] import atlas CLI entrypoint module"
"$python_bin" -c "from cli.atlas_main_cli import main as _m; print('atlas-cli-import-ok')" >/dev/null

if [[ ! -f config.txt ]]; then
  cp config.example.txt config.txt
  echo "[public-smoke] created config.txt from config.example.txt"
fi

echo "[public-smoke] generate demo report"
"$python_bin" -c "import sys; from cli.qol_cli import dispatch; sys.argv = ['atlas', 'demo', '--force', '--status-html']; raise SystemExit(dispatch())" >/dev/null

echo "[public-smoke] complete"
