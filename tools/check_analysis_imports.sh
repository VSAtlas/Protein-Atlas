#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

# Legacy wrapper imports that should not be introduced in new code.
pattern='(from analysis import (run_report|master_schema_export|heatmap_html|manifest_utils|fda_name_map|pocket_second_pass|convert_interactions_to_parquet|summarize_run_errors|dud_eval_(types|log|metrics|schema|labels|discovery|reporting|engine|orchestrate))|import analysis\.(run_report|master_schema_export|heatmap_html|manifest_utils|fda_name_map|pocket_second_pass|convert_interactions_to_parquet|summarize_run_errors|dud_eval_(types|log|metrics|schema|labels|discovery|reporting|engine|orchestrate)))'

if rg -n -P "$pattern" \
  analysis src ml chemdb/tests postrun_hooks.py main.py \
  -g '!analysis/*.py' -g '!**/__pycache__/**'; then
  echo ""
  echo "[import-check] Found legacy analysis wrapper imports. Use canonical modules under analysis.cli, analysis.dud_eval_core, analysis.reporting."
  exit 1
fi

echo "[import-check] OK"
