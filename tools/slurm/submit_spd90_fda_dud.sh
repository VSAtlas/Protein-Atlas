#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RUN_ID="${RUN_ID:-spd90_fda_fda_dud_spr10_40h_$(date +%Y%m%d_%H%M%S)}"
NODES="${NODES:-10}"
WALLTIME="${WALLTIME:-40:00:00}"
PARTITION="${PARTITION:-spr}"
ACCOUNT="${ACCOUNT:-MCB24041}"
JOB_NAME="${JOB_NAME:-atlas_spd90_fda_dud}"
SPD_PDB_FILE="${SPD_PDB_FILE:-analysis/gene_list/spd_targets_dry90_pdbs.txt}"
PREPPED_ROOT="${PREPPED_ROOT:-/scratch/11040/mpg2352/atlas/code/protein_automation/prepped_ligands}"
ALL_DIRS="${ALL_DIRS:-/scratch/11040/mpg2352/atlas/runs/${RUN_ID}}"
DUD_LIBRARY="${DUD_LIBRARY:-fda_dud}"
MAIL_USER="${MAIL_USER:-}"
MAIL_TYPE="${MAIL_TYPE:-BEGIN,END,FAIL}"

cd "${REPO_ROOT}"
if [[ ! -f "${SPD_PDB_FILE}" ]]; then
  echo "Missing SPD_PDB_FILE=${SPD_PDB_FILE}" >&2
  exit 2
fi

ATLAS_BIN="${ATLAS_BIN:-}"
if [[ -z "${ATLAS_BIN}" ]]; then
  ATLAS_BIN="$(command -v atlas || true)"
fi
if [[ -z "${ATLAS_BIN}" && -x "/work2/11040/mpg2352/.mamba/envs/docking-env/bin/atlas" ]]; then
  ATLAS_BIN="/work2/11040/mpg2352/.mamba/envs/docking-env/bin/atlas"
fi
if [[ -z "${ATLAS_BIN}" ]]; then
  echo "Unable to find atlas. Activate docking-env or set ATLAS_BIN=/path/to/atlas." >&2
  exit 127
fi

SPD_PDBS="$(tr '\n' ' ' < "${SPD_PDB_FILE}" | sed 's/[[:space:]]*$//')"

export ATLAS_REPO_ROOT="${REPO_ROOT}"
export ATLAS_PREPPED_ROOT="${PREPPED_ROOT}"
export ATLAS_ALL_DIRS="${ALL_DIRS}"
export ATLAS_SKIP_EXTRACT=1
export TEST_MODE_ENABLE="dud+fda"
export USE_SCORCH=true
export SCORCH_TOP_FRACTION=0.10
export ATLAS_SCORCH_SHARDS_ENABLE="${ATLAS_SCORCH_SHARDS_ENABLE:-1}"
export USE_GNINA="${USE_GNINA:-false}"
export USE_LEDOCK="${USE_LEDOCK:-false}"
export USE_DOCK6="${USE_DOCK6:-false}"

if [[ -n "${MAIL_USER}" ]]; then
  export SBATCH_MAIL_USER="${MAIL_USER}"
  export SBATCH_MAIL_TYPE="${MAIL_TYPE}"
fi

MAIL_ARGS=()
if [[ -n "${MAIL_USER}" ]]; then
  MAIL_ARGS+=(--mail-user "${MAIL_USER}" --mail-type "${MAIL_TYPE}")
fi

"${ATLAS_BIN}" slurm submit \
  --run-id "${RUN_ID}" \
  --array "0-$((NODES - 1))%${NODES}" \
  --cpus-per-task auto \
  --script tools/slurm/run_relocated_spr.sh \
  --finalize-script tools/slurm/run_relocated_finalize.sh \
  --partition "${PARTITION}" \
  --account "${ACCOUNT}" \
  --time "${WALLTIME}" \
  --job-name "${JOB_NAME}" \
  --exclusive-node \
  --execution-mode distributed_combo_hybrid \
  --all-dirs "${ALL_DIRS}" \
  --main-args "${SPD_PDBS} --dud-library ${DUD_LIBRARY} --run-id {run_id}" \
  "${MAIL_ARGS[@]}"
