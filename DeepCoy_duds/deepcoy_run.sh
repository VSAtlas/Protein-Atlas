#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 <ACTIVES_SMI_PATH> <OUTPUT_DIR> [--artifact-dir DIR] [--decoys-per-active N] [--config BASE_CONFIG] [--deepcoy-python PY] [--train-cache PATH] [--valid-cache PATH]" >&2
}

if [ $# -lt 2 ]; then
    usage
    exit 2
fi

ACTIVES_SMI_PATH="$1"
OUTPUT_DIR="$2"
shift 2

DECOYS_PER_ACTIVE=50
BASE_CONFIG=""
DEEPCOY_PYTHON_ARG=""
TRAIN_CACHE_PATH=""
VALID_CACHE_PATH=""
ARTIFACT_DIR_ARG=""

while [ $# -gt 0 ]; do
    case "$1" in
        --decoys-per-active)
            if [ $# -lt 2 ]; then
                echo "ERROR: --decoys-per-active requires a value." >&2
                exit 2
            fi
            DECOYS_PER_ACTIVE="$2"
            shift 2
            ;;
        --config)
            if [ $# -lt 2 ]; then
                echo "ERROR: --config requires a value." >&2
                exit 2
            fi
            BASE_CONFIG="$2"
            shift 2
            ;;
        --artifact-dir)
            if [ $# -lt 2 ]; then
                echo "ERROR: --artifact-dir requires a value." >&2
                exit 2
            fi
            ARTIFACT_DIR_ARG="$2"
            shift 2
            ;;
        --deepcoy-python)
            if [ $# -lt 2 ]; then
                echo "ERROR: --deepcoy-python requires a value." >&2
                exit 2
            fi
            DEEPCOY_PYTHON_ARG="$2"
            shift 2
            ;;
        --train-cache)
            if [ $# -lt 2 ]; then
                echo "ERROR: --train-cache requires a value." >&2
                exit 2
            fi
            TRAIN_CACHE_PATH="$2"
            shift 2
            ;;
        --valid-cache)
            if [ $# -lt 2 ]; then
                echo "ERROR: --valid-cache requires a value." >&2
                exit 2
            fi
            VALID_CACHE_PATH="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 2
            ;;
    esac
done

DEEPCOY_REPO_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
resolve_artifact_dir() {
    local provided="$1"
    if [ -n "$provided" ]; then
        if [[ "$provided" = /* ]]; then
            echo "$provided"
            return
        else
            echo "${DEEPCOY_REPO_PATH}/${provided}"
            return
        fi
    fi
    if [ -n "${DEEPCOY_ARTIFACT_DIR:-}" ]; then
        if [[ "${DEEPCOY_ARTIFACT_DIR}" = /* ]]; then
            echo "${DEEPCOY_ARTIFACT_DIR}"
        else
            echo "${DEEPCOY_REPO_PATH}/${DEEPCOY_ARTIFACT_DIR}"
        fi
        return
    fi
    echo "${DEEPCOY_REPO_PATH}/deepcoy_work/UNKNOWN"
}

ARTIFACT_DIR="$(resolve_artifact_dir "$ARTIFACT_DIR_ARG")"
mkdir -p "$ARTIFACT_DIR"

DEEPCOY_PYTHON_PATH="$DEEPCOY_PYTHON_ARG"
if [ -z "$DEEPCOY_PYTHON_PATH" ]; then
    DEEPCOY_PYTHON_PATH="${DEEPCOY_PYTHON:-}"
fi

if [ -z "$DEEPCOY_PYTHON_PATH" ]; then
    echo "ERROR: DeepCoy python interpreter not set. Use --deepcoy-python or set DEEPCOY_PYTHON." >&2
    exit 2
fi

if [ ! -x "$DEEPCOY_PYTHON_PATH" ]; then
    echo "ERROR: DeepCoy python interpreter not executable: $DEEPCOY_PYTHON_PATH" >&2
    exit 2
fi

if [ -z "$BASE_CONFIG" ]; then
    BASE_CONFIG="${DEEPCOY_REPO_PATH}/deepcoy_config.json"
else
    if [ ! -f "$BASE_CONFIG" ] && [ -f "${DEEPCOY_REPO_PATH}/${BASE_CONFIG}" ]; then
        BASE_CONFIG="${DEEPCOY_REPO_PATH}/${BASE_CONFIG}"
    fi
fi

if [ ! -f "$BASE_CONFIG" ]; then
    echo "ERROR: Base config not found: $BASE_CONFIG" >&2
    exit 2
fi

if [ ! -s "$ACTIVES_SMI_PATH" ]; then
    echo "ERROR: Actives SMILES file not found or empty: $ACTIVES_SMI_PATH" >&2
    exit 2
fi

mkdir -p "$OUTPUT_DIR"

export DEEPCOY_DIR="${DEEPCOY_REPO_PATH}"
export ATLAS_ROOT="${DEEPCOY_REPO_PATH}"

ACTIVES_JSON="${OUTPUT_DIR}/actives_valid.json"
RUNTIME_CONFIG="${ARTIFACT_DIR}/deepcoy_runtime.json"
OUTPUT_SMI="${OUTPUT_DIR}/deepcoy_decoys.smi"
DEFAULT_TRAIN_CACHE="${DEEPCOY_REPO_PATH}/data/cache/molecules_train_zinc.preprocessed.pkl"
DEFAULT_VALID_CACHE="${DEEPCOY_REPO_PATH}/data/cache/molecules_valid_zinc.preprocessed.pkl"

echo "Starting DeepCoy run for: ${ACTIVES_SMI_PATH}"
echo "Repo path detected as: ${DEEPCOY_REPO_PATH}"
echo "Artifact directory: ${ARTIFACT_DIR}"

"${DEEPCOY_PYTHON_PATH}" "${DEEPCOY_REPO_PATH}/smiles_to_valid_json.py" \
    --in-smi "${ACTIVES_SMI_PATH}" \
    --out-json "${ACTIVES_JSON}" \
    --dataset zinc

if [ -z "$TRAIN_CACHE_PATH" ] && [ -f "$DEFAULT_TRAIN_CACHE" ]; then
    TRAIN_CACHE_PATH="$DEFAULT_TRAIN_CACHE"
fi
if [ -z "$VALID_CACHE_PATH" ] && [ -f "$DEFAULT_VALID_CACHE" ]; then
    VALID_CACHE_PATH="$DEFAULT_VALID_CACHE"
fi

RESTORE_PATH="$("${DEEPCOY_PYTHON_PATH}" - "${BASE_CONFIG}" "${ACTIVES_JSON}" "${OUTPUT_SMI}" "${DECOYS_PER_ACTIVE}" "${RUNTIME_CONFIG}" "${TRAIN_CACHE_PATH}" "${VALID_CACHE_PATH}" <<'PY'
import json
import os
import sys

base_path, valid_file, output_name, decoys_per_active, out_path, train_cache, valid_cache = sys.argv[1:]
with open(base_path, "r") as f:
    config = json.load(f)

base_valid_file = config.get("valid_file")
config["generation"] = True
config["valid_file"] = valid_file
config["output_name"] = output_name
config["number_of_generation_per_valid"] = int(decoys_per_active)
config["num_epochs"] = 1
config["epoch_to_generate"] = 1
try:
    with open(valid_file, "r") as vf:
        valid_data = json.load(vf)
    valid_count = len(valid_data) if isinstance(valid_data, list) else 0
except Exception:
    valid_count = 0
base_batch_size = int(config.get("batch_size", 8))
if valid_count > 0:
    config["batch_size"] = max(1, min(base_batch_size, valid_count))
if train_cache and os.path.isfile(train_cache):
    config["train_cache"] = train_cache
if valid_cache and os.path.isfile(valid_cache) and valid_file == base_valid_file:
    config["valid_cache"] = valid_cache

with open(out_path, "w") as f:
    json.dump(config, f, indent=2)

restore = config.get("restore", "")
if restore:
    print(restore)
PY
)"

DEEPCOY_CMD=( "${DEEPCOY_PYTHON_PATH}" "${DEEPCOY_REPO_PATH}/DeepCoy.py" --dataset zinc --config-file "${RUNTIME_CONFIG}" --log_dir "${ARTIFACT_DIR}" )

if [ -n "$RESTORE_PATH" ]; then
    if [[ "$RESTORE_PATH" != /* ]]; then
        RESTORE_PATH="${DEEPCOY_REPO_PATH}/${RESTORE_PATH}"
    fi
    if [ ! -f "$RESTORE_PATH" ]; then
        echo "ERROR: Restore checkpoint not found: $RESTORE_PATH" >&2
        exit 2
    fi
    DEEPCOY_CMD+=( --restore "$RESTORE_PATH" )
fi

( cd "${ARTIFACT_DIR}" && "${DEEPCOY_CMD[@]}" )

if [ ! -s "$OUTPUT_SMI" ]; then
    echo "ERROR: DeepCoy output missing or empty: $OUTPUT_SMI" >&2
    exit 1
fi

if ! find "$ARTIFACT_DIR" -maxdepth 1 -type f \( -name "*_params_*.json" -o -name "*_log_*.json" \) | head -n 1 | grep -q .; then
    echo "ERROR: Expected DeepCoy params/log artifacts not found in ${ARTIFACT_DIR}" >&2
    exit 1
fi

echo "DeepCoy execution completed successfully."
