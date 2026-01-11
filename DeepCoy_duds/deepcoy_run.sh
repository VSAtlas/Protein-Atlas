#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 <ACTIVES_SMI_PATH> <OUTPUT_DIR> [--artifact-dir DIR] [--decoys-per-active N] [--config BASE_CONFIG] [--deepcoy-python PY] [--dataset NAME] [--train-cache PATH] [--valid-cache PATH] [--restrict-data N] [--restore-model PATH] [--random-seed N] [--use-argmax-generation BOOL] [--try-different-starting BOOL] [--num-different-starting N] [--num-samples N]" >&2
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
DATASET="zinc"
TRAIN_CACHE_PATH=""
VALID_CACHE_PATH=""
ARTIFACT_DIR_ARG=""
RESTRICT_DATA=""
RESTORE_MODEL_OVERRIDE=""
RANDOM_SEED=""
USE_ARGMAX_GENERATION=""
TRY_DIFFERENT_STARTING=""
NUM_DIFFERENT_STARTING=""
NUM_SAMPLES=""

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
        --dataset)
            if [ $# -lt 2 ]; then
                echo "ERROR: --dataset requires a value." >&2
                exit 2
            fi
            DATASET="$2"
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
        --restrict-data)
            if [ $# -lt 2 ]; then
                echo "ERROR: --restrict-data requires a value." >&2
                exit 2
            fi
            RESTRICT_DATA="$2"
            shift 2
            ;;
        --restore-model)
            if [ $# -lt 2 ]; then
                echo "ERROR: --restore-model requires a value." >&2
                exit 2
            fi
            RESTORE_MODEL_OVERRIDE="$2"
            shift 2
            ;;
        --random-seed)
            if [ $# -lt 2 ]; then
                echo "ERROR: --random-seed requires a value." >&2
                exit 2
            fi
            RANDOM_SEED="$2"
            shift 2
            ;;
        --use-argmax-generation)
            if [ $# -lt 2 ]; then
                echo "ERROR: --use-argmax-generation requires a value (true/false)." >&2
                exit 2
            fi
            USE_ARGMAX_GENERATION="$2"
            shift 2
            ;;
        --try-different-starting)
            if [ $# -lt 2 ]; then
                echo "ERROR: --try-different-starting requires a value (true/false)." >&2
                exit 2
            fi
            TRY_DIFFERENT_STARTING="$2"
            shift 2
            ;;
        --num-different-starting)
            if [ $# -lt 2 ]; then
                echo "ERROR: --num-different-starting requires a value." >&2
                exit 2
            fi
            NUM_DIFFERENT_STARTING="$2"
            shift 2
            ;;
        --num-samples)
            if [ $# -lt 2 ]; then
                echo "ERROR: --num-samples requires a value." >&2
                exit 2
            fi
            NUM_SAMPLES="$2"
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
echo "Dataset: ${DATASET}"

"${DEEPCOY_PYTHON_PATH}" "${DEEPCOY_REPO_PATH}/smiles_to_valid_json.py" \
    --in-smi "${ACTIVES_SMI_PATH}" \
    --out-json "${ACTIVES_JSON}" \
    --dataset "${DATASET}"

if [ -z "$TRAIN_CACHE_PATH" ] && [ -f "$DEFAULT_TRAIN_CACHE" ]; then
    TRAIN_CACHE_PATH="$DEFAULT_TRAIN_CACHE"
fi
if [ -z "$VALID_CACHE_PATH" ] && [ -f "$DEFAULT_VALID_CACHE" ]; then
    VALID_CACHE_PATH="$DEFAULT_VALID_CACHE"
fi

RESTORE_PATH="$("${DEEPCOY_PYTHON_PATH}" - "${BASE_CONFIG}" "${ACTIVES_JSON}" "${OUTPUT_SMI}" "${DECOYS_PER_ACTIVE}" "${RUNTIME_CONFIG}" "${TRAIN_CACHE_PATH}" "${VALID_CACHE_PATH}" "${RESTORE_MODEL_OVERRIDE}" "${RANDOM_SEED}" "${USE_ARGMAX_GENERATION}" "${TRY_DIFFERENT_STARTING}" "${NUM_DIFFERENT_STARTING}" "${NUM_SAMPLES}" <<'PY'
import json
import os
import sys
from pathlib import Path

(
    base_path,
    valid_file,
    output_name,
    decoys_per_active,
    out_path,
    train_cache,
    valid_cache,
    restore_override,
    random_seed,
    use_argmax_generation,
    try_different_starting,
    num_different_starting,
    num_samples,
) = sys.argv[1:]
with open(base_path, "r") as f:
    config = json.load(f)

def _coerce_bool(val, default):
    if val is None or val == "":
        return default
    val_norm = str(val).strip().lower()
    if val_norm in {"1", "true", "yes", "on"}:
        return True
    if val_norm in {"0", "false", "no", "off"}:
        return False
    return default

def _coerce_int(val, default):
    try:
        return int(val)
    except Exception:
        return default

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

base_dir = Path(base_path).resolve().parent
default_restore = config.get("restore", "")
if restore_override:
    restore_path = Path(restore_override)
    if not restore_path.is_absolute():
        restore_path = base_dir / restore_path
    config["restore"] = str(restore_path)
else:
    if default_restore:
        restore_path = Path(default_restore)
        if not restore_path.is_absolute():
            restore_path = base_dir / restore_path
        config["restore"] = str(restore_path)

seed_default = config.get("random_seed", 0)
if random_seed:
    config["random_seed"] = _coerce_int(random_seed, seed_default)
elif seed_default is not None:
    config["random_seed"] = seed_default

config["use_argmax_generation"] = _coerce_bool(
    use_argmax_generation, config.get("use_argmax_generation", False)
)
config["try_different_starting"] = _coerce_bool(
    try_different_starting, config.get("try_different_starting", True)
)
config["num_different_starting"] = _coerce_int(
    num_different_starting, config.get("num_different_starting", 1)
)
num_samples_default = config.get("number_of_generation_per_valid", config.get("num_samples", 50))
config["num_samples"] = _coerce_int(num_samples, num_samples_default)

with open(out_path, "w") as f:
    json.dump(config, f, indent=2)

restore = config.get("restore", "")
if restore:
    print(restore)
PY
)"

if [ -z "$RESTRICT_DATA" ] && [ -n "${DEEPCOY_RESTRICT_DATA:-}" ]; then
    RESTRICT_DATA="${DEEPCOY_RESTRICT_DATA}"
fi

DEEPCOY_CMD=( "${DEEPCOY_PYTHON_PATH}" "${DEEPCOY_REPO_PATH}/DeepCoy.py" --dataset "${DATASET}" --config-file "${RUNTIME_CONFIG}" --log_dir "${ARTIFACT_DIR}" )

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
if [ -n "$RESTRICT_DATA" ]; then
    DEEPCOY_CMD+=( --restrict_data "$RESTRICT_DATA" )
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
