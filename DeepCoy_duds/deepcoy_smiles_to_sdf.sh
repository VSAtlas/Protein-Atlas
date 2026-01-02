#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 <OUTPUT_DIR>" >&2
}

if [ $# -ne 1 ]; then
    usage
    exit 2
fi

OUTPUT_DIR="$1"
SMILES_FILE="${OUTPUT_DIR}/deepcoy_decoys.smi"
SDF_FILE="${OUTPUT_DIR}/deepcoy_decoys.sdf"
OBABEL_EXE="${OBABEL_EXE:-obabel}"

echo "Converting SMILES in: ${OUTPUT_DIR}"

if [ ! -s "${SMILES_FILE}" ]; then
    echo "ERROR: Output file not found or empty at: ${SMILES_FILE}" >&2
    exit 1
fi

if [[ "${OBABEL_EXE}" == */* ]]; then
    if [ ! -x "${OBABEL_EXE}" ]; then
        echo "ERROR: Open Babel executable not found or not executable: ${OBABEL_EXE}" >&2
        exit 1
    fi
else
    if ! command -v "${OBABEL_EXE}" >/dev/null 2>&1; then
        echo "ERROR: Open Babel executable not found in PATH. Set OBABEL_EXE or install obabel." >&2
        exit 1
    fi
fi

"${OBABEL_EXE}" \
    -i smi "${SMILES_FILE}" \
    -o sdf -O "${SDF_FILE}" \
    --gen2D

if [ ! -s "${SDF_FILE}" ]; then
    echo "ERROR: Conversion completed but output is missing or empty: ${SDF_FILE}" >&2
    exit 1
fi

echo "Success: Created ${SDF_FILE}"
