#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <deepcoy_output_name_or_path> [label]" >&2
  echo "  - If you pass a name (e.g. 'atlas_cah2_control_decoys')," >&2
  echo "    the script will look for tools/DeepCoy/output/<name>.txt" >&2
  echo "  - If you pass a path (e.g. 'output/foo.txt' or '/full/path/foo.txt')," >&2
  echo "    that file will be used directly." >&2
  echo "  - Optional [label] overrides the library name used for folder/file names." >&2
  exit 1
fi

# Directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Atlas root:
#   Default: two levels up from this script (…/atlas)
#   Override: export ATLAS_ROOT=/path/to/atlas_root
ATLAS_ROOT="${ATLAS_ROOT:-"$(cd "$SCRIPT_DIR/../.." && pwd)"}"

DEEPCOY_OUTPUT_DIR="${ATLAS_ROOT}/tools/DeepCoy/output"
BASE_EXTRACT_DIR="${ATLAS_ROOT}/code/protein_automation/extracted_ligands/deepcoy"

INPUT="$1"
LABEL_OVERRIDE="${2:-}"

# Resolve INPUT to an absolute DeepCoy output file path
if [[ "$INPUT" == */* ]]; then
  # Looks like a path (relative or absolute)
  if [[ "$INPUT" = /* ]]; then
    OUT_PATH="$INPUT"
  else
    # strip leading "output/" if present and anchor in tools/DeepCoy/output
    OUT_PATH="${DEEPCOY_OUTPUT_DIR}/${INPUT#output/}"
  fi
else
  # Looks like a bare name, assume tools/DeepCoy/output/<name>.txt
  OUT_PATH="${DEEPCOY_OUTPUT_DIR}/${INPUT}.txt"
fi

if [[ ! -f "$OUT_PATH" ]]; then
  echo "Error: DeepCoy output file not found at '$OUT_PATH'." >&2
  exit 1
fi

# Library name / label
if [[ -n "$LABEL_OVERRIDE" ]]; then
  LABEL="$LABEL_OVERRIDE"
else
  LABEL="$(basename "$OUT_PATH" .txt)"
fi

# Per-library folder:
#   .../extracted_ligands/deepcoy/<LABEL>/
LIB_DIR="${BASE_EXTRACT_DIR}/${LABEL}"
mkdir -p "$LIB_DIR"

DECOY_SMI="${LIB_DIR}/${LABEL}.smi"
DECOY_SDF="${LIB_DIR}/${LABEL}.sdf"

echo "Using DeepCoy output: $OUT_PATH"
echo "Library label:        $LABEL"
echo "Writing decoy SMILES to: $DECOY_SMI"
echo "Writing decoy SDF to:    $DECOY_SDF"

# DeepCoy format: "<active_smiles> <decoy_smiles>"
# We keep only the decoy (2nd column) and assign a name per line.
awk -v label="$LABEL" 'NF>=2 {printf "%s %s_%06d\n", $2, label, NR}' "$OUT_PATH" > "$DECOY_SMI"

# Convert SMILES to SDF with 3D coordinates
obabel -ismi "$DECOY_SMI" -osdf -O "$DECOY_SDF" --gen3D

echo "Done."
echo "Outputs are in: $LIB_DIR"
