#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_ambertools.sh [--prefix DIR] [--package SPEC] [--configure]

Create or update an AmberTools environment for Atlas MMGBSA workflows.
The default prefix is $ATLAS_TOOLS_DIR/envs/ambertools, or
../../tools/envs/ambertools relative to this repository.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
atlas_root="$(cd "$repo_root/../.." && pwd)"
tools_root="${ATLAS_TOOLS_DIR:-$atlas_root/tools}"
prefix="${MMGBSA_AMBERTOOLS_PREFIX:-${AMBERTOOLS_PREFIX:-$tools_root/envs/ambertools}}"
package_spec="${AMBERTOOLS_PACKAGE_SPEC:-ambertools=24.8}"
configure="${ATLAS_CONFIGURE_TOOLS:-0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      prefix="${2:?--prefix requires a value}"
      shift 2
      ;;
    --package)
      package_spec="${2:?--package requires a value}"
      shift 2
      ;;
    --configure)
      configure="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

cd "$repo_root"
mkdir -p "$(dirname "$prefix")"

if command -v micromamba >/dev/null 2>&1; then
  if [[ -d "$prefix" ]]; then
    micromamba install -y -p "$prefix" -c conda-forge "$package_spec"
  else
    micromamba create -y -p "$prefix" -c conda-forge "$package_spec"
  fi
elif command -v conda >/dev/null 2>&1; then
  if [[ -d "$prefix" ]]; then
    conda install -y -p "$prefix" -c conda-forge "$package_spec"
  else
    conda create -y -p "$prefix" -c conda-forge "$package_spec"
  fi
else
  echo "Neither micromamba nor conda was found on PATH." >&2
  exit 2
fi

required_tools=(MMPBSA.py cpptraj tleap antechamber parmchk2 sander MCPB.py)
missing=()
for tool in "${required_tools[@]}"; do
  if [[ ! -x "$prefix/bin/$tool" && ! -f "$prefix/bin/$tool" ]]; then
    missing+=("$tool")
  fi
done
if [[ ${#missing[@]} -gt 0 ]]; then
  echo "AmberTools prefix is missing expected tools: ${missing[*]}" >&2
  echo "Prefix checked: $prefix" >&2
  exit 3
fi

echo "AmberTools prefix: $prefix"
echo "AmberTools package: $package_spec"

if [[ "$configure" == "1" ]]; then
  tools/installers/configure_local_tools.py \
    --set "AMBERTOOLS_PREFIX=$prefix" \
    --set "MMGBSA_AMBERTOOLS_PREFIX=$prefix"
fi
