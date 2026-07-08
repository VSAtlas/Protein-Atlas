#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_adfrsuite.sh [--archive FILE | --existing-root DIR] [--prefix DIR] [--configure]

Install or register a local ADFRsuite/MGLTools tree. Atlas does not redistribute
ADFRsuite or MGLTools; provide your own archive or an already installed root.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
atlas_root="$(cd "$repo_root/../.." && pwd)"
tools_root="${ATLAS_TOOLS_DIR:-$atlas_root/tools}"
prefix="${ADFRSUITE_PREFIX:-$tools_root/adfrsuite}"
archive="${ADFRSUITE_ARCHIVE:-}"
existing_root="${ADFRSUITE_ROOT:-}"
configure="${ATLAS_CONFIGURE_TOOLS:-0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --archive)
      archive="${2:?--archive requires a value}"
      shift 2
      ;;
    --existing-root)
      existing_root="${2:?--existing-root requires a value}"
      shift 2
      ;;
    --prefix)
      prefix="${2:?--prefix requires a value}"
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

if [[ -n "$existing_root" ]]; then
  install_root="$existing_root"
elif [[ -n "$archive" ]]; then
  mkdir -p "$prefix"
  case "$archive" in
    *.tar.gz|*.tgz|*.tar)
      tar -xf "$archive" -C "$prefix"
      ;;
    *.zip)
      command -v unzip >/dev/null 2>&1 || { echo "unzip is required for zip archives" >&2; exit 2; }
      unzip -q "$archive" -d "$prefix"
      ;;
    *)
      echo "Unsupported archive type: $archive" >&2
      exit 2
      ;;
  esac
  install_root="$prefix"
else
  echo "Provide --archive FILE or --existing-root DIR for ADFRsuite/MGLTools." >&2
  echo "This installer is BYOL and does not download restricted binaries." >&2
  exit 2
fi

pythonsh="$(find "$install_root" -maxdepth 5 -type f -name pythonsh -perm -u+x 2>/dev/null | head -n 1 || true)"
prepare_ligand="$(find "$install_root" -maxdepth 7 -type f -name prepare_ligand4.py 2>/dev/null | head -n 1 || true)"
prepare_receptor="$(find "$install_root" -maxdepth 7 -type f -name prepare_receptor4.py 2>/dev/null | head -n 1 || true)"

if [[ -z "$pythonsh" || -z "$prepare_ligand" || -z "$prepare_receptor" ]]; then
  echo "Could not locate pythonsh, prepare_ligand4.py, and prepare_receptor4.py under $install_root" >&2
  exit 2
fi

echo "ADFRsuite/MGLTools root: $install_root"
echo "pythonsh: $pythonsh"
echo "prepare_ligand4.py: $prepare_ligand"
echo "prepare_receptor4.py: $prepare_receptor"

if [[ "$configure" == "1" ]]; then
  tools/installers/configure_local_tools.py \
    --set "MGLTOOLS_PATH=$install_root" \
    --set "MGLTOOLS_PYTHON=$pythonsh" \
    --set "PREPARE_LIGAND_SCRIPT=$prepare_ligand" \
    --set "PREPARE_RECEPTOR_SCRIPT=$prepare_receptor"
fi
