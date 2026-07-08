#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_p2rank.sh [--archive FILE | --url URL | --existing-root DIR] [--prefix DIR] [--configure]

Install or register P2Rank. Pass an explicit archive or URL so publication
workflows can pin the exact upstream release they used.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
atlas_root="$(cd "$repo_root/../.." && pwd)"
tools_root="${ATLAS_TOOLS_DIR:-$atlas_root/tools}"
prefix="${P2RANK_PREFIX:-$tools_root/p2rank}"
archive="${P2RANK_ARCHIVE:-}"
url="${P2RANK_URL:-}"
existing_root="${P2RANK_ROOT:-}"
configure="${ATLAS_CONFIGURE_TOOLS:-0}"
cache_dir="${ATLAS_INSTALLER_CACHE:-$HOME/.cache/atlas_installers}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --archive)
      archive="${2:?--archive requires a value}"
      shift 2
      ;;
    --url)
      url="${2:?--url requires a value}"
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
else
  if [[ -n "$url" && -z "$archive" ]]; then
    mkdir -p "$cache_dir"
    archive="$cache_dir/$(basename "$url")"
    if command -v curl >/dev/null 2>&1; then
      curl -fL "$url" -o "$archive"
    elif command -v wget >/dev/null 2>&1; then
      wget -O "$archive" "$url"
    else
      echo "curl or wget is required to download --url" >&2
      exit 2
    fi
  fi

  if [[ -z "$archive" ]]; then
    if command -v prank >/dev/null 2>&1; then
      install_root="$(dirname "$(dirname "$(command -v prank)")")"
    else
      echo "Provide --archive FILE, --url URL, or --existing-root DIR for P2Rank." >&2
      exit 2
    fi
  else
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
  fi
fi

launcher="$(find "$install_root" -maxdepth 5 \( -type f -o -type l \) -name prank 2>/dev/null | head -n 1 || true)"
if [[ -z "$launcher" && -x "$install_root/prank" ]]; then
  launcher="$install_root/prank"
fi
if [[ -z "$launcher" ]]; then
  echo "Could not locate P2Rank launcher named prank under $install_root" >&2
  exit 2
fi
chmod +x "$launcher" 2>/dev/null || true

echo "P2Rank root: $install_root"
echo "prank launcher: $launcher"

if [[ "$configure" == "1" ]]; then
  tools/installers/configure_local_tools.py --set "P2RANK_PATH=$launcher"
fi
