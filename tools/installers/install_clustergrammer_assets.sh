#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_clustergrammer_assets.sh [--prefix DIR] [--check-only] [--force] [--no-pointer]

Install Clustergrammer browser assets from the pinned NPM metadata under
tools/report_assets/npm. Assets are installed outside the source tree by default
under $ATLAS_TOOLS_DIR/report_assets and exposed to reports through the
repo-local report_assets pointer file.

Flags:
  --prefix DIR    Asset root that will contain clustergrammer/ (default: $ATLAS_TOOLS_DIR/report_assets)
  --check-only    Validate the installed asset bundle and pointer without installing
  --force         Reinstall even when the current asset bundle validates
  --no-pointer    Do not create/update the repo-local report_assets pointer file
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
atlas_root="$(cd "$repo_root/../.." && pwd)"
tools_root="${ATLAS_TOOLS_DIR:-$atlas_root/tools}"
prefix="${ATLAS_REPORT_ASSETS_ROOT:-$tools_root/report_assets}"
npm_metadata_dir="$repo_root/tools/report_assets/npm"
check_only="0"
force="0"
write_pointer="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      prefix="${2:?--prefix requires a value}"
      shift 2
      ;;
    --check-only)
      check_only="1"
      shift
      ;;
    --force)
      force="1"
      shift
      ;;
    --no-pointer)
      write_pointer="0"
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
prefix="$(python -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).expanduser().resolve())' "$prefix")"
asset_root="$prefix/clustergrammer"
build_dir="$prefix/_npm_build"
pointer_path="$repo_root/report_assets"

required_assets=(
  "clustergrammer.js"
  "css/custom.css"
  "lib/css/bootstrap.css"
  "lib/css/font-awesome.min.css"
  "lib/js/bootstrap.min.js"
  "lib/js/d3.js"
  "lib/js/jquery-1.11.2.min.js"
  "lib/js/underscore-min.js"
)

validate_assets() {
  local root="$1"
  local rel=""
  for rel in "${required_assets[@]}"; do
    if [[ ! -s "$root/$rel" ]]; then
      echo "Missing Clustergrammer asset: $root/$rel" >&2
      return 1
    fi
  done
  return 0
}

validate_pointer() {
  if [[ "$write_pointer" != "1" ]]; then
    return 0
  fi
  if [[ -d "$pointer_path" ]]; then
    if validate_assets "$pointer_path/clustergrammer"; then
      return 0
    fi
    echo "Existing report_assets directory does not contain a valid clustergrammer bundle: $pointer_path" >&2
    return 1
  fi
  if [[ ! -f "$pointer_path" ]]; then
    echo "Missing report_assets pointer file: $pointer_path" >&2
    return 1
  fi
  local target_root
  target_root="$(tr -d '\r\n' < "$pointer_path")"
  if [[ -z "$target_root" ]]; then
    echo "Empty report_assets pointer file: $pointer_path" >&2
    return 1
  fi
  validate_assets "$target_root/clustergrammer"
}

if [[ "$check_only" == "1" ]]; then
  validate_assets "$asset_root"
  validate_pointer
  echo "Clustergrammer assets are installed: $asset_root"
  exit 0
fi

if [[ "$force" != "1" ]] && validate_assets "$asset_root" >/dev/null 2>&1; then
  echo "Clustergrammer assets already installed: $asset_root"
else
  command -v npm >/dev/null 2>&1 || {
    echo "npm is required to install Clustergrammer browser assets." >&2
    echo "Run the core Atlas installer first, or install Node.js/npm on PATH." >&2
    exit 2
  }
  [[ -f "$npm_metadata_dir/package-lock.json" ]] || {
    echo "Missing pinned NPM lockfile: $npm_metadata_dir/package-lock.json" >&2
    exit 2
  }

  mkdir -p "$build_dir"
  cp "$npm_metadata_dir/package.json" "$build_dir/package.json"
  cp "$npm_metadata_dir/package-lock.json" "$build_dir/package-lock.json"
  npm ci \
    --omit=dev \
    --ignore-scripts \
    --no-audit \
    --no-fund \
    --cache "$prefix/_npm_cache" \
    --prefix "$build_dir"

  source_root="$build_dir/node_modules/clustergrammer"
  [[ -d "$source_root" ]] || {
    echo "NPM install did not produce $source_root" >&2
    exit 2
  }

  tmp_root="${asset_root}.tmp.$$"
  rm -rf "$tmp_root"
  mkdir -p "$tmp_root"
  cp -R "$source_root/." "$tmp_root/"
  validate_assets "$tmp_root"
  rm -rf "$asset_root"
  mv "$tmp_root" "$asset_root"
  rm -rf "$build_dir"
  echo "Clustergrammer assets installed: $asset_root"
fi

if [[ "$write_pointer" == "1" ]]; then
  if [[ -d "$pointer_path" ]]; then
    validate_assets "$pointer_path/clustergrammer"
    echo "Using existing report_assets directory: $pointer_path"
  else
    printf '%s\n' "$prefix" > "$pointer_path"
    echo "Report assets pointer updated: $pointer_path -> $prefix"
  fi
fi

validate_assets "$asset_root"
validate_pointer
