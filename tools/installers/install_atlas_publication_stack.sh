#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_atlas_publication_stack.sh [FLAGS]

Flags:
  --core             Install the core Atlas environment.
  --adfrsuite        Install/register ADFRsuite or MGLTools. Requires its own flags/env.
  --p2rank           Install/register P2Rank. Requires --p2rank-archive/--p2rank-url or env.
  --scorch           Install/register SCORCH. Requires --scorch-source/--scorch-git-url or env.
  --ambertools       Install/register AmberTools under the Atlas tools root.
  --banana           Install/register BANANA under the Atlas tools root.
  --report-assets    Install Clustergrammer browser assets for offline reports.
  --no-report-assets Skip Clustergrammer browser asset installation.
  --configure        Detect local tools and update config.txt.
  --public-smoke     Run the no-proprietary public smoke/demo report check.
  --no-public-smoke  Skip the default public smoke/demo report check.
  --no-shim          Do not install/update the outer atlas environment shim.
  --all              Run all optional tiers. Restricted/BYOL tiers still need sources.
  --verify-tools     Run atlas --verify-tools after selected tiers.

Pass-through flags:
  --adfrsuite-archive FILE   --adfrsuite-root DIR
  --p2rank-archive FILE      --p2rank-url URL      --p2rank-root DIR
  --scorch-source DIR        --scorch-git-url URL
  --ambertools-prefix DIR     --ambertools-package SPEC
  --banana-source DIR         --banana-existing-root DIR
  --banana-git-url URL        --banana-git-ref REF
  --banana-root DIR           --banana-env-prefix DIR

With no flags, this runs --core --report-assets --public-smoke. The core
installer creates config.txt and auto-detects open/legal tools from the active
environment.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
installer_dir="$repo_root/tools/installers"
run_core="0"
run_adfrsuite="0"
run_p2rank="0"
run_scorch="0"
run_ambertools="0"
run_banana="0"
run_report_assets="0"
run_configure="0"
run_public_smoke="0"
verify_tools="0"
install_shim="${ATLAS_INSTALL_SHIM:-1}"
adfrsuite_args=()
p2rank_args=()
scorch_args=()
ambertools_args=()
banana_args=()

if [[ $# -eq 0 ]]; then
  run_core="1"
  run_report_assets="1"
  run_public_smoke="1"
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --core)
      run_core="1"
      run_report_assets="1"
      shift
      ;;
    --adfrsuite)
      run_adfrsuite="1"
      shift
      ;;
    --p2rank)
      run_p2rank="1"
      shift
      ;;
    --scorch)
      run_scorch="1"
      shift
      ;;
    --ambertools)
      run_ambertools="1"
      shift
      ;;
    --banana)
      run_banana="1"
      shift
      ;;
    --report-assets)
      run_report_assets="1"
      shift
      ;;
    --no-report-assets)
      run_report_assets="0"
      shift
      ;;
    --configure)
      run_configure="1"
      shift
      ;;
    --public-smoke)
      run_report_assets="1"
      run_public_smoke="1"
      shift
      ;;
    --no-public-smoke)
      run_public_smoke="0"
      shift
      ;;
    --no-shim)
      install_shim="0"
      shift
      ;;
    --all)
      run_core="1"
      run_adfrsuite="1"
      run_p2rank="1"
      run_scorch="1"
      run_ambertools="1"
      run_banana="1"
      run_report_assets="1"
      run_configure="1"
      run_public_smoke="1"
      shift
      ;;
    --verify-tools)
      verify_tools="1"
      shift
      ;;
    --adfrsuite-archive)
      adfrsuite_args+=(--archive "${2:?--adfrsuite-archive requires a value}")
      shift 2
      ;;
    --adfrsuite-root)
      adfrsuite_args+=(--existing-root "${2:?--adfrsuite-root requires a value}")
      shift 2
      ;;
    --p2rank-archive)
      p2rank_args+=(--archive "${2:?--p2rank-archive requires a value}")
      shift 2
      ;;
    --p2rank-url)
      p2rank_args+=(--url "${2:?--p2rank-url requires a value}")
      shift 2
      ;;
    --p2rank-root)
      p2rank_args+=(--existing-root "${2:?--p2rank-root requires a value}")
      shift 2
      ;;
    --scorch-source)
      scorch_args+=(--source-dir "${2:?--scorch-source requires a value}")
      shift 2
      ;;
    --scorch-git-url)
      scorch_args+=(--git-url "${2:?--scorch-git-url requires a value}")
      shift 2
      ;;
    --ambertools-prefix)
      ambertools_args+=(--prefix "${2:?--ambertools-prefix requires a value}")
      shift 2
      ;;
    --ambertools-package)
      ambertools_args+=(--package "${2:?--ambertools-package requires a value}")
      shift 2
      ;;
    --banana-source)
      banana_args+=(--source-dir "${2:?--banana-source requires a value}")
      shift 2
      ;;
    --banana-existing-root)
      banana_args+=(--existing-root "${2:?--banana-existing-root requires a value}")
      shift 2
      ;;
    --banana-git-url)
      banana_args+=(--git-url "${2:?--banana-git-url requires a value}")
      shift 2
      ;;
    --banana-git-ref)
      banana_args+=(--git-ref "${2:?--banana-git-ref requires a value}")
      shift 2
      ;;
    --banana-root)
      banana_args+=(--prefix "${2:?--banana-root requires a value}")
      shift 2
      ;;
    --banana-env-prefix)
      banana_args+=(--env-prefix "${2:?--banana-env-prefix requires a value}")
      shift 2
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

run_setup_report() {
  local env_name="${ATLAS_ENV_NAME:-docking-env}"
  if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n "$env_name" atlas setup-report || true
  elif command -v conda >/dev/null 2>&1; then
    conda run -n "$env_name" atlas setup-report || true
  elif command -v atlas >/dev/null 2>&1; then
    atlas setup-report || true
  else
    echo "atlas setup-report unavailable until the Atlas CLI is on PATH"
  fi
}

run_report_assets_installer() {
  local env_name="${ATLAS_ENV_NAME:-docking-env}"
  if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n "$env_name" bash "$installer_dir/install_clustergrammer_assets.sh"
  elif command -v conda >/dev/null 2>&1; then
    conda run -n "$env_name" bash "$installer_dir/install_clustergrammer_assets.sh"
  else
    bash "$installer_dir/install_clustergrammer_assets.sh"
  fi
}

print_success_panel() {
  local env_name="${ATLAS_ENV_NAME:-docking-env}"
  local shim_note="${ATLAS_SHIM_DIR:-auto-selected writable PATH dir, fallback $HOME/.local/bin}"
  local default_tools_root
  local report_assets_root
  default_tools_root="$(cd "$repo_root/../.." && pwd)/tools"
  report_assets_root="${ATLAS_REPORT_ASSETS_ROOT:-${ATLAS_TOOLS_DIR:-$default_tools_root}/report_assets}"
  echo
  echo "Atlas install summary"
  echo "  repo:   $repo_root"
  echo "  env:    $env_name"
  echo "  config: $repo_root/config.txt"
  echo "  assets: $report_assets_root"
  if [[ "$install_shim" == "1" ]]; then
    echo "  shim:   $shim_note"
  fi
  echo
  run_setup_report
  echo
  echo "Ready next steps"
  echo "  atlas smoke public"
  echo "  atlas demo --status-html"
  echo "  atlas init --write-example-inputs"
  echo "  micromamba activate $env_name  # optional; the atlas shim runs one command in the env"
  echo
  echo "Not bundled: proprietary or bring-your-own-license tools."
  echo "Register ADFR/MGLTools, P2Rank, SCORCH, AmberTools, BANANA, or Modeller only if your workflow needs them."
}

if [[ "$run_core" == "1" ]]; then
  ATLAS_INSTALL_SHIM="$install_shim" "$installer_dir/install_atlas_core.sh"
fi
if [[ "$run_adfrsuite" == "1" ]]; then
  "$installer_dir/install_adfrsuite.sh" --configure "${adfrsuite_args[@]}"
fi
if [[ "$run_p2rank" == "1" ]]; then
  "$installer_dir/install_p2rank.sh" --configure "${p2rank_args[@]}"
fi
if [[ "$run_scorch" == "1" ]]; then
  "$installer_dir/install_scorch.sh" --configure "${scorch_args[@]}"
fi
if [[ "$run_ambertools" == "1" ]]; then
  "$installer_dir/install_ambertools.sh" --configure "${ambertools_args[@]}"
fi
if [[ "$run_banana" == "1" ]]; then
  "$installer_dir/install_banana.sh" "${banana_args[@]}"
fi
if [[ "$run_report_assets" == "1" ]]; then
  run_report_assets_installer
fi
if [[ "$run_configure" == "1" ]]; then
  "$installer_dir/configure_local_tools.py"
fi
if [[ "$run_public_smoke" == "1" ]]; then
  if command -v micromamba >/dev/null 2>&1; then
    env_name="${ATLAS_ENV_NAME:-docking-env}"
    micromamba run -n "$env_name" bash tools/public_smoke_check.sh
  elif command -v conda >/dev/null 2>&1; then
    env_name="${ATLAS_ENV_NAME:-docking-env}"
    conda run -n "$env_name" bash tools/public_smoke_check.sh
  else
    bash tools/public_smoke_check.sh
  fi
fi
if [[ "$verify_tools" == "1" ]]; then
  env_name="${ATLAS_ENV_NAME:-docking-env}"
  if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n "$env_name" atlas --verify-tools
  elif command -v conda >/dev/null 2>&1; then
    conda run -n "$env_name" atlas --verify-tools
  else
    atlas --verify-tools
  fi
fi

print_success_panel
