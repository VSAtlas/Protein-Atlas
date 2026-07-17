# -*- coding: utf-8 -*-
# ruff: noqa: E402

"""
Atlas2 multi-stage docking pipeline

Usage:
  atlas [OPTIONS] [PDB_IDS...]

Task-oriented commands:
  atlas init
      Create config.txt if needed and auto-detect local external tool paths.

  atlas demo
      Generate a tiny local report demo under outputs/data/atlas_demo/ without
      running external docking tools.

  atlas setup-report
      Summarize install readiness, detected tool paths, and optional/BYOL gaps.

  atlas new-run --pdb ID --small-library --fast --dry-run
      Validate first real-run inputs and print the command Atlas would launch.

  atlas status [RUN_ID]
      Summarize run lifecycle, stage progress, ETA, Slurm state when available,
      root-cause errors, and key output/report paths.

  atlas status RUN_ID --errors --explain
      Add structured root-cause buckets and plain-language next actions.

  atlas status RUN_ID --html
      Write outputs/data/<RUN_ID>/status.html as a standalone browser dashboard.

  atlas runs
      List recent run IDs, progress, failures, age, and report availability.

  atlas report RUN_ID
      Generate the canonical data/<RUN_ID>/report.html report artifacts.

  atlas targets guide
      Prompt through wetlab target gene-list/query/install setup.

  atlas targets from-genes EGFR ABL1 --out analysis/gene_list/targets.csv
      Query and rank candidate PDB structures from wetlab gene symbols.

  atlas ligands install hmdb
      Download/cache a supported ligand SDF source and prepare its PDBQT library
      with Meeko.

  atlas hmdb --pdb 1ABC --fast
      Convenience form: install HMDB if needed, then run the pipeline with the
      HMDB ligand library selected. `atlas fda ...` does the same for the FDA
      DrugCentral-backed library.

  atlas run [OPTIONS]
      Explicitly run the legacy-compatible docking pipeline. Existing direct
      forms like `atlas --pdb 1BN1 --fast` remain supported.

High-level description:
  Atlas runs the multi-stage docking pipeline for one or more proteins:
    1) Set up logging and run metadata
    2) Extract bound ligands and generate ligand-free PDBs
    3) Prepare receptors (clean PDB + PDBQT), reusing cached prep when possible
    4) Detect active sites and define docking boxes
    5) Prepare and filter ligands from configured libraries
    6) Run multi-stage docking with recentering and fallback heuristics
    7) Validate final poses
    8) Write per-protein score CSVs and update run_manifest.yaml

Run control:
  -h, --help
      Show this help message and exit.

  --run-id RUN_ID
  -run-id RUN_ID
      Explicit run identifier. By default, a timestamp like
      YYYYMMDD_HHMMSS is generated. Required when using -resume.
      Also overrides ATLAS_RUN_ID if set.

  -resume, --resume
      Resume a previous run using its run_id. The pipeline will:
        - Load run_manifest.yaml for that run
        - Reuse the stored configuration snapshot
        - Reconstruct the original argv for parsing
        - Skip proteins that were already marked as completed
      You MUST also pass --run-id <RUN_ID> (or set ATLAS_RUN_ID).

Inputs:
  --configs-dir PATH
      Override CONFIGS_DIR for this run. By default this is derived
      from OVERALL_DIR in the YAML config (e.g. OVERALL_DIR/configs).

  --no-reset-configs
      Do not reset the per-run CONFIG_RUN_DIR. By default, new runs
      reset their run-specific config directory so stale configs do
      not accumulate. In resume mode, configs are never reset.

Protein selection (Specified Proteins Mode):
  You can explicitly choose which PDBs to run via CLI, environment,
  or config. Precedence is: CLI > ENV > CFG. If nothing is specified,
  Atlas2 runs on all .pdb files under INPUT_DIR that do not contain
  '_nolig' in the filename.

  CLI forms (recommended first):
    --pdb ID
    -pdb ID
        Recommended: add one PDB ID per flag (repeatable). Example:
          --pdb 1BN1 --pdb 2OJ9

    --pdbs "ID1,ID2"
    -pdbs "ID1 ID2"
        Compatibility alias: comma- or space-separated list of PDB IDs.

    QoL flags:
      --1BN1, -1BN1
      -1BN1,2OJ9
        Compatibility alias: 4-character alnum short flags are treated as
        PDB IDs, including comma-packed short form; excludes FAST-like flags.

    Bare tokens:
      1BN1
      1BN1.pdb
      1BN1_cleaned.pdb
        If you pass bare arguments that look like 4-char PDB IDs or
        .pdb filenames, they are also treated as target proteins.

  Environment:
    ONLY_PDBS="1BN1 2OJ9"
        Space- or comma-separated list of PDB IDs to run.

  Config (YAML):
    SPECIFIED_PROTEINS:
      - 1BN1
      - 2OJ9
    or
    SPECIFIED_PROTEINS: "1BN1, 2OJ9"

  The final list of PDBs is printed as:
    [config] SPECIFIED_PROTEINS effective=[...]

Docking modes:
Fast mode:
  -fast, --fast, fast
      Enable FAST_MODE. This forces docking exhaustiveness=1 across
      stages for faster but less thorough runs. Helpful for smoke tests
      or CI checks.

Ligand/library modes:
Test FDA library toggle:
  -test-fda, --test-fda
      Switch the default small-molecule library from the full FDA set
      to the small test library configured in TEST_FDA_LIBRARY_SUBDIR
      (default: fda_test_library_10). This is designed for quick,
      lightweight test runs.

Benchmark Mode:
  -bench, --bench
      Enable reproducible BENCHMARK mode. This forces:
        - Only runs on PDBs: bNJS, bOJG, bNNQ
        - DUD test mode enabled (TEST_MODE_ENABLE=dud)
        - Maps these PDBs to specific bench libraries (bench_pur2, etc.)
        - Forces all engines enabled (GNINA, LEDOCK, DOCK6, SCORCH)
        - Forces APO_HOLO_MODE=holo
        - Forces PH_ENSEMBLE=True
      This ensures a standardized, comparable run configuration.
      Writes a 'bench_config.txt' snapshot to the config run dir.
      Cannot be combined with -resume.

Benchmark Mode 2:
  -bench2, --bench2
      Enable benchmark mode using the same PDBs, DUD test mode, and
      bench library mapping as -bench, but with Vina + SCORCH only:
        - Forces GNINA/LEDOCK/DOCK6 disabled; SCORCH enabled
        - Forces APO_HOLO_MODE=holo
        - Forces PH_ENSEMBLE=True
      Writes a 'bench2_config.txt' snapshot to the config run dir.
      Cannot be combined with -resume or -bench.

Benchmark Small:
  -bench-small, --bench-small
      Enable a smaller throughput profile for rapid scheduler iteration:
        - Uses PDBs bNJS, bOJG, bNNQ
        - Uses test_library_20 mapping for those proteins
        - Runs Vina + SCORCH only
        - Forces APO_HOLO_MODE=holo and PH_ENSEMBLE=True
        - Enables distributed combo-chunk assignment mode
      Writes a 'bench_small_config.txt' snapshot to the config run dir.
      Cannot be combined with -resume, -bench, or -bench2.

Benchmark Micro:
  -bench-micro, --bench-micro
      Enable the tiniest throughput profile for rapid scheduler iteration:
        - Uses bOJG by default with test_library_20
        - Caps ligand selection to a few ligands
        - Runs Vina + SCORCH only
        - Intended for cached-prep/local scheduler iteration, not final performance claims
      Cannot be combined with -resume, -bench, -bench2, or -bench-small.

Water Benchmark:
  -water-bench, --water-bench, --waterbench
      Enable a compact DUD-E-style water-policy profile for 3EML, 1UYG,
      and 1XL2 using water_aa2ar, water_hs90a, and water_hivpr libraries.
      Runs Vina + SCORCH only, caps each library at 999 ligands, and reads
      WATER_POLICY from ATLAS_WATER_BENCH_POLICY (default: site_only).
      Use remove_all, site_only, and keep_all runs to compare water effects.
      Cannot be combined with -resume, -bench, -bench2, or -bench-small.

DUD-Only Runtime Mode:
  -dude, --dude
      Force TEST_MODE_ENABLE=dud for this process only (no config file edits).
      This keeps only proteins listed in TEST_LIBRARY_MAP, while still honoring
      explicitly requested proteins from --pdb/--pdbs/ONLY_PDBS/SPECIFIED_PROTEINS.
      Cannot be combined with -bench, -bench2, -bench-small, -water-bench, or -resume.

DUD library selection:
  --dud-library LIBRARY_SUBDIR
      Use one DUD/decoy library subdirectory for every selected protein in this
      run without editing config.txt. This updates TEST_LIBRARY_MAP in memory
      and ensures TEST_MODE_ENABLE includes dud. Example:
        atlas --pdb 1BN1 --dud-library aa2ar
      For target-specific DUD-E mappings, keep using TEST_LIBRARY_MAP.

Single-ligand mode:
  --single PATTERN
      Enable SINGLE_LIGAND mode and restrict docking to a single ligand
      whose name contains PATTERN. Lookup policy is fixed in code:
        - search order: fda_library, per_protein
        - global search disabled
        - prefix matching disabled
        - manifest-only mode enabled
      When combined with multiple PDBs,
      the pipeline will dock the same ligand (if found) across all
      selected proteins in parallel.

pH/APO/HOLO options:
pH ensemble and pH-ligand modes:
  --ph-ligand-mode MODE
      Control how ligands are selected for pH-ensemble runs. Precedence:
        CLI (--ph-ligand-mode) > ENV (PH_LIGAND_MODE) > CFG (PH_LIGAND_MODE)
      Typical values:
        off       - standard behavior (no special pH filtering)
        context_window   - only dock ligands whose pH microstate is relevant
                    to the current receptor microstate/context


  Related environment variables:
    PH_ENSEMBLE (bool)
        Turn pH-ensemble mode on/off (overrides config).
    PH_SCOPE (string)
        Limit pH-ensemble to a subset of residues or region. If unset,
        defaults to config behavior.
    PH_RADIUS (float)
        Radius (Å) for pH-ensemble context. Non-positive values fall
        back to a large default.

No-library docking mode:
  --no-docking
      Enable NO_LIBRARY_DOCKING=True. In this mode, Atlas2 will:
        - Run receptor preparation
        - Run control redocking
        - Run pocket detection and ligand planning
      but it will SKIP docking for DUD/FDA libraries. This is useful
      when you only want control validation and planned ligand lists.

Debug/developer options:
Tool verification:
  --verify-tools
      Run startup tool verification using fixed "prefer_config" resolution
      (configured path first, PATH fallback second). This verification is
      CLI-only and does not auto-run from config defaults.
      Core required tools must resolve; optional checks (for example,
      SCORCH script/env) are reported but do not fail verification.
      When used, Atlas exits after verification.

Doctor diagnostics:
  doctor
  --doctor
      Run lightweight environment diagnostics and exit. This reports
      active Python executable, detected env name, key Python imports,
      and external tool path resolution.

Effective config:
  --print-effective-config
      Print the fully resolved runtime config as pretty JSON and exit.
      Resolution uses the same precedence as normal startup.

  -rebuild
      Rebuild/repair per-library `_manifest.json` files for all direct
      libraries under `prepped_ligands`, then exit without docking.

  --rebuild-root PATH
      Optional base repo path for `-rebuild` autodiscovery (for relocated
      installs). Example: /path/to/protein_automation (uses
      /path/to/protein_automation/prepped_ligands).

Reporting:
Artifact retention:
  --retain
      Force artifact retention on for this run.

  --noretain
      Force artifact retention off for this run.

  --retainmode MODE
      Force artifact retention mode. Allowed:
        rerun_safe | minimal_disk
      Cannot be combined with --noretain.

Apo/holo and variants:
  The apo/holo mode and variant list are resolved from the YAML
  configuration (APO_HOLO_MODE) and possibly env variables. This
  script will:
    - Iterate over all variants
    - Set APO_HOLO_VARIANT in the environment for each variant
    - Route all paths via the path_router module
  There is currently no direct CLI flag to override APO_HOLO_MODE; use
  the YAML config for that.

Environment variables (summary):
  ATLAS_RUN_ID
      Alternate way to supply the run ID. Used if --run-id is not
      provided. main() will still derive a timestamp if neither is set.

  ATLAS_DISTRIBUTED_MODE
      Set to "slurm_array" to enable distributed array sharding.
      Requires SLURM_ARRAY_TASK_ID and SLURM_ARRAY_TASK_COUNT (or
      ATLAS_DIST_TASK_ID/ATLAS_DIST_TASK_COUNT). In multi-task mode,
      Atlas uses combo-chunk distribution by default.

  ATLAS_DISABLE_COMBO_CHUNKS
      Optional distributed compatibility switch. Set to true/1/on to
      force legacy per-PDB sharding instead of combo-chunk distribution.

  PH_ENSEMBLE, PH_SCOPE, PH_RADIUS
      See pH ensemble section above.

  PH_LIGAND_MODE
      Fallback for --ph-ligand-mode if the CLI flag is not used.

  NO_LIBRARY_DOCKING
      Fallback for --no-docking when CLI is not used.

  ARTIFACT_RETENTION
      Artifact retention mode:
        rerun_safe | minimal_disk | off
      Precedence: CLI (--retain/--noretain/--retainmode)
                  > ENV (ARTIFACT_RETENTION)
                  > CFG (ARTIFACT_RETENTION)
                  > default (rerun_safe)

  SINGLE_LIGAND
      Fallback pattern for --single when no CLI arg is given.

  LIBRARY_MANIFEST_BUILD_ON_SCAN
      Boolean controlling whether per-library manifest files are built
      during directory scans.

Logging and outputs:
  - Logs:
      * A top-level main log is created under ./logs as:
          logs/main_<RUN_ID>.log
        In Slurm-array distributed mode:
          logs/main_<RUN_ID>_a<ARRAY_JOB_ID>_t<TASK_ID>.log
        Stdout/stderr are tee'd into that log.

  - Manifests:
      * run_manifest.yaml is created under manifests/<run-id> and tracks
        run configuration, apo/holo/water decisions, and per-protein
        status (start/success/failure).

  - Failed proteins:
      * Per-PDB failure logs are placed under:
          <OVERALL_DIR>/failed/<RUN_ID>__<PDB>__<variant>.log

  - Docking and configs:
      * Receptors, docking configs, and pose/score outputs are routed
        by path_router into subdirectories under:
          INPUT_DIR, PROCESSED_PDBS_DIR, DOCKED_DIR, PREPPED_LIGANDS_DIR,
          CONFIGS_DIR, etc., as defined in the YAML config.

Examples:
  # Standard run on all PDBs in INPUT_DIR
  atlas

  # Fast run on two specific PDBs
  atlas --pdb 1BN1 --pdb 2OJ9 --fast

  # Use small FDA test library and fast mode
  atlas -test-fda --fast

  # Resume a previous run by ID
  atlas -resume --run-id 20251205_225826

  # Single-ligand docking across multiple proteins
  atlas --pdbs "1BN1,2OJ9" --single imatinib

"""
from __future__ import annotations

# Canonical implementation for the legacy-compatible docking pipeline.
# Root main.py is a compatibility bridge only; keep new CLI/runtime seams in src/cli.

import json
import logging
import os
import sys
import time
import warnings
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Dict, Mapping

# Ensure local imports resolve when running from a source checkout without PYTHONPATH.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_SRC_ROOT = _REPO_ROOT / "src"
if str(_SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(_SRC_ROOT))
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from cli.cli_utils import (
    _cli_has,
    _cli_val,
)
from cli.local_scheduler_runtime import (
    run_global_scheduler,
    run_multi_pdb_single_ligand,
    run_serial_scheduler,
)
from cli.run_bootstrap import bootstrap_pre_variant_state, select_pdb_files_for_run
from cli.run_execution import (
    VariantExecutionContext,
    VariantExecutionDeps,
    VariantExecutionResult,
    VariantExecutionSharedState,
    run_variant_execution,
)
from cli.run_lifecycle import (
    finalize_deferred_retention,
    finalize_post_run,
    finalize_scheduler_and_retention,
)
from cli.process_lifecycle import process_lifecycle
from cli.run_context import (
    ConfigDict,
    _apply_resume_config_from_snapshot,
    _prepare_run_logfile,
    _resolve_run_id,
)
from cli.resume_manifest_repair import auto_repair_resume_manifest
from cli.run_process_one import (
    ProcessOneContext,
    ProcessOneDeps,
    ProcessOneSharedState,
    build_process_one_runner,
)
from cli import run_profiles as _run_profiles
from cli import run_tool_verification as _run_tool_verification
from cli import doctor as _doctor
from cli.distributed_chunk_planner import (
    _build_combo_work_items,
    _chunk_ligand_key,
    _load_post_scored_ligand_keys,
    _load_scored_ligand_keys_from_summary,
    _normalize_ph_tag_token,
    _resolve_combo_docking_summary_csv,
    _resolve_combo_output_dir,
    _resolve_combo_post_consensus_csv,
    _resolve_global_scheduler_plan,
    _update_combo_coverage_snapshot,
    _verify_chunk_combo_outputs,
)
from cli.planner_manifest import maybe_run_manifest_rebuild_only
from cli.debug_fs import install_debug_makedirs
from docking.fallback_recenter import load_recenter_params as get_recenter_params
from docking.global_scheduler import (
    acquire_global_slot,
)
from cli.run_bootstrap_helpers import smoke_emit_config_demo
from config.normalize import _to_bool
from config.runtime_config import load_inputs, validate_config
from cli.runtime_logging import (
    _tee_stdio_to,
    bootstrap_root_logging,
)
from cli.main_script_entrypoint import run_with_email_notification
from path_router.path_router import Paths as RouterPaths
import post_docking.mmgbsa.mmgbsa_pipeline as mmgbsa_pipeline
from cli.postrun_hooks_runtime import (
    _maybe_run_artifact_retention_for_combo,
    _maybe_run_scorch_rescore_for_pdb,
)
from post_docking.mmgbsa.mmgbsa_pipeline import (
    make_mmgbsa_trajectory,
    run_implicit_md,
    run_mmgbsa,
)
from cli.run_manifest_runtime import (
    apply_pocket_detection_events,
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
    update_manifest_for_protein_success,
)

_BOOTSTRAPPED = False


def run_distributed_combo_variant(*args: Any, **kwargs: Any) -> Any:
    from cli.distributed_chunk_runtime import run_distributed_combo_variant as impl

    return impl(*args, **kwargs)


def process_one_protein(*args: Any, **kwargs: Any) -> Any:
    from docking.docking import process_one_protein as impl

    return impl(*args, **kwargs)


def _bootstrap_runtime() -> None:
    """Apply runtime-only side effects once before pipeline execution."""
    global _BOOTSTRAPPED
    if _BOOTSTRAPPED:
        return

    warnings.filterwarnings(
        "ignore",
        message=(
            "to-Python converter for boost::shared_ptr<RDKit::FilterHierarchyMatcher> "
            "already registered; second conversion method ignored."
        ),
        category=RuntimeWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message=r".*joblib will operate in serial mode.*",
        category=UserWarning,
    )

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")

    install_debug_makedirs()
    _BOOTSTRAPPED = True


HELP_TEXT = __doc__ or "Atlas2 multi-stage docking pipeline"

def _print_effective_config_and_exit(cfg: ConfigDict) -> None:
    payload: dict[str, Any] = {}
    for key in sorted(cfg.keys()):
        k = str(key)
        if k.startswith("_"):
            continue
        payload[k] = cfg[key]
    print(json.dumps(payload, indent=2, sort_keys=True, default=str))
    sys.exit(0)


# Planner and distributed chunk helper functions moved to
# src/cli/distributed_chunk_planner.py.


@dataclass(frozen=True)
class _MainCliFlags:
    is_resume: bool
    cli_run_id: str | None
    is_bench: bool
    is_bench2: bool
    is_bench_micro: bool
    is_bench_small: bool
    is_water_bench: bool
    is_dude: bool
    cli_dud_library: str | None

    @property
    def bench_family(self) -> bool:
        return (
            self.is_bench
            or self.is_bench2
            or self.is_bench_micro
            or self.is_bench_small
            or self.is_water_bench
        )


def _apply_doctor_profile_overrides(cfg: ConfigDict, argv: list[str]) -> None:
    if _run_profiles._dude_enabled(argv):
        _run_profiles._apply_dude_overrides(cfg)
    if _run_profiles._bench_enabled(argv):
        _run_profiles._apply_bench_overrides(cfg)
    elif _run_profiles._bench2_enabled(argv):
        _run_profiles._apply_bench2_overrides(cfg)
    elif _run_profiles._bench_micro_enabled(argv):
        _run_profiles._apply_bench_micro_overrides(cfg)
    elif _run_profiles._bench_small_enabled(argv):
        _run_profiles._apply_bench_small_overrides(cfg)
    elif _run_profiles._water_bench_enabled(argv):
        _run_profiles._apply_water_bench_overrides(cfg)


def _maybe_exit_for_main_early_commands(argv: list[str]) -> None:
    wants_doctor = _doctor.cli_wants_doctor(argv)
    if not wants_doctor and (_cli_has(argv, "-h") or _cli_has(argv, "--help")):
        print(HELP_TEXT)
        sys.exit(0)

    if _cli_has(argv, "--print-effective-config"):
        cfg = ConfigDict(load_inputs())
        _print_effective_config_and_exit(cfg)

    if wants_doctor:
        cfg = ConfigDict(load_inputs())
        _apply_doctor_profile_overrides(cfg, argv)
        sys.exit(_doctor.run_doctor(cfg, argv))

    if _cli_has(argv, "--verify-tools"):
        cfg = ConfigDict(load_inputs())
        _run_tool_verification._verify_tools_if_requested(cfg, argv)
        sys.exit(0)


def _parse_main_cli_flags(argv: list[str]) -> _MainCliFlags:
    try:
        cli_dud_library = _run_profiles._dud_library_value(argv)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    return _MainCliFlags(
        is_resume=_cli_has(argv, "-resume") or _cli_has(argv, "--resume"),
        cli_run_id=_cli_val(argv, "--run-id") or _cli_val(argv, "-run-id"),
        is_bench=_run_profiles._bench_enabled(argv),
        is_bench2=_run_profiles._bench2_enabled(argv),
        is_bench_micro=_run_profiles._bench_micro_enabled(argv),
        is_bench_small=_run_profiles._bench_small_enabled(argv),
        is_water_bench=_run_profiles._water_bench_enabled(argv),
        is_dude=_run_profiles._dude_enabled(argv),
        cli_dud_library=cli_dud_library,
    )


_BENCH_MODE_ATTRS: tuple[tuple[str, str], ...] = (
    ("is_bench", "-bench"),
    ("is_bench2", "-bench2"),
    ("is_bench_micro", "-bench-micro"),
    ("is_bench_small", "-bench-small"),
    ("is_water_bench", "-water-bench"),
)


def _exit_cli_error(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(2)


def _validate_bench_mode_exclusivity(flags: _MainCliFlags) -> None:
    if sum(1 for attr, _ in _BENCH_MODE_ATTRS if getattr(flags, attr)) > 1:
        _exit_cli_error(
            "-bench, -bench2, -bench-micro, -bench-small, and -water-bench are mutually exclusive"
        )


def _validate_resume_compatibility(flags: _MainCliFlags) -> None:
    if not flags.is_resume:
        return
    for attr, label in _BENCH_MODE_ATTRS:
        if getattr(flags, attr):
            _exit_cli_error(f"{label} cannot be used with -resume")
    if not flags.cli_run_id:
        _exit_cli_error("-resume requires --run-id <RUN_ID>")


def _validate_dude_compatibility(flags: _MainCliFlags) -> None:
    if flags.is_dude and flags.bench_family:
        _exit_cli_error(
            "-dude cannot be used with -bench, -bench2, -bench-micro, -bench-small, or -water-bench"
        )
    if flags.cli_dud_library and flags.bench_family:
        _exit_cli_error(
            "--dud-library cannot be used with -bench, -bench2, -bench-micro, -bench-small, or -water-bench"
        )
    if flags.is_dude and flags.is_resume:
        _exit_cli_error("-dude cannot be used with -resume")


def _validate_main_cli_flags(flags: _MainCliFlags) -> None:
    _validate_bench_mode_exclusivity(flags)
    _validate_resume_compatibility(flags)
    _validate_dude_compatibility(flags)


def _apply_run_profile_overrides(cfg: ConfigDict, flags: _MainCliFlags) -> None:
    if flags.is_dude:
        _run_profiles._apply_dude_overrides(cfg)
        logging.info("[dude] DUD-only runtime mode ENABLED. Overrides applied.")
    if flags.is_bench:
        _run_profiles._apply_bench_overrides(cfg)
        logging.info("[bench] benchmark mode ENABLED. Overrides applied.")
    elif flags.is_bench2:
        _run_profiles._apply_bench2_overrides(cfg)
        logging.info("[bench2] benchmark2 mode ENABLED. Overrides applied (vina + scorch).")
    elif flags.is_bench_micro:
        _run_profiles._apply_bench_micro_overrides(cfg)
        logging.info(
            "[bench-micro] micro throughput profile ENABLED. Overrides applied (vina + scorch)."
        )
    elif flags.is_bench_small:
        _run_profiles._apply_bench_small_overrides(cfg)
        logging.info(
            "[bench-small] small throughput profile ENABLED. Overrides applied (vina + scorch)."
        )
    elif flags.is_water_bench:
        _run_profiles._apply_water_bench_overrides(cfg)
        logging.info(
            "[water-bench] water-policy DUD-E profile ENABLED. Overrides applied (vina + scorch)."
        )


def _run_main_pipeline_core(
    *,
    cfg: ConfigDict,
    run_id: str,
    argv_for_parsing: list[str],
    flags: _MainCliFlags,
    log_path: str,
    resume_protein_ids: list[str],
    completed_combo_lookup: set[tuple[str, str, str]],
) -> None:
    bootstrap_state = bootstrap_pre_variant_state(
        cfg=cfg,
        run_id=run_id,
        argv_for_parsing=argv_for_parsing,
        is_resume=flags.is_resume,
        is_bench=flags.is_bench,
        is_bench2=flags.is_bench2,
        is_bench_small=flags.is_bench_small or flags.is_bench_micro,
        resume_protein_ids=resume_protein_ids,
        completed_combo_lookup=completed_combo_lookup,
        log_path=log_path,
        bench_test_library_map=_run_profiles.BENCH_TEST_LIBRARY_MAP,
        bench_small_test_library_map=_run_profiles.BENCH_SMALL_TEST_LIBRARY_MAP,
        bench_small_target_ligands=_run_profiles.BENCH_SMALL_TARGET_LIGANDS,
        resolve_artifact_retention_mode=_run_profiles._resolve_artifact_retention_mode,
        select_pdb_files_for_run=select_pdb_files_for_run,
        get_recenter_params=get_recenter_params,
    )
    cfg = bootstrap_state.cfg
    argv_for_parsing = bootstrap_state.argv_for_parsing
    dist_ctx = bootstrap_state.dist_ctx
    completed_combo_lookup = bootstrap_state.completed_combo_lookup
    params = bootstrap_state.params
    stages = bootstrap_state.stages
    pdb_files = bootstrap_state.pdb_files
    tokens = bootstrap_state.tokens
    test_mode = bootstrap_state.test_mode
    dist_combo_chunk_mode = bootstrap_state.dist_combo_chunk_mode
    dist_combo_hybrid_mode = bootstrap_state.dist_combo_hybrid_mode
    execution_mode_decision = bootstrap_state.execution_mode_decision
    plan_only = bootstrap_state.plan_only
    mode = bootstrap_state.mode
    variants = bootstrap_state.variants
    start = bootstrap_state.start
    global_start = bootstrap_state.global_start
    failed_root = bootstrap_state.failed_root
    failed_entries = bootstrap_state.failed_entries
    run_scope_completion_times = bootstrap_state.run_scope_completion_times
    run_chunk_rebalance_count = bootstrap_state.run_chunk_rebalance_count
    recorded_terminal_chunk_ids = bootstrap_state.recorded_terminal_chunk_ids
    distributed_assigned_pdb_ids = bootstrap_state.distributed_assigned_pdb_ids
    retention_lock = bootstrap_state.retention_lock
    pending_retention_combos = bootstrap_state.pending_retention_combos
    pending_coverage_refresh = bootstrap_state.pending_coverage_refresh
    pending_retention_lock = bootstrap_state.pending_retention_lock
    scorch_queue_service = bootstrap_state.scorch_queue_service

    execution_context = VariantExecutionContext(
        run_id=str(run_id),
        cfg=cfg,
        variants=list(variants),
        pdb_files=list(pdb_files),
        tokens=list(tokens),
        mode=str(mode),
        is_resume=bool(flags.is_resume),
        dist_ctx=dist_ctx,
        dist_combo_chunk_mode=bool(dist_combo_chunk_mode),
        dist_combo_hybrid_mode=bool(dist_combo_hybrid_mode),
        execution_mode_decision=execution_mode_decision,
        stages=stages,
        params=params,
        global_start=float(global_start),
    )
    execution_shared_state = VariantExecutionSharedState(
        failed_root=str(failed_root),
        completed_combo_lookup=completed_combo_lookup,
        failed_entries=failed_entries,
        run_scope_completion_times=run_scope_completion_times,
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        distributed_assigned_pdb_ids=distributed_assigned_pdb_ids,
        run_chunk_rebalance_count=int(run_chunk_rebalance_count),
        retention_lock=retention_lock,
        pending_retention_combos=pending_retention_combos,
        pending_coverage_refresh=pending_coverage_refresh,
        pending_retention_lock=pending_retention_lock,
        scorch_queue_service=scorch_queue_service,
    )
    execution_deps = VariantExecutionDeps(
        to_bool=_to_bool,
        build_combo_work_items=_build_combo_work_items,
        resolve_global_scheduler_plan=_resolve_global_scheduler_plan,
        build_process_one_runner=build_process_one_runner,
        process_one_context_type=ProcessOneContext,
        process_one_shared_state_type=ProcessOneSharedState,
        process_one_deps_type=ProcessOneDeps,
        normalize_ph_tag_token=_normalize_ph_tag_token,
        chunk_ligand_key=_chunk_ligand_key,
        process_one_protein=process_one_protein,
        update_manifest_for_protein_start=update_manifest_for_protein_start,
        update_manifest_for_protein_success=update_manifest_for_protein_success,
        update_manifest_for_protein_failure=update_manifest_for_protein_failure,
        verify_chunk_combo_outputs=_verify_chunk_combo_outputs,
        resolve_combo_output_dir=_resolve_combo_output_dir,
        resolve_combo_docking_summary_csv=_resolve_combo_docking_summary_csv,
        load_scored_ligand_keys_from_summary=_load_scored_ligand_keys_from_summary,
        update_combo_coverage_snapshot=_update_combo_coverage_snapshot,
        resolve_combo_post_consensus_csv=_resolve_combo_post_consensus_csv,
        load_post_scored_ligand_keys=_load_post_scored_ligand_keys,
        acquire_global_slot=acquire_global_slot,
        maybe_run_scorch_rescore_for_pdb=_maybe_run_scorch_rescore_for_pdb,
        maybe_run_artifact_retention_for_combo=_maybe_run_artifact_retention_for_combo,
        run_distributed_combo_variant=run_distributed_combo_variant,
        run_multi_pdb_single_ligand=run_multi_pdb_single_ligand,
        run_global_scheduler=run_global_scheduler,
        run_serial_scheduler=run_serial_scheduler,
        apply_pocket_detection_events=apply_pocket_detection_events,
    )
    execution_result: VariantExecutionResult = run_variant_execution(
        context=execution_context,
        shared_state=execution_shared_state,
        deps=execution_deps,
    )
    run_chunk_rebalance_count = int(execution_result.run_chunk_rebalance_count)

    scorch_queue_failed, retention_barrier_ok = finalize_scheduler_and_retention(
        cfg=cfg,
        run_id=str(run_id),
        global_start=float(global_start),
        run_scope_completion_times=run_scope_completion_times,
        run_chunk_rebalance_count=int(run_chunk_rebalance_count),
        dist_ctx=dist_ctx,
        scorch_queue_service=scorch_queue_service,
        pending_coverage_refresh=pending_coverage_refresh,
        pending_retention_combos=pending_retention_combos,
        retention_lock=retention_lock,
    )
    elapsed_min = (time.time() - start) / 60.0
    print(f"\nAll proteins processed in {elapsed_min:.2f} minutes.")
    if scorch_queue_failed > 0:
        print(f"\nSCORCH queue reported {scorch_queue_failed} failed queued jobs.")

    if failed_entries:
        print("\nThe following proteins failed. See per-PDB logs under:", failed_root)
        for pdb_id, label, failed_log_path, exc_type, exc_msg in failed_entries:
            print(
                f"  - {pdb_id} ({label}): {exc_type} â€” {exc_msg}\n"
                f"      log: {failed_log_path}"
            )
    else:
        print("\nNo proteins recorded as failed.")

    should_stop = finalize_post_run(
        cfg=cfg,
        run_id=str(run_id),
        global_start=float(global_start),
        failed_entries=failed_entries,
        plan_only=bool(plan_only),
        variants=list(variants),
        pdb_files=list(pdb_files),
        test_mode=str(test_mode),
        dist_ctx=dist_ctx,
        dist_combo_chunk_mode=bool(dist_combo_chunk_mode),
        distributed_assigned_pdb_ids=distributed_assigned_pdb_ids,
        is_bench=bool(flags.is_bench),
        is_bench2=bool(flags.is_bench2),
        is_bench_small=bool(flags.is_bench_small or flags.is_bench_micro),
    )
    if plan_only:
        sys.exit(0)
    if scorch_queue_failed > 0:
        sys.exit(1)
    if should_stop:
        return
    finalize_deferred_retention(
        cfg=cfg,
        run_id=str(run_id),
        retention_barrier_ok=bool(retention_barrier_ok),
        pending_retention_combos=pending_retention_combos,
        retention_lock=retention_lock,
    )


# >>> PATHS CLASS START
Paths = RouterPaths
# >>> PATHS CLASS END


# --- Single-ligand ---
def _mmgbsa_five_replicate_runner(
    topo: Dict[str, object],
    out_dir: Path,
    cfg: Mapping[str, Any],
    force: bool,
    md_enabled: bool,
    stage_dir: str,
    ligand_stem: str,
    pdb_id: str,
    variant_dir: str,
    ph_label: str,
    run_id: str,
    logger: logging.Logger,
) -> Dict[str, object]:
    # Refresh mmgbsa_pipeline globals so test monkeypatches on main.* apply.
    mmgbsa_pipeline.run_implicit_md = run_implicit_md
    mmgbsa_pipeline.run_mmgbsa = run_mmgbsa
    mmgbsa_pipeline.make_mmgbsa_trajectory = make_mmgbsa_trajectory
    return mmgbsa_pipeline._mmgbsa_five_replicate_runner(
        topo=topo,
        out_dir=out_dir,
        cfg=cfg,
        force=force,
        md_enabled=md_enabled,
        stage_dir=stage_dir,
        ligand_stem=ligand_stem,
        pdb_id=pdb_id,
        variant_dir=variant_dir,
        ph_label=ph_label,
        run_id=run_id,
        logger=logger,
    )


# ======================
# Program entry point
# ======================
def main() -> None:
    _bootstrap_runtime()
    _maybe_exit_for_main_early_commands(sys.argv)

    print("MODELLER is working with license.")
    flags = _parse_main_cli_flags(sys.argv)
    _validate_main_cli_flags(flags)

    if flags.bench_family:
        os.environ["PH_ENSEMBLE"] = "1"

    run_id = _resolve_run_id(sys.argv)
    argv_for_parsing = list(sys.argv)
    resume_protein_ids: list[str] = []
    completed_combo_lookup: set[tuple[str, str, str]] = set()
    os.environ["ATLAS_RUN_ID"] = run_id
    log_path = _prepare_run_logfile(run_id)
    os.environ["ATLAS_LOG_FILE"] = log_path
    _tee_stdio_to(log_path)
    cfg = ConfigDict(load_inputs())
    cfg["_DUD_LIBRARY_OVERRIDE"] = flags.cli_dud_library

    _apply_run_profile_overrides(cfg, flags)

    bootstrap_root_logging(cfg, log_path)
    logging.info("[probe.root] root-logger INFO now visible")
    print(f"[run] log_file={log_path} run_id={run_id}")

    if flags.is_resume:
        auto_repair_resume_manifest(cfg, run_id)
        snap_dict = _apply_resume_config_from_snapshot(cfg, run_id)
        cfg = ConfigDict(snap_dict)
        cfg["_DUD_LIBRARY_OVERRIDE"] = flags.cli_dud_library
        cfg["RUN_ID"] = run_id
    else:
        cfg["RUN_ID"] = run_id

    validate_config(cfg)
    if maybe_run_manifest_rebuild_only(cfg, argv_for_parsing):
        return
    _run_tool_verification._verify_tools_if_requested(cfg, argv_for_parsing)

    with process_lifecycle(run_id=str(run_id), is_resume=bool(flags.is_resume)):
        _run_main_pipeline_core(
            cfg=cfg,
            run_id=str(run_id),
            argv_for_parsing=argv_for_parsing,
            flags=flags,
            log_path=log_path,
            resume_protein_ids=resume_protein_ids,
            completed_combo_lookup=completed_combo_lookup,
        )


def _run_with_email_notification() -> None:
    run_with_email_notification(
        bootstrap_runtime=_bootstrap_runtime,
        smoke_emit_config_demo=smoke_emit_config_demo,
        main_func=main,
    )


if __name__ == "__main__":
    _run_with_email_notification()
