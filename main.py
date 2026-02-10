# -*- coding: utf-8 -*-
from __future__ import annotations

HELP_TEXT = """
Atlas2 multi-stage docking pipeline

Usage:
  python main.py [OPTIONS] [PDB_IDS...]

High-level description:
  This script runs the Atlas2 pipeline for one or more proteins:
    1) Set up logging and run metadata
    2) Extract bound ligands and generate ligand-free PDBs
    3) Prepare receptors (clean PDB + PDBQT), reusing cached prep when possible
    4) Detect active sites and define docking boxes
    5) Prepare and filter ligands from configured libraries
    6) Run multi-stage docking with recentering and fallback heuristics
    7) Validate final poses and optionally capture screenshots
    8) Write per-protein score CSVs and update run_manifest.yaml

Core options:
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

Config directory and manifests:
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

  CLI forms:
    --pdb ID
    -pdb ID
        Add a single PDB ID (repeatable). Example:
          --pdb 1BN1 --pdb 2OJ9

    --pdbs "ID1,ID2"
    -pdbs "ID1 ID2"
        Comma- or space-separated list of PDB IDs.

    QoL flags:
      --1BN1, -1BN1
        Any 4-character alphanumeric token used as a flag (--XXXX or -XXXX)
        is treated as a PDB ID, except special tokens like FAST.

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

Fast mode:
  -fast, --fast, fast
      Enable FAST_MODE. This forces docking exhaustiveness=1 across
      stages for faster but less thorough runs. Helpful for smoke tests
      or CI checks.

Test FDA library toggle:
  -test-fda, --test-fda
      Switch the default small-molecule library from the full FDA set
      to the small test library configured in TEST_FDA_LIBRARY_SUBDIR
      (default: fda_test_library_10). This is designed for quick,
      lightweight test runs.

Calibrator refresh:
  --force-calibrator, -force-calibrator
      Force regenerating calibrator ligands/labels by clearing any cached
      calibrator artifacts before pocket evaluation.

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

DUD-Only Runtime Mode:
  -dude, --dude
      Force TEST_MODE_ENABLE=dud for this process only (no config file edits).
      This keeps only proteins listed in TEST_LIBRARY_MAP, while still honoring
      explicitly requested proteins from --pdb/--pdbs/ONLY_PDBS/SPECIFIED_PROTEINS.
      Cannot be combined with -bench, -bench2, or -resume.

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
  python main.py

  # Fast run on two specific PDBs
  python main.py --pdb 1BN1 --pdb 2OJ9 --fast

  # Use small FDA test library and fast mode
  python main.py -test-fda --fast

  # Resume a previous run by ID
  python main.py -resume --run-id 20251205_225826

  # Single-ligand docking across multiple proteins
  python main.py --pdbs "1BN1,2OJ9" --single imatinib

"""
from pdb_fixer import get_atom_rules
import sys
import logging
import time
import os
import shutil
import shlex
import subprocess
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Any, Mapping
from tqdm import tqdm
import traceback
from logging_topics import (
    _tee_stdio_to,
    bootstrap_root_logging,
    ensure_file_handler,
)
from run_manifest import (
    apply_pocket_detection_events,
    finalize_run_manifest,
    init_run_manifest,
    load_run_manifest,
    update_manifest_for_config_hash,
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
    update_manifest_for_protein_success,
    update_manifest_for_run_config,
    update_manifest_for_scheduled_proteins,
)

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
from input_and_export_functions import (
    load_inputs,
    validate_config,
    define_docking_stages,
    _to_bool,
    init_config_run_dir,
)
from cli.cli_utils import (
    _cli_val,
    _cli_has,
    _parse_single_from_cli,
    _parse_fast_flag,
    _parse_specified_proteins,
    _norm_pdb_id,
)
from cli.run_context import (
    _resolve_run_id,
    _prepare_run_logfile,
    _apply_resume_config_from_snapshot,
    ConfigDict,
)
from docking.docking_vina import emit_vina_config as _emit_vina_config_impl
from docking.docking_ligands import _resolve_test_mode
from docking.library_mode import _coerce_test_map, parse_test_libraries
from docking.docking import process_one_protein
from path_router.path_router import (
    make_paths,
    Paths as RouterPaths,
    receptor_file,
    docked_dir,
    config_dir as router_config_dir,
    run_logs_dir,
)
from prep_ligands.deepcoy_integration import (
    apply_deepcoy_cli_overrides,
    apply_deepcoy_defaults,
)
import post_docking.mmgbsa.mmgbsa_pipeline as mmgbsa_pipeline
from post_docking.mmgbsa.mmgbsa_pipeline import (
    _maybe_run_mmgbsa_for_pdb,
    make_mmgbsa_trajectory,
    run_implicit_md,
    run_mmgbsa,
)
from docking.fallback_recenter import RecenterParams
from apo_holo_mode import (
    resolve_apo_holo_mode,
    _debug_normalize_mode_token,
)
from debug_fs import install_debug_makedirs
from postrun_hooks import (
    _maybe_run_dud_eval,
    _maybe_run_scorch_rescore_for_pdb,
    _log_rescore_verification,
    _maybe_run_master_schema_export,
    _maybe_run_report_generation,
    _maybe_run_artifact_retention_for_pdb,
)

# Install debug wrappers for Path.mkdir and os.makedirs at import time,
# preserving the previous behavior.
install_debug_makedirs()


# --- Benchmark Constants & Helpers ---
BENCH_PDB_IDS = ["bNJS", "bOJG", "bNNQ"]
BENCH_TEST_LIBRARY_MAP = {
    "bNJS": "bench_pur2",
    "bOJG": "bench_mk01",
    "bNNQ": "bench_fabp4",
}


def _bench_enabled(argv: list[str]) -> bool:
    """Check if -bench or --bench is present in arguments."""
    return _cli_has(argv, "-bench") or _cli_has(argv, "--bench")


def _bench2_enabled(argv: list[str]) -> bool:
    """Check if -bench2 or --bench2 is present in arguments."""
    return _cli_has(argv, "-bench2") or _cli_has(argv, "--bench2")


def _dude_enabled(argv: list[str]) -> bool:
    """Check if -dude or --dude is present in arguments."""
    return _cli_has(argv, "-dude") or _cli_has(argv, "--dude")


def _normalize_artifact_retention_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "rerunsafe", "rerun_safe"}:
        return "rerun_safe"
    if token in {"minimaldisk", "minimal_disk"}:
        return "minimal_disk"
    if token in {"off", "none", "false", "0", "disabled", "disable"}:
        return "off"
    return "rerun_safe"


def _resolve_artifact_retention_mode(
    cfg: Mapping[str, Any], argv: list[str]
) -> tuple[str, str]:
    cli_retain = _cli_has(argv, "--retain")
    cli_noretain = _cli_has(argv, "--noretain")
    cli_has_mode = _cli_has(argv, "--retainmode")
    cli_mode_raw = _cli_val(argv, "--retainmode")
    if cli_has_mode and not cli_mode_raw:
        print("ERROR: --retainmode requires a value", file=sys.stderr)
        sys.exit(2)

    if cli_noretain and (cli_retain or cli_mode_raw):
        print(
            "ERROR: --noretain cannot be combined with --retain or --retainmode",
            file=sys.stderr,
        )
        sys.exit(2)

    if cli_mode_raw:
        mode = _normalize_artifact_retention_mode(cli_mode_raw)
        if mode == "off":
            print(
                "ERROR: --retainmode must be rerun_safe or minimal_disk",
                file=sys.stderr,
            )
            sys.exit(2)
        return mode, "CLI(--retainmode)"

    if cli_noretain:
        return "off", "CLI(--noretain)"

    if cli_retain:
        return "rerun_safe", "CLI(--retain)"

    env_raw = os.environ.get("ARTIFACT_RETENTION")
    if env_raw is not None and str(env_raw).strip():
        return _normalize_artifact_retention_mode(env_raw), "ENV(ARTIFACT_RETENTION)"

    cfg_raw = cfg.get("ARTIFACT_RETENTION")
    if cfg_raw is not None and str(cfg_raw).strip():
        return _normalize_artifact_retention_mode(cfg_raw), "CFG(ARTIFACT_RETENTION)"

    return "rerun_safe", "DEFAULT"


def _apply_bench_overrides(cfg: ConfigDict) -> None:
    """Force benchmark configuration settings."""
    cfg["TEST_MODE_ENABLE"] = "dud"
    
    # Enable all engines
    cfg["USE_GNINA"] = True
    cfg["USE_LEDOCK"] = True
    cfg["USE_DOCK6"] = True
    cfg["USE_SCORCH"] = True
    
    # Force Holo + pH Ensemble
    cfg["APO_HOLO_MODE"] = "holo"
    cfg["PH_ENSEMBLE"] = True
    
    # Restrict proteins
    cfg["SPECIFIED_PROTEINS"] = ",".join(BENCH_PDB_IDS)
    
    # Merge library map
    current_map = _coerce_test_map(cfg.get("TEST_LIBRARY_MAP", {}))
    # _coerce_test_map returns a dict-like object or dict. Ensure it's a dict.
    if hasattr(current_map, "copy"):
        new_map = dict(current_map)
    else:
        new_map = {}
        
    for k, v in BENCH_TEST_LIBRARY_MAP.items():
        new_map[k] = v
        
    cfg["TEST_LIBRARY_MAP"] = new_map


def _apply_bench2_overrides(cfg: ConfigDict) -> None:
    """Force benchmark2 configuration settings (Vina + SCORCH only)."""
    _apply_bench_overrides(cfg)
    cfg["USE_GNINA"] = False
    cfg["USE_LEDOCK"] = False
    cfg["USE_DOCK6"] = False
    cfg["USE_SCORCH"] = True


def _apply_dude_overrides(cfg: ConfigDict) -> None:
    """Force DUD-only runtime mode without touching on-disk config."""
    cfg["TEST_MODE_ENABLE"] = "dud"


def select_pdb_files_for_run(
    cfg: ConfigDict,
    argv_for_parsing: list[str],
    *,
    is_resume: bool = False,
    resume_protein_ids: list[str] | None = None,
) -> list[str]:
    cfg.setdefault("SPECIFIED_PROTEINS", "")
    requested_ids, _ = _parse_specified_proteins(argv_for_parsing, cfg)
    cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] = requested_ids
    print(
        f"[config] SPECIFIED_PROTEINS effective={requested_ids} (precedence: CLI>ENV>CFG)"
    )

    pdb_files = [
        f
        for f in os.listdir(cfg["INPUT_DIR"])
        if f.lower().endswith(".pdb") and "_nolig" not in f.lower()
    ]

    id_index: dict[str, str] = {}
    for f in pdb_files:
        base = os.path.splitext(f)[0].replace("_cleaned", "")
        nid = _norm_pdb_id(base)
        if nid:
            id_index.setdefault(nid, f)

    req = list(cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS", []) or [])
    if is_resume and not req and resume_protein_ids:
        req = list(resume_protein_ids)
    if req:
        hits = [nid for nid in req if nid in id_index]
        miss = [nid for nid in req if nid not in id_index]

        print(
            f"[filter.proteins] mode=on requested={len(req)} present={len(hits)} missing={len(miss)} ? {hits}"
        )
        for m in miss:
            print(
                f"WARNING: requested PDB '{m}' not found under INPUT_DIR={cfg['INPUT_DIR']} or was excluded (_nolig)."
            )

        if not hits:
            print(
                "ERROR: No requested proteins found. Exiting with status 2 to avoid a no-op run."
            )
            sys.exit(2)

        pdb_files = [id_index[nid] for nid in hits]
        print("Selected proteins (Specified Proteins Mode): " + ", ".join(hits))
    else:
        print(
            f"[filter.proteins] mode=off requested=0 present={len(pdb_files)} missing=0 ? []"
        )

    tokens = parse_test_libraries(cfg)
    raw_map = cfg.get("TEST_LIBRARY_MAP", {})
    test_map = _coerce_test_map(raw_map)
    try:
        cfg["_TEST_LIBRARY_CANONICAL"] = {
            str(k).upper(): str(v) for k, v in getattr(test_map, "items", lambda: [])()
        }
    except Exception:
        cfg["_TEST_LIBRARY_CANONICAL"] = {}

    if "dud" in tokens:
        test_keys = set()
        for k in getattr(test_map, "keys", lambda: [])():
            nid = _norm_pdb_id(str(k))
            if nid:
                test_keys.add(nid)

        if test_keys:
            specified_raw = cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS") or []
            specified_keys = set()
            for s in specified_raw:
                nid = _norm_pdb_id(str(s))
                if nid:
                    specified_keys.add(nid)

            kept, skipped = [], []
            for f in pdb_files:
                nid = _norm_pdb_id(f)
                if nid and (nid in test_keys or nid in specified_keys):
                    kept.append(f)
                else:
                    skipped.append(f)

            if skipped:
                print(
                    f"[test-mode] Enabled tokens={'+'.join(tokens)}; restricting to "
                    f"{len(kept)} PDBs from TEST_LIBRARY_MAP keys"
                    + (" (plus specified proteins)." if specified_keys else ".")
                )
                for s in skipped:
                    print(f"[test-mode] Skipping {s} (not in TEST_LIBRARY_MAP).")

            pdb_files = kept
        else:
            print(
                "[test-mode] TEST_LIBRARY_MAP empty/invalid; no extra filtering applied."
            )

    return pdb_files


def _log_cfg_emit_path_check(pdb_id, receptor_path, variant, legacy):
    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    variant_label = variant_token or "legacy"
    contains_variant = bool(variant_token and variant_token in str(receptor_path))
    logging.info(
        "[cfg.emit.check] pdb=%s variant=%s receptor_path=%s path_contains_variant=%s",
        pdb_id,
        variant_label,
        receptor_path,
        contains_variant,
    )
    if variant_token and not contains_variant and not legacy:
        logging.warning(
            "[variant.mismatch] expected_variant=%s wrote_legacy_path=%s action=fail_ci",
            variant_token,
            receptor_path,
        )


def emit_vina_config(
    cfg,
    pdb_id,
    receptor_pdbqt,
    center,
    box_size,
    ligand_path,
    stage_name,
    stage_info,
    cpu_per_job,
    logger=None,
    *,
    variant=None,
    ph_token=None,
    legacy=False,
    skip_manifest_if_exists=False,
):
    if not skip_manifest_if_exists:
        result = _emit_vina_config_impl(
            cfg,
            pdb_id,
            receptor_pdbqt,
            center,
            box_size,
            ligand_path,
            stage_name,
            stage_info,
            cpu_per_job,
            logger,
            variant=variant,
            ph_token=ph_token,
            legacy=legacy,
        )
    else:
        result = _emit_vina_config_skip_manifest(
            cfg,
            pdb_id,
            receptor_pdbqt,
            center,
            box_size,
            ligand_path,
            stage_name,
            stage_info,
            cpu_per_job,
            logger,
            variant=variant,
            ph_token=ph_token,
            legacy=legacy,
        )
    _log_cfg_emit_path_check(pdb_id, receptor_pdbqt, variant, legacy)
    return result


def _emit_vina_config_skip_manifest(
    cfg,
    pdb_id,
    receptor_pdbqt,
    center,
    box_size,
    ligand_path,
    stage_name,
    stage_info,
    cpu_per_job,
    logger,
    *,
    variant=None,
    ph_token=None,
    legacy=False,
):
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    ph_label = (str(ph_token).strip() or None) if ph_token is not None else None
    legacy_mode = bool(legacy)

    lig_base = Path(ligand_path).stem
    run_id = cfg["RUN_ID"]

    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    cfg_dir.mkdir(parents=True, exist_ok=True)

    cfg_path = cfg_dir / f"{lig_base}_{stage_name}.txt"

    stage_root = docked_dir(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    stage_root.mkdir(parents=True, exist_ok=True)
    out_dir = stage_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{lig_base}_{stage_name}.pdbqt"

    expected_receptor = receptor_file(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    receptor_exists = expected_receptor.exists()
    receptor_for_config = str(expected_receptor)

    variant_display = variant_token or "None"
    ph_display = ph_label or "None"
    breadcrumb = (
        "[cfg.emit] run=%s pdb=%s stage=%s variant=%s ph=%s\n"
        "           cfg_dir=%s receptor=%s out_root=%s"
    )
    breadcrumb_args = (
        run_id,
        pdb_id,
        stage_name,
        variant_display,
        ph_display,
        str(cfg_dir),
        receptor_for_config,
        str(stage_root),
    )
    if logger:
        logger.info(breadcrumb, *breadcrumb_args)
    else:
        print(breadcrumb % breadcrumb_args)

    if not receptor_exists:
        msg = (
            f"[router.error] missing receptor for pdb={pdb_id} variant={variant_display} "
            f"ph={ph_display} -> {expected_receptor}"
        )
        if logger:
            logger.error(msg)
        else:
            print(msg)

    lines = [
        f"receptor = {receptor_for_config}",
        f"ligand   = {ligand_path}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {box_size[0]:.3f}",
        f"size_y   = {box_size[1]:.3f}",
        f"size_z   = {box_size[2]:.3f}",
        f"cpu      = {int(cpu_per_job)}",
        f"exhaustiveness = {int(stage_info.get('exhaustiveness', 8))}",
        f"energy_range   = {int(stage_info.get('energy_range', 4))}",
        f"num_modes      = {int(stage_info.get('num_modes', 4))}",
        f"verbosity      = {int(stage_info.get('verbosity', 0))}",
        f"out = {out_path}",
    ]

    if logger:
        cx, cy, cz = center
        sx, sy, sz = box_size
        lig_name = os.path.basename(str(ligand_path))
        logger.info(
            "[vina.cfg] lig=%s center=(%.3f,%.3f,%.3f) size=(%.1f,%.1f,%.1f)",
            lig_name,
            cx,
            cy,
            cz,
            sx,
            sy,
            sz,
        )

    payload = ("\n".join(lines)).encode("utf-8")
    overwrite = cfg_path.exists()

    tmp = cfg_path.with_suffix(".part")
    with open(tmp, "wb") as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, cfg_path)

    emit_msg = (
        "[cfg.emit] run=%s pdb=%s variant=%s ph=%s stage=%s ligand=%s "
        "cfg_dir=%s docked_root=%s path=%s overwrite=%s bytes=%d"
    )
    emit_args = (
        run_id,
        pdb_id,
        variant_display,
        ph_display,
        stage_name,
        lig_base,
        str(cfg_dir),
        str(stage_root),
        str(cfg_path),
        str(overwrite).lower(),
        len(payload),
    )
    if logger:
        logger.info(emit_msg, *emit_args)
    else:
        print(emit_msg % emit_args)

    return str(cfg_path), str(out_path)


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


def get_recenter_params(cfg: Dict) -> RecenterParams:
    """Load recenter parameters from config with safe defaults."""
    return RecenterParams(
        EARLY_RECENTER_RATIO=float(cfg.get("EARLY_RECENTER_RATIO", 0.70)),
        EARLY_RECENTER_MIN_EVAL=int(cfg.get("EARLY_RECENTER_MIN_EVAL", 10)),
        EARLY_RECENTER_FAR_A=float(cfg.get("EARLY_RECENTER_FAR_A", 15.0)),
        EARLY_RECENTER_MEDIAN_A=float(cfg.get("EARLY_RECENTER_MEDIAN_A", 10.0)),
        ALLOW_BOX_EXPAND=bool(cfg.get("ALLOW_BOX_EXPAND", True)),
        MAX_RECENTER_ATTEMPTS=int(cfg.get("MAX_RECENTER_ATTEMPTS", 1)),
    )


# >>> MAKE_PATHS SHIM START
# make_paths is imported from path_router above (legacy helper removed).
# >>> MAKE_PATHS SHIM END


def _smoke_emit_config_demo() -> None:
    """Emit a small config to exercise router paths in isolation."""
    smoke_log = logging.getLogger("smoke")
    old_variant = os.environ.get("APO_HOLO_VARIANT")
    old_run_env = os.environ.get("ATLAS_RUN_ID")
    try:
        base_cfg = ConfigDict(load_inputs())
    except Exception as exc:
        smoke_log.warning("[smoke.emit.skip] reason=%s", exc)
        return

    try:
        cfg = ConfigDict(base_cfg.copy())
        repo_root = Path(__file__).resolve().parent
        smoke_root = repo_root / "analysis" / "_smoke"
        overrides = {
            "OVERALL_DIR": smoke_root,
            "INPUT_DIR": smoke_root / "input_pdbs",
            "OUTPUT_DIR": smoke_root / "processed_pdbs",
            "DOCKED_DIR": smoke_root / "docked",
            "PREPPED_LIGANDS_DIR": smoke_root / "prepped_ligands",
            "LIGANDS_MOL2_DIR": smoke_root / "ligands_mol2",
            "CONFIGS_DIR": smoke_root / "configs",
        }
        for key, path_value in overrides.items():
            cfg[key] = str(path_value)
            Path(path_value).mkdir(parents=True, exist_ok=True)

        cfg["RUN_ID"] = "smoke_demo"
        cfg["RESET_CONFIGS"] = False
        init_config_run_dir(cfg, run_id=cfg["RUN_ID"], reset=False, logger=smoke_log)

        mode, variants = resolve_apo_holo_mode(cfg)
        smoke_log.info("[smoke.emit] mode=%s variants=%s", mode, variants)

        pdb_id = "3CS9"
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        lig_dir = paths.prepped_ligands_dir
        lig_dir.mkdir(parents=True, exist_ok=True)
        lig_path = lig_dir / "smoke_ligand.pdbqt"
        if not lig_path.exists():
            lig_path.write_text("SMOKE", encoding="utf-8")

        ph_label = "pH6_7"
        os.environ["APO_HOLO_VARIANT"] = "APO"

        receptor_path = receptor_file(
            paths.pdb_id, variant="APO", ph_tag=ph_label, legacy=False
        )
        receptor_path.parent.mkdir(parents=True, exist_ok=True)
        if not receptor_path.exists():
            receptor_path.write_text("RECEPTOR", encoding="utf-8")

        stage_info = {
            "name": "smoke_stage",
            "exhaustiveness": 8,
            "num_modes": 9,
            "verbosity": 0,
        }
        conf_path, out_path = emit_vina_config(
            cfg,
            paths.pdb_id,
            str(receptor_path),
            (0.0, 0.0, 0.0),
            (20.0, 20.0, 20.0),
            str(lig_path),
            stage_info["name"],
            stage_info,
            1,
            smoke_log,
            variant="APO",
            ph_token=ph_label,
            legacy=False,
        )
        smoke_log.info("[smoke.emit.done] config=%s out=%s", conf_path, out_path)
    except Exception as exc:
        smoke_log.warning("[smoke.emit.skip] reason=%s", exc)
    finally:
        if old_run_env is None:
            os.environ.pop("ATLAS_RUN_ID", None)
        else:
            os.environ["ATLAS_RUN_ID"] = old_run_env
        if old_variant is None:
            os.environ.pop("APO_HOLO_VARIANT", None)
        else:
            os.environ["APO_HOLO_VARIANT"] = old_variant


# ======================
# Program entry point
# ======================
def main() -> None:
    if _cli_has(sys.argv, "-h") or _cli_has(sys.argv, "--help"):
        print(HELP_TEXT)
        sys.exit(0)

    print("MODELLER is working with license.")
    is_resume = _cli_has(sys.argv, "-resume") or _cli_has(sys.argv, "--resume")
    cli_run_id = _cli_val(sys.argv, "--run-id") or _cli_val(sys.argv, "-run-id")
    is_bench = _bench_enabled(sys.argv)
    is_bench2 = _bench2_enabled(sys.argv)
    is_dude = _dude_enabled(sys.argv)

    if is_bench and is_bench2:
        print("ERROR: -bench and -bench2 cannot be combined", file=sys.stderr)
        sys.exit(2)

    if is_resume and is_bench:
        print("ERROR: -bench cannot be used with -resume", file=sys.stderr)
        sys.exit(2)

    if is_resume and is_bench2:
        print("ERROR: -bench2 cannot be used with -resume", file=sys.stderr)
        sys.exit(2)

    if is_dude and (is_bench or is_bench2):
        print(
            "ERROR: -dude cannot be used with -bench or -bench2",
            file=sys.stderr,
        )
        sys.exit(2)

    if is_dude and is_resume:
        print("ERROR: -dude cannot be used with -resume", file=sys.stderr)
        sys.exit(2)

    if is_resume and not cli_run_id:
        print("ERROR: -resume requires --run-id <RUN_ID>", file=sys.stderr)
        sys.exit(2)

    # Force PH_ENSEMBLE env var early if bench mode to bypass env override logic later
    if is_bench or is_bench2:
        os.environ["PH_ENSEMBLE"] = "1"

    run_id = _resolve_run_id(sys.argv)
    argv_for_parsing = list(sys.argv)
    resume_manifest = None
    resume_protein_ids: list[str] = []
    completed_lookup: dict[tuple[str, str], bool] = {}
    os.environ["ATLAS_RUN_ID"] = run_id
    log_path = _prepare_run_logfile(run_id)
    os.environ["ATLAS_LOG_FILE"] = log_path
    _tee_stdio_to(log_path)
    cfg = ConfigDict(load_inputs())

    if is_dude:
        _apply_dude_overrides(cfg)
        logging.info("[dude] DUD-only runtime mode ENABLED. Overrides applied.")

    if is_bench:
        _apply_bench_overrides(cfg)
        logging.info("[bench] benchmark mode ENABLED. Overrides applied.")
    elif is_bench2:
        _apply_bench2_overrides(cfg)
        logging.info("[bench2] benchmark2 mode ENABLED. Overrides applied (vina + scorch).")

    bootstrap_root_logging(cfg, log_path)
    logging.info("[probe.root] root-logger INFO now visible")
    print(f"[run] log_file={log_path} run_id={run_id}")

    if is_resume:
        snap_dict = _apply_resume_config_from_snapshot(cfg, run_id)
        cfg = ConfigDict(snap_dict)
        cfg["RUN_ID"] = run_id
    else:
        cfg["RUN_ID"] = run_id

    validate_config(cfg)

    rules = get_atom_rules()
    alias_sets = getattr(rules, "alias_sets", None)
    waters_set = set(getattr(rules, "waters", set()))
    if not waters_set and alias_sets is not None:
        waters_set = {
            str(tok).strip().upper()
            for tok in getattr(alias_sets, "waters", set())
            if str(tok).strip()
        }
    cofactors_set = {
        str(tok).strip().upper()
        for tok in getattr(rules, "cofactors", set())
        if str(tok).strip()
    }
    elements_set = {
        str(tok).strip().upper()
        for tok in getattr(rules, "elem_tokens_canonical", set())
        if str(tok).strip()
    }
    logging.info(
        "[aliases.summary] mode=%s keep={waters:%d, cofactors:%d, elements:%d}",
        getattr(rules, "policy_mode", "LEGACY"),
        len(waters_set),
        len(cofactors_set),
        len(elements_set),
    )

    cfg.setdefault("PH_ENSEMBLE", False)
    cfg.setdefault("PH_RADIUS", 1000000.0)

    env_ph_flag = os.environ.get("PH_ENSEMBLE")
    if env_ph_flag is not None:
        cfg.PH_ENSEMBLE = _to_bool(env_ph_flag)
    else:
        cfg.PH_ENSEMBLE = _to_bool(cfg.PH_ENSEMBLE)

    scope_env = os.environ.get("PH_SCOPE")
    if scope_env is not None:
        cfg.PH_SCOPE = scope_env.strip()
    elif "PH_SCOPE" in cfg:
        cfg.PH_SCOPE = cfg["PH_SCOPE"]

    radius_env = os.environ.get("PH_RADIUS")
    if radius_env is not None:
        try:
            cfg.PH_RADIUS = float(radius_env)
        except Exception:
            cfg.PH_RADIUS = 1000000.0
    else:
        try:
            cfg.PH_RADIUS = float(cfg.PH_RADIUS)
        except Exception:
            cfg.PH_RADIUS = 100000000.0
    if cfg.PH_RADIUS <= 0:
        cfg.PH_RADIUS = 10000000.0

    log = logging.getLogger("ph_ensemble")
    log.info(
        "[ph_ensemble.mode] enabled=%s scope=%s radius=%s",
        cfg.PH_ENSEMBLE,
        getattr(cfg, "PH_SCOPE", "auto"),
        getattr(cfg, "PH_RADIUS", 1000000.0),
    )

    if is_resume:
        resume_manifest = load_run_manifest(cfg, run_id)
        if resume_manifest is None:
            print(
                f"ERROR: -resume requested but run_manifest.yaml for run_id={run_id} "
                f"could not be loaded.",
                file=sys.stderr,
            )
            sys.exit(2)

        command_section = resume_manifest.get("command") or {}
        stored_argv = command_section.get("argv")
        if stored_argv:
            try:
                parsed_args = shlex.split(str(stored_argv))
                if parsed_args:
                    argv_for_parsing = parsed_args
                    logging.info("[resume] using stored argv=%s", parsed_args)
            except Exception:
                logging.warning(
                    "[resume] Failed to parse stored argv=%r",
                    stored_argv,
                    exc_info=True,
                )

        listed = command_section.get("pdb_list") or []
        if isinstance(listed, (list, tuple)):
            resume_protein_ids.extend(
                str(x).strip().upper() for x in listed if str(x).strip()
            )

        if not resume_protein_ids:
            summary_section = resume_manifest.get("summary") or {}
            summary_list = summary_section.get("total_protein_list") or []
            if isinstance(summary_list, (list, tuple)):
                resume_protein_ids.extend(
                    str(x).strip().upper() for x in summary_list if str(x).strip()
                )

        proteins = resume_manifest.get("proteins") or {}
        aggregated: dict[tuple[str, str], bool] = {}
        if isinstance(proteins, Mapping):
            for key, entry in proteins.items():
                try:
                    raw_key = str(key)
                    parts = raw_key.split("|", 2)
                    if len(parts) != 3:
                        continue
                    pdb_part, variant_part, _ph = parts
                    pdb_token = (pdb_part or "").strip().upper()
                    if not pdb_token:
                        continue
                    variant_token = (
                        variant_part or "legacy"
                    ).strip().lower() or "legacy"
                    ckey = (pdb_token, variant_token)
                    if not resume_protein_ids:
                        resume_protein_ids.append(pdb_token)
                    aggregated.setdefault(ckey, True)
                    status = (entry.get("status") or "").strip().lower()
                    if status != "completed":
                        aggregated[ckey] = False
                except Exception:
                    logging.warning(
                        "[resume.lookup.skip] key=%r entry=%r",
                        key,
                        entry,
                        exc_info=True,
                    )

        completed_lookup = {k: True for k, v in aggregated.items() if v}
        if resume_protein_ids:
            seen_ids = set()
            deduped = []
            for pid in resume_protein_ids:
                token = pid.strip().upper()
                if not token or token in seen_ids:
                    continue
                seen_ids.add(token)
                deduped.append(token)
            resume_protein_ids = deduped
        logging.info(
            "[resume] loaded manifest for run_id=%s; completed protein entries=%d",
            run_id,
            len(completed_lookup),
        )

    # --- PH-ligand mode (CLI > ENV > CFG) ---
    cfg.setdefault("PH_LIGAND_MODE", "off")
    cli_ph_mode = _cli_val(argv_for_parsing, "--ph-ligand-mode")
    env_ph_mode = os.environ.get("PH_LIGAND_MODE")
    cfg_ph_mode = str(cfg.get("PH_LIGAND_MODE", "off"))

    effective_ph_mode = next(
        (
            m
            for m in (cli_ph_mode, env_ph_mode, cfg_ph_mode)
            if m is not None and str(m).strip() != ""
        ),
        "off",
    )
    cfg["PH_LIGAND_MODE"] = effective_ph_mode

    log.info(
        "[ph_ligand.mode] effective=%r (source=%s)",
        effective_ph_mode,
        "CLI" if cli_ph_mode else "ENV" if env_ph_mode else "CFG",
    )

    # --- Per-run configs (RUN_DIR) ---
    cfg.setdefault("CONFIGS_DIR", str(Path(cfg["OVERALL_DIR"]) / "configs"))
    cfg.setdefault("RESET_CONFIGS", True)

    # CLI > ENV > CFG
    cli_cfg_dir = _cli_val(argv_for_parsing, "--configs-dir")
    cli_no_reset = _cli_has(argv_for_parsing, "--no-reset-configs")

    if cli_cfg_dir:
        cfg["CONFIGS_DIR"] = cli_cfg_dir
    cfg["RUN_ID"] = run_id
    if is_resume:
        cfg["RESET_CONFIGS"] = False
    elif cli_no_reset:
        cfg["RESET_CONFIGS"] = False

    init_config_run_dir(
        cfg,
        run_id=cfg.get("RUN_ID"),
        reset=cfg.get("RESET_CONFIGS"),
        logger=logging.getLogger("run"),
    )
    run_log_dir = run_logs_dir(cfg)
    pipeline_log_path = run_log_dir / "pipeline.log"
    ensure_file_handler(logging.getLogger(), pipeline_log_path, cfg, level=logging.INFO)
    logging.info("[run.log] pipeline_log=%s", pipeline_log_path)

    if is_bench and "CONFIG_RUN_DIR" in cfg:
        try:
            bench_snap_path = Path(cfg["CONFIG_RUN_DIR"]) / "bench_config.txt"
            lines = [
                "# Benchmark Configuration Snapshot",
                f"TEST_MODE_ENABLE={cfg.get('TEST_MODE_ENABLE')}",
                "USE_GNINA=true",
                "USE_LEDOCK=true",
                "USE_DOCK6=true",
                "USE_SCORCH=true",
                "APO_HOLO_MODE=holo",
                "PH_ENSEMBLE=true",
                f"BENCH_PDBS={cfg.get('SPECIFIED_PROTEINS')}",
            ]
            for k, v in BENCH_TEST_LIBRARY_MAP.items():
                lines.append(f"TEST_LIBRARY_MAP.{k}={v}")
            
            bench_snap_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            logging.info("[bench] wrote snapshot config to %s", bench_snap_path)
        except Exception:
            logging.warning("[bench] failed to write snapshot config", exc_info=True)

    if is_bench2 and "CONFIG_RUN_DIR" in cfg:
        try:
            bench2_snap_path = Path(cfg["CONFIG_RUN_DIR"]) / "bench2_config.txt"
            lines = [
                "# Benchmark2 Configuration Snapshot",
                f"TEST_MODE_ENABLE={cfg.get('TEST_MODE_ENABLE')}",
                "USE_GNINA=false",
                "USE_LEDOCK=false",
                "USE_DOCK6=false",
                "USE_SCORCH=true",
                "APO_HOLO_MODE=holo",
                "PH_ENSEMBLE=true",
                f"BENCH_PDBS={cfg.get('SPECIFIED_PROTEINS')}",
            ]
            for k, v in BENCH_TEST_LIBRARY_MAP.items():
                lines.append(f"TEST_LIBRARY_MAP.{k}={v}")

            bench2_snap_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            logging.info("[bench2] wrote snapshot config to %s", bench2_snap_path)
        except Exception:
            logging.warning("[bench2] failed to write snapshot config", exc_info=True)

    print(f"[cfg.run] run_id={cfg['RUN_ID']} run_dir={cfg['CONFIG_RUN_DIR']}")
    print(f"[ph.mode] PH_ENSEMBLE={cfg.get('PH_ENSEMBLE', False)}")

    global_start = time.time()

    # Best-effort manifest initialization
    if not is_resume:
        try:
            init_run_manifest(cfg, run_id, sys.argv, log_path)
        except Exception:
            logging.warning("Failed to initialize run_manifest.yaml", exc_info=True)
        else:
            try:
                update_manifest_for_config_hash(cfg, run_id)
            except Exception:
                logging.warning(
                    "[run-manifest] Failed to record config hash run_id=%s",
                    run_id,
                    exc_info=True,
                )
    else:
        logging.info("[resume] Skipping manifest initialization for run_id=%s", run_id)

    # --- Single-ligand config (ported) ---------------------------------------
    cfg.setdefault("SINGLE_LIGAND", "")
    cfg.setdefault("SINGLE_LIGAND_SUGGESTIONS", 5)
    cfg.setdefault("ALLOW_FDA_FALLBACK", False)
    cfg.setdefault("LIBRARY_MANIFEST_FILENAME", "_manifest.json")
    # Default + env override for building manifests during fallback scans
    cfg.setdefault("LIBRARY_MANIFEST_BUILD_ON_SCAN", True)
    env_build_flag = os.environ.get("LIBRARY_MANIFEST_BUILD_ON_SCAN")
    if env_build_flag is not None:
        try:
            cfg["LIBRARY_MANIFEST_BUILD_ON_SCAN"] = _to_bool(env_build_flag)
        except Exception:
            # If parsing fails, keep the config/default
            pass

    logging.getLogger("lib-manifest").info(
        "[lib-manifest.scan] build_on_scan=%s",
        str(bool(cfg.get("LIBRARY_MANIFEST_BUILD_ON_SCAN", True))).lower(),
    )
    cfg.setdefault(
        "FDA_MAPPING_CSV", str(Path(__file__).with_name("fda_mapping_from_pdbqt.csv"))
    )

    # CLI > ENV > CFG precedence
    cli_single = _parse_single_from_cli(argv_for_parsing)
    env_single = os.environ.get("SINGLE_LIGAND", "").strip()
    cfg_single = str(cfg.get("SINGLE_LIGAND", "")).strip()

    effective_single = next((x for x in (cli_single, env_single, cfg_single) if x), "")
    cfg["_EFFECTIVE_SINGLE_LIGAND"] = effective_single
    if effective_single:
        print(
            f"[config] SINGLE_LIGAND effective='{effective_single}' "
            f"(order=CLI>{'ENV' if env_single else ''}>{'CFG' if cfg_single else ''})"
        )

    # --- Library subfolder selection -----------------------------------
    cfg.setdefault("LIBRARY_SUBDIR_DEFAULT", "fda_library")
    cfg.setdefault("HMDB_LIBRARY_SUBDIR", "hmdb")
    cfg.setdefault("TEST_MODE_ENABLE", "off")
    # Accept dict or JSON-ish string
    if "TEST_LIBRARY_MAP" not in cfg:
        cfg["TEST_LIBRARY_MAP"] = {}

    # DeepCoy defaults (decoy autogen remains off unless toggled)
    apply_deepcoy_defaults(cfg)

    # Small FDA test library override.
    # When called with: python main.py -test -test-fda -fast
    # we want to use a tiny FDA test library instead of the full FDA set.
    #
    # By default this points at "fda_test_library_10", which should correspond to:
    #   prepped_ligands/fda_test_library_10/
    #   extracted_ligands/fda_test_library_10/
    cfg.setdefault("TEST_FDA_LIBRARY_SUBDIR", "fda_test_library_10")

    if _cli_has(argv_for_parsing, "-test-fda") or _cli_has(
        argv_for_parsing, "--test-fda"
    ):
        cfg["LIBRARY_SUBDIR_DEFAULT"] = cfg.get(
            "TEST_FDA_LIBRARY_SUBDIR",
            "fda_test_library_10",
        )
        print(
            f"[config] TEST_FDA_LIBRARY enabled: "
            f"LIBRARY_SUBDIR_DEFAULT={cfg['LIBRARY_SUBDIR_DEFAULT']}"
        )
    # --- Fast mode: force exhaustiveness=1 everywhere ---
    cfg["FAST_MODE"] = _parse_fast_flag(argv_for_parsing) or bool(
        cfg.get("FAST_MODE", False)
    )
    apply_deepcoy_cli_overrides(cfg, argv_for_parsing)
    if cfg["FAST_MODE"]:
        print("[config] FAST_MODE effective=True (exhaustiveness=1)")

    force_calibrator = _cli_has(argv_for_parsing, "--force-calibrator") or _cli_has(
        argv_for_parsing, "-force-calibrator"
    )
    cfg_force_raw = cfg.get("FORCE_CALIBRATOR", False)
    try:
        cfg_force_bool = _to_bool(cfg_force_raw)
    except Exception:
        cfg_force_bool = bool(cfg_force_raw)
    cfg["FORCE_CALIBRATOR"] = bool(force_calibrator or cfg_force_bool)
    if cfg["FORCE_CALIBRATOR"]:
        logging.info("[calibrator.force] enabled=true")

    # --- Control consensus (multi-engine control docking) -------------
    cfg.setdefault("CONTROL_CONSENSUS", False)
    control_consensus = False
    try:
        control_consensus = _to_bool(cfg.get("CONTROL_CONSENSUS", False))
    except Exception:
        control_consensus = bool(cfg.get("CONTROL_CONSENSUS", False))

    if _cli_has(argv_for_parsing, "-control-consensus") or _cli_has(
        argv_for_parsing, "--control-consensus"
    ):
        control_consensus = True

    cfg["CONTROL_CONSENSUS"] = bool(control_consensus)
    if cfg["CONTROL_CONSENSUS"]:
        print(
            "[config] CONTROL_CONSENSUS effective=True (multi-engine control docking enabled)"
        )

    # --- No-library docking mode (controls-only + ligand planning) ------
    #
    # If enabled, we still run receptor prep + control redocking +
    # pocket detection + ligand selection/enumeration, but we skip the
    # actual docking calls for DUD/FDA libraries.
    #
    cfg.setdefault("NO_LIBRARY_DOCKING", False)

    cli_no_dock = _cli_has(argv_for_parsing, "--no-docking")
    env_no_dock = os.environ.get("NO_LIBRARY_DOCKING")

    effective_no_dock = False
    if cli_no_dock:
        effective_no_dock = True
    elif env_no_dock is not None:
        try:
            effective_no_dock = _to_bool(env_no_dock)
        except Exception:
            effective_no_dock = False
    else:
        try:
            effective_no_dock = bool(cfg.get("NO_LIBRARY_DOCKING", False))
        except Exception:
            effective_no_dock = False

    cfg["NO_LIBRARY_DOCKING"] = effective_no_dock
    if effective_no_dock:
        print(
            "[config] NO_LIBRARY_DOCKING=True "
            "(controls-only; skip DUD/FDA docking, but still enumerate ligands)"
        )

    # --- Artifact retention mode (CLI > ENV > CFG > default) ---
    retention_mode, retention_source = _resolve_artifact_retention_mode(cfg, sys.argv)
    cfg["ARTIFACT_RETENTION"] = retention_mode
    print(
        f"[config] ARTIFACT_RETENTION effective={retention_mode} source={retention_source}"
    )

    # --- Center selection knobs (safe defaults) ---
    cfg.setdefault(
        "CENTER_MODE", "control-first"
    )  # ["control-first","hybrid","library-first"]
    cfg.setdefault("CONTROL_BLACKLIST", "GOL,EDO,PG4,MPD,ACT,SO4,PO4,CL,NA,CA")
    cfg.setdefault("CONTROL_MIN_HEAVY_ATOMS", 10)
    cfg.setdefault(
        "CONTROL_ANCHOR_MIN_VALID_RATE", 0.10
    )  # if current cluster has control hits + =10% valid, anchor
    cfg.setdefault("ALLOW_SWITCH_FROM_CONTROL", True)
    cfg.setdefault("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False)
    cfg.setdefault(
        "SWITCH_AWAY_FROM_CONTROL_MIN_BOOST", 2.5
    )  # kcal/mol median boost needed to leave control
    # Optional lock score gate (kcal/mol). Use a large positive number (or remove) to lock on RMSD alone.
    cfg.setdefault("CONTROL_LOCK_SCORE_MAX", -6.0)
    cfg.setdefault("CONTROL_LOCK_MIN_HITS", 1)  # require = this many validated controls
    cfg.setdefault(
        "CONTROL_LOCK_CENTER_MAX_DIST", 4.0
    )  # A; control centroid must be within this of center

    # clustering + switching thresholds
    cfg.setdefault("CLUSTER_EPS_ANG", 3.5)
    cfg.setdefault("SWITCH_VALID_RATE_MIN", 0.40)
    cfg.setdefault("SWITCH_SCORE_IMPROVE_MIN", 1.5)
    cfg.setdefault("SWITCH_LE_GAIN_MIN", 0.02)
    cfg.setdefault("SWITCH_SCORE_THRESHOLD", 0.70)
    cfg.setdefault("SWITCH_SCORE_HYSTERESIS", 0.50)

    # Hard cap on global switching (early recenter, empty-stage fallback, selector promotions)
    cfg.setdefault("MAX_GLOBAL_CENTER_SWITCHES", 2)

    # Critical defaults
    cfg.setdefault("THREADS_PER_VINA", 1)
    cfg.setdefault("RMSD_FILTER_ANG", 2.0)
    cfg.setdefault("RMSD_MAX_MODELS", 3)

    # --- safe defaults ---
    cfg.setdefault("RETRY_NEAR_MISS", True)
    cfg.setdefault("RETRY_EXHAUST_MULT", 2)
    cfg.setdefault("RETRY_DIST_THRESH", 7.0)
    cfg.setdefault("RETRY_SCORE_THRESH", -7.5)

    cfg.setdefault("ADAPTIVE_SHRINK_ENABLE", True)
    cfg.setdefault("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0)
    cfg.setdefault("ADAPTIVE_SHRINK_DEC", 4.0)
    cfg.setdefault("ADAPTIVE_SHRINK_MIN_BOX", 14.0)

    cfg.setdefault("CHECKPOINT_ENABLE", True)
    cfg.setdefault("PAINS_NAMES", "")
    cfg.setdefault("RECEPTOR_SANITY_CHECK", True)
    cfg.setdefault("AUDIT_JSON", True)

    # --- logging/noise controls ---
    cfg.setdefault("QUIET_CONSOLE", False)  # console shows WARN+ only; file keeps DEBUG
    cfg.setdefault("VINA_VERBOSITY", 2)  # 0=minimal, 1=normal, 2=verbose
    cfg.setdefault(
        "FILTER_VINA_STDOUT", False
    )  # reserved if we need extra filtering later

    # RMSD PARAMETERS
    cfg.setdefault("SELF_RMSD_MAX_ANG", 2.0)
    cfg.setdefault(
        "SELF_RMSD_REQUIRE_FOR_CONTROLS", True
    )  # reserved for future stricter gating
    cfg.setdefault(
        "EARLY_EXIT_MAX_MODELS", 3
    )  # validate at most N poses, stop on first PASS
    cfg.setdefault(
        "MAX_RETRY_SECONDS_PER_LIGAND", 300
    )  # wall-clock for retries/validation per ligand

    params = get_recenter_params(cfg)
    stages = define_docking_stages(cfg.get("DOCKING_MODE", "discovery").lower())
    print("current docking mode is ", cfg.get("DOCKING_MODE"))

    pdb_files = select_pdb_files_for_run(
        cfg,
        argv_for_parsing,
        is_resume=is_resume,
        resume_protein_ids=resume_protein_ids,
    )
    tokens = parse_test_libraries(cfg)
    test_mode = _resolve_test_mode(cfg)

    print("Working directory:", os.getcwd())
    print("Loaded config keys:", list(cfg.keys()))
    print(f"Proteins queued: {len(pdb_files)}")
    if not is_resume:
        try:
            scheduled_ids: list[str] = []
            for f in pdb_files:
                nid = _norm_pdb_id(f)
                if nid:
                    scheduled_ids.append(nid.upper())
            update_manifest_for_scheduled_proteins(cfg, run_id, scheduled_ids)
        except Exception:
            print(
                "[run-manifest] WARNING: failed to record scheduled proteins in manifest",
                file=sys.stderr,
            )

    start = time.time()
    plan_only = os.environ.get("A2_PLAN_ONLY") == "1"
    mode, variants = resolve_apo_holo_mode(cfg)
    router_legacy = mode == "legacy"
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))
    cfg_raw_mode = cfg.get("APO_HOLO_MODE")
    logging.info(
        "[apo-holo] cfg_token_raw=%r resolved_mode=%s variants=%s",
        cfg_raw_mode,
        mode,
        variants,
    )
    logging.info(
        "[apo-holo] router_legacy=%s ph_enabled=%s",
        router_legacy,
        ph_enabled,
    )
    logging.info(
        "[apo-holo] normalized_mode_token=%s",
        _debug_normalize_mode_token(cfg_raw_mode),
    )
    cfg["_ROUTER_LEGACY"] = router_legacy
    cfg["_RESOLVED_APO_HOLO_MODE"] = mode
    logging.info(
        "[apo-holo] resolved mode=%s variants=%s cfg_token=%r",
        mode,
        variants,
        cfg_raw_mode,
    )
    try:
        update_manifest_for_run_config(cfg, run_id)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record run-level apo/holo + water policy",
            exc_info=True,
        )

    # Where to write per-PDB failure logs
    overall_dir = cfg.get("OVERALL_DIR", ".")
    failed_root = os.path.join(overall_dir, "failed")
    os.makedirs(failed_root, exist_ok=True)
    logging.info("[apo-holo] failed log directory: %s", failed_root)

    failed_entries = []  # (pdb_id, variant_label, log_path, exc_type, exc_msg)
    retention_lock = threading.Lock()

    for variant in variants:
        # Make variant visible to any module still reading env (legacy compatibility)
        os.environ["APO_HOLO_VARIANT"] = "" if variant is None else str(variant).upper()
        label = "legacy" if variant is None else str(variant).lower()
        logging.info(
            "[apo-holo] start_variant mode=%s variant_label=%s env_token=%s",
            mode,
            label,
            os.environ.get("APO_HOLO_VARIANT", ""),
        )

        with tqdm(
            total=len(pdb_files),
            desc=f"Processing Proteins ({label})",
            unit="protein",
            position=0,
            dynamic_ncols=True,
            mininterval=0.2,
            leave=True,
            file=sys.stdout,
        ) as bar:
            cfg_v = cfg  # no per-variant mutation; variant is propagated via APO_HOLO_VARIANT env

            single_ligand_mode = bool(cfg.get("_EFFECTIVE_SINGLE_LIGAND"))
            multi_pdb_single_ligand = single_ligand_mode and len(pdb_files) > 1
            cpu = int(cfg.get("CPU", os.cpu_count() or 1))
            max_pdb_workers = max(1, min(cpu, len(pdb_files))) if pdb_files else 1

            def _process_one(pdb_file: str, cfg_for_pdb):
                pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()
                pdb_start = time.time()

                if is_resume:
                    resume_key = (pdb_id, label)
                    if completed_lookup.get(resume_key):
                        logging.info(
                            "[resume.skip] pdb_id=%s variant=%s already completed; skipping.",
                            pdb_id,
                            label,
                        )
                        return

                try:
                    library_name = None
                    try:
                        if "dud" in tokens:
                            lib_map = (
                                cfg_for_pdb.get("_TEST_LIBRARY_CANONICAL", {}) or {}
                            )
                            library_name = lib_map.get(pdb_id)
                        if not library_name:
                            library_name = cfg_for_pdb.get("LIBRARY_SUBDIR_DEFAULT")
                    except Exception:
                        library_name = cfg_for_pdb.get("LIBRARY_SUBDIR_DEFAULT")

                    try:
                        update_manifest_for_protein_start(
                            cfg_for_pdb, run_id, pdb_id, label, library_name
                        )
                    except Exception:
                        logging.warning(
                            "Failed to update run_manifest for start of %s (%s)",
                            pdb_id,
                            label,
                            exc_info=True,
                        )

                    # Main per-PDB work
                    process_one_protein(cfg_for_pdb, pdb_file, stages, params)
                    # Run SCORCH per protein immediately after docking completion.
                    # The rescoring script handles per-protein idempotency using _DONE sentinels.
                    try:
                        _maybe_run_scorch_rescore_for_pdb(
                            cfg_for_pdb,
                            run_id,
                            pdb_id,
                            variant=None if variant is None else str(variant),
                            verbose=False,
                        )
                    except Exception:
                        logging.warning(
                            "[scorch-rescore.hook] action=skip reason=unexpected_exception pdb_id=%s variant=%s",
                            pdb_id,
                            label,
                            exc_info=True,
                        )

                    pdb_elapsed = time.time() - pdb_start
                    try:
                        update_manifest_for_protein_success(
                            cfg_for_pdb, run_id, pdb_id, label, pdb_elapsed
                        )
                    except Exception:
                        logging.warning(
                            "Failed to update run_manifest for success of %s (%s)",
                            pdb_id,
                            label,
                            exc_info=True,
                        )

                    try:
                        _maybe_run_artifact_retention_for_pdb(
                            cfg_for_pdb,
                            run_id,
                            pdb_id,
                            lock=retention_lock,
                        )
                    except Exception:
                        logging.warning(
                            "[artifact-retention.hook] action=skip reason=unexpected_exception pdb_id=%s variant=%s",
                            pdb_id,
                            label,
                            exc_info=True,
                        )

                except Exception as exc:
                    # Per-PDB failure handling
                    exc_type = type(exc).__name__
                    exc_msg = str(exc)
                    traceback_str = traceback.format_exc()

                    fail_log_path = (
                        Path(failed_root) / f"{run_id}__{pdb_id}__{label}.log"
                    )
                    tmp_fail_log_path = fail_log_path.with_suffix(
                        fail_log_path.suffix + ".part"
                    )

                    # Write a dedicated failure log for this PDB+variant
                    with tmp_fail_log_path.open("w", encoding="utf-8") as fh:
                        fh.write(
                            f"[FAILED PDB]\n"
                            f"  pdb_id        = {pdb_id}\n"
                            f"  variant_label = {label}\n"
                            f"  mode          = {mode}\n"
                            f"  pdb_file      = {pdb_file} run_id={run_id}\n"
                            f"  exception     = {exc_type}: {exc_msg}\n\n"
                            f"[TRACEBACK]\n"
                            f"{traceback_str}\n"
                        )
                        traceback.print_exc(file=fh)
                    os.replace(tmp_fail_log_path, fail_log_path)

                    logging.error(
                        "[apo-holo] FAILED pdb_id=%s variant_label=%s; "
                        "see failure log at %s",
                        pdb_id,
                        label,
                        fail_log_path,
                    )

                    failed_entries.append(
                        (pdb_id, label, str(fail_log_path), exc_type, exc_msg)
                    )

                    try:
                        update_manifest_for_protein_failure(
                            cfg_for_pdb, run_id, pdb_id, label, fail_log_path
                        )
                    except Exception:
                        logging.warning(
                            "Failed to update run_manifest for failure of %s (%s)",
                            pdb_id,
                            label,
                            exc_info=True,
                        )

            if multi_pdb_single_ligand:
                logging.info(
                    "[main.parallel.single_ligand] mode=%s n_pdb=%d max_workers=%d",
                    label.upper(),
                    len(pdb_files),
                    max_pdb_workers,
                )

                def _cfg_for_pdb() -> ConfigDict:
                    cfg_local = cfg_v.copy()
                    cfg_local["CPU"] = 1
                    return cfg_local

                with ThreadPoolExecutor(max_workers=max_pdb_workers) as pool:
                    future_map = {
                        pool.submit(_process_one, pdb_file, _cfg_for_pdb()): pdb_file
                        for pdb_file in pdb_files
                    }
                    for fut in as_completed(future_map):
                        pdb_file = future_map[fut]
                        try:
                            fut.result()
                        except Exception as exc:
                            logging.exception(
                                "[main.parallel.single_ligand.error] pdb=%s error=%s",
                                pdb_file,
                                exc,
                            )
                        finally:
                            bar.update(1)
            else:
                for pdb_file in pdb_files:
                    try:
                        _process_one(pdb_file, cfg_v)
                    finally:
                        # Always advance the progress bar, even if this PDB failed
                        bar.update(1)

            try:
                apply_pocket_detection_events(cfg, run_id)
            except Exception:
                logging.warning(
                    "[run-manifest.pocket_detection.events.apply] run_id=%s action=skip",
                    run_id,
                    exc_info=True,
                )

    elapsed_min = (time.time() - start) / 60.0
    print(f"\nAll proteins processed in {elapsed_min:.2f} minutes.")

    if failed_entries:
        print("\nThe following proteins failed. See per-PDB logs under:", failed_root)
        for pdb_id, label, log_path, exc_type, exc_msg in failed_entries:
            print(
                f"  - {pdb_id} ({label}): {exc_type} â€” {exc_msg}\n"
                f"      log: {log_path}"
            )
    else:
        print("\nNo proteins recorded as failed.")

    try:
        finalize_run_manifest(
            cfg, run_id, start_time=global_start, failed_entries=failed_entries
        )
    except Exception:
        logging.warning("Failed to finalize run_manifest.yaml", exc_info=True)

    if plan_only:
        sys.exit(0)

    try:
        failed_lookup = {(entry[0], entry[1]) for entry in failed_entries}
        for variant in variants:
            label = "legacy" if variant is None else str(variant).lower()
            mmgbsa_variant_token: str | None = (
                None if variant is None else str(variant).upper()
            )
            legacy_mode = variant is None
            for pdb_file in pdb_files:
                pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()
                if (pdb_id, label) in failed_lookup:
                    continue
                _maybe_run_mmgbsa_for_pdb(
                    cfg,
                    pdb_file,
                    pdb_id,
                    mmgbsa_variant_token,
                    run_id,
                    test_mode,
                    legacy_mode,
                )
    except Exception:
        if _to_bool(cfg.get("MMGBSA_STRICT", False)):
            raise
        logging.warning(
            "[mmgbsa.pipeline] action=skip reason=unexpected_exception", exc_info=True
        )

    try:
        _maybe_run_dud_eval(cfg, run_id, pdb_files)
    except Exception:
        logging.warning(
            "[dud-eval.invoke] action=skip reason=unexpected_exception", exc_info=True
        )

    try:
        _maybe_run_master_schema_export(cfg, run_id)
    except Exception:
        logging.warning(
            "[master-export.invoke] action=skip reason=unexpected_exception",
            exc_info=True,
        )

    try:
        _maybe_run_report_generation(cfg, run_id)
    except Exception:
        logging.warning(
            "[report.invoke] action=skip reason=unexpected_exception", exc_info=True
        )

    try:
        _log_rescore_verification(run_id)
    except Exception:
        logging.warning(
            "[post-check.invoke] action=skip reason=unexpected_exception", exc_info=True
        )


def _send_run_email(status: int, start_time: str, end_time: str) -> None:
    """
    Best-effort email notification using the system `mail` command.
    Must never raise, so it is safe to call from finally blocks.
    """
    try:
        # Hostname: use uname if available (Linux/Unix), else fall back.
        try:
            host = os.uname().nodename
        except AttributeError:
            host = "unknown-host"

        subject = f"Atlas run exited with status {status} on {host}"

        # Mirror your shell script body as closely as possible
        cmd_line = "python main.py " + " ".join(sys.argv[1:])
        body_lines = [
            f"Atlas run finished on host: {host}",
            f"Start time: {start_time}",
            f"End time:   {end_time}",
            f"Exit status: {status}",
            "",
            "Command:",
            cmd_line,
            "",
        ]
        body = "\n".join(body_lines)

        mail_bin = shutil.which("mail")
        if not mail_bin:
            logging.info("[notify] mail command unavailable; skipping email")
            return

        try:
            proc = subprocess.Popen(
                [mail_bin, "-s", subject, "mpg2352@utexas.edu"],
                stdin=subprocess.PIPE,
                text=True,
            )
        except Exception as exc:
            logging.warning("[notify] failed to spawn mail command err=%s", exc)
            return

        try:
            proc.communicate(body, timeout=30)
        except Exception as exc:
            logging.warning("[notify] mail command failed err=%s", exc)
            try:
                proc.kill()
            except Exception:
                pass

    except Exception:
        # Absolute last-resort guard: never let notification kill the process
        try:
            logging.exception(
                "[notify] unexpected error while building notification email"
            )
        except Exception:
            pass


def _run_with_email_notification() -> None:
    """
    Wrapper used ONLY when main.py is invoked as a script.

    - Calls _smoke_emit_config_demo() and main() in the same order as before.
    - Preserves all exit codes (SystemExit, unhandled exceptions).
    - Always attempts to send an email in a finally block.
    """
    # Import-time already brought in `time`, `os`, `sys`, `logging`, etc.
    start_time = time.strftime("%Y-%m-%d %H:%M:%S")
    status: int = 0

    try:
        _smoke_emit_config_demo()
        main()
        # If main returns normally, status 0
        status = 0
    except SystemExit as exc:
        # Preserve the original exit code from sys.exit()
        code = exc.code
        status = code if isinstance(code, int) else 1
        raise
    except BaseException:
        # KeyboardInterrupt and other errors â†’ non-zero status
        status = 1
        raise
    finally:
        end_time = time.strftime("%Y-%m-%d %H:%M:%S")
        try:
            _send_run_email(status=status, start_time=start_time, end_time=end_time)
        except Exception:
            # Never let notification interfere with the original exit behavior
            try:
                logging.exception("[notify] email wrapper raised unexpectedly")
            except Exception:
                pass


if __name__ == "__main__":
    _run_with_email_notification()
