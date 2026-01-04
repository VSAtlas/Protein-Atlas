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

Single-ligand mode:
  --single PATTERN
      Enable SINGLE_LIGAND mode and restrict docking to a single ligand
      whose name contains PATTERN. The search order is controlled by:
        SINGLE_LIGAND_SEARCH_ORDER
        SINGLE_LIGAND_ALLOW_PREFIX
        SINGLE_LIGAND_MANIFEST_ONLY
      as defined in the YAML config. When combined with multiple PDBs,
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
          <OVERALL_DIR>/failed/<PDB>.<variant>.log

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
from activesite import extract_and_remove_ligands, get_atom_rules
import sys, hashlib, re, logging, json, time, os, shutil, re, shlex, subprocess, math, random
import csv
import statistics
from dataclasses import dataclass, field
import datetime
import yaml
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, Any, Mapping, Sequence
import numpy as np
from tqdm import tqdm
import traceback
from logging_topics import (
    _tee_stdio_to,
    make_protein_logger,
    bootstrap_root_logging,
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
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
from input_and_export_functions import (
    load_inputs,
    validate_config,
    define_docking_stages,
    write_score_summary_to_csv,
    extract_best_score,
    _to_bool,
    init_config_run_dir,
)
from docking_vina import emit_vina_config as _emit_vina_config_impl
from docking import (
    norm,
    _fingerprint_stage,
    select_ligands_for_next,
    final_pose_validation_and_screenshots,
    RetryManager,
    run_one_stage,
    _coerce_test_map,
    _resolve_test_mode,
    process_one_protein,
)
from prep_for_mmgbsa import prep_mmgbsa_from_sdfs
import protein_prep_mmgbsa as ppm
from protein_prep_mmgbsa import (
    prep_mmgbsa_receptor,
    run_tleap,
    write_leap_for_ligands,
)
from mmgbsa_trajectory import make_mmgbsa_trajectory, run_implicit_md
from run_mmgbsa import run_mmgbsa, parse_mmpbsa_delta_total, write_aggregated_mmpbsa_results
from path_router import (
    expand_variants,
    make_paths,
    Paths as RouterPaths,
    receptor_file,
    docked_dir,
    config_dir as router_config_dir,
    ph_ensemble_dir as router_ph_ensemble_dir,
)
from fallback_recenter import RecenterParams
from apo_holo_mode import (
    resolve_apo_holo_mode,
    _debug_normalize_mode_token,
    _variant_receptor_path,
    file_sha1,
    delete_variant_trees,
    _record_apo_holo_usage,
    _record_apo_holo_decision,
)
import automate_protein_prep as protein_prep
from single_ligand_index import (
    _ensure_single_ligand_index,
    _resolve_single_ligand,
)
from debug_fs import install_debug_makedirs

# Install debug wrappers for Path.mkdir and os.makedirs at import time,
# preserving the previous behavior.
install_debug_makedirs()


def _resolve_run_id(argv: list[str]) -> str:
    cli_run_id = _cli_val(argv, "--run-id") or _cli_val(argv, "-run-id")
    env_run_id = (os.environ.get("ATLAS_RUN_ID") or "").strip()
    if cli_run_id:
        return cli_run_id
    if env_run_id and env_run_id.lower() not in {"smoke_demo", "main_smoke_demo"}:
        return env_run_id
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def _prepare_run_logfile(run_id: str) -> str:
    logs_dir = Path.cwd() / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    return str(logs_dir / f"main_{run_id}.log")




class ConfigDict(dict):
    __slots__ = ()
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc
    def __setattr__(self, key, value):
        self[key] = value
    def copy(self):
        return ConfigDict(super().copy())


def _apply_resume_config_from_snapshot(
    cfg: Mapping[str, Any],
    run_id: str,
) -> dict:
    """
    For resume mode:

    - Load run_manifest.yaml for run_id.
    - Find the per-run config snapshot (paths.run_dir + paths.config_file).
    - Load it as a dict.
    - Return that dict so the caller can use it as the canonical cfg.

    If anything goes wrong, this logs a WARNING and returns the original cfg as a plain dict.
    """
    logger = logging.getLogger("resume.config")
    try:
        manifest = load_run_manifest(cfg, run_id)
        if not manifest:
            logger.warning(
                "[resume.config] manifest missing or empty for run_id=%s; "
                "falling back to current cfg.",
                run_id,
            )
            return dict(cfg)

        paths = manifest.get("paths") or {}
        run_dir = paths.get("run_dir")
        cfg_file_name = paths.get("config_file", "run_config.yaml")

        if not run_dir:
            logger.warning(
                "[resume.config] run_dir missing in manifest paths for run_id=%s; "
                "falling back to current cfg.",
                run_id,
            )
            return dict(cfg)

        cfg_path = Path(run_dir) / cfg_file_name
        if not cfg_path.exists():
            logger.warning(
                "[resume.config] snapshot config file missing at %s for run_id=%s; "
                "falling back to current cfg.",
                cfg_path,
                run_id,
            )
            return dict(cfg)

        with cfg_path.open("r", encoding="utf-8") as fh:
            snap = yaml.safe_load(fh) or {}
        if not isinstance(snap, dict):
            logger.warning(
                "[resume.config] snapshot config not a mapping at %s; falling back.",
                cfg_path,
            )
            return dict(cfg)

        cmd = manifest.get("command") or {}
        manifest_hash = cmd.get("config_hash")
        if manifest_hash:
            try:
                from run_manifest import compute_config_hash  # avoid import cycles

                snap_hash = compute_config_hash(snap)
                if snap_hash != manifest_hash:
                    logger.warning(
                        "[resume.config] config_hash mismatch for run_id=%s: "
                        "manifest=%s snapshot=%s",
                        run_id,
                        manifest_hash,
                        snap_hash,
                    )
            except Exception:
                logger.warning(
                    "[resume.config] unable to compute config_hash for snapshot",
                    exc_info=True,
                )

        logger.info(
            "[resume.config] loaded snapshot config from %s for run_id=%s",
            cfg_path,
            run_id,
        )
        return snap
    except Exception:
        logger.warning(
            "[resume.config] unexpected error loading snapshot config; falling back",
            exc_info=True,
        )
        return dict(cfg)


def _log_cfg_emit_path_check(pdb_id, receptor_path, variant, legacy):
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
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

    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
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
def _parse_single_from_cli(argv) -> str:
    """
    Minimal CLI parser for: --single <pattern>
    Returns the pattern string or "" if not provided.
    """
    try:
        if "--single" in argv:
            i = argv.index("--single")
            if i + 1 < len(argv) and not argv[i+1].startswith("-"):
                return argv[i+1]
    except Exception:
        pass
    return ""
def _parse_fast_flag(argv) -> bool:
    """Return True if argv includes fast/-fast/--fast (case-insensitive)."""
    try:
        return any(tok.lower().lstrip("-") == "fast" for tok in argv)
    except Exception:
        return False


def _parse_force_deepcoy_flag(argv) -> bool:
    """
    Return True if argv includes -force-deepcoy/--force-deepcoy (case-insensitive).
    """
    try:
        return any(_normalize_flag_name(tok) == "force-deepcoy" for tok in argv if tok.startswith("-"))
    except Exception:
        return False
def _iter_pdbqt_dirfirst(root: Path, allowed_subdirs: Optional[set[str]] = None):
    """
    Yield .pdbqt files with a directory-first strategy:
      - list files directly under `root`
      - then list files under first-level subdirs, optionally restricted by `allowed_subdirs`
    Falls back to rglob if listing fails (robustness over speed).
    """
    try:
        if not root or not root.exists():
            return
        # files at root
        for p in root.glob("*.pdbqt"):
            yield p
        # first-level subdirs (dir-name filter first, then files)
        for d in root.iterdir():
            if not d.is_dir():
                continue
            if allowed_subdirs is not None and d.name not in allowed_subdirs:
                continue
            for p in d.glob("*.pdbqt"):
                yield p
    except Exception:
        # robust fallback
        for p in root.rglob("*.pdbqt"):
            yield p

def _normalize_flag_name(tok: str) -> str:
    """
    Normalize a CLI flag token by stripping leading dashes and
    lowercasing. This lets '-flag' and '--flag' be interchangeable.
    """
    return str(tok).lstrip("-").lower()

def _cli_val(argv, flag):
    """
    Return the value following a flag, treating '-flag' and '--flag'
    as equivalent. Example: _cli_val(sys.argv, "--run-id") will find
    values from either '-run-id' or '--run-id'.
    """
    try:
        target = _normalize_flag_name(flag)
        for i, tok in enumerate(argv):
            # Only consider tokens that look like flags
            if not tok.startswith("-"):
                continue
            if _normalize_flag_name(tok) == target:
                if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                    return argv[i + 1]
    except Exception:
        pass
    return None

def _cli_has(argv, flag):
    """
    Return True if a flag is present, treating '-flag' and '--flag'
    as equivalent.
    """
    try:
        target = _normalize_flag_name(flag)
        for tok in argv:
            if not tok.startswith("-"):
                continue
            if _normalize_flag_name(tok) == target:
                return True
        return False
    except Exception:
        return False

# --- Specified Proteins Mode helpers (NEW) -----------------------------------
from typing import Iterable

def _norm_pdb_id(token: str) -> Optional[str]:
    """
    Normalize a user token to a 4-char PDB ID (uppercase).
    Accepts bare IDs (2HYY), QoL flags (--2HYY or -2HYY), and filenames (2HYY.pdb).
    Returns None if it cannot produce a 4-char alnum ID.
    """
    if not token:
        return None
    t = str(token).strip()
    # Strip any leading dashes (one or two)
    while t.startswith("-"):
        t = t[1:]
    t = os.path.basename(t)
    if t.lower().endswith(".pdb"):
        t = t[:-4]
    t = t.replace("_cleaned", "")
    t = t.upper()
    if len(t) >= 4:
        cand = t[:4]
        return cand if cand.isalnum() else None
    return None


def _split_ids(s: str) -> list[str]:
    """Split a comma/whitespace separated string into normalized 4-char IDs."""
    if not s:
        return []
    parts = s.replace(",", " ").split()
    out = []
    for p in parts:
        nid = _norm_pdb_id(p)
        if nid:
            out.append(nid)
    return out

def _dedupe_order(seq: Iterable[str]) -> list[str]:
    """De-duplicate while preserving first-seen order."""
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out



from typing import Iterable

def _parse_specified_proteins(argv, cfg) -> tuple[list[str], str]:
    """
    Resolve requested PDB IDs with precedence CLI > ENV > CFG.

    CLI:
      --pdb  XIAP        (repeatable)
      -pdb   XIAP        (short form, repeatable)
      --pdbs "XIAP,1BN1" (comma/space separated)
      -pdbs  "XIAP 1BN1" (short form)
      --XIAP, -XIAP      (QoL: any --<4char> alnum)
      XIAP, xiap.pdb     (bare tokens)

    ENV:
      ONLY_PDBS="XIAP 1BN1"

    CFG:
      SPECIFIED_PROTEINS: JSON list or string "XIAP, 1BN1"

    Returns: (normalized_ids, source or "")
    """
    cli_ids: list[str] = []
    consumed_value_idx: set[int] = set()

    # --- --pdb / -pdb (repeatable) ---
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in ("--pdb", "-pdb"):
            if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                nid = _norm_pdb_id(argv[i + 1])
                if nid:
                    cli_ids.append(nid)
                    consumed_value_idx.add(i + 1)
            i += 2
            continue
        i += 1

    # --- --pdbs / -pdbs "XIAP,1BN1" ---
    for flag in ("--pdbs", "-pdbs"):
        try:
            if flag in argv:
                j = argv.index(flag)
                if j + 1 < len(argv) and not argv[j + 1].startswith("-"):
                    cli_ids.extend(_split_ids(argv[j + 1]))
                    consumed_value_idx.add(j + 1)
        except Exception:
            # be robust to weird argv
            pass

    # --- QoL: --XIAP / -XIAP style (exact length, 4-char alnum) ---
    for tok in argv:
        low = tok.lower()
        # don't treat fast/-fast/--fast as a PDB short-form token
        if low in ("fast", "-fast", "--fast"):
            continue
        if (tok.startswith("--") and len(tok) == 6) or (tok.startswith("-") and len(tok) == 5):
            nid = _norm_pdb_id(tok)
            if nid:
                cli_ids.append(nid)

    # --- Bare tokens: XIAP, xiap.pdb, XIAP_cleaned.pdb ---
    #
    # We skip:
    #   - argv[0] (script name)
    #   - tokens we've already consumed as values to known flags
    #   - anything starting with '-' (flags)
    #
    for idx, tok in enumerate(argv[1:], start=1):
        if idx in consumed_value_idx:
            continue
        if tok.startswith("-"):
            continue

        base = os.path.basename(tok)

        # Strip common suffix patterns
        lower = base.lower()
        if lower.endswith("_cleaned.pdb"):
            core = base[:-len("_cleaned.pdb")]
        elif lower.endswith(".pdb"):
            core = base[:-4]
        else:
            core = base

        core = core.strip()
        if not core:
            continue

        # Require exactly 4 alnum chars to avoid grabbing e.g. run-id strings
        if len(core) == 4 and core.isalnum():
            nid = _norm_pdb_id(core)
            if nid:
                cli_ids.append(nid)

    if cli_ids:
        return _dedupe_order(cli_ids), "CLI"

    # --- ENV ---
    env_val = os.environ.get("ONLY_PDBS", "").strip()
    if env_val:
        return _dedupe_order(_split_ids(env_val)), "ENV"

    # --- CFG ---
    val = cfg.get("SPECIFIED_PROTEINS", "")
    if isinstance(val, list):
        cfg_ids = [_norm_pdb_id(x) for x in val]
        cfg_ids = [x for x in cfg_ids if x]
        if cfg_ids:
            return _dedupe_order(cfg_ids), "CFG"
        return [], ""
    s = str(val or "").strip()
    if s:
        return _dedupe_order(_split_ids(s)), "CFG"

    return [], ""




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
            "OUTPUT_LIGANDS_DIR": smoke_root / "prepped_ligands",
            "PREPPED_LIGANDS_ROOT": smoke_root / "prepped_ligands",
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

        receptor_path = receptor_file(paths.pdb_id, variant="APO", ph_tag=ph_label, legacy=False)
        receptor_path.parent.mkdir(parents=True, exist_ok=True)
        if not receptor_path.exists():
            receptor_path.write_text("RECEPTOR", encoding="utf-8")

        stage_info = {"name": "smoke_stage", "exhaustiveness": 8, "num_modes": 9, "verbosity": 0}
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
    if is_resume and not cli_run_id:
        print("ERROR: -resume requires --run-id <RUN_ID>", file=sys.stderr)
        sys.exit(2)
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
        waters_set = {str(tok).strip().upper() for tok in getattr(alias_sets, "waters", set()) if str(tok).strip()}
    cofactors_set = {str(tok).strip().upper() for tok in getattr(rules, "cofactors", set()) if str(tok).strip()}
    elements_set = {str(tok).strip().upper() for tok in getattr(rules, "elem_tokens_canonical", set()) if str(tok).strip()}
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
    log.info("[ph_ensemble.mode] enabled=%s scope=%s radius=%s", cfg.PH_ENSEMBLE, getattr(cfg, "PH_SCOPE", "auto"), getattr(cfg, "PH_RADIUS", 1000000.0))

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
                    "[resume] Failed to parse stored argv=%r", stored_argv, exc_info=True
                )

        listed = command_section.get("pdb_list") or []
        if isinstance(listed, (list, tuple)):
            resume_protein_ids.extend(str(x).strip().upper() for x in listed if str(x).strip())

        if not resume_protein_ids:
            summary_section = resume_manifest.get("summary") or {}
            summary_list = summary_section.get("total_protein_list") or []
            if isinstance(summary_list, (list, tuple)):
                resume_protein_ids.extend(str(x).strip().upper() for x in summary_list if str(x).strip())

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
                    variant_token = (variant_part or "legacy").strip().lower() or "legacy"
                    ckey = (pdb_token, variant_token)
                    if not resume_protein_ids:
                        resume_protein_ids.append(pdb_token)
                    aggregated.setdefault(ckey, True)
                    status = (entry.get("status") or "").strip().lower()
                    if status != "completed":
                        aggregated[ckey] = False
                except Exception:
                    logging.warning(
                        "[resume.lookup.skip] key=%r entry=%r", key, entry, exc_info=True
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

    init_config_run_dir(cfg, run_id=cfg.get("RUN_ID"), reset=cfg.get("RESET_CONFIGS"),
                        logger=logging.getLogger("run"))
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
    cfg.setdefault("SINGLE_LIGAND_SEARCH_ORDER", "fda_library,per_protein,global")
    cfg.setdefault("SINGLE_LIGAND_ALLOW_PREFIX", False)
    cfg.setdefault("SINGLE_LIGAND_MANIFEST_ONLY", True)
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
    cfg.setdefault("FDA_MAPPING_CSV", str(Path(__file__).with_name("fda_mapping_from_pdbqt.csv")))

    # CLI > ENV > CFG precedence
    cli_single = _parse_single_from_cli(argv_for_parsing)
    env_single = os.environ.get("SINGLE_LIGAND", "").strip()
    cfg_single = str(cfg.get("SINGLE_LIGAND", "")).strip()

    effective_single = next((x for x in (cli_single, env_single, cfg_single) if x), "")
    cfg["_EFFECTIVE_SINGLE_LIGAND"] = effective_single
    if effective_single:
        print(f"[config] SINGLE_LIGAND effective='{effective_single}' "
              f"(order=CLI>{'ENV' if env_single else ''}>{'CFG' if cfg_single else ''})")

    # --- Library subfolder selection -----------------------------------
    cfg.setdefault("LIBRARY_SUBDIR_DEFAULT", "fda_library")
    cfg.setdefault("HMDB_LIBRARY_SUBDIR", "hmdb")
    cfg.setdefault("TEST_MODE_ENABLE", "off")
    # Accept dict or JSON-ish string
    if "TEST_LIBRARY_MAP" not in cfg:
        cfg["TEST_LIBRARY_MAP"] = {}

    # DeepCoy defaults (decoy autogen remains off unless toggled)
    cfg.setdefault("DEEPCOY_DUDS_DIR", "/home/michael/atlas/code/protein_automation/DeepCoy_duds")
    cfg.setdefault("DEEPCOY_PYTHON", "/home/michael/atlas/anaconda3/envs/DeepCoy-env/bin/python")
    cfg.setdefault("DEEPCOY_INPUT_PDB_DIR", "/home/michael/atlas/code/protein_automation/input_pdbs")
    cfg.setdefault(
        "DEEPCOY_OUT_ROOT",
        os.path.join(
            cfg.get("LIGAND_EXTRACTED_DIR", "/home/michael/atlas/code/protein_automation/extracted_ligands"),
            "deepcoy",
        ),
    )
    cfg.setdefault("DEEPCOY_WORK_ROOT", os.path.join(cfg["DEEPCOY_DUDS_DIR"], "deepcoy_work"))
    cfg.setdefault("DEEPCOY_ENABLE_AUTOGEN_SDF", "on")
    cfg.setdefault("DEEPCOY_ENABLE_AUTOGEN_PDBQT", "on")
    cfg.setdefault("DEEPCOY_FALLBACK_SMILES", "")
    cfg.setdefault("DEEPCOY_RESTRICT_DATA", 0)
    cfg.setdefault("DEEPCOY_PREPPED_SUBDIR", "deepcoy")
    cfg.setdefault("DEEPCOY_LIGPREP_FORCE", "off")
    cfg.setdefault("DEEPCOY_USE_AS_DUD_LIBRARY", "on")
    cfg.setdefault("DEEPCOY_FORCE", "off")
    env_deepcoy_toggle = os.environ.get("DEEPCOY_ENABLE_AUTOGEN_SDF")
    if env_deepcoy_toggle is not None:
        cfg["DEEPCOY_ENABLE_AUTOGEN_SDF"] = env_deepcoy_toggle

    # Small FDA test library override.
    # When called with: python main.py -test -test-fda -fast
    # we want to use a tiny FDA test library instead of the full FDA set.
    #
    # By default this points at "fda_test_library_10", which should correspond to:
    #   prepped_ligands/fda_test_library_10/
    #   extracted_ligands/fda_test_library_10/
    cfg.setdefault("TEST_FDA_LIBRARY_SUBDIR", "fda_test_library_10")

    if _cli_has(argv_for_parsing, "-test-fda") or _cli_has(argv_for_parsing, "--test-fda"):
        cfg["LIBRARY_SUBDIR_DEFAULT"] = cfg.get(
            "TEST_FDA_LIBRARY_SUBDIR",
            "fda_test_library_10",
        )
        print(
            f"[config] TEST_FDA_LIBRARY enabled: "
            f"LIBRARY_SUBDIR_DEFAULT={cfg['LIBRARY_SUBDIR_DEFAULT']}"
        )
    # --- Specified Proteins Mode ---------------------------------------
    cfg.setdefault("SPECIFIED_PROTEINS", "")
    requested_ids, _sel_src = _parse_specified_proteins(argv_for_parsing, cfg)
    cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] = requested_ids
    print(f"[config] SPECIFIED_PROTEINS effective={requested_ids} (precedence: CLI>ENV>CFG)")
    # --- Fast mode: force exhaustiveness=1 everywhere ---
    cfg["FAST_MODE"] = _parse_fast_flag(argv_for_parsing) or bool(cfg.get("FAST_MODE", False))
    if _parse_force_deepcoy_flag(argv_for_parsing):
        cfg["DEEPCOY_FORCE"] = "on"
    if cfg["FAST_MODE"]:
        print("[config] FAST_MODE effective=True (exhaustiveness=1)")

    # --- Control consensus (multi-engine control docking) -------------
    cfg.setdefault("CONTROL_CONSENSUS", False)
    control_consensus = False
    try:
        control_consensus = _to_bool(cfg.get("CONTROL_CONSENSUS", False))
    except Exception:
        control_consensus = bool(cfg.get("CONTROL_CONSENSUS", False))

    if _cli_has(argv_for_parsing, "-control-consensus") or _cli_has(argv_for_parsing, "--control-consensus"):
        control_consensus = True

    cfg["CONTROL_CONSENSUS"] = bool(control_consensus)
    if cfg["CONTROL_CONSENSUS"]:
        print("[config] CONTROL_CONSENSUS effective=True (multi-engine control docking enabled)")

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

    # --- Center selection knobs (safe defaults) ---
    cfg.setdefault("CENTER_MODE", "control-first")  # ["control-first","hybrid","library-first"]
    cfg.setdefault("CONTROL_BLACKLIST", "GOL,EDO,PG4,MPD,ACT,SO4,PO4,CL,NA,CA")
    cfg.setdefault("CONTROL_MIN_HEAVY_ATOMS", 10)
    cfg.setdefault("CONTROL_ANCHOR_MIN_VALID_RATE", 0.10)  # if current cluster has control hits + =10% valid, anchor
    cfg.setdefault("ALLOW_SWITCH_FROM_CONTROL", True)
    cfg.setdefault("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False)
    cfg.setdefault("SWITCH_AWAY_FROM_CONTROL_MIN_BOOST", 2.5)  # kcal/mol median boost needed to leave control
    # Optional lock score gate (kcal/mol). Use a large positive number (or remove) to lock on RMSD alone.
    cfg.setdefault("CONTROL_LOCK_SCORE_MAX", -6.0)
    cfg.setdefault("CONTROL_LOCK_MIN_HITS", 1)  # require = this many validated controls
    cfg.setdefault("CONTROL_LOCK_CENTER_MAX_DIST", 4.0)  # A; control centroid must be within this of center

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
    cfg.setdefault("QUIET_CONSOLE", False)   # console shows WARN+ only; file keeps DEBUG
    cfg.setdefault("VINA_VERBOSITY", 2)     # 0=minimal, 1=normal, 2=verbose
    cfg.setdefault("FILTER_VINA_STDOUT", False)  # reserved if we need extra filtering later

    #RMSD PARAMETERS
    cfg.setdefault("SELF_RMSD_MAX_ANG", 2.0)
    cfg.setdefault("SELF_RMSD_REQUIRE_FOR_CONTROLS", True)  # reserved for future stricter gating
    cfg.setdefault("EARLY_EXIT_MAX_MODELS", 3)  # validate at most N poses, stop on first PASS
    cfg.setdefault("MAX_RETRY_SECONDS_PER_LIGAND", 300)  # wall-clock for retries/validation per ligand

    params = get_recenter_params(cfg)
    stages = define_docking_stages(cfg.get("DOCKING_MODE", "discovery").lower())
    print("current docking mode is ", cfg.get("DOCKING_MODE"))

    # Discover all candidate PDB files (unchanged default behavior)
    pdb_files = [
        f for f in os.listdir(cfg["INPUT_DIR"])
        if f.lower().endswith(".pdb") and "_nolig" not in f.lower()
    ]

    # Build an index: PDBID (4-char, upper) -> filename
    id_index: dict[str, str] = {}
    for f in pdb_files:
        base = os.path.splitext(f)[0].replace("_cleaned", "")
        nid = _norm_pdb_id(base)
        if nid:
            # preserve first occurrence to retain directory order
            id_index.setdefault(nid, f)

    req = list(cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS", []) or [])
    if is_resume and not req and resume_protein_ids:
        req = list(resume_protein_ids)
    if req:
        # Compute present/missing and apply filter in user-specified order
        hits = [nid for nid in req if nid in id_index]
        miss = [nid for nid in req if nid not in id_index]

        print(f"[filter.proteins] mode=on requested={len(req)} present={len(hits)} missing={len(miss)} ? {hits}")
        for m in miss:
            print(f"WARNING: requested PDB '{m}' not found under INPUT_DIR={cfg['INPUT_DIR']} or was excluded (_nolig).")

        if not hits:
            print("ERROR: No requested proteins found. Exiting with status 2 to avoid a no-op run.")
            sys.exit(2)

        # Restrict queue to the selected files, preserving user order
        pdb_files = [id_index[nid] for nid in hits]
        print("Selected proteins (Specified Proteins Mode): " + ", ".join(hits))
    else:
        print(f"[filter.proteins] mode=off requested=0 present={len(pdb_files)} missing=0 ? []")



    # --- Test-mode protein filter: keep only PDBs listed in TEST_LIBRARY_MAP ---
    test_mode = _resolve_test_mode(cfg)
    if test_mode != "off":
        raw_map = cfg.get("TEST_LIBRARY_MAP", {})
        test_map = _coerce_test_map(raw_map)
        try:
            cfg["_TEST_LIBRARY_CANONICAL"] = {
                str(k).upper(): str(v) for k, v in getattr(test_map, "items", lambda: [])()
            }
        except Exception:
            cfg["_TEST_LIBRARY_CANONICAL"] = {}

        # Normalize all TEST_LIBRARY_MAP keys to canonical 4-char uppercase PDB IDs.
        # This makes matching robust to case and minor suffix differences.
        test_keys = set()
        for k in getattr(test_map, "keys", lambda: [])():
            nid = _norm_pdb_id(str(k))
            if nid:
                test_keys.add(nid)

        if test_keys:
            # Always honor explicitly specified proteins even in test mode.
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
                    f"[test-mode] Enabled mode={test_mode}; restricting to "
                    f"{len(kept)} PDBs from TEST_LIBRARY_MAP keys"
                    + (" (plus specified proteins)." if specified_keys else ".")
                )
                for s in skipped:
                    print(f"[test-mode] Skipping {s} (not in TEST_LIBRARY_MAP).")

            pdb_files = kept
        else:
            print("[test-mode] TEST_LIBRARY_MAP empty/invalid; no extra filtering applied.")
    else:
        cfg["_TEST_LIBRARY_CANONICAL"] = {}


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
    router_legacy = (mode == "legacy")
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
                        if test_mode != "off":
                            lib_map = cfg_for_pdb.get("_TEST_LIBRARY_CANONICAL", {}) or {}
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

                except Exception as exc:
                    # Per-PDB failure handling
                    exc_type = type(exc).__name__
                    exc_msg = str(exc)
                    traceback_str = traceback.format_exc()

                    fail_log_path = Path(failed_root) / f"{pdb_id}.{label}.log"

                    # Write a dedicated failure log for this PDB+variant
                    with fail_log_path.open("w", encoding="utf-8") as fh:
                        fh.write(
                            f"[FAILED PDB]\n"
                            f"  pdb_id        = {pdb_id}\n"
                            f"  variant_label = {label}\n"
                            f"  mode          = {mode}\n"
                            f"  pdb_file      = {pdb_file}\n"
                            f"  exception     = {exc_type}: {exc_msg}\n\n"
                            f"[TRACEBACK]\n"
                            f"{traceback_str}\n"
                        )
                        traceback.print_exc(file=fh)

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
                    cfg_local["MAX_PARALLEL_JOBS"] = 1
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
        _maybe_run_scorch_rescore(cfg, run_id, verbose=False)
    except Exception:
        logging.warning("[scorch-rescore.hook] action=skip reason=unexpected_exception", exc_info=True)

    try:
        failed_lookup = {(entry[0], entry[1]) for entry in failed_entries}
        for variant in variants:
            label = "legacy" if variant is None else str(variant).lower()
            variant_token = None if variant is None else str(variant).upper()
            legacy_mode = variant is None
            for pdb_file in pdb_files:
                pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()
                if (pdb_id, label) in failed_lookup:
                    continue
                _maybe_run_mmgbsa_for_pdb(
                    cfg,
                    pdb_file,
                    pdb_id,
                    variant_token,
                    run_id,
                    test_mode,
                    legacy_mode,
                )
    except Exception:
        if _to_bool(cfg.get("MMGBSA_STRICT", False)):
            raise
        logging.warning("[mmgbsa.pipeline] action=skip reason=unexpected_exception", exc_info=True)

    try:
        _maybe_run_dud_eval(cfg, run_id, pdb_files)
    except Exception:
        logging.warning("[dud-eval.invoke] action=skip reason=unexpected_exception", exc_info=True)

    try:
        _log_rescore_verification(run_id)
    except Exception:
        logging.warning("[post-check.invoke] action=skip reason=unexpected_exception", exc_info=True)


def _maybe_run_dud_eval(cfg: Mapping[str, Any], run_id: str, pdb_files: list[str]) -> None:
    """
    Best-effort post-run DUD evaluator. Never raises.
    """
    logger = logging.getLogger("dud-eval")
    if not run_id:
        logger.info("[dud-eval.skip] reason=missing_run_id")
        return

    test_mode = str(cfg.get("TEST_MODE_ENABLE", "off")).lower()
    if "dud" not in test_mode:
        logger.info("[dud-eval.skip] reason=test_mode_off test_mode=%s", test_mode)
        return
    if cfg.get("NO_LIBRARY_DOCKING"):
        logger.info("[dud-eval.skip] reason=no_library_docking")
        return

    repo_root = Path(__file__).resolve().parent
    preferred = repo_root / "analysis" / "dud_eval.py"
    fallback = repo_root / "dud_eval.py"
    dud_eval_path = preferred if preferred.exists() else fallback
    if not dud_eval_path.exists():
        logger.warning("[dud-eval.skip] reason=script_missing path=%s", dud_eval_path)
        return

    cmd = [
        sys.executable,
        str(dud_eval_path),
        "--run-id",
        str(run_id),
        "--docked-root",
        "docked",
        "--out-dir",
        "analysis/dud_eval",
    ]

    if len(pdb_files) == 1:
        pdb_id = _norm_pdb_id(pdb_files[0])
        if pdb_id:
            cmd.extend(["--pdb-id", pdb_id])
            logger.info("[dud-eval.filter] pdb_id=%s", pdb_id)

    logger.info("[dud-eval.run] cmd=%s", shlex.join(cmd))
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        logger.warning(
            "[dud-eval.fail] run_id=%s returncode=%s", run_id, result.returncode
        )
    else:
        logger.info(
            "[dud-eval.done] run_id=%s returncode=%s out_root=%s",
            run_id,
            result.returncode,
            "analysis/dud_eval",
        )


def _maybe_run_scorch_rescore(cfg: Mapping[str, Any], run_id: str, verbose: bool = False) -> None:
    """
    Best-effort post-run SCORCH rescoring. Never raises.
    """
    logger = logging.getLogger("scorch-rescore-hook")
    if not run_id:
        logger.info("[scorch-rescore.skip] reason=missing_run_id")
        return

    repo_root = Path(__file__).resolve().parent
    script_path = repo_root / "rescoring_scorch.py"
    if not script_path.exists():
        logger.warning("[scorch-rescore.skip] reason=missing_script path=%s", script_path)
        return

    def _as_int(val: Any, fallback: int) -> int:
        try:
            return int(val)
        except Exception:
            return fallback

    total_cpu = _as_int(cfg.get("CPU") or cfg.get("MAX_PARALLEL_JOBS") or (os.cpu_count() or 1), os.cpu_count() or 1)
    max_jobs_cfg = _as_int(cfg.get("MAX_PARALLEL_JOBS") or total_cpu, total_cpu)
    total_cpu = max(1, total_cpu)
    max_jobs_cfg = max(1, max_jobs_cfg)

    jobs = max(1, min(max_jobs_cfg, total_cpu))
    threads = max(1, total_cpu // jobs)

    scorch_jobs = cfg.get("SCORCH_JOBS")
    scorch_threads = cfg.get("SCORCH_THREADS")
    if scorch_jobs is not None:
        jobs = max(1, min(_as_int(scorch_jobs, jobs), total_cpu))
        threads = max(1, total_cpu // jobs)
    if scorch_threads is not None:
        threads = max(1, min(_as_int(scorch_threads, threads), total_cpu))
        if jobs * threads > total_cpu:
            threads = max(1, total_cpu // jobs)

    cmd = [
        sys.executable,
        str(script_path),
        "--run-id",
        run_id,
        "--repo-root",
        str(repo_root),
        "--threads",
        str(threads),
        "--jobs",
        str(jobs),
    ]
    if verbose or logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        cmd.append("--verbose")

    logger.info(
        "[scorch-rescore.invoke] run_id=%s jobs=%d threads=%d total_cpu=%d cmd=%s",
        run_id,
        jobs,
        threads,
        total_cpu,
        shlex.join(str(c) for c in cmd),
    )
    try:
        proc = subprocess.run(cmd, cwd=str(repo_root), check=False)
    except Exception:
        logger.warning("[scorch-rescore.invoke] action=skip reason=execution_failed", exc_info=True)
        return

    if proc.returncode != 0:
        logger.warning(
            "[scorch-rescore.result] status=failed run_id=%s returncode=%d",
            run_id,
            proc.returncode,
        )
    else:
        logger.info("[scorch-rescore.result] status=ok run_id=%s", run_id)


def _log_rescore_verification(run_id: str) -> None:
    """
    Lightweight best-effort check for consensus and SCORCH outputs.
    """
    logger = logging.getLogger("post-run-checks")
    if not run_id:
        logger.info("[post-check.skip] reason=missing_run_id")
        return

    repo_root = Path(__file__).resolve().parent
    docked_root = repo_root / "docked" / run_id
    post_root = repo_root / "post_docked" / run_id

    consensus_paths = [p for p in docked_root.rglob("consensus_docking_scores.csv") if p.is_file()]
    scorch_paths = [p for p in post_root.rglob("scorch_scores_all.csv") if p.is_file()]

    if consensus_paths:
        logger.info("[post-check.consensus] found=%d sample=%s", len(consensus_paths), consensus_paths[0])
    else:
        logger.warning("[post-check.consensus] reason=missing_files run_id=%s", run_id)

    if scorch_paths:
        logger.info("[post-check.scorch] found=%d sample=%s", len(scorch_paths), scorch_paths[0])
    else:
        logger.warning("[post-check.scorch] reason=missing_files run_id=%s", run_id)


def _mmgbsa_find_vina_config(
    run_id: str,
    pdb_id: str,
    stage_dir: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
) -> Optional[Path]:
    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage_dir,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    if not cfg_dir.exists():
        return None
    txt_files = sorted(cfg_dir.glob("*.txt"))
    return txt_files[0] if txt_files else None


def _mmgbsa_parse_vina_config(cfg_path: Path) -> tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]]]:
    center_vals: dict[str, float] = {}
    size_vals: dict[str, float] = {}
    try:
        lines = cfg_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None, None

    for line in lines:
        if "=" not in line:
            continue
        key, raw_val = line.split("=", 1)
        key = key.strip().lower()
        val_text = raw_val.strip()
        try:
            val = float(val_text)
        except Exception:
            continue
        if key in {"center_x", "center_y", "center_z"}:
            center_vals[key] = val
        elif key in {"size_x", "size_y", "size_z"}:
            size_vals[key] = val

    if len(center_vals) == 3:
        center = (center_vals["center_x"], center_vals["center_y"], center_vals["center_z"])
    else:
        center = None

    if len(size_vals) == 3:
        size = (size_vals["size_x"], size_vals["size_y"], size_vals["size_z"])
    else:
        size = None

    return center, size


def _mmgbsa_resolve_center_radius(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    stage_dir: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
) -> tuple[Optional[Tuple[float, float, float]], float, Optional[Path]]:
    fallback_radius = float(cfg.get("MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK", 6.0))
    cfg_path = _mmgbsa_find_vina_config(run_id, pdb_id, stage_dir, variant_token, ph_label, legacy_mode)
    if not cfg_path:
        return None, fallback_radius, None

    center, size = _mmgbsa_parse_vina_config(cfg_path)
    if center is None:
        return None, fallback_radius, cfg_path

    radius = fallback_radius
    if size is not None:
        radius = max(size) / 2.0
    return center, radius, cfg_path


def _mmgbsa_resolve_receptor_pdb(
    cfg: Mapping[str, Any],
    pdb_id: str,
    pdb_file: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
) -> Optional[Path]:
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=pdb_file)
    candidates: List[Path] = []

    if ph_label:
        ensemble_dir = router_ph_ensemble_dir(pdb_id, variant=variant_token, legacy=legacy_mode)
        tag = f"{pdb_id}_{ph_label}"
        candidates.append(ensemble_dir / f"{tag}.withH.pdb")
        candidates.append(ensemble_dir / f"{tag}.pdb")

    candidates.append(paths.receptor_cleaned_pdb(variant_token))
    candidates.append(paths.receptor_dir(variant_token) / f"{pdb_id}.pdb")

    for cand in candidates:
        if cand.exists():
            return cand
    return None


def _mmgbsa_normalize_ligand_base(name: str) -> str:
    s = str(name or "").strip()
    if not s:
        return ""
    s = re.sub(r"\.(pdbqt|mol2|sdf)$", "", s, flags=re.IGNORECASE)
    s = s.replace(".sanitized", "")
    return s


def _mmgbsa_sdf_base_candidates(stem: str, stage_dir_label: str) -> List[str]:
    candidates: set[str] = set()
    base = stem.replace(".sanitized", "")
    candidates.add(base)

    if stage_dir_label:
        for sep in ("_", "__"):
            suffix = f"{sep}{stage_dir_label}"
            if base.endswith(suffix):
                candidates.add(base[: -len(suffix)])

    for pattern in (
        r"(_+gnina_stage\d+)$",
        r"(_+ledock_stage\d+)$",
        r"(_+dock6_stage\d+)$",
        r"(_+stage\d+)$",
    ):
        stripped = re.sub(pattern, "", base, flags=re.IGNORECASE)
        if stripped != base:
            candidates.add(stripped)

    cleaned = set()
    for cand in candidates:
        cleaned.add(re.sub(r"_+$", "", cand))
    return sorted(cleaned)


def _mmgbsa_default_pose_group_regexes() -> List[str]:
    return [
        r"(_pose\d+)$",
        r"(_rank\d+)$",
        r"(_conf\d+)$",
        r"(_model\d+)$",
        r"(_p\d+)$",
        r"(_stage\d+)(_pose\d+)$",
        r"(_gnina_stage\d+)(_pose\d+)$",
        r"(_ledock_stage\d+)(_pose\d+)$",
        r"(_dock6_stage\d+)(_pose\d+)$",
    ]


def _mmgbsa_pose_group_regexes(cfg: Mapping[str, Any], logger: logging.Logger) -> List[re.Pattern]:
    raw = str(cfg.get("MMGBSA_POSE_GROUP_REGEXES", "") or "").strip()
    patterns: List[str] = []
    if raw:
        for part in raw.split(","):
            part = part.strip()
            if part:
                patterns.append(part)
    if not patterns:
        patterns = _mmgbsa_default_pose_group_regexes()

    compiled: List[re.Pattern] = []
    for pat in patterns:
        try:
            compiled.append(re.compile(pat, flags=re.IGNORECASE))
        except re.error as exc:
            logger.warning(
                "[mmgbsa.pipeline] selection=pose_regex_invalid regex=%s err=%s",
                pat,
                exc,
            )
    if not compiled:
        compiled = [re.compile(pat, flags=re.IGNORECASE) for pat in _mmgbsa_default_pose_group_regexes()]
    return compiled


def _mmgbsa_pose_candidates(stem: str, stage_dir_label: str, regexes: List[re.Pattern]) -> List[str]:
    base = stem.replace(".sanitized", "")
    seeds: set[str] = {base}
    for regex in regexes:
        stripped = regex.sub("", base)
        if stripped != base:
            seeds.add(stripped)

    candidates: set[str] = set()
    for seed in seeds:
        for cand in _mmgbsa_sdf_base_candidates(seed, stage_dir_label):
            candidates.add(cand)
        for regex in regexes:
            stripped = regex.sub("", seed)
            if stripped != seed:
                for cand in _mmgbsa_sdf_base_candidates(stripped, stage_dir_label):
                    candidates.add(cand)

    cleaned: set[str] = set()
    for cand in candidates:
        cleaned.add(re.sub(r"_+$", "", cand))
    return sorted([c for c in cleaned if c])


def _mmgbsa_canonical_ligand_id(
    stem: str,
    stage_dir_label: str,
    regexes: List[re.Pattern],
    logger: logging.Logger,
) -> str:
    candidates = _mmgbsa_pose_candidates(stem, stage_dir_label, regexes)
    if not candidates:
        return _mmgbsa_normalize_ligand_base(stem)

    min_len = min(len(cand) for cand in candidates)
    shortest = sorted([cand for cand in candidates if len(cand) == min_len])
    chosen = shortest[0]
    if len(shortest) > 1:
        logger.warning(
            "[mmgbsa.pipeline] selection=ligand_id_collision stem=%s candidates=%s chosen=%s",
            stem,
            shortest,
            chosen,
        )
    return chosen


def _mmgbsa_group_pose_sdfs(
    sdfs: List[Path],
    stage_dir_label: str,
    regexes: List[re.Pattern],
    logger: logging.Logger,
) -> Dict[str, List[Path]]:
    grouped: Dict[str, List[Path]] = {}
    for sdf in sdfs:
        ligand_id = _mmgbsa_canonical_ligand_id(sdf.stem, stage_dir_label, regexes, logger)
        if not ligand_id:
            ligand_id = sdf.stem
        grouped.setdefault(ligand_id, []).append(sdf)
    return grouped


def _mmgbsa_pose_numeric_rank(stem: str) -> Optional[int]:
    primary = re.search(r"(?:pose|rank|conf|model|p)(\d+)$", stem, flags=re.IGNORECASE)
    if primary:
        try:
            return int(primary.group(1))
        except Exception:
            return None

    trailing = re.search(r"(\d+)$", stem)
    if trailing:
        try:
            return int(trailing.group(1))
        except Exception:
            return None
    return None


def _mmgbsa_pose_sort_key(path: Path, mode: str) -> Tuple[int, object, str]:
    stem = path.stem.lower()
    if mode == "name_numeric":
        rank = _mmgbsa_pose_numeric_rank(path.stem)
        if rank is not None:
            return (0, rank, stem)
        return (1, stem, stem)
    return (0, stem, stem)


def _mmgbsa_index_sdfs(sdfs: List[Path], stage_dir_label: str) -> Dict[str, Path]:
    mapping: Dict[str, Path] = {}
    for sdf in sdfs:
        for key in _mmgbsa_sdf_base_candidates(sdf.stem, stage_dir_label):
            if key and key not in mapping:
                mapping[key] = sdf
    return mapping


def _mmgbsa_load_reranked_bases(csv_path: Path, logger: logging.Logger) -> List[str]:
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return []

    rows: List[Tuple[Optional[int], int, str]] = []
    try:
        with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle)
            for idx, row in enumerate(reader):
                lig_base = row.get("ligand_base") or _mmgbsa_normalize_ligand_base(row.get("ligand", ""))
                if not lig_base:
                    continue
                rank_val = None
                rank_raw = str(row.get("final_rank", "")).strip()
                if rank_raw:
                    try:
                        rank_val = int(float(rank_raw))
                    except Exception:
                        rank_val = None
                rows.append((rank_val, idx, lig_base))
    except Exception as exc:
        logger.warning(
            "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=csv_read_error path=%s err=%s",
            csv_path,
            exc,
        )
        return []

    if not rows:
        return []

    rows.sort(key=lambda x: (x[0] if x[0] is not None else 1_000_000_000, x[1]))
    seen: set[str] = set()
    ordered: List[str] = []
    for _, _, base in rows:
        if base in seen:
            continue
        seen.add(base)
        ordered.append(base)
    return ordered


def _mmgbsa_select_sdfs(
    stage_dir: Path,
    max_ligands: Optional[int],
    poses_per_ligand: int,
    pose_sort_mode: str,
    pose_group_regexes: List[re.Pattern],
    reranked_csv: Optional[Path] = None,
    top_pct: Optional[float] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[Path], bool, bool, List[str]]:
    sdfs = sorted(stage_dir.glob("*.sdf"))
    if not sdfs:
        return [], False, False, []

    if logger is None:
        logger = logging.getLogger("mmgbsa.pipeline")

    used_reranked = False
    used_actives = False
    stage_label = stage_dir.name
    pose_groups = _mmgbsa_group_pose_sdfs(sdfs, stage_label, pose_group_regexes, logger)
    selected_group_map = pose_groups
    selected_ligands: List[str] = []

    if reranked_csv is not None and top_pct is not None and top_pct > 0:
        if not reranked_csv.exists():
            logger.warning(
                "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=missing_reranked_csv path=%s",
                reranked_csv,
            )
        else:
            ordered_bases = _mmgbsa_load_reranked_bases(reranked_csv, logger)
            if ordered_bases:
                pct_val = float(top_pct)
                if pct_val > 100.0:
                    pct_val = 100.0
                if pct_val <= 0.0:
                    pct_val = 0.0
                ordered_ligands: List[str] = []
                seen: set[str] = set()
                for base in ordered_bases:
                    ligand_id = _mmgbsa_canonical_ligand_id(base, stage_label, pose_group_regexes, logger)
                    if ligand_id in seen or ligand_id not in pose_groups:
                        continue
                    seen.add(ligand_id)
                    ordered_ligands.append(ligand_id)

                total = len(ordered_ligands)
                target = max(1, int(math.ceil(total * pct_val / 100.0))) if total else 0
                if ordered_ligands and target:
                    selected_ligands = ordered_ligands[:target]
                    used_reranked = True
                    logger.info(
                        "[mmgbsa.pipeline] selection=reranked_top_pct pct=%.3g total=%d target=%d selected=%d path=%s",
                        pct_val,
                        total,
                        target,
                        len(selected_ligands),
                        reranked_csv,
                    )
                else:
                    logger.warning(
                        "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=no_matching_sdfs path=%s stage_dir=%s",
                        reranked_csv,
                        stage_dir,
                    )
            else:
                logger.warning(
                    "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=empty_csv path=%s",
                    reranked_csv,
                )

    if not selected_ligands:
        actives = [p for p in sdfs if "actives_final" in p.name]
        if actives:
            used_actives = True
            selected_group_map = _mmgbsa_group_pose_sdfs(actives, stage_label, pose_group_regexes, logger)
            selected_ligands = sorted(selected_group_map.keys())
        else:
            selected_group_map = pose_groups
            selected_ligands = sorted(pose_groups.keys())

    if max_ligands is not None and max_ligands > 0:
        selected_ligands = selected_ligands[:max_ligands]

    selected_sdfs: List[Path] = []
    for ligand_id in selected_ligands:
        group = selected_group_map.get(ligand_id, [])
        if not group:
            continue
        sorted_group = sorted(group, key=lambda p: _mmgbsa_pose_sort_key(p, pose_sort_mode))
        if poses_per_ligand > 0:
            sorted_group = sorted_group[:poses_per_ligand]
        selected_sdfs.extend(sorted_group)

    logger.info(
        "[mmgbsa.pipeline] selection=poses stage_dir=%s ligand_ids=%d poses_per_ligand=%d poses_selected=%d reranked=%s actives_only=%s",
        stage_label,
        len(selected_ligands),
        poses_per_ligand,
        len(selected_sdfs),
        used_reranked,
        used_actives,
    )
    return selected_sdfs, used_reranked, used_actives, selected_ligands


def _mmgbsa_parse_delta_total(csv_path: Path) -> Optional[float]:
    if not csv_path.exists():
        return None
    try:
        lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None

    in_delta = False
    header = None
    frame_values: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("DELTA Energy Terms"):
            in_delta = True
            header = None
            continue
        if not in_delta:
            continue
        if header is None:
            header = [t.strip() for t in stripped.split(",")]
            continue
        values = [t.strip() for t in stripped.split(",")]
        if not header:
            return None
        try:
            idx = header.index("DELTA TOTAL")
        except ValueError:
            return None
        if idx >= len(values):
            continue
        try:
            frame_values.append(float(values[idx]))
        except Exception:
            continue
    if not frame_values:
        return None
    return float(frame_values[0])


def _mmgbsa_parse_delta_frames(csv_path: Path) -> Dict[str, object]:
    result = {
        "frame_values": [],
        "n_frames": 0,
        "frame_mean": None,
        "frame_median": None,
        "frame_sd": None,
        "frame_min": None,
        "frame_max": None,
    }
    if not csv_path.exists():
        return result
    try:
        lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return result

    in_delta = False
    header = None
    frames: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("DELTA Energy Terms"):
            in_delta = True
            header = None
            continue
        if not in_delta:
            continue
        if header is None:
            header = [t.strip() for t in stripped.split(",")]
            continue
        values = [t.strip() for t in stripped.split(",")]
        if not header:
            break
        try:
            idx = header.index("DELTA TOTAL")
        except ValueError:
            break
        if idx >= len(values):
            continue
        try:
            frames.append(float(values[idx]))
        except Exception:
            continue

    if not frames:
        return result

    result["frame_values"] = frames
    result["n_frames"] = len(frames)
    try:
        result["frame_mean"] = float(statistics.mean(frames))
    except Exception:
        result["frame_mean"] = None
    try:
        result["frame_median"] = float(statistics.median(frames))
    except Exception:
        result["frame_median"] = None
    try:
        result["frame_sd"] = float(statistics.pstdev(frames))
    except Exception:
        result["frame_sd"] = None
    try:
        result["frame_min"] = float(min(frames))
        result["frame_max"] = float(max(frames))
    except Exception:
        pass
    return result


def _mmgbsa_write_replicate_summary(summary_path: Path, rows: List[Dict[str, object]], mean: float, sd: float) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["replicate", "seed", "score", "ok", "notes"]
    tmp_path = summary_path.with_suffix(summary_path.suffix + ".part")
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        writer.writerow({"replicate": "mean", "score": f"{mean:.6g}", "ok": ""})
        writer.writerow({"replicate": "stddev", "score": f"{sd:.6g}", "ok": ""})
    os.replace(tmp_path, summary_path)


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
    seeds = _mmgbsa_five_rep_seeds(run_id, pdb_id, variant_dir, ph_label, stage_dir, ligand_stem)
    out_dat_name = str(cfg.get("MMGBSA_MMPBSA_OUT_DAT", "FINAL_RESULTS_MMPBSA.dat") or "FINAL_RESULTS_MMPBSA.dat")
    out_csv_name = str(cfg.get("MMGBSA_MMPBSA_OUT_CSV", "FINAL_RESULTS_MMPBSA.csv") or "FINAL_RESULTS_MMPBSA.csv")
    agg_dat = out_dir / out_dat_name
    agg_csv = out_dir / out_csv_name

    if not force and agg_dat.exists() and agg_csv.exists():
        existing_delta = parse_mmpbsa_delta_total(agg_csv)
        if existing_delta is not None:
            return {
                "ok": True,
                "delta_total": existing_delta,
                "results_csv": str(agg_csv),
                "results_dat": str(agg_dat),
                "replicates_ok": 5,
                "replicates_total": 5,
                "summary_path": str(out_dir / "mmgbsa_replicates_summary.csv"),
                "notes": "cached",
            }

    rep_rows: List[Dict[str, object]] = []
    scores: List[float] = []
    summary_path = out_dir / "mmgbsa_replicates_summary.csv"

    if not md_enabled:
        logger.warning(
            "[mmgbsa.replicate] MD_FIVE_REPLICATE enabled but MD disabled; running single MMGBSA and cloning results"
        )
        rep_dir = out_dir / "rep1"
        traj_result = make_mmgbsa_trajectory(
            complex_prmtop=topo["complex_prmtop"],
            complex_inpcrd=topo["complex_inpcrd"],
            out_dir=str(rep_dir),
            cfg=cfg,
            force=force,
            run_cpptraj=_to_bool(cfg.get("MMGBSA_CPPTRAJ_RUN", True)),
        )
        traj_path = Path(traj_result.get("trajout_path") or "")
        mmpbsa_dir = rep_dir / "mmpbsa"
        mmpbsa_dir.mkdir(parents=True, exist_ok=True)
        mm_res = run_mmgbsa(
            complex_prmtop=topo["complex_prmtop"],
            receptor_prmtop=topo["receptor_prmtop"],
            ligand_prmtop=topo["ligand_prmtop"],
            trajectory_path=str(traj_path),
            work_dir=str(mmpbsa_dir),
            cfg=cfg,
            force=force,
            run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
        )
        rep_csv = Path(mm_res.get("out_csv", ""))
        score = parse_mmpbsa_delta_total(rep_csv)
        rep_ok = score is not None and rep_csv.exists()
        rep_rows.append({"replicate": 1, "seed": seeds[0], "score": score if score is not None else "", "ok": rep_ok, "notes": ""})
        if not rep_ok:
            _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=0.0, sd=0.0)
            return {"ok": False, "notes": "rep1_failed", "replicates_ok": 0, "replicates_total": 5}
        scores.append(float(score))
        # clone outputs for other reps
        for idx in range(2, 6):
            clone_dir = out_dir / f"rep{idx}" / "mmpbsa"
            clone_dir.mkdir(parents=True, exist_ok=True)
            for src in (rep_csv, Path(mm_res.get("out_dat", "")), Path(mm_res.get("log_path", "")), Path(mm_res.get("input_path", ""))):
                if not src:
                    continue
                if not src.exists():
                    continue
                dest = clone_dir / src.name
                if dest.exists() and not force:
                    continue
                shutil.copy2(src, dest)
            rep_rows.append({"replicate": idx, "seed": seeds[idx - 1], "score": score, "ok": True, "notes": "cloned"})
            scores.append(float(score))
    else:
        for rep_idx, seed in enumerate(seeds, start=1):
            rep_dir = out_dir / f"rep{rep_idx}"
            md_result = run_implicit_md(
                complex_prmtop=topo["complex_prmtop"],
                complex_inpcrd=topo["complex_inpcrd"],
                out_dir=str(out_dir),
                cfg=cfg,
                replicate_index=rep_idx,
                seed=int(seed),
                force=force,
                run=True,
            )
            traj_path = Path(md_result.get("traj_path") or "")
            rep_ok = False
            score = None
            notes = ""
            if not md_result.get("ok"):
                notes = "md_failed"
            elif not traj_path.exists() or traj_path.stat().st_size == 0:
                notes = "missing_trajectory"
            else:
                mmpbsa_dir = rep_dir / "mmpbsa"
                mmpbsa_dir.mkdir(parents=True, exist_ok=True)
                mm_res = run_mmgbsa(
                    complex_prmtop=topo["complex_prmtop"],
                    receptor_prmtop=topo["receptor_prmtop"],
                    ligand_prmtop=topo["ligand_prmtop"],
                    trajectory_path=str(traj_path),
                    work_dir=str(mmpbsa_dir),
                    cfg=cfg,
                    force=force,
                    run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
                )
                rep_csv = Path(mm_res.get("out_csv", ""))
                score = parse_mmpbsa_delta_total(rep_csv)
                if score is not None and rep_csv.exists():
                    rep_ok = True
                else:
                    notes = "missing_outputs"
            rep_rows.append({"replicate": rep_idx, "seed": seed, "score": score if score is not None else "", "ok": rep_ok, "notes": notes})
            if rep_ok and score is not None:
                scores.append(float(score))

    if len(scores) != 5 or any(not r.get("ok") for r in rep_rows):
        _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=0.0, sd=0.0)
        return {
            "ok": False,
            "notes": "replicate_failed",
            "replicates_ok": len(scores),
            "replicates_total": 5,
            "summary_path": str(summary_path),
        }

    mean_val = float(statistics.mean(scores))
    try:
        sd_val = float(statistics.pstdev(scores))
    except Exception:
        sd_val = 0.0

    _mmgbsa_write_replicate_summary(summary_path, rep_rows, mean=mean_val, sd=sd_val)
    agg_paths = write_aggregated_mmpbsa_results(work_dir=out_dir, mean_score=mean_val, std_score=sd_val, cfg=cfg, force=force)

    logger.info(
        "[mmgbsa.aggregate] n=5 mean=%.6g sd=%.6g out_csv=%s out_dat=%s",
        mean_val,
        sd_val,
        agg_paths.get("out_csv"),
        agg_paths.get("out_dat"),
    )

    return {
        "ok": True,
        "delta_total": mean_val,
        "results_csv": agg_paths.get("out_csv", ""),
        "results_dat": agg_paths.get("out_dat", ""),
        "replicates_ok": 5,
        "replicates_total": 5,
        "summary_path": str(summary_path),
        "notes": "",
        "frame_mean": mean_val,
        "frame_sd": sd_val,
    }


def _mmgbsa_trimmed_mean(values: List[float], trim_fraction: float = 0.1) -> float:
    if not values:
        raise ValueError("no values for trimmed mean")
    if trim_fraction <= 0:
        return float(statistics.mean(values))
    sorted_vals = sorted(values)
    trim_n = int(len(sorted_vals) * trim_fraction)
    max_trim = max(0, (len(sorted_vals) - 1) // 2)
    trim_n = max(0, min(trim_n, max_trim))
    if trim_n > 0:
        trimmed = sorted_vals[trim_n:-trim_n] or sorted_vals
    else:
        trimmed = sorted_vals
    return float(statistics.mean(trimmed))


def _mmgbsa_apply_agg(
    values: List[float],
    method: str,
    logger: Optional[logging.Logger] = None,
    context: str = "mmgbsa",
) -> Optional[float]:
    if not values:
        return None
    method_norm = str(method or "mean").strip().lower()
    try:
        if method_norm == "median":
            return float(statistics.median(values))
        if method_norm == "min":
            return float(min(values))
        if method_norm == "trimmed_mean":
            return _mmgbsa_trimmed_mean(values, trim_fraction=0.1)
        return float(statistics.mean(values))
    except Exception as exc:
        if logger:
            logger.warning(
                "[mmgbsa.pipeline] aggregate_failed context=%s method=%s err=%s",
                context,
                method_norm,
                exc,
            )
    return None


def _mmgbsa_frame_aggregate(csv_path: Path, cfg: Mapping[str, Any], logger: logging.Logger) -> Dict[str, object]:
    parsed = _mmgbsa_parse_delta_frames(csv_path)
    frames: List[float] = list(parsed.get("frame_values") or [])
    total_frames = len(frames)
    try:
        frame_limit = int(cfg.get("MMGBSA_FRAME_LIMIT", 0))
    except Exception:
        frame_limit = 0
    if frame_limit > 0 and frames:
        frames = frames[:frame_limit]
    frames_used = len(frames)
    frame_method = str(cfg.get("MMGBSA_FRAME_AGG", "mean") or "mean").strip().lower()

    agg_delta = _mmgbsa_apply_agg(frames, frame_method, logger=logger, context="frame")

    subset_mean = None
    subset_median = None
    subset_sd = None
    try:
        subset_mean = float(statistics.mean(frames))
    except Exception:
        subset_mean = None
    try:
        subset_median = float(statistics.median(frames))
    except Exception:
        subset_median = None
    try:
        subset_sd = float(statistics.pstdev(frames))
    except Exception:
        subset_sd = None

    if frames_used:
        logger.info(
            "[mmgbsa.frame] csv=%s frames_total=%d frames_used=%d limit=%d method=%s delta=%s mean=%s sd=%s",
            csv_path,
            total_frames,
            frames_used,
            frame_limit,
            frame_method,
            f"{agg_delta:.6g}" if agg_delta is not None else "",
            f"{subset_mean:.6g}" if subset_mean is not None else "",
            f"{subset_sd:.6g}" if subset_sd is not None else "",
        )
    else:
        logger.warning(
            "[mmgbsa.frame] csv=%s frames_total=%d frames_used=%d limit=%d method=%s reason=missing_frames",
            csv_path,
            total_frames,
            frames_used,
            frame_limit,
            frame_method,
        )

    return {
        "delta_total": agg_delta,
        "frame_agg_method": frame_method,
        "frames_used": frames_used,
        "frames_total": total_frames,
        "frame_mean": subset_mean,
        "frame_median": subset_median,
        "frame_sd": subset_sd,
        "frame_values": frames,
        "frame_limit": frame_limit,
        "notes": "" if agg_delta is not None else "missing_delta",
    }


def _mmgbsa_append_summary(summary_path: Path, row: Dict[str, object]) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    exists = summary_path.exists() and summary_path.stat().st_size > 0
    fieldnames = [
        "stage_dir",
        "ligand_stem",
        "delta_total",
        "results_csv",
        "ok",
        "notes",
        "frame_agg_method",
        "frames_used",
        "frames_total",
        "frame_mean",
        "frame_sd",
        "frame_median",
        "rep_agg_method",
        "replicates_ok",
        "replicates_total",
    ]
    with summary_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", restval="")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def _mmgbsa_write_pose_aggregate(
    summary_path: Path,
    out_path: Path,
    stage_dir_label: str,
    pose_group_regexes: List[re.Pattern],
    agg_method: str,
    logger: logging.Logger,
) -> None:
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        logger.info("[mmgbsa.pipeline] aggregate=skip reason=missing_summary path=%s", summary_path)
        return

    groups: Dict[Tuple[str, str], Dict[str, object]] = {}
    try:
        with summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                stage_dir = str(row.get("stage_dir") or stage_dir_label or "").strip()
                ligand_stem = str(row.get("ligand_stem") or "").strip()
                if not ligand_stem:
                    continue
                ligand_id = _mmgbsa_canonical_ligand_id(ligand_stem, stage_dir, pose_group_regexes, logger)
                key = (stage_dir, ligand_id)
                entry = groups.setdefault(key, {"total": 0, "ok": []})
                entry["total"] = int(entry.get("total", 0)) + 1

                ok_val = _to_bool(row.get("ok", False))
                delta_val = None
                if ok_val:
                    raw = str(row.get("delta_total", "")).strip()
                    if raw:
                        try:
                            delta_val = float(raw)
                        except Exception:
                            delta_val = None
                if ok_val and delta_val is not None:
                    entry["ok"].append((delta_val, ligand_stem))
    except Exception as exc:
        logger.warning(
            "[mmgbsa.pipeline] aggregate=skip reason=read_error path=%s err=%s",
            summary_path,
            exc,
        )
        return

    agg_method_norm = str(agg_method or "min").strip().lower()
    if agg_method_norm not in ("min", "median"):
        logger.warning(
            "[mmgbsa.pipeline] aggregate=unknown_method method=%s fallback=min",
            agg_method_norm,
        )
        agg_method_norm = "min"

    rows_out: List[Dict[str, object]] = []
    for (stage_dir, ligand_id), entry in sorted(groups.items()):
        ok_list = list(entry.get("ok", []))
        n_total = int(entry.get("total", 0))
        n_ok = len(ok_list)
        agg_delta = ""
        best_pose = ""
        included: List[str] = []

        if n_ok:
            deltas = [val for val, _ in ok_list]
            if agg_method_norm == "median":
                agg_val = float(statistics.median(deltas))
                best = min(ok_list, key=lambda x: (abs(x[0] - agg_val), x[0], x[1]))
            else:
                agg_val = min(deltas)
                best = min(ok_list, key=lambda x: (x[0], x[1]))
            agg_delta = f"{agg_val:.6g}"
            best_pose = best[1]
            included = [stem for _, stem in ok_list]

        rows_out.append(
            {
                "stage_dir": stage_dir,
                "ligand_id": ligand_id,
                "n_poses_total": n_total,
                "n_poses_ok": n_ok,
                "agg_method": agg_method_norm,
                "agg_delta": agg_delta,
                "best_pose_stem": best_pose,
                "pose_stems_included": ";".join(included),
            }
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    fieldnames = [
        "stage_dir",
        "ligand_id",
        "n_poses_total",
        "n_poses_ok",
        "agg_method",
        "agg_delta",
        "best_pose_stem",
        "pose_stems_included",
    ]
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)
    os.replace(tmp_path, out_path)
    logger.info(
        "[mmgbsa.pipeline] aggregate=pose_delta path=%s rows=%d method=%s",
        out_path,
        len(rows_out),
        agg_method_norm,
    )


def _mmgbsa_md_seed_list(cfg: Mapping[str, Any], n_reps: int, logger: logging.Logger) -> List[int]:
    try:
        base_seed = int(cfg.get("MMGBSA_MD_BASE_SEED", 12345))
    except Exception:
        base_seed = 12345

    mode = str(cfg.get("MMGBSA_MD_SEED_MODE", "increment") or "increment").strip().lower()
    if n_reps < 1:
        n_reps = 1

    seeds: List[int] = []
    if mode == "random":
        rng = random.Random(base_seed)
        for _ in range(n_reps):
            seeds.append(rng.randint(1, 2_147_483_647))
    else:
        for idx in range(n_reps):
            seeds.append(max(1, base_seed + idx))

    logger.info(
        "[mmgbsa.pipeline] md_seeds mode=%s base=%s reps=%d",
        mode,
        base_seed,
        n_reps,
    )
    return seeds


def _mmgbsa_five_rep_seeds(run_id: str, pdb_id: str, variant: str, ph_label: str, stage_dir: str, ligand: str) -> List[int]:
    context = f"{run_id}|{pdb_id}|{variant}|{ph_label}|{stage_dir}|{ligand}"
    base_seed = int(hashlib.md5(context.encode("utf-8")).hexdigest()[:8], 16)
    seeds: List[int] = []
    for idx in range(5):
        seeds.append(max(1, base_seed + idx * 10007))
    return seeds


def _mmgbsa_md_aggregate(deltas: List[float], method: str) -> Optional[float]:
    return _mmgbsa_apply_agg(deltas, method, logger=logging.getLogger("mmgbsa.pipeline"), context="replicate")


_MMGBSA_DEPRECATED_MD_KEYS = ("MMGBSA_MD_RUN", "MMGBSA_TRAJ_MODE")
_MMGBSA_DEPRECATED_MD_WARNED = False


def _mmgbsa_effective_md_config(cfg: Mapping[str, Any], logger: Optional[logging.Logger]) -> Dict[str, object]:
    global _MMGBSA_DEPRECATED_MD_WARNED
    deprecated_keys = []
    for key in _MMGBSA_DEPRECATED_MD_KEYS:
        if key in cfg and str(cfg.get(key)).strip() != "":
            deprecated_keys.append(key)

    md_enabled_raw = cfg.get("MMGBSA_MD_ENABLED", None)
    md_enabled_set = "MMGBSA_MD_ENABLED" in cfg and str(md_enabled_raw).strip() != ""
    if md_enabled_set:
        md_enabled = _to_bool(md_enabled_raw)
    else:
        traj_mode = str(cfg.get("MMGBSA_TRAJ_MODE", "") or "").strip().upper()
        md_enabled = traj_mode == "IMPLICIT_MD" or _to_bool(cfg.get("MMGBSA_MD_RUN", False))

    if deprecated_keys and not _MMGBSA_DEPRECATED_MD_WARNED and logger is not None:
        logger.warning(
            "[mmgbsa.pipeline] deprecated_keys=%s msg=Deprecated MMGBSA keys detected; please use MMGBSA_MD_ENABLED",
            ",".join(sorted(deprecated_keys)),
        )
        _MMGBSA_DEPRECATED_MD_WARNED = True

    return {
        "md_enabled": md_enabled,
        "effective_mode": "IMPLICIT_MD" if md_enabled else "ONEFRAME",
        "deprecated_keys": deprecated_keys,
    }


def _mmgbsa_apply_env_overrides(cfg: Mapping[str, Any], keys: Sequence[str]) -> Dict[str, Any]:
    updated = dict(cfg)
    for key in keys:
        value = os.environ.get(key)
        if value is None or str(value).strip() == "":
            continue
        updated[key] = value
    return updated


def _maybe_run_mmgbsa_for_pdb(
    cfg: Mapping[str, Any],
    pdb_file: str,
    pdb_id: str,
    variant_token: Optional[str],
    run_id: str,
    test_mode: str,
    legacy_mode: bool,
) -> None:
    logger = logging.getLogger("mmgbsa.pipeline")
    cfg = _mmgbsa_apply_env_overrides(
        cfg,
        [
            "MMGBSA_MAX_LIGANDS",
            "MMGBSA_RERANKED_TOP_PCT",
            "MMGBSA_TRAJ_MODE",
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES",
            "MMGBSA_MD_ENABLED",
            "MMGBSA_MD_RUN",
            "MMGBSA_MD_ENGINE",
            "MMGBSA_MD_NREPLICATES",
            "MMGBSA_MD_SEED_MODE",
            "MMGBSA_MD_BASE_SEED",
            "MMGBSA_MD_IGB",
            "MMGBSA_MD_SALTCON",
            "MMGBSA_MD_RESTRAIN_PROTEIN_HEAVY",
            "MMGBSA_MD_RESTRAINT_WT",
            "MMGBSA_MD_RESTRAINT_MASK",
            "MMGBSA_MD_DT_PS",
            "MMGBSA_MD_HEAT_PS",
            "MMGBSA_MD_EQUIL_PS",
            "MMGBSA_MD_PROD_PS",
            "MMGBSA_MD_TEMP0",
            "MMGBSA_MD_NTT",
            "MMGBSA_MD_GAMMA_LN",
            "MMGBSA_MD_NTC",
            "MMGBSA_MD_NTF",
            "MMGBSA_MD_FRAME_STRIDE_PS",
            "MMGBSA_MD_TRAJ_FORMAT",
            "MMGBSA_MD_TRAJ_NAME",
            "MMGBSA_MD_REP_AGG",
            "MMGBSA_MD_COPY_BEST_REPLICATE",
            "MMGBSA_FRAME_AGG",
            "MMGBSA_FRAME_LIMIT",
        ],
    )
    if not _to_bool(cfg.get("MMGBSA_ENABLED", False)):
        logger.info("[mmgbsa.pipeline] action=skip reason=disabled")
        return

    variant_dir = (variant_token or "legacy").upper()
    post_root = Path(cfg.get("OVERALL_DIR", ".")) / "post_docked" / run_id / pdb_id
    post_base = post_root if legacy_mode else post_root / variant_dir
    if not post_base.exists():
        logger.info("[mmgbsa.pipeline] action=skip reason=missing_post_docked path=%s", post_base)
        return

    try:
        max_ligands_val = int(cfg.get("MMGBSA_MAX_LIGANDS", 1))
    except Exception:
        max_ligands_val = 1
    if max_ligands_val <= 0:
        max_ligands_val = None

    try:
        top_pct_val = float(cfg.get("MMGBSA_RERANKED_TOP_PCT", 0))
    except Exception:
        top_pct_val = 0.0
    if top_pct_val <= 0:
        top_pct_val = None

    try:
        poses_per_ligand = int(cfg.get("MMGBSA_POSES_PER_LIGAND", 1))
    except Exception:
        poses_per_ligand = 1
    if poses_per_ligand <= 0:
        poses_per_ligand = 1

    pose_sort_mode = str(cfg.get("MMGBSA_POSE_SORT_MODE", "name_numeric") or "name_numeric").strip().lower()
    pose_group_regexes = _mmgbsa_pose_group_regexes(cfg, logger)

    agg_enabled = _to_bool(cfg.get("MMGBSA_AGGREGATE_PER_LIGAND", False))
    agg_method = str(cfg.get("MMGBSA_AGG_METHOD", "min") or "min").strip().lower()
    agg_output_name = str(cfg.get("MMGBSA_AGG_OUTPUT_CSV", "mmgbsa_pose_aggregate.csv") or "mmgbsa_pose_aggregate.csv")

    md_cfg = _mmgbsa_effective_md_config(cfg, logger)
    md_enabled = _to_bool(md_cfg.get("md_enabled", False))
    effective_mode = str(md_cfg.get("effective_mode", "ONEFRAME") or "ONEFRAME")
    five_reps = _to_bool(cfg.get("MD_FIVE_REPLICATE", False))
    try:
        md_reps = int(cfg.get("MMGBSA_MD_NREPLICATES", 1))
    except Exception:
        md_reps = 1
    if md_reps < 1:
        md_reps = 1
    md_rep_agg = str(cfg.get("MMGBSA_MD_REP_AGG", "mean") or "mean").strip().lower()
    md_copy_best = _to_bool(cfg.get("MMGBSA_MD_COPY_BEST_REPLICATE", False))
    frame_agg_method = str(cfg.get("MMGBSA_FRAME_AGG", "mean") or "mean").strip().lower()
    try:
        frame_limit = int(cfg.get("MMGBSA_FRAME_LIMIT", 0))
    except Exception:
        frame_limit = 0
    if frame_limit < 0:
        frame_limit = 0

    force = _to_bool(cfg.get("MMGBSA_FORCE", False))
    strict = _to_bool(cfg.get("MMGBSA_STRICT", False))
    stage_dir_name = str(cfg.get("MMGBSA_INPUT_STAGE_DIR", "stage1") or "stage1")
    strip_all_h_for_leap = _to_bool(cfg.get("MMGBSA_TLEAP_STRIP_ALL_H", True))
    map_hoh_to_wat = _to_bool(cfg.get("MMGBSA_TLEAP_MAP_HOH_TO_WAT", True))
    water_model = str(cfg.get("MMGBSA_TLEAP_WATER_MODEL", "tip3p") or "tip3p").strip().lower()
    if water_model not in {"tip3p"}:
        raise ValueError(f"MMGBSA_TLEAP_WATER_MODEL supports tip3p only (got {water_model})")
    mmgbsa_logger = ppm._get_logger()

    test_override: Dict[str, object] = {}
    if test_mode != "off" and top_pct_val is None:
        test_stage = stage_dir_name
        candidate_dir = post_base / "pH7_0" / test_stage
        if candidate_dir.is_dir():
            actives = sorted(candidate_dir.glob("actives_final*.sdf"))
            if actives:
                test_override = {"ph_dir": candidate_dir.parent, "sdfs": [actives[0]], "stage_dir": test_stage}
            else:
                sdfs = sorted(candidate_dir.glob("*.sdf"))
                if sdfs:
                    test_override = {"ph_dir": candidate_dir.parent, "sdfs": [sdfs[0]], "stage_dir": test_stage}

    if test_override:
        ph_dirs = [test_override["ph_dir"]]
    else:
        ph_dirs = sorted([p for p in post_base.iterdir() if p.is_dir()])

    if not ph_dirs:
        logger.info("[mmgbsa.pipeline] action=skip reason=no_ph_dirs path=%s", post_base)
        return

    amber_prefix = cfg.get("MMGBSA_AMBERTOOLS_PREFIX") or cfg.get("AMBERTOOLS_PREFIX")
    amber_prefix = str(amber_prefix).strip() if amber_prefix else None

    for ph_dir in ph_dirs:
        ph_label = ph_dir.name
        if test_override:
            sdfs = list(test_override["sdfs"])
            stage_dir_label = str(test_override["stage_dir"])
            pose_groups = _mmgbsa_group_pose_sdfs(sdfs, stage_dir_label, pose_group_regexes, logger)
            selected_ligands = sorted(pose_groups.keys())
            logger.info(
                "[mmgbsa.pipeline] selection=poses stage_dir=%s ligand_ids=%d poses_per_ligand=%d poses_selected=%d reranked=%s actives_only=%s",
                stage_dir_label,
                len(selected_ligands),
                poses_per_ligand,
                len(sdfs),
                False,
                False,
            )
        else:
            stage_dir_label = stage_dir_name
            stage_dir = ph_dir / stage_dir_label
            reranked_csv = ph_dir / "consensus_reranked_scorch.csv"
            sdfs, _, _, selected_ligands = _mmgbsa_select_sdfs(
                stage_dir,
                max_ligands_val,
                poses_per_ligand,
                pose_sort_mode,
                pose_group_regexes,
                reranked_csv=reranked_csv,
                top_pct=top_pct_val,
                logger=logger,
            )

        if not sdfs:
            logger.info(
                "[mmgbsa.pipeline] action=skip reason=no_sdfs pdb=%s variant=%s ph=%s stage=%s",
                pdb_id,
                variant_dir,
                ph_label,
                stage_dir_label,
            )
            continue

        expected_frames = None
        if md_enabled:
            try:
                prod_ps = float(cfg.get("MMGBSA_MD_PROD_PS", 100.0))
            except Exception:
                prod_ps = 100.0
            try:
                stride_ps = float(cfg.get("MMGBSA_MD_FRAME_STRIDE_PS", 2.0))
            except Exception:
                stride_ps = 2.0
            if stride_ps > 0:
                expected_frames = max(1, int(prod_ps / stride_ps))

        mmpbsa_interval = 1
        for key in ("MMGBSA_GENERAL_INTERVAL", "MMGBSA_MMPBSA_INTERVAL"):
            val = cfg.get(key, None)
            if val is None or str(val).strip() == "":
                continue
            try:
                mmpbsa_interval = int(val)
            except Exception:
                mmpbsa_interval = 1
            break

        try:
            mmpbsa_igb = int(cfg.get("MMGBSA_GB_IGB", 5))
        except Exception:
            mmpbsa_igb = 5
        try:
            mmpbsa_saltcon = float(cfg.get("MMGBSA_GB_SALTCON", 0.150))
        except Exception:
            mmpbsa_saltcon = 0.150
        logger.info(
            "[mmgbsa.pipeline] effective_cfg mode=%s md_enabled=%s expected_frames_per_rep=%s replicates=%d mmpbsa_interval=%d igb=%d saltcon=%.3f",
            effective_mode,
            md_enabled,
            expected_frames if md_enabled else "n/a",
            md_reps,
            mmpbsa_interval,
            mmpbsa_igb,
            mmpbsa_saltcon,
        )

        center, radius, cfg_path = _mmgbsa_resolve_center_radius(
            cfg,
            run_id,
            pdb_id,
            stage_dir_label,
            variant_token,
            ph_label,
            legacy_mode,
        )
        if center is None:
            logger.warning(
                "[mmgbsa.pipeline] action=skip reason=missing_center pdb=%s variant=%s ph=%s cfg=%s",
                pdb_id,
                variant_dir,
                ph_label,
                cfg_path,
            )
            continue

        receptor_candidates = [
            ph_dir / stage_dir_label / "receptor.pdb",
            ph_dir / "receptor.pdb",
        ]
        receptor_pdb = next((cand for cand in receptor_candidates if cand.exists()), None)
        if receptor_pdb is None:
            receptor_pdb = _mmgbsa_resolve_receptor_pdb(
                cfg,
                pdb_id,
                pdb_file,
                variant_token,
                ph_label,
                legacy_mode,
            )

        if receptor_pdb is None:
            logger.warning(
                "[mmgbsa.pipeline] action=skip reason=missing_receptor pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_dir,
                ph_label,
            )
            continue

        receptor_input = receptor_pdb
        if "post_docked" not in receptor_pdb.parts:
            receptor_input_dir = ph_dir / "mmgbsa_receptor_inputs"
            receptor_input_dir.mkdir(parents=True, exist_ok=True)
            receptor_input = receptor_input_dir / receptor_pdb.name
            if not receptor_input.exists() or force:
                try:
                    shutil.copy2(receptor_pdb, receptor_input)
                except Exception as exc:
                    logger.error(
                        "[mmgbsa.pipeline] action=skip reason=receptor_copy_failed pdb=%s variant=%s ph=%s err=%s",
                        pdb_id,
                        variant_dir,
                        ph_label,
                        exc,
                    )
                    if strict:
                        raise
                    continue
            logger.info(
                "[mmgbsa.pipeline] receptor_materialized source=%s dest=%s",
                receptor_pdb,
                receptor_input,
            )

        logger.info(
            "[mmgbsa.pipeline] start pdb=%s variant=%s ph=%s ligands=%d",
            pdb_id,
            variant_dir,
            ph_label,
            len(sdfs),
        )

        try:
            receptor_result = prep_mmgbsa_receptor(
                pdb_path=str(receptor_input),
                runid=run_id,
                center=center,
                radius=radius,
                force=force,
            )
        except Exception as exc:
            logger.error(
                "[mmgbsa.pipeline] action=skip reason=receptor_prep_failed pdb=%s variant=%s ph=%s err=%s",
                pdb_id,
                variant_dir,
                ph_label,
                exc,
            )
            if strict:
                raise
            continue

        receptor_for_leap_path, strip_info = ppm._strip_receptor_h_for_leap(
            Path(receptor_result["output_path"]),
            strip_all_h_for_leap,
        )
        hoh_residue_count = strip_info.get("hoh_residue_count", 0)
        ppm._log_leap_prep(
            mmgbsa_logger,
            "INFO",
            {
                "strip_all_h": strip_all_h_for_leap,
                "input": receptor_result["output_path"],
                "output": str(receptor_for_leap_path),
                "removed_H": strip_info.get("removed_h", 0),
                "hoh_residues": hoh_residue_count,
            },
        )
        ppm._log_leap_prep(
            mmgbsa_logger,
            "INFO",
            {"map_hoh_to_wat": map_hoh_to_wat, "water_model": water_model, "hoh_residues": hoh_residue_count},
            "water_mapping",
        )

        mmgbsa_dir = ph_dir / "mmgbsa"
        summary_path = mmgbsa_dir / "mmgbsa_results_summary.csv"
        work_root = mmgbsa_dir / "work"

        prep_results = prep_mmgbsa_from_sdfs(
            [str(p) for p in sdfs],
            cfg=cfg,
            max_ligands=0,
            force=force,
            amber_prefix=amber_prefix,
        )

        ligand_entries: List[dict] = []
        for res in prep_results:
            if res.get("error"):
                _mmgbsa_append_summary(
                    summary_path,
                    {
                        "stage_dir": res.get("stage_dir", stage_dir_label),
                        "ligand_stem": Path(res.get("sdf_path", "ligand")).stem,
                        "delta_total": "",
                        "results_csv": "",
                        "ok": False,
                        "notes": res.get("error", "prep_failed"),
                    },
                )
                continue
            ligand_entries.append(
                {
                    "stage_dir": res.get("stage_dir", stage_dir_label),
                    "ligand_base": Path(res.get("sdf_path", "ligand")).stem,
                    "mol2_path": res.get("mol2_path"),
                    "frcmod_path": res.get("frcmod_path"),
                }
            )

        if not ligand_entries:
            logger.info(
                "[mmgbsa.pipeline] action=skip reason=ligand_prep_failed pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_dir,
                ph_label,
            )
            continue

        topo_results = write_leap_for_ligands(
            receptor_pdb_path=str(receptor_for_leap_path),
            ligand_mol2_frcmod_pairs=ligand_entries,
            out_dir_base=str(work_root),
            force=force or _to_bool(cfg.get("MMGBSA_TLEAP_FORCE", False)),
            water_model=water_model,
            map_hoh_to_wat=map_hoh_to_wat,
            hoh_residue_count=hoh_residue_count,
        )

        run_tleap_flag = _to_bool(cfg.get("MMGBSA_TLEAP_ENABLED", True)) and _to_bool(
            cfg.get("MMGBSA_TLEAP_RUN", True)
        )

        for topo in topo_results:
            ligand_stem = topo.get("ligand_base")
            stage_dir = topo.get("stage_dir")
            out_dir = Path(topo.get("output_dir", work_root))
            notes = ""
            ok = False
            results_csv = ""
            delta_total = ""
            summary_frames: Dict[str, object] = {
                "frame_agg_method": frame_agg_method,
                "frames_used": "",
                "frame_mean": "",
                "frame_sd": "",
                "frames_total": "",
                "frame_median": "",
            }
            replicates_total = md_reps if md_enabled else 1
            replicates_ok = 0
            rep_agg_method_used = md_rep_agg if md_enabled else "single"

            ppm._log_leap(
                mmgbsa_logger,
                "INFO",
                {
                    "stage_dir": stage_dir,
                    "ligand": ligand_stem,
                    "receptor_pdb_used": str(receptor_for_leap_path),
                    "build_leap_path": topo.get("leap_file"),
                    "tleap_log_path": out_dir / "tleap.log",
                    "tleap_run": run_tleap_flag,
                    "skip": topo.get("skip_tleap"),
                },
                "ligand_topology",
            )

            try:
                if run_tleap_flag and not topo.get("skip_tleap"):
                    run_tleap(topo["leap_file"], str(out_dir))

                required = {
                    "complex_prmtop": topo.get("complex_prmtop"),
                    "complex_inpcrd": topo.get("complex_inpcrd"),
                    "receptor_prmtop": topo.get("receptor_prmtop"),
                    "ligand_prmtop": topo.get("ligand_prmtop"),
                }
                missing = []
                for key, value in required.items():
                    if not value:
                        missing.append(key)
                        continue
                    path = Path(value)
                    if not path.exists() or path.stat().st_size == 0:
                        missing.append(key)

                if missing:
                    notes = f"missing_topology:{','.join(missing)}"
                    _mmgbsa_append_summary(
                        summary_path,
                        {
                            "stage_dir": stage_dir,
                            "ligand_stem": ligand_stem,
                            "delta_total": "",
                            "results_csv": "",
                            "ok": False,
                            "notes": notes,
                        },
                    )
                    continue

                if five_reps:
                    rep_result = _mmgbsa_five_replicate_runner(
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
                    ok = bool(rep_result.get("ok"))
                    delta_val = rep_result.get("delta_total")
                    if delta_val is not None:
                        delta_total = f"{float(delta_val):.6g}"
                    results_csv = rep_result.get("results_csv", "")
                    notes = rep_result.get("notes", "")
                    summary_frames.update(
                        {
                            "frame_mean": rep_result.get("frame_mean", ""),
                            "frame_sd": rep_result.get("frame_sd", ""),
                        }
                    )
                    replicates_ok = int(rep_result.get("replicates_ok", 0))
                    replicates_total = int(rep_result.get("replicates_total", 5))
                    _mmgbsa_append_summary(
                        summary_path,
                        {
                            "stage_dir": stage_dir,
                            "ligand_stem": ligand_stem,
                            "delta_total": delta_total,
                            "results_csv": results_csv,
                            "ok": ok,
                            "notes": notes or "",
                            "frame_agg_method": summary_frames.get("frame_agg_method", "") or "",
                            "frames_used": summary_frames.get("frames_used", "") or "",
                            "frames_total": summary_frames.get("frames_total", "") or "",
                            "frame_mean": summary_frames.get("frame_mean", "") or "",
                            "frame_sd": summary_frames.get("frame_sd", "") or "",
                            "frame_median": summary_frames.get("frame_median", "") or "",
                            "rep_agg_method": rep_agg_method_used,
                            "replicates_ok": replicates_ok,
                            "replicates_total": replicates_total,
                        },
                    )
                    continue

                if not md_enabled:
                    traj_result = make_mmgbsa_trajectory(
                        complex_prmtop=topo["complex_prmtop"],
                        complex_inpcrd=topo["complex_inpcrd"],
                        out_dir=str(out_dir),
                        cfg=cfg,
                        force=force,
                        run_cpptraj=_to_bool(cfg.get("MMGBSA_CPPTRAJ_RUN", True)),
                    )
                    traj_path = Path(traj_result.get("trajout_path") or "")
                    if not traj_path.exists() or traj_path.stat().st_size == 0:
                        default_traj = str(cfg.get("MMGBSA_DEFAULT_TRAJ_NAME", "mdcrd") or "mdcrd")
                        candidate = out_dir / default_traj
                        if candidate.exists() and candidate.stat().st_size > 0:
                            traj_path = candidate
                        else:
                            notes = "missing_trajectory"
                            _mmgbsa_append_summary(
                                summary_path,
                                {
                                    "stage_dir": stage_dir,
                                    "ligand_stem": ligand_stem,
                                    "delta_total": "",
                                    "results_csv": "",
                                    "ok": False,
                                    "notes": notes,
                                },
                            )
                            continue

                    mmpbsa_result = run_mmgbsa(
                        complex_prmtop=topo["complex_prmtop"],
                        receptor_prmtop=topo["receptor_prmtop"],
                        ligand_prmtop=topo["ligand_prmtop"],
                        trajectory_path=str(traj_path),
                        work_dir=str(out_dir),
                        cfg=cfg,
                        force=force,
                        run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
                    )

                    results_csv = mmpbsa_result.get("out_csv", "")
                    results_dat = mmpbsa_result.get("out_dat", "")
                    if not mmpbsa_result.get("enabled", True):
                        notes = "mmpbsa_disabled"
                    else:
                        if results_csv and results_dat and Path(results_csv).exists() and Path(results_dat).exists():
                            frame_meta = _mmgbsa_frame_aggregate(Path(results_csv), cfg, logger)
                            delta_val = frame_meta.get("delta_total")
                            if delta_val is not None:
                                delta_total = f"{float(delta_val):.6g}"
                                ok = True
                            else:
                                notes = frame_meta.get("notes", "missing_delta") or "missing_delta"
                            summary_frames.update(
                                {
                                    "frame_agg_method": frame_meta.get("frame_agg_method", ""),
                                    "frames_used": frame_meta.get("frames_used", ""),
                                    "frame_mean": frame_meta.get("frame_mean", ""),
                                    "frame_sd": frame_meta.get("frame_sd", ""),
                                    "frames_total": frame_meta.get("frames_total", ""),
                                    "frame_median": frame_meta.get("frame_median", ""),
                                }
                            )
                            replicates_ok = 1 if ok else 0
                        else:
                            notes = "missing_outputs"
                else:
                    rep_seeds = _mmgbsa_md_seed_list(cfg, md_reps, logger)
                    rep_results: List[Dict[str, object]] = []
                    for rep_idx, seed in enumerate(rep_seeds, start=1):
                        rep_dir = out_dir / f"rep{rep_idx}"
                        md_result = run_implicit_md(
                            complex_prmtop=topo["complex_prmtop"],
                            complex_inpcrd=topo["complex_inpcrd"],
                            out_dir=str(out_dir),
                            cfg=cfg,
                            replicate_index=rep_idx,
                            seed=int(seed),
                            force=force,
                            run=True,
                        )
                        traj_path = Path(md_result.get("traj_path") or "")
                        rep_notes = ""
                        rep_ok = False
                        rep_delta = None
                        rep_frames: Dict[str, object] = {
                            "frame_agg_method": frame_agg_method,
                            "frames_used": "",
                            "frames_total": "",
                            "frame_mean": "",
                            "frame_sd": "",
                            "frame_median": "",
                        }
                        rep_csv = ""
                        rep_dat = ""
                        rep_log = ""

                        if not md_result.get("ok"):
                            rep_notes = "md_failed"
                        elif not traj_path.exists() or traj_path.stat().st_size == 0:
                            rep_notes = "missing_trajectory"
                        else:
                            try:
                                mmpbsa_result = run_mmgbsa(
                                    complex_prmtop=topo["complex_prmtop"],
                                    receptor_prmtop=topo["receptor_prmtop"],
                                    ligand_prmtop=topo["ligand_prmtop"],
                                    trajectory_path=str(traj_path),
                                    work_dir=str(rep_dir),
                                    cfg=cfg,
                                    force=force,
                                    run=_to_bool(cfg.get("MMGBSA_MMPBSA_RUN", True)),
                                )
                                rep_csv = mmpbsa_result.get("out_csv", "")
                                rep_dat = mmpbsa_result.get("out_dat", "")
                                rep_log = mmpbsa_result.get("log_path", "")
                                if not mmpbsa_result.get("enabled", True):
                                    rep_notes = "mmpbsa_disabled"
                                elif rep_csv and rep_dat and Path(rep_csv).exists() and Path(rep_dat).exists():
                                    frame_meta = _mmgbsa_frame_aggregate(Path(rep_csv), cfg, logger)
                                    delta_val = frame_meta.get("delta_total")
                                    rep_frames.update(
                                        {
                                            "frame_agg_method": frame_meta.get("frame_agg_method", frame_agg_method),
                                            "frames_used": frame_meta.get("frames_used", ""),
                                            "frames_total": frame_meta.get("frames_total", ""),
                                            "frame_mean": frame_meta.get("frame_mean", ""),
                                            "frame_sd": frame_meta.get("frame_sd", ""),
                                            "frame_median": frame_meta.get("frame_median", ""),
                                        }
                                    )
                                    if delta_val is not None:
                                        rep_delta = float(delta_val)
                                        rep_ok = True
                                    else:
                                        rep_notes = frame_meta.get("notes", "missing_delta") or "missing_delta"
                                else:
                                    rep_notes = "missing_outputs"
                            except Exception as exc:
                                rep_notes = f"mmpbsa_error:{type(exc).__name__}"
                                logger.warning(
                                    "[mmgbsa.pipeline] md_replicate_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s rep=%d err=%s",
                                    pdb_id,
                                    variant_dir,
                                    ph_label,
                                    stage_dir,
                                    ligand_stem,
                                    rep_idx,
                                    exc,
                                )
                                if strict:
                                    raise

                        rep_results.append(
                            {
                                "replicate": rep_idx,
                                "seed": seed,
                                "ok": rep_ok,
                                "delta_total": rep_delta,
                                "results_csv": rep_csv,
                                "results_dat": rep_dat,
                                "log_path": rep_log,
                                "traj_path": str(traj_path),
                                "work_dir": str(rep_dir),
                                "notes": rep_notes,
                                "frames_used": rep_frames.get("frames_used", ""),
                                "frames_total": rep_frames.get("frames_total", ""),
                                "frame_mean": rep_frames.get("frame_mean", ""),
                                "frame_sd": rep_frames.get("frame_sd", ""),
                                "frame_agg_method": rep_frames.get("frame_agg_method", frame_agg_method),
                                "frame_median": rep_frames.get("frame_median", ""),
                            }
                        )

                    ok_reps = [r for r in rep_results if r.get("ok")]
                    n_ok = len(ok_reps)
                    replicates_ok = n_ok
                    agg_delta_val = _mmgbsa_md_aggregate(
                        [float(r["delta_total"]) for r in ok_reps if r.get("delta_total") is not None],
                        md_rep_agg,
                    )
                    if agg_delta_val is not None:
                        delta_total = f"{agg_delta_val:.6g}"
                        ok = True
                    else:
                        ok = False

                    best_rep = None
                    if ok_reps:
                        if md_rep_agg == "min":
                            best_rep = min(ok_reps, key=lambda r: (r.get("delta_total", 0), r.get("replicate", 0)))
                        elif md_rep_agg == "median":
                            best_rep = min(
                                ok_reps,
                                key=lambda r: (abs(float(r.get("delta_total", 0)) - float(agg_delta_val or 0)), r.get("replicate", 0)),
                            )
                        else:
                            best_rep = min(
                                ok_reps,
                                key=lambda r: (abs(float(r.get("delta_total", 0)) - float(agg_delta_val or 0)), r.get("replicate", 0)),
                            )

                    if best_rep:
                        results_csv = str(best_rep.get("results_csv", ""))
                        summary_frames.update(
                            {
                                "frame_agg_method": best_rep.get("frame_agg_method", frame_agg_method),
                                "frames_used": best_rep.get("frames_used", ""),
                                "frame_mean": best_rep.get("frame_mean", ""),
                                "frame_sd": best_rep.get("frame_sd", ""),
                                "frames_total": best_rep.get("frames_total", ""),
                                "frame_median": best_rep.get("frame_median", ""),
                            }
                        )
                    notes = f"md_reps_ok={n_ok}/{len(rep_results)} rep_agg={md_rep_agg} frame_agg={frame_agg_method}"

                    replicate_summary = {
                        "ok": ok,
                        "n_reps_total": len(rep_results),
                        "n_reps_ok": n_ok,
                        "agg_method": md_rep_agg,
                        "agg_delta": agg_delta_val,
                        "best_replicate": best_rep,
                        "replicates": rep_results,
                    }
                    summary_json = out_dir / "mmgbsa_replicate_summary.json"
                    try:
                        tmp_json = summary_json.with_suffix(".json.part")
                        with tmp_json.open("w", encoding="utf-8") as handle:
                            json.dump(replicate_summary, handle, indent=2, sort_keys=True)
                        os.replace(tmp_json, summary_json)
                    except Exception as exc:
                        logger.warning(
                            "[mmgbsa.pipeline] md_replicate_summary_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s err=%s",
                            pdb_id,
                            variant_dir,
                            ph_label,
                            stage_dir,
                            ligand_stem,
                            exc,
                        )

                    if md_copy_best and best_rep:
                        for filename in ("FINAL_RESULTS_MMPBSA.dat", "FINAL_RESULTS_MMPBSA.csv", "mmpbsa.log", "mmpbsa.in"):
                            src = Path(best_rep.get("work_dir", "")) / filename
                            dest = out_dir / filename
                            if not src.exists():
                                continue
                            if dest.exists() and not force:
                                continue
                            try:
                                shutil.copy2(src, dest)
                            except Exception as exc:
                                logger.warning(
                                    "[mmgbsa.pipeline] md_copy_best_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s file=%s err=%s",
                                    pdb_id,
                                    variant_dir,
                                    ph_label,
                                    stage_dir,
                                    ligand_stem,
                                    filename,
                                    exc,
                                )

            except Exception as exc:
                notes = f"error:{type(exc).__name__}"
                ok = False
                logger.error(
                    "[mmgbsa.pipeline] ligand_failed pdb=%s variant=%s ph=%s stage=%s ligand=%s err=%s",
                    pdb_id,
                    variant_dir,
                    ph_label,
                    stage_dir,
                    ligand_stem,
                    exc,
                )
                if strict:
                    raise

            _mmgbsa_append_summary(
                summary_path,
                {
                    "stage_dir": stage_dir,
                    "ligand_stem": ligand_stem,
                    "delta_total": delta_total,
                    "results_csv": results_csv,
                    "ok": ok,
                    "notes": notes,
                    "frame_agg_method": summary_frames.get("frame_agg_method", "") or "",
                    "frames_used": summary_frames.get("frames_used", "") or "",
                    "frames_total": summary_frames.get("frames_total", "") or "",
                    "frame_mean": summary_frames.get("frame_mean", "") or "",
                    "frame_sd": summary_frames.get("frame_sd", "") or "",
                    "frame_median": summary_frames.get("frame_median", "") or "",
                    "rep_agg_method": rep_agg_method_used,
                    "replicates_ok": replicates_ok,
                    "replicates_total": replicates_total,
                },
            )

        if agg_enabled:
            agg_path = mmgbsa_dir / agg_output_name
            _mmgbsa_write_pose_aggregate(
                summary_path=summary_path,
                out_path=agg_path,
                stage_dir_label=stage_dir_label,
                pose_group_regexes=pose_group_regexes,
                agg_method=agg_method,
                logger=logger,
            )

        logger.info(
            "[mmgbsa.pipeline] done pdb=%s variant=%s ph=%s summary=%s",
            pdb_id,
            variant_dir,
            ph_label,
            summary_path,
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
            logging.exception("[notify] unexpected error while building notification email")
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
