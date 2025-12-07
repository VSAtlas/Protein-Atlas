# -*- coding: utf-8 -*-
# High-level pipeline for multi-stage docking.
#
# Phases per protein:
# 1) Setup & logging
# 2) Extract ligands and generate a ligand-free PDB
# 3) Protein preparation (clean PDB + receptor PDBQT), reuse if cached
# 4) Active-site detection (center, box size)
# 5) Ligand preparation & filtering
# 6) Multi-stage docking with early/fallback recenter heuristics + CenterSelector
# 7) Final pose validation & optional screenshots
# 8) Write per-protein score CSV

from __future__ import annotations
from activesite import extract_and_remove_ligands, get_atom_rules
import sys, hashlib, re, logging, json, time, os, shutil, re
from dataclasses import dataclass, field
import datetime
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, Any, Mapping
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
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
    update_manifest_for_protein_success,
    update_manifest_for_scheduled_proteins,
)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
from input_and_export_functions import (
    load_inputs, validate_config, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, emit_vina_config as _emit_vina_config_impl, _to_bool, init_config_run_dir
)
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
from path_router import (
    expand_variants,
    make_paths,
    Paths as RouterPaths,
    receptor_file,
    docked_dir,
    config_dir as router_config_dir,
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
from debug_fs import install_debug_makedirs, collapse_sanitized_names_for_cfg

# Install debug wrappers for Path.mkdir and os.makedirs at import time,
# preserving the previous behavior.
install_debug_makedirs()


def _resolve_run_id(argv: list[str]) -> str:
    cli_run_id = _cli_val(argv, "--run-id")
    env_run_id = (os.environ.get("ATLAS_RUN_ID") or "").strip()
    if cli_run_id:
        return cli_run_id
    if env_run_id:
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



def _cli_val(argv, flag):
    try:
        if flag in argv:
            i = argv.index(flag)
            if i + 1 < len(argv) and not argv[i+1].startswith("-"):
                return argv[i+1]
    except Exception:
        pass
    return None

def _cli_has(argv, flag):
    try:
        return flag in argv
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
        if old_variant is None:
            os.environ.pop("APO_HOLO_VARIANT", None)
        else:
            os.environ["APO_HOLO_VARIANT"] = old_variant


# ======================
# Program entry point
# ======================
def main() -> None:
    print("MODELLER is working with license.")
    run_id = _resolve_run_id(sys.argv)
    os.environ["ATLAS_RUN_ID"] = run_id
    log_path = _prepare_run_logfile(run_id)
    os.environ["ATLAS_LOG_FILE"] = log_path
    _tee_stdio_to(log_path)
    cfg = ConfigDict(load_inputs())
    bootstrap_root_logging(cfg, log_path)
    logging.info("[probe.root] root-logger INFO now visible")
    print(f"[run] log_file={log_path} run_id={run_id}")
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
    cfg.setdefault("PH_RADIUS", 10.0)

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
            cfg.PH_RADIUS = 10.0
    else:
        try:
            cfg.PH_RADIUS = float(cfg.PH_RADIUS)
        except Exception:
            cfg.PH_RADIUS = 10.0
    if cfg.PH_RADIUS <= 0:
        cfg.PH_RADIUS = 10.0

    log = logging.getLogger("ph_ensemble")
    log.info("[ph_ensemble.mode] enabled=%s scope=%s radius=%s", cfg.PH_ENSEMBLE, getattr(cfg, "PH_SCOPE", "auto"), getattr(cfg, "PH_RADIUS", 10.0))
    # --- PH-ligand mode (CLI > ENV > CFG) ---
    cfg.setdefault("PH_LIGAND_MODE", "off")
    cli_ph_mode = _cli_val(sys.argv, "--ph-ligand-mode")
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
    cli_cfg_dir = _cli_val(sys.argv, "--configs-dir")
    cli_no_reset = _cli_has(sys.argv, "--no-reset-configs")

    if cli_cfg_dir:  cfg["CONFIGS_DIR"] = cli_cfg_dir
    cfg["RUN_ID"] = run_id
    if cli_no_reset: cfg["RESET_CONFIGS"] = False

    init_config_run_dir(cfg, run_id=cfg.get("RUN_ID"), reset=cfg.get("RESET_CONFIGS"),
                        logger=logging.getLogger("run"))
    print(f"[cfg.run] run_id={cfg['RUN_ID']} run_dir={cfg['CONFIG_RUN_DIR']}")
    print(f"[ph.mode] PH_ENSEMBLE={cfg.get('PH_ENSEMBLE', False)}")

    global_start = time.time()

    # Best-effort manifest initialization
    try:
        init_run_manifest(cfg, run_id, sys.argv, log_path)
    except Exception:
        logging.warning("Failed to initialize run_manifest.yaml", exc_info=True)

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
    cli_single = _parse_single_from_cli(sys.argv)
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

    # Small FDA test library override.
    # When called with: python main.py -test -test-fda -fast
    # we want to use a tiny FDA test library instead of the full FDA set.
    #
    # By default this points at "fda_test_library_10", which should correspond to:
    #   prepped_ligands/fda_test_library_10/
    #   extracted_ligands/fda_test_library_10/
    cfg.setdefault("TEST_FDA_LIBRARY_SUBDIR", "fda_test_library_10")

    if _cli_has(sys.argv, "-test-fda") or _cli_has(sys.argv, "--test-fda"):
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
    requested_ids, _sel_src = _parse_specified_proteins(sys.argv, cfg)
    cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] = requested_ids
    print(f"[config] SPECIFIED_PROTEINS effective={requested_ids} (precedence: CLI>ENV>CFG)")
    # --- Fast mode: force exhaustiveness=1 everywhere ---
    cfg["FAST_MODE"] = _parse_fast_flag(sys.argv) or bool(cfg.get("FAST_MODE", False))
    if cfg["FAST_MODE"]:
        print("[config] FAST_MODE effective=True (exhaustiveness=1)")

    # --- No-library docking mode (controls-only + ligand planning) ------
    #
    # If enabled, we still run receptor prep + control redocking +
    # pocket detection + ligand selection/enumeration, but we skip the
    # actual docking calls for DUD/FDA libraries.
    #
    cfg.setdefault("NO_LIBRARY_DOCKING", False)

    cli_no_dock = _cli_has(sys.argv, "--no-docking")
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
            kept, skipped = [], []
            for f in pdb_files:
                nid = _norm_pdb_id(f)
                if nid and nid in test_keys:
                    kept.append(f)
                else:
                    skipped.append(f)

            if skipped:
                print(
                    f"[test-mode] Enabled mode={test_mode}; restricting to "
                    f"{len(kept)} PDBs from TEST_LIBRARY_MAP keys."
                )
                for s in skipped:
                    print(f"[test-mode] Skipping {s} (not in TEST_LIBRARY_MAP).")

            pdb_files = kept
        else:
            print("[test-mode] TEST_LIBRARY_MAP empty/invalid; no extra filtering applied.")
    else:
        cfg["_TEST_LIBRARY_CANONICAL"] = {}


    # Normalize any repeated '.sanitized' tokens in ligand filenames
    sanitize_logger = logging.getLogger("sanitize")
    collapse_sanitized_names_for_cfg(cfg, logger=sanitize_logger)

    print("Working directory:", os.getcwd())
    print("Loaded config keys:", list(cfg.keys()))
    print(f"Proteins queued: {len(pdb_files)}")
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
                f"  - {pdb_id} ({label}): {exc_type} — {exc_msg}\n"
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

        # Use the same `mail` CLI you already tested in your bash wrapper
        try:
            pipe = os.popen(f'mail -s "{subject}" mpg2352@utexas.edu', "w")
            try:
                pipe.write(body)
            finally:
                pipe.close()
        except Exception:
            # If mail fails, log it but never break the run
            try:
                logging.exception("[notify] failed to send mail notification")
            except Exception:
                # Logging itself should not be able to kill the run
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
        # KeyboardInterrupt and other errors → non-zero status
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
