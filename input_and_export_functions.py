import csv
import math
import os
import re
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path
from distutils.util import strtobool
from typing import Any, Dict, Optional

import pandas as pd
import sitecustomize  # noqa: F401  # ensure HOME is writable for micromamba/pytest sandboxes
from path_router import make_paths

# -------------------------
# OS guards (Windows-only)
# -------------------------
IS_WINDOWS = sys.platform.startswith("win")
try:
    if IS_WINDOWS:
        import win32api  # type: ignore
        import win32file  # type: ignore
    else:
        win32api = None
        win32file = None
except Exception:
    win32api = None
    win32file = None


# -------------------------
# Small helpers
# -------------------------
def _to_bool(x):
    try:
        return bool(strtobool(str(x)))
    except Exception:
        return False


def _to_int(x, default=None):
    try:
        return int(x)
    except Exception:
        return default


def _to_float(x, default=None):
    try:
        return float(x)
    except Exception:
        return default


def which_or_exists(candidates):
    """
    Return the first path that exists (if absolute),
    or the first found on PATH via shutil.which.
    """
    for c in candidates:
        p = Path(str(c))
        if p.is_absolute() and p.exists():
            return str(p)
        w = shutil.which(p.name)
        if w:
            return w
    return None


# -------------------------
# Defaults for BRCF layout
# -------------------------
def _default_base_dir() -> Path:
    # env wins; else directory of this file
    return Path(
        os.environ.get("PROTEIN_AUTOMATION_DIR", Path(__file__).resolve().parent)
    )


def _default_dirs(base: Path) -> Dict[str, str]:
    """
    Opinionated defaults that match BRCF layout:
      /stor/home/<user>/atlas/code/protein_automation/{subdirs}
    """
    return {
        "OVERALL_DIR": str(base),
        "INPUT_DIR": str(base / "input_pdbs"),
        "PROTEIN_DIR": str(base / "pdbqts"),
        "LIGAND_DIR": str(base / "input_ligands"),
        "LIGAND_EXTRACTED_DIR": str(base / "extracted_ligands"),
        "LIGANDS_MOL2_DIR": str(base / "ligands_mol2"),
        "OUTPUT_LIGANDS_DIR": str(base / "prepped_ligands"),
        "OUTPUT_DIR": str(base / "processed_pdbs"),
        "PDBQT_DIR": str(base / "pdbqts"),
        "DOCKED_DIR": str(base / "docked"),
        # p2rank
        "P2RANK_OUTPUT_DIR": str(base / "p2rank_out"),
        # cleanup script (kept in repo)
        "PHENIX_CLEAN_SCRIPT": str(base / "phenix_clean.py"),
    }


def _default_tools() -> Dict[str, str]:
    """
    Tool defaults prefer environment & PATH.
    """
    base = _default_base_dir()
    scorch_root = base.parent.parent / "tools" / "SCORCH"
    if (scorch_root / "scorch.py").exists():
        scorch_script_default = str(scorch_root / "scorch.py")
    else:
        scorch_script_default = which_or_exists(["scorch.py"])

    return {
        "VINA_PATH": which_or_exists(["vina"]),
        "VINA_EXE": which_or_exists(["vina"]),
        "GNINA_EXE": which_or_exists(["gnina"]),
        "OPENBABEL_PATH": which_or_exists(["obabel"]),
        "PYMOL_PATH": which_or_exists(["pymol"]),
        "REDUCE_EXE": which_or_exists(["reduce"]),
        # prefer env MGLTOOLS_* if set; otherwise try ADFRsuite pythonsh on PATH
        "MGLTOOLS_PYTHON": os.environ.get("MGL_PYTHON")
        or which_or_exists(["pythonsh"]),
        "PREPARE_LIGAND_SCRIPT": os.environ.get("PREPARE_LIGAND_SCRIPT")
        or which_or_exists(["prepare_ligand4.py"]),
        "PREPARE_RECEPTOR_SCRIPT": os.environ.get("PREPARE_RECEPTOR_SCRIPT")
        or which_or_exists(["prepare_receptor4.py"]),
        # p2rank: either absolute prank or found on PATH
        "P2RANK_PATH": shutil.which("prank") or "prank",
        # SCORCH: optional rescoring stage, default to sibling tools checkout
        "SCORCH_SCRIPT": os.environ.get("SCORCH_SCRIPT") or scorch_script_default,
        "SCORCH_ENV": os.environ.get("SCORCH_ENV") or "scorch-env",
    }


def _default_runtime() -> Dict[str, Any]:
    return {
        "CPU_ONLY": True,
        "CPU": os.cpu_count() or 8,
        "MAX_PARALLEL_JOBS": max(1, (os.cpu_count() or 8) // 2),
        "FORCE_REPROCESS": False,
        "DOCKING_MODE": "discovery",
        "QUIET_CONSOLE": False,
        # docking/recenter knobs
        "EARLY_RECENTER_RATIO": 0.70,
        "EARLY_RECENTER_MIN_EVAL": 10,
        "EARLY_RECENTER_FAR_A": 15.0,
        "EARLY_RECENTER_MEDIAN_A": 10.0,
        "ALLOW_BOX_EXPAND": True,
        "MAX_RECENTER_ATTEMPTS": 3,
        # pocket evaluation (optional)
        "POCKET_EVAL": False,
        "POCKET_EVAL_MAX_CALIBRATORS": 0,
        "POCKET_EVAL_SEED": 0,
        "POCKET_EVAL_FOLDS": 5,
        # control-centering knobs
        "CONTROL_CENTER_POLICY": "best_redock",
        "CONTROL_CENTER_CLOSE_MAX_A": 8.0,
        "CTRL_REDOCK_EXHAUSTIVENESS": 32,
        "CTRL_REDOCK_NMODES": 9,
        # --- benchmark policy knobs ---
        "BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY": True,
        "BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING": False,
        # --- MMGBSA receptor prep ---
        "MMGBSA_STRIP_METALS": True,
        "MMGBSA_WATER_POLICY": "ACTIVE_SITE",
        "MMGBSA_WATER_KEEP_RADIUS_A": 6.0,
        "MMGBSA_WATER_USE_ALIASES": True,
        "MMGBSA_METAL_USE_ALIASES": True,
        "MMGBSA_TOPOLOGY_PREP_ENABLED": True,
        "MMGBSA_TLEAP_RUN": True,
        "MMGBSA_TOPOLOGY_DIRNAME": "mmgbsa_topologies",
        "MMGBSA_AMBERTOOLS_PREFIX": "",
        "MMGBSA_CPPTRAJ_ENABLED": True,
        "MMGBSA_CPPTRAJ_RUN": True,
        "MMGBSA_TRAJOUT_NAME": "mdcrd",
        "MMGBSA_TRAJOUT_FORMAT": "mdcrd",
        "MMGBSA_TRAJIN_SOURCE": "INPCRD",
        "MMGBSA_TRAJIN_PATH": "",
        "MMGBSA_TRAJIN_FORMAT": "inpcrd",
        "MMGBSA_TRAJ_STARTFRAME": 1,
        "MMGBSA_TRAJ_ENDFRAME": 1,
        "MMGBSA_TRAJ_INTERVAL": 1,
        "MMGBSA_MMPBSA_ENABLED": True,
        "MMGBSA_MMPBSA_RUN": True,
        "MMGBSA_MMPBSA_STARTFRAME": 1,
        "MMGBSA_MMPBSA_ENDFRAME": 1,
        "MMGBSA_MMPBSA_INTERVAL": 1,
        "MMGBSA_MMPBSA_VERBOSE": 2,
        "MMGBSA_GB_IGB": 5,
        "MMGBSA_GB_SALTCON": 0.150,
        "MMGBSA_MMPBSA_INPUT_NAME": "mmpbsa.in",
        "MMGBSA_MMPBSA_LOG_NAME": "mmpbsa.log",
        "MMGBSA_MMPBSA_OUT_DAT": "FINAL_RESULTS_MMPBSA.dat",
        "MMGBSA_MMPBSA_OUT_CSV": "FINAL_RESULTS_MMPBSA.csv",
        "MMGBSA_DEFAULT_TRAJ_NAME": "mdcrd",
        "MMGBSA_ENABLED": True,
        "MMGBSA_KEEP_WATERS": True,
        "MMGBSA_WATER_KEEP_RADIUS": 6.0,
        "MMGBSA_KEEP_METALS": False,
        "MMGBSA_METAL_RETAIN_TOKENS": "",
        "MMGBSA_WATER_RETAIN_TOKENS": "",
        "MMGBSA_RECEPTOR_FORCE": False,
        "MMGBSA_FORCE": False,
        "MMGBSA_STRICT": False,
        "MMGBSA_LIGAND_AT": "gaff2",
        "MMGBSA_LIGAND_CHARGE_METHOD": "bcc",
        "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD": "bcc",
        "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD": "gas",
        "MMGBSA_LIGAND_NOMINAL_NET_CHARGE": 0,
        "MMGBSA_LIGAND_BCC_CHARGE_SWEEP": "-1,1",
        "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2": False,
        "MMGBSA_LIGAND_FORCE": False,
        "MMGBSA_RDKit_VALIDATE": True,
        "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD": 1,
        "MMGBSA_LIGAND_NET_CHARGE": 0,
        "MMGBSA_LIGAND_SQM_LEVEL": 2,
        "MMGBSA_INPUT_STAGE_DIR": "stage1",
        "MMGBSA_MAX_LIGANDS": 1,
        "MMGBSA_RERANKED_TOP_PCT": 0.0,
        "MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK": 6.0,
        "MMGBSA_TLEAP_ENABLED": True,
        "MMGBSA_TLEAP_FORCE": False,
        "MMGBSA_GENERAL_STARTFRAME": 1,
        "MMGBSA_GENERAL_ENDFRAME": 1,
        "MMGBSA_GENERAL_INTERVAL": 1,
        "MMGBSA_GENERAL_VERBOSE": 2,
    }


# -------------------------
# Config loading & validation
# -------------------------
_ALLOWED_ENV_OVERRIDES = {
    "VINA_EXE",
    "VINA_PATH",
    "GNINA_EXE",
    "OPENBABEL_PATH",
    "MGLTOOLS_PYTHON",
    "PREPARE_LIGAND_SCRIPT",
    "PREPARE_RECEPTOR_SCRIPT",
    "PYMOL_PATH",
    "P2RANK_PATH",
    "PHENIX_DIR",
    "PHENIX_LIB_PATH",
    "PHENIX_CLEAN_SCRIPT",
    "INPUT_DIR",
    "OUTPUT_DIR",
    "PDBQT_DIR",
    "DOCKED_DIR",
    "LIGAND_DIR",
    "LIGAND_EXTRACTED_DIR",
    "LIGANDS_MOL2_DIR",
    "OUTPUT_LIGANDS_DIR",
    "P2RANK_OUTPUT_DIR",
    "CPU",
    "CPU_ONLY",
    "MAX_PARALLEL_JOBS",
    "DOCKING_MODE",
    "REDUCE_EXE",
    "USE_MEEKO",
    "OVERALL_DIR",
    "PROTEIN_DIR",
    "DEEPCOY_DECOYS_PER_ACTIVE",
    "DEEPCOY_PYTHON",
    "DEEPCOY_CHUNK_SIZE",
    "DEEPCOY_BASE_SEED",
    "DEEPCOY_SEED_PER_CHUNK",
    "DEEPCOY_USE_ARGMAX_GENERATION",
    "DEEPCOY_TRY_DIFFERENT_STARTING",
    "DEEPCOY_NUM_DIFFERENT_STARTING",
    "DEEPCOY_NUM_SAMPLES",
    "DEEPCOY_KEEP_CHUNKS",
    "DEEPCOY_ACTIVE_SOURCES",
    "DEEPCOY_SOURCE_AUDIT",
    "DEEPCOY_SOURCE_AUDIT_MAX_LINES",
    "DEEPCOY_CHEMBL_MAX_PAGES",
    "DEEPCOY_ACTIVE_POTENCY_CUTOFF_NM",
    "DEEPCOY_ACTIVE_POTENCY_KEEP_UNKNOWN",
    # ---  allow ENV override for the knobs ---
    "BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY",
    "BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING",
    "CONTROL_CENTER_POLICY",
    "CONTROL_CENTER_CLOSE_MAX_A",
    "CTRL_REDOCK_EXHAUSTIVENESS",
    "CTRL_REDOCK_NMODES",
    # logging/topic gates (opt-in; safe to ignore if unset)
    "LOG_TOPICS",
    "LOG_LEVEL_FILE",
    "LOG_LEVEL_CONSOLE",
    # GNINA follow-up toggle
    "USE_GNINA",
    # LeDock follow-up toggle
    "USE_LEDOCK",
    # DOCK6 follow-up toggle
    "USE_DOCK6",
    # SCORCH rescoring knobs
    "SCORCH_SCRIPT",
    "SCORCH_ENV",
    # MMGBSA receptor prep knobs
    "MMGBSA_STRIP_METALS",
    "MMGBSA_WATER_POLICY",
    "MMGBSA_WATER_KEEP_RADIUS_A",
    "MMGBSA_WATER_USE_ALIASES",
    "MMGBSA_METAL_USE_ALIASES",
    "MMGBSA_TOPOLOGY_PREP_ENABLED",
    "MMGBSA_TLEAP_RUN",
    "MMGBSA_TOPOLOGY_DIRNAME",
    "MMGBSA_AMBERTOOLS_PREFIX",
    "MMGBSA_CPPTRAJ_ENABLED",
    "MMGBSA_CPPTRAJ_RUN",
    "MMGBSA_TRAJOUT_NAME",
    "MMGBSA_TRAJOUT_FORMAT",
    "MMGBSA_TRAJIN_SOURCE",
    "MMGBSA_TRAJIN_PATH",
    "MMGBSA_TRAJIN_FORMAT",
    "MMGBSA_TRAJ_STARTFRAME",
    "MMGBSA_TRAJ_ENDFRAME",
    "MMGBSA_TRAJ_INTERVAL",
    "MMGBSA_MMPBSA_ENABLED",
    "MMGBSA_MMPBSA_RUN",
    "MMGBSA_MMPBSA_STARTFRAME",
    "MMGBSA_MMPBSA_ENDFRAME",
    "MMGBSA_MMPBSA_INTERVAL",
    "MMGBSA_MMPBSA_VERBOSE",
    "MMGBSA_GB_IGB",
    "MMGBSA_GB_SALTCON",
    "MMGBSA_MMPBSA_INPUT_NAME",
    "MMGBSA_MMPBSA_LOG_NAME",
    "MMGBSA_MMPBSA_OUT_DAT",
    "MMGBSA_MMPBSA_OUT_CSV",
    "MMGBSA_DEFAULT_TRAJ_NAME",
    "MMGBSA_ENABLED",
    "MMGBSA_KEEP_WATERS",
    "MMGBSA_WATER_KEEP_RADIUS",
    "MMGBSA_KEEP_METALS",
    "MMGBSA_METAL_RETAIN_TOKENS",
    "MMGBSA_WATER_RETAIN_TOKENS",
    "MMGBSA_RECEPTOR_FORCE",
    "MMGBSA_FORCE",
    "MMGBSA_STRICT",
    "MMGBSA_LIGAND_AT",
    "MMGBSA_LIGAND_CHARGE_METHOD",
    "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD",
    "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD",
    "MMGBSA_LIGAND_NOMINAL_NET_CHARGE",
    "MMGBSA_LIGAND_BCC_CHARGE_SWEEP",
    "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2",
    "MMGBSA_LIGAND_FORCE",
    "MMGBSA_RDKit_VALIDATE",
    "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD",
    "MMGBSA_LIGAND_NET_CHARGE",
    "MMGBSA_LIGAND_SQM_LEVEL",
    "MMGBSA_INPUT_STAGE_DIR",
    "MMGBSA_MAX_LIGANDS",
    "MMGBSA_RERANKED_TOP_PCT",
    "MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK",
    "MMGBSA_TLEAP_ENABLED",
    "MMGBSA_TLEAP_FORCE",
    "MMGBSA_GENERAL_STARTFRAME",
    "MMGBSA_GENERAL_ENDFRAME",
    "MMGBSA_GENERAL_INTERVAL",
    "MMGBSA_GENERAL_VERBOSE",
}


def _extract_brace_block(text: str, start_idx: int, open_char="{", close_char="}"):
    """Return the brace-balanced substring starting at the first open_char after start_idx."""
    i = text.find(open_char, start_idx)
    if i == -1:
        return None
    level = 0
    j = i
    while j < len(text):
        c = text[j]
        if c == open_char:
            level += 1
        elif c == close_char:
            level -= 1
            if level == 0:
                return text[i : j + 1]
        j += 1
    return None  # unbalanced


def _strip_inline_comment(value: str) -> str:
    """Strip inline comments starting with #, except when inside quotes."""
    in_single = False
    in_double = False
    escaped = False
    for idx, ch in enumerate(value):
        if escaped:
            escaped = False
            continue
        if ch == "\\":
            escaped = True
            continue
        if ch == "'" and not in_double:
            in_single = not in_single
            continue
        if ch == '"' and not in_single:
            in_double = not in_double
            continue
        if ch == "#" and not in_single and not in_double:
            return value[:idx]
    return value


def _parse_kv_config(path: Path) -> Dict[str, str]:
    """
    Lightweight key=value parser; allows comments starting with '#'.
    """
    out = {}
    if not path.exists():
        return out
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                val = _strip_inline_comment(v).strip()
                out[k.strip().upper()] = val
    return out


def _expand_vars_in_value(val: str, cfg_now: Dict[str, Any]) -> str:
    """
    Expand simple config variables inside string values.

    Supported forms:
      - ${VAR_NAME}   (preferred, "professional" style)
      - {OVERALL_DIR} (legacy)
      - $OVERALL_DIR  (legacy)

    For ${VAR_NAME}, we look up VAR_NAME (uppercased) in the current cfg dict.
    We intentionally keep this modest, not a full shell-style expansion.
    """
    if not isinstance(val, str):
        return val

    # 1) Generic ${VAR_NAME} expansion using keys from cfg_now
    def _repl(match: re.Match) -> str:
        key = match.group(1) or ""
        key_upper = key.upper()
        if key_upper in cfg_now and cfg_now[key_upper] is not None:
            return str(cfg_now[key_upper])
        # If we don't know the key, leave the original text unchanged
        return match.group(0)

    # Replace all ${VAR_NAME} occurrences
    val = re.sub(r"\$\{([A-Za-z0-9_]+)\}", _repl, val)

    # 2) Backwards-compatible OVERALL_DIR shorthands
    od = str(cfg_now.get("OVERALL_DIR", ""))
    if od:
        val = val.replace("{OVERALL_DIR}", od)
        val = val.replace("$OVERALL_DIR", od)

    return val


def _expand_all_vars(cfg: Dict[str, Any]) -> Dict[str, Any]:
    # single pass is sufficient for our use: expand OVERALL_DIR in the other paths
    out = dict(cfg)
    for k, v in list(out.items()):
        if isinstance(v, str):
            out[k] = _expand_vars_in_value(v, out)
    return out


def load_config(
    config_path: str = "config.txt", base_dir: Path | None = None
) -> Dict[str, Any]:
    """
    Load configuration with precedence:
      1) defaults (BRCF-aware)
      2) config.txt (key=value, optional)
      3) environment variables (allowed set only)
    Then expand {OVERALL_DIR}/$OVERALL_DIR inside string values.
    """
    base = base_dir or _default_base_dir()

    cfg: Dict[str, Any] = {}
    cfg.update(_default_dirs(base))
    cfg.update(_default_tools())
    cfg.update(_default_runtime())

    # overlay config file
    file_cfg = _parse_kv_config(Path(config_path))

    # --- preserve multi-line TEST_LIBRARY_MAP block ---
    try:
        raw = file_cfg.get("TEST_LIBRARY_MAP")
        if isinstance(raw, str):
            s = raw.strip()
            # If it looks like a dict but was truncated to one line, re-extract the full brace block
            if s.startswith("{") and not s.endswith("}"):
                txt = Path(config_path).read_text(encoding="utf-8", errors="ignore")
                key_idx = txt.find("TEST_LIBRARY_MAP")
                if key_idx != -1:
                    block = _extract_brace_block(txt, key_idx, "{", "}")
                    if block:
                        file_cfg["TEST_LIBRARY_MAP"] = block

    except Exception:
        # Non-fatal: if anything goes wrong, keep the original single-line value
        pass

    cfg.update(file_cfg)
    # Preserve raw file-sourced config (before env overrides) for precedence checks.
    cfg["_FILE_CFG"] = dict(file_cfg)

    # overlay env vars (uppercased keys only)
    for k, v in os.environ.items():
        K = k.upper()
        if (K in cfg or K in _ALLOWED_ENV_OVERRIDES) and v:
            cfg[K] = v

    # type coercion (before expansion is fine)
    for k in [
        "CPU_ONLY",
        "FORCE_REPROCESS",
        "ALLOW_BOX_EXPAND",
        "QUIET_CONSOLE",
        "USE_MEEKO",
        "USE_DOCK6",
        "use_dock6",
        "DEEPCOY_SEED_PER_CHUNK",
        "DEEPCOY_KEEP_CHUNKS",
        "DEEPCOY_USE_ARGMAX_GENERATION",
        "DEEPCOY_TRY_DIFFERENT_STARTING",
        "DEEPCOY_SOURCE_AUDIT",
        "DEEPCOY_ACTIVE_POTENCY_KEEP_UNKNOWN",
        "MMGBSA_STRIP_METALS",
        "MMGBSA_WATER_USE_ALIASES",
        "MMGBSA_METAL_USE_ALIASES",
        "MMGBSA_TOPOLOGY_PREP_ENABLED",
        "MMGBSA_TLEAP_RUN",
        "MMGBSA_CPPTRAJ_ENABLED",
        "MMGBSA_CPPTRAJ_RUN",
        "MMGBSA_MMPBSA_ENABLED",
        "MMGBSA_MMPBSA_RUN",
        "MMGBSA_ENABLED",
        "MMGBSA_KEEP_WATERS",
        "MMGBSA_KEEP_METALS",
        "MMGBSA_RECEPTOR_FORCE",
        "MMGBSA_FORCE",
        "MMGBSA_STRICT",
        "MMGBSA_TLEAP_ENABLED",
        "MMGBSA_TLEAP_FORCE",
        "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2",
        "MMGBSA_LIGAND_FORCE",
        "MMGBSA_RDKit_VALIDATE",
        "POCKET_EVAL",
    ]:
        if k in cfg:
            cfg[k] = _to_bool(cfg[k])

    if ("USE_DOCK6" in cfg) or ("use_dock6" in cfg):
        dock6_flag = _to_bool(cfg.get("USE_DOCK6", cfg.get("use_dock6")))
        cfg["USE_DOCK6"] = dock6_flag
        cfg["use_dock6"] = dock6_flag

    for k in [
        "MAX_PARALLEL_JOBS",
        "CPU",
        "MAX_RECENTER_ATTEMPTS",
        "EARLY_RECENTER_MIN_EVAL",
        "CTRL_REDOCK_EXHAUSTIVENESS",
        "CTRL_REDOCK_NMODES",
        "DEEPCOY_DECOYS_PER_ACTIVE",
        "DEEPCOY_CHUNK_SIZE",
        "DEEPCOY_BASE_SEED",
        "DEEPCOY_NUM_DIFFERENT_STARTING",
        "DEEPCOY_NUM_SAMPLES",
        "DEEPCOY_SOURCE_AUDIT_MAX_LINES",
        "DEEPCOY_CHEMBL_MAX_PAGES",
        "DEEPCOY_ACTIVE_POTENCY_CUTOFF_NM",
        "MMGBSA_TRAJ_STARTFRAME",
        "MMGBSA_TRAJ_ENDFRAME",
        "MMGBSA_TRAJ_INTERVAL",
        "MMGBSA_MMPBSA_STARTFRAME",
        "MMGBSA_MMPBSA_ENDFRAME",
        "MMGBSA_MMPBSA_INTERVAL",
        "MMGBSA_MMPBSA_VERBOSE",
        "MMGBSA_GB_IGB",
        "MMGBSA_LIGAND_SQM_LEVEL",
        "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD",
        "MMGBSA_MAX_LIGANDS",
        "MMGBSA_RERANKED_TOP_PCT",
        "MMGBSA_GENERAL_STARTFRAME",
        "MMGBSA_GENERAL_ENDFRAME",
        "MMGBSA_GENERAL_INTERVAL",
        "MMGBSA_GENERAL_VERBOSE",
        "POCKET_EVAL_MAX_CALIBRATORS",
        "POCKET_EVAL_SEED",
        "POCKET_EVAL_FOLDS",
    ]:
        if k in cfg:
            cfg[k] = _to_int(cfg[k], cfg[k])

    for k in [
        "EARLY_RECENTER_RATIO",
        "EARLY_RECENTER_FAR_A",
        "EARLY_RECENTER_MEDIAN_A",
        "CONTROL_CENTER_CLOSE_MAX_A",
        "MMGBSA_WATER_KEEP_RADIUS_A",
        "MMGBSA_WATER_KEEP_RADIUS",
        "MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK",
        "MMGBSA_RERANKED_TOP_PCT",
        "MMGBSA_GB_SALTCON",
    ]:
        if k in cfg:
            cfg[k] = _to_float(cfg[k], cfg[k])

    # normalize mode
    cfg["DOCKING_MODE"] = str(cfg.get("DOCKING_MODE", "discovery")).lower()

    # now expand {OVERALL_DIR}/$OVERALL_DIR appearances
    cfg = _expand_all_vars(cfg)

    return cfg


def validate_config(cfg: Dict[str, Any]):
    """
    Require OVERALL_DIR and INPUT_DIR to exist.
    Auto-create typical output directories if missing.
    """
    required_keys = [
        "OUTPUT_DIR",
        "INPUT_DIR",
        "PDBQT_DIR",
        "DOCKED_DIR",
        "OVERALL_DIR",
        "MGLTOOLS_PYTHON",
        "PREPARE_RECEPTOR_SCRIPT",
        "VINA_EXE",
        "MAX_PARALLEL_JOBS",
        "PYMOL_PATH",
    ]
    missing = [k for k in required_keys if not cfg.get(k)]
    if missing:
        raise ValueError(f"Missing required config keys or values: {missing}")

    # Must exist
    for path_key in ["OVERALL_DIR", "INPUT_DIR"]:
        p = Path(cfg[path_key])
        if not p.exists():
            raise FileNotFoundError(f"Required path does not exist: {p}")

    # Create-if-missing for outputs/caches
    create_keys = [
        "OUTPUT_DIR",
        "PDBQT_DIR",
        "DOCKED_DIR",
        "OUTPUT_LIGANDS_DIR",
        "LIGAND_EXTRACTED_DIR",
        "LIGANDS_MOL2_DIR",
        "P2RANK_OUTPUT_DIR",
    ]
    for path_key in create_keys:
        p = Path(cfg[path_key])
        try:
            p.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise FileNotFoundError(
                f"Could not create directory for {path_key}: {p} ({e})"
            )


def init_config_run_dir(cfg, run_id=None, reset=None, logger=None):
    root = Path(cfg.get("CONFIGS_DIR", Path(cfg["OVERALL_DIR"]) / "configs"))
    rid = run_id or cfg.get("RUN_ID")
    if not rid:
        rid = time.strftime("%Y%m%d_%H%M%S")
    cfg["RUN_ID"] = rid
    run_dir = root / rid
    cfg["CONFIG_RUN_DIR"] = str(run_dir)

    reset = bool(cfg.get("RESET_CONFIGS", True)) if reset is None else bool(reset)
    removed = 0
    if reset and run_dir.exists():
        for p in run_dir.rglob("*"):
            removed += 1
        shutil.rmtree(run_dir, ignore_errors=False)
    run_dir.mkdir(parents=True, exist_ok=True)

    msg = f"[cfg.reset] run_dir={run_dir} removed={removed}"
    (logger.info(msg) if logger else print(msg))


# -------------------------
# P2Rank helpers (optional)
# -------------------------
def p2rank_out_dir(cfg: Dict[str, Any]) -> Path:
    """
    Central place to define where P2Rank writes results.
    Default: OVERALL_DIR/p2rank_out, but user may override P2RANK_OUTPUT_DIR.
    """
    out = Path(cfg.get("P2RANK_OUTPUT_DIR") or Path(cfg["OVERALL_DIR"]) / "p2rank_out")
    out.mkdir(parents=True, exist_ok=True)
    return out


# -------------------------
# Docking stages
# -------------------------
def define_docking_stages(mode="discovery"):
    if mode == "discovery":
        return [
            {"name": "stage1", "num_modes": 1, "energy_range": 2, "exhaustiveness": 2},
            {"name": "stage2", "num_modes": 3, "energy_range": 3, "exhaustiveness": 4},
            {"name": "stage3", "num_modes": 5, "energy_range": 4, "exhaustiveness": 6},
            {"name": "stage4", "num_modes": 9, "energy_range": 6, "exhaustiveness": 8},
            {
                "name": "stage5",
                "num_modes": 20,
                "energy_range": 9,
                "exhaustiveness": 20,
            },
        ]
    elif mode == "polypharmacology":
        return [
            {"name": "stage1", "num_modes": 3, "energy_range": 2, "exhaustiveness": 4},
            {
                "name": "stage2",
                "num_modes": 10,
                "energy_range": 6,
                "exhaustiveness": 12,
            },
            {
                "name": "stage3",
                "num_modes": 20,
                "energy_range": 9,
                "exhaustiveness": 24,
            },
        ]
    else:
        raise ValueError(f"Unknown docking mode: {mode}")


# -------------------------
# Score I/O
# -------------------------
def write_score_summary_to_csv(
    score_history, output_path="docking_score_summary.csv", run_id=None, variant=None
):
    ligand_stage_pattern = re.compile(r"^(.*?)(_stage\d+)?\.pdbqt$", re.IGNORECASE)
    ligand_scores = defaultdict(dict)

    for stage, stage_scores in score_history.items():
        for path, score in stage_scores.items():
            ligand_filename = os.path.basename(path)
            match = ligand_stage_pattern.match(ligand_filename)
            if match:
                ligand_core = f"{match.group(1)}.pdbqt"
                ligand_scores[ligand_core][stage] = score

    all_ligands = sorted(ligand_scores.keys())
    all_stages = sorted(score_history.keys())

    run_id_value = "" if run_id is None else str(run_id)
    variant_value = (variant or "").strip()
    include_variant = bool(variant is not None and variant_value)

    header = ["run_id"]
    if include_variant:
        header.append("variant")
    header.append("Ligand")
    header.extend(all_stages)
    rows = []
    for ligand in all_ligands:
        row = [run_id_value]
        if include_variant:
            row.append(variant_value)
        row.append(ligand)
        for stage in all_stages:
            score = ligand_scores[ligand].get(stage, "")
            if isinstance(score, (float, int)):
                row.append(f"{score:.2f}")
            else:
                row.append(str(score) if score else "")
        rows.append(row)

    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(rows)

    print(f"\nScore summary written to: {output_path}")


def extract_best_score(docked_pdbqt_path):
    best_score = None
    with open(docked_pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.strip().split()
                score = float(parts[3])
                if best_score is None or score < best_score:
                    best_score = score
    return best_score


def extract_gnina_scores(docked_pdbqt_path: str) -> Dict[str, Optional[float]]:
    """
    Parse GNINA REMARK lines from a docked PDBQT.

    Returns keys:
      - minimized_affinity_kcal
      - cnn_score
      - cnn_affinity_pK
    Values are float or None if not present/parseable.
    """
    metrics = {
        "minimized_affinity_kcal": None,
        "cnn_score": None,
        "cnn_affinity_pK": None,
    }
    try:
        with open(docked_pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.startswith("REMARK"):
                    continue
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                token = parts[1].lower()
                try:
                    val = float(parts[-1])
                except Exception:
                    val = None
                if token == "minimizedaffinity":
                    metrics["minimized_affinity_kcal"] = val
                elif token == "cnnscore":
                    metrics["cnn_score"] = val
                elif token == "cnnaffinity":
                    metrics["cnn_affinity_pK"] = val
    except Exception:
        pass
    return metrics


# -------------------------
# Vina config writer, uses shim(lazy yeah)
# -------------------------
def generate_config(
    output_dir,
    pdb_id,
    receptor_pdbqt,
    center,
    box_size,
    ligand_path,
    stage,
    stage_info,
    cpu_per_job,
    docked_dir=None,
):
    # Legacy shim: derive a minimal cfg for older callers
    from docking.docking_vina import emit_vina_config

    cfg = {
        "OVERALL_DIR": output_dir,
        "CONFIGS_DIR": os.path.join(output_dir, "configs"),
        "DOCKED_DIR": docked_dir or os.path.join(output_dir, "docked"),
        "RUN_ID": "legacy",
        "RESET_CONFIGS": False,
    }
    cfg["CONFIG_RUN_DIR"] = os.path.join(cfg["CONFIGS_DIR"], cfg["RUN_ID"])
    Path(cfg["CONFIG_RUN_DIR"]).mkdir(parents=True, exist_ok=True)
    if not str(ligand_path).lower().endswith(".pdbqt"):
        raise ValueError(f"Ligand is not a .pdbqt file: {ligand_path}")
    return emit_vina_config(
        cfg,
        pdb_id,
        receptor_pdbqt,
        center,
        box_size,
        ligand_path,
        stage,
        stage_info,
        cpu_per_job,
        logger=None,
        legacy=True,
    )


# -------------------------
# Score bookkeeping
# -------------------------
def record_score(score_history, stage_name, ligand, score: Any, valid, reason=None):
    def coerce_score(x):
        if isinstance(x, (int, float)):
            return float(x)
        if isinstance(x, str) and x.lower().endswith(".pdbqt"):
            try:
                return extract_best_score(x)
            except Exception:
                return None
        try:
            return float(x)
        except Exception:
            return None

    s = coerce_score(score)
    score_history[stage_name][ligand] = {
        "score": s,
        "valid": bool(valid),
        "reason": reason,
    }


def score_key(item):
    """item = (ligand, rec). Sort by numeric score; invalid or None go to bottom."""
    _, rec = item
    s = rec.get("score")
    if s is None:
        return math.inf
    try:
        s = float(s)
    except (TypeError, ValueError):
        return math.inf
    return s if rec.get("valid", False) else s + 1e-6


def annotate_fda_long_csv_with_t_scores_vs_decoys(
    cfg: Dict[str, Any],
    pdb_id: str,
    ph_label: Optional[str] = None,
    logger=None,
) -> Optional[str]:
    """
    Post-processing helper used by docking.py when TEST_MODE_ENABLE includes DUD
    and we are running the FDA subrun.

    It reads dud_docking_score_long.csv to compute mean/std of best decoy scores,
    then annotates docking_score_long.csv for FDA ligands with a t_vs_decoys column.

    Returns the path to the updated FDA long CSV, or None if skipped.
    """
    try:
        from dud_eval import (
            compute_decoy_stats_from_long_csv,
            guess_ligfile_col,
            guess_score_col,
        )
    except Exception as e:
        if logger:
            logger.warning("[t-score.skip] pdb_id=%s reason=import_error %s", pdb_id, e)
        return None

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    ph_token = (ph_label or "").strip() or None
    variant_root = Path(paths.docked_variant_root(var, ph_token))

    dud_csv = variant_root / "dud_docking_score_long.csv"
    fda_csv = variant_root / "docking_score_long.csv"

    if not dud_csv.exists() or not fda_csv.exists():
        if logger:
            logger.info(
                "[t-score.skip] pdb_id=%s ph=%s reason=missing_csv dud=%s fda=%s",
                pdb_id,
                ph_label or "base",
                str(dud_csv),
                str(fda_csv),
            )
        return None

    mu, sigma, n_decoys = compute_decoy_stats_from_long_csv(dud_csv)
    if (
        not n_decoys
        or not math.isfinite(mu)
        or not math.isfinite(sigma)
        or sigma == 0.0
    ):
        if logger:
            logger.info(
                "[t-score.skip] pdb_id=%s ph=%s reason=degenerate_stats n=%s mu=%s sigma=%s",
                pdb_id,
                ph_label or "base",
                n_decoys,
                mu,
                sigma,
            )
        return None

    df = pd.read_csv(fda_csv)

    lig_col = guess_ligfile_col(df, None)
    score_col = guess_score_col(df, None)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    # Best score per ligand
    best = df.groupby(lig_col, as_index=False).agg(best_score=(score_col, "min"))
    best["t_vs_decoys"] = (mu - best["best_score"]) / sigma
    t_map = dict(zip(best[lig_col], best["t_vs_decoys"]))

    df["t_vs_decoys"] = df[lig_col].map(t_map)

    df.to_csv(fda_csv, index=False)

    if logger:
        logger.info(
            "[t-score.ok] pdb_id=%s ph=%s n_decoys=%s mean=%.3f std=%.3f out=%s",
            pdb_id,
            ph_label or "base",
            n_decoys,
            mu,
            sigma,
            str(fda_csv),
        )

    return str(fda_csv)


# -------------------------
# Common path builder per protein
# -------------------------
# >>> BUILD_PATHS SHIM START
def build_paths_for_protein(cfg, base_id, pdb_file, variant=None, ph_token=None):
    p = make_paths(cfg, base_id=base_id, pdb_file=pdb_file)
    return {
        "pdb_id": base_id.upper(),
        "pdb_path": str(p.input_pdb_path),
        "nolig_pdb_path": str(p.nolig_pdb_path),
        "ligand_output_dir": str(p.ligand_output_dir),
        "ligands_mol2_dir": str(p.ligands_mol2_dir),
        "prepped_ligands_dir": str(p.prepped_ligands_dir),
        "cleaned_pdb_path": str(p.receptor_cleaned_pdb(variant)),
        "receptor_pdbqt_path": str(p.receptor_pdbqt(variant, ph_token=ph_token)),
    }


# >>> BUILD_PATHS SHIM END
# -------------------------
# Backward-compat shims
# -------------------------


def load_inputs():
    """Legacy name used by existing scripts. Now just calls load_config()."""
    return load_config()


def get_default_config(prompt: bool = False):
    """
    Legacy helper kept for compatibility.
    simply return the fully merged config (defaults + config.txt + env).
    """
    return load_config()
