# -*- coding: utf-8 -*-
"""
Automated protein preparation pipeline (Windows/WSL-friendly)

This module prepares a receptor structure from a raw PDB by performing:
  1) Alternate-conformation filtering (altLoc) with deterministic selection
  2) Removal of nonstandard residues (policy-aware keep/drop of cofactors/metals/waters)
  3) Element-column fixes to ensure PDB compliance (via activesite/YAML)
  4) Optional loop/residue completion via MODELLER
  5) Sanity checks for incomplete residues
  6) Phenix cleaning passes (auto-detect Linux Phenix, Windows Phenix, or skip—WSL-safe)
  7) Hydrogen sanity cleanup (CONECT- and geometry-based)
  8) Chain validation (ensure at least one CA-containing chain)
  9) Protonation via Reduce with automatic fallbacks (temp retry/OpenBabel)
 10) Final polish via phenix.pdbtools when available; otherwise degrade gracefully
 11) Receptor PDBQT preparation (Meeko preferred; ADT fallback)

Design notes
------------
• Single source of truth for element handling: **activesite** helpers (YAML-driven). No local element inference.
• WSL-aware external calls.
• Logging instead of prints; short, numbered steps for easier debugging.
• Never touch PDBQT with any PDB element-fixing logic.

Prerequisites (recommended)
---------------------------
• Phenix (Linux: phenix.pdbtools; Windows: phenix.pdbtools.bat, phenix.python.bat)
• Reduce (reduce.exe or reduce) — often bundled with Phenix
• MODELLER (licensed) — optional but recommended for loop filling
• Open Babel (obabel) — fallback for hydrogen addition
• Biopython — for PDB parsing
• MGLTools (for prepare_receptor4.py) and/or Meeko (mk_prepare_receptor)
"""
from __future__ import annotations

import json
import platform
from collections import Counter, defaultdict
from math import sqrt
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union
import os, sys, shutil, subprocess, logging, re
# >>> PATHS IMPORT START
from path_router import make_paths
# >>> PATHS IMPORT END
from propka_wire import pdb2pqr_protonate

from installation import load_config
from logger_setup import setup_logger
import activesite as _activesite_mod
from activesite import (
    fix_element_columns_in_file,
    scan_helium_counts,
    assert_no_helium_in_hydrogen_names,
    get_atom_rules,
    fix_pdb_elements,
    scan_helium_counts_with_hits,
    format_pdb_atom_debug,
    rules_version,
)
ALIASES = get_atom_rules()
RULES = ALIASES.__dict__ if hasattr(ALIASES, "__dict__") else dict(ALIASES)
_HE_POSTWRITE_VERBOSE = (os.environ.get("HELIUM_POSTWRITE_VERBOSE", "1") != "0")

def load_aliases():
    return get_atom_rules()
def _flatten_semicolons(items):
    out = []
    for item in items or []:
        # items might be lines like "A; B; C"
        parts = [p.strip().upper() for p in str(item).split(";") if p.strip()]
        out.extend(parts)
    return out

_RETAIN = set(_flatten_semicolons(RULES.get("retain_in_receptor_resnames", [])))
_WATER_NAMES = {w for w in _RETAIN if w in {"HOH","WAT","DOD","H2O","TIP","TIP3","SOL"}}
# --- Canonical element tokens (upper-cased 1�2 letter symbols) ---
def _canonize_two_letter(xs) -> set[str]:
    out = set()
    for line in (xs or []):
        for tok in str(line).split(";"):
            t = tok.strip()
            if t:
                out.add(t.upper())
    return out

_ES      = RULES.get("element_sets", {}) or {}
_ONE     = {str(s).upper() for s in (_ES.get("one_letter_elements") or [])}
_TWORAW  = _ES.get("two_letter_elements", []) or []
_TWO     = _canonize_two_letter(_TWORAW)
# All allowed element tokens
_ELEM_CANON = _ONE | _TWO
_MEEKO_DROP_IONS = set(_flatten_semicolons(RULES.get("meeko_drop_free_ions", []))) or {"NA", "K", "LI"}

_SALT_RESNAMES = {"NA", "K", "CL", "BR", "I"}
_METAL_RESNAMES = {"MG", "MN", "FE", "ZN", "CU", "CO", "NI", "CA"}
_IONS_CFG_CACHE: tuple[set[str], str] | None = None
_IONS_CFG_LOGGED = False


def _resolve_variant_token(cfg: Optional[dict] = None, override: Optional[str] = None) -> Optional[str]:
    token = (override or "").strip().upper() if override else ""
    if token in {"APO", "HOLO"}:
        return token
    env_token = (os.environ.get("APO_HOLO_VARIANT") or "").strip().upper()
    if env_token in {"APO", "HOLO"}:
        return env_token
    if cfg is not None:
        cfg_token = str(cfg.get("_CURRENT_VARIANT", "")).strip().upper()
        if cfg_token in {"APO", "HOLO"}:
            return cfg_token
    return None


def _normalize_ion_policy(cfg: Optional[dict]) -> str:
    raw = "by_variant"
    if cfg is not None:
        raw = str(cfg.get("ION_STRIP_POLICY", "by_variant")).strip().lower() or "by_variant"
    if raw not in {"by_variant", "always_strip", "never_strip"}:
        logging.warning("[ions.cfg] unsupported_policy=%s fallback=by_variant", raw)
        return "by_variant"
    return raw


def _load_retain_allowlist(cfg: Optional[dict]) -> tuple[set[str], str]:
    global _IONS_CFG_CACHE, _IONS_CFG_LOGGED
    if _IONS_CFG_CACHE is not None:
        allow, source = _IONS_CFG_CACHE
        if not _IONS_CFG_LOGGED:
            logging.info("[ions.cfg] retain_in_receptor_resnames=%s source=%s", ",".join(sorted(allow)), source)
            _IONS_CFG_LOGGED = True
        return allow, source

    candidate = None
    if cfg and "retain_in_receptor_resnames" in cfg:
        candidate = cfg.get("retain_in_receptor_resnames")
    elif "retain_in_receptor_resnames" in config:
        candidate = config.get("retain_in_receptor_resnames")

    allow_items: Iterable[str] | None = None
    source = "default"

    if isinstance(candidate, (list, tuple, set)):
        allow_items = list(candidate)
        source = "inline"
    elif isinstance(candidate, str) and candidate.strip():
        text = candidate.strip()
        parsed: Iterable[str] | None = None
        if text.startswith("[") and text.endswith("]"):
            try:
                loaded = json.loads(text)
                if isinstance(loaded, list):
                    parsed = loaded
                    source = "inline"
            except Exception as exc:
                logging.warning("[ions.cfg] inline_json_parse_failed=%s err=%s", text[:40], exc)
        if parsed is None:
            yaml_path = Path(text).expanduser()
            if yaml_path.is_file():
                try:
                    with open(yaml_path, "r", encoding="utf-8") as fh:
                        data = _activesite_mod.yaml.safe_load(fh) or []
                    if isinstance(data, dict):
                        payload = data.get("retain_in_receptor_resnames")
                        if isinstance(payload, list):
                            data = payload
                    if isinstance(data, list):
                        parsed = data
                        source = str(yaml_path)
                    else:
                        raise TypeError("yaml_payload_not_list")
                except Exception as exc:
                    logging.warning("[ions.cfg] yaml_load_failed path=%s err=%s", yaml_path, exc)
            if parsed is None:
                tokens = [tok.strip() for tok in text.replace(";", ",").split(",") if tok.strip()]
                if tokens:
                    parsed = tokens
                    source = "inline"
        allow_items = parsed

    if allow_items is None:
        allow_items = RULES.get("retain_in_receptor_resnames", [])
        source = "default"

    allow_set = {tok.upper() for tok in _flatten_semicolons(allow_items)}
    _IONS_CFG_CACHE = (allow_set, source)
    if not _IONS_CFG_LOGGED:
        logging.info(
            "[ions.cfg] retain_in_receptor_resnames=%s source=%s",
            ",".join(sorted(allow_set)),
            source,
        )
        _IONS_CFG_LOGGED = True
    return allow_set, source


def _bucket_counts(counter: Counter) -> dict[str, int]:
    keys = ["ZN", "MG", "NA", "K", "CA", "MN", "FE", "CL"]
    out = {k: 0 for k in keys}
    other = 0
    for resn, count in counter.items():
        token = resn.upper()
        if token in out:
            out[token] += count
        else:
            other += count
    out["OTHER"] = other
    return out


def _format_counts(counter: Counter) -> str:
    bucketed = _bucket_counts(counter)
    keys = ["ZN", "MG", "NA", "K", "CA", "MN", "FE", "CL", "OTHER"]
    parts = [f"{k}={bucketed.get(k, 0)}" for k in keys]
    return " ".join(parts)


def _distance_from_center(line: str, center: Optional[tuple[float, float, float]]) -> float | None:
    if center is None:
        return None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except Exception:
        return None
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]
    return sqrt(dx * dx + dy * dy + dz * dz)


def _maybe_strip_ions(
    pdb_path: Union[str, Path],
    cfg: Optional[dict] = None,
    *,
    variant: Optional[str] = None,
    pocket_center: Optional[tuple[float, float, float]] = None,
    extra_keep: Optional[Iterable[str]] = None,
) -> int:
    path = Path(pdb_path)
    if not path.exists():
        logging.warning("[ions] skip_missing file=%s", path)
        return 0

    allow_tokens, _ = _load_retain_allowlist(cfg)
    allow_set = {str(tok).strip().upper() for tok in allow_tokens if str(tok).strip()}
    extra_set: set[str] = set()
    if extra_keep:
        for token in extra_keep:
            if token:
                extra_set.add(str(token).strip().upper())
    holo_keep_tokens = allow_set | extra_set
    apo_keep_tokens = set(extra_set)

    policy = _normalize_ion_policy(cfg)
    variant_token = _resolve_variant_token(cfg, variant)
    variant_label = variant_token or "legacy"

    radius_cfg = 6.0
    if cfg is not None:
        try:
            radius_cfg = float(cfg.get("HOLO_SALT_STRIP_RADIUS", 6.0) or 0.0)
        except Exception:
            radius_cfg = 6.0
    radius_cfg = max(0.0, radius_cfg)
    radius = radius_cfg if (policy == "by_variant" and variant_token == "HOLO") else 0.0
    radius_term = f"{radius:.2f}" if radius > 0.0 else "none"

    # [ions] stage=clean instrumentation
    logging.info(
        "[ions.policy] stage=clean variant=%s policy=%s salts_radius=%s allowlist=%d file=%s",
        variant_label,
        policy,
        radius_term,
        len(allow_set),
        path,
    )

    try:
        text = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as exc:
        logging.warning("[ions] read_failed file=%s err=%s", path, exc)
        return 0

    residues: dict[tuple[str, str, str, str], list[tuple[int, str]]] = {}
    for idx, line in enumerate(text):
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        key = (line[21], line[22:26], line[26], resname)
        residues.setdefault(key, []).append((idx, line))

    before = Counter()
    stripped = Counter()
    kept_counter = Counter()
    category_totals: dict[str, dict[str, int]] = {
        "metal": {"kept": 0, "stripped": 0},
        "salt": {"kept": 0, "stripped": 0},
        "other": {"kept": 0, "stripped": 0},
    }

    warn_missing_center = False
    remove_indices: set[int] = set()

    for key, atoms in residues.items():
        if len(atoms) != 1:
            continue
        chain, resi, icode, resname = key
        res_token = resname.upper()
        line_idx, line = atoms[0]
        before[res_token] += 1

        remove = False
        reason = ""
        elem = (line[76:78].strip() or res_token).upper()

        is_metal = res_token in _METAL_RESNAMES or elem in _METAL_RESNAMES
        is_salt = res_token in _SALT_RESNAMES or elem in _SALT_RESNAMES
        category = "metal" if is_metal else "salt" if is_salt else "other"

        if policy == "never_strip":
            remove = False
        elif policy == "always_strip":
            if res_token in holo_keep_tokens:
                remove = False
            else:
                remove = True
                reason = "policy_always_strip"
        else:  # policy == by_variant
            if variant_token == "HOLO":
                if res_token in holo_keep_tokens or is_metal:
                    remove = False
                elif is_salt:
                    if radius > 0.0:
                        dist = _distance_from_center(line, pocket_center)
                        if dist is None:
                            warn_missing_center = True
                            remove = False
                        elif dist >= radius:
                            remove = True
                            reason = "salt_far"
                        else:
                            remove = False
                    else:
                        remove = False
                else:
                    remove = False
            else:
                if res_token in apo_keep_tokens:
                    remove = False
                    reason = "apo_keep_override"
                else:
                    if is_metal:
                        remove = True
                        reason = "apo_strip_metal"
                    elif is_salt:
                        remove = True
                        reason = "apo_strip_salt"
                    else:
                        remove = True
                        reason = "apo_strip_other"

        if remove:
            stripped[res_token] += 1
            remove_indices.add(line_idx)
            logging.debug(
                "[ions.remove] elem=%s resname=%s serial=%s chain=%s resi=%s reason=%s",
                elem,
                res_token,
                line[6:11].strip(),
                chain.strip() or "-",
                (resi or "0").strip() or "0",
                reason,
            )
            category_totals[category]["stripped"] += 1
        else:
            kept_counter[res_token] += 1
            category_totals[category]["kept"] += 1

    logging.info(
        "[ions.counts.before] stage=clean variant=%s file=%s detail=%s",
        variant_label,
        path,
        _format_counts(before),
    )
    kept = before - stripped
    logging.info(
        "[ions.counts.after] stage=clean variant=%s file=%s detail=%s",
        variant_label,
        path,
        _format_counts(kept),
    )
    total_kept = sum(kept_counter.values())
    total_stripped = sum(stripped.values())
    logging.info(
        "[ions.summary] kept=%d stripped=%d metals_kept=%d metals_stripped=%d salts_kept=%d salts_stripped=%d",
        total_kept,
        total_stripped,
        category_totals["metal"]["kept"],
        category_totals["metal"]["stripped"],
        category_totals["salt"]["kept"],
        category_totals["salt"]["stripped"],
    )

    if warn_missing_center and radius > 0.0 and variant_token == "HOLO":
        logging.warning(
            "[ions] holo_salt_radius_set_but_no_center action=keep_salts radius=%.2f",
            radius,
        )

    if not remove_indices:
        logging.debug("[ions] no_monoatomic_hits remove=0")
        return 0

    for idx in sorted(remove_indices):
        text[idx] = None  # type: ignore

    rewritten = [ln for ln in text if ln is not None]
    try:
        path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    except Exception as exc:
        logging.warning("[ions] write_failed file=%s err=%s", path, exc)
        return 0

    return len(remove_indices)


def _is_element_token(sym):
    s = str(sym).strip().upper()
    return (len(s) in (1,2)) and (s in _ELEM_CANON)

# Treat as “ion-like” if it looks like an element token and is in the retain set
def _is_retained_ion(resname: str) -> bool:
    r = (resname or "").upper()
    return r in _RETAIN and _is_element_token(r)


# =============================
# Configuration helpers
# =============================
# Single unified config load for this module:
# Priority: environment overrides > config.txt keys > legacy key aliases > defaults.
try:
    from input_and_export_functions import load_config, validate_config
    _CFG = load_config("config.txt")
    validate_config(_CFG)
except Exception:
    # Fallback for older setups
    from installation import load_config as _legacy_load_config
    _CFG = _legacy_load_config()

# Keep both names so existing call-sites continue to work, ik this is lazy
config = _CFG
_CONFIG = _CFG

def _cfg(key: str, default: str = "", legacy_key: str | None = None) -> str:
    """env > config[key] > config[legacy_key] > default (all coerced to str)"""
    v = os.environ.get(key)
    if v not in (None, ""):
        return v
    if key in _CONFIG and str(_CONFIG.get(key)) != "":
        return str(_CONFIG.get(key))
    if legacy_key and legacy_key in _CONFIG and str(_CONFIG.get(legacy_key)) != "":
        return str(_CONFIG.get(legacy_key))
    return default

PHENIX_DIR          = _cfg("PHENIX_DIR", "", "phenix_dir")
PHENIX_LIB_PATH     = _cfg("PHENIX_LIB_PATH", "", "phenix_lib_path")
PHENIX_CLEAN_SCRIPT = _cfg("PHENIX_CLEAN_SCRIPT", "", "phenix_clean_script")
INPUT_DIR           = _cfg("INPUT_DIR", ".")

MGLTOOLS_PYTHON          = _cfg("MGLTOOLS_PYTHON", "")
PREPARE_RECEPTOR_SCRIPT  = _cfg("PREPARE_RECEPTOR_SCRIPT", "")
USE_MEEKO                = str(_cfg("USE_MEEKO", "")).lower() in ("1", "true", "yes")


# prefer explicit env var; else fall back to obabel on PATH
OPENBABEL_PATH = os.environ.get("OPENBABEL_PATH") or shutil.which("obabel") or ""

# Windows .bat wrappers (if using Windows Phenix)
PDBTOOLS_BAT     = str(Path(PHENIX_DIR) / "phenix.pdbtools.bat") if PHENIX_DIR else ""
MOLPROBITY_BAT   = str(Path(PHENIX_DIR) / "phenix.molprobity.bat") if PHENIX_DIR else ""
PHENIX_PYTHON_BAT= str(Path(PHENIX_DIR) / "phenix.python.bat")   if PHENIX_DIR else ""

def _pick_reduce_exe() -> str:
    """
    Choose the actual 'reduce' binary robustly.
    Precedence:
      1) $REDUCE_EXE (env) or config value
      2) Repo-relative local build: <repo>/tools/reduce/reduce_src/reduce
      3) Phenix conda_base/bin/reduce (next to PHENIX_DIR)
      4) reduce on PATH
    """
    # 1) Explicit env/config
    explicit = (os.environ.get("REDUCE_EXE") or os.environ.get("REDUCE_BIN") or
                _cfg("REDUCE_EXE", "", "reduce_exe"))
    if explicit and Path(explicit).exists():
        return explicit

    # 2) Repo-relative (portable)
    repo_root = Path(__file__).resolve().parent / "tools" / "reduce" / "reduce_src" / "reduce"
    if repo_root.exists() and os.access(str(repo_root), os.X_OK):
        return str(repo_root)

    # 3) Phenix conda-base candidate near PHENIX_DIR
    if PHENIX_DIR:
        phenix_root = Path(PHENIX_DIR).resolve().parent
        cb = phenix_root / "conda_base" / "bin" / "reduce"
        if cb.exists() and os.access(str(cb), os.X_OK):
            return str(cb)

    # 4) PATH fallback
    which = shutil.which("reduce") or shutil.which("reduce.exe")
    return which or "reduce"


REDUCE_EXE = _pick_reduce_exe()
logging.info("Using Reduce at: %s", REDUCE_EXE)

# Default HET dict (repo-relative) if not provided by env/config
DEFAULT_HET = Path(__file__).resolve().parent / "tools" / "reduce" / "reduce_wwPDB_het_dict.txt"
def _het_dict_path() -> str | None:
    hd = os.environ.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT", "")
    if hd and Path(hd).is_file():
        return hd
    return str(DEFAULT_HET) if DEFAULT_HET.is_file() else None


import shlex
from pathlib import Path
from typing import List, Tuple

def _persist_subproc(tag: str, cmd: List[str], cp: "subprocess.CompletedProcess",
                     outdir: Path, receptor_pdbqt: Path) -> None:
    """
    Save command, stdout, stderr to <work>/<tag>.* and emit a one-line summary with rc and receptor size.
    """
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    (outdir / f"{tag}.cmd.txt").write_text(" ".join(shlex.quote(x) for x in cmd), encoding="utf-8", errors="ignore")
    (outdir / f"{tag}.stdout.txt").write_text(cp.stdout or "", encoding="utf-8", errors="ignore")
    (outdir / f"{tag}.stderr.txt").write_text(cp.stderr or "", encoding="utf-8", errors="ignore")
    exists = receptor_pdbqt.exists()
    size = receptor_pdbqt.stat().st_size if exists else 0
    logging.info("[receptor-attempt %s] rc=%s exists=%s size=%d stderr=%s",
                 tag, cp.returncode, exists, size, str(outdir / f"{tag}.stderr.txt"))

def _first_last_lines(p: Path, n: int = 50) -> Tuple[list, list]:
    try:
        lines = (p.read_text(encoding="utf-8", errors="ignore")).splitlines()
        return lines[:n], lines[-n:]
    except Exception:
        return [], []

def _ok_receptor_file(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


# --- Chain-prune config helpers (reads via _cfg) ---
def _cfg_bool(key: str, default: bool = False) -> bool:
    v = (_cfg(key, str(int(default))) or "").strip().lower()
    return v in ("1","true","yes","y","on")

def _cfg_float(key: str, default: float) -> float:
    try:
        return float(_cfg(key, str(default)))
    except Exception:
        return float(default)

def _cfg_int(key: str, default: int) -> int:
    try:
        return int(float(_cfg(key, str(default))))
    except Exception:
        return int(default)

def _cfg_chain_keep_list() -> set[str]:
    raw = (_cfg("CHAIN_KEEP_LIST", "") or "")
    toks = [t.strip() for t in raw.replace(";",",").split(",") if t.strip()]
    return {t if len(t)==1 else t[:1] for t in toks}






def _post_write_element_guard(step: str, pdb_path: str | Path) -> None:
    """
    Normalize element columns after each receptor write step and log helium status.
    """
    p = Path(pdb_path)
    try:
        fixed = fix_element_columns_in_file(p, dst_path=p, rewrite_atoms=False)
    except Exception as e:
        logging.warning(
            "[helium] stage=post_write step=%s file=%s He->H=? residual_He=? note=elemfix_error:%s",
            step, p, e
        )
        return

    txt = p.read_text(encoding="utf-8", errors="ignore")
    he_count = scan_helium_counts(txt)
    logging.info(
        "[helium] stage=post_write step=%s file=%s He->H=%s residual_He=%d",
        step, p, fixed if isinstance(fixed, int) else -1, he_count
    )


def _meeko_preflight_or_fail(pdb_input: str | Path, work_dir: str | Path) -> Path:
    """
    Ensure the exact PDB Meeko will read contains no 'He'.
    Saves copy to work/meeko_input_pre_sanitize.pdb for forensics.
    """
    src = Path(pdb_input).resolve()
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    try:
        fix_element_columns_in_file(src, dst_path=src, rewrite_atoms=False)
    except Exception as e:
        logging.warning("[helium] stage=preflight file=%s note=elemfix_exception:%s", src, e)

    txt = src.read_text(encoding="utf-8", errors="ignore")
    he_count, first_hits = scan_helium_counts_with_hits(txt, max_hits=5)

    # Save what Meeko will actually see
    prefile = work / "meeko_input_pre_sanitize.pdb"
    prefile.write_text(txt, encoding="utf-8")

    if he_count > 0:
        logging.error(
            "[helium] stage=preflight file=%s He_count=%d rules=%s",
            src, he_count, rules_version(),
        )
        for hit in first_hits:
            logging.error("[helium] offender %s", hit)
        raise RuntimeError("helium_preflight_failed")

    return src




def _helium_postwrite_counter(step_name: str, pdb_path: str | Path) -> None:
    """
    Run the element fixer on the freshly written receptor PDB and emit:
      [helium] stage=post_write step=<name> file=<path> He->H=<n> residual_He=<n>
    Fail fast if residual_He > 0.
    """
    p = Path(pdb_path)
    try:
        before = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        before = ""

    # Count BEFORE
    he_before = scan_helium_counts(before.splitlines()) if before else 0

    # Fix elements in-place on the exact file we pass forward
    fix_element_columns_in_file(p, p, rewrite_atoms=True)

    try:
        after = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        after = ""

    # Count AFTER and compute delta
    he_after = scan_helium_counts(after.splitlines()) if after else 0
    delta = max(0, he_before - he_after)

    if _HE_POSTWRITE_VERBOSE:
        logging.info(
            "[helium] stage=post_write step=%s file=%s He->H=%d residual_He=%d",
            step_name, str(p), delta, he_after
        )

    if he_after > 0:
        try:
            bad = [ln for ln in after.splitlines()
                   if (" He" in ln) or (len(ln) >= 78 and ln[76:78].strip() == "HE")][:3]
            logging.warning("[helium] residual examples: %r", bad)
        except Exception:
            pass
        # Fail fast so we can see which step leaked helium
        raise RuntimeError("helium_residual_post_write")





try:
    from rdkit import Chem  # type: ignore
    _HAS_RDKIT = True
except Exception:
    _HAS_RDKIT = False

# =============================
# OS / WSL helpers
# =============================

def _is_wsl() -> bool:
    try:
        return "microsoft" in platform.release().lower() or "WSL_INTEROP" in os.environ
    except Exception:
        return False

def _bin_on_path(name: str) -> bool:
    return shutil.which(name) is not None

def _win_path(p: Union[str, Path]) -> str:
    s = str(p)
    if not _is_wsl():
        return s
    try:
        out = subprocess.check_output(["wslpath", "-w", s], text=True).strip()
        return out
    except Exception:
        return s

def _powershell(cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run([
        "powershell.exe", "-NoProfile", "-NonInteractive", "-Command", cmd
    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)





def _as_path(p) -> Path:
    return p if isinstance(p, Path) else Path(p)
def collapse_sanitized_once(p: Union[str, Path]) -> Path:
    """
    Return a Path with a sanitized basename (single pass; no filesystem changes).
    Keeps directory the same, replaces non [A-Za-z0-9._-] with underscores.
    """
    p = Path(p)
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", p.name)
    return p.with_name(safe)


def compute_control_centroids(ligands_dir: Union[str, Path]) -> list[tuple[float, float, float]]:
    """Return one centroid per *.pdb control ligand under ligands_dir."""
    ligands_dir = _as_path(ligands_dir)
    pts: list[tuple[float, float, float]] = []
    if not ligands_dir.exists():
        return pts
    for p in sorted(ligands_dir.glob("*.pdb")):
        try:
            n = 0
            sx = sy = sz = 0.0
            with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if not ln.startswith(("ATOM  ", "HETATM")):
                        continue
                    # drop hydrogens by element column (77–78)
                    if len(ln) >= 78 and ln[76:78].strip().upper() == "H":
                        continue
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                    sx += x; sy += y; sz += z; n += 1
            if n > 0:
                pts.append((sx / n, sy / n, sz / n))
        except Exception:
            # best-effort; ignore malformed ligand files
            pass
    return pts

def filter_waters_near_points(src_pdb: Union[str, Path],
                              dst_pdb: Union[str, Path],
                              points: list[tuple[float, float, float]],
                              radius_A: float) -> int:
    """
    Copy src_pdb → dst_pdb keeping all non-water records and only water residues
    with any atom within radius_A of any reference point. Returns #water residues kept.

    Water 3-letter residue names are loaded from YAML if available, else a fallback set.
    Optional: set env WATER_NAMES_YAML to a file path. YAML may be a list (["HOH","WAT",...])
    or a dict containing any of these keys: water_resnames, water, waters, solvent_water, water_aliases.
    """
    import os
    R2 = float(radius_A) * float(radius_A)

    # --- Load water residue names ---
    water3: set[str] = {
        # Common crystallographic water names
        "HOH", "WAT", "OH2", "DOD", "D2O",
        # MD water models that sometimes leak into PDB-like files
        "TIP", "TP3", "TP4", "TP5", "TIP3", "TIP4", "TIP5", "SOL", "W",
        # protonation variants
        "H3O", "OH-", "OHO", "H2O",
    }
    try:
        import yaml  # PyYAML
        # prefer explicit env; else try a conventional chemdb path beside this file
        default_yaml = os.path.join(os.path.dirname(__file__), "chemdb", "water_names.yaml")
        yaml_path = os.environ.get("WATER_NAMES_YAML", default_yaml)
        if os.path.exists(yaml_path):
            with open(yaml_path, "r", encoding="utf-8") as yf:
                y = yaml.safe_load(yf)
            candidates: list[str] = []
            if isinstance(y, dict):
                for k in ("water_resnames", "water", "waters", "solvent_water", "water_aliases"):
                    v = y.get(k, [])
                    if isinstance(v, str):
                        candidates.append(v)
                    elif isinstance(v, (list, tuple, set)):
                        candidates.extend(v)
            elif isinstance(y, (list, tuple, set)):
                candidates = list(y)
            # Normalize to 3-char PDB residue names
            for itm in candidates:
                try:
                    water3.add(str(itm).strip().upper()[:3])
                except Exception:
                    pass
    except Exception:
        # YAML absent or parse error → fall back silently
        pass

    keep: set[tuple[str, str]] = set()  # (chain, resseq+icode)

    # --- First pass: decide which water residues to keep ---
    try:
        with open(src_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                resn = ln[17:20].strip().upper()
                if resn not in water3:
                    continue
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                except Exception:
                    continue
                for (cx, cy, cz) in points:
                    dx = x - cx
                    dy = y - cy
                    dz = z - cz
                    if (dx * dx + dy * dy + dz * dz) <= R2:
                        keep.add((ln[21], ln[22:27]))  # chain, resseq+icode
                        break
    except Exception:
        # If anything goes wrong during reading, just fall back to keeping none.
        pass

    # --- Second pass: write out everything but only keep selected waters ---
    with open(src_pdb, "r", encoding="utf-8", errors="ignore") as fh, \
         open(dst_pdb, "w", encoding="utf-8") as out:
        for ln in fh:
            if ln.startswith(("ATOM  ", "HETATM")):
                resn = ln[17:20].strip().upper()
                if resn in water3:
                    key = (ln[21], ln[22:27])
                    if key in keep:
                        out.write(ln)
                    # else drop water line
                else:
                    out.write(ln)
            else:
                out.write(ln)

    return len(keep)




def count_waters_within(pdb_path: Union[str, Path],
                        point_xyz: tuple[float, float, float],
                        radius_A: float) -> int:
    """Count distinct HOH residues within radius_A of point_xyz."""
    R2 = float(radius_A) * float(radius_A)
    seen: set[tuple[str, str]] = set()
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")) or ln[17:20] != "HOH":
                    continue
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                except Exception:
                    continue
                dx = x - point_xyz[0]; dy = y - point_xyz[1]; dz = z - point_xyz[2]
                if (dx*dx + dy*dy + dz*dz) <= R2:
                    seen.add((ln[21], ln[22:27]))
    except Exception:
        return 0
    return len(seen)





def _cfg_env_or_default(key: str, default: Optional[str] = None) -> Optional[str]:
    """Prefer env var, then config.txt (sibling of this file), else default."""
    v = os.environ.get(key)
    if v:
        return v
    try:
        root = Path(__file__).resolve().parent
        cfg = root / "config.txt"
        if cfg.is_file():
            for line in cfg.read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, val = line.split("=", 1)
                if k.strip() == key:
                    return val.strip()
    except Exception:
        pass
    return default

def _canon_base(output_root: Path, pdb_id: str) -> Path:
    """Canonical per-protein base: processed_pdbs/<PDB>"""
    output_root = _as_path(output_root)
    return (output_root / pdb_id.upper()).resolve()

def _merge_dir(src: Path, dst: Path) -> None:
    """Merge src directory into dst; remove src after moving."""
    src, dst = src.resolve(), dst.resolve()
    if not src.exists():
        return
    dst.mkdir(parents=True, exist_ok=True)
    for root, dirs, files in os.walk(src):
        r = Path(root)
        rel = r.relative_to(src)
        (dst / rel).mkdir(parents=True, exist_ok=True)
        for d in dirs:
            (dst / rel / d).mkdir(parents=True, exist_ok=True)
        for f in files:
            s = r / f
            t = (dst / rel / f)
            if t.exists():
                try:
                    if s.stat().st_size == t.stat().st_size:
                        continue
                except Exception:
                    pass
            shutil.move(str(s), str(t))
    try:
        shutil.rmtree(src)
    except Exception:
        pass

def fold_legacy_layout(pdb_id: str, output_root) -> None:
    """
    Extended: migrate uppercase legacy dirs into canonical tree:
      <PDB>_NOLIG            -> processed_pdbs/<PDB>/nolig
      <PDB>_CLEANED_LIGANDS  -> processed_pdbs/<PDB>/ligands_raw
      <PDB>_nolig(.pdb)      -> processed_pdbs/<PDB>/nolig/<PDB>_nolig_phenix_clean.pdb
    """
    try:
        root = _as_path(output_root).resolve()
        pdb_idU = pdb_id.upper()
        base = _canon_base(root, pdb_idU)
        (base / "nolig").mkdir(parents=True, exist_ok=True)
        (base / "ligands_raw").mkdir(parents=True, exist_ok=True)

        legacy_dirs = [
            (root / f"{pdb_id}_nolig",                   base / "nolig"),
            (root / f"{pdb_id.lower()}_nolig",           base / "nolig"),
            (root / f"{pdb_idU}_NOLIG",                  base / "nolig"),
            (root / f"{pdb_id}_cleaned_ligands",         base / "ligands_raw"),
            (root / f"{pdb_id.lower()}_cleaned_ligands", base / "ligands_raw"),
            (root / f"{pdb_idU}_CLEANED_LIGANDS",        base / "ligands_raw"),
        ]
        for src, dst in legacy_dirs:
            if src.exists():
                logging.info("Migrating legacy directory %s -> %s", src, dst)
                _merge_dir(src, dst)

        candidates_files = [
            (root / f"{pdb_id}_nolig.pdb",         base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
            (root / f"{pdb_id.lower()}_nolig.pdb", base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
        ]
        for src, dst in candidates_files:
            if src.exists() and not dst.exists():
                logging.info("Moving legacy file %s -> %s", src, dst)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))
    except Exception as e:
        logging.warning("fold_legacy_layout (extended) failed for %s: %s", pdb_id, e)


# =============================
# Ligand pristine reference (optional; for RMSD downstream)
# =============================

def _write_pristine_reference(pdb_lig_path: Path) -> None:
    """Write pristine ligand copy next to ligands_raw (SDF + light atom map JSON)."""
    ref_dir = pdb_lig_path.parent.parent / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    ref_sdf = ref_dir / (pdb_lig_path.stem + ".sdf")
    ref_json = ref_dir / (pdb_lig_path.stem + ".map.json")

    if _HAS_RDKIT:
        m = Chem.MolFromPDBFile(str(pdb_lig_path), sanitize=False, removeHs=False)
        if m is not None:
            w = Chem.SDWriter(str(ref_sdf)); w.write(m); w.close()
            amap = []
            with open(pdb_lig_path, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    if ln.startswith(("ATOM","HETATM")):
                        amap.append({
                            "name": ln[12:16].strip(),
                            "element": (ln[76:78].strip() or ln[12:16].strip()[:1].upper())
                        })
            json.dump({"atoms": amap}, open(ref_json, "w"), indent=2)
            return

    # Fallback minimal SDF stub
    with open(ref_sdf, "w", encoding="utf-8") as out:
        out.write(f"{pdb_lig_path.stem}\n  -Pristine-\n\n")
        out.write("$$$$\n")

# =============================
# Canonical per-protein directory layout
# =============================

def canon_paths(pdb_id: str, output_root: Union[str, Path]) -> Dict[str, Path]:
    root = Path(output_root).resolve()
    base = root / pdb_id.upper()
    return {
        "protein_root": base,
        "raw":          base / "raw",
        "work":         base / "work",
        "ligands_raw":  base / "ligands_raw",
        "nolig":        base / "nolig",
        "receptor":     base / "receptor",
    }

# =============================
# Ligand extraction (single source lives here)
# =============================

def extract_ligands_from_filtered(filtered_pdb: Union[str, Path], out_dir: Union[str, Path]) -> List[Path]:
    """Extract non-water HETATM residues into individual PDBs (RES_CHAINRESI.pdb),
    and run text-level element repair on each (YAML rules)."""
    outd = Path(out_dir)
    outd.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    cur_key = None
    bucket: List[str] = []

    def flush():
        nonlocal bucket, cur_key, written
        if not bucket or cur_key is None:
            return
        resname, chain, resseq, icode = cur_key
        name = f"{resname}_{chain}{int(resseq)}"
        outp = outd / f"{name}.pdb"
        with open(outp, "w") as w:
            w.write(f"REMARK Extracted {name}\n")
            for ln in bucket:
                w.write(ln)
            w.write("TER\nEND\n")

        # normalize once (no-op if not needed)
        norm = collapse_sanitized_once(outp)
        if norm.name != outp.name:
            try:
                outp.replace(norm)
                outp = norm
            except Exception:
                pass

        # YAML-backed element repair (PDB only)
        try:
            fix_element_columns_in_file(outp, outp)
        except Exception as e:
            logging.warning("Element-fix skipped for %s: %s", outp.name, e)
        written.append(outp)
        bucket = []
        cur_key = None
        try:
            _write_pristine_reference(outp)
        except Exception as e:
            logging.warning("Could not write pristine reference for %s: %s", outp.name, e)

    with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith("HETATM"):
                continue
            resname = ln[17:20].strip().upper()
            if resname in _WATER_NAMES:
                continue  # never extract waters
            if resname in _RETAIN:
                continue  # don't extract cofactors/metals/ions you keep with protein
            chain = ln[21]
            resseq = ln[22:26].strip() or "0"
            icode = ln[26]
            key = (resname, chain, resseq, icode)
            if cur_key is None:
                cur_key = key
            if key != cur_key:
                flush()
                cur_key = key
            bucket.append(ln)
    flush()
    logging.info("Extracted %d ligand residues to %s", len(written), out_dir)
    return written



def element_fix_all_in_dir(dir_path: Union[str, Path], rewrite_atoms: bool = False) -> int:
    """
    Run the text-level element column fixer on every *.pdb under dir_path.
    Returns the number of files rewritten.
    """
    d = Path(dir_path)
    if not d.exists():
        return 0
    n = 0
    for p in d.rglob("*.pdb"):
        try:
            fix_element_columns_in_file(p, dst_path=p, rewrite_atoms=rewrite_atoms)
            n += 1
        except Exception as e:
            logging.warning("element_fix_all_in_dir skipped %s: %s", p, e)
    logging.info("Element column sweep fixed %d PDB files under %s", n, d)
    return n

def expose_ligand_intermediates_for_debug(src_dir: Union[str, Path],
                                          link_dir: Union[str, Path]) -> None:
    """
    Create/refresh a symlink 'link_dir' -> 'src_dir' for easy browsing.
    On Windows, tries directory junction fallback if symlink fails.
    Quietly skips if the source doesn't exist yet.
    """
    src = Path(src_dir).resolve()
    dst = Path(link_dir)

    # Avoid noise if source not present yet during early prep
    if not src.exists():
        logging.debug("intermediates: skip symlink, source missing: %s", src)
        return

    # Ensure link's parent exists
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logging.warning("intermediates: could not ensure parent dir for %s: %s", dst, e)
        return

    try:
        # Remove existing link/dir
        if dst.is_symlink() or dst.exists():
            try:
                if dst.is_symlink():
                    dst.unlink()
                else:
                    shutil.rmtree(dst)
            except Exception:
                pass

        # Create link (with Windows junction fallback)
        if os.name == "nt":
            try:
                os.symlink(src, dst, target_is_directory=True)
            except OSError:
                cmd = ['cmd', '/c', 'mklink', '/J', str(dst), str(src)]
                subprocess.run(cmd, check=True)
        else:
            os.symlink(src, dst, target_is_directory=True)

        logging.info("Exposed ligand intermediates: %s -> %s", dst, src)
    except Exception as e:
        logging.warning("intermediates: could not create symlink %s -> %s: %s", dst, src, e)



# =============================
# Element & Hydrogen utilities (PDB only)
# =============================

def hydrogenation_status(pdb_path: Union[str, Path]) -> Tuple[str, int, int, float]:
    h = heavy = 0
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            el = line[76:78].strip().upper()
            if el == "H":
                h += 1
            else:
                heavy += 1
    ratio = h / max(heavy, 1)
    if h == 0:
        return "NO_H", h, heavy, ratio
    if ratio < 0.20:
        return "SUSPECT_LOW_H", h, heavy, ratio
    return "HAS_H", h, heavy, ratio


def file_contains_hydrogens(pdb_path: Union[str, Path]) -> bool:
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                el = line[76:78].strip().upper()
                if el in {"H", "D", "T"}:
                    return True
    except Exception as e:
        logging.warning("Could not read %s to check for H atoms: %s", pdb_path, e)
    return False



def strip_monoatomic_ions_inplace(
    pdb_path: Union[str, Path],
    keep_resnames: Optional[Set[str]] = None,
    *,
    cfg: Optional[dict] = None,
    variant: Optional[str] = None,
    pocket_center: Optional[tuple[float, float, float]] = None,
    force_policy: Optional[str] = None,
) -> int:
    """Legacy wrapper that now delegates to the policy-aware ion stripping."""

    base_cfg = cfg if cfg is not None else config
    if isinstance(base_cfg, dict):
        cfg_obj: dict = dict(base_cfg)
    else:
        cfg_obj = {}

    if force_policy:
        cfg_obj = cfg_obj or {}
        cfg_obj["ION_STRIP_POLICY"] = force_policy

    extra_keep = set(keep_resnames or []) or None

    removed = _maybe_strip_ions(
        pdb_path,
        cfg=cfg_obj,
        variant=variant,
        pocket_center=pocket_center,
        extra_keep=extra_keep,
    )

    if removed:
        logging.info("Stripped %d monoatomic ions from %s", removed, pdb_path)

    return removed or 0

# Compatibility alias for older callers
def _strip_monoatomic_ions_inplace(*args, **kwargs):
    return strip_monoatomic_ions_inplace(*args, **kwargs)




def detect_catalytic_metals(pdb_path: Union[str, Path]) -> Set[str]:
    """
    Return the set of retained ion-like resnames present in the PDB (e.g., ZN, MG).
    Uses YAML retain list + element tokens. For diagnostics only.
    """
    RETAIN = set()
    for item in RULES.get("retain_in_receptor_resnames", []):
        for t in str(item).split(";"):
            tok = t.strip().upper()
            if tok:
                RETAIN.add(tok)

    ionic = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith("HETATM"):
                continue
            resname = ln[17:20].strip().upper()
            if (resname in RETAIN) and _is_element_token(resname):
                ionic.add(resname)

    if ionic:
        logging.info("Catalytic/retained ions present: %s", sorted(ionic))
    return ionic

def compare_ion_presence_between_pdb_and_pdbqt(clean_pdb: Union[str, Path],
                                                receptor_pdbqt: Union[str, Path]) -> None:
    """
    Log a warning if an ion present in the cleaned PDB is missing in receptor PDBQT.
    Uses YAML retain list + element tokens as "ions" definition.
    """
    def ions_in_pdb(p: Union[str, Path]) -> Set[Tuple[str, str]]:
        ions = set()
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith("HETATM"):
                    continue
                resname = ln[17:20].strip().upper()
                chain = ln[21]
                resi  = ln[22:26].strip()
                if _is_element_token(resname):
                    ions.add((resname, f"{chain}:{resi}"))
        return ions

    def ions_in_pdbqt(p: Union[str, Path]) -> Set[Tuple[str, str]]:
        ions = set()
        if not Path(p).exists():
            return ions
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                resname = ln[17:20].strip().upper()
                chain = ln[21]
                resi  = ln[22:26].strip()
                if _is_element_token(resname):
                    ions.add((resname, f"{chain}:{resi}"))
        return ions

    pdb_ions   = ions_in_pdb(clean_pdb)
    pdbqt_ions = ions_in_pdbqt(receptor_pdbqt)
    missing = pdb_ions - pdbqt_ions
    if missing:
        logging.warning("Ions present in PDB but not in receptor PDBQT: %s", sorted(missing))

def log_possible_metal_mislabels(pdb_path: Union[str, Path]) -> None:
    """
    Heuristic: flag HET residues whose resname is not an element token, but
    whose element column shows an element token and the residue has very few atoms.
    """
    suspects = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        by_res = {}
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                key = (ln[21], ln[22:26], ln[26], ln[17:20].strip().upper())
                by_res.setdefault(key, []).append(ln)

    for (chain, resi, icode, resname), atms in by_res.items():
        if len(atms) > 4:
            continue
        if _is_element_token(resname):
            continue
        elem_tokens = {ln[76:78].strip().upper() for ln in atms if len(ln) >= 78}
        if any(_is_element_token(e) for e in elem_tokens):
            suspects.append((resname, f"{chain}:{resi}", sorted(elem_tokens)))

    if suspects:
        logging.warning("Possible metal mislabels (check resname vs element cols): %s", suspects)


# =============================
# Structural cleaning utilities
# =============================

def filter_altlocs(pdb_input_path: Union[str, Path], pdb_output_path: Union[str, Path]) -> None:
    """Filter alternate locations (altLoc) deterministically and log decisions."""
    lines: List[str] = []
    atoms = {}
    removed_count = 0

    with open(pdb_input_path, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                atom_name = line[12:16]
                altLoc = line[16]
                chainID = line[21]
                resSeq = line[22:26].strip()
                iCode = line[26]
                key = (chainID, resSeq, iCode, atom_name.strip())
                atoms.setdefault(key, {})[altLoc] = line
            else:
                lines.append(line)

    filtered_atoms: List[str] = []
    for key, altloc_dict in atoms.items():
        altLocs = list(altloc_dict.keys())
        if len(altLocs) > 1:
            logging.info("AltLocs found for %s: %s", key, altLocs)
        if ' ' in altloc_dict:
            selected = ' '
        elif 'A' in altloc_dict:
            selected = 'A'
        else:
            selected = sorted(altloc_dict.keys())[0]
        if len(altloc_dict) > 1:
            removed = [alt for alt in altLocs if alt != selected]
            removed_count += len(removed)
            logging.info("Keeping altLoc '%s' for atom %s, removed %s", selected, key, removed)
        filtered_atoms.append(altloc_dict[selected])

    def _int_safe(s: str) -> int:
        s = s.strip()
        return int(s) if s and s.lstrip("-").isdigit() else 0

    filtered_atoms.sort(key=lambda l: (l[21], _int_safe(l[22:26]), l[26], l[12:16].strip()))

    with open(pdb_output_path, 'w', encoding="utf-8") as f:
        for line in lines:
            f.write(line)
        for atom_line in filtered_atoms:
            f.write(atom_line)
    _post_write_element_guard("altloc_filter", pdb_output_path)
    logging.info("Filtered altLocs in %s → %s (removed %d alternates)", pdb_input_path, pdb_output_path, removed_count)


def build_missing_loops(input_pdb: Union[str, Path], output_dir: Union[str, Path]) -> str:
    """Fill missing loops/residues using MODELLER; return output PDB path."""
    try:
        from modeller import environ, log
        from modeller.scripts import complete_pdb
    except Exception as e:
        logging.warning("MODELLER not available; skipping loop completion: %s", e)
        return str(input_pdb)

    output_pdb = os.path.join(str(output_dir), "modeller_filled.pdb")
    log.none()
    logging.info("Running MODELLER to complete missing parts of %s", input_pdb)

    try:
        env = environ()
        env.io.hetatm = True
        env.io.water = True
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        # MODELLER loop building
        mdl = complete_pdb(env, str(input_pdb))
        mdl.write(file=output_pdb)
        _post_write_element_guard("MODELLER", output_pdb)

        if os.path.exists(output_pdb):
            logging.info("MODELLER filled PDB saved to %s", output_pdb)
            return output_pdb
        logging.warning("MODELLER did not produce expected output: %s", output_pdb)
        return str(input_pdb)
    except Exception as e:
        logging.warning("MODELLER failed on %s: %s", input_pdb, e)
        return str(input_pdb)


def filter_invalid_chains(pdb_path: Union[str, Path], output_path: Union[str, Path]) -> None:
    """Remove entire chains lacking backbone atoms (CA/N/C/O)."""
    chains: Dict[str, List[str]] = defaultdict(list)
    valid_chains: Set[str] = set()

    with open(pdb_path, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                chain_id = line[21]
                atom_name = line[12:16].strip()
                chains[chain_id].append(line)
                if atom_name in {"CA", "N", "C", "O"}:
                    valid_chains.add(chain_id)
            else:
                chains["HEADER"].append(line)

    with open(output_path, 'w', encoding="utf-8") as f:
        for chain_id in chains:
            if chain_id == "HEADER" or chain_id in valid_chains:
                f.writelines(chains[chain_id])
            else:
                logging.warning("Skipping invalid chain '%s' (no CA atoms)", chain_id)


def _has_backbone_atoms(lines):
    req = {"N","CA","C","O"}
    seen = set()
    for ln in lines:
        if not ln.startswith(("ATOM","HETATM")): continue
        name = ln[12:16].strip()
        if name in req: seen.add(name)
    return req.issubset(seen)

def _group_by_chain(lines):
    chains = {}
    for ln in lines:
        if not ln.startswith(("ATOM","HETATM")): continue
        ch = ln[21]
        chains.setdefault(ch, []).append(ln)
    return chains

def _prune_chains_conservative(lines, keep_chains):
    out=[]
    for ln in lines:
        if ln.startswith(("ATOM","HETATM")) and ln[21] not in keep_chains:
            continue
        out.append(ln)
    return out

def _chains_to_keep(lines, pocket_center=None, r=12.0):
    chains = _group_by_chain(lines)
    keep=set()
    # rule 1: chain has full backbone atoms somewhere
    for ch, seg in chains.items():
        if _has_backbone_atoms(seg): keep.add(ch)
    # rule 2: optional geometric proximity if center known
    if pocket_center:
        near=set()
        x0,y0,z0 = pocket_center
        for ch, seg in chains.items():
            for ln in seg:
                if not ln.startswith(("ATOM","HETATM")): continue
                try:
                    x=float(ln[30:38]); y=float(ln[38:46]); z=float(ln[46:54])
                except Exception: continue
                if (x-x0)**2+(y-y0)**2+(z-z0)**2 <= r*r:
                    near.add(ch); break
        if near: keep = keep & near
    return keep or set(chains.keys())






# --- Hydrogen cleanup (geometry + CONECT) ---

def _parse_xyz(line: str):
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except Exception:
        return None


def conect_coverage(pdb_path: Union[str, Path]) -> float:
    atom_ids, conect_ids = set(), set()
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_ids.add(line[6:11].strip())
            elif line.startswith("CONECT"):
                parts = line.split()
                conect_ids.update(parts[1:])
    return len(conect_ids & atom_ids) / max(len(atom_ids), 1)


def remove_unbonded_atoms(pdb_path: Union[str, Path]) -> None:
    bonded_atoms: Set[str] = set()
    all_atoms: List[str] = []
    with open(pdb_path, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("CONECT"):
                parts = line.split()
                for atom_serial in parts[1:]:
                    bonded_atoms.add(atom_serial)
            elif line.startswith(("ATOM", "HETATM")):
                all_atoms.append(line)

    filtered: List[str] = []
    for line in all_atoms:
        atom_serial = line[6:11].strip()
        if atom_serial in bonded_atoms or line[76:78].strip() != 'H':
            filtered.append(line)
        else:
            logging.info("Removed unbonded hydrogen: %s", line.strip())

    with open(pdb_path, 'w', encoding="utf-8") as f:
        f.writelines(filtered)


def remove_implausible_hydrogens_by_distance(pdb_path: Union[str, Path]) -> None:
    atoms: List[Tuple[str, str, Optional[Tuple[float, float, float]]]] = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    atoms.append((line, "UNK", (None, None, None)))
                    continue
                el = line[76:78].strip() or line[12:16].strip()[:1]
                atoms.append((line, el.upper(), (x, y, z)))

    kept: List[str] = []
    heavy_coords = [a[2] for a in atoms if a[1] != "H" and a[2][0] is not None]

    def near_heavy(coord: Optional[Tuple[float, float, float]]) -> bool:
        if coord[0] is None:
            return True
        x, y, z = coord
        for X, Y, Z in heavy_coords:
            dx = x - X; dy = y - Y; dz = z - Z
            if (dx*dx + dy*dy + dz*dz) <= (1.35 * 1.35):
                return True
        return False

    for line, el, coord in atoms:
        if el != "H" or near_heavy(coord):
            kept.append(line)
        else:
            logging.info("Removed implausible H: %s", line.strip())

    with open(pdb_path, "w", encoding="utf-8") as out:
        out.writelines(kept)


def clean_hydrogens(pdb_path: Union[str, Path], use_conect_if_reliable: bool = True, conect_min_cov: float = 0.6) -> None:
    cov = conect_coverage(pdb_path) if use_conect_if_reliable else 0.0
    if cov >= conect_min_cov:
        remove_unbonded_atoms(pdb_path)
    remove_implausible_hydrogens_by_distance(pdb_path)
    _post_write_element_guard("hydrogen_cleanup", pdb_path)


# =============================
# Cofactors / waters policy
# =============================

def _is_metal(resname: str) -> bool:
    """Treat as 'metal/ion' iff YAML retain list contains this resname AND it looks like an element token."""
    rn = (resname or "").upper()
    return (rn in _RETAIN) and _is_element_token(rn)

def _cofactor_policy_keep(resname: str) -> bool:
    """
    YAML-only policy:
      • Keep anything listed under 'retain_in_receptor_resnames'
      • Keep standard residues (handled elsewhere)
      • Water handling is done in strip_nonstandard_residues(), so return False here for waters.
      • Everything else → remove
    """
    rn = (resname or "").upper()
    if rn in _RETAIN:
        return True
    if rn in _WATER_NAMES:
        return False
    return False


def strip_nonstandard_residues(input_pdb: Union[str, Path], output_pdb: Union[str, Path]) -> Tuple[int, str]:
    """
    Remove nonstandard residues while keeping what the YAML says to keep
    (retain_in_receptor_resnames) and standard amino acids.
    Water handling still honors your numeric policy (radius/B-factor) but
    the water names themselves come from the YAML.
    """
    standard_residues = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
        "HID","HIE","HIP","SEC","PYL","MSE",
    }

    removed: Set[str] = set()
    kept_lines: List[str] = []

    # Estimate pocket center from any YAML-retained cofactors/metals present
    cofm_xyz: List[Tuple[float, float, float]] = []
    with open(input_pdb, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper()
            if resname in _RETAIN:
                xyz = _parse_xyz(line)
                if xyz: cofm_xyz.append(xyz)

    pocket_center: Optional[Tuple[float, float, float]] = None
    if cofm_xyz:
        from statistics import fmean
        xs, ys, zs = zip(*cofm_xyz)
        pocket_center = (fmean(xs), fmean(ys), fmean(zs))

    water_policy = (config.get("WATER_POLICY", "site_only") or "site_only").lower()
    water_radius = float(config.get("WATER_SITE_RADIUS_ANG", 6.0))
    water_bmax   = float(config.get("WATER_MAX_BFACTOR", 60.0))

    with open(input_pdb, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM  "):
                resname = line[17:20].strip().upper()
                if resname in standard_residues:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
                continue

            if line.startswith("HETATM"):
                resname = line[17:20].strip().upper()

                # Water names come from YAML (no hardcoded list here)
                if resname in _WATER_NAMES:
                    if water_policy == "keep_all":
                        kept_lines.append(line); continue
                    if water_policy == "remove_all":
                        removed.add(resname);      continue
                    # site_only / auto: keep if near pocket center and not too mobile
                    if pocket_center is not None:
                        xyz = _parse_xyz(line)
                        keep = False
                        if xyz:
                            dx = xyz[0]-pocket_center[0]; dy = xyz[1]-pocket_center[1]; dz = xyz[2]-pocket_center[2]
                            dist2 = dx*dx + dy*dy + dz*dz
                            if dist2 <= water_radius*water_radius:
                                try:
                                    b = float(line[60:66]); keep = (b <= water_bmax)
                                except Exception:
                                    keep = True
                        if keep:
                            kept_lines.append(line)
                        else:
                            removed.add(resname)
                    else:
                        # no pocket estimate: default to remove unless keep_all
                        removed.add(resname)
                    continue

                # Retain anything listed in YAML retain block (cofactors, metals, ions, etc.)
                if resname in _RETAIN:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
                continue

            # non-coordinate records pass through
            kept_lines.append(line)

    with open(output_pdb, 'w', encoding="utf-8") as f:
        f.writelines(kept_lines)
    _post_write_element_guard("strip_nonstandard", output_pdb)

    logging.info("Removed nonstandard residues (YAML-driven): %s", sorted(removed))
    if pocket_center:
        logging.info("Estimated pocket center: (%.2f, %.2f, %.2f)", *pocket_center)
    return len(removed), str(output_pdb)


def quick_element_histogram(pdb_path: Union[str, Path]) -> None:
    cnt = Counter()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM","HETATM")):
                el = ln[76:78].strip().upper()
                cnt[el or ""] += 1
    logging.info("[Elem histogram %s] %s", os.path.basename(str(pdb_path)), dict(sorted(cnt.items())))

    
def assert_no_metal_in_peptidic(pdb_path: Union[str, Path]) -> None:
    # Build a peptide-like name set from YAML
    pep_name_lines = RULES["element_sets"].get("peptide_like_names", [])
    peptidey = set()
    for line in pep_name_lines:
        for tok in str(line).split(";"):
            t = tok.strip().upper()
            if t:
                peptidey.add(t)

    # PTM whitelist: YAML override if present; else default PTR/SEP/TPO
    ptm_yaml = set(_flatten_semicolons(RULES.get("element_sets", {}).get("ptm_resnames", [])))
    ptm_resnames = ptm_yaml or {"PTR", "SEP", "TPO"}

    # Element tokens considered "ionic" from YAML context (retain + element list)
    ionic_tokens = {r for r in _RETAIN if _is_element_token(r)}

    res_atoms = defaultdict(list)
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ","HETATM")):
                key = (ln[21], ln[22:26], ln[26], ln[17:20].strip().upper())
                res_atoms[key].append(ln)

    offenders = []
    for key, lines in res_atoms.items():
        chain, resi, icode, resname = key
        # Skip known PTMs entirely
        if resname in ptm_resnames:
            continue
        names = [ln[12:16].strip().upper() for ln in lines]
        if not names:
            continue
        hits = sum((n in peptidey) or (n[:2] in {"OE","NE","OD","ND","SD"}) for n in names)
        pep_like = hits >= max(4, 0.6*len(names))
        if not pep_like:
            continue
        # if any atom's element is an ionic token (from YAML), flag
        if any((ln[76:78].strip().upper() in ionic_tokens) for ln in lines):
            offenders.append(key)

    # Compact info line to make the warning greppable/noisy only when meaningful
    logging.info("[peptide-ion check] offenders=%s peptidey_n=%d ionic_tokens_n=%d",
                 offenders, len(peptidey), len(ionic_tokens))

    if offenders:
        logging.warning("Peptide-like residues contain ionic elements (check labeling): %s", offenders)


#DELETE CHAINS HELPERS
def detect_pocket_center_from_ligands(filtered_pdb: Union[str, Path],
                                      ligands_dir: Union[str, Path]
) -> Optional[Tuple[float,float,float]]:
    """
    Prefer center of extracted ligands in ligands_raw/; fallback to YAML-retained cofactors/metals
    present in filtered_pdb. Returns (x,y,z) or None.
    """
    from statistics import fmean
    pts: list[Tuple[float,float,float]] = []

    # 1) Extracted ligands (*.pdb) under ligands_dir
    ligd = Path(ligands_dir)
    for p in ligd.glob("*.pdb"):
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    if ln.startswith(("ATOM","HETATM")):
                        el = ln[76:78].strip().upper()
                        if el == "H":
                            continue
                        try:
                            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                            pts.append((x,y,z))
                        except Exception:
                            pass
        except Exception:
            pass

    # 2) Fallback: retained cofactors/metals from YAML in filtered_pdb
    if not pts:
        with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith("HETATM"):
                    continue
                resname = ln[17:20].strip().upper()
                if resname in _RETAIN:
                    try:
                        x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                        pts.append((x,y,z))
                    except Exception:
                        pass

    if not pts:
        return None
    xs, ys, zs = zip(*pts)
    return (fmean(xs), fmean(ys), fmean(zs))
def score_chain_contacts(pdb_path: Union[str, Path], keep_chains: set[str]) -> Dict[str, int]:
    """
    Return heavy-atom pair counts within CHAIN_CONTACT_DIST_ANG between each non-kept chain
    and the *union* of kept chains.
    """
    dist = _cfg_float("CHAIN_CONTACT_DIST_ANG", 5.0)
    atoms_by_chain: dict[str, list[Tuple[float,float,float]]] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith(("ATOM  ","HETATM")):
                continue
            el = ln[76:78].strip().upper()
            if el == "H":
                continue
            c = ln[21]
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except Exception:
                continue
            atoms_by_chain.setdefault(c, []).append((x,y,z))

    kept_pts = []
    for kc in keep_chains:
        kept_pts.extend(atoms_by_chain.get(kc, []))

    out: Dict[str,int] = {}
    if not kept_pts:
        return {c: 0 for c in atoms_by_chain}  # nothing to compare to

    d2 = dist*dist
    for c, pts in atoms_by_chain.items():
        if c in keep_chains:
            continue
        cnt = 0
        # simple O(N*M); early stop not necessary but harmless
        for x,y,z in pts:
            for X,Y,Z in kept_pts:
                dx = x-X; dy = y-Y; dz = z-Z
                if (dx*dx + dy*dy + dz*dz) <= d2:
                    cnt += 1
        out[c] = cnt
    return out
def select_chains_to_keep(filtered_pdb: Union[str, Path],
                          ligands_dir: Union[str, Path],
                          cfg=None
) -> set[str]:
    """
    Decide chains to keep using pocket proximity, contact counts, keep-list, and safety rails.
    Emits a compact decision table when CHAIN_LOG_DECISIONS is enabled.
    """
    radius = _cfg_float("CHAIN_POCKET_RADIUS_ANG", 10.0)
    min_contacts = _cfg_int("CHAIN_MIN_CONTACTS", 200)
    keep_list = _cfg_chain_keep_list()
    log_dec = _cfg_bool("CHAIN_LOG_DECISIONS", True)

    pocket = detect_pocket_center_from_ligands(filtered_pdb, ligands_dir)

    # Build per-chain stats
    chains: dict[str, dict] = {}
    all_chains: set[str] = set()
    ca_counts: dict[str,int] = {}
    min_dists: dict[str,float] = {}

    def _dist2(pt, xyz):
        dx = pt[0]-xyz[0]; dy = pt[1]-xyz[1]; dz = pt[2]-xyz[2]
        return dx*dx + dy*dy + dz*dz

    with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ","HETATM")):
                c = ln[21]
                all_chains.add(c)
                if ln.startswith("ATOM  ") and ln[12:16].strip() == "CA":
                    ca_counts[c] = ca_counts.get(c,0)+1
                if pocket is not None:
                    try:
                        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                        d2 = _dist2(pocket, xyz)
                        md = min_dists.get(c, float("inf"))
                        if d2 < md:
                            min_dists[c] = d2
                    except Exception:
                        pass

    if not all_chains:
        return set()  # nothing to decide

    # Initial kept set
    kept: set[str] = set(keep_list)

    # Always keep blank-chain records conservatively
    if " " in all_chains:
        kept.add(" ")

    # Pocket proximity
    if pocket is not None:
        r2 = radius*radius
        for c in all_chains:
            if min_dists.get(c, float("inf")) <= r2:
                kept.add(c)

    # If we still have no pocket & no keep_list -> abort prune
    if (pocket is None) and not keep_list:
        if log_dec:
            logging.info("[chain_prune] skipped reason=no_pocket_center or keep_list")
        return set()  # signal: do not prune

    # Contact expansion to capture interfaces
    contacts = score_chain_contacts(filtered_pdb, kept)
    for c, cnt in contacts.items():
        if cnt >= min_contacts:
            kept.add(c)

    # Safety rails
    if not kept:
        if log_dec:
            logging.info("[chain_prune] conservative_abort reason=kept_empty")
        return set()
    largest_ca_chain = max(ca_counts, key=lambda c: ca_counts.get(c,0)) if ca_counts else None
    drop_set = {c for c in all_chains if c not in kept}
    if len(drop_set) > len(all_chains) / 2:
        if log_dec:
            logging.warning("[chain_prune] conservative_abort reason=drop_gt_50pct all=%s drop=%s",
                            sorted(all_chains), sorted(drop_set))
        return set()
    if largest_ca_chain and (largest_ca_chain in drop_set):
        if log_dec:
            logging.warning("[chain_prune] conservative_abort reason=largest_CA_chain_would_be_dropped largest=%s",
                            largest_ca_chain)
        return set()

    # Decision table
    if log_dec:
        rows = []
        for c in sorted(all_chains):
            rows.append({
                "chain": c,
                "CA_count": ca_counts.get(c, 0),
                "near_pocket": ("yes" if (pocket is not None and min_dists.get(c,float("inf")) <= (radius*radius)) else "no"),
                "min_dist": (0.0 if pocket is None else (min_dists.get(c,float("inf"))**0.5)),
                "contact_count_to_kept": contacts.get(c, 0),
                "decision": ("keep" if c in kept else "drop"),
                "reason": ("whitelist" if c in keep_list else
                           "near_pocket" if (pocket is not None and min_dists.get(c,float("inf")) <= (radius*radius)) else
                           "interface_contacts" if contacts.get(c,0) >= min_contacts else
                           "far_and_sparse")
            })
        # Emit compact table
        hdr = "# chain  CA  near  min_d  contacts  keep  reason"
        logging.info(hdr)
        for r in rows:
            logging.info("  %-5s  %-3d %-5s %6.2f    %-7d %-4s  %s",
                         r["chain"], r["CA_count"], r["near_pocket"], r["min_dist"],
                         r["contact_count_to_kept"], ("yes" if r["decision"]=="keep" else "no"),
                         r["reason"])
        kept_ids = "".join(sorted(kept)) or "-"
        dropped_ids = "".join(sorted(drop_set)) or "-"
        logging.info("[chain_prune] kept=%s dropped=%s reason=see_table", kept_ids, dropped_ids)

    return kept
def prune_to_chains(input_pdb: Union[str, Path],
                    kept_chains: set[str],
                    output_pdb: Union[str, Path]) -> None:
    """
    Write only coordinate records belonging to kept_chains; copy all non-coordinate lines through.
    """
    kept = set(kept_chains or set())
    with open(input_pdb, "r", encoding="utf-8", errors="ignore") as f, \
         open(output_pdb, "w", encoding="utf-8") as w:
        for ln in f:
            if ln.startswith(("ATOM  ","HETATM")):
                c = ln[21]
                if (c in kept):
                    w.write(ln)
            else:
                w.write(ln)
    _post_write_element_guard("chain_pruned", output_pdb)











def _maybe_get_target_ph() -> float | None:
    # Priority: explicit env, then context_ph json drop, else None.
    v = os.environ.get("TARGET_PH", "").strip()
    if v:
        try: return float(v)
        except Exception: pass
    # allow a sidecar json dumped by context_ph: <pdb_id>.ph.json with {"target_pH": 7.8}
    try:
        pdb_base = os.path.splitext(os.path.basename(input_pdb_path))[0]
        sidecar = Path(input_pdb_path).parent / f"{pdb_base}.ph.json"
        if sidecar.exists():
            import json
            data = json.loads(sidecar.read_text())
            t = data.get("target_pH", None)
            if isinstance(t, (int, float)): return float(t)
    except Exception:
        pass
    return None

def _protonate_with_pdb2pqr_if_available(nolig_pdb_path: str, out_dir: Path, logger) -> tuple[Path, bool, str|None]:
    """
    Try PDB2PQR at the *pipeline pH* if available; fall back to the input PDB.
    Returns: (pdb_for_reduce, used_pdb2pqr, propka_log_path_or_None)
    """
    # Choose a pH: use context_ph if you already resolved one before calling this,
    # otherwise default to 7.0 (harmless; you can feed in your target later).
    try:
        from context_ph import select_ph_values_for_protonation
        phs = select_ph_values_for_protonation(nolig_pdb_path)
        target_ph = float(phs[0]) if phs else 7.0
    except Exception:
        target_ph = 7.0

    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    pdb_from_p2p, pk_log = pdb2pqr_protonate(nolig_pdb_path, target_ph, out_dir)
    if pdb_from_p2p:
        logger.info("PDB2PQR succeeded at pH %.2f -> %s", target_ph, pdb_from_p2p)
        return Path(pdb_from_p2p), True, pk_log
    else:
        logger.warning("PDB2PQR unavailable/failed; Reduce will build hydrogens.")
        return Path(nolig_pdb_path), False, None






# ============================
# End-to-end Cleaning Pipeline
# =============================
def clean_pdb(pdb_file: Union[str, Path], output_root: Union[str, Path]) -> Optional[str]:
    """Run the full cleaning pipeline and return path to final cleaned PDB (receptor)."""
    output_root = str(output_root)
    Path(output_root).mkdir(parents=True, exist_ok=True)

    raw_stem = os.path.splitext(os.path.basename(str(pdb_file)))[0]
    pdb_id = re.sub(r"(_nolig(_cleaned)?|_cleaned)$", "", raw_stem, flags=re.I).upper()
    logging.info("[prep.id] clean_pdb stem=%s -> base_id=%s", raw_stem, pdb_id)
    paths = canon_paths(pdb_id, output_root)
    logging.info("[prep.paths] protein_root=%s receptor=%s nolig=%s work=%s",
                 paths["protein_root"], paths["receptor"], paths["nolig"], paths["work"])


    logging.info("[proteinprep] entering clean_pdb pdb_file=%s output_root=%s", pdb_file, output_root)

    for d in ["protein_root", "raw", "work", "ligands_raw", "nolig", "receptor"]:
        paths[d].mkdir(parents=True, exist_ok=True)

    # (1) Working copy → raw/
    working_pdb = paths["raw"] / f"{pdb_id}_working.pdb"
    shutil.copyfile(str(pdb_file), working_pdb)
    _helium_postwrite_counter("copy_working", working_pdb)
    # (2) AltLoc filtering → raw/filtered.pdb
    filtered_pdb = paths["raw"] / f"{pdb_id}_filtered.pdb"
    filter_altlocs(working_pdb, filtered_pdb)

    # (2a) EARLY text-level element fix (YAML-driven), before any heavy tools
    try:
        fix_element_columns_in_file(filtered_pdb, filtered_pdb, rewrite_atoms=True)
        _helium_postwrite_counter("elemfix_filtered", filtered_pdb)
        logging.info("Early text-level element fix applied to %s", filtered_pdb)
    except Exception as e:
        logging.warning("Early text-level element fix skipped for %s: %s", filtered_pdb, e)

    # (3) Extract ligands now (controls live here), with YAML element repair per-file
    _ = extract_ligands_from_filtered(filtered_pdb, paths["ligands_raw"])

    if os.environ.get("EARLY_CHAIN_PRUNE", "1").lower() not in {"0", "false", "no"}:
        try:
            with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as fh:
                lines = fh.readlines()
            keep = _chains_to_keep(lines, pocket_center=None)  # center can be wired later
            orig = {ln[21] for ln in lines if ln.startswith(("ATOM", "HETATM"))}
            if keep and keep != orig:
                pruned = _prune_chains_conservative(lines, keep)
                with open(filtered_pdb, "w", encoding="utf-8") as out:
                    out.writelines(pruned)
                logging.info("[chains] early-pruned chains keep=%s drop=%s", "".join(sorted(keep)),
                             "".join(sorted(orig - keep)))
        except Exception as e:
            logging.warning("[chains] early prune skipped: %s", e)
    # (3a) Optional early chain-prune (conservative, pocket-aware)
    source_for_strip = filtered_pdb
    if _cfg_bool("CHAIN_PRUNE", False):
        try:
            kept = select_chains_to_keep(filtered_pdb, paths["ligands_raw"], config)
            if kept:
                pruned_filtered = paths["raw"] / f"{pdb_id}_filtered_pruned.pdb"
                prune_to_chains(filtered_pdb, kept, pruned_filtered)
                logging.info("[chain_prune] using pruned source for step (4): %s", pruned_filtered)
                source_for_strip = pruned_filtered
            else:
                logging.info("[chain_prune] not applied (kept empty or conservative_abort); using unpruned file")
        except Exception as e:
            logging.warning("[chain_prune] skipped due to exception: %s", e)

    # (4) Strip nonstandard from protein (policy aware) → work/stripped.pdb
    stripped_pdb = paths["work"] / f"{pdb_id}_stripped.pdb"
    removed_count, _out = strip_nonstandard_residues(source_for_strip, stripped_pdb)
    logging.info("Removed %d nonstandard residue lines.", removed_count)

    # (5) Element fix → MODELLER → element fix again (PDB only)
    elemfix_pdb = paths["work"] / f"{pdb_id}_elemfix.pdb"
    
    fix_pdb_elements(stripped_pdb, elemfix_pdb)
    _helium_postwrite_counter("elemfix_before_modeller", elemfix_pdb)
    loop_fixed_pdb = build_missing_loops(elemfix_pdb, paths["work"])
    fix_pdb_elements(loop_fixed_pdb, loop_fixed_pdb)
    _helium_postwrite_counter("elemfix_after_modeller", loop_fixed_pdb)

    # MODELLER (detect whether a new file was actually produced)
    modeller_ok = (
            os.path.basename(loop_fixed_pdb) == "modeller_filled.pdb"
            and os.path.isfile(loop_fixed_pdb)
    )


    receptor_pdb = paths["receptor"] / f"{pdb_id}_cleaned.pdb"
    # (6) Optional external Phenix polish (non-fatal if missing)
    phenix_ok = False
    _use_phenix = str(config.get("use_phenix", config.get("USE_PHENIX", "false"))).strip().lower() in ("1", "true",
                                                                                                       "yes")
    if _use_phenix:        # Water policy & radius
        _remove_waters = str(_cfg_env_or_default("REMOVE_WATERS", "true")).strip().lower() in ("1", "true", "yes")
        _policy = (_cfg_env_or_default("WATER_KEEP_POLICY", "none") or "none").strip().lower()
        _keep_R = float(_cfg_env_or_default("KEEP_WATERS_WITHIN_A", "6.0") or 6.0)

        # Default: feed Phenix the loop-fixed input
        _phenix_in = loop_fixed_pdb

        if _remove_waters and _policy != "none":
            # Policy active: derive reference points
            ref_pts: list[tuple[float, float, float]] = []
            # Prefer chosen center when available; here we’re early, so fall back to control centroids
            try:
                ref_pts = compute_control_centroids(paths["ligands_raw"])
            except Exception:
                ref_pts = []
            if ref_pts:
                _phenix_in = paths["work"] / f"{pdb_id}_prefiltered_waters.pdb"
                kept = filter_waters_near_points(loop_fixed_pdb, _phenix_in, ref_pts, _keep_R)
                logging.info("[waters] policy=%s kept=%d within %.1f Å of %d centers",
                             _policy, kept, _keep_R, len(ref_pts))
            else:
                logging.info("[waters] policy=%s but no reference points found; skipping prefilter", _policy)

        # Blanket removal only when policy is 'none'
        _phenix_remove = bool(_remove_waters and _policy == "none")
        phenix_ok = run_phenix_pdbtools(input_pdb=_phenix_in, output_pdb=receptor_pdb, remove_waters=_phenix_remove)

        # (6b) Dry-run sanity: count HOH within 8 Å of control-centroid center in final receptor
        try:
            ref_pts = compute_control_centroids(paths["ligands_raw"])
            center0 = None
            if ref_pts:
                # quick average as an approximate center for the dry-run note
                cx = sum(p[0] for p in ref_pts) / len(ref_pts)
                cy = sum(p[1] for p in ref_pts) / len(ref_pts)
                cz = sum(p[2] for p in ref_pts) / len(ref_pts)
                center0 = (cx, cy, cz)
            if center0:
                kept8 = count_waters_within(receptor_pdb, center0, 8.0)
                logging.info("[waters] dry-run kept_within_8A=%d center=(%.2f,%.2f,%.2f) file=%s",
                             kept8, center0[0], center0[1], center0[2], receptor_pdb)
        except Exception as _e:
            logging.debug("[waters] dry-run check skipped: %s", _e)

    if not phenix_ok:
        shutil.copyfile(loop_fixed_pdb, receptor_pdb)





    _helium_postwrite_counter("phenix_or_copy_receptor", receptor_pdb)
    
    # choose the file to pass downstream
    pdb_for_reduce = loop_fixed_pdb if modeller_ok else elemfix_pdb

    reduce_deferred = True

    logging.info(
        "[proteinprep] steps: Reduce=%s Phenix=%s MODELLER=%s",
        "deferred" if reduce_deferred else "applied",
        str(phenix_ok),
        str(modeller_ok),
    )
    if modeller_ok:
        sz = os.path.getsize(loop_fixed_pdb)
        logging.info("[proteinprep] modeller_out=%s size=%d", loop_fixed_pdb, sz)
    else:
        logging.info("[proteinprep] modeller_out=none (kept %s)", elemfix_pdb)

    if phenix_ok:
        sz = os.path.getsize(receptor_pdb)
        logging.info("[proteinprep] phenix_applied_to=%s size=%d", receptor_pdb, sz)

    # (7) Hydrogen cleanup & chain validation
    debulked_pdb = paths["work"] / f"{pdb_id}_debulked.pdb"
    shutil.copyfile(receptor_pdb, debulked_pdb)
    clean_hydrogens(debulked_pdb, use_conect_if_reliable=True, conect_min_cov=0.6)

    chain_validated_pdb = paths["work"] / f"{pdb_id}_validated.pdb"
    filter_invalid_chains(debulked_pdb, chain_validated_pdb)
    _helium_postwrite_counter("chain_validate", chain_validated_pdb)

    try:
        _txt_before = Path(chain_validated_pdb).read_text(encoding="utf-8", errors="ignore")
    except Exception:
        _txt_before = ""

    try:
        # Rewrite PDB element columns (77–78) using the unified rules (also for ATOM when rewrite_atoms=True)
        fix_element_columns_in_file(chain_validated_pdb, chain_validated_pdb, rewrite_atoms=True)

        # Secondary invariant: if any H-named atom still carries He, fix and summarize
        _txt_after = Path(chain_validated_pdb).read_text(encoding="utf-8", errors="ignore")
        _fixed_text, _nname = assert_no_helium_in_hydrogen_names(_txt_after)
        if _nname > 0:
            Path(chain_validated_pdb).write_text(_fixed_text, encoding="utf-8")

        # Grep-friendly one-liner with He→H delta
        _before = scan_helium_counts(_txt_before)
        _after = scan_helium_counts(Path(chain_validated_pdb).read_text(encoding="utf-8", errors="ignore"))
        _delta = max(0, _before - _after)
        logging.info(f"[elem-fix] file={Path(chain_validated_pdb).name} stage=preflight He->H={_delta}")
    except Exception as _e:
        logging.warning(f"[elements] receptor preflight failed for {Path(chain_validated_pdb).name}: {_e}")
        
    # (8) Protonation (Reduce when safe; else Open Babel fallback)
    def _present_resnames(pdb_path: Path) -> set[str]:
        res = set()
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith(("ATOM  ", "HETATM")):
                    res.add(ln[17:20].strip().upper())
        return res

    present_resnames = _present_resnames(chain_validated_pdb)
    NUC_LIKE = set(_flatten_semicolons(RULES.get("nucleotide_like_resnames", [])))
    use_reduce = not any(r in NUC_LIKE for r in present_resnames)
    # ---- Prefer PDB2PQR (with PROPKA) at pipeline pH; fall back to Reduce ----
    # This uses the helper already defined earlier in this file.
    pdb_for_reduce, used_pdb2pqr, pk_log = _protonate_with_pdb2pqr_if_available(
        str(chain_validated_pdb),         # protonate the validated, ligand-free coordinates
        str(paths["work"]),               # write PROPKA/PDB2PQR artifacts into the work directory
        logging
    )

    # If PDB2PQR succeeded, pdb_for_reduce now has hydrogens and titration states.
    # assign_protonation_states() will detect H presence and run Reduce WITHOUT -BUILD,
    # i.e., do flips/cleanup only. If PDB2PQR failed, pdb_for_reduce == chain_validated_pdb
    # and Reduce will run with -BUILD as needed.

    if not use_reduce:
        logging.info("[protonation] Skipping Reduce due to detected nucleotides; using OpenBabel path.")

    reduced_pdb = paths["work"] / f"{pdb_id}_reduced.pdb"
    assign_protonation_states(
        pdb_for_reduce,
        reduced_pdb,
        reduce_exe=REDUCE_EXE if use_reduce else None,  # skip Reduce for nucleotide cofactors
    )

    _helium_postwrite_counter("reduce_or_fallback", reduced_pdb)
    print(f"[proteinprep] Reduce/alt_protonation wrote={Path(reduced_pdb).is_file()} -> {reduced_pdb}")

    # (9) Final element fix and sanity on the protonated file
    fix_pdb_elements(reduced_pdb)
    _helium_postwrite_counter("elemfix_after_reduce", reduced_pdb)

    quick_element_histogram(reduced_pdb)

    assert_no_metal_in_peptidic(reduced_pdb)

    # (10) Move to receptor and re-fix (post-step edits)
    shutil.copyfile(reduced_pdb, receptor_pdb)
    _helium_postwrite_counter("promote_receptor_copy", receptor_pdb)

    fix_pdb_elements(receptor_pdb)
    _helium_postwrite_counter("elemfix_final_receptor", receptor_pdb)

    quick_element_histogram(receptor_pdb)
    assert file_contains_hydrogens(receptor_pdb), f"[FATAL] Cleaned file lost hydrogens: {receptor_pdb}"

    logging.info("Cleaned receptor: %s", receptor_pdb)
    print(f"[proteinprep] cleaned receptor exists={Path(receptor_pdb).is_file()} -> {receptor_pdb}")
    return str(receptor_pdb)


from collections import defaultdict

def _classify_and_rename_histidines(pdb_in: Union[str, Path], pdb_out: Union[str, Path]) -> None:
    """
    Inspect each HIS residue's side-chain hydrogens and rename:
      - HD1 only  -> HID
      - HE2 only  -> HIE
      - both      -> HIP
      - neither   -> keep HIS (caller may rewrite later)
    """
    pdb_in, pdb_out = str(pdb_in), str(pdb_out)

    residues = defaultdict(list)
    other_lines = []
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                chain = ln[21]
                resseq = ln[22:26]
                icode  = ln[26]
                resn   = ln[17:20]
                key = (chain, resseq, icode, resn)
                residues[key].append(ln)
            else:
                other_lines.append(ln)

    out_lines = []
    for key, atoms in residues.items():
        chain, resseq, icode, resn = key
        resname = resn.strip()
        if resname != "HIS":
            out_lines.extend(atoms)
            continue

        names = {ln[12:16].strip().upper() for ln in atoms}
        has_hd1 = "HD1" in names  # proton on ND1
        has_he2 = "HE2" in names  # proton on NE2

        if has_hd1 and has_he2:
            new = "HIP"
        elif has_hd1 and not has_he2:
            new = "HID"
        elif has_he2 and not has_hd1:
            new = "HIE"
        else:
            new = None

        if new:
            for ln in atoms:
                if ln.startswith(("ATOM  ", "HETATM")):
                    out_lines.append(ln[:17] + f"{new:>3}" + ln[20:])
                else:
                    out_lines.append(ln)
        else:
            out_lines.extend(atoms)

    with open(pdb_out, "w", encoding="utf-8") as w:
        for ln in other_lines:
            w.write(ln)
        for ln in out_lines:
            w.write(ln)


def _rewrite_his_default(pdb_in: Union[str, Path], pdb_out: Union[str, Path], default: str = "HIE") -> None:
    pdb_in, pdb_out = str(pdb_in), str(pdb_out)
    default = default.upper()
    assert default in {"HIE","HID","HIP"}
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as f, \
         open(pdb_out, "w", encoding="utf-8") as w:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")) and ln[17:20] == "HIS":
                w.write(ln[:17] + default + ln[20:])
            else:
                w.write(ln)

# =============================
# External Tools (Meeko/ADT, Phenix, OpenBabel, Reduce)
# =============================
def _drop_free_ions_for_meeko(
    pdb_in: str | Path,
    pdb_out: str | Path,
    banlist: set[str] | None = None,
    *,
    cfg: Optional[dict] = None,
    variant: Optional[str] = None,
    pocket_center: Optional[tuple[float, float, float]] = None,
) -> int:
    """
    Remove HETATM entries for simple ions that Meeko chokes on (e.g., Na, K, Li).
    Writes to pdb_out. Returns #atoms dropped.
    Never drops ions explicitly retained by YAML (retain_in_receptor_resnames) that are elemental tokens.
    Set MEEKO_DROP_VERBOSE=1 to log each dropped residue position.
    """
    cfg_obj: Optional[dict] = None
    if isinstance(cfg, dict):
        cfg_obj = cfg
    elif isinstance(config, dict):
        cfg_obj = config

    allow_tokens, _ = _load_retain_allowlist(cfg_obj)
    allow_set = {str(tok).strip().upper() for tok in allow_tokens if str(tok).strip()}
    retain_ions = {r for r in _RETAIN if _is_element_token(r)}
    base_ban = set(banlist) if banlist else _MEEKO_DROP_IONS

    policy = _normalize_ion_policy(cfg_obj)
    variant_token = _resolve_variant_token(cfg_obj, variant)
    variant_label = variant_token or "legacy"

    radius_cfg = 0.0
    if cfg_obj is not None:
        try:
            radius_cfg = float(cfg_obj.get("HOLO_SALT_STRIP_RADIUS", 0.0) or 0.0)
        except Exception:
            radius_cfg = 0.0
    radius_cfg = max(0.0, radius_cfg)
    radius = radius_cfg if (policy == "by_variant" and variant_token == "HOLO") else 0.0
    radius_term = f"{radius:.2f}" if radius > 0.0 else "none"

    drop_metals = False
    drop_salts = False
    if policy == "always_strip":
        drop_metals = True
        drop_salts = True
    elif policy == "by_variant":
        if variant_token == "HOLO":
            drop_metals = False
            drop_salts = radius > 0.0
        else:
            drop_metals = True
            drop_salts = True
    elif policy == "never_strip":
        drop_metals = False
        drop_salts = False

    if drop_salts and radius > 0.0 and pocket_center is None and variant_token == "HOLO":
        logging.warning(
            "[ions] holo_salt_radius_set_but_no_center stage=meeko action=keep_salts radius=%.2f",
            radius,
        )
        drop_salts = False

    logging.info(
        "[ions.policy] stage=meeko variant=%s drop_metals=%s drop_salts=%s radius=%s allowlist=%d",
        variant_label,
        drop_metals,
        drop_salts,
        radius_term,
        len(allow_set),
    )

    verbose = (os.environ.get("MEEKO_DROP_VERBOSE", "0") not in ("0", "false", "False"))

    dropped = 0
    kept_lines: list[str] = []
    warn_missing_center = False
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if not ln.startswith("HETATM"):
                kept_lines.append(ln)
                continue
            res = ln[17:20].strip().upper()
            elem = (ln[76:78].strip() or res).upper()
            token = res or elem
            if token in allow_set or elem in allow_set:
                kept_lines.append(ln)
                continue

            is_metal = res in _METAL_RESNAMES or elem in _METAL_RESNAMES
            is_salt = res in _SALT_RESNAMES or elem in _SALT_RESNAMES or res in base_ban or elem in base_ban

            remove = False
            reason = ""
            if policy == "never_strip":
                remove = False
            elif is_metal and drop_metals and (token not in retain_ions):
                remove = True
                reason = "meeko_strip_metal"
            elif is_salt and drop_salts:
                if radius > 0.0:
                    dist = _distance_from_center(ln, pocket_center)
                    if dist is None:
                        warn_missing_center = True
                        remove = False
                    elif dist >= radius:
                        remove = True
                        reason = "meeko_salt_far"
                elif (res in base_ban or elem in base_ban) and policy != "never_strip":
                    remove = True
                    reason = "meeko_salt_policy"
            elif (res in base_ban or elem in base_ban) and policy == "always_strip":
                remove = True
                reason = "meeko_policy"

            if remove:
                dropped += 1
                if verbose:
                    resi = ln[22:26].strip()
                    chain = ln[21]
                    logging.info(
                        "[meeko drop] res=%s chain=%s resi=%s elem=%s reason=%s",
                        res,
                        chain,
                        resi,
                        elem or res,
                        reason or "policy",
                    )
                continue

            kept_lines.append(ln)

    if warn_missing_center and radius > 0.0 and variant_token == "HOLO":
        logging.warning(
            "[ions] holo_salt_radius_set_but_no_center stage=meeko action=keep_salts radius=%.2f",
            radius,
        )

    Path(pdb_out).write_text("".join(kept_lines), encoding="utf-8")
    logging.info(
        "[ions.touch] stage=meeko variant=%s dropped=%d kept=%d file=%s",
        variant_label,
        dropped,
        len(kept_lines),
        pdb_out,
    )
    if dropped:
        logging.warning(
            "Meeko pre-sanitize: dropped %d monoatomics (ban=%s).",
            dropped,
            ",".join(sorted(base_ban - retain_ions)),
        )
    return dropped



def run_prepare_receptor(input_pdb: Union[str, Path], output_pdbqt: Union[str, Path], cfg: dict) -> bool:
    """
    Robust Meeko/ADT wrapper with:
      â€¢ Early check for HIS tautomer ambiguity on the modern (--read_pdb) call
      â€¢ Auto -n mapping for all truly ambiguous HIS (no HD1/HE2)
      â€¢ Optional -a (allow_bad_res) retry on template-mismatch
      â€¢ Heavy-handed HIS?<default> fallback
      â€¢ ADT fallback
    """
    import re
    from tempfile import NamedTemporaryFile

    input_pdb = str(input_pdb)
    output_pdbqt = str(output_pdbqt)
    variant_token = _resolve_variant_token(cfg, cfg.get("_CURRENT_VARIANT"))
    variant_label = variant_token or "legacy"
    # Persist all attempt logs here (processed_pdbs/<PDB>/work) NOTE DIFF FROM _WORK_DIR
    work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    # --- Helium preflight on the EXACT file we’re about to feed into Meeko pipeline
    _work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    skip_meeko = False
    try:
        input_pdb = str(_meeko_preflight_or_fail(input_pdb, _work_dir))
    except RuntimeError as e:
        if str(e) == "helium_preflight_failed":
            logging.error("Helium persists pre-Meeko; skipping Meeko and trying ADT fallback.")
            skip_meeko = True
        else:
            raise


    his_default = str(cfg.get("HIS_DEFAULT", "HIE")).upper()
    if his_default not in {"HIE", "HID", "HIP"}:
        his_default = "HIE"

    # improved truthiness parsing
    allow_bad_res = str(cfg.get("MEEKO_ALLOW_BAD_RES", "true")).lower() in ("1", "true", "yes")
    default_altloc = (cfg.get("MEEKO_DEFAULT_ALTLOC") or "").strip()  # e.g. "A" or ""

    def _meeko_cmd() -> list[str]:
        """
        Build a Meeko CLI command that is portable across machines/environments.
        Prefer invoking by module with the current interpreter to avoid PATH/script shims.
        """
        import sys
        # First, try module-based invocation (preferred, shim-free)
        try:
            import meeko  # noqa: F401
            return [sys.executable, "-m", "meeko.cli.mk_prepare_receptor"]
        except Exception:
            pass

        # Fallbacks (still avoid PATH hassles as much as possible)
        import shutil, sysconfig, os
        # 1) PATH (with/without .py)
        for name in ("mk_prepare_receptor", "mk_prepare_receptor.py"):
            exe = shutil.which(name)
            if exe:
                return [exe]
        # 2) Scripts dir of the current interpreter
        try:
            scripts = sysconfig.get_path("scripts")
            cand = os.path.join(scripts, "mk_prepare_receptor.py")
            if os.path.exists(cand):
                return [sys.executable, cand]  # run via python for consistency
        except Exception:
            pass

        raise FileNotFoundError("Meeko CLI not found via module or script.")

    def _run(cmd: List[str]) -> subprocess.CompletedProcess:
        logging.info("Meeko: %s", " ".join(cmd))
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        (logging.info if r.returncode == 0 else logging.warning)(
            "rc=%s\nSTDOUT:\n%s\nSTDERR:\n%s", r.returncode, (r.stdout or ""), (r.stderr or "")
        )
        return r

    def _build_his_override_mapping(pdb_path: str) -> str:
        """Return a single -n mapping like 'A:27,A:94=HIE' for HIS with no HD1/HE2."""
        from collections import defaultdict
        by_res = defaultdict(set)
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                if ln[17:20] != "HIS":
                    continue
                chain = ln[21]
                resi = ln[22:26].strip()
                aname = ln[12:16].strip().upper()
                by_res[(chain, resi)].add(aname)
        targets = []
        for (chain, resi), names in by_res.items():
            if ("HD1" not in names) and ("HE2" not in names):
                try:
                    rnum = int(resi)
                except Exception:
                    rnum = int(resi.strip() or "0")
                targets.append(f"{chain}:{rnum}")
        return (",".join(sorted(targets)) + f"={his_default}") if targets else ""

    def _modern_meeko(pdb_path: str,
                      extra_flags: Optional[List[str]] = None,
                      tag: Optional[str] = None,
                      work_dir_for_tag: Optional[Path] = None) -> subprocess.CompletedProcess:
        try:
            cmd = _meeko_cmd() + ["--read_pdb", pdb_path, "-p", output_pdbqt]
        except FileNotFoundError as e:
            from types import SimpleNamespace
            cp = SimpleNamespace(returncode=127, stdout="", stderr=str(e))
            if tag and work_dir_for_tag:
                _persist_subproc(tag, ["<meeko-not-found>"], cp, work_dir_for_tag, Path(output_pdbqt))
            return cp
        flags = extra_flags or []
        cp = _run(cmd + flags)
        if tag and work_dir_for_tag:
            _persist_subproc(tag, cmd + flags, cp, work_dir_for_tag, Path(output_pdbqt))
            # Inline empty-file hint for fast triage
            try:
                if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
                    head, tail = _first_last_lines(work_dir_for_tag / f"{tag}.stderr.txt")
                    logging.warning("[%s] receptor PDBQT is empty. head=%r tail=%r", tag, head, tail)
            except Exception:
                pass
        return cp

    def _legacy_meeko(pdb_path: str,
                      extra_flags: Optional[List[str]] = None,
                      tag: Optional[str] = None,
                      work_dir_for_tag: Optional[Path] = None) -> bool:
        try:
            base = _meeko_cmd()
        except FileNotFoundError:
            if tag and work_dir_for_tag:
                from types import SimpleNamespace
                _persist_subproc(tag, ["<meeko-not-found>"],
                                 SimpleNamespace(returncode=127, stdout="", stderr="mk_prepare_receptor not found"),
                                 work_dir_for_tag, Path(output_pdbqt))
            return False
        flags = extra_flags or []

        # try '-i' form
        cmd_i = base + ["-i", pdb_path, "-p", output_pdbqt] + flags
        r1 = _run(cmd_i)
        if tag and work_dir_for_tag:
            _persist_subproc(f"{tag}_i", cmd_i, r1, work_dir_for_tag, Path(output_pdbqt))
        if r1.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0:
            return True

        # try '-r/-o' form
        cmd_r = base + ["-r", pdb_path, "-o", output_pdbqt] + flags
        r2 = _run(cmd_r)
        if tag and work_dir_for_tag:
            _persist_subproc(f"{tag}_ro", cmd_r, r2, work_dir_for_tag, Path(output_pdbqt))
            try:
                if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
                    head, tail = _first_last_lines(work_dir_for_tag / f"{tag}_ro.stderr.txt")
                    logging.warning("[%s_ro] receptor PDBQT is empty. head=%r tail=%r", tag, head, tail)
            except Exception:
                pass
        return r2.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0

    # --- Step 0: classify HIS by existing hydrogens (best-case, no coord change)
    with NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp1:
        tmp1_path = tmp1.name
    try:
        _classify_and_rename_histidines(input_pdb, tmp1_path)  # leaves non-diagnostic HIS as HIS
    except Exception as e:
        logging.warning("HIS classify/rename step failed (continuing with original): %s", e)
        tmp1_path = input_pdb

    if _cfg_bool("MEEKO_DROP_FREE_IONS", True):
        logging.info(
            "[ions.touch] stage=meeko-pre variant=%s action=pre_sanitize file=%s",
            variant_label,
            tmp1_path,
        )
        # Never drop YAML-retained elemental ions (banlist computed inside)
        _drop_free_ions_for_meeko(
            tmp1_path,
            tmp1_path,
            cfg=cfg,
            variant=variant_token,
            pocket_center=None,
        )

    # Preflight again on the HIS-classified file we’re actually handing to Meeko now
    if not skip_meeko:
        try:
            _meeko_preflight_or_fail(tmp1_path, Path(output_pdbqt).resolve().parent.parent / "work")
        except RuntimeError as e:
            if str(e) == "helium_preflight_failed":
                logging.error("Helium persists after HIS/ion tweaks; skipping Meeko and trying ADT fallback.")
                skip_meeko = True
            else:
                raise
    # define once near the top of run_prepare_receptor, right after output_pdbqt:
    work_dir = Path(output_pdbqt).resolve().parent.parent / "work"

    # Build explicit commands for each attempt
    try:
        base_meeko = _meeko_cmd()
    except FileNotFoundError as e:
        base_meeko = None
        logging.warning("Meeko CLI not found at prebuild stage: %s", e)

    # Meeko modern (--read_pdb) uses the HIS-classified, ion-sanitized tmp1_path
    if base_meeko:
        meeko_cmd = base_meeko + ["--read_pdb", tmp1_path, "-p", output_pdbqt]
        # Legacy single-shot form (we still capture rc/file size for logging)
        meeko_cmd_legacy = base_meeko + ["-r", tmp1_path, "-o", output_pdbqt]
    else:
        meeko_cmd = None
        meeko_cmd_legacy = None
    # Only attempt ADT if BOTH are explicitly configured and present
    if not (MGLTOOLS_PYTHON and Path(MGLTOOLS_PYTHON).is_file() and os.access(MGLTOOLS_PYTHON, os.X_OK) and
            PREPARE_RECEPTOR_SCRIPT and Path(PREPARE_RECEPTOR_SCRIPT).is_file()):
        logging.info("[receptor] ADT fallback disabled (missing MGLTOOLS_PYTHON or PREPARE_RECEPTOR_SCRIPT)")
        return False  # keep modern Meeko as the only path unless ADT is truly configured

    # --- Step 3: ADT prepare_receptor4 fallback (only if both keys are valid)
    adt_ok = False
    mgltools_python = MGLTOOLS_PYTHON
    prepare_script = PREPARE_RECEPTOR_SCRIPT
    if (mgltools_python and Path(mgltools_python).is_file() and os.access(mgltools_python, os.X_OK)
            and prepare_script and Path(prepare_script).is_file()):
        allow_tokens, _ = _load_retain_allowlist(cfg)
        logging.info(
            "[ions.policy] stage=adt variant=%s drop_metals=%s drop_salts=%s radius=none allowlist=%d",
            variant_label,
            False,
            False,
            len({str(tok).strip().upper() for tok in allow_tokens if str(tok).strip()}),
        )
        adt_cmd = [mgltools_python, prepare_script, "-r", tmp1_path, "-o", output_pdbqt,
                   "-A", "none", "-U", "nphs_lps_nonstdres"]
        cp = subprocess.run(adt_cmd, capture_output=True, text=True)
        _persist_subproc("adt_prepare_receptor4", adt_cmd, cp, work_dir, Path(output_pdbqt))
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            logging.info("Prepared receptor with ADT prepare_receptor4.py.")
            return True
        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(work_dir / "adt_prepare_receptor4.stderr.txt")
            logging.warning("[adt_prepare_receptor4] receptor PDBQT is empty. head=%r tail=%r", head, tail)
    else:
        logging.info("[receptor] ADT fallback disabled (MGLTOOLS_PYTHON/PREPARE_RECEPTOR_SCRIPT not set)")
    return False

    # --- Step 1: Modern Meeko attempt
    if not skip_meeko and meeko_cmd:
        cp = subprocess.run(meeko_cmd, capture_output=True, text=True)
        _persist_subproc("meeko_modern", meeko_cmd, cp, work_dir, Path(output_pdbqt))
    else:
        # If Meeko is skipped or unavailable, synthesize a "failed" result object
        from types import SimpleNamespace
        cp = SimpleNamespace(returncode=127, stdout="", stderr=("meeko_skipped" if skip_meeko else "meeko_not_found"))

    if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
        logging.info("Prepared receptor PDBQT with modern Meeko.")
        # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
        try:
            # Compute present/missing sets in-line (keeps compare logic decoupled from warn-only helper)
            def _ions_in_pdb(p):
                s = set()
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    for _ln in fh:
                        if not _ln.startswith("HETATM"): continue
                        _res = _ln[17:20].strip().upper();
                        _c = _ln[21];
                        _i = _ln[22:26].strip()
                        if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                return s

            def _ions_in_pdbqt(p):
                s = set()
                if not Path(p).exists(): return s
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    for _ln in fh:
                        if not _ln.startswith(("ATOM  ", "HETATM")): continue
                        _res = _ln[17:20].strip().upper();
                        _c = _ln[21];
                        _i = _ln[22:26].strip()
                        if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                return s

            _pdb_ions = _ions_in_pdb(input_pdb)
            _pdbqt_ions = _ions_in_pdbqt(output_pdbqt)
            _missing = _pdb_ions - _pdbqt_ions
            if _pdb_ions or _missing:
                logging.warning("[ion diff] present_pdb=%s missing_in_pdbqt=%s",
                                sorted(_pdb_ions), sorted(_missing))
            # If any missing ion is YAML-retained, emit a targeted warning
            _retained = {r for r in _RETAIN if _is_element_token(r)}
            _lost_retained = sorted([x for x in _missing if x[0] in _retained])
            if _lost_retained:
                logging.warning("[ion lost] %s", _lost_retained)
        except Exception as _e:
            logging.warning("[ion diff] skipped note=%s", _e)
        return True

    if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
        head, tail = _first_last_lines(work_dir / "meeko_modern.stderr.txt")
        logging.warning("[meeko_modern] receptor PDBQT is empty. head=%r tail=%r", head, tail)

    # Decide if we need an -n retry (histidine tie) based on se0
    se0 = cp.stderr or ""
    his_tie = ("tied for fewest missing h" in se0.lower()) and ("hie" in se0.lower() and "hid" in se0.lower())
    mapping = ""
    if his_tie:
        mapping = _build_his_override_mapping(tmp1_path)
        if not mapping:
            m = re.search(r"residue_key='([A-Za-z]):(\d+)'", se0)
            if m:
                mapping = f"{m.group(1)}:{int(m.group(2))}={his_default}"

    # --- Step 1b: Retry with -n if needed
    if mapping and meeko_cmd:
        meeko_cmd_with_n = meeko_cmd + ["-n", mapping]
        cp = subprocess.run(meeko_cmd_with_n, capture_output=True, text=True)
        _persist_subproc("meeko_retry_nmap", meeko_cmd_with_n, cp, work_dir, Path(output_pdbqt))
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            logging.info("Prepared receptor after HIS -n mapping.")
            # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
            try:
                def _ions_in_pdb(p):
                    s = set()
                    with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                        for _ln in fh:
                            if not _ln.startswith("HETATM"): continue
                            _res = _ln[17:20].strip().upper();
                            _c = _ln[21];
                            _i = _ln[22:26].strip()
                            if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                    return s

                def _ions_in_pdbqt(p):
                    s = set()
                    if not Path(p).exists(): return s
                    with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                        for _ln in fh:
                            if not _ln.startswith(("ATOM  ", "HETATM")): continue
                            _res = _ln[17:20].strip().upper();
                            _c = _ln[21];
                            _i = _ln[22:26].strip()
                            if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                    return s

                _pdb_ions = _ions_in_pdb(input_pdb)
                _pdbqt_ions = _ions_in_pdbqt(output_pdbqt)
                _missing = _pdb_ions - _pdbqt_ions
                if _pdb_ions or _missing:
                    logging.warning("[ion diff] present_pdb=%s missing_in_pdbqt=%s",
                                    sorted(_pdb_ions), sorted(_missing))
                _retained = {r for r in _RETAIN if _is_element_token(r)}
                _lost_retained = sorted([x for x in _missing if x[0] in _retained])
                if _lost_retained:
                    logging.warning("[ion lost] %s", _lost_retained)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            return True

        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(work_dir / "meeko_retry_nmap.stderr.txt")
            logging.warning("[meeko_retry_nmap] receptor PDBQT is empty. head=%r tail=%r", head, tail)

    # --- Step 1c: Retry with -a allow_bad_res (and default_altloc if hinted)
    templ_fail = ("no template matched for residue_key" in se0.lower()) or (
                "template matching failed" in se0.lower()) or ("template matched failed" in se0.lower())
    extra = []
    if meeko_cmd and allow_bad_res and (
            templ_fail or "allow_bad_res" in se0.lower() or "recommendations" in se0.lower()):
        extra = ["-a"]
        if (not default_altloc) and ("default_altloc" in se0.lower() or "altloc" in se0.lower()):
            default_altloc = "A"
            logging.warning("Assuming --default_altloc A based on Meeko hint.")
        if default_altloc:
            extra += ["--default_altloc", default_altloc]
        meeko_cmd_allow_bad = meeko_cmd + extra
        cp = subprocess.run(meeko_cmd_allow_bad, capture_output=True, text=True)
        _persist_subproc("meeko_retry_allow_bad_res", meeko_cmd_allow_bad, cp, work_dir, Path(output_pdbqt))
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            logging.info("Prepared receptor after -a allow_bad_res.")
            # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
            try:
                def _ions_in_pdb(p):
                    s = set()
                    with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                        for _ln in fh:
                            if not _ln.startswith("HETATM"): continue
                            _res = _ln[17:20].strip().upper();
                            _c = _ln[21];
                            _i = _ln[22:26].strip()
                            if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                    return s

                def _ions_in_pdbqt(p):
                    s = set()
                    if not Path(p).exists(): return s
                    with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                        for _ln in fh:
                            if not _ln.startswith(("ATOM  ", "HETATM")): continue
                            _res = _ln[17:20].strip().upper();
                            _c = _ln[21];
                            _i = _ln[22:26].strip()
                            if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                    return s

                _pdb_ions = _ions_in_pdb(input_pdb)
                _pdbqt_ions = _ions_in_pdbqt(output_pdbqt)
                _missing = _pdb_ions - _pdbqt_ions
                if _pdb_ions or _missing:
                    logging.warning("[ion diff] present_pdb=%s missing_in_pdbqt=%s",
                                    sorted(_pdb_ions), sorted(_missing))
                _retained = {r for r in _RETAIN if _is_element_token(r)}
                _lost_retained = sorted([x for x in _missing if x[0] in _retained])
                if _lost_retained:
                    logging.warning("[ion lost] %s", _lost_retained)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            return True

        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(work_dir / "meeko_retry_allow_bad_res.stderr.txt")
            logging.warning("[meeko_retry_allow_bad_res] receptor PDBQT is empty. head=%r tail=%r", head, tail)
    # --- Step 2: Legacy Meeko
    if not skip_meeko and meeko_cmd_legacy:
        cp = subprocess.run(meeko_cmd_legacy, capture_output=True, text=True)
        _persist_subproc("meeko_legacy", meeko_cmd_legacy, cp, work_dir, Path(output_pdbqt))
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            logging.info("Prepared receptor with legacy Meeko.")
            # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
            try:
                def _ions_in_pdb(p):
                    s = set()
                    with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                        for _ln in fh:
                            if not _ln.startswith("HETATM"): continue
                            _res = _ln[17:20].strip().upper();
                            _c = _ln[21];
                            _i = _ln[22:26].strip()
                            if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                    return s

                def _ions_in_pdbqt(p):
                    s = set()
                    if not Path(p).exists(): return s
                    with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                        for _ln in fh:
                            if not _ln.startswith(("ATOM  ", "HETATM")): continue
                            _res = _ln[17:20].strip().upper();
                            _c = _ln[21];
                            _i = _ln[22:26].strip()
                            if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                    return s

                _pdb_ions = _ions_in_pdb(input_pdb)
                _pdbqt_ions = _ions_in_pdbqt(output_pdbqt)
                _missing = _pdb_ions - _pdbqt_ions
                if _pdb_ions or _missing:
                    logging.warning("[ion diff] present_pdb=%s missing_in_pdbqt=%s",
                                    sorted(_pdb_ions), sorted(_missing))
                _retained = {r for r in _RETAIN if _is_element_token(r)}
                _lost_retained = sorted([x for x in _missing if x[0] in _retained])
                if _lost_retained:
                    logging.warning("[ion lost] %s", _lost_retained)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            return True

        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(work_dir / "meeko_legacy.stderr.txt")
            logging.warning("[meeko_legacy] receptor PDBQT is empty. head=%r tail=%r", head, tail)

    # --- Step 3: ADT prepare_receptor4 fallback
    cp = subprocess.run(adt_cmd, capture_output=True, text=True)
    _persist_subproc("adt_prepare_receptor4", adt_cmd, cp, work_dir, Path(output_pdbqt))
    if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
        logging.info("Prepared receptor with ADT prepare_receptor4.py.")
        # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
        try:
            def _ions_in_pdb(p):
                s = set()
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    for _ln in fh:
                        if not _ln.startswith("HETATM"): continue
                        _res = _ln[17:20].strip().upper();
                        _c = _ln[21];
                        _i = _ln[22:26].strip()
                        if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                return s

            def _ions_in_pdbqt(p):
                s = set()
                if not Path(p).exists(): return s
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    for _ln in fh:
                        if not _ln.startswith(("ATOM  ", "HETATM")): continue
                        _res = _ln[17:20].strip().upper();
                        _c = _ln[21];
                        _i = _ln[22:26].strip()
                        if _is_element_token(_res): s.add((_res, f"{_c}:{_i}"))
                return s

            _pdb_ions = _ions_in_pdb(input_pdb)
            _pdbqt_ions = _ions_in_pdbqt(output_pdbqt)
            _missing = _pdb_ions - _pdbqt_ions
            if _pdb_ions or _missing:
                logging.warning("[ion diff] present_pdb=%s missing_in_pdbqt=%s",
                                sorted(_pdb_ions), sorted(_missing))
            _retained = {r for r in _RETAIN if _is_element_token(r)}
            _lost_retained = sorted([x for x in _missing if x[0] in _retained])
            if _lost_retained:
                logging.warning("[ion lost] %s", _lost_retained)
        except Exception as _e:
            logging.warning("[ion diff] skipped note=%s", _e)
        return True

    if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
        head, tail = _first_last_lines(work_dir / "adt_prepare_receptor4.stderr.txt")
        logging.warning("[adt_prepare_receptor4] receptor PDBQT is empty. head=%r tail=%r", head, tail)



def _phenix_detect() -> tuple[str, list[str] | None, dict]:
    """
    Return (mode, cmd_list, env) for phenix.pdbtools detection.
    mode ∈ {'linux-path','linux-explicit','wsl-bat','unavailable'}
    """
    env = dict(os.environ)
    if PHENIX_LIB_PATH:
        env["PYTHONPATH"] = PHENIX_LIB_PATH

    # 1) On PATH (native Linux installs)
    if shutil.which("phenix.pdbtools"):
        return ("linux-path", ["phenix.pdbtools"], env)

    # 2) Explicit Linux path via PHENIX_DIR (your case)
    if PHENIX_DIR:
        explicit = Path(PHENIX_DIR) / "phenix.pdbtools"
        if explicit.exists() and os.access(str(explicit), os.X_OK):
            return ("linux-explicit", [str(explicit)], env)

    # 3) Windows .bat (usable from WSL via PowerShell)
    if PDBTOOLS_BAT and Path(PDBTOOLS_BAT).exists():
        return ("wsl-bat", [PDBTOOLS_BAT], env)

    return ("unavailable", None, env)


def _run_and_log(cmd: list[str], *, env: dict, check: bool = False) -> subprocess.CompletedProcess:
    logging.info("[phenix] exec: %s", " ".join(cmd))
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if cp.stdout:
        logging.debug("[phenix][stdout]\n%s", cp.stdout)
    if cp.stderr:
        logging.debug("[phenix][stderr]\n%s", cp.stderr)
    if check and cp.returncode != 0:
        raise subprocess.CalledProcessError(cp.returncode, cmd, cp.stdout, cp.stderr)
    return cp


def run_phenix_pdbtools(input_pdb: Union[str, Path],
                        output_pdb: Union[str, Path],
                        remove_waters: bool = True) -> bool:
    """
    Run Phenix pdbtools if available (Linux, or Windows via PowerShell under WSL).
    Return True if successful. Always log branch as: [phenix] mode=...
    """
    inp  = os.path.abspath(str(input_pdb)).replace("\\", "/")
    outp = os.path.abspath(str(output_pdb)).replace("\\", "/")

    mode, base_cmd, env = _phenix_detect()
    logging.info("[phenix] mode=%s", mode)

    if mode == "unavailable" or not base_cmd:
        logging.warning("Phenix not available; skipping pdbtools polish")
        return False

    cmd = base_cmd + [inp, f"output.file_name={outp}"]
    if remove_waters:
        cmd.append('remove="resname HOH"')

    try:
        cp = _run_and_log(cmd, env=env, check=True)
        if os.path.exists(outp) and os.path.getsize(outp) > 0:
            logging.info("phenix.pdbtools wrote %s", outp)
            return True
        logging.warning("phenix.pdbtools finished but did not create expected output: %s", outp)
        return False
    except Exception as e:
        logging.error("phenix.pdbtools failed: %s", e)
        # Diagnostic: phenix.python import iotbx (same env)
        py_cmd = None
        if PHENIX_DIR and (Path(PHENIX_DIR) / "phenix.python").exists():
            py_cmd = [str(Path(PHENIX_DIR) / "phenix.python")]
        elif shutil.which("phenix.python"):
            py_cmd = ["phenix.python"]
        elif PHENIX_PYTHON_BAT and Path(PHENIX_PYTHON_BAT).exists():
            py_cmd = [PHENIX_PYTHON_BAT]
        if py_cmd:
            try:
                _run_and_log(py_cmd + ["-c", "from iotbx import pdb; print('iotbx_ok')"], env=env, check=False)
            except Exception:
                pass
        return False


def run_windows_phenix_clean_script(loop_fixed_pdb: Union[str, Path], nolig_dir: Union[str, Path]) -> int:
    """Call your PHENIX_CLEAN_SCRIPT using Windows Phenix Python from WSL. Returns returncode."""
    if not PHENIX_PYTHON_BAT or not Path(PHENIX_PYTHON_BAT).exists() or not PHENIX_CLEAN_SCRIPT:
        return 127
    win_script = _win_path(PHENIX_CLEAN_SCRIPT)
    win_in = _win_path(loop_fixed_pdb)
    win_outdir = _win_path(nolig_dir)
    ps = (
        "$ErrorActionPreference='Stop';"
        f"if (!(Test-Path '{win_outdir}')) {{ New-Item -ItemType Directory -Force -Path '{win_outdir}' | Out-Null }};"
        f"Push-Location '{win_outdir}';"
        f"& '{_win_path(PHENIX_PYTHON_BAT)}' '{win_script}' '{win_in}' '{win_outdir}';"
        "Pop-Location;"
    )
    r = _powershell(ps)
    return r.returncode

def run_openbabel_add_h(input_pdb: Union[str, Path], output_pdb: Union[str, Path]) -> None:
    obabel = OPENBABEL_PATH or shutil.which("obabel")
    if not obabel or not shutil.which(Path(obabel).name):
        raise RuntimeError("Open Babel not found. Install or set OPENBABEL_PATH.")
    cmd = [obabel, str(input_pdb), "-O", str(output_pdb), "-h"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Open Babel failed:\n{r.stderr}")
    logging.info("Open Babel hydrogenation succeeded")
from typing import Union, Optional
from pathlib import Path
import shutil, subprocess, logging, os


def assign_protonation_states(input_pdb: Union[str, Path],
                              output_pdb: Union[str, Path],
                              reduce_exe: Optional[str] = None) -> str:
    import uuid
    from pathlib import Path
    # --- toggles ---
    _REDUCE_RETRY   = str(os.environ.get("REDUCE_RETRY", "1")).lower() not in {"0","false","no"}
    _FALLBACK_ADDH  = str(os.environ.get("FALLBACK_ADDH", "1")).lower() not in {"0","false","no"}

    input_pdb  = os.path.abspath(str(input_pdb)).replace("\\", "/")
    output_pdb = os.path.abspath(str(output_pdb)).replace("\\", "/")

    # Optional precheck: very low CONECT coverage is a strong predictor of Reduce complaints.
    try:
        cc = conect_coverage(input_pdb)
        if cc < 0.05:
            logging.info("[reduce precheck] low_conect=%.2f file=%s", cc, input_pdb)
            try:
                # Idempotent, YAML-driven element repair
                fix_pdb_elements(input_pdb)
            except Exception as _e:
                logging.warning("[reduce precheck] elemfix_skip note=%s", _e)
    except Exception as _e:
        logging.warning("[reduce precheck] coverage_skip note=%s", _e)

    has_h = hydrogenation_status(input_pdb)[0] != "NO_H"
    # If input already has H (e.g., from PDB2PQR), do flip/cleanup only; else do a full H build.
    reduce_flags = ["-FLIP", "-Quiet"] if has_h else ["-BUILD", "-Quiet"]

    exe = reduce_exe
    exe_dir = os.path.dirname(exe) if exe else None

    def run_reduce(in_pdb: str, stage_name: str) -> str:
        env = dict(os.environ)
        het = _het_dict_path()
        if het:
            env["REDUCE_HET_DICT"] = het
        with open(output_pdb, "w", encoding="utf-8") as out:
            cp = subprocess.run([exe] + reduce_flags + [in_pdb], stdout=out, stderr=subprocess.PIPE, text=True, env=env)
        # Persist stderr and count atoms written
        try:
            (Path(output_pdb).parent / f"{stage_name}.stderr.txt").write_text(cp.stderr or "", encoding="utf-8")
        except Exception:
            pass
        wrote_atoms = 0
        try:
            with open(output_pdb, "r", encoding="utf-8", errors="ignore") as fh:
                wrote_atoms = sum(1 for ln in fh if ln.startswith(("ATOM  ", "HETATM")))
        except Exception:
            wrote_atoms = 0
        logging.warning("[reduce] stage=%s rc=%s flags=%s wrote_atoms=%d", stage_name, cp.returncode,
                        " ".join(reduce_flags), wrote_atoms)
        return output_pdb

    status_before, h0, hv0, r0 = hydrogenation_status(input_pdb)
    logging.info("[H-Scan before] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_before, h0, hv0, r0)
    if exe is None:
        # Caller requested to skip Reduce (e.g., nucleotides). Go straight to fallback.
        if _FALLBACK_ADDH:
            logging.info("[protonate] skipping Reduce by request; using OpenBabel fallback")
            run_openbabel_add_h(input_pdb, output_pdb)
        else:
            # Keep as-is if fallback is disabled
            shutil.copy(input_pdb, output_pdb)
            logging.warning("[protonate] Reduce skipped and FALLBACK_ADDH=0; copying input")
        # Continue into the existing post-guard checks (size/H-Scan) below
        sz = Path(output_pdb).stat().st_size if Path(output_pdb).exists() else 0
        if sz == 0:
            raise RuntimeError("protonation_empty_output")

        #  strict ATOM/HETATM guard (Reduce can write tiny, header-only files)
        def _file_has_atoms(p: str) -> bool:
            try:
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    return any(ln.startswith(("ATOM  ", "HETATM")) for ln in fh)
            except Exception:
                return False

        if not _file_has_atoms(output_pdb):
            logging.error("[protonate] wrote_atoms=0 after Reduce; attempting OpenBabel fallback")
            try:
                tmp_babel = output_pdb + ".babel.pdb"
                run_openbabel_add_h(input_pdb, tmp_babel)
                shutil.move(tmp_babel, output_pdb)
                logging.warning("[fallback] openbabel applied (post-Reduce)")
            except Exception as e:
                logging.error("[fallback] openbabel failed: %s; copying input→output", str(e)[:200])
                shutil.copy(input_pdb, output_pdb)
        status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
        logging.info("[H-Scan after ] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_after, h1, hv1, r1)
        return output_pdb

    # --- First attempt (default flags) ---
    try:
        run_reduce(input_pdb, "Reduce#1")
    except Exception as e:
        logging.warning("[reduce#1 fail] exe=%s het_dict=%s has_h=%s note=%s",
                        exe, os.environ.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT",""),
                        str(has_h), str(e)[:200])
        completed = False
        # Quick element repair before a retry
        try:
            fix_pdb_elements(input_pdb)
        except Exception:
            pass

        # ---  retry ladder ---
        retried_ok = False
        if _REDUCE_RETRY:
            try:
                # Second attempt: add -noflip (Reduce sometimes flips/complains on tricky H networks)
                reduce_flags_noflip = (["-noflip"] + reduce_flags) if "-noflip" not in reduce_flags else reduce_flags
                def run_reduce_noflip(in_pdb: str, stage_name: str) -> str:
                    env = dict(os.environ)
                    het = env.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT", "")
                    if het:
                        env["REDUCE_HET_DICT"] = het
                    with open(output_pdb, "w", encoding="utf-8") as out:
                        cp = subprocess.run([exe] + reduce_flags_noflip + [in_pdb],
                                            stdout=out, stderr=subprocess.PIPE, text=True,
                                            cwd=exe_dir, env=env)
                    (logging.info if cp.returncode == 0 else logging.warning)(
                        "[reduce] stage=%s rc=%s flags=%r in=%s", stage_name, cp.returncode, reduce_flags_noflip, in_pdb
                    )
                    if cp.returncode != 0:
                        raise RuntimeError(f"{stage_name} reduce failed: {cp.stderr.strip()}")
                    return output_pdb

                run_reduce_noflip(input_pdb, "Reduce#2(noflip)")
                retried_ok = True
                completed = True
            except Exception as e2:
                logging.warning("[reduce#2 fail] noflip note=%s", str(e2)[:200])

        if not retried_ok:
            # --- H-add fallback path (OpenBabel) ---
            if _FALLBACK_ADDH:
                try:
                    # Guard element columns first to avoid He/column drift in fallback
                    try:
                        fix_pdb_elements(input_pdb)
                    except Exception as _e_fix:
                        logging.warning("[element] guard before fallback failed: %s", _e_fix)
                    logging.info("[protonate fallback] using=OpenBabel in=%s out=%s", input_pdb, output_pdb)
                    run_openbabel_add_h(input_pdb, output_pdb)
                    logging.warning("[reduce] fallback H-add applied")
                    completed = True

                except Exception as babel_error:
                    logging.error("OpenBabel fallback failed: %s", babel_error)
                    # Last resort: if input already had H, keep them; else raise
                    if has_h:
                        shutil.copy(input_pdb, output_pdb)
                        logging.warning("Reduce+fallback failed; keeping existing hydrogens.")
                    else:
                        raise
        if not completed:
            try:
                # Repair element columns (YAML-driven fixer under the hood)
                fix_pdb_elements(input_pdb)
            except Exception:
                pass
            try:
                tmp_in = os.path.splitext(input_pdb)[0] + f"_retry_{uuid.uuid4().hex}.pdb"
                shutil.copy(input_pdb, tmp_in)
                try:
                    run_reduce(tmp_in, "Reduce#3(temp)")
                finally:
                    try:
                        os.remove(tmp_in)
                    except Exception:
                        pass
            except Exception as e2:
                logging.error("Reduce temp retry failed: %s", e2)
                try:
                    if has_h:
                        shutil.copy(input_pdb, output_pdb)
                        logging.warning("Reduce failed; keeping existing hydrogens (no rebuild).")
                    else:
                        logging.info("[protonate fallback] using=OpenBabel in=%s out=%s", input_pdb, output_pdb)
                        run_openbabel_add_h(input_pdb, output_pdb)
                        logging.info("Open Babel used to add hydrogens.")
                except Exception as babel_error:
                    logging.error("OpenBabel fallback failed: %s", babel_error)
                    shutil.copy(input_pdb, output_pdb)
                    logging.warning("Hydrogenation skipped; copied input to output.")

    # Guard: the path we will scan must exist and be non-empty
    sz = Path(output_pdb).stat().st_size if Path(output_pdb).exists() else 0
    if sz == 0:
        raise RuntimeError("protonation_empty_output")

    status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
    logging.info("[H-Scan after ] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_after, h1, hv1, r1)

    if status_after == "NO_H":
        logging.warning("Reduce produced no hydrogens; trying Open Babel fallback.")
        try:
            logging.info("[protonate fallback] using=OpenBabel in=%s out=%s", input_pdb, output_pdb)
            tmp_babel = output_pdb + ".babel.pdb"
            run_openbabel_add_h(input_pdb, tmp_babel)
            shutil.move(tmp_babel, output_pdb)
            status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
            logging.info("[H-Scan babel] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_after, h1, hv1, r1)
        except Exception as e:
            logging.error("OpenBabel fallback failed after Reduce: %s", e)

        if hydrogenation_status(output_pdb)[0] == "NO_H":
            if str(os.environ.get("ALLOW_NO_HYDROGENS", "")).lower() in ("1","true","yes"):
                logging.warning("Continuing with NO_H due to ALLOW_NO_HYDROGENS env override.")
                # Still emit result line for grep
                try:
                    sz = Path(output_pdb).stat().st_size if Path(output_pdb).exists() else 0
                except Exception:
                    sz = 0
                logging.info("[protonate result] status=%s H=%d Heavy=%d ratio=%.2f size=%d",
                             "NO_H", h1, hv1, r1, sz)
                return output_pdb
            raise RuntimeError("Protonation produced no hydrogens.")

    # Final greppable result line
    try:
        sz = Path(output_pdb).stat().st_size if Path(output_pdb).exists() else 0
    except Exception:
        sz = 0
    logging.info("[protonate result] status=%s H=%d Heavy=%d ratio=%.2f size=%d",
                 status_after, h1, hv1, r1, sz)

    return output_pdb




def run_molprobity_validate(pdb_path: Union[str, Path], work_dir: Optional[Union[str, Path]] = None) -> int:
    """
    Try running MolProbity validation via 'phenix.molprobity' if available.
    Returns the process return code (0 means success). Non-fatal if missing.
    """
    exe = shutil.which("phenix.molprobity")
    if not exe:
        logging.info("MolProbity not available; skipping validation.")
        return 127
    work = Path(work_dir or Path(pdb_path).with_suffix(".molprobity")).resolve()
    work.mkdir(parents=True, exist_ok=True)
    cmd = [exe, str(Path(pdb_path).resolve())]
    logging.info("Running MolProbity: %s", " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(work), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        logging.warning("MolProbity failed (rc=%s)\nSTDERR:\n%s", r.returncode, (r.stderr or ""))
    else:
        logging.info("MolProbity finished; results under %s", work)
    return r.returncode
# =============================
# Module Entrypoint
# =============================

def main(pdb_filename: str, output_dir: Union[str, Path] = r"./processed_pdbs"):
    """High-level wrapper: clean a PDB and prepare the receptor PDBQT."""
    try:
        # Resolve input path robustly:
        # 1) If the provided path already points to a file (absolute or relative), use it as-is.
        # 2) Otherwise, fall back to INPUT_DIR/pdb_filename.
        cand = Path(pdb_filename)
        if cand.is_file():
            pdb_path = str(cand.resolve())
        else:
            pdb_path = str(Path(_cfg("INPUT_DIR", ".")).joinpath(pdb_filename))

        logging.info("[prep] argv pdb_filename=%s resolved=%s output_dir=%s cwd=%s",
                     pdb_filename, pdb_path, str(output_dir), os.getcwd())

        if not os.path.isfile(pdb_path):
            logging.error("ERROR: File does not exist: %s", pdb_path)
            raise FileNotFoundError(pdb_path)


        # Run cleaning → returns the final cleaned receptor PDB path
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        cleaned_pdb = clean_pdb(pdb_path, output_dir)
        if not cleaned_pdb:
            logging.error("ERROR: Cleaning failed for %s", pdb_filename)
            raise RuntimeError(f"cleaning_failed: {pdb_filename}")

        # Prepare receptor PDBQT next to the cleaned tree
        raw_stem = Path(pdb_filename).stem
        pdb_id = re.sub(r"(_nolig(_cleaned)?|_cleaned)$", "", raw_stem, flags=re.I).upper()
        logging.info("[prep.id] main stem=%s -> base_id=%s", raw_stem, pdb_id)
        # >>> LEGACY PATHS PATCH START
        legacy_paths = canon_paths(pdb_id, output_dir)
        logging.info("[prep.paths] protein_root=%s receptor=%s nolig=%s work=%s",
                     legacy_paths["protein_root"], legacy_paths["receptor"],
                     legacy_paths["nolig"], legacy_paths["work"])
        # >>> LEGACY PATHS PATCH END
        # >>> PATHS INIT START
        paths = make_paths(config, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        # >>> PATHS INIT END
        # >>> RECEPTOR PATHS PATCH START
        cleaned_pdb_out = str(paths.receptor_cleaned_pdb(None))
        receptor_pdbqt_out = str(paths.receptor_pdbqt(None, ph_token=None))
        # >>> RECEPTOR PATHS PATCH END
        # >>> RECEPTOR OUTPUT PATCH START
        cleaned_target = Path(cleaned_pdb_out)
        cleaned_target.parent.mkdir(parents=True, exist_ok=True)
        if Path(cleaned_pdb).resolve() != cleaned_target.resolve():
            shutil.copyfile(cleaned_pdb, str(cleaned_target))
        cleaned_pdb = str(cleaned_target)
        receptor_target = Path(receptor_pdbqt_out)
        receptor_target.parent.mkdir(parents=True, exist_ok=True)
        output_pdbqt = str(receptor_target)
        # >>> RECEPTOR OUTPUT PATCH END
        fix_element_columns_in_file(cleaned_pdb, cleaned_pdb, rewrite_atoms=True)

        if not run_prepare_receptor(cleaned_pdb, output_pdbqt, config):
            logging.error("ERROR: Failed to prepare receptor PDBQT for %s", pdb_id)
            try:
                rp = Path(output_pdbqt)
                exists = rp.exists()
                size = rp.stat().st_size if exists else 0
            except Exception:
                exists = False
                size = 0
            logging.info("[receptor-summary]\n"
                         "cleaned_pdb=%s\n"
                         "receptor_pdbqt=%s exists=%s size=%d\n"
                         "meeko_attempts=(see work/*.cmd.txt | *.stderr.txt)",
                         cleaned_pdb, output_pdbqt, exists, size)
            raise RuntimeError(f"receptor_pdbqt_failed: {pdb_id}")


        # Success path summary
        try:
            rp = Path(output_pdbqt)
            exists = rp.exists()
            size = rp.stat().st_size if exists else 0
        except Exception:
            exists = False
            size = 0
        logging.info("[receptor-summary]\n"
                     "cleaned_pdb=%s\n"
                     "receptor_pdbqt=%s exists=%s size=%d\n"
                     "meeko_attempts=(see work/*.cmd.txt | *.stderr.txt)",
                     cleaned_pdb, output_pdbqt, exists, size)
        
        logging.info("[prep.return] cleaned=%s receptor_pdbqt=%s", cleaned_pdb, output_pdbqt)
        logging.info("Prepared receptor PDBQT: %s", output_pdbqt)
        return cleaned_pdb, output_pdbqt

    except Exception as e:
        logging.exception("[FATAL] automate_protein_prep.main() failed: %s", e)
        raise



