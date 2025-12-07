import os, sys, shutil, json
from pathlib import Path
from distutils.util import strtobool
from typing import Any, Dict, Optional
import re, csv, math
import pandas as pd
from collections import defaultdict
# >>> PATHS IMPORT START
from path_router import (
    make_paths,
    config_dir as router_config_dir,
    config_file as router_config_file,
    docked_dir as router_docked_dir,
    receptor_file as router_receptor_file,
)
# >>> PATHS IMPORT END

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
    return Path(os.environ.get("PROTEIN_AUTOMATION_DIR", Path(__file__).resolve().parent))

def _default_dirs(base: Path) -> Dict[str, str]:
    """
    Opinionated defaults that match BRCF layout:
      /stor/home/<user>/atlas/code/protein_automation/{subdirs}
    """
    return {
        "OVERALL_DIR":            str(base),
        "INPUT_DIR":              str(base / "input_pdbs"),
        "PROTEIN_DIR":            str(base / "pdbqts"),
        "LIGAND_DIR":             str(base / "input_ligands"),
        "LIGAND_EXTRACTED_DIR":   str(base / "extracted_ligands"),
        "LIGANDS_MOL2_DIR":       str(base / "ligands_mol2"),
        "OUTPUT_LIGANDS_DIR":     str(base / "prepped_ligands"),
        "OUTPUT_DIR":             str(base / "processed_pdbs"),
        "PDBQT_DIR":              str(base / "pdbqts"),
        "DOCKED_DIR":             str(base / "docked"),
        # p2rank
        "P2RANK_OUTPUT_DIR":      str(base / "p2rank_out"),
        # cleanup script (kept in repo)
        "PHENIX_CLEAN_SCRIPT":    str(base / "phenix_clean.py"),
    }

def _default_tools() -> Dict[str, str]:
    """
    Tool defaults prefer environment & PATH.
    """
    return {
        "VINA_PATH": which_or_exists(["vina"]),
        "VINA_EXE": which_or_exists(["vina"]),
        "OPENBABEL_PATH": which_or_exists(["obabel"]),
        "PYMOL_PATH": which_or_exists(["pymol"]),
        "REDUCE_EXE": which_or_exists(["reduce"]),
        # prefer env MGLTOOLS_* if set; otherwise try ADFRsuite pythonsh on PATH
        "MGLTOOLS_PYTHON": os.environ.get("MGL_PYTHON") or which_or_exists(["pythonsh"]),
        "PREPARE_LIGAND_SCRIPT": os.environ.get("PREPARE_LIGAND_SCRIPT") or which_or_exists(["prepare_ligand4.py"]),
        "PREPARE_RECEPTOR_SCRIPT": os.environ.get("PREPARE_RECEPTOR_SCRIPT") or which_or_exists(["prepare_receptor4.py"]),
        # p2rank: either absolute prank or found on PATH
        "P2RANK_PATH": shutil.which("prank") or "prank",
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
        # control-centering knobs
        "CONTROL_CENTER_POLICY": "best_redock",
        "CONTROL_CENTER_CLOSE_MAX_A": 8.0,
        "CTRL_REDOCK_EXHAUSTIVENESS": 64,
        "CTRL_REDOCK_NMODES": 9,

        # --- benchmark policy knobs ---
        "BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY": True,
        "BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING": False,
    }


# -------------------------
# Config loading & validation
# -------------------------
_ALLOWED_ENV_OVERRIDES = {
    "VINA_EXE","VINA_PATH","OPENBABEL_PATH","MGLTOOLS_PYTHON",
    "PREPARE_LIGAND_SCRIPT","PREPARE_RECEPTOR_SCRIPT","PYMOL_PATH",
    "P2RANK_PATH","PHENIX_DIR","PHENIX_LIB_PATH","PHENIX_CLEAN_SCRIPT",
    "INPUT_DIR","OUTPUT_DIR","PDBQT_DIR","DOCKED_DIR","LIGAND_DIR",
    "LIGAND_EXTRACTED_DIR","LIGANDS_MOL2_DIR","OUTPUT_LIGANDS_DIR",
    "P2RANK_OUTPUT_DIR",
    "CPU","CPU_ONLY","MAX_PARALLEL_JOBS","DOCKING_MODE","REDUCE_EXE",
    "USE_MEEKO","OVERALL_DIR","PROTEIN_DIR",

    # ---  allow ENV override for the knobs ---
    "BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY",
    "BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING",
    "CONTROL_CENTER_POLICY",
    "CONTROL_CENTER_CLOSE_MAX_A",
    "CTRL_REDOCK_EXHAUSTIVENESS",
    "CTRL_REDOCK_NMODES",
    # logging/topic gates (opt-in; safe to ignore if unset)
    "LOG_TOPICS", "LOG_LEVEL_FILE", "LOG_LEVEL_CONSOLE",

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
                return text[i:j+1]
        j += 1
    return None  # unbalanced
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
                out[k.strip().upper()] = v.strip()
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

def load_config(config_path: str = "config.txt", base_dir: Path | None = None) -> Dict[str, Any]:
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

    # overlay env vars (uppercased keys only)
    for k, v in os.environ.items():
        K = k.upper()
        if K in _ALLOWED_ENV_OVERRIDES:
            cfg[K] = v

    # type coercion (before expansion is fine)
    for k in ["CPU_ONLY","FORCE_REPROCESS","ALLOW_BOX_EXPAND","QUIET_CONSOLE","USE_MEEKO"]:
        if k in cfg:
            cfg[k] = _to_bool(cfg[k])

    for k in ["MAX_PARALLEL_JOBS","CPU","MAX_RECENTER_ATTEMPTS","EARLY_RECENTER_MIN_EVAL","CTRL_REDOCK_EXHAUSTIVENESS","CTRL_REDOCK_NMODES"]:
        if k in cfg:
            cfg[k] = _to_int(cfg[k], cfg[k])


    for k in ["EARLY_RECENTER_RATIO","EARLY_RECENTER_FAR_A","EARLY_RECENTER_MEDIAN_A","CONTROL_CENTER_CLOSE_MAX_A"]:
        if k in cfg:
            cfg[k] = _to_float(cfg[k], cfg[k])


    # normalize mode
    cfg["DOCKING_MODE"] = str(cfg.get("DOCKING_MODE","discovery")).lower()

    # now expand {OVERALL_DIR}/$OVERALL_DIR appearances
    cfg = _expand_all_vars(cfg)

    return cfg

def validate_config(cfg: Dict[str, Any]):
    """
    Require OVERALL_DIR and INPUT_DIR to exist.
    Auto-create typical output directories if missing.
    """
    required_keys = [
        "OUTPUT_DIR","INPUT_DIR","PDBQT_DIR","DOCKED_DIR","OVERALL_DIR",
        "MGLTOOLS_PYTHON","PREPARE_RECEPTOR_SCRIPT","VINA_EXE","MAX_PARALLEL_JOBS",
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
        "OUTPUT_DIR","PDBQT_DIR","DOCKED_DIR","OUTPUT_LIGANDS_DIR",
        "LIGAND_EXTRACTED_DIR","LIGANDS_MOL2_DIR","P2RANK_OUTPUT_DIR",
    ]
    for path_key in create_keys:
        p = Path(cfg[path_key])
        try:
            p.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise FileNotFoundError(f"Could not create directory for {path_key}: {p} ({e})")





import time, logging, tempfile

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

def emit_vina_config(
    cfg: Dict[str, Any],
    pdb_id: str,
    receptor_pdbqt: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    ligand_path: str,
    stage_name: str,
    stage_info: Dict[str, Any],
    cpu_per_job: int,
    logger: logging.Logger | None = None,
    *,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    legacy: bool = False,
):
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
    ph_label = (str(ph_token).strip() or None) if ph_token is not None else None
    legacy_mode = bool(legacy)

    from pathlib import Path

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

    manifest_path = router_config_file(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        name="vina.json",
        legacy=legacy_mode,
    )

    cfg_path = cfg_dir / f"{lig_base}_{stage_name}.txt"

    stage_root = router_docked_dir(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_dir = stage_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{lig_base}_{stage_name}.pdbqt"

    expected_receptor = router_receptor_file(
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

    if "seed" in stage_info:
        lines.append(f"seed = {int(stage_info['seed'])}")
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

    manifest_data: Dict[str, Any]
    entries_map: Dict[str, Dict[str, Any]]
    if manifest_path.exists():
        try:
            manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            manifest_data = {}
    else:
        manifest_data = {}

    entries = manifest_data.get("entries") if isinstance(manifest_data, dict) else None
    entries_map = {}
    if isinstance(entries, list):
        for item in entries:
            if isinstance(item, dict):
                lig = str(item.get("ligand", ""))
                if lig:
                    entries_map[lig] = item

    entry = {
        "ligand": lig_base,
        "config": str(cfg_path),
        "out": str(out_path),
        "receptor": receptor_for_config,
    }
    entries_map[lig_base] = entry

    manifest_data = {
        "run_id": run_id,
        "pdb_id": pdb_id,
        "stage": stage_name,
        "variant": variant_token,
        "ph": ph_label,
        "legacy": legacy_mode,
        "entries": [entries_map[k] for k in sorted(entries_map.keys())],
    }

    manifest_tmp = manifest_path.with_suffix(".part")
    try:
        manifest_tmp.parent.mkdir(parents=True, exist_ok=True)
        with open(manifest_tmp, "w", encoding="utf-8") as fh:
            json.dump(manifest_data, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(manifest_tmp, manifest_path)
    except FileNotFoundError:
        # Rare race when another process cleans up temp files; non-fatal for ctrl_redock.
        if logger:
            logger.warning(
                "[cfg.emit] manifest_tmp missing; skipping manifest update tmp=%s dest=%s stage=%s",
                str(manifest_tmp),
                str(manifest_path),
                stage_name,
            )
    except Exception:
        # Any manifest issue should not block docking; log and continue.
        if logger:
            logger.exception(
                "[cfg.emit] manifest update failed; continuing without manifest stage=%s path=%s",
                stage_name,
                str(manifest_path),
            )

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
            {"name": "stage1", "num_modes": 1,  "energy_range": 2, "exhaustiveness": 2},
            {"name": "stage2", "num_modes": 3,  "energy_range": 3, "exhaustiveness": 4},
            {"name": "stage3", "num_modes": 5,  "energy_range": 4, "exhaustiveness": 6},
            {"name": "stage4", "num_modes": 9,  "energy_range": 6, "exhaustiveness": 8},
            {"name": "stage5", "num_modes": 20, "energy_range": 9, "exhaustiveness": 20},
        ]
    elif mode == "polypharmacology":
        return [
            {"name": "stage1", "num_modes": 3,  "energy_range": 2, "exhaustiveness": 4},
            {"name": "stage2", "num_modes": 10, "energy_range": 6, "exhaustiveness": 12},
            {"name": "stage3", "num_modes": 20, "energy_range": 9, "exhaustiveness": 24},
        ]
    else:
        raise ValueError(f"Unknown docking mode: {mode}")


# -------------------------
# Score I/O
# -------------------------
def write_score_summary_to_csv(score_history, output_path="docking_score_summary.csv", run_id=None, variant=None):
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


# -------------------------
# Vina config writer, uses shim(lazy yeah)
# -------------------------
def generate_config(output_dir, pdb_id, receptor_pdbqt, center, box_size, ligand_path, stage, stage_info, cpu_per_job, docked_dir=None):
    # Legacy shim: derive a minimal cfg for older callers
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

def write_scores_csv(cfg, pdb_id, score_history):
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    variant_root = paths.docked_variant_root(var)
    os.makedirs(variant_root, exist_ok=True)
    csv_output_path = os.path.join(variant_root, "docking_score_summary.csv")

    run_id_value = str(cfg.get("RUN_ID") or "")

    flat_history = {}
    for stage_name, stage_map in score_history.items():
        flat_history[stage_name] = {}
        for lig, rec in stage_map.items():
            s = rec.get("score", None)
            if rec.get("valid", False):
                flat_history[stage_name][lig] = s if s is not None else ""
            else:
                flat_history[stage_name][lig] = (f"{s:.2f} (invalid)" if isinstance(s, (int, float)) else "(invalid)")

    write_score_summary_to_csv(
        flat_history,
        output_path=csv_output_path,
        run_id=run_id_value,
        variant=var,
    )
    return csv_output_path


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
    best = (
        df.groupby(lig_col, as_index=False)
          .agg(best_score=(score_col, "min"))
    )
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
        "pdb_path":             str(p.input_pdb_path),
        "nolig_pdb_path":       str(p.nolig_pdb_path),
        "ligand_output_dir":    str(p.ligand_output_dir),
        "ligands_mol2_dir":     str(p.ligands_mol2_dir),
        "prepped_ligands_dir":  str(p.prepped_ligands_dir),
        "cleaned_pdb_path":     str(p.receptor_cleaned_pdb(variant)),
        "receptor_pdbqt_path":  str(p.receptor_pdbqt(variant, ph_token=ph_token)),
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
