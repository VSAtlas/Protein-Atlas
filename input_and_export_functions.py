import os, sys, shutil
from pathlib import Path
from distutils.util import strtobool
from typing import Any, Dict
import re, csv, math
from collections import defaultdict
# >>> PATHS IMPORT START
from path_router import make_paths
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
    Expand {OVERALL_DIR} and $OVERALL_DIR occurrences inside config values.
    We intentionally keep it simple, not a full env expansion.
    """
    if not isinstance(val, str):
        return val
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
):
    # Determine active variant from the per-pass env that main.py sets
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_label = (cfg.get("_ACTIVE_PH_LABEL") or "").strip() or None

    # ---  define lig_base from ligand_path (fixes NameError) ---
    from pathlib import Path
    lig_base = Path(ligand_path).stem

    # Variant-aware config dir (prefer path_router helpers; fallback to old layout)
    run_id = cfg["RUN_ID"]
    if hasattr(paths, "configs_stage_dir"):
        conf_dir = paths.configs_stage_dir(run_id, var, stage_name, ph_label)
    else:
        # Fallback: configs/<RUN_ID>/<PDB>/<VARIANT>/<stage> (omit VARIANT if None)
        conf_dir = Path(cfg["CONFIG_RUN_DIR"]) / pdb_id
        conf_dir = (conf_dir / var / stage_name) if var else (conf_dir / stage_name)
    conf_dir.mkdir(parents=True, exist_ok=True)

    # Variant-aware Vina output dir (prefer helper; fallback to old layout)
    if hasattr(paths, "docked_stage_dir"):
        out_dir = paths.docked_stage_dir(var, stage_name, ph_label)
    else:
        # Fallback: docked/<PDB>/<VARIANT>/<stage> (omit VARIANT if None)
        root = paths.docked_pdb_root()  # expected to be a Path-like
        out_dir = (root / var / stage_name) if var else (root / stage_name)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / f"{lig_base}_{stage_name}.pdbqt"

    if logger:
        logger.debug(
            "[emit_vina_config] variant=%r stage=%s lig_base=%s conf_dir=%s out_path=%s",
            var, stage_name, lig_base, str(conf_dir), str(out_path)
        )


    lines = [
        f"receptor = {receptor_pdbqt}",
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
    # --- AUDIT: compact Vina config trace ---
    if logger:
        cx, cy, cz = center
        sx, sy, sz = box_size
        lig_name = os.path.basename(str(ligand_path))
        logger.info("[vina.cfg] lig=%s center=(%.3f,%.3f,%.3f) size=(%.1f,%.1f,%.1f)",
                    lig_name, cx, cy, cz, sx, sy, sz)
    cfg_path = conf_dir / f"{lig_base}_{stage_name}.txt"
    payload = ("\n".join(lines)).encode("utf-8")
    overwrite = cfg_path.exists()


    # atomic write
    tmp = cfg_path.with_suffix(".part")
    with open(tmp, "wb") as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, cfg_path)

    msg = f"[cfg.emit] pdb={pdb_id} stage={stage_name} ligand={lig_base} path={cfg_path} overwrite={str(overwrite).lower()} bytes={len(payload)}"
    (logger.info(msg) if logger else print(msg))

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
    return emit_vina_config(cfg, pdb_id, receptor_pdbqt, center, box_size, ligand_path, stage, stage_info, cpu_per_job, logger=None)





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
