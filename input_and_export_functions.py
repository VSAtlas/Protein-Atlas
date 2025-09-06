import os, sys, logging, csv
from pathlib import Path
from collections import defaultdict
from distutils.util import strtobool

IS_WINDOWS = sys.platform.startswith("win")

# Windows-only modules (pywin32). Guard them so Linux/WSL can import this file.
try:
    if IS_WINDOWS:
        import win32api      # type: ignore
        import win32file     # type: ignore
    else:
        win32api = None
        win32file = None
except Exception:
    win32api = None
    win32file = None


def list_local_drives():
    """Return a list of drive/mount roots on this OS."""
    if IS_WINDOWS and win32api:
        return [d for d in win32api.GetLogicalDriveStrings().split("\000") if d]
    # Linux/macOS: collect mount roots (WSL exposes Windows drives under /mnt/*)
    drives = []
    for base in ("/mnt", "/media"):
        if os.path.isdir(base):
            for name in os.listdir(base):
                p = os.path.join(base, name)
                if os.path.ismount(p):
                    # mimic 'C:\\' style ending with separator for callers that expect it
                    drives.append(p if p.endswith(os.sep) else p + os.sep)
    # always include home as a last-resort root
    drives.append(str(Path.home()))
    return drives

def get_short_path(path: str) -> str:
    """Windows short path if available; otherwise the original path."""
    if IS_WINDOWS and win32api:
        try:
            return win32api.GetShortPathName(path)
        except Exception:
            return path
    return path

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

# === Tool Locator ===
def find_tool_on_any_drive(possible_subpaths):
    """Search through all fixed drives for given relative subpaths."""
    drives = list_local_drives()
    for drive in drives:
        for subpath in possible_subpaths:
            candidate = Path(drive) / subpath
            if candidate.exists():
                return str(candidate.resolve())
    return None

# === Config Loader (non-interactive) ===
def get_default_config(prompt=True):
    from pathlib import Path
    import os, shutil, sys

    NONINTERACTIVE = (not prompt) or (not sys.stdin.isatty())

    def which_or_exists(candidates):
        """Return first existing absolute path or which() result for a name."""
        for c in candidates:
            p = Path(c)
            if p.is_absolute() and p.exists():
                return str(p)
            w = shutil.which(Path(c).name)
            if w:
                return w
        return None

    def get_or_prompt(name, candidates):
        found = find_tool_on_any_drive(candidates) or which_or_exists(candidates)
        if found:
            return found
        if NONINTERACTIVE:
            return f"{name.upper().replace(' ','_')}_NOT_FOUND"
        try:
            manual = input(f"[!] {name} not found. Please paste full path or leave blank to skip: ").strip()
        except EOFError:
            manual = ""
        return manual if manual else f"{name.upper().replace(' ','_')}_NOT_FOUND"

    base_dir = Path(__file__).resolve().parent

    vina_path = get_or_prompt("AutoDock Vina", [
        "/home/michael/miniconda3/envs/docking-env/bin/vina", "vina"
    ])

    phenix_dir = get_or_prompt("Phenix bin", [
        "/mnt/e/phenix/Library/bin"
    ])
    phenix_lib_path = get_or_prompt("Phenix site-packages", [
        "/mnt/e/phenix/Lib/site-packages"
    ])

    p2rank_path = get_or_prompt("P2Rank", ["/mnt/e/p2rank_2.5", "p2rank_2.5"])

    # Prefer CONDA MGL/ADT scripts if present
    mgltools_python = get_or_prompt("MGLTools Python", [
        "/home/michael/miniconda3/envs/docking-env/bin/python"
    ])
    prepare_ligand_script = get_or_prompt("prepare_ligand4.py", [
        "/home/michael/miniconda3/envs/docking-env/bin/prepare_ligand4.py", "prepare_ligand4.py"
    ])
    prepare_receptor_script = get_or_prompt("prepare_receptor4.py", [
        "/home/michael/miniconda3/envs/docking-env/bin/prepare_receptor4.py", "prepare_receptor4.py"
    ])

    openbabel_path = get_or_prompt("OpenBabel", [
        "/home/michael/miniconda3/envs/docking-env/bin/obabel", "obabel"
    ])

    reduce_exe = get_or_prompt("Reduce", [
        "/mnt/e/phenix/Library/bin/reduce.exe", "reduce"
    ])

    pymol_path = get_or_prompt("PyMOL", [
        "/home/michael/miniconda3/envs/docking-env/bin/pymol", "pymol"
    ])

    config = {
        "OVERALL_DIR": str(base_dir),
        "INPUT_DIR": str(base_dir / "input_pdbs"),
        "PROTEIN_DIR": str(base_dir / "pdbqts"),
        "LIGAND_DIR": str(base_dir / "input_ligands"),
        "LIGAND_EXTRACTED_DIR": str(base_dir / "extracted_ligands"),
        "LIGANDS_MOL2_DIR": str(base_dir / "ligands_mol2"),
        "OUTPUT_LIGANDS_DIR": str(base_dir / "prepped_ligands"),
        "OUTPUT_DIR": str(base_dir / "processed_pdbs"),
        "PDBQT_DIR": str(base_dir / "pdbqts"),
        "DOCKED_DIR": str(base_dir / "docked"),

        "PHENIX_CLEAN_SCRIPT": str(base_dir / "phenix_clean.py"),

        "VINA_PATH": vina_path,
        "VINA_EXE": vina_path,

        "PHENIX_DIR": phenix_dir,
        "PHENIX_LIB_PATH": phenix_lib_path,

        "P2RANK_PATH": p2rank_path,

        "MGLTOOLS_PATH": "/home/michael/miniconda3/envs/docking-env",  # logical home
        "MGLTOOLS_PYTHON": mgltools_python,
        "MGLTOOLS_DIR": mgltools_python,
        "PREPARE_LIGAND_SCRIPT": prepare_ligand_script,
        "PREPARE_RECEPTOR_SCRIPT": prepare_receptor_script,

        "OPENBABEL_PATH": openbabel_path,
        "REDUCE_EXE": reduce_exe,

        "CPU_ONLY": True,
        "CPU": 12,
        "MAX_PARALLEL_JOBS": 6,
        "FORCE_REPROCESS": False,
        "PYMOL_PATH": pymol_path,
        "DOCKING_MODE": "discovery",
    }
    return {k.upper(): v for k, v in config.items()}


def load_inputs():
    import os
    config_path = Path("config.txt")
    parsed = {}

    if config_path.exists():
        with open(config_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "=" in line:
                    k, v = line.split("=", 1)
                    parsed[k.strip().upper()] = v.strip()

    # Normalize common aliases from config
    if "PYMOL_EXE" in parsed and "PYMOL_PATH" not in parsed:
        parsed["PYMOL_PATH"] = parsed["PYMOL_EXE"]

    # Overlay environment variables (from .env.wsl)
    for k, v in os.environ.items():
        K = k.upper()
        if K in {
            "VINA_EXE","VINA_PATH","OPENBABEL_PATH","MGLTOOLS_PYTHON",
            "PREPARE_LIGAND_SCRIPT","PREPARE_RECEPTOR_SCRIPT","PYMOL_PATH",
            "P2RANK_PATH","PHENIX_DIR","PHENIX_LIB_PATH","PHENIX_CLEAN_SCRIPT",
            "INPUT_DIR","OUTPUT_DIR","PDBQT_DIR","DOCKED_DIR","LIGAND_DIR",
            "LIGAND_EXTRACTED_DIR","LIGANDS_MOL2_DIR","OUTPUT_LIGANDS_DIR",
            "CPU","CPU_ONLY","MAX_PARALLEL_JOBS","DOCKING_MODE","REDUCE_EXE","USE_MEEKO",  
        }:
            parsed[K] = v

    default_cfg = get_default_config(prompt=False)
    merged_cfg = {**default_cfg, **parsed}  # user/env overrides defaults

    # Coerce types
    for k in ["CPU_ONLY","FORCE_REPROCESS","ALLOW_BOX_EXPAND",
              "QUIET_CONSOLE","RECEPTOR_SANITY_CHECK","CHECKPOINT_ENABLE",
              "FILTER_VINA_STDOUT","USE_MEEKO"]:
        if k in merged_cfg:
            merged_cfg[k] = _to_bool(merged_cfg[k])

    for k in ["MAX_PARALLEL_JOBS","CPU","MAX_RECENTER_ATTEMPTS","EARLY_RECENTER_MIN_EVAL"]:
        if k in merged_cfg:
            merged_cfg[k] = _to_int(merged_cfg[k], merged_cfg[k])

    for k in ["EARLY_RECENTER_RATIO","EARLY_RECENTER_FAR_A","EARLY_RECENTER_MEDIAN_A"]:
        if k in merged_cfg:
            merged_cfg[k] = _to_float(merged_cfg[k], merged_cfg[k])

    merged_cfg["DOCKING_MODE"] = str(merged_cfg.get("DOCKING_MODE","discovery")).lower()
    return merged_cfg


# === Config Validation ===

def validate_config(cfg):
    required_keys = [
        "OUTPUT_DIR", "INPUT_DIR", "MGLTOOLS_DIR", "MGLTOOLS_PYTHON",
        "PREPARE_RECEPTOR_SCRIPT", "VINA_EXE", "MAX_PARALLEL_JOBS",
        "PDBQT_DIR", "DOCKED_DIR", "OVERALL_DIR", "DOCKING_MODE","PYMOL_PATH",
    ]

    missing = [key for key in required_keys if key not in cfg]
    if missing:
        raise ValueError(f"Missing required config keys: {missing}")

    for path_key in ["INPUT_DIR", "OUTPUT_DIR", "MGLTOOLS_DIR", "PDBQT_DIR", "DOCKED_DIR", "OVERALL_DIR"]:
        path = Path(cfg[path_key])
        if not path.exists():
            raise FileNotFoundError(f"Required path does not exist: {path}")

# === Docking Stage Definitions ===

def define_docking_stages(mode="discovery"):
    if mode == "discovery":
        return [
            {"name": "stage1", "num_modes": 1, "energy_range": 2, "exhaustiveness": 2},
            {"name": "stage2", "num_modes": 3, "energy_range": 3, "exhaustiveness": 4},
            {"name": "stage3", "num_modes": 5, "energy_range": 4, "exhaustiveness": 6},
            {"name": "stage4", "num_modes": 9, "energy_range": 6, "exhaustiveness": 8},
            {"name": "stage5", "num_modes": 20, "energy_range": 9, "exhaustiveness": 20},
        ]
    elif mode == "polypharmacology":
        return [
            {"name": "stage1", "num_modes": 3, "energy_range": 2, "exhaustiveness": 4},
            {"name": "stage2", "num_modes": 10, "energy_range": 6, "exhaustiveness": 12},
            {"name": "stage3", "num_modes": 20, "energy_range": 9, "exhaustiveness": 24},
        ]
    else:
        raise ValueError(f"Unknown docking mode: {mode}")

# === Docking Score Output ===
import re

def write_score_summary_to_csv(score_history, output_path="docking_score_summary.csv"):
    # Normalize ligand names (remove _stageX)
    ligand_stage_pattern = re.compile(r"^(.*?)(_stage\d+)?\.pdbqt$", re.IGNORECASE)

    ligand_scores = defaultdict(dict)

    for stage, stage_scores in score_history.items():
        for path, score in stage_scores.items():
            ligand_filename = os.path.basename(path)
            match = ligand_stage_pattern.match(ligand_filename)
            if match:
                ligand_core = f"{match.group(1)}.pdbqt"
                ligand_scores[ligand_core][stage] = score

    # Get sorted list of unique ligands and stages
    all_ligands = sorted(ligand_scores.keys())
    all_stages = sorted(score_history.keys())

    header = ["Ligand"] + all_stages
    rows = []

    for ligand in all_ligands:
        row = [ligand]
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


# === Extract Best Docking Score ===

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

# === Docking Config Writer ===
def generate_config(output_dir, pdb_id, receptor_pdbqt, center, box_size, ligand_path, stage, stage_info, cpu_per_job):
    config_dir = os.path.join(output_dir, "configs", pdb_id, stage)
    os.makedirs(config_dir, exist_ok=True)

    ligand_name = os.path.splitext(os.path.basename(ligand_path))[0]
    out_path = os.path.join(output_dir, "docked", pdb_id, stage, f"{ligand_name}_{stage}.pdbqt")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    config_lines = [
        f"receptor = {receptor_pdbqt}",
        f"ligand = {ligand_path}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x = {box_size[0]:.3f}",
        f"size_y = {box_size[1]:.3f}",
        f"size_z = {box_size[2]:.3f}",
        f"num_modes = {stage_info['num_modes']}",
        f"energy_range = {stage_info['energy_range']}",
        f"exhaustiveness = {stage_info['exhaustiveness']}",
        f"out = {out_path}",
        f"cpu = {cpu_per_job}"
    ]

    config_path = os.path.join(config_dir, f"{ligand_name}.txt")
    with open(config_path, "w") as f:
        f.write("\n".join(config_lines))
    return config_path, out_path

import re
import math

from typing import Any
def record_score(score_history, stage_name, ligand, score: Any, valid, reason=None):
    def coerce_score(x):
        # numeric already?
        if isinstance(x, (int, float)):
            return float(x)
        # PDBQT path? try to read best score
        if isinstance(x, str) and x.lower().endswith(".pdbqt"):
            try:
                return extract_best_score(x)
            except Exception:
                return None
        # numeric string?
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

    # Keep native Vina ordering (more negative = better)
    # If you want invalid to always be worse than same-number valid, add a small bump:
    return s if rec.get("valid", False) else s + 1e-6

def write_scores_csv(cfg, pdb_id, score_history):
    """
    Flatten score_history and write docking_score_summary.csv into DOCKED_DIR/<pdb>/.
    """
    protein_dock_dir = os.path.join(cfg["DOCKED_DIR"], pdb_id)
    os.makedirs(protein_dock_dir, exist_ok=True)
    csv_output_path = os.path.join(protein_dock_dir, "docking_score_summary.csv")

    flat_history = {}
    for stage_name, stage_map in score_history.items():
        flat_history[stage_name] = {}
        for lig, rec in stage_map.items():
            s = rec.get("score", None)
            if rec.get("valid", False):
                flat_history[stage_name][lig] = s if s is not None else ""
            else:
                flat_history[stage_name][lig] = (f"{s:.2f} (invalid)" if isinstance(s, (int, float)) else "(invalid)")

    write_score_summary_to_csv(flat_history, output_path=csv_output_path)
    return csv_output_path

def build_paths_for_protein(cfg, base_id, pdb_file):
    """
    Build common paths/dirs for a protein and ensure required folders exist.
    Returns a dict with input/output paths.
    """
    pdb_id = base_id
    ligand_output_dir = Path(cfg["OUTPUT_DIR"]) / pdb_id / f"{pdb_id}_cleaned_ligands"
    ligands_mol2_dir = Path(cfg["LIGANDS_MOL2_DIR"]) / pdb_id
    prepped_ligands_dir = Path(cfg["OUTPUT_LIGANDS_DIR"]) / pdb_id
    protein_dir = Path(cfg["OUTPUT_DIR"]) / f"{base_id}_nolig"
    cleaned_pdb_path = protein_dir / f"{base_id}_nolig_cleaned.pdb"
    receptor_pdbqt_path = Path(cfg["PDBQT_DIR"]) / f"{base_id}_receptor.pdbqt"
    pdb_path = os.path.join(cfg["INPUT_DIR"], pdb_file)
    nolig_pdb_path = os.path.join(cfg["OUTPUT_DIR"], f"{base_id}_nolig.pdb")

    ligand_output_dir.mkdir(parents=True, exist_ok=True)
    prepped_ligands_dir.mkdir(parents=True, exist_ok=True)
    protein_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)

    return {
        "pdb_id": pdb_id,
        "pdb_path": pdb_path,
        "nolig_pdb_path": nolig_pdb_path,
        "ligand_output_dir": ligand_output_dir,
        "ligands_mol2_dir": ligands_mol2_dir,
        "prepped_ligands_dir": prepped_ligands_dir,
        "cleaned_pdb_path": cleaned_pdb_path,
        "receptor_pdbqt_path": receptor_pdbqt_path,
    }

