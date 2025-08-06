import os
import csv
from pathlib import Path
import win32api  # Requires: pip install pywin32
import win32file

# === Tool Locator ===
def find_tool_on_any_drive(possible_subpaths):
    """Search through all fixed drives for given relative subpaths."""
    drives = [
        d for d in win32api.GetLogicalDriveStrings().split('\000')
        if d and win32file.GetDriveType(d) == win32file.DRIVE_FIXED
    ]
    for drive in drives:
        for subpath in possible_subpaths:
            candidate = Path(drive) / subpath
            if candidate.exists():
                return str(candidate.resolve())
    return None

# === Config Loader === currently overwrites the actual config AAAAAAA
def get_default_config():
    from pathlib import Path
    import os

    user_home = Path.home()
    base_dir = Path(__file__).resolve().parent

    def get_or_prompt(name, candidates):
        path = find_tool_on_any_drive(candidates)
        if path:
            return path
        manual = input(f"[!] {name} not found. Please paste full path or leave blank to skip: ").strip()
        return manual if manual else f"{name.lower()}_not_found"

    vina_path = get_or_prompt("AutoDock Vina", [
        "AutoDock-Vina-1.2.7/vina_1.2.7_win.exe",
        "Users/Public/AutoDock-Vina-1.2.7/vina_1.2.7_win.exe",
        "vina_1.2.7_win.exe"
    ])

    phenix_dir = get_or_prompt("Phenix bin", [
        "phenix-1.21.2-5419/Library/bin",
        "Program Files/Phenix/Library/bin",
        "phenix/Library/bin",
        "Library/bin",
        "Phenix/Library/bin",  # Capitalization may matter
        "E:/phenix/Library/bin"
    ])

    phenix_lib_path = get_or_prompt("Phenix site-packages", [
        "phenix-1.21.2-5419/Lib/site-packages",
        "phenix/Lib/site-packages",
        "Program Files/Phenix/Lib/site-packages",
        "Lib/site-packages"
    ])

    p2rank_path = get_or_prompt("P2Rank", [
        "p2rank_2.5",
        "Users/Public/p2rank_2.5",
        "Program Files/p2rank_2.5"
    ])

    mgltools_path = get_or_prompt("MGLTools", [
        "MGLTools-1.5.7",
        "Program Files (x86)/MGLTools-1.5.7",
        "Programs/MGLTools-1.5.7"
    ])

    openbabel_path = get_or_prompt("OpenBabel", [
        "OpenBabel-3.1.1/obabel.exe",
        "Program Files (x86)/OpenBabel-3.1.1/obabel.exe",
        "Program Files/OpenBabel-3.1.1/obabel.exe"
    ])

    reduce_exe = get_or_prompt("Reduce", [
        "miniconda3/envs/docking-env/Library/bin/reduce.exe",
        "Library/bin/reduce.exe",
        "Program Files/reduce.exe",
        "reduce.exe",
        "phenix/Library/bin/reduce.exe",
        "E:/phenix/Library/bin/reduce.exe",
    ])
    pymol_exe = get_or_prompt("PyMOL Python", [
        "E:/pymol/python.exe",
        "E:/Program Files/PyMOL/python.exe",
        "C:/Program Files/PyMOL/python.exe",
        "python.exe",
    ])

    mgltools_python = str(Path(mgltools_path) / "python.exe") if "not_found" not in mgltools_path else "MGLTOOLS_PYTHON_NOT_FOUND"
    prepare_receptor_script = str(Path(mgltools_path) / "Lib/site-packages/AutoDockTools/Utilities24/prepare_receptor4.py") if "not_found" not in mgltools_path else "PREPARE_RECEPTOR_SCRIPT_NOT_FOUND"

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
        "MGLTOOLS_PATH": mgltools_path,
        "MGLTOOLS_PYTHON": mgltools_python,
        "PREPARE_RECEPTOR_SCRIPT": prepare_receptor_script,
        "MGLTOOLS_DIR": mgltools_path,
        "OPENBABEL_PATH": openbabel_path,
        "REDUCE_EXE": reduce_exe,
        "CPU_ONLY": True,
        "CPU": 12,
        "MAX_PARALLEL_JOBS": 6,
        "FORCE_REPROCESS": False,
        "PYMOL_PATH": pymol_exe,
        "DOCKING_MODE": "discovery",
    }

    # Ensure all keys are uppercase
    return {k.upper(): v for k, v in config.items()}


def load_inputs():
    config_path = Path("config.txt")
    if config_path.exists():
        config = {}
        with open(config_path) as f:
            for line in f:
                if "=" in line:
                    k, v = line.strip().split("=", 1)
                    config[k.strip().upper()] = v.strip()
        # Fill in missing keys from default
        default_cfg = get_default_config()
        config["MAX_PARALLEL_JOBS"] = int(config["MAX_PARALLEL_JOBS"])
        merged_cfg = {**default_cfg, **config}
        config["DOCKING_MODE"] = config.get("DOCKING_MODE", "discovery").lower()
        if "MAX_PARALLEL_JOBS" in merged_cfg:
            merged_cfg["MAX_PARALLEL_JOBS"] = int(merged_cfg["MAX_PARALLEL_JOBS"])
        return merged_cfg
    else:
        return get_default_config()


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
            {"name": "stage1", "num_modes": 5, "energy_range": 4, "exhaustiveness": 8},
            {"name": "stage2", "num_modes": 10, "energy_range": 6, "exhaustiveness": 12},
            {"name": "stage3", "num_modes": 20, "energy_range": 9, "exhaustiveness": 24},
        ]
    else:
        raise ValueError(f"Unknown docking mode: {mode}")

# === Docking Score Output ===

def write_score_summary_to_csv(score_history, output_path="docking_score_summary.csv"):
    all_ligands = set()
    for stage_scores in score_history.values():
        all_ligands.update(os.path.basename(lig) for lig in stage_scores)

    all_ligands = sorted(all_ligands)
    stages = sorted(score_history.keys())
    header = ["Ligand"] + stages

    rows = []
    for ligand in all_ligands:
        row = [ligand]
        for stage in stages:
            stage_scores = score_history.get(stage, {})
            matched = [score for path, score in stage_scores.items() if os.path.basename(path) == ligand]

            if matched:
                value = matched[0]
                if isinstance(value, (float, int)):
                    row.append(f"{value:.2f}")
                else:
                    row.append(str(value))  # fall back to string if it's not numeric
            else:
                row.append("")
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
