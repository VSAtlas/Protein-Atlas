from pathlib import Path
import subprocess
import sys
import ctypes
import os
from ctypes import wintypes, create_unicode_buffer
from pathlib import Path
import subprocess
import logging
from activesite import fix_pdb_elements

os.environ["BABEL_DATADIR"] = "E:/OpenBabel-3.1.1/data"
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "SEC", "PYL", "MSE"
}

def prep_ligands_from_pdb(ligand_output_dir: Path, ligands_mol2_dir: Path, prepped_ligands_dir: Path):
    logging.info("Starting ligand preparation from PDB files")

    config = read_config()
    mgltools_python = config.get("MGLTOOLS_PYTHON")
    mgltools_path = config.get("MGLTOOLS_PATH")
    obabel_exe = config.get("OPENBABEL_PATH")

    if not mgltools_python or not mgltools_path or not obabel_exe:
        raise RuntimeError("Missing paths in config.txt: MGLTOOLS_PYTHON, MGLTOOLS_PATH, OPENBABEL_PATH")

    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = Path(mgltools_path) / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))
    obabel_exe_short = get_short_path_name(obabel_exe)

    pdb_files = list(ligand_output_dir.glob("*.pdb"))
    logging.info(f"Found {len(pdb_files)} PDB ligand file(s)")

    for pdb_file in pdb_files:
        logging.info(f"Processing: {pdb_file.name}")

        # skip standard amino acids
        residue_name = pdb_file.stem.split("_")[0].upper()
        if residue_name in STANDARD_AMINO_ACIDS:
            logging.info(f"Skipping standard amino acid residue: {pdb_file.name}")
            continue

        # Optional: Skip ligands with very few atoms
        with open(pdb_file, "r") as f:
            atom_lines = [line for line in f if line.startswith("HETATM") or line.startswith("ATOM")]
        if len(atom_lines) < 5:
            logging.info(f"Skipping tiny ligand (fewer than 5 atoms): {pdb_file.name}")
            continue

        # Fix malformed element columns before Open Babel reads it
        fix_pdb_elements(str(pdb_file))

        mol2_file = ligands_mol2_dir / f"{pdb_file.stem}.mol2"
        mol2_file.parent.mkdir(parents=True, exist_ok=True)
        if mol2_file.suffix.lower() != ".mol2":
            logging.warning(f"Unexpected file extension for {mol2_file.name}, skipping.")
            continue

        obabel_cmd = [
            obabel_exe_short, "-ipdb", str(pdb_file),
            "--addh", "--gen3d", "--partialcharge", "gasteiger",
            "-omol2", "-O", str(mol2_file)
        ]

        try:
            subprocess.run(obabel_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Open Babel failed for {pdb_file.name}:\n{e}")
            continue

        # Now that Open Babel has run, check file validity
        if not mol2_file.exists() or mol2_file.stat().st_size < 100:
            logging.warning(f"Mol2 file is empty or malformed, skipping: {mol2_file.name}")
            continue
        from rdkit import Chem

        #  check mol2 file using RDKit
        rdkit_mol = Chem.MolFromMol2File(str(mol2_file), sanitize=False)
        if rdkit_mol is None:
            logging.warning(f"RDKit could not parse MOL2, skipping: {mol2_file.name}")
            continue
        try:
            Chem.SanitizeMol(rdkit_mol)
            logging.info(f"RDKit sanitization succeeded: {mol2_file.name}")
        except Exception as e:
            logging.warning(f"RDKit sanitization failed for {mol2_file.name}: {e}")
            logging.warning(f"Proceeding with unsanitized ligand for docking (may be invalid)")
            # Optional: move to "warned" folder instead of quarantining
            warned_dir = ligands_mol2_dir / "warned"
            warned_dir.mkdir(exist_ok=True)
            (mol2_file).replace(warned_dir / mol2_file.name)
            mol2_file = warned_dir / mol2_file.name

        pdbqt_path = prepped_ligands_dir / f"{pdb_file.stem}.pdbqt"
        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
        prepare_cmd = [
            mgltools_python_short,
            prepare_script_short,
            "-l", get_short_path_name(str(mol2_file.resolve())),
            "-o", get_short_path_name(str(pdbqt_path.resolve())),
            "-A", "checkhydrogens"
        ]

        logging.info(f"Preparing ligand: {' '.join(prepare_cmd)}")
        try:
            result = subprocess.run(
                prepare_cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=str(mol2_file.parent)
            )
            logging.info(result.stdout)
            logging.info(f"Created: {pdbqt_path.name}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to prepare {mol2_file.name}:\n{e.stderr}")

def read_config(path="config.txt"):
    config = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                config[key.strip()] = value.strip()
    return config

def get_short_path_name(long_name):
    if sys.platform != 'win32':
        return long_name

    _GetShortPathNameW = ctypes.windll.kernel32.GetShortPathNameW
    _GetShortPathNameW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
    _GetShortPathNameW.restype = wintypes.DWORD

    output_buf_size = 260
    output_buf = create_unicode_buffer(output_buf_size)
    ret = _GetShortPathNameW(long_name, output_buf, output_buf_size)

    if ret == 0 or ret > output_buf_size:
        #print(f"⚠️ Warning: Failed to get short path for {long_name}")
        return long_name
    return output_buf.value

def convert_sdf_to_mol2_single(sdf_path: Path, mol2_output_dir: Path, obabel_exe: str) -> Path:
    mol2_output_dir.mkdir(parents=True, exist_ok=True)
    output_file = mol2_output_dir / f"{sdf_path.stem}.mol2"

    command = [
        obabel_exe,
        "-isdf", str(sdf_path),
        "--gen3d",
        "-omol2",
        "-O", str(output_file)
    ]

    print("🔄 Running Open Babel:", " ".join(map(str, command)))
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(f"Open Babel stderr:\n{result.stderr}")
        return output_file
    except subprocess.CalledProcessError as e:
        print(f" Open Babel conversion error:\n{e.stderr}")
        return None


# Im splitting the ligand files up into individual files because vina only
#accepts individual ligands for the docking function
def split_multi_mol2(multi_mol2_path: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(multi_mol2_path, "r") as f:
        content = f.read()

    molecules = content.split("@<TRIPOS>MOLECULE")[1:]  # skip leading blank chunk
    paths = []

    for i, mol_text in enumerate(molecules, 1):
        mol_data = "@<TRIPOS>MOLECULE" + mol_text
        out_file = output_dir / f"{multi_mol2_path.stem}_{i:04d}.mol2"
        with open(out_file, "w") as out_f:
            out_f.write(mol_data)
        paths.append(out_file)

    print(f"Split {len(paths)} molecules from {multi_mol2_path.name}")
    return paths


#literally the same thing as prep ligands from pdbs but for sdf files
def prep_ligands_with_mgltools():
    print("Starting ligand preparation")

    config = read_config()

    ligand_extracted_dir = Path(config.get("ligand_extracted_dir", "")).resolve()
    ligands_mol2_dir = Path(config.get("ligands_mol2_dir", "")).resolve()
    output_ligands_dir = Path(config.get("output_ligands_dir", "")).resolve()
    output_ligands_dir.mkdir(parents=True, exist_ok=True)

    mgltools_python = config.get("mgltools_python")
    mgltools_path = config.get("mgltools_path")
    obabel_exe = config.get("openbabel_path")

    if not mgltools_python or not mgltools_path or not obabel_exe:
        raise RuntimeError("Missing paths in config.txt: mgltools_python, mgltools_path, openbabel_path")

    # Short paths required for MGLTools compatibility
    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = Path(mgltools_path) / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    if not prepare_script.exists():
        raise FileNotFoundError(f"prepare_ligand4.py not found at {prepare_script}")
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))
    obabel_exe_short = get_short_path_name(obabel_exe)

    sdf_files = list(ligand_extracted_dir.glob("*.sdf"))
    print(f"🔍 Found {len(sdf_files)} SDF file(s)")

    if not sdf_files:
        return

    for sdf_file in sdf_files:
        print(f"\nProcessing: {sdf_file.name}")
        sdf_abs = sdf_file.resolve()

        # Step 1: Convert to multi-mol2
        mol2_multi_file = convert_sdf_to_mol2_single(sdf_abs, ligands_mol2_dir, obabel_exe_short)
        if not mol2_multi_file:
            continue

        # Step 2: Split multi-mol2
        mol2_files = split_multi_mol2(mol2_multi_file, ligands_mol2_dir)
        if not mol2_files:
            continue

        # Step 3: Prepare with MGLTools
        for mol2_file in mol2_files:
            pdbqt_path = output_ligands_dir / f"{mol2_file.stem}.pdbqt"
            mol2_short = mol2_file.name  # just the filename since cwd will be set
            pdbqt_short = get_short_path_name(str(pdbqt_path.resolve()))

            command = [
                mgltools_python_short,
                prepare_script_short,
                "-l", mol2_short,
                "-o", pdbqt_short,
                "-A", "hydrogens"
            ]

            print("🧪 Preparing ligand:", " ".join(command))
            try:
                result = subprocess.run(
                    command,
                    check=True,
                    capture_output=True,
                    text=True,
                    cwd=str(mol2_file.parent)  # ensure prepare_ligand4.py runs in correct folder
                )
                print(result.stdout)
                print(f"Created: {pdbqt_path.name}")
            except subprocess.CalledProcessError as e:
                print(f"Failed to prepare {mol2_file.name}:\n{e.stderr}")

from pathlib import Path

def is_valid_ligand(path: Path, log_dir: Path) -> bool:
    """
    Minimal check that ensures the .pdbqt is not malformed.
    Returns True if basic structure seems intact.
    """
    try:
        with open(path, 'r') as f:
            lines = f.readlines()

        atom_lines = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        torsion_lines = [line for line in lines if line.startswith("TORSDOF")]

        if not atom_lines:
            raise ValueError("No ATOM or HETATM lines found.")
        if not torsion_lines:
            raise ValueError("No torsion info (TORSDOF) found.")

        return True

    except Exception as e:
        malformed_log = Path(log_dir) / "malformed_ligands.txt"
        with open(malformed_log, "a") as f:
            f.write(f"{path.name} - PDBQT validation failed: {e}\n")
        return False


if __name__ == "__main__":
    prep_ligands_with_mgltools()
