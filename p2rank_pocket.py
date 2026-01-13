import os
import csv
import subprocess
import logging

from installation import load_config

logger = logging.getLogger(__name__)

# Load config
config = load_config()
P2RANK_DIR = config.get("P2RANK_PATH")


def get_box_from_p2rank_csv(pdb_file):
    pdb_file = os.path.abspath(pdb_file)
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    pred_dir = os.path.join(P2RANK_DIR, "test_output", f"predict_{pdb_name}")
    pred_file = os.path.join(pred_dir, f"{pdb_name}.pdb_predictions.csv")

    logging.info(f"Running P2Rank for: {pdb_file}")

    # --- Choose P2Rank launcher deterministically (avoid PRANK-MSA on PATH) ---
    from shutil import which

    pr_base = None
    pr_exe = None
    pr_jar = None

    # Note: P2RANK_DIR is loaded from config (same value you call P2RANK_PATH in cfg)
    if P2RANK_DIR:
        cfg = P2RANK_DIR
        if os.path.isfile(cfg) and cfg.lower().endswith(".jar"):
            pr_jar = cfg
            pr_base = os.path.dirname(
                os.path.dirname(cfg)
            )  # .../bin/p2rank.jar -> install root
            exe_candidate = os.path.join(pr_base, "bin", "prank")
            if os.path.isfile(exe_candidate):
                pr_exe = exe_candidate
        elif os.path.isdir(cfg):
            pr_base = cfg
            exe_candidate = os.path.join(cfg, "bin", "prank")
            jar_candidate = os.path.join(cfg, "bin", "p2rank.jar")
            if os.path.isfile(exe_candidate):
                pr_exe = exe_candidate
            if os.path.isfile(jar_candidate):
                pr_jar = jar_candidate

    # Preference:
    # 1) Configured bin/prank (best: sets classpath)
    # 2) Configured jar (java -jar)
    # 3) PATH 'prank' *only if* it looks like P2Rank (not PRANK-MSA)
    run_cmd = None

    if pr_exe:
        run_cmd = [pr_exe, "predict"]
        logging.info("P2Rank launcher: bin/prank (configured)")
    elif pr_jar:
        run_cmd = ["java", "-Xmx4G", "-jar", pr_jar, "predict"]
        logging.info("P2Rank launcher: java -jar (configured)")
    else:
        pr_path = which("prank")
        if pr_path:
            # Heuristic guard: PRANK-MSA prints 'prank v.' and lacks 'predict' help
            try:
                out = subprocess.run(
                    [pr_path, "-version"], capture_output=True, text=True
                )
                banner = (out.stdout + out.stderr).lower()
                if "p2rank" in banner or "predict" in banner:
                    run_cmd = [pr_path, "predict"]
                    logging.info("P2Rank launcher: PATH prank (guarded OK)")
                else:
                    logging.error(
                        "Found '%s' on PATH, but it is PRANK (MSA), not P2Rank. "
                        "Set P2RANK_PATH to your P2Rank install dir.",
                        pr_path,
                    )
                    return None, None
            except Exception as e:
                logging.error(
                    "Unable to validate PATH prank (%s). Set P2RANK_PATH to the P2Rank install.",
                    e,
                )
                return None, None
        else:
            logging.error(
                "P2Rank not found. Set P2RANK_PATH to the install dir (with bin/prank) or bin/p2rank.jar."
            )
            return None, None

    # Always write to a known output folder next to the input PDB
    out_dir = os.path.join(os.path.dirname(os.path.abspath(pdb_file)), "_p2rank")
    os.makedirs(out_dir, exist_ok=True)

    try:
        # Explicit -o ensures we know exactly where predictions land
        cmd = run_cmd + ["-f", pdb_file, "-o", out_dir]
        logging.info("P2Rank cmd: %s", " ".join(cmd))
        subprocess.run(cmd, check=True, shell=False)
        logging.info(f"P2Rank ran successfully for {pdb_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"P2Rank failed: {e}")
        return None, None

    # Prefer the explicit output location; fall back to legacy locations if needed
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    pred_file = os.path.join(out_dir, f"{pdb_name}.pdb_predictions.csv")
    if not os.path.isfile(pred_file):
        # legacy fallback: current working dir default
        legacy = os.path.join(
            "test_output", f"predict_{pdb_name}", f"{pdb_name}.pdb_predictions.csv"
        )
        alt = os.path.join(
            pr_base or "",
            "test_output",
            f"predict_{pdb_name}",
            f"{pdb_name}.pdb_predictions.csv",
        )
        for probe in (legacy, alt):
            if probe and os.path.isfile(probe):
                pred_file = probe
                break

    if not os.path.isfile(pred_file):
        logging.warning(f"Prediction file not created: {pred_file}")
        return None, None
    with open(pred_file, "r", newline="") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [field.strip() for field in reader.fieldnames]

        for i, row in enumerate(reader):
            print(f"Row {i}: {row}")
            break  # for debugging

        # Rewind the file and parse again to get the top pocket
        f.seek(0)
        reader = csv.DictReader(f)
        reader.fieldnames = [field.strip() for field in reader.fieldnames]

        top_pocket = None
        for row in reader:
            try:
                rank_val = int(row.get("rank", "").strip())
                if rank_val == 1:
                    top_pocket = row
                    break
            except ValueError:
                logging.warning(f"Could not convert rank to int: {row.get('rank')}")

    if not top_pocket:
        logging.warning("No pocket found with rank 1")
        return None, None

    try:
        center = (
            float(top_pocket.get("center_x", "").strip()),
            float(top_pocket.get("center_y", "").strip()),
            float(top_pocket.get("center_z", "").strip()),
        )
        box_dim = float(top_pocket.get("surf_atoms", "20.0").strip())
        MAX_BOX_SIZE = 40.0
        box_size = (box_dim, box_dim, box_dim)
        clamped_box_size = tuple(min(dim, MAX_BOX_SIZE) for dim in box_size)
        print(f"[DEBUG] Returning center={center}, box_size={clamped_box_size}")
        return center, clamped_box_size
    except (KeyError, ValueError) as e:
        logging.error(f"Error parsing P2Rank pocket fields: {e}")
        return None, None


def detect_pocket(cleaned_pdb, logger):
    """
    Detect pocket center and box size from cleaned PDB.
    Returns (center, box_size) or (None, None) if failed.
    """
    center, box_size = get_box_from_p2rank_csv(cleaned_pdb)
    if center is None:
        logger.warning("Active-site detection failed.")
    return center, box_size
