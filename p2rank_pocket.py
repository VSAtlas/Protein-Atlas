import os
import csv
import subprocess
import logging
import shlex
from shutil import which

from installation import load_config

logger = logging.getLogger(__name__)

# Load config
config = load_config()
P2RANK_DIR = config.get("P2RANK_PATH")
P2RANK_JAVA_EXE = str(config.get("P2RANK_JAVA_EXE", "java") or "java")
P2RANK_JAVA_OPTS = str(config.get("P2RANK_JAVA_OPTS", "-Xmx4G") or "-Xmx4G")


def _java_base_cmd():
    cmd = [P2RANK_JAVA_EXE]
    opts = P2RANK_JAVA_OPTS.strip()
    if opts:
        cmd.extend(shlex.split(opts))
    return cmd


def _resolve_p2rank_install(p2rank_path):
    """
    Resolve configured P2Rank install artifacts.
    Returns (prank_exe, p2rank_jar, install_base_dir, lib_glob_or_none).
    """
    if not p2rank_path:
        return None, None, None, None

    raw = os.path.expanduser(str(p2rank_path).strip())
    if not raw:
        return None, None, None, None

    if os.path.isfile(raw):
        lower = raw.lower()
        if lower.endswith(".jar"):
            jar_path = raw
            bin_dir = os.path.dirname(jar_path)
            install_base = os.path.dirname(bin_dir)
            lib_glob = os.path.join(bin_dir, "lib", "*")
            prank_exe = os.path.join(bin_dir, "prank")
            if not os.path.isfile(prank_exe):
                prank_exe = None
            return prank_exe, jar_path, install_base, lib_glob
        if os.path.basename(lower) in {"prank", "prank.bat", "prank.cmd"}:
            return raw, None, os.path.dirname(os.path.dirname(raw)), None
        return None, None, None, None

    if os.path.isdir(raw):
        if os.path.basename(raw).lower() == "bin":
            bin_dir = raw
            install_base = os.path.dirname(raw)
        else:
            install_base = raw
            bin_dir = os.path.join(raw, "bin")

        prank_exe = os.path.join(bin_dir, "prank")
        if not os.path.isfile(prank_exe):
            prank_exe = None
        jar_path = os.path.join(bin_dir, "p2rank.jar")
        if not os.path.isfile(jar_path):
            jar_path = None
        lib_glob = os.path.join(bin_dir, "lib", "*")
        return prank_exe, jar_path, install_base, lib_glob

    return None, None, None, None


def _launch_attempts_from_config():
    prank_exe, jar_path, install_base, lib_glob = _resolve_p2rank_install(P2RANK_DIR)
    attempts = []

    if prank_exe:
        attempts.append(("configured_prank", [prank_exe, "predict"], install_base))

    if jar_path:
        java_cmd = _java_base_cmd()
        if os.path.isdir(os.path.dirname(lib_glob or "")):
            classpath = os.pathsep.join([jar_path, lib_glob])
            attempts.append(
                (
                    "configured_classpath",
                    java_cmd
                    + ["-cp", classpath, "cz.siret.prank.program.Main", "predict"],
                    install_base,
                )
            )
        attempts.append(
            ("configured_jar", java_cmd + ["-jar", jar_path, "predict"], install_base)
        )

    return attempts


def _path_prank_attempt():
    pr_path = which("prank")
    if not pr_path:
        return None

    try:
        out = subprocess.run([pr_path, "-version"], capture_output=True, text=True)
        banner = (out.stdout + out.stderr).lower()
        if "p2rank" in banner or "predict" in banner:
            return ("path_prank", [pr_path, "predict"], None)
    except Exception:
        return None

    logging.error(
        "Found '%s' on PATH, but it looks like PRANK (MSA), not P2Rank. "
        "Set P2RANK_PATH to your P2Rank install dir.",
        pr_path,
    )
    return None


def get_box_from_p2rank_csv(pdb_file):
    pdb_file = os.path.abspath(pdb_file)
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

    logging.info(f"Running P2Rank for: {pdb_file}")

    launch_attempts = _launch_attempts_from_config()
    if not launch_attempts:
        path_attempt = _path_prank_attempt()
        if path_attempt:
            launch_attempts = [path_attempt]

    if not launch_attempts:
        logging.error(
            "P2Rank not found. Set P2RANK_PATH to install root/bin dir/bin/prank/bin/p2rank.jar."
        )
        return None, None

    # Always write to a known output folder next to the input PDB
    out_dir = os.path.join(os.path.dirname(os.path.abspath(pdb_file)), "_p2rank")
    os.makedirs(out_dir, exist_ok=True)

    pr_base = None
    failure_summaries = []
    ran_ok = False
    for mode, run_cmd, attempt_base in launch_attempts:
        cmd = run_cmd + ["-f", pdb_file, "-o", out_dir]
        logging.info("P2Rank launcher: %s", mode)
        logging.info("P2Rank cmd: %s", " ".join(cmd))
        pr_base = attempt_base or pr_base

        try:
            completed = subprocess.run(cmd, check=False, shell=False, capture_output=True, text=True)
        except Exception as e:
            failure_summaries.append(f"{mode}: exec_error={e}")
            continue

        if completed.returncode == 0:
            logging.info(f"P2Rank ran successfully for {pdb_file}")
            ran_ok = True
            break

        err_blob = "\n".join(
            x for x in [(completed.stderr or "").strip(), (completed.stdout or "").strip()] if x
        )
        if "NoClassDefFoundError: groovy/lang/GroovyObject" in err_blob:
            logging.warning(
                "P2Rank launcher '%s' missing Groovy runtime on classpath; trying next launcher.",
                mode,
            )
        failure_summaries.append(f"{mode}: returncode={completed.returncode}")

    if not ran_ok:
        logging.error("P2Rank failed after %d launcher attempt(s): %s", len(launch_attempts), "; ".join(failure_summaries))
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
