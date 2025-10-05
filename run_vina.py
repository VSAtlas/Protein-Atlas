# run_vina.py
import os
import re
import subprocess
from pathlib import Path
from typing import Tuple, Optional, List

from input_and_export_functions import extract_best_score

# ----------------------------
# Selection helper (kept here to match your current imports)
# ----------------------------
def select_for_next_stage(docking_mode, i, stages, scores, logger):
    """
    Select a fraction of ligands to carry forward to the next stage based on mode and rank.
    Returns the list of next-stage ligands.
    """
    if not scores:
        return []

    percentages = {
        "discovery": [1.0, 0.1, 0.01, 0.001, 0.001],
        "polypharmacology": [1.0, 0.05, 0.005],
    }.get(docking_mode, [1.0] * len(stages))

    pct = percentages[i + 1] if i + 1 < len(percentages) else 0.01
    num_to_select = int(len(scores) * pct)
    if num_to_select < 1:
        logger.warning(f"Percentage {pct*100:.5f}% yielded <1 ligand. Using best-scoring ligand.")
        num_to_select = 1

    top_ligands = sorted(
        ((l, s) for l, s in scores.items() if isinstance(s, (int, float))),
        key=lambda x: x[1]
    )[:num_to_select]
    next_list = [l for l, _ in top_ligands if l]
    logger.info(f"Selected top {len(next_list)} ligands ({pct * 100:.5f}%) for next stage.")
    return next_list


# ----------------------------
# Vina stdout/stderr filtering
# ----------------------------
_ERR_PAT = re.compile(
    r"(error|failed|fatal|exception|segmentation|not found|cannot|invalid|timeout)",
    re.IGNORECASE,
)

def _should_filter_stdout() -> bool:
    v = os.environ.get("FILTER_VINA_STDOUT", "").strip().lower()
    return v in {"1", "true", "yes", "on"}

def _maybe_print_useful_lines(stdout: str, stderr: str) -> None:
    """
    Print only lines that look like actual problems. This hides banners/progress noise.
    """
    if not stdout and not stderr:
        return
    for stream, label in ((stdout, "stdout"), (stderr, "stderr")):
        if not stream:
            continue
        for ln in stream.splitlines():
            if _ERR_PAT.search(ln):
                # Surface only meaningful lines
                print(f"[vina:{label}] {ln}")

def _get_timeout() -> Optional[int]:
    try:
        t = int(os.environ.get("VINA_TIMEOUT_SEC", "0"))
        return t if t > 0 else None
    except Exception:
        return None


# ----------------------------
# Run a single docking task
# ----------------------------
def run_docking_task(vina_exe: str, config_path: str, ligand_name: str, out_path: str) -> Tuple[str, Optional[float]]:
    """
    Runs Vina with a prepared config file.
    - Captures stdout/stderr.
    - Optionally filters out banner/progress noise (FILTER_VINA_STDOUT=1).
    - Optional timeout via VINA_TIMEOUT_SEC env var.
    - Returns (ligand_name, best_score or None).
    """
    try:
        if not Path(vina_exe).exists():
            raise FileNotFoundError(f"[!] AutoDock Vina not found at: {vina_exe}")
        if not Path(config_path).exists():
            raise FileNotFoundError(f"[!] Vina config not found at: {config_path}")

        # Run vina
        timeout = _get_timeout()
        proc = subprocess.run(
            [vina_exe, "--config", config_path],
            text=True,
            capture_output=True,   # capture so we can filter noise
            timeout=timeout,
            check=False,           # handle return codes ourselves to still extract scores when possible
        )

        # Optionally print only interesting lines
        if _should_filter_stdout():
            _maybe_print_useful_lines(proc.stdout, proc.stderr)
        else:
            # Fall back to original behavior: show full stderr only if nonzero exit
            if proc.returncode != 0 and proc.stderr:
                print(proc.stderr)

        # If Vina failed, still try to parse score (file might exist/contain partial output)
        if proc.returncode != 0:
            # Surface a compact failure line
            msg = proc.stderr.strip().splitlines()[-1] if proc.stderr else f"Return code {proc.returncode}"
            print(f"Docking failed for {ligand_name}: {msg}")

        # Extract best score (if pose file exists)
        score = extract_best_score(out_path)
        return ligand_name, score

    except subprocess.TimeoutExpired:
        print(f"Docking timed out for {ligand_name} (>{_get_timeout()}s)")
        return ligand_name, None
    except FileNotFoundError as e:
        print(f"File error for {ligand_name}: {e}")
        return ligand_name, None
    except Exception as e:
        print(f"Docking crashed for {ligand_name}: {e}")
        return ligand_name, None


# ----------------------------
# Pose helpers
# ----------------------------
def extract_models(pdbqt_path: str) -> List[List[str]]:
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    models: List[List[str]] = []
    current: List[str] = []
    for line in lines:
        if line.startswith("MODEL"):
            current = [line]
        elif line.startswith("ENDMDL"):
            current.append(line)
            models.append(current[:])
        else:
            current.append(line)
    return models

def validate_all_poses(
    pdbqt_path: str,
    receptor_pdbqt: str,
    center,
    surface_coords,
    validate_fn,
):
    pdbqt_path = Path(pdbqt_path)
    receptor_pdbqt = Path(receptor_pdbqt)
    models = extract_models(str(pdbqt_path))
    best_valid_score: Optional[float] = None
    best_valid_model: Optional[Path] = None
    temp_paths: List[Path] = []

    try:
        for i, model_lines in enumerate(models):
            temp_path = Path(str(pdbqt_path).replace(".pdbqt", f"_model{i + 1}.pdbqt"))
            with open(temp_path, "w", encoding="utf-8") as f:
                f.writelines(model_lines)
            temp_paths.append(temp_path)

            result = validate_fn(
                protein_pdbqt=str(receptor_pdbqt),
                ligand_pdbqt=str(temp_path),
                pocket_center=center,
                clash_threshold=2.0,
                CLASH_TOLERANCE=3,
                DIST_THRESHOLD_SURFACE=6.0,
                DIST_THRESHOLD_CENTROID=4.5,
                surface_atom_coords=surface_coords,
            )

            if result.get("valid", False):
                score = extract_best_score(str(temp_path))
                if best_valid_score is None or (score is not None and score < best_valid_score):
                    best_valid_score = score
                    best_valid_model = temp_path

        return best_valid_model, best_valid_score

    finally:
        # Clean up intermediates
        for p in temp_paths:
            try:
                if p.exists():
                    p.unlink()
            except Exception:
                pass
