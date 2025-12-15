# run_vina.py
import os
import re
import subprocess
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Any, Iterable
import logging, json
_log = logging.getLogger("vina")

from input_and_export_functions import extract_best_score
try:
    from input_and_export_functions import load_config, validate_config
    _CFG = load_config("config.txt")
    validate_config(_CFG)
except Exception:
    _CFG = {}





import io

from path_router import (
    make_paths,
    config_dir as router_config_dir,
    docked_dir as router_docked_dir,
)

_SCORE_LINE = re.compile(r"REMARK\s+VINA\s+RESULT[:\s]+(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_MODEL_START = re.compile(r"^\s*MODEL\b", re.IGNORECASE)

def _extract_best_score_robust(pdbqt_path: str) -> tuple[float|None, int]:
    """
    Robust score scan that tolerates spacing/case/line-endings and multi-MODEL files.
    Returns (best_energy, n_models_seen).
    """
    p = Path(pdbqt_path)
    if not p.exists() or p.stat().st_size < 32:
        return None, 0
    best = None
    n_models = 0
    try:
        with open(p, "r", encoding="utf-8", errors="ignore") as fh:
            buf = []
            saw_model_tag = False
            for ln in fh:
                if _MODEL_START.match(ln):
                    # new model begins
                    if buf:
                        n_models += 1
                    buf = [ln]
                    saw_model_tag = True
                else:
                    buf.append(ln)
                # Examine every line for a score marker
                m = _SCORE_LINE.search(ln)
                if m:
                    try:
                        e = float(m.group(1))
                        if (best is None) or (e < best):
                            best = e
                    except Exception:
                        pass
            # count last unclosed model if we saw any MODEL tag
            if saw_model_tag:
                n_models += 1
    except Exception:
        return None, n_models
    # If no explicit MODEL tags, treat file as single model (n_models=1)
    if n_models == 0:
        n_models = 1
    return best, n_models

    
    
# ----------------------------
# Selection helper 
# ----------------------------
def select_for_next_stage(docking_mode, i, stages, scores, logger):
    """
    Select a fraction of ligands to carry forward to the next stage based on mode and rank.
    Returns the list of next-stage ligands.
    """
    if not scores:
        return []

    def _parse_pcts(val, default_list):
        if isinstance(val, (list, tuple)): return [float(x) for x in val]
        if isinstance(val, str) and val.strip():
            try:
                return [float(x) for x in val.split(",")]
            except Exception:
                return default_list
        return default_list

    disc_default = [1.0, 0.1, 0.01, 0.001, 0.001]
    poly_default = [1.0, 0.05, 0.005]
    disc = _parse_pcts(_CFG.get("DISCOVERY_SELECTION_PCTS", disc_default), disc_default)
    poly = _parse_pcts(_CFG.get("POLYPHARM_SELECTION_PCTS", poly_default), poly_default)
    percentages = {"discovery": disc, "polypharmacology": poly}.get(docking_mode, [1.0] * len(stages))
    pct = percentages[i + 1] if i + 1 < len(percentages) else 0.0
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
    v = os.environ.get("FILTER_VINA_STDOUT")
    if v is not None:
        return v.strip().lower() in {"1","true","yes","on"}
    return str(_CFG.get("FILTER_VINA_STDOUT", "false")).strip().lower() in {"1","true","yes","on"}


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
                _log.warning("[vina.emit] %s", ln)

def _get_timeout() -> Optional[int]:
    v = os.environ.get("VINA_TIMEOUT_SEC")
    if v not in (None, ""):
        try:
            t = int(v);  return t if t > 0 else None
        except Exception:
            return None
    try:
        t = int(_CFG.get("VINA_TIMEOUT_SEC", 0))
        return t if t > 0 else None
    except Exception:
        return None


# ----------------------------
# Run a single docking task
# ----------------------------
def run_docking_task(vina_exe: str, config_path: str, ligand_name: str, out_path: str):
    """
    Runs Vina with a prepared config file and extracts the best score.
    Returns (ligand_name, best_score or None).
    """
    try:
        from shutil import which
        from pathlib import Path

        if not Path(config_path).exists():
            raise FileNotFoundError(f"[!] Vina config not found at: {config_path}")

        # Resolve vina_exe: absolute path OR on PATH
        resolved_exe = vina_exe
        if not Path(resolved_exe).exists():
            found = which(vina_exe)
            if not found:
                raise FileNotFoundError(f"[!] AutoDock Vina not found (VINA_EXE='{vina_exe}', not a file and not on PATH)")
            resolved_exe = found

        Path(out_path).parent.mkdir(parents=True, exist_ok=True)

        timeout = _get_timeout()

        proc = subprocess.run(
            [resolved_exe, "--config", config_path, "--out", out_path],
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )


        if _should_filter_stdout():
            _maybe_print_useful_lines(proc.stdout, proc.stderr)
        elif proc.returncode != 0 and proc.stderr:
            _log.warning("[vina.emit] %s", proc.stderr.strip())

        if proc.returncode != 0:
            msg = proc.stderr.strip().splitlines()[-1] if proc.stderr else f"Return code {proc.returncode}"
            _log.warning("[vina.emit] Docking failed for %s: %s", ligand_name, msg)


        # Primary parse using project helper
        score = extract_best_score(out_path)

        # Hardened fallback: tolerant scan across all MODEL blocks
        if score is None:
            robust_score, n_models = _extract_best_score_robust(out_path)
            score = robust_score
            if score is None:
                _log.info("[parser] %s no Vina score found in %s (models_in_file=%d)", ligand_name, out_path, n_models)

        return ligand_name, score
    
    except subprocess.TimeoutExpired:
        _log.error("[vina.emit] Docking timed out for %s (> %ss)", ligand_name, _get_timeout())
        return ligand_name, None
    except FileNotFoundError as e:
        _log.error("[vina.emit] File error for %s: %s", ligand_name, e)
        return ligand_name, None
    except Exception as e:
        _log.error("[vina.emit] Docking crashed for %s: %s", ligand_name, e)
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
    clash_thr = float(_CFG.get("POSE_CLASH_THRESHOLD_A", 2.0))
    clash_tol = int(_CFG.get("POSE_CLASH_TOLERANCE", 3))
    dmax_surf = float(_CFG.get("POSE_DIST_SURFACE_MAX_A", 6.0))
    dmax_cent = float(_CFG.get("POSE_DIST_CENTROID_MAX_A", 4.5))

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
                clash_threshold=clash_thr,
                CLASH_TOLERANCE=clash_tol,
                DIST_THRESHOLD_SURFACE=dmax_surf,
                DIST_THRESHOLD_CENTROID=dmax_cent,
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


# ----------------------------
# Path routing helpers
# ----------------------------
def resolve_stage_and_config_dirs(
    cfg,
    pdb_id: str,
    run_id: str,
    stage: str,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    legacy: bool = False,
) -> Tuple[Path, Path]:
    """Return (stage_dir, config_dir) using the centralized path router."""

    if not pdb_id:
        raise ValueError("pdb_id must be provided")

    # >>> PATHS INIT START
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # >>> PATHS INIT END

    # >>> DOCKED PATHS PATCH START
    stage_root = router_docked_dir(pdb_id, variant=variant, ph_tag=ph_token, legacy=legacy)
    stage_dir = stage_root / stage
    # >>> DOCKED PATHS PATCH END

    # >>> CONFIG PATHS PATCH START
    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage,
        variant=variant,
        ph_tag=ph_token,
        legacy=legacy,
    )
    # >>> CONFIG PATHS PATCH END

    return stage_dir, cfg_dir


MANIFEST_NAME = "vina.json"


def _load_manifest(manifest_path: Path) -> Dict[str, Any]:
    try:
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        if isinstance(data, dict):
            return data
    except Exception:
        pass
    return {}


def iter_stage_jobs(
    cfg: Dict[str, Any],
    run_id: str,
    pdb_id: str,
    stage: str,
    variants: Iterable[Optional[str]],
    ph_tags: Optional[Iterable[Optional[str]]],
    legacy: bool,
) -> List[Dict[str, Any]]:
    jobs: List[Dict[str, Any]] = []
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_list = list(variants) or [None]
    ph_iterable = list(ph_tags) if ph_tags is not None else [None]

    for variant in variant_list:
        for ph_token in ph_iterable:
            cfg_dir = router_config_dir(
                run_id,
                pdb_id,
                stage,
                variant=variant,
                ph_tag=ph_token,
                legacy=legacy,
            )
            manifest_path = cfg_dir / MANIFEST_NAME
            if not manifest_path.exists():
                continue
            payload = _load_manifest(manifest_path)
            entries = payload.get("entries") if isinstance(payload, dict) else None
            if not isinstance(entries, list):
                continue
            for item in entries:
                if not isinstance(item, dict):
                    continue
                config_path = item.get("config") or item.get("config_path")
                out_path = item.get("out") or item.get("out_path")
                ligand = item.get("ligand")
                receptor = item.get("receptor")
                if not (config_path and out_path and ligand):
                    continue
                job = {
                    "pdb_id": pdb_id,
                    "stage": stage,
                    "variant": variant,
                    "ph": ph_token,
                    "config": config_path,
                    "out": out_path,
                    "ligand": ligand,
                    "receptor": receptor,
                }
                jobs.append(job)
    return jobs


def dispatch_job(vina_exe: str, job: Dict[str, Any]):
    variant_display = job.get("variant") or "None"
    ph_display = job.get("ph") or "None"
    cfg_path = job["config"]
    ligand = job.get("ligand") or Path(cfg_path).stem
    out_path = job.get("out")
    _log.info(
        "[run.dispatch] pdb=%s stage=%s variant=%s ph=%s cfg=%s",
        job.get("pdb_id"),
        job.get("stage"),
        variant_display,
        ph_display,
        cfg_path,
    )
    return run_docking_task(vina_exe, cfg_path, ligand, out_path)
