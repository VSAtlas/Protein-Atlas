# run_vina.py
import os
import re
import subprocess
import logging
import json
from uuid import uuid4
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Any, Iterable

from docking.docking_utils import write_failure_marker
from docking.score_io import (
    VINA_MODEL_START_RE,
    extract_best_vina_score,
    extract_best_vina_score_from_lines,
)
from config.runtime_config import load_config, validate_config
from path_router.path_router import (
    make_paths,
    config_dir as router_config_dir,
    docked_dir as router_docked_dir,
)

_CFG = load_config()
validate_config(_CFG)

_log = logging.getLogger("vina")

_FS_ERRORS: tuple[type[BaseException], ...] = (OSError,)
"""Expected failures writing markers, unlinking temps, reading manifests."""


def _write_failure_marker_safe(
    target: Path,
    reason: str,
    stdout_tail: Optional[str] = None,
    stderr_tail: Optional[str] = None,
) -> None:
    try:
        write_failure_marker(target, reason, stdout_tail, stderr_tail)
    except _FS_ERRORS as exc:
        _log.warning(
            "[vina.emit] failure_marker_unwritable path=%s docking_reason=%s err=%s",
            target,
            reason,
            exc,
        )


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
        if isinstance(val, (list, tuple)):
            try:
                return [float(x) for x in val]
            except (TypeError, ValueError):
                return default_list
        if isinstance(val, str) and val.strip():
            try:
                return [float(x) for x in val.split(",")]
            except ValueError:
                return default_list
        return default_list

    disc_default = [1.0, 0.1, 0.01, 0.001, 0.001]
    poly_default = [1.0, 0.05, 0.005]
    disc = _parse_pcts(_CFG.get("DISCOVERY_SELECTION_PCTS", disc_default), disc_default)
    poly = _parse_pcts(_CFG.get("POLYPHARM_SELECTION_PCTS", poly_default), poly_default)
    percentages = {"discovery": disc, "polypharmacology": poly}.get(
        docking_mode, [1.0] * len(stages)
    )
    pct = percentages[i + 1] if i + 1 < len(percentages) else 0.0
    num_to_select = int(len(scores) * pct)
    if num_to_select < 1:
        logger.warning(
            f"Percentage {pct * 100:.5f}% yielded <1 ligand. Using best-scoring ligand."
        )
        num_to_select = 1

    top_ligands = sorted(
        (
            (lig_name, score)
            for lig_name, score in scores.items()
            if isinstance(score, (int, float))
        ),
        key=lambda x: x[1],
    )[:num_to_select]
    next_list = [lig_name for lig_name, _score in top_ligands if lig_name]
    logger.info(
        f"Selected top {len(next_list)} ligands ({pct * 100:.5f}%) for next stage."
    )
    return next_list


# ----------------------------
# Vina stdout/stderr filtering
# ----------------------------
_ERR_PAT = re.compile(
    r"(error|failed|fatal|exception|segmentation|not found|cannot|invalid|timeout)",
    re.IGNORECASE,
)
_CPU_CFG_RE = re.compile(r"^\s*cpu\s*=.*$", re.IGNORECASE | re.MULTILINE)


def _should_filter_stdout() -> bool:
    v = os.environ.get("FILTER_VINA_STDOUT")
    if v is not None:
        return v.strip().lower() in {"1", "true", "yes", "on"}
    return str(_CFG.get("FILTER_VINA_STDOUT", "false")).strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


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
    if isinstance(v, str) and v:
        try:
            t = int(v)
            return t if t > 0 else None
        except ValueError:
            return None
    try:
        timeout_value = _CFG.get("VINA_TIMEOUT_SEC", 0)
        t = int(timeout_value if timeout_value not in (None, "") else 0)
        return t if t > 0 else None
    except (TypeError, ValueError):
        return None


def _summarize_vina_failure(stderr: str, returncode: int) -> str:
    lines = [ln.strip() for ln in str(stderr or "").splitlines() if ln.strip()]
    if not lines:
        return f"Return code {returncode}"
    for ln in lines:
        low = ln.lower()
        if "cannot be specified more than once" in low and "--cpu" in low:
            return "Command line parse error: option '--cpu' cannot be specified more than once"
        if "command line parse error" in low:
            return ln
    return lines[-1]


# ----------------------------
# Run a single docking task
# ----------------------------
def run_docking_task(
    vina_exe: str,
    config_path: str,
    ligand_name: str,
    out_path: str,
    *,
    write_failure_marker_flag: bool = False,
    cpu_override: int | None = None,
):
    """
    Runs Vina with a prepared config file and extracts the best score.
    Returns (ligand_name, best_score or None).
    """
    temp_cfg_path: Optional[Path] = None
    timeout_sec: Optional[int] = None
    try:
        from shutil import which

        out_path_obj = Path(out_path)
        if not Path(config_path).exists():
            raise FileNotFoundError(f"[!] Vina config not found at: {config_path}")

        # Resolve vina_exe: absolute path OR on PATH
        resolved_exe = vina_exe
        if not Path(resolved_exe).exists():
            found = which(vina_exe)
            if not found:
                raise FileNotFoundError(
                    f"[!] AutoDock Vina not found (VINA_EXE='{vina_exe}', not a file and not on PATH)"
                )
            resolved_exe = found

        Path(out_path).parent.mkdir(parents=True, exist_ok=True)

        timeout_sec = _get_timeout()

        config_for_run = str(config_path)
        append_cli_cpu = False
        cpu_n: Optional[int] = None
        if cpu_override is not None:
            try:
                cpu_n = max(1, int(cpu_override))
            except (TypeError, ValueError):
                cpu_n = None
        if cpu_n is not None:
            cfg_text = Path(config_path).read_text(encoding="utf-8", errors="ignore")
            if _CPU_CFG_RE.search(cfg_text):
                # Keep exactly one cpu setting in rewritten configs to avoid parse ambiguity.
                without_cpu = _CPU_CFG_RE.sub("", cfg_text)
                overridden = f"{without_cpu.rstrip()}\ncpu = {cpu_n}\n"
                if overridden != cfg_text:
                    temp_cfg_path = Path(config_path).with_name(
                        f"{Path(config_path).stem}.cpu{cpu_n}.{os.getpid()}.{uuid4().hex}.tmp"
                    )
                    temp_cfg_path.write_text(overridden, encoding="utf-8")
                    config_for_run = str(temp_cfg_path)
            else:
                append_cli_cpu = True

        cmd = [resolved_exe, "--config", config_for_run]
        if append_cli_cpu and cpu_n is not None:
            cmd.extend(["--cpu", str(cpu_n)])
        cmd.extend(["--out", out_path])

        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            timeout=timeout_sec,
            check=False,
        )

        stdout_tail = proc.stdout or ""
        stderr_tail = proc.stderr or ""
        failure_reason = None

        if _should_filter_stdout():
            _maybe_print_useful_lines(proc.stdout, proc.stderr)
        elif proc.returncode != 0 and proc.stderr:
            _log.warning("[vina.emit] %s", proc.stderr.strip())

        if proc.returncode != 0:
            msg = _summarize_vina_failure(proc.stderr or "", proc.returncode)
            _log.warning("[vina.emit] Docking failed for %s: %s", ligand_name, msg)
            failure_reason = msg

        if not out_path_obj.exists():
            msg = failure_reason or "missing_output"
            _log.warning(
                "[vina.emit] Docking produced no output for %s (config=%s)",
                ligand_name,
                config_path,
            )
            if write_failure_marker_flag:
                _write_failure_marker_safe(
                    out_path_obj, msg, stdout_tail, stderr_tail
                )
            return ligand_name, None

        score, n_models = extract_best_vina_score(out_path)
        if score is None:
            _log.info(
                "[parser] %s no Vina score found in %s (models_in_file=%d)",
                ligand_name,
                out_path,
                n_models,
            )
            failure_reason = failure_reason or "no_score"

        if write_failure_marker_flag and failure_reason:
            _write_failure_marker_safe(
                out_path_obj,
                failure_reason,
                stdout_tail,
                stderr_tail,
            )

        return ligand_name, score

    except subprocess.TimeoutExpired:
        to = timeout_sec if timeout_sec is not None else _get_timeout()
        _log.error(
            "[vina.emit] Docking timed out for %s (limit=%ss)",
            ligand_name,
            to if to is not None else "unset",
        )
        if write_failure_marker_flag:
            _write_failure_marker_safe(Path(out_path), "timeout")
        return ligand_name, None
    except FileNotFoundError as e:
        _log.error("[vina.emit] File error for %s: %s", ligand_name, e)
        if write_failure_marker_flag:
            _write_failure_marker_safe(Path(out_path), f"file_error:{e}")
        return ligand_name, None
    except Exception as e:
        _log.exception("[vina.emit] Docking crashed for %s", ligand_name)
        if write_failure_marker_flag:
            _write_failure_marker_safe(Path(out_path), f"exception:{e}")
        return ligand_name, None
    finally:
        if temp_cfg_path is not None:
            try:
                temp_cfg_path.unlink(missing_ok=True)
            except _FS_ERRORS as exc:
                _log.debug(
                    "[vina.emit] temp_cfg_unlink_failed path=%s err=%s",
                    temp_cfg_path,
                    exc,
                )


# ----------------------------
# Pose helpers
# ----------------------------
def extract_models(pdbqt_path: str) -> List[List[str]]:
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    models: List[List[str]] = []
    current: List[str] = []
    for line in lines:
        if VINA_MODEL_START_RE.match(line):
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
    pdbqt_path_obj = Path(pdbqt_path)
    receptor_pdbqt_obj = Path(receptor_pdbqt)
    models = extract_models(str(pdbqt_path_obj))
    best_valid_score: Optional[float] = None
    best_valid_model: Optional[Path] = None
    temp_paths: List[Path] = []
    clash_thr = float(_CFG.get("POSE_CLASH_THRESHOLD_A", 2.0))
    clash_tol = int(_CFG.get("POSE_CLASH_TOLERANCE", 3))
    dmax_surf = float(_CFG.get("POSE_DIST_SURFACE_MAX_A", 6.0))
    dmax_cent = float(_CFG.get("POSE_DIST_CENTROID_MAX_A", 4.5))

    try:
        for i, model_lines in enumerate(models):
            temp_path = Path(str(pdbqt_path_obj).replace(".pdbqt", f"_model{i + 1}.pdbqt"))
            with open(temp_path, "w", encoding="utf-8") as f:
                f.writelines(model_lines)
            temp_paths.append(temp_path)

            result = validate_fn(
                protein_pdbqt=str(receptor_pdbqt_obj),
                ligand_pdbqt=str(temp_path),
                pocket_center=center,
                clash_threshold=clash_thr,
                CLASH_TOLERANCE=clash_tol,
                DIST_THRESHOLD_SURFACE=dmax_surf,
                DIST_THRESHOLD_CENTROID=dmax_cent,
                surface_atom_coords=surface_coords,
            )

            if result.get("valid", False):
                score = extract_best_vina_score_from_lines(model_lines)
                if best_valid_score is None or (
                    score is not None and score < best_valid_score
                ):
                    best_valid_score = score
                    best_valid_model = temp_path

        return best_valid_model, best_valid_score

    finally:
        # Clean up intermediates
        for p in temp_paths:
            try:
                if p.exists():
                    p.unlink()
            except _FS_ERRORS as exc:
                _log.debug(
                    "[vina.pose] temp_model_unlink_failed path=%s err=%s",
                    p,
                    exc,
                )


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

    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    stage_root = router_docked_dir(
        pdb_id, variant=variant, ph_tag=ph_token, legacy=legacy
    )
    stage_dir = stage_root / stage

    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage,
        variant=variant,
        ph_tag=ph_token,
        legacy=legacy,
    )

    return stage_dir, cfg_dir


MANIFEST_NAME = "vina.json"


def _load_manifest(manifest_path: Path) -> Dict[str, Any]:
    try:
        raw = manifest_path.read_text(encoding="utf-8")
    except _FS_ERRORS as exc:
        _log.debug(
            "[vina.manifest] read_failed path=%s err=%s", manifest_path, exc
        )
        return {}
    except UnicodeDecodeError as exc:
        _log.warning(
            "[vina.manifest] decode_failed path=%s err=%s",
            manifest_path,
            exc,
        )
        return {}
    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        _log.warning(
            "[vina.manifest] json_invalid path=%s err=%s", manifest_path, exc
        )
        return {}
    if isinstance(data, dict):
        return data
    _log.warning(
        "[vina.manifest] unexpected_shape path=%s type=%s",
        manifest_path,
        type(data).__name__,
    )
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
    cfg_path = str(job["config"])
    ligand = str(job.get("ligand") or Path(cfg_path).stem)
    out_path_value = job.get("out")
    out_path = str(out_path_value) if out_path_value else ""
    _log.info(
        "[run.dispatch] pdb=%s stage=%s variant=%s ph=%s cfg=%s",
        job.get("pdb_id"),
        job.get("stage"),
        variant_display,
        ph_display,
        cfg_path,
    )
    return run_docking_task(vina_exe, cfg_path, ligand, out_path)
