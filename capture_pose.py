from __future__ import annotations

"""
PyMOL helpers for native and docked-pose screenshots + simple selection utilities.

Features
- pick_control_and_nearest_rdk(): choose a control and a "nearest-in-score" RDK ligand
- _render_three_views_with_pymol(): save front/side/top PNGs for receptor + ligands
- _render_native_on_original_pdb(): show native ligands from the original PDB
- write_multiview_pml(): create a .pml with stored scenes (front/side/top)
- render_pml_headless(): optional CLI fallback to render PNGs without pymol2
- _safe_open_csv_for_write(): Windows-friendly writer when CSV is locked
"""

from string import Template
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import subprocess
import time
import os, sys, json, statistics, tempfile, threading
# ------------------------
# Exclusions: common ions/solvents/cofactors to hide in "native" views
# ------------------------
EXCLUDE_HET_IDS = {
    "HOH", "WAT", "NA", "K", "CL", "MG", "MN", "CA", "ZN", "FE", "CO", "CU", "NI", "MO",
    "SO4", "PO4", "ACT", "ACE", "IPH", "FMT", "BME", "MPD", "DMS", "IPA", "IMD", "DTT",
    "TRS", "MES", "HEP", "CIT", "TAR", "TLA", "GLY", "EDO", "GOL", "PEG",
    "GLC", "GAL", "MAN", "NAG", "BMA", "FUC", "TRE", "BGC", "BOG",
    "HEM", "FAD", "FMN", "NAD", "NAP", "NADH", "SAM", "SAH",
}


# >>> PATHS IMPORT START
from path_router import make_paths
# >>> PATHS IMPORT END
from concurrent.futures import ProcessPoolExecutor, as_completed
try:
    from input_and_export_functions import load_config, validate_config
    _cfg = load_config("config.txt") or {}
    try:
        validate_config(_cfg)
    except Exception:
        pass
except Exception:
    _cfg = {}

def _cfg_bool(name: str, *, env: str | None = None, default: bool = False) -> bool:
    v = _cfg.get(name, None)
    if v is None and env:
        v = os.environ.get(env)
    if v is None:
        return default
    s = str(v).strip().lower()
    return s in ("1", "true", "yes", "on")

def _cfg_int(name: str, *, env: str | None = None, default: int = 0) -> int:
    v = _cfg.get(name, None)
    if v is None and env:
        v = os.environ.get(env)
    try:
        return int(v)
    except Exception:
        return default

def _cfg_float(name: str, default: float) -> float:
    try:
        return float(_cfg.get(name, default))
    except Exception:
        return default


# >>> RENDER PATHS PATCH START
def _resolve_receptor_and_outprefix(cfg, pdb_id, variant=None, ph_token=None, tag="renders", *, receptor_kind: str = "pdbqt"):
    p = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    if str(receptor_kind).lower() == "cleaned":
        receptor = p.receptor_cleaned_pdb(variant)
    else:
        receptor = p.receptor_pdbqt(variant, ph_token)
    outdir   = p.docked_variant_root(variant) / tag
    outdir.mkdir(parents=True, exist_ok=True)
    return p, str(receptor), outdir
# >>> RENDER PATHS PATCH END


def _normalize_candidate_path(candidate, base_dir: Optional[Path]) -> Optional[Path]:
    if candidate is None:
        return None
    cand_str = str(candidate).strip()
    if not cand_str:
        return None
    cand_path = Path(cand_str).expanduser()
    if cand_path.is_absolute() or base_dir is None:
        return cand_path.resolve()
    return (base_dir / cand_path).resolve()


def _default_outprefix_name(candidates: Sequence[Optional[str]], fallback: str) -> str:
    for cand in candidates:
        if not cand:
            continue
        stem = Path(str(cand)).stem
        if stem:
            return stem
    return fallback
_DEFER_MODE = _cfg_bool("DEFER_PYMOL", env="DEFER_PYMOL", default=False)
_QUEUE_PATH = (
    _cfg.get("PYMOL_DEFER_QUEUE")
    or os.environ.get("PYMOL_DEFER_QUEUE")
    or str(Path(os.environ.get("OVERALL_DIR", Path.cwd())) / "deferred_pymol_jobs.jsonl")
)
_Q_LOCK = threading.Lock()

def set_defer_mode(enable: bool, queue_path: Optional[str] = None):
    """Toggle deferral globally (preferred entry point for callers)."""
    global _DEFER_MODE, _QUEUE_PATH
    _DEFER_MODE = bool(enable)
    if queue_path:
        _QUEUE_PATH = str(queue_path)

def _enqueue_job(kind: str, payload: dict):
    Path(_QUEUE_PATH).parent.mkdir(parents=True, exist_ok=True)
    rec = {"kind": kind, **payload}
    line = json.dumps(rec, ensure_ascii=False)
    with _Q_LOCK:
        with open(_QUEUE_PATH, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")

def _iter_jobs(path: str):
    p = Path(path)
    if not p.is_file():
        return
    with open(p, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            try:
                yield json.loads(ln)
            except Exception:
                continue

def replay_deferred_jobs(queue_path: Optional[str] = None, max_workers: int = 2):
    qp = str(queue_path or _QUEUE_PATH)
    jobs = list(_iter_jobs(qp)) or []
    if not jobs:
        print("[capture_pose] no deferred PyMOL jobs to replay.")
        return 0

    print(f"[capture_pose] replaying {len(jobs)} deferred PyMOL jobs | workers={max_workers}")

    global _DEFER_MODE
    _prev = _DEFER_MODE
    _DEFER_MODE = False
    try:
        def _run(job):
            k = job.get("kind")
            if k == "three":
                return _render_three_views_with_pymol(
                    job.get("receptor_path"), job.get("ligand_paths_and_colors", []),
                    job.get("outprefix"),
                    cfg=job.get("cfg"),
                    pdb_id=job.get("pdb_id"),
                    variant=job.get("variant"),
                    ph_token=job.get("ph_token"),
                    tag=job.get("tag", "renders"),
                    label_top_n_res=job.get("label_top_n_res", 5),
                    label_cutoff=job.get("label_cutoff", 5.0),
                    viewport=tuple(job.get("viewport", (192, 144))),
                    hide_receptor=bool(job.get("hide_receptor", False)),
                )
            elif k == "native":
                return _render_native_on_original_pdb(
                    job.get("original_pdb"), job.get("outprefix"),
                    cfg=job.get("cfg"),
                    pdb_id=job.get("pdb_id"),
                    variant=job.get("variant"),
                    ph_token=job.get("ph_token"),
                    tag=job.get("tag", "renders"),
                    exclude_resns=job.get("exclude_resns", []),
                    viewport=tuple(job.get("viewport", (192, 144))),
                )

        workers = max(1, int(max_workers))
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(_run, jobs))
    finally:
        _DEFER_MODE = _prev
        # Clear queue after successful replay
        try:
            Path(qp).unlink()
        except Exception:
            pass
    print("[capture_pose] deferred PyMOL jobs complete.")
    return len(jobs)

def replay_deferred_jobs_mp(queue_path: Optional[str] = None,
                            max_workers: int = 99,
                            mode: str = "cli",
                            slow_ms: int = 15) -> int:
    """
    Multi-process renderer that replays JSONL queue without sharing a PyMOL session.
    Uses a single persistent pool to avoid per-batch spawn overhead.
    - mode: "cli" (pymol -cq) preferred; "pymol2" fallback per process.
    - slow_ms: only used for summary stats (no automatic backoff here).
    """
    qp = str(queue_path or _QUEUE_PATH)
    jobs = list(_iter_jobs(qp)) or []
    if not jobs:
        print(f"[render-replay] workers=0 jobs=0 mode={mode}")
        return 0

    global _DEFER_MODE
    prev = _DEFER_MODE
    _DEFER_MODE = False
    workers = max(1, int(max_workers))
    print(f"[render-replay] workers={workers} jobs={len(jobs)} mode={mode}")

    ok_count = fail_count = 0
    durations: List[int] = []

    # Use a stable mapping so we can print tags with the right job
    from concurrent.futures import ProcessPoolExecutor, as_completed
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fut2job = {ex.submit(_run_one_job_mp, job, mode): job for job in jobs}
        for fut in as_completed(fut2job):
            j = fut2job[fut]
            kind = j.get("kind", "?")
            tag = Path(j.get("outprefix", "?")).name
            try:
                ok, t_ms, err = fut.result()
            except Exception as e:
                ok, t_ms, err = False, 0, str(e)
            durations.append(t_ms)
            print(f"[render-done] {kind}:{tag} t_ms={t_ms} ok={int(ok)}"
                  + (f" err={str(err).splitlines()[0][:80]}" if err else ""))
            if ok:
                ok_count += 1
            else:
                fail_count += 1

    # Clear queue
    try:
        Path(qp).unlink()
    except Exception:
        pass
    _DEFER_MODE = prev

    # Summary
    if durations:
        avg = int(sum(durations) / len(durations))
        import statistics
        p50 = int(statistics.median(durations))
        p90 = int(sorted(durations)[max(0, int(0.9 * len(durations)) - 1)])
    else:
        avg = p50 = p90 = 0
    print(f"[render-summary] jobs={len(jobs)} ok={ok_count} fail={fail_count} avg_ms={avg} p50={p50} p90={p90} workers_used<={workers}")
    return ok_count





# ============================================================
# Utilities
# ============================================================

def _with_pymol():
    """Return PyMOL class if available, else None with a single clear message."""
    try:
        from pymol2 import PyMOL  # noqa: F401
        return PyMOL
    except Exception as e:
        print(f"[capture_pose] PyMOL not available; skipping renders: {e}")
        return None


def _safe_open_csv_for_write(target_path: Path, retries: int = 3, delay: float = 0.3):
    """
    Try to open `target_path` for writing. If it's locked (e.g., opened in Excel),
    retry a few times; then fall back to an alternate file:
      <stem>__alt1.csv, __alt2.csv, ...
    Returns: (file_handle, path_used)
    """
    target_path = Path(target_path)
    target_path.parent.mkdir(parents=True, exist_ok=True)

    last_err = None
    for _ in range(max(1, retries)):
        try:
            f = open(target_path, "w", newline="", encoding="utf-8")
            return f, target_path
        except OSError as e:
            last_err = e
            time.sleep(max(0.0, delay))

    # Fallback to alternate filename if still locked
    for k in range(1, 100):
        alt = target_path.with_name(f"{target_path.stem}__alt{k}{target_path.suffix}")
        try:
            f = open(alt, "w", newline="", encoding="utf-8")
            return f, alt
        except OSError:
            continue

    # If everything failed, re-raise the last error
    raise last_err if last_err else OSError(f"Unable to open {target_path}")


def _prefer_best_pdb(pose_path: str) -> str:
    """
    If a sibling '<stem>.best.pdb' exists next to a .pdbqt, prefer that for PyMOL GUI.
    Falls back to the input path otherwise.
    """
    p = Path(pose_path or "")
    if not p:
        return pose_path
    cand = p.with_name(p.stem + ".best.pdb")
    return str(cand) if cand.is_file() else pose_path


# ============================================================
# Selection helper used by your benchmark flow
# ============================================================

def pick_control_and_nearest_rdk(
    results_for_stage: Dict[str, Dict],
    raw_docked: Dict[str, str],
    control_stems_lower: Sequence[str],
) -> Tuple[Optional[str], Optional[str]]:
    """
    Returns (best_control_ligand_path, nearest_rdk_ligand_path).
    - results_for_stage: {lig_path: {"score": float, ...}}
    - raw_docked: {lig_path: best_pose_path.pdbqt}
    - control_stems_lower: set/list of native ligand stems, lowercase
    """
    ctrl_set = set(control_stems_lower)

    # best control (lowest score)
    best_ctrl, best_ctrl_score = None, float("inf")
    for lig, rec in results_for_stage.items():
        stem = Path(lig).stem.split("_stage")[0].lower()
        sc = rec.get("score")
        if stem in ctrl_set and isinstance(sc, (int, float)) and raw_docked.get(lig):
            if sc < best_ctrl_score:
                best_ctrl, best_ctrl_score = lig, sc

    # if control exists, choose RDK closest in score; else lowest-scoring RDK
    nearest_rdk, metric = None, float("inf")
    if best_ctrl is not None:
        for lig, rec in results_for_stage.items():
            stem = Path(lig).stem.split("_stage")[0].lower()
            sc = rec.get("score")
            if stem not in ctrl_set and isinstance(sc, (int, float)) and raw_docked.get(lig):
                d = abs(sc - best_ctrl_score)
                if d < metric:
                    nearest_rdk, metric = lig, d
    else:
        for lig, rec in results_for_stage.items():
            stem = Path(lig).stem.split("_stage")[0].lower()
            sc = rec.get("score")
            if stem not in ctrl_set and isinstance(sc, (int, float)) and raw_docked.get(lig):
                if sc < metric:
                    nearest_rdk, metric = lig, sc

    return best_ctrl, nearest_rdk


# ============================================================
# PyMOL rendering
# ============================================================

def _render_three_views_with_pymol(
    receptor_path: str | None,
    ligand_paths_and_colors: list[tuple[str, str, str]],
    outprefix: str | Path | None,
    *,
    cfg: Optional[dict] = None,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    tag: str = "renders",
    label_top_n_res: int = None,
    label_cutoff: float = None,
    viewport: tuple[int, int] | None = None,
    hide_receptor: bool = False,
) -> list[str]:
    """
    Saves <outprefix>_[front|side|top].png.
    Receptor shown as transparent surface; ligands as sticks; labels top-N closest residues (CA) within cutoff Å.
    If hide_receptor is True, receptor is not drawn (pair-only views).
    """
    # Resolve defaults from cfg once if not provided by caller
    if label_top_n_res is None:
        label_top_n_res = _cfg_int("LABEL_TOP_N_RES", default=5)
    if label_cutoff is None:
        label_cutoff = _cfg_float("LABEL_CUTOFF_A", 5.0)
    if viewport is None:
        viewport = ( _cfg_int("VIEWPORT_W", default=640), _cfg_int("VIEWPORT_H", default=480) )
    surface_transparency = _cfg_float("SURFACE_TRANSPARENCY", 0.30)

    cfg_data = cfg if cfg is not None else (_cfg or None)
    resolved_outprefix = None
    resolved_receptor = None
    router_parent: Optional[Path] = None
    outdir: Optional[Path] = None

    if cfg_data and pdb_id:
        try:
            _paths_obj, default_receptor, outdir = _resolve_receptor_and_outprefix(
                cfg_data, pdb_id, variant=variant, ph_token=ph_token, tag=tag, receptor_kind="cleaned"
            )
        except Exception:
            outdir = None
        else:
            router_parent = Path(default_receptor).parent
            resolved_receptor = _normalize_candidate_path(receptor_path, router_parent) if receptor_path else None
            if resolved_receptor is None:
                resolved_receptor = Path(default_receptor).resolve()

            outprefix_candidate = _normalize_candidate_path(outprefix, outdir)
            if outprefix_candidate is None:
                fallback_name = _default_outprefix_name(
                    [lp for lp, _, _ in ligand_paths_and_colors] + [obj for _, obj, _ in ligand_paths_and_colors],
                    fallback=f"{pdb_id}_render",
                )
                outprefix_candidate = (outdir / fallback_name).resolve()
            resolved_outprefix = outprefix_candidate
            resolved_outprefix.parent.mkdir(parents=True, exist_ok=True)

    if resolved_receptor is None:
        resolved_receptor = _normalize_candidate_path(receptor_path, router_parent)
        if resolved_receptor is None and receptor_path:
            resolved_receptor = Path(str(receptor_path)).expanduser().resolve()

    if resolved_outprefix is None:
        resolved_outprefix = _normalize_candidate_path(outprefix, outdir)
        if resolved_outprefix is None and outprefix:
            resolved_outprefix = Path(str(outprefix)).expanduser().resolve()
        if resolved_outprefix is not None:
            resolved_outprefix.parent.mkdir(parents=True, exist_ok=True)

    if resolved_receptor is None or not str(resolved_receptor):
        print("[capture_pose] Receptor path could not be resolved; skipping renders.")
        return []

    receptor_path = str(resolved_receptor)
    outprefix = str(resolved_outprefix) if resolved_outprefix is not None else str(outprefix or "")

    if not outprefix:
        print("[capture_pose] Output prefix could not be resolved; skipping renders.")
        return []

    if _DEFER_MODE:
        _enqueue_job("three", {
            "receptor_path": str(receptor_path),
            "ligand_paths_and_colors": [(str(p), str(n), str(c)) for (p, n, c) in ligand_paths_and_colors],
            "outprefix": str(outprefix),
            "cfg": cfg_data if cfg is not None else None,
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_token": ph_token,
            "tag": tag,
            "label_top_n_res": int(label_top_n_res),
            "label_cutoff": float(label_cutoff),
            "viewport": list(viewport),
            "hide_receptor": bool(hide_receptor),
        })
        return

    PyMOL = _with_pymol()
    if PyMOL is None:
        return

    if not Path(receptor_path).is_file():
        print(f"[capture_pose] Receptor missing: {receptor_path}")
        return

    with PyMOL() as pm:
        cmd = pm.cmd
        cmd.reinitialize()

        # Receptor (uniform gray, semi-transparent)
        cmd.load(receptor_path, "receptor")
        cmd.hide("everything", "receptor")
        if not hide_receptor:
            cmd.show("surface", "receptor")
            cmd.set_color("gray90", [230, 230, 230])
            cmd.color("gray90", "receptor")
            cmd.set("transparency", 0.30, "receptor")

        # Ligands (keep colors as provided)
        lig_objects: List[str] = []
        for lig_path, obj_name, color in ligand_paths_and_colors:
            if not lig_path or not Path(lig_path).is_file():
                continue
            lig_path_eff = _prefer_best_pdb(lig_path)
            if not Path(lig_path_eff).is_file():
                print(f"[capture_pose] missing ligand file: {lig_path_eff} (original: {lig_path})")
                continue

            cmd.load(lig_path_eff, obj_name)
            cmd.show("sticks", obj_name)
            cmd.color(color, obj_name)
            lig_objects.append(obj_name)

        lig_union = " or ".join(lig_objects) if lig_objects else "receptor"

        # Label nearest residues (sticks colored gray to avoid color noise)
        if lig_objects:
            cmd.select("active_site_all", f"receptor within {label_cutoff} of ({lig_union})")

            lig_model = cmd.get_model(lig_union)
            lig_coords = [(a.coord[0], a.coord[1], a.coord[2]) for a in lig_model.atom]

            def _min_dist_to_lig(x, y, z, coords):
                if not coords:
                    return float("inf")
                dx = x - coords[0][0]
                dy = y - coords[0][1]
                dz = z - coords[0][2]
                best = (dx * dx + dy * dy + dz * dz) ** 0.5
                for (lx, ly, lz) in coords[1:]:
                    dx = x - lx; dy = y - ly; dz = z - lz
                    d = (dx * dx + dy * dy + dz * dz) ** 0.5
                    if d < best:
                        best = d
                return best

            distances = []
            sel_ca = "active_site_all and name CA and (alt '' or alt A)"
            ca_model = cmd.get_model(sel_ca)
            for a in ca_model.atom:
                d = _min_dist_to_lig(a.coord[0], a.coord[1], a.coord[2], lig_coords)
                distances.append((a.model, a.chain, a.resi, a.resn, d))

            top_res = sorted(distances, key=lambda t: t[4])[:label_top_n_res] if distances else []
            if top_res:
                top_sel = " or ".join([f"receptor and chain {c} and resi {r}" for _, c, r, _, _ in top_res])
                cmd.select("top_site", top_sel)
                cmd.show("sticks", "top_site")
                cmd.color("gray", "top_site")
                cmd.label("top_site and name CA and (alt '' or alt A)", "resn + '-' + resi")

        # Views (unchanged performance settings)
        focus_sel = (lig_union if lig_objects else "receptor")
        cmd.zoom(focus_sel, 10)
        cmd.viewport(*viewport)
        cmd.set("antialias", 2)
        cmd.set("ray_opaque_background", 0)

        cmd.sync(); cmd.refresh()
        cmd.png(f"{outprefix}_front.png", ray=0)

        cmd.sync(); cmd.refresh()
        cmd.turn("y", 90)
        cmd.png(f"{outprefix}_side.png", ray=0)

        cmd.sync(); cmd.refresh()
        cmd.turn("x", 90)
        cmd.png(f"{outprefix}_top.png", ray=0)





def _render_native_on_original_pdb(
    original_pdb: str | None,
    outprefix: str | Path | None,
    *,
    cfg: Optional[dict] = None,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    tag: str = "renders",
    exclude_resns: Sequence[str] = tuple(EXCLUDE_HET_IDS),
    viewport: Tuple[int, int] = (192, 144),
) -> None:
    """
    Render original PDB with native ligand(s): polymer surface + HETATM (minus excludes) as sticks.
    Saves <outprefix>_[front|side|top].png.
    """
    cfg_data = cfg if cfg is not None else (_cfg or None)
    resolved_outprefix = None
    resolved_original = None
    outdir: Optional[Path] = None
    paths_obj = None

    if cfg_data and pdb_id:
        try:
            paths_obj, _default_receptor, outdir = _resolve_receptor_and_outprefix(
                cfg_data, pdb_id, variant=variant, ph_token=ph_token, tag=tag, receptor_kind="cleaned"
            )
        except Exception:
            outdir = None
        else:
            resolved_outprefix = _normalize_candidate_path(outprefix, outdir)
            if resolved_outprefix is None:
                fallback = _default_outprefix_name([original_pdb], fallback=f"{pdb_id}_native")
                resolved_outprefix = (outdir / fallback).resolve()
            resolved_outprefix.parent.mkdir(parents=True, exist_ok=True)

    if resolved_outprefix is None and outprefix:
        resolved_outprefix = _normalize_candidate_path(outprefix, outdir)
        if resolved_outprefix is None:
            resolved_outprefix = Path(str(outprefix)).expanduser().resolve()
        resolved_outprefix.parent.mkdir(parents=True, exist_ok=True)

    if paths_obj is not None:
        base_dirs = [paths_obj.root_pdb_dir, paths_obj.processed_root, paths_obj.input_root]
    else:
        base_dirs = []

    if original_pdb:
        resolved_original = _normalize_candidate_path(original_pdb, None)
        if resolved_original is None or not resolved_original.is_file():
            for base in base_dirs:
                cand = _normalize_candidate_path(original_pdb, base)
                if cand is not None and cand.is_file():
                    resolved_original = cand
                    break
            if resolved_original is None:
                resolved_original = Path(str(original_pdb)).expanduser().resolve()

    outprefix_value = str(resolved_outprefix) if resolved_outprefix is not None else str(outprefix or "")
    if not outprefix_value:
        print("[capture_pose] Output prefix could not be resolved; skipping native render.")
        return

    if _DEFER_MODE:
        _enqueue_job("native", {
            "original_pdb": str(resolved_original) if resolved_original else str(original_pdb or ""),
            "outprefix": outprefix_value,
            "cfg": cfg_data if cfg is not None else None,
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_token": ph_token,
            "tag": tag,
            "exclude_resns": list(exclude_resns or []),
            "viewport": list(viewport),
        })
        return

    PyMOL = _with_pymol()
    if PyMOL is None:
        return

    original_path = resolved_original if resolved_original is not None else (Path(str(original_pdb)).expanduser() if original_pdb else None)
    if not original_path or not Path(original_path).is_file():
        print(f"[capture_pose] Original PDB missing: {original_pdb}")
        return

    with PyMOL() as pm:
        cmd = pm.cmd
        cmd.reinitialize()

        cmd.load(str(original_path), "orig")
        cmd.hide("everything")

        # Protein polymer: uniform gray + transparency; no auto-coloring
        cmd.show("surface", "orig and polymer")
        cmd.set_color("gray90", [230, 230, 230])
        cmd.color("gray90", "orig and polymer")
        cmd.set("transparency", 0.30, "orig and polymer")

        # Native ligands (keep visible & colored as before)
        excl = "+".join(sorted(set(exclude_resns or [])))
        cmd.select("native_lig", f"(hetatm and not polymer and not solvent) and not resn {excl}")
        if cmd.count_atoms("native_lig") > 0:
            cmd.show("sticks", "native_lig")
            cmd.color("green", "native_lig")

            # Nearby residues: keep sticks but gray to avoid visual noise
            cmd.select("near_native", "orig within 5 of native_lig and polymer.protein")
            cmd.show("sticks", "near_native")
            cmd.color("gray", "near_native")
            cmd.label("near_native and name CA", "resn + '-' + resi")

            focus_sel = "native_lig or near_native"
        else:
            focus_sel = "orig and polymer"

        # Views (unchanged)
        cmd.zoom(focus_sel, 10)
        cmd.viewport(*viewport)
        cmd.set("antialias", 2)
        cmd.set("ray_opaque_background", 0)
        cmd.sync(); cmd.refresh()
        cmd.png(f"{outprefix_value}_front.png", ray=0)
        cmd.sync(); cmd.refresh()
        cmd.turn("y", 90)
        cmd.png(f"{outprefix_value}_side.png", ray=0)
        cmd.sync(); cmd.refresh()
        cmd.turn("x", 90)
        cmd.png(f"{outprefix_value}_top.png", ray=0)


# ============================================================
# PML writer + launch/fallback
# ============================================================
from pathlib import Path
from typing import Tuple
from string import Template
def write_multiview_pml(
    receptor_path: str | None,
    control_path: str,
    rdk_path: str,
    out_pml: Path,
    label_top_n_res: int = 5,
    label_cutoff: float = 5.0,
    outprefix: str | Path = "",
    viewport: Tuple[int, int] = (192, 144),
    *,
    cfg: Optional[dict] = None,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    tag: str = "renders",
) -> Path:
    cfg_data = cfg if cfg is not None else (_cfg or None)
    resolved_outprefix: Optional[Path] = None
    resolved_receptor: Optional[Path] = None
    outdir: Optional[Path] = None

    if cfg_data and pdb_id:
        try:
            _paths_obj, default_receptor, outdir = _resolve_receptor_and_outprefix(
                cfg_data, pdb_id, variant=variant, ph_token=ph_token, tag=tag, receptor_kind="cleaned"
            )
        except Exception:
            outdir = None
        else:
            resolved_receptor = _normalize_candidate_path(receptor_path, Path(default_receptor).parent) if receptor_path else None
            if resolved_receptor is None:
                resolved_receptor = Path(default_receptor).resolve()

            resolved_outprefix = _normalize_candidate_path(outprefix, outdir)
            if resolved_outprefix is None:
                fallback = _default_outprefix_name([control_path, rdk_path], fallback=f"{pdb_id}_render")
                resolved_outprefix = (outdir / fallback).resolve()
            resolved_outprefix.parent.mkdir(parents=True, exist_ok=True)

    if resolved_receptor is None and receptor_path:
        resolved_receptor = _normalize_candidate_path(receptor_path, None)
        if resolved_receptor is None:
            resolved_receptor = Path(str(receptor_path)).expanduser().resolve()

    if resolved_outprefix is None:
        resolved_outprefix = _normalize_candidate_path(outprefix, outdir)
        if resolved_outprefix is None and outprefix:
            resolved_outprefix = Path(str(outprefix)).expanduser().resolve()
        if resolved_outprefix is not None:
            resolved_outprefix.parent.mkdir(parents=True, exist_ok=True)

    out_pml = Path(out_pml)
    if outdir is not None and not out_pml.is_absolute():
        out_pml = (outdir / out_pml).resolve()
    out_pml.parent.mkdir(parents=True, exist_ok=True)

    receptor_path = str(resolved_receptor) if resolved_receptor is not None else str(receptor_path)
    outprefix = str(resolved_outprefix) if resolved_outprefix is not None else str(outprefix)

    def _posix(p: str) -> str:
        return (Path(p).resolve().as_posix() if p else "")

    receptor_posix = _posix(receptor_path)
    control_posix  = _posix(control_path) if control_path else ""
    rdk_posix      = _posix(rdk_path) if rdk_path else ""

    opref = Path(outprefix) if outprefix else out_pml.with_suffix("")
    if not opref.is_absolute():
        opref = opref.resolve()
    outprefix = str(opref)
    out_front = (opref.with_name(opref.name + "_front")).resolve().as_posix()
    out_side  = (opref.with_name(opref.name + "_side")).resolve().as_posix()
    out_top   = (opref.with_name(opref.name + "_top")).resolve().as_posix()
    Path(out_front).parent.mkdir(parents=True, exist_ok=True)

    parts = []
    if control_posix: parts.append("control")
    if rdk_posix:     parts.append("rdk")
    lig_union = " or ".join(parts) if parts else "receptor"

    tpl = Template(r"""
reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 1

set async_builds, off
set defer_builds_mode, 3
feedback disable, all, results
feedback disable, all, actions

load $RECEPTOR, receptor
hide everything, receptor
show surface, receptor
set transparency, 0.30
# No custom set_color needed; gray90 is built-in
color gray90, receptor
set surface_color, gray90, receptor
set transparency, 0.30, receptor
set two_sided_lighting, on
set ambient, 0.4
set specular, 0.2
set spec_power, 150
set light_count, 8
set depth_cue, on
set fog_start, 0.6


$LOAD_CONTROL
$LOAD_RDK

# Re-assert uniform gray on protein after all loads (belt-and-suspenders)
color gray90, receptor and surface

select lig_union, ($LIG_UNION)

# Label top-N nearby residues (sticks colored gray)
select near_res, (receptor within $CUTOFF of lig_union) and polymer.protein and name CA and (alt '' or alt A)
python
TOPN = $TOPN
import math
from pymol import cmd
lig_model = cmd.get_model("lig_union")
lig_coords = [a.coord for a in lig_model.atom]
def _mindist(x,y,z,coords):
    if not coords: return float("inf")
    best = float("inf")
    for lx,ly,lz in coords:
        dx=x-lx; dy=y-ly; dz=z-lz
        d=dx*dx+dy*dy+dz*dz
        if d<best: best=d
    return math.sqrt(best)
dlist=[]
for (obj_name, atom_index) in cmd.index("near_res"):
    m=cmd.get_model(f"{obj_name} and index {atom_index}")
    if not m.atom: continue
    a=m.atom[0]
    dist=_mindist(a.coord[0],a.coord[1],a.coord[2],lig_coords)
    dlist.append((a.model,a.chain,a.resi,dist))
dlist.sort(key=lambda t:t[3])
keep,seen=[],set()
for (m,c,i,_) in dlist:
    key=(m,c,i)
    if key in seen: continue
    seen.add(key); keep.append(key)
    if len(keep)>=TOPN: break
sel=" or ".join([f"{m}//{c}/{i}" for (m,c,i) in keep])
if sel:
    cmd.select("top_site", sel)
    cmd.show("sticks", "top_site")
    cmd.color("gray", "top_site")
    cmd.label("top_site and name CA and (alt '' or alt A)", "resn + '-' + resi")
python end

orient lig_union
zoom lig_union, 10
scene front, store
turn y, 90
scene side, store
turn x, 90
scene top, store

viewport $W, $H
png $OUT_FRONT, ray=0
turn y, 90
png $OUT_SIDE, ray=0
turn x, 90
png $OUT_TOP, ray=0
""")

    load_control = ""
    if control_posix:
        load_control = '\n'.join([
            f'load {control_posix}, control',
            'show sticks, control',
            'color blue, control',
        ])
    load_rdk = ""
    if rdk_posix:
        load_rdk = '\n'.join([
            f'load {rdk_posix}, rdk',
            'show sticks, rdk',
            'color orange, rdk',
        ])

    w, h = int(viewport[0]), int(viewport[1])
    pml_text = tpl.substitute(
        RECEPTOR=receptor_posix,
        LOAD_CONTROL=load_control,
        LOAD_RDK=load_rdk,
        LIG_UNION=lig_union,
        CUTOFF=f"{label_cutoff:.2f}",
        TOPN=str(label_top_n_res),
        W=str(w), H=str(h),
        OUT_FRONT=out_front,
        OUT_SIDE=out_side,
        OUT_TOP=out_top,
    )
    out_pml.write_text(pml_text, encoding="utf-8")
    return out_pml


def write_native_pml(original_pdb: str, outprefix: Path, exclude_resns: Sequence[str], viewport=(192,144)) -> Path:
    outprefix = Path(outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)
    excl = " ".join(sorted(set(exclude_resns or [])))
    tpl = Template(r"""
reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 1

set async_builds, off
set defer_builds_mode, 3
feedback disable, all, results
feedback disable, all, actions

load "$PDB", orig
hide everything
show surface, orig and polymer
set transparency, 0.30, orig and polymer
set_color gray90, [230,230,230]
color gray90, orig and polymer
python
python end
select native_lig, (hetatm and not polymer and not solvent) and not resn $EXCL
if (count_atoms("native_lig")>0) {
    show sticks, native_lig
    color green, native_lig
    select near_native, orig within 5 of native_lig and polymer.protein
    show sticks, near_native
    color gray, near_native
    label near_native and name CA, resn + "-" + resi
    zoom native_lig or near_native, 10
} else {
    zoom orig and polymer, 10
}
viewport $W, $H
png "$OUT_front", ray=0
turn y, 90
png "$OUT_side", ray=0
turn x, 90
png "$OUT_top",  ray=0
""")
    w, h = int(viewport[0]), int(viewport[1])
    pml = outprefix.with_suffix(".native.pml")
    pml.write_text(tpl.substitute(
        PDB=Path(original_pdb).resolve().as_posix(),
        EXCL=",".join(exclude_resns or []),
        W=str(w), H=str(h),
        OUT_front=(str(outprefix) + "_front.png").replace("\\","/"),
        OUT_side=(str(outprefix) + "_side.png").replace("\\","/"),
        OUT_top =(str(outprefix) + "_top.png").replace("\\","/"),
    ), encoding="utf-8")
    return pml
def _job_to_pml(job: dict, tmpdir: Path) -> Optional[Path]:
    kind = job.get("kind")
    if kind == "three":
        # Map list into control/rdk slots by name if available (fallback = rdk only)
        ctrl_path, rdk_path = "", ""
        for p, obj, _col in job.get("ligand_paths_and_colors", []):
            name = (obj or "").lower()
            if "control" in name and not ctrl_path:
                ctrl_path = p
            elif "rdk" in name and not rdk_path:
                rdk_path = p
        # if only one ligand, treat it as rdk
        if not ctrl_path and not rdk_path and job.get("ligand_paths_and_colors"):
            rdk_path = job["ligand_paths_and_colors"][0][0]
        out_pml = tmpdir / (Path(job["outprefix"]).name + ".pml")
        return write_multiview_pml(
            receptor_path=job["receptor_path"],
            control_path=ctrl_path or "",
            rdk_path=rdk_path or "",
            out_pml=out_pml,
            label_top_n_res=int(job.get("label_top_n_res", 5)),
            label_cutoff=float(job.get("label_cutoff", 5.0)),
            outprefix=str(job.get("outprefix", "")),
            viewport=tuple(job.get("viewport", (192, 144))),
        )
    elif kind == "native":
        outprefix = Path(job["outprefix"])
        return write_native_pml(
            original_pdb=job["original_pdb"],
            outprefix=outprefix,
            exclude_resns=job.get("exclude_resns", []),
            viewport=tuple(job.get("viewport", (192,144))),
        )
    return None
def _run_one_job_mp(job: dict, mode: str = "cli") -> tuple:
    """
    Returns: (ok:bool, t_ms:int, err:str|None)
    """
    t0 = time.time()
    err_first = None
    try:
        if mode == "cli":
            with tempfile.TemporaryDirectory() as td:
                pml = _job_to_pml(job, Path(td))
                if pml is None:
                    raise RuntimeError("unsupported job")
                rc = render_pml_headless(pml)
                ok = (rc == 0)
        else:
            # process-local pymol2 fallback (still one job per process)
            k = job.get("kind")
            if k == "three":
                _render_three_views_with_pymol(
                    job.get("receptor_path"), job.get("ligand_paths_and_colors", []), job.get("outprefix"),
                    cfg=job.get("cfg"),
                    pdb_id=job.get("pdb_id"),
                    variant=job.get("variant"),
                    ph_token=job.get("ph_token"),
                    tag=job.get("tag", "renders"),
                    label_top_n_res=job.get("label_top_n_res", 5),
                    label_cutoff=job.get("label_cutoff", 5.0),
                    viewport=tuple(job.get("viewport", (192,144))),
                    hide_receptor=bool(job.get("hide_receptor", False)),
                )
                ok = True
            elif k == "native":
                _render_native_on_original_pdb(
                    job.get("original_pdb"), job.get("outprefix"),
                    cfg=job.get("cfg"),
                    pdb_id=job.get("pdb_id"),
                    variant=job.get("variant"),
                    ph_token=job.get("ph_token"),
                    tag=job.get("tag", "renders"),
                    exclude_resns=job.get("exclude_resns", []),
                    viewport=tuple(job.get("viewport", (192,144))),
                )
                ok = True
            else:
                ok = False
        t_ms = int((time.time() - t0) * 1000)
        return (ok, t_ms, err_first)
    except Exception as e:
        t_ms = int((time.time() - t0) * 1000)
        return (False, t_ms, str(e))

def launch_pymol_with_pml(pml_path: Path, pymol_exe: Optional[str] = None) -> None:
    """
    Launch PyMOL GUI with the given .pml (non-blocking).
    If pymol_exe is None, tries 'pymol' on PATH.
    """
    exe = pymol_exe or "pymol"
    try:
        subprocess.Popen([exe, str(pml_path)], shell=False)
    except FileNotFoundError:
        print(f"[capture_pose] Could not find PyMOL executable '{exe}'. Open manually: {pml_path}")


def render_pml_headless(pml_path: Path, pymol_exe: Optional[str] = None) -> int:
    """
    Render a .pml without pymol2 using PyMOL CLI.
    The .pml must include PNG commands; returns process returncode.
    """
    exe = pymol_exe or "pymol"
    try:
        import shlex, subprocess
        pml_abs = Path(pml_path).resolve()
        cmd = [exe, "-cq", str(pml_abs)]
        print("[render-cli] " + " ".join(shlex.quote(x) for x in cmd))
        return subprocess.run(cmd, check=False).returncode
    except FileNotFoundError:
        print(f"[capture_pose] PyMOL CLI not found: '{exe}'")
        return 127



# ============================================================
# Simple CLI (one-ligand) for quick testing
# ============================================================
def _cli_render_active_site(
    receptor: str,
    ligand: str,
    outprefix: str,
    top_n_residues: int = 5,
    proximity_cutoff: float = 5.0,
    viewport: Tuple[int, int] = (192, 144),
    **kwargs,
) -> None:
    """
    Headless-friendly renderer (NO RAY TRACING):
      - Loads receptor/ligand
      - Computes ligand centroid (heavy atoms)
      - Finds receptor CA atoms within cutoff; ranks nearest
      - Shows ligand sticks, highlights residues, saves 3 PNGs
      - Applies optional color-blind–friendly palette and colors

    Optional kwargs (all optional; ignored if missing):
      - palette_defs: Dict[str, List[int]]
      - protein_color: str
      - protein_transparency: float
      - ligand_color: str
      - pocket_color: str
    """
    import os
    from typing import List, Dict
    try:
        from pymol import cmd
    except Exception as e:
        raise RuntimeError("PyMOL (pymol2) is required to render.") from e

    rec_obj = "receptor"
    lig_obj = "ligand"

    cmd.reinitialize()

    # Load
    cmd.load(receptor, rec_obj)
    cmd.load(ligand, lig_obj)

    # Viewport unchanged
    try:
        w, h = int(viewport[0]), int(viewport[1])
        cmd.viewport(w, h)
    except Exception:
        pass

    # Apply palette definitions first (so names are usable)
    palette_defs = kwargs.get("palette_defs") or {}
    if isinstance(palette_defs, dict):
        for name, rgb in palette_defs.items():
            if isinstance(rgb, (list, tuple)) and len(rgb) == 3:
                cmd.set_color(str(name), [int(rgb[0]), int(rgb[1]), int(rgb[2])])

    # Ensure receptor is uniformly gray and semi-transparent by default
    cmd.hide("everything", rec_obj)
    cmd.show("surface", rec_obj)
    cmd.set_color("gray90", [230, 230, 230])
    cmd.color("gray90", rec_obj)
    cmd.set("transparency", 0.30, rec_obj)


    protein_transparency = kwargs.get("protein_transparency")
    if protein_transparency is not None:
        try:
            cmd.set("transparency", float(protein_transparency), rec_obj)
        except Exception:
            pass

    # Ligand color (keep if provided; otherwise default elsewhere)
    lig_color = kwargs.get("ligand_color")
    if lig_color:
        cmd.color(str(lig_color), lig_obj)

    pocket_color = kwargs.get("pocket_color")
    if pocket_color:
        try:
            cmd.color(str(pocket_color), "pocket")
        except Exception:
            pass

    # Sticks for ligand; protein already surface
    cmd.show("sticks", lig_obj)

    # Centroid + labels (unchanged)
    cmd.select("lig_heavy_tmp", f"({lig_obj}) and not elem H")
    coords = []
    cmd.iterate_state(1, "lig_heavy_tmp", "coords.append([x,y,z])", space={"coords": coords})
    if coords:
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        cmd.pseudoatom("lig_centroid", pos=[cx, cy, cz])
    else:
        cmd.orient(lig_obj)

    cutoff = float(proximity_cutoff)
    cmd.select("near_CA", f"byres ({rec_obj} and name CA within {cutoff} of {lig_obj})")
    cmd.label("near_CA and name CA", '"%s-%s" % (resn, resi)')

    cmd.orient(lig_obj)

    base = outprefix
    if os.path.isdir(outprefix) or outprefix.endswith(os.sep):
        base = os.path.join(outprefix, "top_pose")
    os.makedirs(os.path.dirname(base) or ".", exist_ok=True)

    cmd.sync(); cmd.refresh()
    cmd.png(base + "_front.png", ray=0, width=w, height=h)
    cmd.sync(); cmd.refresh()
    cmd.turn("y", 90)
    cmd.png(base + "_side.png", ray=0, width=w, height=h)
    cmd.sync(); cmd.refresh()
    cmd.turn("x", 90)
    cmd.png(base + "_top.png", ray=0, width=w, height=h)

    try:
        cmd.delete("lig_heavy_tmp")
        cmd.delete("lig_centroid")
    except Exception:
        pass



def capture_pose(receptor_path, ligand_path, out_path_or_dir, **kwargs):
    """
    Backward-compatible default: ligand stays magenta, existing behavior unchanged.

    Optional controls (only applied if provided):
      - ligand_color: explicit PyMOL color name or custom color name
      - ligand_role: "reference" | "candidate"   (maps to blue/orange under okabe_ito)
      - color_scheme: "okabe_ito"                 (defines safe palette names)
      - palette_defs: dict like {"orangeOI":[230,159,0], ...} (override/extend)
      - protein_color: e.g., "gray90"
      - protein_transparency: float 0..1
      - pocket_color: e.g., "tealOI"
    """
    import os
    from pathlib import Path

    # Built-in color-blind–friendly palette (only used if requested)
    base_palettes = {
        "okabe_ito": {
            "custom_defs": {
                "orangeOI": [230, 159,   0],  # candidate
                "blueOI":   [  0, 114, 178],  # reference
                "tealOI":   [  0, 158, 115],
                "yellowOI": [240, 228,  66],
                "vermOI":   [213,  94,   0],
                "purpleOI": [204, 121, 167],
            },
            "ligand_colors": {
                "reference": "blueOI",
                "candidate": "orangeOI",
            },
        }
    }

    # Did the caller request any custom color behavior?
    wants_custom = any(k in kwargs for k in (
        "ligand_color", "ligand_role", "color_scheme", "palette_defs",
        "protein_color", "protein_transparency", "pocket_color"
    ))

    # Defaults that preserve current behavior
    chosen_ligand_color = "magenta"

    # Resolve custom ligand color if requested
    if wants_custom:
        scheme_key = kwargs.get("color_scheme")
        scheme = base_palettes.get(scheme_key, {}) if scheme_key else {}
        role = (kwargs.get("ligand_role") or "candidate").lower()
        # Direct override wins
        if "ligand_color" in kwargs:
            chosen_ligand_color = kwargs["ligand_color"]
        # Otherwise pick from scheme if provided
        elif scheme and "ligand_colors" in scheme:
            chosen_ligand_color = scheme["ligand_colors"].get(role, scheme["ligand_colors"].get("candidate", "magenta"))

    # Compute output prefix if a directory is given
    outprefix = out_path_or_dir
    try:
        if os.path.isdir(out_path_or_dir):
            outprefix = str(Path(out_path_or_dir) / Path(ligand_path).stem)
    except Exception:
        pass

    # Build payload for deferred renderer
    payload = {
        "receptor_path": str(receptor_path),
        "ligand_paths_and_colors": [(str(ligand_path), "ligand", chosen_ligand_color)],
        "outprefix": str(outprefix),
        "label_top_n_res": int(kwargs.get("top_n_residues", 6)),
        "label_cutoff": float(kwargs.get("proximity_cutoff", 4.5)),
        "viewport": list(kwargs.get("viewport", (192, 144))),
        "hide_receptor": False,
    }

    # Attach extras only if custom requested
    if wants_custom:
        # Merge built-in palette with user overrides
        palette_defs = {}
        if kwargs.get("color_scheme") in base_palettes:
            palette_defs.update(base_palettes[kwargs["color_scheme"]]["custom_defs"])
        palette_defs.update(kwargs.get("palette_defs", {}))
        if palette_defs:
            payload["palette_defs"] = palette_defs
        if "protein_color" in kwargs:
            payload["protein_color"] = kwargs["protein_color"]
        if "protein_transparency" in kwargs:
            payload["protein_transparency"] = float(kwargs["protein_transparency"])
        if "pocket_color" in kwargs:
            payload["pocket_color"] = kwargs["pocket_color"]

    # Defer path (unchanged API)
    if _DEFER_MODE:
        _enqueue_job("three", payload)
        return

    # Immediate path
    outpath_final = (str(outprefix) if not os.path.isdir(out_path_or_dir)
                     else str(Path(out_path_or_dir) / Path(ligand_path).stem))

    appearance_kwargs = {}
    if wants_custom:
        if "protein_color" in payload:
            appearance_kwargs["protein_color"] = payload["protein_color"]
        if "protein_transparency" in payload:
            appearance_kwargs["protein_transparency"] = payload["protein_transparency"]
        if "palette_defs" in payload:
            appearance_kwargs["palette_defs"] = payload["palette_defs"]
        if "pocket_color" in payload:
            appearance_kwargs["pocket_color"] = payload["pocket_color"]
        appearance_kwargs["ligand_color"] = chosen_ligand_color

    try:
        return _cli_render_active_site(
            str(receptor_path),
            str(ligand_path),
            outpath_final,
            top_n_residues=kwargs.get("top_n_residues", 6),
            proximity_cutoff=kwargs.get("proximity_cutoff", 4.5),
            viewport=kwargs.get("viewport", (192, 144)),
            **appearance_kwargs
        )
    except TypeError:
        # Fallback if old signature: ignore appearance kwargs
        return _cli_render_active_site(
            str(receptor_path),
            str(ligand_path),
            outpath_final,
            top_n_residues=kwargs.get("top_n_residues", 6),
            proximity_cutoff=kwargs.get("proximity_cutoff", 4.5),
            viewport=kwargs.get("viewport", (192, 144)),
        )

