"""PML generation and deferred replay helpers for pose capture."""

from __future__ import annotations

import subprocess
import tempfile
import time
from pathlib import Path
from string import Template
from typing import List, Optional, Sequence, Tuple

from docking import capture_pose_common as common
from docking.capture_pose_queue import _iter_jobs
from docking.capture_pose_render import _render_native_on_original_pdb, _render_three_views_with_pymol

def _dispatch_render_job(job: dict) -> None:
    kind = job.get("kind")
    if kind == "three":
        _render_three_views_with_pymol(
            job.get("receptor_path"),
            job.get("ligand_paths_and_colors", []),
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
        return
    if kind == "native":
        _render_native_on_original_pdb(
            job.get("original_pdb"),
            job.get("outprefix"),
            cfg=job.get("cfg"),
            pdb_id=job.get("pdb_id"),
            variant=job.get("variant"),
            ph_token=job.get("ph_token"),
            tag=job.get("tag", "renders"),
            exclude_resns=job.get("exclude_resns", []),
            viewport=tuple(job.get("viewport", (192, 144))),
        )
        return
    raise RuntimeError(f"unsupported job kind: {kind}")


def _clear_queue(path: str) -> None:
    try:
        Path(path).unlink()
    except Exception:
        pass


def replay_deferred_jobs(queue_path: Optional[str] = None, max_workers: int = 2):
    qp = str(queue_path or common._QUEUE_PATH)
    jobs = list(_iter_jobs(qp)) or []
    if not jobs:
        print("[capture_pose] no deferred PyMOL jobs to replay.")
        return 0

    print(f"[capture_pose] replaying {len(jobs)} deferred PyMOL jobs | workers={max_workers}")

    prev = common._DEFER_MODE
    common._DEFER_MODE = False
    try:
        workers = max(1, int(max_workers))
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(_dispatch_render_job, jobs))
    finally:
        common._DEFER_MODE = prev
        _clear_queue(qp)
    print("[capture_pose] deferred PyMOL jobs complete.")
    return len(jobs)

def _duration_summary(durations: List[int]) -> tuple[int, int, int]:
    if not durations:
        return 0, 0, 0
    avg = int(sum(durations) / len(durations))
    import statistics

    p50 = int(statistics.median(durations))
    p90 = int(sorted(durations)[max(0, int(0.9 * len(durations)) - 1)])
    return avg, p50, p90


def replay_deferred_jobs_mp(
    queue_path: Optional[str] = None,
    max_workers: int = 99,
    mode: str = "cli",
    slow_ms: int = 15,
) -> int:
    """
    Multi-process renderer that replays JSONL queue without sharing a PyMOL session.
    Uses a single persistent pool to avoid per-batch spawn overhead.
    - mode: "cli" (pymol -cq) preferred; "pymol2" fallback per process.
    - slow_ms: only used for summary stats (no automatic backoff here).
    """
    qp = str(queue_path or common._QUEUE_PATH)
    jobs = list(_iter_jobs(qp)) or []
    if not jobs:
        print(f"[render-replay] workers=0 jobs=0 mode={mode}")
        return 0

    prev = common._DEFER_MODE
    common._DEFER_MODE = False
    workers = max(1, int(max_workers))
    print(f"[render-replay] workers={workers} jobs={len(jobs)} mode={mode}")

    ok_count = fail_count = 0
    durations: List[int] = []

    from concurrent.futures import ProcessPoolExecutor, as_completed

    with ProcessPoolExecutor(max_workers=workers) as ex:
        fut2job = {ex.submit(_run_one_job_mp, job, mode): job for job in jobs}
        for fut in as_completed(fut2job):
            job = fut2job[fut]
            tag = Path(job.get("outprefix", "?")).name
            kind = job.get("kind", "?")
            try:
                ok, t_ms, err = fut.result()
            except Exception as exc:
                ok, t_ms, err = False, 0, str(exc)
            durations.append(t_ms)
            print(
                f"[render-done] {kind}:{tag} t_ms={t_ms} ok={int(ok)}"
                + (f" err={str(err).splitlines()[0][:80]}" if err else "")
            )
            if ok:
                ok_count += 1
            else:
                fail_count += 1

    _clear_queue(qp)
    common._DEFER_MODE = prev
    avg, p50, p90 = _duration_summary(durations)
    print(
        f"[render-summary] jobs={len(jobs)} ok={ok_count} fail={fail_count} avg_ms={avg} p50={p50} p90={p90} workers_used<={workers}"
    )
    return ok_count

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
    cfg_data = common._cfg_data_or_default(cfg)
    default_receptor = None
    outdir: Optional[Path] = None
    if cfg_data and pdb_id:
        try:
            _paths_obj, default_receptor, outdir = common._resolve_receptor_and_outprefix(
                cfg_data,
                pdb_id,
                variant=variant,
                ph_token=ph_token,
                tag=tag,
                receptor_kind="cleaned",
            )
        except Exception:
            outdir = None
            default_receptor = None

    resolved_receptor = common._resolve_receptor_path(receptor_path, default_receptor)
    resolved_outprefix = common._resolve_output_prefix(
        outprefix,
        outdir,
        fallback_candidates=[control_path, rdk_path],
        fallback_name=f"{pdb_id}_render" if pdb_id else "render",
    )

    out_pml = Path(out_pml)
    if outdir is not None and not out_pml.is_absolute():
        out_pml = (outdir / out_pml).resolve()
    out_pml.parent.mkdir(parents=True, exist_ok=True)

    receptor_value = str(resolved_receptor) if resolved_receptor is not None else str(receptor_path or "")
    outprefix_value = str(resolved_outprefix) if resolved_outprefix is not None else str(outprefix)

    def _posix(path_value: str) -> str:
        return Path(path_value).resolve().as_posix() if path_value else ""

    receptor_posix = _posix(receptor_value)
    control_posix = _posix(control_path) if control_path else ""
    rdk_posix = _posix(rdk_path) if rdk_path else ""

    opref = Path(outprefix_value) if outprefix_value else out_pml.with_suffix("")
    if not opref.is_absolute():
        opref = opref.resolve()
    Path(opref).parent.mkdir(parents=True, exist_ok=True)
    out_front = (opref.with_name(opref.name + "_front")).resolve().as_posix()
    out_side = (opref.with_name(opref.name + "_side")).resolve().as_posix()
    out_top = (opref.with_name(opref.name + "_top")).resolve().as_posix()

    lig_parts = [name for name, value in (("control", control_posix), ("rdk", rdk_posix)) if value]
    lig_union = " or ".join(lig_parts) if lig_parts else "receptor"

    tpl = Template(
        r"""
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
color gray90, receptor and surface
select lig_union, ($LIG_UNION)
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
"""
    )

    def _load_block(path_value: str, obj_name: str, color: str) -> str:
        if not path_value:
            return ""
        return "\n".join(
            [
                f"load {path_value}, {obj_name}",
                f"show sticks, {obj_name}",
                f"color {color}, {obj_name}",
            ]
        )

    w, h = int(viewport[0]), int(viewport[1])
    out_pml.write_text(
        tpl.substitute(
            RECEPTOR=receptor_posix,
            LOAD_CONTROL=_load_block(control_posix, "control", "blue"),
            LOAD_RDK=_load_block(rdk_posix, "rdk", "orange"),
            LIG_UNION=lig_union,
            CUTOFF=f"{label_cutoff:.2f}",
            TOPN=str(label_top_n_res),
            W=str(w),
            H=str(h),
            OUT_FRONT=out_front,
            OUT_SIDE=out_side,
            OUT_TOP=out_top,
        ),
        encoding="utf-8",
    )
    return out_pml

def write_native_pml(
    original_pdb: str,
    outprefix: Path,
    exclude_resns: Sequence[str],
    viewport=(192, 144),
) -> Path:
    outprefix = Path(outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)
    tpl = Template(
        r"""
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
"""
    )
    w, h = int(viewport[0]), int(viewport[1])
    pml = outprefix.with_suffix(".native.pml")
    pml.write_text(
        tpl.substitute(
            PDB=Path(original_pdb).resolve().as_posix(),
            EXCL=",".join(exclude_resns or []),
            W=str(w),
            H=str(h),
            OUT_front=(str(outprefix) + "_front.png").replace("\\", "/"),
            OUT_side=(str(outprefix) + "_side.png").replace("\\", "/"),
            OUT_top=(str(outprefix) + "_top.png").replace("\\", "/"),
        ),
        encoding="utf-8",
    )
    return pml

def _control_and_rdk_paths(job: dict) -> tuple[str, str]:
    ctrl_path, rdk_path = "", ""
    for path_value, obj_name, _color in job.get("ligand_paths_and_colors", []):
        lowered = str(obj_name or "").lower()
        if "control" in lowered and not ctrl_path:
            ctrl_path = path_value
        elif "rdk" in lowered and not rdk_path:
            rdk_path = path_value
    if not ctrl_path and not rdk_path and job.get("ligand_paths_and_colors"):
        rdk_path = job["ligand_paths_and_colors"][0][0]
    return ctrl_path, rdk_path


def _job_to_pml(job: dict, tmpdir: Path) -> Optional[Path]:
    kind = job.get("kind")
    if kind == "three":
        ctrl_path, rdk_path = _control_and_rdk_paths(job)
        return write_multiview_pml(
            receptor_path=job["receptor_path"],
            control_path=ctrl_path,
            rdk_path=rdk_path,
            out_pml=tmpdir / (Path(job["outprefix"]).name + ".pml"),
            label_top_n_res=int(job.get("label_top_n_res", 5)),
            label_cutoff=float(job.get("label_cutoff", 5.0)),
            outprefix=str(job.get("outprefix", "")),
            viewport=tuple(job.get("viewport", (192, 144))),
        )
    if kind == "native":
        return write_native_pml(
            original_pdb=job["original_pdb"],
            outprefix=Path(job["outprefix"]),
            exclude_resns=job.get("exclude_resns", []),
            viewport=tuple(job.get("viewport", (192, 144))),
        )
    return None

def _run_job_via_cli(job: dict) -> bool:
    with tempfile.TemporaryDirectory() as td:
        pml = _job_to_pml(job, Path(td))
        if pml is None:
            raise RuntimeError("unsupported job")
        return render_pml_headless(pml) == 0


def _run_one_job_mp(job: dict, mode: str = "cli") -> tuple:
    """
    Returns: (ok:bool, t_ms:int, err:str|None)
    """
    t0 = time.time()
    try:
        if mode == "cli":
            ok = _run_job_via_cli(job)
        else:
            _dispatch_render_job(job)
            ok = True
        return ok, int((time.time() - t0) * 1000), None
    except Exception as exc:
        return False, int((time.time() - t0) * 1000), str(exc)

def launch_pymol_with_pml(pml_path: Path, pymol_exe: Optional[str] = None) -> None:
    """
    Launch PyMOL GUI with the given .pml (non-blocking).
    If pymol_exe is None, tries 'pymol' on PATH.
    """
    exe = pymol_exe or "pymol"
    try:
        subprocess.Popen([exe, str(pml_path)], shell=False)
    except FileNotFoundError:
        print(
            f"[capture_pose] Could not find PyMOL executable '{exe}'. Open manually: {pml_path}"
        )

def render_pml_headless(pml_path: Path, pymol_exe: Optional[str] = None) -> int:
    """
    Render a .pml without pymol2 using PyMOL CLI.
    The .pml must include PNG commands; returns process returncode.
    """
    exe = pymol_exe or "pymol"
    try:
        import shlex
        import subprocess

        pml_abs = Path(pml_path).resolve()
        cmd = [exe, "-cq", str(pml_abs)]
        print("[render-cli] " + " ".join(shlex.quote(x) for x in cmd))
        return subprocess.run(cmd, check=False).returncode
    except FileNotFoundError:
        print(f"[capture_pose] PyMOL CLI not found: '{exe}'")
        return 127
