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
    receptor_path: str,
    ligand_paths_and_colors: List[Tuple[str, str, str]],  # [(path, object_name, color), ...]
    outprefix: str,
    label_top_n_res: int = 5,
    label_cutoff: float = 5.0,
    viewport: Tuple[int, int] = (1200, 900),
    hide_receptor: bool = False,
) -> None:
    """
    Saves <outprefix>_[front|side|top].png.
    Receptor shown as transparent surface; ligands as sticks; labels top-N closest residues (CA) within cutoff Å.
    If hide_receptor is True, receptor is not drawn (pair-only views).
    """
    PyMOL = _with_pymol()
    if PyMOL is None:
        return

    if not Path(receptor_path).is_file():
        print(f"[capture_pose] Receptor missing: {receptor_path}")
        return

    with PyMOL() as pm:
        cmd = pm.cmd
        cmd.reinitialize()

        # Receptor
        cmd.load(receptor_path, "receptor")
        cmd.hide("everything", "receptor")
        if not hide_receptor:
            cmd.show("surface", "receptor")
            cmd.set("transparency", 0.30, "receptor")
            cmd.color("slate", "receptor")

        # Ligands
        lig_objects: List[str] = []
        for lig_path, obj_name, color in ligand_paths_and_colors:
            if not lig_path or not Path(lig_path).is_file():
                continue
            cmd.load(lig_path, obj_name)
            cmd.show("sticks", obj_name)
            cmd.color(color, obj_name)
            lig_objects.append(obj_name)

        lig_union = " or ".join(lig_objects) if lig_objects else "receptor"

        # Label nearest residues (rank by CA distance)
        if lig_objects:
            cmd.select("active_site_all", f"receptor within {label_cutoff} of ({lig_union})")

            # Precompute ligand atom coordinates once to avoid get_distance on multi-atom selections
            lig_model = cmd.get_model(lig_union)
            lig_coords = [(a.coord[0], a.coord[1], a.coord[2]) for a in lig_model.atom]

            def _min_dist_to_lig(x, y, z, coords):
                if not coords:
                    return float("inf")
                # Euclidean distance to nearest ligand atom
                dx = x - coords[0][0]
                dy = y - coords[0][1]
                dz = z - coords[0][2]
                best = (dx * dx + dy * dy + dz * dz) ** 0.5
                for (lx, ly, lz) in coords[1:]:
                    dx = x - lx;
                    dy = y - ly;
                    dz = z - lz
                    d = (dx * dx + dy * dy + dz * dz) ** 0.5
                    if d < best:
                        best = d
                return best

            distances = []
            sel_ca = "active_site_all and name CA and (alt '' or alt A)"

            # Pull CA atoms (with coords) directly and compute min distance to ligand atoms
            ca_model = cmd.get_model(sel_ca)
            for a in ca_model.atom:
                d = _min_dist_to_lig(a.coord[0], a.coord[1], a.coord[2], lig_coords)
                distances.append((a.model, a.chain, a.resi, a.resn, d))

            top_res = sorted(distances, key=lambda t: t[4])[:label_top_n_res] if distances else []
            if top_res:
                top_sel = " or ".join([f"receptor and chain {c} and resi {r}" for _, c, r, _, _ in top_res])
                cmd.select("top_site", top_sel)
                cmd.show("sticks", "top_site")
                cmd.color("cyan", "top_site")
                cmd.label("top_site and name CA and (alt '' or alt A)", "resn + resi")

        # Views
        cmd.zoom(lig_union, 10)
        cmd.viewport(*viewport)
        cmd.set("antialias", 2)
        cmd.set("ray_opaque_background", 0)

        cmd.png(f"{outprefix}_front.png", ray=1)
        cmd.turn("y", 90)
        cmd.png(f"{outprefix}_side.png", ray=1)
        cmd.turn("x", 90)
        cmd.png(f"{outprefix}_top.png", ray=1)


def _render_native_on_original_pdb(
    original_pdb: str,
    outprefix: str,
    exclude_resns: Sequence[str] = tuple(EXCLUDE_HET_IDS),
    viewport: Tuple[int, int] = (1200, 900),
) -> None:
    """
    Render original PDB with native ligand(s): polymer surface + HETATM (minus excludes) as sticks.
    Saves <outprefix>_[front|side|top].png.
    """
    PyMOL = _with_pymol()
    if PyMOL is None:
        return

    if not Path(original_pdb).is_file():
        print(f"[capture_pose] Original PDB missing: {original_pdb}")
        return

    with PyMOL() as pm:
        cmd = pm.cmd
        cmd.reinitialize()

        cmd.load(original_pdb, "orig")
        cmd.hide("everything")

        # Polymer & coloring
        cmd.show("surface", "orig and polymer")
        cmd.set("transparency", 0.30, "orig and polymer")
        cmd.color("slate", "orig and polymer")
        try:
            cmd.util.cbag("orig and polymer")
        except Exception:
            pass

        # Native ligands (exclude common junk)
        excl = "+".join(sorted(set(exclude_resns or [])))
        cmd.select("native_lig", f"(hetatm and not polymer and not solvent) and not resn {excl}")
        if cmd.count_atoms("native_lig") > 0:
            cmd.show("sticks", "native_lig")
            cmd.color("green", "native_lig")

            # Nearby residues (sticks on residues; labels only on CA)
            cmd.select("near_native", "orig within 5 of native_lig and polymer.protein")
            cmd.show("sticks", "near_native")
            cmd.color("cyan", "near_native")
            cmd.label("near_native and name CA", "resn + resi")

            focus_sel = "native_lig or near_native"
        else:
            focus_sel = "orig and polymer"

        # Views
        cmd.zoom(focus_sel, 10)
        cmd.viewport(*viewport)
        cmd.set("antialias", 2)
        cmd.set("ray_opaque_background", 0)

        cmd.png(f"{outprefix}_front.png", ray=1)
        cmd.turn("y", 90)
        cmd.png(f"{outprefix}_side.png", ray=1)
        cmd.turn("x", 90)
        cmd.png(f"{outprefix}_top.png", ray=1)


# ============================================================
# PML writer + launch/fallback
# ============================================================

def write_multiview_pml(
    receptor_path: str,
    control_path: str,
    rdk_path: str,
    out_pml: Path,
    label_top_n_res: int = 5,
    label_cutoff: float = 5.0,
) -> Path:
    """
    Writes a PyMOL .pml that:
      - loads receptor (transparent surface)
      - optionally loads CONTROL (green sticks) and RDK (magenta sticks)
      - labels top-N nearby residues within cutoff Å
      - stores three scenes: front, side, top
    """
    out_pml.parent.mkdir(parents=True, exist_ok=True)

    control_path = _prefer_best_pdb(control_path) if control_path else ""
    rdk_path = _prefer_best_pdb(rdk_path) if rdk_path else ""

    # Use forward slashes so PyMOL on Windows is happy
    def as_posix_or_empty(path: str) -> str:
        return (Path(path).resolve().as_posix() if path else "")

    receptor_posix = as_posix_or_empty(receptor_path)
    control_posix = as_posix_or_empty(control_path)
    rdk_posix = as_posix_or_empty(rdk_path)

    # Build the ligand union expression used for zoom/labels
    lig_parts = []
    if control_posix:
        lig_parts.append("control")
    if rdk_posix:
        lig_parts.append("rdk")
    lig_union = " or ".join(lig_parts) if lig_parts else "receptor"

    tpl = Template(r"""
reinitialize
bg_color white
set ray_opaque_background, 0
set antialias, 2

load "$RECEPTOR", receptor
hide everything, receptor
show surface, receptor
set transparency, 0.30, receptor
color slate, receptor
python
try:
    cmd.util.cbag("receptor and polymer")
except Exception:
    pass
python end

# Ligands (conditionally present)
$LOAD_CONTROL
$LOAD_RDK

# Union for focus
select lig_union, ($LIG_UNION)

# Label nearby residues (build a ranked list in a small Python block)
select near_res, (receptor within $CUTOFF of lig_union) and polymer.protein and name CA and (alt '' or alt A)
python
TOPN = $TOPN
import math

# Cache ligand atom coordinates once (avoid get_distance with multi-atom selections)
lig_model = cmd.get_model("lig_union")
lig_coords = [(a.coord[0], a.coord[1], a.coord[2]) for a in lig_model.atom]

def min_dist_to_lig(x, y, z, coords):
    if not coords:
        return float("inf")
    best = float("inf")
    for (lx, ly, lz) in coords:
        dx = x - lx; dy = y - ly; dz = z - lz
        d = math.sqrt(dx*dx + dy*dy + dz*dz)
        if d < best:
            best = d
    return best

dlist = []
# Collect distances to CA atoms near lig_union (alt '' or A only per selection above)
for (model_id, atom_index) in cmd.index("near_res"):
    m = cmd.get_model("near_res and index %d" % atom_index)
    if not m.atom:
        continue
    a = m.atom[0]
    dist = min_dist_to_lig(a.coord[0], a.coord[1], a.coord[2], lig_coords)
    dlist.append((a.model, a.chain, a.resi, dist))

# Sort by distance and keep unique residues up to TOPN
dlist.sort(key=lambda t: t[3])
keep = []
seen = set()
for (m,c,i,_) in dlist:
    key = (m,c,i)
    if key in seen:
        continue
    seen.add(key)
    keep.append(key)
    if len(keep) >= TOPN:
        break

sel = " or ".join(["%s//%s/%s" % (m,c,i) for (m,c,i) in keep])
if sel:
    cmd.select("top_site", sel)
    cmd.show("sticks", "top_site")
    cmd.color("cyan", "top_site")
    cmd.label("top_site and name CA and (alt '' or alt A)", "resn + resi")
python end

# Views & scenes
orient lig_union
zoom lig_union, 10
scene front, store
turn y, 90
scene side, store
turn x, 90
scene top, store

set scene_buttons, on
""")

    load_control = ""
    if control_posix:
        load_control = '\n'.join([
            f'load "{control_posix}", control',
            'show sticks, control',
            'color green, control',
        ])

    load_rdk = ""
    if rdk_posix:
        load_rdk = '\n'.join([
            f'load "{rdk_posix}", rdk',
            'show sticks, rdk',
            'color magenta, rdk',
        ])

    pml_text = tpl.substitute(
        RECEPTOR=receptor_posix,
        LOAD_CONTROL=load_control,
        LOAD_RDK=load_rdk,
        LIG_UNION=lig_union,
        CUTOFF=f"{label_cutoff:.2f}",
        TOPN=str(label_top_n_res),
    )

    out_pml.write_text(pml_text, encoding="utf-8")
    return out_pml


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
        # -cq = console (no GUI) + quiet
        return subprocess.run([exe, "-cq", str(pml_path)], check=False).returncode
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
    viewport: Tuple[int, int] = (800, 600),
) -> None:
    """
    Original one-ligand script:
    - ligand orange sticks
    - active-site residues within cutoff Å
    - rank by CA distance; label top-N
    - receptor transparent surface (slate)
    - save front/side/top PNGs
    """
    PyMOL = _with_pymol()
    if PyMOL is None:
        return

    if not Path(receptor).is_file():
        print(f"[capture_pose] Receptor missing: {receptor}")
        return
    if not Path(ligand).is_file():
        print(f"[capture_pose] Ligand missing: {ligand}")
        return

    with PyMOL() as pymol:
        cmd = pymol.cmd

        cmd.load(receptor, "receptor")
        cmd.load(ligand, "ligand")

        cmd.hide("everything")
        cmd.color("orange", "ligand")
        cmd.show("sticks", "ligand")

        # Active site: residues within cutoff Å
        cmd.select("active_site_all", f"receptor within {proximity_cutoff} of ligand")

        # Rank residues by proximity to ligand CA atoms (no temp measurements)
        distances: List[Tuple[str, str, str, str, float]] = []
        cmd.iterate(
            "active_site_all and name CA",
            "distances.append((model, chain, resi, resn, cmd.get_distance('ligand', f'{model}//{chain}/{resi}/CA')))",
            space={"distances": distances, "cmd": cmd},
        )

        # Sort and select top N
        top_residues = sorted(distances, key=lambda x: x[4])[:top_n_residues] if distances else []
        top_sel = " or ".join([f"receptor and chain {c} and resi {r}" for _, c, r, _, _ in top_residues])

        if top_sel:
            cmd.select("top_site", top_sel)
            cmd.show("sticks", "top_site")
            cmd.color("cyan", "top_site")
            cmd.label("top_site and name CA", "resn + resi")

        # Full receptor as transparent surface
        cmd.show("surface", "receptor")
        cmd.set("transparency", 0.3, "receptor")
        cmd.color("slate", "receptor")

        # View and save
        cmd.zoom("ligand or top_site", 10)
        cmd.viewport(*viewport)
        cmd.set("ray_opaque_background", 0)

        cmd.png(f"{outprefix}_front.png", ray=1)
        cmd.turn("y", 90)
        cmd.png(f"{outprefix}_side.png", ray=1)
        cmd.turn("x", 90)
        cmd.png(f"{outprefix}_top.png", ray=1)


# ============================================================
# __main__
# ============================================================

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("Usage: python capture_pose.py <receptor.pdb|pdbqt> <ligand.pdb|pdbqt> <outprefix>")
        sys.exit(2)
    _cli_render_active_site(sys.argv[1], sys.argv[2], sys.argv[3])