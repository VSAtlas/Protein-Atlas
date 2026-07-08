"""PyMOL-backed pose rendering helpers."""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from docking import capture_pose_common as common
from docking.capture_pose_queue import _enqueue_job

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
            if (
                stem not in ctrl_set
                and isinstance(sc, (int, float))
                and raw_docked.get(lig)
            ):
                d = abs(sc - best_ctrl_score)
                if d < metric:
                    nearest_rdk, metric = lig, d
    else:
        for lig, rec in results_for_stage.items():
            stem = Path(lig).stem.split("_stage")[0].lower()
            sc = rec.get("score")
            if (
                stem not in ctrl_set
                and isinstance(sc, (int, float))
                and raw_docked.get(lig)
            ):
                if sc < metric:
                    nearest_rdk, metric = lig, sc

    return best_ctrl, nearest_rdk

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
    label_top_n_res: int | None = None,
    label_cutoff: float | None = None,
    viewport: tuple[int, int] | None = None,
    hide_receptor: bool = False,
) -> list[str]:
    """
    Saves <outprefix>_[front|side|top].png.
    Receptor shown as transparent surface; ligands as sticks; labels top-N closest residues (CA) within cutoff Å.
    If hide_receptor is True, receptor is not drawn (pair-only views).
    """
    if label_top_n_res is None:
        label_top_n_res = common._cfg_int("LABEL_TOP_N_RES", default=5)
    if label_cutoff is None:
        label_cutoff = common._cfg_float("LABEL_CUTOFF_A", 5.0)
    if viewport is None:
        viewport = (
            common._cfg_int("VIEWPORT_W", default=640),
            common._cfg_int("VIEWPORT_H", default=480),
        )

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
        fallback_candidates=[lp for lp, _, _ in ligand_paths_and_colors]
        + [obj for _, obj, _ in ligand_paths_and_colors],
        fallback_name=f"{pdb_id}_render" if pdb_id else "render",
    )

    if resolved_receptor is None or not str(resolved_receptor):
        print("[capture_pose] Receptor path could not be resolved; skipping renders.")
        return []
    if resolved_outprefix is None:
        print("[capture_pose] Output prefix could not be resolved; skipping renders.")
        return []

    receptor_path = str(resolved_receptor)
    outprefix = str(resolved_outprefix)

    if common._DEFER_MODE:
        _enqueue_job(
            "three",
            {
                "receptor_path": receptor_path,
                "ligand_paths_and_colors": [
                    (str(p), str(n), str(c)) for (p, n, c) in ligand_paths_and_colors
                ],
                "outprefix": outprefix,
                "cfg": cfg_data if cfg is not None else None,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_token": ph_token,
                "tag": tag,
                "label_top_n_res": int(label_top_n_res),
                "label_cutoff": float(label_cutoff),
                "viewport": list(viewport),
                "hide_receptor": bool(hide_receptor),
            },
        )
        return []

    PyMOL = _with_pymol()
    if PyMOL is None:
        return []

    if not Path(receptor_path).is_file():
        print(f"[capture_pose] Receptor missing: {receptor_path}")
        return []

    with PyMOL() as pm:
        cmd = pm.cmd
        cmd.reinitialize()

        cmd.load(receptor_path, "receptor")
        cmd.hide("everything", "receptor")
        if not hide_receptor:
            cmd.show("surface", "receptor")
            cmd.set_color("gray90", [230, 230, 230])
            cmd.color("gray90", "receptor")
            cmd.set("transparency", 0.30, "receptor")

        lig_objects: List[str] = []
        for lig_path, obj_name, color in ligand_paths_and_colors:
            if not lig_path or not Path(lig_path).is_file():
                continue
            lig_path_eff = _prefer_best_pdb(lig_path)
            if not Path(lig_path_eff).is_file():
                print(
                    f"[capture_pose] missing ligand file: {lig_path_eff} (original: {lig_path})"
                )
                continue

            cmd.load(lig_path_eff, obj_name)
            cmd.show("sticks", obj_name)
            cmd.color(color, obj_name)
            lig_objects.append(obj_name)

        lig_union = " or ".join(lig_objects) if lig_objects else "receptor"

        if lig_objects:
            cmd.select(
                "active_site_all", f"receptor within {label_cutoff} of ({lig_union})"
            )

            lig_model = cmd.get_model(lig_union)
            lig_coords = [(a.coord[0], a.coord[1], a.coord[2]) for a in lig_model.atom]

            def _min_dist_to_lig(x, y, z, coords):
                if not coords:
                    return float("inf")
                dx = x - coords[0][0]
                dy = y - coords[0][1]
                dz = z - coords[0][2]
                best = (dx * dx + dy * dy + dz * dz) ** 0.5
                for lx, ly, lz in coords[1:]:
                    dx = x - lx
                    dy = y - ly
                    dz = z - lz
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

            top_res = (
                sorted(distances, key=lambda t: t[4])[:label_top_n_res]
                if distances
                else []
            )
            if top_res:
                top_sel = " or ".join(
                    [
                        f"receptor and chain {c} and resi {r}"
                        for _, c, r, _, _ in top_res
                    ]
                )
                cmd.select("top_site", top_sel)
                cmd.show("sticks", "top_site")
                cmd.color("gray", "top_site")
                cmd.label(
                    "top_site and name CA and (alt '' or alt A)", "resn + '-' + resi"
                )

        focus_sel = lig_union if lig_objects else "receptor"
        cmd.zoom(focus_sel, 10)
        cmd.viewport(*viewport)
        cmd.set("antialias", 2)
        cmd.set("ray_opaque_background", 0)

        cmd.sync()
        cmd.refresh()
        cmd.png(f"{outprefix}_front.png", ray=0)

        cmd.sync()
        cmd.refresh()
        cmd.turn("y", 90)
        cmd.png(f"{outprefix}_side.png", ray=0)

        cmd.sync()
        cmd.refresh()
        cmd.turn("x", 90)
        cmd.png(f"{outprefix}_top.png", ray=0)

        return [
            f"{outprefix}_front.png",
            f"{outprefix}_side.png",
            f"{outprefix}_top.png",
        ]

def _render_native_on_original_pdb(
    original_pdb: str | None,
    outprefix: str | Path | None,
    *,
    cfg: Optional[dict] = None,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    tag: str = "renders",
    exclude_resns: Sequence[str] = tuple(common.EXCLUDE_HET_IDS),
    viewport: Tuple[int, int] = (192, 144),
) -> None:
    """
    Render original PDB with native ligand(s): polymer surface + HETATM (minus excludes) as sticks.
    Saves <outprefix>_[front|side|top].png.
    """
    cfg_data = common._cfg_data_or_default(cfg)
    paths_obj = None
    outdir: Optional[Path] = None
    if cfg_data and pdb_id:
        try:
            paths_obj, _default_receptor, outdir = common._resolve_receptor_and_outprefix(
                cfg_data,
                pdb_id,
                variant=variant,
                ph_token=ph_token,
                tag=tag,
                receptor_kind="cleaned",
            )
        except Exception:
            paths_obj = None
            outdir = None

    resolved_outprefix = common._resolve_output_prefix(
        outprefix,
        outdir,
        fallback_candidates=[original_pdb],
        fallback_name=f"{pdb_id}_native" if pdb_id else "native",
    )
    resolved_original = common._resolve_original_path(original_pdb, paths_obj)

    if resolved_outprefix is None:
        print("[capture_pose] Output prefix could not be resolved; skipping native render.")
        return

    outprefix_value = str(resolved_outprefix)

    if common._DEFER_MODE:
        _enqueue_job(
            "native",
            {
                "original_pdb": str(resolved_original) if resolved_original else str(original_pdb or ""),
                "outprefix": outprefix_value,
                "cfg": cfg_data if cfg is not None else None,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_token": ph_token,
                "tag": tag,
                "exclude_resns": list(exclude_resns or []),
                "viewport": list(viewport),
            },
        )
        return

    PyMOL = _with_pymol()
    if PyMOL is None:
        return

    if not resolved_original or not resolved_original.is_file():
        print(f"[capture_pose] Original PDB missing: {original_pdb}")
        return

    with PyMOL() as pm:
        cmd = pm.cmd
        cmd.reinitialize()

        cmd.load(str(resolved_original), "orig")
        cmd.hide("everything")
        cmd.show("surface", "orig and polymer")
        cmd.set_color("gray90", [230, 230, 230])
        cmd.color("gray90", "orig and polymer")
        cmd.set("transparency", 0.30, "orig and polymer")

        excluded = "+".join(sorted(set(exclude_resns or [])))
        cmd.select(
            "native_lig",
            f"(hetatm and not polymer and not solvent) and not resn {excluded}",
        )
        if cmd.count_atoms("native_lig") > 0:
            cmd.show("sticks", "native_lig")
            cmd.color("green", "native_lig")

            cmd.select("near_native", "orig within 5 of native_lig and polymer.protein")
            cmd.show("sticks", "near_native")
            cmd.color("gray", "near_native")
            cmd.label("near_native and name CA", "resn + '-' + resi")

            focus_sel = "native_lig or near_native"
        else:
            focus_sel = "orig and polymer"

        cmd.zoom(focus_sel, 10)
        cmd.viewport(*viewport)
        cmd.set("antialias", 2)
        cmd.set("ray_opaque_background", 0)
        cmd.sync()
        cmd.refresh()
        cmd.png(f"{outprefix_value}_front.png", ray=0)
        cmd.sync()
        cmd.refresh()
        cmd.turn("y", 90)
        cmd.png(f"{outprefix_value}_side.png", ray=0)
        cmd.sync()
        cmd.refresh()
        cmd.turn("x", 90)
        cmd.png(f"{outprefix_value}_top.png", ray=0)

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

    try:
        from pymol import cmd
    except Exception as e:
        raise RuntimeError(
            "PyMOL is required to render poses. Install with: "
            "conda install -c conda-forge pymol-open-source"
        ) from e

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

    # Keep receptor hidden until the ligand-centered pocket selection is defined.
    cmd.hide("everything", rec_obj)
    cmd.set_color("gray90", [230, 230, 230])
    cmd.color("gray90", rec_obj)
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.55)
    try:
        cmd.set("surface_quality", 1)
        cmd.set("surface_smooth_edges", 1)
        cmd.set("mesh_width", 0.25)
    except Exception:
        pass

    protein_transparency = kwargs.get("protein_transparency")
    receptor_transparency = 0.65
    if protein_transparency is not None:
        try:
            receptor_transparency = float(protein_transparency)
        except Exception:
            pass

    # Ligand color (keep if provided; otherwise default elsewhere)
    lig_color = kwargs.get("ligand_color")
    if lig_color:
        cmd.color(str(lig_color), lig_obj)

    hide_ligand = bool(kwargs.get("hide_ligand") or kwargs.get("without_ligand"))
    if hide_ligand:
        cmd.hide("everything", lig_obj)
    else:
        # Sticks for ligand; protein already surface
        cmd.show("sticks", lig_obj)

    # Centroid + labels (unchanged)
    cmd.select("lig_heavy_tmp", f"({lig_obj}) and not elem H")
    coords: list[list[float]] = []
    cmd.iterate_state(
        1, "lig_heavy_tmp", "coords.append([x,y,z])", space={"coords": coords}
    )
    centroid_pos: tuple[float, float, float] | None = None
    if coords:
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        centroid_pos = (cx, cy, cz)
        cmd.pseudoatom("lig_centroid", pos=[cx, cy, cz])
        cmd.hide("everything", "lig_centroid")
    else:
        cmd.orient(lig_obj)

    cutoff = float(proximity_cutoff)
    cmd.select("active_site_all", f"byres ({rec_obj} within {cutoff} of {lig_obj})")
    try:
        if int(cmd.count_atoms("active_site_all")) <= 0:
            cmd.select("active_site_all", rec_obj)
    except Exception:
        cmd.select("active_site_all", rec_obj)
    pocket_color = str(kwargs.get("pocket_color") or "gray90")
    protein_style = str(kwargs.get("protein_style") or "cartoon").strip().lower()
    carve_cutoff = float(kwargs.get("surface_carve_cutoff", 5.0))
    carve_normal_cutoff = float(kwargs.get("surface_carve_normal_cutoff", -0.1))
    carved_surface = False
    if protein_style in {"surface", "both"}:
        try:
            cmd.set("surface_carve_selection", "lig_heavy_tmp", rec_obj)
            cmd.set("surface_carve_cutoff", carve_cutoff, rec_obj)
            cmd.set("surface_carve_normal_cutoff", carve_normal_cutoff, rec_obj)
            cmd.show("surface", rec_obj)
            cmd.color(pocket_color, rec_obj)
            cmd.set("transparency", receptor_transparency, rec_obj)
            carved_surface = True
        except Exception:
            carved_surface = False
    if protein_style == "cartoon":
        cmd.show("cartoon", rec_obj)
        cmd.color(pocket_color, rec_obj)
        try:
            cmd.set("cartoon_transparency", 0.20, rec_obj)
        except Exception:
            pass
    elif protein_style == "lines":
        cmd.show("lines", "active_site_all")
        cmd.color(pocket_color, "active_site_all")
    elif protein_style == "mesh":
        cmd.show("mesh", "active_site_all")
        cmd.color(pocket_color, "active_site_all")
    elif protein_style == "both":
        cmd.show("mesh", "active_site_all")
        cmd.color(pocket_color, "active_site_all")
    elif protein_style != "none" and not carved_surface:
        cmd.show("surface", "active_site_all")
        cmd.color(pocket_color, "active_site_all")
        cmd.set("transparency", receptor_transparency, "active_site_all")

    cmd.select("near_CA", "active_site_all and name CA")
    label_selection = "near_CA"
    try:
        label_count = int(top_n_residues)
        if label_count <= 0:
            label_selection = "none"
        elif centroid_pos is not None:
            cx, cy, cz = centroid_pos
            ca_atoms = sorted(
                cmd.get_model("near_CA").atom,
                key=lambda atom: (
                    (float(atom.coord[0]) - cx) ** 2
                    + (float(atom.coord[1]) - cy) ** 2
                    + (float(atom.coord[2]) - cz) ** 2
                ),
            )
            cmd.select("label_CA", "none")
            for atom in ca_atoms[:label_count]:
                cmd.select("label_CA", f"label_CA or ({rec_obj} and index {int(atom.index)})")
            label_selection = "label_CA"
    except Exception:
        label_selection = "near_CA"
    cmd.label(f"{label_selection} and name CA", '"%s-%s" % (resn, resi)')

    context_mode = str(kwargs.get("view_context") or "full").strip().lower()
    if context_mode == "pocket":
        cmd.select("view_context", f"active_site_all or {lig_obj}")
    else:
        cmd.select("view_context", f"{rec_obj} or {lig_obj}")
    zoom_buffer = float(kwargs.get("zoom_buffer", 4.0))
    cmd.orient("view_context")
    cmd.zoom("view_context", buffer=zoom_buffer, complete=1)
    clip_slab = float(kwargs.get("clip_slab", 28.0))
    if clip_slab > 0:
        try:
            cmd.clip("slab", clip_slab, "view_context")
        except Exception:
            pass

    base = outprefix
    if os.path.isdir(outprefix) or outprefix.endswith(os.sep):
        base = os.path.join(outprefix, "top_pose")
    os.makedirs(os.path.dirname(base) or ".", exist_ok=True)
    ray = int(bool(kwargs.get("ray", False)))

    cmd.sync()
    cmd.refresh()
    cmd.png(base + "_front.png", ray=ray, width=w, height=h)
    cmd.sync()
    cmd.refresh()
    cmd.turn("y", 90)
    cmd.png(base + "_side.png", ray=ray, width=w, height=h)
    cmd.sync()
    cmd.refresh()
    cmd.turn("x", 90)
    cmd.png(base + "_top.png", ray=ray, width=w, height=h)

    try:
        cmd.delete("lig_heavy_tmp")
        cmd.delete("lig_centroid")
        cmd.delete("active_site_all")
        cmd.delete("view_context")
        cmd.delete("label_CA")
    except Exception:
        pass
