"""Public pose-capture API composed from focused helpers."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from docking import capture_pose_common as common
from docking.capture_pose_queue import _enqueue_job, set_defer_mode
from docking.capture_pose_render import (
    _cli_render_active_site,
    _render_native_on_original_pdb,
    _render_three_views_with_pymol,
    _safe_open_csv_for_write,
    pick_control_and_nearest_rdk,
)
from docking.capture_pose_replay import (
    launch_pymol_with_pml,
    render_pml_headless,
    replay_deferred_jobs,
    replay_deferred_jobs_mp,
    write_multiview_pml,
    write_native_pml,
)

_COLOR_SCHEMES = {
    "okabe_ito": {
        "custom_defs": {
            "orangeOI": [230, 159, 0],
            "blueOI": [0, 114, 178],
            "tealOI": [0, 158, 115],
            "yellowOI": [240, 228, 66],
            "vermOI": [213, 94, 0],
            "purpleOI": [204, 121, 167],
        },
        "ligand_colors": {
            "reference": "blueOI",
            "candidate": "orangeOI",
        },
    }
}
_APPEARANCE_KEYS = (
    "ligand_color",
    "ligand_role",
    "color_scheme",
    "palette_defs",
    "protein_color",
    "protein_transparency",
    "pocket_color",
    "hide_ligand",
)


def _capture_outprefix(ligand_path: Any, out_path_or_dir: Any) -> str:
    try:
        if os.path.isdir(out_path_or_dir):
            return str(Path(out_path_or_dir) / Path(str(ligand_path)).stem)
    except Exception:
        pass
    return str(out_path_or_dir)


def _capture_scheme(kwargs: dict[str, Any]) -> dict[str, Any]:
    scheme_key = kwargs.get("color_scheme")
    if not scheme_key:
        return {}
    scheme = _COLOR_SCHEMES.get(str(scheme_key), {})
    return scheme if isinstance(scheme, dict) else {}


def _capture_ligand_color(kwargs: dict[str, Any]) -> str:
    if "ligand_color" in kwargs:
        return str(kwargs["ligand_color"])
    scheme = _capture_scheme(kwargs)
    ligand_colors = scheme.get("ligand_colors")
    if isinstance(ligand_colors, dict):
        role = str(kwargs.get("ligand_role") or "candidate").lower()
        return str(ligand_colors.get(role, ligand_colors.get("candidate", "magenta")))
    return "magenta"


def _capture_palette_defs(kwargs: dict[str, Any]) -> dict[str, Any]:
    palette_defs: dict[str, Any] = {}
    scheme = _capture_scheme(kwargs)
    custom_defs = scheme.get("custom_defs")
    if isinstance(custom_defs, dict):
        palette_defs.update(custom_defs)
    user_defs = kwargs.get("palette_defs")
    if isinstance(user_defs, dict):
        palette_defs.update(user_defs)
    return palette_defs


def _build_capture_payload(
    receptor_path: Any,
    ligand_path: Any,
    outprefix: str,
    kwargs: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    wants_custom = any(key in kwargs for key in _APPEARANCE_KEYS)
    ligand_color = _capture_ligand_color(kwargs) if wants_custom else "magenta"
    payload: dict[str, Any] = {
        "receptor_path": str(receptor_path),
        "ligand_paths_and_colors": [(str(ligand_path), "ligand", ligand_color)],
        "outprefix": outprefix,
        "label_top_n_res": int(kwargs.get("top_n_residues", 6)),
        "label_cutoff": float(kwargs.get("proximity_cutoff", 4.5)),
        "viewport": list(kwargs.get("viewport", (192, 144))),
        "hide_receptor": False,
    }
    appearance_kwargs: dict[str, Any] = {}
    if not wants_custom:
        return payload, appearance_kwargs

    palette_defs = _capture_palette_defs(kwargs)
    if palette_defs:
        payload["palette_defs"] = palette_defs
        appearance_kwargs["palette_defs"] = palette_defs
    for key in ("protein_color", "pocket_color", "protein_style", "view_context"):
        if key in kwargs:
            payload[key] = kwargs[key]
            appearance_kwargs[key] = kwargs[key]
    if "protein_transparency" in kwargs:
        transparency = float(kwargs["protein_transparency"])
        payload["protein_transparency"] = transparency
        appearance_kwargs["protein_transparency"] = transparency
    if "hide_ligand" in kwargs:
        hide_ligand = bool(kwargs["hide_ligand"])
        payload["hide_ligand"] = hide_ligand
        appearance_kwargs["hide_ligand"] = hide_ligand
    for key in ("surface_carve_cutoff", "surface_carve_normal_cutoff", "clip_slab", "ray"):
        if key in kwargs:
            payload[key] = kwargs[key]
    appearance_kwargs["ligand_color"] = ligand_color
    return payload, appearance_kwargs


def _render_capture_immediate(
    receptor_path: Any,
    ligand_path: Any,
    outprefix: str,
    kwargs: dict[str, Any],
    appearance_kwargs: dict[str, Any],
):
    base_kwargs = {
        "top_n_residues": kwargs.get("top_n_residues", 6),
        "proximity_cutoff": kwargs.get("proximity_cutoff", 4.5),
        "viewport": kwargs.get("viewport", (192, 144)),
        "zoom_buffer": kwargs.get("zoom_buffer", 4.0),
        "protein_style": kwargs.get("protein_style", "cartoon"),
        "surface_carve_cutoff": kwargs.get("surface_carve_cutoff", 5.0),
        "surface_carve_normal_cutoff": kwargs.get("surface_carve_normal_cutoff", -0.1),
        "clip_slab": kwargs.get("clip_slab", 0.0),
        "view_context": kwargs.get("view_context", "full"),
        "ray": kwargs.get("ray", False),
    }
    try:
        return _cli_render_active_site(
            str(receptor_path),
            str(ligand_path),
            outprefix,
            **base_kwargs,
            **appearance_kwargs,
        )
    except TypeError:
        return _cli_render_active_site(
            str(receptor_path),
            str(ligand_path),
            outprefix,
            **base_kwargs,
        )


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
      - protein_style: "cartoon", "lines", "mesh", "surface", "both", or "none"
      - view_context: "full" or "pocket"
      - hide_ligand: bool, use ligand for centering/selection but hide it in PNGs
      - surface_carve_cutoff: float, carved surface distance around the ligand
      - surface_carve_normal_cutoff: float, hide surface triangles not facing ligand
      - clip_slab: float, camera slab thickness around the pocket context
      - ray: bool, ray trace output PNGs
    """
    outprefix = _capture_outprefix(ligand_path, out_path_or_dir)
    payload, appearance_kwargs = _build_capture_payload(
        receptor_path,
        ligand_path,
        outprefix,
        kwargs,
    )
    if common._DEFER_MODE:
        _enqueue_job("three", payload)
        return None
    return _render_capture_immediate(
        receptor_path,
        ligand_path,
        outprefix,
        kwargs,
        appearance_kwargs,
    )


__all__ = [
    "_cli_render_active_site",
    "_render_native_on_original_pdb",
    "_render_three_views_with_pymol",
    "_safe_open_csv_for_write",
    "capture_pose",
    "launch_pymol_with_pml",
    "pick_control_and_nearest_rdk",
    "render_pml_headless",
    "replay_deferred_jobs",
    "replay_deferred_jobs_mp",
    "set_defer_mode",
    "write_multiview_pml",
    "write_native_pml",
]
