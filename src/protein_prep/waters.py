"""Water and centroid helper utilities used by protein prep."""

from __future__ import annotations

from pathlib import Path
from typing import Union

from protein_prep.os_utils import _as_path
from protein_prep.pdb_records import line_xyz

def compute_control_centroids(
    ligands_dir: Union[str, Path],
) -> list[tuple[float, float, float]]:
    """Return one centroid per *.pdb control ligand under ligands_dir."""
    ligands_path = _as_path(ligands_dir)
    pts: list[tuple[float, float, float]] = []
    if not ligands_path.exists():
        return pts
    for p in sorted(ligands_path.glob("*.pdb")):
        try:
            n = 0
            sx = sy = sz = 0.0
            with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if not ln.startswith(("ATOM  ", "HETATM")):
                        continue
                    # drop hydrogens by element column (77–78)
                    if len(ln) >= 78 and ln[76:78].strip().upper() == "H":
                        continue
                    xyz = line_xyz(ln)
                    if xyz is None:
                        continue
                    x, y, z = xyz
                    sx += x
                    sy += y
                    sz += z
                    n += 1
            if n > 0:
                pts.append((sx / n, sy / n, sz / n))
        except Exception:
            # best-effort; ignore malformed ligand files
            pass
    return pts

def filter_waters_near_points(
    src_pdb: Union[str, Path],
    dst_pdb: Union[str, Path],
    points: list[tuple[float, float, float]],
    radius_A: float,
) -> int:
    """
    Copy src_pdb → dst_pdb keeping all non-water records and only water residues
    with any atom within radius_A of any reference point. Returns #water residues kept.

    Water 3-letter residue names are loaded from YAML if available, else a fallback set.
    Optional: set env WATER_NAMES_YAML to a file path. YAML may be a list (["HOH","WAT",...])
    or a dict containing any of these keys: water_resnames, water, waters, solvent_water, water_aliases.
    """
    import os

    R2 = float(radius_A) * float(radius_A)

    # --- Load water residue names ---
    water3: set[str] = {
        # Common crystallographic water names
        "HOH",
        "WAT",
        "OH2",
        "DOD",
        "D2O",
        # MD water models that sometimes leak into PDB-like files
        "TIP",
        "TP3",
        "TP4",
        "TP5",
        "TIP3",
        "TIP4",
        "TIP5",
        "SOL",
        "W",
        # protonation variants
        "H3O",
        "OH-",
        "OHO",
        "H2O",
    }
    try:
        import yaml  # type: ignore[import-untyped]  # PyYAML

        # prefer explicit env; else try a conventional chemdb path beside this file
        default_yaml = os.path.join(
            os.path.dirname(__file__), "chemdb", "water_names.yaml"
        )
        yaml_path = os.environ.get("WATER_NAMES_YAML", default_yaml)
        if os.path.exists(yaml_path):
            with open(yaml_path, "r", encoding="utf-8") as yf:
                y = yaml.safe_load(yf)
            candidates: list[str] = []
            if isinstance(y, dict):
                for k in (
                    "water_resnames",
                    "water",
                    "waters",
                    "solvent_water",
                    "water_aliases",
                ):
                    v = y.get(k, [])
                    if isinstance(v, str):
                        candidates.append(v)
                    elif isinstance(v, (list, tuple, set)):
                        candidates.extend(v)
            elif isinstance(y, (list, tuple, set)):
                candidates = list(y)
            # Normalize to 3-char PDB residue names
            for itm in candidates:
                try:
                    water3.add(str(itm).strip().upper()[:3])
                except Exception:
                    pass
    except Exception:
        # YAML absent or parse error → fall back silently
        pass

    keep: set[tuple[str, str]] = set()  # (chain, resseq+icode)

    # --- First pass: decide which water residues to keep ---
    try:
        with open(src_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                resn = ln[17:20].strip().upper()
                if resn not in water3:
                    continue
                xyz = line_xyz(ln)
                if xyz is None:
                    continue
                x, y, z = xyz
                for cx, cy, cz in points:
                    dx = x - cx
                    dy = y - cy
                    dz = z - cz
                    if (dx * dx + dy * dy + dz * dz) <= R2:
                        keep.add((ln[21], ln[22:27]))  # chain, resseq+icode
                        break
    except Exception:
        # If anything goes wrong during reading, just fall back to keeping none.
        pass

    # --- Second pass: write out everything but only keep selected waters ---
    with (
        open(src_pdb, "r", encoding="utf-8", errors="ignore") as fh,
        open(dst_pdb, "w", encoding="utf-8") as out,
    ):
        for ln in fh:
            if ln.startswith(("ATOM  ", "HETATM")):
                resn = ln[17:20].strip().upper()
                if resn in water3:
                    key = (ln[21], ln[22:27])
                    if key in keep:
                        out.write(ln)
                    # else drop water line
                else:
                    out.write(ln)
            else:
                out.write(ln)

    return len(keep)
