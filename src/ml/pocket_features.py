from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


POCKET_FEATURE_COLUMNS = (
    "pocket_frac_polar",
    "pocket_frac_hydrophobic",
    "pocket_frac_aromatic",
    "pocket_frac_charged_pos",
    "pocket_frac_charged_neg",
    "pocket_net_charge_proxy",
    "pocket_rgyr",
    "pocket_bbox_vol",
    "pocket_pca_ratio1",
    "pocket_pca_ratio2",
    "pocket_metal_zn",
    "pocket_metal_fe",
    "pocket_metal_mg",
    "pocket_metal_ca",
    "pocket_metal_mn",
    "pocket_metal_cu",
    "pocket_metal_co",
    "pocket_metal_ni",
)

_POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
_HYDROPHOBIC = {"ALA", "VAL", "ILE", "LEU", "MET", "PRO"}
_AROMATIC = {"PHE", "TRP", "TYR", "HIS"}
_POS = {"LYS", "ARG", "HIS"}
_NEG = {"ASP", "GLU"}
_METALS = {"ZN", "FE", "MG", "CA", "MN", "CU", "CO", "NI"}


def load_pdb_atoms(pdb_path: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if not Path(pdb_path).exists():
        return pd.DataFrame(
            columns=["resname", "chain", "resid", "atom_name", "x", "y", "z", "element"]
        )
    with Path(pdb_path).open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            resname = (line[17:20] or "").strip().upper()
            chain = (line[21:22] or "").strip()
            resid = (line[22:27] or "").strip()
            atom_name = (line[12:16] or "").strip().upper()
            element = (line[76:78] or "").strip().upper()
            if not element and atom_name:
                element = "".join(ch for ch in atom_name if ch.isalpha())[:2].upper()
            rows.append(
                {
                    "resname": resname,
                    "chain": chain,
                    "resid": resid,
                    "atom_name": atom_name,
                    "x": x,
                    "y": y,
                    "z": z,
                    "element": element,
                }
            )
    return pd.DataFrame(rows)


def pocket_residue_composition(atom_table: pd.DataFrame) -> dict[str, float]:
    if atom_table.empty:
        return {
            "pocket_frac_polar": 0.0,
            "pocket_frac_hydrophobic": 0.0,
            "pocket_frac_aromatic": 0.0,
            "pocket_frac_charged_pos": 0.0,
            "pocket_frac_charged_neg": 0.0,
            "pocket_net_charge_proxy": 0.0,
        }
    residues = (
        atom_table[["resname", "chain", "resid"]]
        .drop_duplicates()
        .assign(resname=lambda x: x["resname"].astype(str).str.upper())
    )
    residues = residues.loc[
        ~residues["resname"].isin(_METALS)
        & residues["resname"].str.match(r"^[A-Z]{3}$", na=False)
    ].copy()
    if residues.empty:
        return {
            "pocket_frac_polar": 0.0,
            "pocket_frac_hydrophobic": 0.0,
            "pocket_frac_aromatic": 0.0,
            "pocket_frac_charged_pos": 0.0,
            "pocket_frac_charged_neg": 0.0,
            "pocket_net_charge_proxy": 0.0,
        }
    n = float(max(1, len(residues)))
    resnames = residues["resname"].tolist()
    n_polar = 0
    n_hydrophobic = 0
    n_aromatic = 0
    n_pos = 0
    n_neg = 0
    for resname in resnames:
        if resname in _POS:
            n_pos += 1
        elif resname in _NEG:
            n_neg += 1
        elif resname in _AROMATIC:
            n_aromatic += 1
        elif resname in _POLAR:
            n_polar += 1
        elif resname in _HYDROPHOBIC:
            n_hydrophobic += 1
    return {
        "pocket_frac_polar": float(n_polar / n),
        "pocket_frac_hydrophobic": float(n_hydrophobic / n),
        "pocket_frac_aromatic": float(n_aromatic / n),
        "pocket_frac_charged_pos": float(n_pos / n),
        "pocket_frac_charged_neg": float(n_neg / n),
        "pocket_net_charge_proxy": float(n_pos - n_neg),
    }


def pocket_geometry_features(atom_table: pd.DataFrame) -> dict[str, float]:
    if atom_table.empty:
        return {
            "pocket_rgyr": 0.0,
            "pocket_bbox_vol": 0.0,
            "pocket_pca_ratio1": 0.0,
            "pocket_pca_ratio2": 0.0,
        }
    coords = atom_table[["x", "y", "z"]].to_numpy(dtype=float)
    center = np.mean(coords, axis=0)
    centered = coords - center
    rgyr = float(np.sqrt(np.mean(np.sum(centered**2, axis=1))))
    mins = np.min(coords, axis=0)
    maxs = np.max(coords, axis=0)
    span = np.maximum(maxs - mins, 1e-6)
    bbox_vol = float(span[0] * span[1] * span[2])

    cov = np.cov(centered.T)
    evals = np.linalg.eigvalsh(cov)
    evals = np.sort(np.clip(np.asarray(evals, dtype=float), 1e-12, None))[::-1]
    total = float(np.sum(evals))
    if total <= 0:
        ratio1 = 0.0
        ratio2 = 0.0
    else:
        ratio1 = float(evals[0] / total)
        ratio2 = float(evals[1] / total) if evals.size > 1 else 0.0

    return {
        "pocket_rgyr": rgyr,
        "pocket_bbox_vol": bbox_vol,
        "pocket_pca_ratio1": ratio1,
        "pocket_pca_ratio2": ratio2,
    }


def metal_context_features(atom_table: pd.DataFrame) -> dict[str, float]:
    elements = atom_table.get("element", pd.Series(dtype=object)).astype(str).str.upper()
    values: dict[str, float] = {}
    for metal in sorted(_METALS):
        values[f"pocket_metal_{metal.lower()}"] = float((elements == metal).sum())
    return values


def _cache_key_for_path(path: Path) -> str:
    payload = str(Path(path).resolve()).encode("utf-8")
    return hashlib.sha1(payload).hexdigest()  # noqa: S324 - non-crypto cache key


def compute_pocket_features(
    pocket_pdb_path: Path,
    *,
    cache_root: Path | None = None,
) -> dict[str, float]:
    cache_dir = Path(cache_root or Path(__file__).resolve().parent / "cache" / "pocket_features")
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = _cache_key_for_path(pocket_pdb_path)
    cache_path = cache_dir / f"{key}.json"
    if cache_path.exists():
        try:
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            if isinstance(payload, dict):
                return {k: float(payload.get(k, 0.0)) for k in POCKET_FEATURE_COLUMNS}
        except Exception:
            pass

    atoms = load_pdb_atoms(pocket_pdb_path)
    features = {}
    features.update(pocket_residue_composition(atoms))
    features.update(pocket_geometry_features(atoms))
    features.update(metal_context_features(atoms))
    normalized = {k: float(features.get(k, 0.0)) for k in POCKET_FEATURE_COLUMNS}
    cache_path.write_text(json.dumps(normalized, indent=2), encoding="utf-8")
    return normalized
