from __future__ import annotations

import json
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


WATER_RESNAMES = {"HOH", "WAT", "DOD", "SOL", "TIP", "TIP3"}
ION_RESNAMES = {
    "AG",
    "AL",
    "BA",
    "BR",
    "CA",
    "CD",
    "CL",
    "CO",
    "CS",
    "CU",
    "FE",
    "HG",
    "I",
    "K",
    "LI",
    "MG",
    "MN",
    "NA",
    "NI",
    "PB",
    "RB",
    "SR",
    "ZN",
}
COMMON_BUFFER_RESNAMES = {
    "ACT",
    "ACE",
    "ACN",
    "BME",
    "BOG",
    "DMS",
    "EDO",
    "EPE",
    "GOL",
    "IPA",
    "MPD",
    "PEG",
    "PE4",
    "PO4",
    "SO4",
    "TRS",
}
GLYCAN_RESNAMES = {
    "A2G",
    "A2M",
    "BMA",
    "FUC",
    "GAL",
    "GLC",
    "MAN",
    "NAG",
    "NDG",
    "NNN",
    "SIA",
    "SIR",
}
STANDARD_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


@dataclass(frozen=True)
class AtomRecord:
    line: str
    record: str
    atom_name: str
    resname: str
    chain: str
    resseq: str
    icode: str
    xyz: tuple[float, float, float]

    @property
    def residue_key(self) -> tuple[str, str, str, str]:
        return (self.chain, self.resseq, self.icode, self.resname)


def _read_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def _parse_xyz(line: str) -> tuple[float, float, float] | None:
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except ValueError:
        return None


def _parse_atom_line(line: str) -> AtomRecord | None:
    record = line[:6].strip()
    if record not in {"ATOM", "HETATM"}:
        return None
    xyz = _parse_xyz(line)
    if xyz is None:
        return None
    return AtomRecord(
        line=line.rstrip("\n"),
        record=record,
        atom_name=line[12:16].strip(),
        resname=line[17:20].strip().upper(),
        chain=line[21].strip(),
        resseq=line[22:26].strip(),
        icode=line[26].strip(),
        xyz=xyz,
    )


def _read_atoms(path: Path) -> list[AtomRecord]:
    atoms: list[AtomRecord] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            atom = _parse_atom_line(line)
            if atom is not None:
                atoms.append(atom)
    return atoms


def _distance2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b, strict=False))


def _centroid(atoms: list[AtomRecord]) -> tuple[float, float, float] | None:
    if not atoms:
        return None
    n = float(len(atoms))
    return (
        sum(atom.xyz[0] for atom in atoms) / n,
        sum(atom.xyz[1] for atom in atoms) / n,
        sum(atom.xyz[2] for atom in atoms) / n,
    )


def _heavy_atom_count(atoms: list[AtomRecord]) -> int:
    count = 0
    for atom in atoms:
        element = atom.line[76:78].strip().upper() or atom.atom_name[:1].upper()
        if element not in {"H", "D"}:
            count += 1
    return count


def _is_excluded_ligand_resname(resname: str, *, allow_glycans: bool) -> bool:
    if resname in WATER_RESNAMES or resname in ION_RESNAMES:
        return True
    if resname in COMMON_BUFFER_RESNAMES or resname in STANDARD_AA:
        return True
    return bool(resname in GLYCAN_RESNAMES and not allow_glycans)


def _choose_ligand_center(
    atoms: list[AtomRecord],
    *,
    min_heavy_atoms: int,
    allow_glycans: bool,
) -> tuple[tuple[float, float, float], str, str] | None:
    ligands: dict[tuple[str, str, str, str], list[AtomRecord]] = defaultdict(list)
    for atom in atoms:
        if atom.record != "HETATM":
            continue
        if _is_excluded_ligand_resname(atom.resname, allow_glycans=allow_glycans):
            continue
        ligands[atom.residue_key].append(atom)
    candidates: list[tuple[int, tuple[str, str, str, str], list[AtomRecord]]] = []
    for key, ligand_atoms in ligands.items():
        heavy_atoms = _heavy_atom_count(ligand_atoms)
        if heavy_atoms >= min_heavy_atoms:
            candidates.append((heavy_atoms, key, ligand_atoms))
    if not candidates:
        return None
    heavy_atoms, key, ligand_atoms = sorted(candidates, key=lambda item: item[0], reverse=True)[0]
    center = _centroid(ligand_atoms)
    if center is None:
        return None
    chain, resseq, icode, resname = key
    ligand_id = f"{resname}_{chain or '_'}{resseq}{icode or ''}"
    return center, ligand_id, f"cocrystal_ligand_heavy_atoms={heavy_atoms}"


def _geometric_center(atoms: list[AtomRecord]) -> tuple[float, float, float] | None:
    protein = [atom for atom in atoms if atom.record == "ATOM"]
    return _centroid(protein)


def _download_pdb(pdb_id: str, out_path: Path, *, timeout_sec: int = 60) -> bool:
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=timeout_sec) as response:
            payload = response.read()
    except Exception:
        return False
    if not payload.startswith((b"HEADER", b"TITLE", b"ATOM", b"HETATM")):
        return False
    out_path.write_bytes(payload)
    return True


def _ensure_raw_pdb(
    pdb_id: str,
    raw_root: Path,
    *,
    allow_download: bool,
) -> tuple[Path | None, str]:
    raw_path = raw_root / f"{pdb_id}.pdb"
    if raw_path.exists():
        return raw_path, "local_raw_pdb"
    if allow_download and _download_pdb(pdb_id, raw_path):
        return raw_path, "rcsb_download"
    return None, "missing_raw_pdb"


def _pdb_ids_from_pair_table(pair_table: pd.DataFrame) -> list[str]:
    if "pdb_id" not in pair_table.columns:
        raise ValueError("pair table must contain pdb_id")
    return sorted(pair_table["pdb_id"].dropna().astype(str).str.upper().unique())


def _candidate_receptor_paths(pdb_id: str, search_dirs: list[Path]) -> list[Path]:
    names = [
        f"{pdb_id}_cleaned.pdb",
        f"{pdb_id}.pdb",
        "receptor.pdb",
    ]
    candidates: list[Path] = []
    for root in search_dirs:
        candidates.extend(root / name for name in names)
        candidates.extend(root / pdb_id / name for name in names)
        candidates.extend(root / pdb_id.upper() / name for name in names)
        # Common Atlas receptor-prep layouts without a repository-wide recursive glob.
        for variant in ("HOLO", "APO"):
            candidates.extend(root / pdb_id / variant / "receptor" / "ph_ensemble" / name for name in names)
            candidates.extend(root / pdb_id.upper() / variant / "receptor" / "ph_ensemble" / name for name in names)
    paths = []
    seen: set[Path] = set()
    for path in candidates:
        if path.exists():
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                paths.append(resolved)
    return paths


def _select_pocket_atoms(
    atoms: list[AtomRecord],
    center: tuple[float, float, float],
    *,
    radius_a: float,
    include_hetatm: bool,
) -> list[AtomRecord]:
    radius2 = radius_a * radius_a
    residue_keys = {
        atom.residue_key
        for atom in atoms
        if atom.record == "ATOM" and _distance2(atom.xyz, center) <= radius2
    }
    selected = [
        atom
        for atom in atoms
        if atom.record == "ATOM" and atom.residue_key in residue_keys
    ]
    if include_hetatm:
        selected.extend(
            atom
            for atom in atoms
            if atom.record == "HETATM"
            and not _is_excluded_ligand_resname(atom.resname, allow_glycans=False)
            and _distance2(atom.xyz, center) <= radius2
        )
    return selected


def _write_pocket(path: Path, atoms: list[AtomRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for atom in atoms:
            handle.write(f"{atom.line}\n")
        handle.write("END\n")


def _round_center(center: tuple[float, float, float] | None) -> list[float] | None:
    if center is None:
        return None
    return [round(float(value), 3) for value in center]


def _ligand_center_from_atoms(
    atoms: list[AtomRecord],
    *,
    min_ligand_heavy_atoms: int,
    allow_glycan_centers: bool,
) -> tuple[tuple[float, float, float], str, str] | None:
    return _choose_ligand_center(
        atoms,
        min_heavy_atoms=min_ligand_heavy_atoms,
        allow_glycans=allow_glycan_centers,
    )


def build_banana_pocket_map(
    pair_table_path: str | Path,
    out_path: str | Path,
    *,
    raw_pdb_dir: str | Path = "data/pilotstudy/banana_pockets/raw_pdbs",
    pocket_dir: str | Path = "data/pilotstudy/banana_pockets/pockets",
    receptor_search_dirs: list[str | Path] | None = None,
    radius_a: float = 8.0,
    min_ligand_heavy_atoms: int = 6,
    allow_download: bool = True,
    allow_glycan_centers: bool = False,
    allow_geometric_fallback: bool = True,
    include_hetatm: bool = False,
    require_complete: bool = True,
) -> pd.DataFrame:
    """Build one BANANA pocket-residue PDB per PDB ID in a pair table.

    The preferred center is the centroid of a non-buffer, non-ion co-crystal
    ligand in the source PDB. A geometric fallback is explicit and should only
    be used to keep large runs complete while preserving provenance.
    """

    pair_table = _read_table(pair_table_path)
    pdb_ids = _pdb_ids_from_pair_table(pair_table)
    raw_root = Path(raw_pdb_dir)
    pocket_root = Path(pocket_dir)
    search_roots = [Path(path) for path in (receptor_search_dirs or ["outputs/processed_pdbs"])]
    rows: list[dict[str, Any]] = []
    for pdb_id in pdb_ids:
        receptor_path: Path | None = None
        source = "missing"
        existing = _candidate_receptor_paths(pdb_id, search_roots)
        if existing:
            receptor_path = existing[0]
            source = "local_receptor"
        else:
            raw_path, raw_source = _ensure_raw_pdb(
                pdb_id,
                raw_root,
                allow_download=allow_download,
            )
            if raw_path is not None:
                receptor_path = raw_path
                source = raw_source
        if receptor_path is None or not receptor_path.exists():
            rows.append({"pdb_id": pdb_id, "status": "missing_receptor_pdb"})
            continue

        atoms = _read_atoms(receptor_path)
        ligand = _choose_ligand_center(
            atoms,
            min_heavy_atoms=min_ligand_heavy_atoms,
            allow_glycans=allow_glycan_centers,
        )
        center: tuple[float, float, float] | None = None
        center_source = ""
        ligand_id = ""
        center_detail = ""
        center_reference_pdb = str(receptor_path)
        if ligand is not None:
            center, ligand_id, center_detail = ligand
            center_source = "cocrystal_ligand"
        elif source == "local_receptor":
            raw_path, raw_source = _ensure_raw_pdb(
                pdb_id,
                raw_root,
                allow_download=allow_download,
            )
            if raw_path is not None:
                raw_ligand = _ligand_center_from_atoms(
                    _read_atoms(raw_path),
                    min_ligand_heavy_atoms=min_ligand_heavy_atoms,
                    allow_glycan_centers=allow_glycan_centers,
                )
                if raw_ligand is not None:
                    center, ligand_id, center_detail = raw_ligand
                    center_source = f"raw_pdb_cocrystal_ligand:{raw_source}"
                    center_reference_pdb = str(raw_path)
        if center is None and allow_geometric_fallback:
            center = _geometric_center(atoms)
            center_source = "protein_geometric_fallback"
            center_detail = "no_suitable_cocrystal_ligand"
        if center is None:
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "status": "missing_pocket_center",
                    "receptor_pdb": str(receptor_path),
                    "receptor_source": source,
                }
            )
            continue

        pocket_atoms = _select_pocket_atoms(
            atoms,
            center,
            radius_a=radius_a,
            include_hetatm=include_hetatm,
        )
        if not pocket_atoms:
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "status": "empty_pocket",
                    "receptor_pdb": str(receptor_path),
                    "receptor_source": source,
                    "center_source": center_source,
                    "center_xyz": json.dumps(_round_center(center)),
                }
            )
            continue
        pocket_path = pocket_root / f"{pdb_id}_banana_pocket_r{radius_a:g}.pdb"
        _write_pocket(pocket_path, pocket_atoms)
        rows.append(
            {
                "pdb_id": pdb_id,
                "pocket_pdb": str(pocket_path),
                "status": "ready",
                "receptor_pdb": str(receptor_path),
                "receptor_source": source,
                "center_reference_pdb": center_reference_pdb,
                "center_source": center_source,
                "center_ligand_id": ligand_id,
                "center_detail": center_detail,
                "center_x": center[0],
                "center_y": center[1],
                "center_z": center[2],
                "radius_a": float(radius_a),
                "pocket_atom_count": int(len(pocket_atoms)),
                "pocket_residue_count": int(len({atom.residue_key for atom in pocket_atoms})),
            }
        )
    out = pd.DataFrame(rows)
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    expected_pdbs = int(len(pdb_ids))
    ready_pdbs = int((out.get("status", pd.Series(dtype=str)) == "ready").sum())
    summary = {
        "rows": int(len(out)),
        "expected_pdbs": expected_pdbs,
        "ready_pdbs": ready_pdbs,
        "status_counts": {
            str(key): int(value)
            for key, value in out.get("status", pd.Series(dtype=str)).value_counts().to_dict().items()
        },
        "center_source_counts": {
            str(key): int(value)
            for key, value in out.get("center_source", pd.Series(dtype=str)).fillna("").value_counts().to_dict().items()
            if key
        },
        "radius_a": float(radius_a),
        "min_ligand_heavy_atoms": int(min_ligand_heavy_atoms),
        "allow_glycan_centers": bool(allow_glycan_centers),
        "allow_geometric_fallback": bool(allow_geometric_fallback),
    }
    path.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    if require_complete and ready_pdbs != expected_pdbs:
        missing = out.loc[out["status"].astype(str).ne("ready"), "pdb_id"].astype(str).tolist()
        raise RuntimeError(f"BANANA pocket map incomplete: {len(missing)} missing PDBs: {', '.join(missing)}")
    return out
