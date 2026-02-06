"""Metal-site donor coordination auditing helpers."""

from __future__ import annotations

import json
import logging
import os
from math import sqrt
from pathlib import Path
from typing import Collection, Mapping, Optional, Sequence

from pdb_fixer import load_canonical_metals, load_canonical_waters


def _parse_atoms_from_pdb_like_lines(
    lines: list[str], canonical_metals: Collection[str] | None = None
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Parse ATOM/HETATM-like records and return (all_atoms, metal_atoms)."""
    parsed_atoms: list[dict[str, object]] = []
    metal_atoms: list[dict[str, object]] = []
    canon_metals = set(canonical_metals or [])
    for ln in lines:
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        try:
            x = float(ln[30:38])
            y = float(ln[38:46])
            z = float(ln[46:54])
        except Exception:
            continue
        record = ln[:6].strip().upper()
        resname = ln[17:20].strip().upper()
        chain = ln[21:22]
        resseq = ln[22:26].strip()
        icode = ln[26:27].strip()
        atom_name = ln[12:16].strip()
        element = ln[76:78].strip().upper() or atom_name[:2].strip().upper()
        atom_info = {
            "record": record,
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "atom_name": atom_name,
            "coords": (x, y, z),
            "element": element,
        }
        parsed_atoms.append(atom_info)
        is_metal = (element in canon_metals) or (resname in canon_metals)
        if record == "HETATM" and is_metal:
            metal_atoms.append(atom_info)
    return parsed_atoms, metal_atoms


def _find_metal_donors(
    metal: Mapping[str, object],
    parsed_atoms: Sequence[Mapping[str, object]],
    canonical_waters: Collection[str] | None = None,
) -> list[dict[str, object]]:
    """Identify N/O/S donors within 3.2 Å of a metal center."""
    donors: list[dict[str, object]] = []
    coords = metal.get("coords")
    if not coords:
        return donors
    try:
        mx, my, mz = coords  # type: ignore[misc]
    except Exception:
        return donors
    water_tokens = {str(w).strip().upper() for w in (canonical_waters or [])}
    for atom in parsed_atoms:
        atom_coords = atom.get("coords")
        if not atom_coords:
            continue
        try:
            ax, ay, az = atom_coords  # type: ignore[misc]
        except Exception:
            continue
        element = str(atom.get("element", "")).upper()
        if element not in {"N", "O", "S"}:
            continue
        dx = mx - ax
        dy = my - ay
        dz = mz - az
        dist = sqrt(dx * dx + dy * dy + dz * dz)
        if dist > 3.2:
            continue
        record = str(atom.get("record", "")).upper()
        resname = str(atom.get("resname", "")).upper()
        chain = str(atom.get("chain", ""))
        resseq = str(atom.get("resseq", ""))
        atom_name = str(atom.get("atom_name", ""))
        if record == "ATOM":
            category = "protein"
        elif resname in water_tokens:
            category = "water"
        else:
            category = "ligand"
        donors.append(
            {
                "category": category,
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "atom_name": atom_name,
                "distance": round(dist, 3),
            }
        )
    return donors


def run_metal_site_audit(
    *,
    pdb_id: str,
    router_paths,
    input_pdb_path: str,
    receptor_pdb_path: str | None,
    receptor_pdbqt_path: str | None,
    center: Optional[tuple[float, float, float]] = None,
    variant_label: str | None = None,
    ph_label: str | None = None,
) -> None:
    """Build a single metal_site_audit.json comparing input vs receptor donors."""

    try:
        canonical_metals = set(load_canonical_metals(None))
        canonical_waters = set(load_canonical_waters(None))
    except Exception as exc:
        logging.warning(
            "[holo.metal_audit] skip reason=canonical_load_failed pdb=%s err=%s",
            pdb_id,
            exc,
        )
        return

    try:
        input_lines = (
            Path(input_pdb_path)
            .read_text(encoding="utf-8", errors="ignore")
            .splitlines()
        )
    except Exception as exc:
        logging.warning(
            "[holo.metal_audit] skip reason=input_read_failed pdb=%s err=%s",
            pdb_id,
            exc,
        )
        return

    parsed_atoms_pre, metal_atoms_pre = _parse_atoms_from_pdb_like_lines(
        input_lines, canonical_metals
    )

    def _metal_key(atom: Mapping[str, object]) -> tuple[str, str, str, str]:
        return (
            str(atom.get("element", "")).upper(),
            str(atom.get("chain", "")),
            str(atom.get("resseq", "")),
            str(atom.get("atom_name", "")),
        )

    def _coords_list(atom: Mapping[str, object] | None) -> list[float] | None:
        if not atom:
            return None
        coords = atom.get("coords")
        if not coords:
            return None
        try:
            ax, ay, az = coords  # type: ignore[misc]
        except Exception:
            return None
        return [float(ax), float(ay), float(az)]

    metal_pre_by_key: dict[tuple[str, str, str, str], Mapping[str, object]] = {}
    donors_input_map: dict[tuple[str, str, str, str], list[dict[str, object]]] = {}
    for metal in metal_atoms_pre:
        key = _metal_key(metal)
        metal_pre_by_key[key] = metal
        donors_input_map[key] = _find_metal_donors(
            metal, parsed_atoms_pre, canonical_waters
        )

    donors_receptor_map: dict[tuple[str, str, str, str], list[dict[str, object]]] = {}
    metal_post_by_key: dict[tuple[str, str, str, str], Mapping[str, object]] = {}
    if receptor_pdbqt_path:
        pdbqt_lines: list[str] = []
        try:
            pdbqt_lines = (
                Path(receptor_pdbqt_path)
                .read_text(encoding="utf-8", errors="ignore")
                .splitlines()
            )
        except Exception as exc:
            logging.warning(
                "[holo.metal_audit] skip reason=receptor_read_failed pdb=%s file=%s err=%s",
                pdb_id,
                receptor_pdbqt_path,
                exc,
            )
        if pdbqt_lines:
            parsed_atoms_post, metal_atoms_post = _parse_atoms_from_pdb_like_lines(
                pdbqt_lines, canonical_metals
            )
            for metal in metal_atoms_post:
                key = _metal_key(metal)
                metal_post_by_key[key] = metal
                donors_receptor_map[key] = _find_metal_donors(
                    metal, parsed_atoms_post, canonical_waters
                )

    all_keys = sorted(set(metal_pre_by_key) | set(metal_post_by_key))

    def _donor_key(entry: Mapping[str, object]) -> tuple[str, str, str, str, str]:
        return (
            str(entry.get("category", "")),
            str(entry.get("resname", "")),
            str(entry.get("chain", "")),
            str(entry.get("resseq", "")),
            str(entry.get("atom_name", "")),
        )

    metals_payload: list[dict[str, object]] = []
    for key in all_keys:
        donors_input = donors_input_map.get(key, [])
        donors_receptor = donors_receptor_map.get(key, [])
        input_set = {_donor_key(d) for d in donors_input}
        receptor_set = {_donor_key(d) for d in donors_receptor}
        lost_keys = input_set - receptor_set
        gained_keys = receptor_set - input_set
        lost_donors = [d for d in donors_input if _donor_key(d) in lost_keys]
        gained_donors = [d for d in donors_receptor if _donor_key(d) in gained_keys]

        metal_source = metal_pre_by_key.get(key) or metal_post_by_key.get(key) or {}
        resname = str(metal_source.get("resname", ""))
        chain = str(metal_source.get("chain", ""))
        resseq = str(metal_source.get("resseq", ""))
        icode = str(metal_source.get("icode", ""))
        atom_name = str(metal_source.get("atom_name", ""))
        element = str(metal_source.get("element", key[0] if key else ""))
        coords_input = _coords_list(metal_pre_by_key.get(key))
        coords_receptor = _coords_list(metal_post_by_key.get(key))
        metals_payload.append(
            {
                "id": f"{element or 'UNK'}_{chain or '-'}_{resseq or '-'}_{atom_name or '-'}",
                "element": element,
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "atom_name": atom_name,
                "coords_input": coords_input,
                "coords_receptor": coords_receptor,
                "donors_input": donors_input,
                "donors_receptor": donors_receptor,
                "lost_donors": lost_donors,
                "gained_donors": gained_donors,
            }
        )

    total_input = 0
    total_receptor = 0
    per_metal_summary: list[dict[str, object]] = []
    for metal in metals_payload:
        donors_input = metal.get("donors_input", []) or []
        donors_receptor = metal.get("donors_receptor", []) or []
        input_count = len(donors_input)
        receptor_count = len(donors_receptor)
        total_input += input_count
        total_receptor += receptor_count
        per_metal_summary.append(
            {
                "id": metal.get("id"),
                "input_count": input_count,
                "receptor_count": receptor_count,
                "delta": input_count - receptor_count,
            }
        )

    coordination_summary = {
        "total_input_count": total_input,
        "total_receptor_count": total_receptor,
        "total_delta": total_input - total_receptor,
        "per_metal": per_metal_summary,
    }

    variant_for_path = variant_label or "HOLO"
    if router_paths is None:
        logging.warning(
            "[holo.metal_audit] skip reason=router_missing pdb=%s variant=%s",
            pdb_id,
            variant_for_path,
        )
        return

    try:
        docked_variant_dir = router_paths.docked_variant_root(
            variant_for_path, ph_label=ph_label
        )
        docked_variant_dir.mkdir(parents=True, exist_ok=True)
    except Exception as exc:
        logging.warning(
            "[holo.metal_audit] skip reason=path_resolve_failed pdb=%s err=%s",
            pdb_id,
            exc,
        )
        return

    payload = {
        "pdb_id": pdb_id,
        "variant": variant_label or variant_for_path,
        "ph_label": ph_label,
        "source_files": {
            "input_pdb": input_pdb_path,
            "receptor_pdb": receptor_pdb_path,
            "receptor_pdbqt": receptor_pdbqt_path,
        },
        "metals": metals_payload,
        "coordination_summary": coordination_summary,
    }

    json_path = docked_variant_dir / "metal_site_audit.json"
    tmp_path = json_path.with_suffix(".part")
    try:
        tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        os.replace(tmp_path, json_path)
        logging.info(
            "[holo.metal_audit] pdb=%s metals=%d file=%s",
            pdb_id,
            len(metals_payload),
            json_path,
        )
    except Exception as exc:
        logging.warning(
            "[holo.metal_audit] write_failed pdb=%s file=%s err=%s",
            pdb_id,
            json_path,
            exc,
        )
