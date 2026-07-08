"""Metal-site donor coordination auditing helpers."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from math import sqrt
from pathlib import Path
from typing import Collection, Mapping, Optional, Sequence

from protein_prep.metal.autodock4zn import (
    build_autodock4zn_plan,
    summarize_autodock4zn_plan,
    write_autodock4zn_plan,
)
from protein_prep.metal.chemistry import audit_metal_chemistry, summarize_metal_chemistry
from protein_prep.metal.parameterization import (
    build_mcpb_parameterization_plan,
    summarize_parameterization_plan,
    write_parameterization_plan,
)
from protein_prep.pdb_fixer_runtime import load_canonical_metals, load_canonical_waters


def _write_text_atomic(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            delete=False,
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".part",
            encoding="utf-8",
        ) as handle:
            tmp_path = handle.name
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
        tmp_path = None
    finally:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except FileNotFoundError:
                pass


def _coords3(value: object) -> tuple[float, float, float] | None:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        return None
    try:
        x, y, z = value
        return float(x), float(y), float(z)
    except Exception:
        return None


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
        serial = ln[6:11].strip()
        element = ln[76:78].strip().upper() or atom_name[:2].strip().upper()
        atom_info: dict[str, object] = {
            "serial": serial,
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
    coords = _coords3(metal.get("coords"))
    if coords is None:
        return donors
    mx, my, mz = coords
    water_tokens = {str(w).strip().upper() for w in (canonical_waters or [])}
    for atom in parsed_atoms:
        atom_coords = _coords3(atom.get("coords"))
        if atom_coords is None:
            continue
        ax, ay, az = atom_coords
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
                "element": element,
                "coords": [round(ax, 3), round(ay, 3), round(az, 3)],
                "distance": round(dist, 3),
            }
        )
    return donors


def _preferred_metal_chemistry(
    receptor: Mapping[str, object],
    cleaned: Mapping[str, object],
    input_row: Mapping[str, object],
) -> Mapping[str, object]:
    for candidate in (receptor, cleaned, input_row):
        try:
            if int(str(candidate.get("coordination_number", 0) or 0)) > 0:
                return candidate
        except Exception:
            continue
    return receptor or cleaned or input_row


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
        coords = _coords3(atom.get("coords"))
        if coords is None:
            return None
        ax, ay, az = coords
        return [float(ax), float(ay), float(az)]

    metal_pre_by_key: dict[tuple[str, str, str, str], Mapping[str, object]] = {}
    donors_input_map: dict[tuple[str, str, str, str], list[dict[str, object]]] = {}
    for metal in metal_atoms_pre:
        key = _metal_key(metal)
        metal_pre_by_key[key] = metal
        donors_input_map[key] = _find_metal_donors(
            metal, parsed_atoms_pre, canonical_waters
        )

    def _stage_maps(
        stage_path: str | None, stage_name: str
    ) -> tuple[
        dict[tuple[str, str, str, str], Mapping[str, object]],
        dict[tuple[str, str, str, str], list[dict[str, object]]],
    ]:
        metal_by_key: dict[tuple[str, str, str, str], Mapping[str, object]] = {}
        donors_by_key: dict[tuple[str, str, str, str], list[dict[str, object]]] = {}
        if not stage_path:
            return metal_by_key, donors_by_key
        stage_lines: list[str] = []
        try:
            stage_lines = (
                Path(stage_path)
                .read_text(encoding="utf-8", errors="ignore")
                .splitlines()
            )
        except Exception as exc:
            logging.warning(
                "[holo.metal_audit] stage=%s read_failed pdb=%s file=%s err=%s",
                stage_name,
                pdb_id,
                stage_path,
                exc,
            )
            return metal_by_key, donors_by_key
        if not stage_lines:
            return metal_by_key, donors_by_key
        parsed_atoms_stage, metal_atoms_stage = _parse_atoms_from_pdb_like_lines(
            stage_lines, canonical_metals
        )
        for metal in metal_atoms_stage:
            key = _metal_key(metal)
            metal_by_key[key] = metal
            donors_by_key[key] = _find_metal_donors(
                metal, parsed_atoms_stage, canonical_waters
            )
        return metal_by_key, donors_by_key

    metal_cleaned_by_key, donors_cleaned_map = _stage_maps(
        receptor_pdb_path, "cleaned_pdb"
    )
    metal_post_by_key, donors_receptor_map = _stage_maps(
        receptor_pdbqt_path, "receptor_pdbqt"
    )

    all_keys = sorted(
        set(metal_pre_by_key) | set(metal_cleaned_by_key) | set(metal_post_by_key)
    )

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
        donors_cleaned = donors_cleaned_map.get(key, [])
        donors_receptor = donors_receptor_map.get(key, [])
        input_set = {_donor_key(d) for d in donors_input}
        cleaned_set = {_donor_key(d) for d in donors_cleaned}
        receptor_set = {_donor_key(d) for d in donors_receptor}
        lost_cleaned_keys = input_set - cleaned_set
        gained_cleaned_keys = cleaned_set - input_set
        lost_keys = input_set - receptor_set
        gained_keys = receptor_set - input_set
        lost_cleaned_donors = [
            d for d in donors_input if _donor_key(d) in lost_cleaned_keys
        ]
        gained_cleaned_donors = [
            d for d in donors_cleaned if _donor_key(d) in gained_cleaned_keys
        ]
        lost_donors = [d for d in donors_input if _donor_key(d) in lost_keys]
        gained_donors = [d for d in donors_receptor if _donor_key(d) in gained_keys]

        metal_source = (
            metal_pre_by_key.get(key)
            or metal_cleaned_by_key.get(key)
            or metal_post_by_key.get(key)
            or {}
        )
        resname = str(metal_source.get("resname", ""))
        chain = str(metal_source.get("chain", ""))
        resseq = str(metal_source.get("resseq", ""))
        icode = str(metal_source.get("icode", ""))
        atom_name = str(metal_source.get("atom_name", ""))
        serial = str(metal_source.get("serial", ""))
        element = str(metal_source.get("element", key[0] if key else ""))
        coords_input = _coords_list(metal_pre_by_key.get(key))
        coords_cleaned = _coords_list(metal_cleaned_by_key.get(key))
        coords_receptor = _coords_list(metal_post_by_key.get(key))
        chemistry_input = audit_metal_chemistry(
            element=element,
            metal_xyz=_coords3(metal_pre_by_key.get(key, {}).get("coords")),
            donors=donors_input,
        )
        chemistry_cleaned = audit_metal_chemistry(
            element=element,
            metal_xyz=_coords3(metal_cleaned_by_key.get(key, {}).get("coords")),
            donors=donors_cleaned,
        )
        chemistry_receptor = audit_metal_chemistry(
            element=element,
            metal_xyz=_coords3(metal_post_by_key.get(key, {}).get("coords")),
            donors=donors_receptor,
        )
        metals_payload.append(
            {
                "id": f"{element or 'UNK'}_{chain or '-'}_{resseq or '-'}_{atom_name or '-'}",
                "element": element,
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "atom_name": atom_name,
                "serial": serial,
                "coords_input": coords_input,
                "coords_cleaned_pdb": coords_cleaned,
                "coords_receptor": coords_receptor,
                "donors_input": donors_input,
                "donors_cleaned_pdb": donors_cleaned,
                "donors_receptor": donors_receptor,
                "lost_donors_cleaned_pdb": lost_cleaned_donors,
                "gained_donors_cleaned_pdb": gained_cleaned_donors,
                "lost_donors": lost_donors,
                "gained_donors": gained_donors,
                "chemistry_input": chemistry_input,
                "chemistry_cleaned_pdb": chemistry_cleaned,
                "chemistry_receptor": chemistry_receptor,
                "chemistry": _preferred_metal_chemistry(
                    chemistry_receptor,
                    chemistry_cleaned,
                    chemistry_input,
                ),
            }
        )

    total_input = 0
    total_cleaned = 0
    total_receptor = 0
    per_metal_summary: list[dict[str, object]] = []
    for metal in metals_payload:
        raw_donors_input = metal.get("donors_input", []) or []
        raw_donors_cleaned = metal.get("donors_cleaned_pdb", []) or []
        raw_donors_receptor = metal.get("donors_receptor", []) or []
        donor_input_rows: list[dict[str, object]]
        donor_cleaned_rows: list[dict[str, object]]
        donor_receptor_rows: list[dict[str, object]]
        donor_input_rows = raw_donors_input if isinstance(raw_donors_input, list) else []
        donor_cleaned_rows = (
            raw_donors_cleaned if isinstance(raw_donors_cleaned, list) else []
        )
        donor_receptor_rows = (
            raw_donors_receptor if isinstance(raw_donors_receptor, list) else []
        )
        input_count = len(donor_input_rows)
        cleaned_count = len(donor_cleaned_rows)
        receptor_count = len(donor_receptor_rows)
        total_input += input_count
        total_cleaned += cleaned_count
        total_receptor += receptor_count
        per_metal_summary.append(
            {
                "id": metal.get("id"),
                "input_count": input_count,
                "cleaned_pdb_count": cleaned_count,
                "receptor_count": receptor_count,
                "cleaned_pdb_delta": input_count - cleaned_count,
                "pdbqt_delta": input_count - receptor_count,
                "delta": input_count - receptor_count,
            }
        )

    coordination_summary = {
        "total_input_count": total_input,
        "total_cleaned_pdb_count": total_cleaned,
        "total_receptor_count": total_receptor,
        "total_cleaned_pdb_delta": total_input - total_cleaned,
        "total_pdbqt_delta": total_input - total_receptor,
        "total_delta": total_input - total_receptor,
        "per_metal": per_metal_summary,
    }
    chemistry_rows: list[Mapping[str, object]] = []
    for metal in metals_payload:
        chemistry = metal.get("chemistry", {})
        if isinstance(chemistry, Mapping):
            chemistry_row = dict(chemistry)
            chemistry_row.update(
                {
                    "id": metal.get("id", ""),
                    "atom_serial": metal.get("serial", ""),
                    "resname": metal.get("resname", ""),
                    "chain": metal.get("chain", ""),
                    "resseq": metal.get("resseq", ""),
                    "atom_name": metal.get("atom_name", ""),
                    "coords": metal.get("coords_cleaned_pdb")
                    or metal.get("coords_input")
                    or metal.get("coords_receptor"),
                    "donors": metal.get("donors_cleaned_pdb")
                    or metal.get("donors_input")
                    or metal.get("donors_receptor")
                    or [],
                }
            )
            chemistry_rows.append(chemistry_row)
    chemistry_summary = summarize_metal_chemistry(chemistry_rows)

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
        "chemistry_summary": chemistry_summary,
    }
    ad4zn_plan_path = docked_variant_dir / "autodock4zn_plan.json"
    ad4zn_plan = build_autodock4zn_plan(
        receptor_pdbqt=Path(receptor_pdbqt_path) if receptor_pdbqt_path else None,
        metal_rows=chemistry_rows,
        work_dir=docked_variant_dir / "autodock4zn",
        center=center,
    )
    ad4zn_plan["plan_path"] = str(ad4zn_plan_path)
    parameterization_plan_path = docked_variant_dir / "metal_parameterization_plan.json"
    parameterization_plan = build_mcpb_parameterization_plan(
        receptor_pdb=Path(receptor_pdb_path) if receptor_pdb_path else None,
        metal_rows=chemistry_rows,
        work_dir=docked_variant_dir / "metal_parameterization",
    )
    parameterization_plan["plan_path"] = str(parameterization_plan_path)
    payload["autodock4zn_summary"] = summarize_autodock4zn_plan(ad4zn_plan)
    payload["autodock4zn_plan"] = ad4zn_plan
    payload["metal_parameterization_summary"] = summarize_parameterization_plan(
        parameterization_plan
    )
    payload["metal_parameterization_plan"] = parameterization_plan

    json_path = docked_variant_dir / "metal_site_audit.json"
    try:
        write_autodock4zn_plan(ad4zn_plan_path, ad4zn_plan)
        write_parameterization_plan(parameterization_plan_path, parameterization_plan)
        _write_text_atomic(json_path, json.dumps(payload, indent=2))
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
