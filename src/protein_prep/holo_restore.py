"""HOLO restore helpers."""

from __future__ import annotations

import logging
import os
from collections import defaultdict
from math import sqrt
from pathlib import Path
from typing import Any, Mapping, Optional, Tuple

from protein_prep.pdb_fixer_runtime import (
    fix_element_columns_in_file,
    load_canonical_cofactors,
    load_canonical_metals,
    load_canonical_waters,
)
from path_router import make_paths
from protein_prep.metal_site_audit import (
    _find_metal_donors,
    _parse_atoms_from_pdb_like_lines,
)
from protein_prep.prep_utils import _resolve_variant_token


def _hydrate_legacy_globals() -> None:
    """Compatibility shim retained for older tests/callers."""
    return None


def _coords3(value: object) -> tuple[float, float, float] | None:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        return None
    try:
        x, y, z = value
        return float(x), float(y), float(z)
    except Exception:
        return None


def _holo_restore_from_input_if_needed(
    pdb_id: str,
    cleaned_pdb: str,
    output_pdbqt: str,
    config: Mapping[str, Any],
    center: Optional[Tuple[float, float, float]] = None,
    box_size: Optional[Tuple[float, float, float]] = None,
) -> Tuple[int, int, bool]:
    """
    HOLO-only: re-add missing metals/cofactors from input_pdbs/<PDB>.pdb into the
    final cleaned PDB, normalize element columns, and regenerate PDBQT once.

    Returns: (metals_added, cofactors_added, regenerated_pdbqt)
    """
    logging.info(
        "[holo.box.debug] pdb=%s center=%s box_size=%s",
        pdb_id,
        center,
        box_size,
    )
    try:
        logging.info("[holo.restore] precheck pdb=%s", pdb_id)

        # Resolve APO/HOLO mode (strict)
        mode = (
            (os.environ.get("APO_HOLO_MODE") or config.get("APO_HOLO_MODE") or "")
            .strip()
            .lower()
        )
        mode_norm = {
            "": "legacy",
            "none": "legacy",
            "null": "legacy",
            "false": "legacy",
            "0": "legacy",
        }.get(mode, mode)
        if mode_norm in {"legacy", "apo"}:
            logging.warning("[holo.restore] skip reason=mode=%s", mode_norm)
            return (0, 0, False)

        # Resolve variant
        variant_token = (_resolve_variant_token(config) or "").strip().upper()
        if variant_token != "HOLO":
            logging.warning(
                "[holo.restore] skip reason=variant=%s", variant_token or "None"
            )
            return (0, 0, False)

        # Locate original input PDB via existing router
        router_paths = None
        try:
            router_paths = make_paths(config, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            input_pdb = str(router_paths.input_pdb_path)
        except Exception:
            input_root = config.get("INPUT_DIR", "input_pdbs")
            input_pdb = str(Path(input_root) / f"{pdb_id}.pdb")

        target_pdb = str(cleaned_pdb)

        # Load canonical sets from aliases.yaml (no hard-coded lists)
        try:
            # NOTE: pass None so loaders use chemdb/aliases.yaml via activesite._load_default_alias_cfg
            canonical_metals = set(load_canonical_metals(None))
            canonical_cofactors = set(load_canonical_cofactors(None))
            canonical_waters = set(load_canonical_waters(None))
        except Exception as e:
            logging.warning(
                "[holo.restore] skip reason=cannot_load_canonical_sets err=%s", e
            )
            return (0, 0, False)

        logging.info(
            "[holo.restore.debug] canonical_sizes metals=%d cofactors=%d waters=%d",
            len(canonical_metals),
            len(canonical_cofactors),
            len(canonical_waters),
        )

        # Scan original input PDB for candidate HETATMs
        metals_found = 0
        cofactors_found = 0
        candidates: list[str] = []
        try:
            with open(input_pdb, "r", encoding="utf-8", errors="ignore") as f:
                input_lines = f.readlines()
        except Exception as e:
            logging.warning("[holo.restore] skip reason=missing_input err=%s", e)
            return (0, 0, False)

        hetatm_total = sum(1 for ln in input_lines if ln.startswith("HETATM"))
        logging.info(
            "[holo.restore.debug] input_hetatm total=%d file=%s",
            hetatm_total,
            input_pdb,
        )

        # --- Metal coordination scan for ligand restoration (JSON emitted later) ---
        parsed_atoms_pre, metal_atoms_pre = _parse_atoms_from_pdb_like_lines(
            input_lines, canonical_metals
        )
        parsed_atoms = parsed_atoms_pre
        donor_icode_lookup: dict[tuple[str, str, str, str], str] = {}
        for atom in parsed_atoms_pre:
            resname = str(atom.get("resname", "")).upper()
            chain = str(atom.get("chain", ""))
            resseq = str(atom.get("resseq", ""))
            atom_name = str(atom.get("atom_name", ""))
            icode = str(atom.get("icode", ""))
            donor_icode_lookup[(resname, chain, resseq, atom_name)] = icode

        coord_ligand_residues: set[tuple[str, str, str, str]] = set()
        core_water_residues: set[tuple[str, str, str]] = set()
        for metal in metal_atoms_pre:
            donors = _find_metal_donors(metal, parsed_atoms_pre, canonical_waters)
            coords_val = metal.get("coords")
            if coords_val and center is not None and len(center) == 3:
                try:
                    coords = _coords3(coords_val)
                    if coords is None:
                        raise ValueError("invalid coords")
                    mx, my, mz = coords
                    cx, cy, cz = center
                    dx = mx - cx
                    dy = my - cy
                    dz = mz - cz
                    m_dist = sqrt(dx * dx + dy * dy + dz * dz)
                    dist_term = f"{m_dist:.3f}"
                except Exception:
                    dist_term = "none"
            else:
                dist_term = "none"
            logging.info(
                "[holo.box.metal] pdb=%s metal=%s chain=%s resseq=%s center_dist=%s",
                pdb_id,
                metal.get("element", ""),
                metal.get("chain", ""),
                metal.get("resseq", ""),
                dist_term,
            )
            for donor in donors:
                resname = str(donor.get("resname", "")).upper()
                chain = str(donor.get("chain", ""))
                resseq = str(donor.get("resseq", ""))
                atom_name = str(donor.get("atom_name", ""))
                icode = donor_icode_lookup.get((resname, chain, resseq, atom_name), "")
                if donor.get("category") == "water":
                    core_water_residues.add((chain, resseq, icode))
                    continue
                if donor.get("category") != "ligand":
                    continue
                coord_ligand_residues.add((resname, chain, resseq, icode))

        def _point_in_box(pt, center_val, box_val):
            if center_val is None or box_val is None:
                return None
            try:
                cx, cy, cz = center_val
                sx, sy, sz = box_val
            except Exception:
                return None
            x, y, z = pt
            return (
                abs(x - cx) <= sx / 2.0
                and abs(y - cy) <= sy / 2.0
                and abs(z - cz) <= sz / 2.0
            )

        ligand_box_class: dict[tuple[str, str, str, str], str] = {}
        residue_atoms: dict[
            tuple[str, str, str, str], list[dict[str, object]]
        ] = defaultdict(list)
        for atom in parsed_atoms:
            if str(atom.get("record", "")).upper() != "HETATM":
                continue
            key = (
                str(atom.get("resname", "")).upper(),
                str(atom.get("chain", "")),
                str(atom.get("resseq", "")),
                str(atom.get("icode", "")),
            )
            residue_atoms[key].append(atom)

        for (
            resname,
            chain_id,
            resseq,
            icode,
        ), residue_atoms_list in residue_atoms.items():
            if (
                resname in canonical_metals
                or resname in canonical_waters
                or resname in canonical_cofactors
            ):
                continue
            classification = "out_of_box"
            has_box_info = False
            for atom in residue_atoms_list:
                coords = _coords3(atom.get("coords"))
                if coords is None:
                    continue
                inside = _point_in_box(coords, center, box_size)
                if inside is None:
                    classification = "no_box"
                    break
                has_box_info = True
                if inside:
                    classification = "in_box"
                    break
            if classification == "out_of_box" and not has_box_info:
                classification = "no_box"
            ligand_box_class[(resname, chain_id, resseq, icode)] = classification
            logging.info(
                "[holo.ligand.box] pdb=%s resname=%s chain=%s resseq=%s classification=%s num_atoms=%d",
                pdb_id,
                resname,
                chain_id,
                resseq,
                classification,
                len(residue_atoms_list),
            )

        # Build residue-present keys from the cleaned PDB to ensure idempotency
        present_keys = set()
        try:
            with open(target_pdb, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    if not ln.startswith("HETATM"):
                        continue
                    resname = ln[17:20].strip().upper()
                    chain = ln[21:22]
                    resseq = ln[22:26].strip()
                    icode = ln[26:27].strip()
                    present_keys.add((resname, chain, resseq, icode))
        except Exception as e:
            logging.warning("[holo.restore] skip reason=cleaned_unreadable err=%s", e)
            return (0, 0, False)

        # Classify metals/cofactors; treat water donors separately
        for ln in input_lines:
            if not ln.startswith("HETATM"):
                continue
            resname = ln[17:20].strip().upper()
            element = ln[76:78].strip().upper()
            chain = ln[21:22]
            resseq = ln[22:26].strip()
            icode = ln[26:27].strip()
            key = (resname, chain, resseq, icode)
            water_key = (chain, resseq, icode)
            if resname in canonical_waters:
                if water_key not in core_water_residues:
                    continue
                if key not in present_keys:
                    candidates.append(ln)
                    logging.info(
                        "[holo.water.restore] pdb=%s resname=%s chain=%s resseq=%s reason=metal_donor",
                        pdb_id,
                        resname,
                        chain,
                        resseq,
                    )
                continue
            classification = ligand_box_class.get(key, "no_box")
            is_metal = (element in canonical_metals) or (resname in canonical_metals)
            is_cofac = (resname in canonical_cofactors) and not is_metal
            is_coord_ligand = (key in coord_ligand_residues) or is_cofac
            if is_metal:
                metals_found += 1
            if is_cofac:
                cofactors_found += 1
            if is_metal:
                if key not in present_keys:
                    candidates.append(ln)
                continue
            if not is_coord_ligand:
                continue
            if classification == "in_box":
                logging.info(
                    "[holo.ligand.restore.skip] pdb=%s resname=%s chain=%s resseq=%s classification=%s reason=in_box",
                    pdb_id,
                    resname,
                    chain,
                    resseq,
                    classification,
                )
                continue
            if key not in present_keys:
                candidates.append(ln)
                logging.info(
                    "[holo.ligand.restore.add] pdb=%s resname=%s chain=%s resseq=%s classification=%s",
                    pdb_id,
                    resname,
                    chain,
                    resseq,
                    classification,
                )

        logging.info(
            "[holo.restore] source=%s target=%s metals_found=%d cofactors_found=%d",
            input_pdb,
            target_pdb,
            metals_found,
            cofactors_found,
        )

        if not candidates:
            logging.info("[holo.restore] metals_added=0 cofactors_added=0")
            logging.info("[holo.restore] regenerating_pdbqt=false out=%s", output_pdbqt)
            return (0, 0, False)

        # Append missing HETATMs (strip trailing END/TER/ENDMDL first), then ensure END
        try:
            with open(target_pdb, "r", encoding="utf-8", errors="ignore") as f:
                out_lines = f.read().splitlines()
            while out_lines and out_lines[-1].strip() in {"END", "ENDMDL", "TER"}:
                out_lines.pop()
            with open(target_pdb, "w", encoding="utf-8") as w:
                if out_lines:
                    w.write("\n".join(out_lines) + "\n")
                for ln in candidates:
                    w.write(ln.rstrip() + "\n")
                w.write("END\n")
        except Exception as e:
            logging.warning("[holo.restore] skip reason=append_failed err=%s", e)
            return (0, 0, False)

        # Normalize element columns after append
        try:
            fix_element_columns_in_file(target_pdb)
        except Exception as e:
            logging.warning("[holo.restore] skip reason=elemfix_failed err=%s", e)
            return (0, 0, False)

        metals_added = sum(
            1
            for ln in candidates
            if (ln[76:78].strip().upper() in canonical_metals)
            or (ln[17:20].strip().upper() in canonical_metals)
        )
        cofactors_added = sum(
            1
            for ln in candidates
            if (ln[17:20].strip().upper() in canonical_cofactors)
            and (ln[76:78].strip().upper() not in canonical_metals)
        )

        logging.info(
            "[holo.restore] metals_added=%d cofactors_added=%d",
            metals_added,
            cofactors_added,
        )

        # We only set regenerated_pdbqt=True here; the actual regeneration is done at the call site
        regenerated_pdbqt = bool(candidates)
        return (metals_added, cofactors_added, regenerated_pdbqt)

    except Exception as exc:
        logging.warning("[holo.restore] skip reason=unexpected err=%s", exc)
        return (0, 0, False)
