import json
import os
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.PDB import NeighborSearch, PDBParser

from protein_prep import pdb_fixer_runtime as pdb_fixer

logger = logging.getLogger(__name__)


def extract_and_remove_ligands(pdb_path, output_cleaned_pdb, ligands_dir):
    os.makedirs(ligands_dir, exist_ok=True)
    logging.info(
        "[extract] in=%s out=%s ldir=%s", pdb_path, output_cleaned_pdb, ligands_dir
    )

    ligands = defaultdict(list)
    ligand_coords = []  # collect all ligand atom coords for box calculation
    rules = pdb_fixer.get_atom_rules()
    retained_resnames = rules.retain_resnames  # already uppercased
    # [ions] retention counters
    metal_tokens = {"ZN", "MG", "MN", "FE", "CA", "CU", "CO", "NI", "HG"}
    variant_label = (os.environ.get("APO_HOLO_MODE") or "").strip().upper() or "legacy"
    kept_metals = 0
    stripped_metals = 0
    stripped_detail: List[str] = []

    # Do NOT retain common cryos/buffers: GOL/EDO/PG4/MPD/ACT/TRS/PO4/PEG → they’ll be extracted

    with open(pdb_path, "r") as infile, open(output_cleaned_pdb, "w") as outfile:
        for line in infile:
            if line.startswith("HETATM"):
                resname_raw = line[17:20]
                resname = (
                    resname_raw.strip().upper()
                )  # normalize for set membership & filenames
                chain = line[21]
                resnum = line[22:26].strip()
                elem = (line[76:78].strip() or resname).upper()
                if resname not in retained_resnames:
                    ligands[(chain, resname, resnum)].append(line)
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append((x, y, z))
                    except ValueError:
                        logging.warning(
                            f"Invalid ligand coordinates in line: {line.strip()}"
                        )
                    if elem in metal_tokens or resname in metal_tokens:
                        stripped_metals += 1
                        stripped_detail.append(
                            f"{resname}:{chain or '-'}:{resnum or '?'}"
                        )
                        logging.info(
                            "[ions.drop] reason=not_in_retain resname=%s chain=%s resSeq=%s element=%s variant=%s",
                            resname,
                            chain.strip() or "-",
                            resnum or "?",
                            elem or "?",
                            variant_label,
                        )
                    continue
                if elem in metal_tokens or resname in metal_tokens:
                    kept_metals += 1
            outfile.write(line)

    # Write separate ligand files (with element repair)
    for (chain, resname, resnum), lines in ligands.items():
        ligand_fname = os.path.join(ligands_dir, f"{resname}_{chain}{resnum}.pdb")
        fixed = pdb_fixer.fix_ligand_element_columns(lines)
        with open(ligand_fname, "w") as lf:
            lf.writelines(fixed)
        logging.info(
            f"Saved ligand {resname} {chain}{resnum} to {ligand_fname} (elements repaired)"
        )

    logging.info(f"Ligands extracted and removed from {pdb_path}.")
    logging.info(
        "[ions.summary] stage=extract_and_remove_ligands variant=%s kept=%d stripped=%d",
        variant_label,
        kept_metals,
        stripped_metals,
    )
    logging.info(
        "[ions.diff.input→cleaned] variant=%s lost=%s",
        variant_label,
        ",".join(stripped_detail) if stripped_detail else "none",
    )
    return ligands, ligand_coords


def rank_ligands_by_atom_count(ligands_dict):
    return sorted(ligands_dict.items(), key=lambda item: len(item[1]), reverse=True)


def calculate_ligand_protein_contacts(
    protein_pdb_path,
    ligand_lines,
    distance_cutoff=4.0,
):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_pdb_path)
    model = structure[0]

    # Collect all protein atoms
    protein_atoms = [
        atom for atom in model.get_atoms() if atom.get_parent().get_id()[0] == " "
    ]

    # Parse ligand atoms from lines
    ligand_atoms = []
    for line in ligand_lines:
        if line.startswith(("HETATM", "ATOM")):
            # Simple PDB atom line parsing for coords
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            # Create a dummy atom-like object
            class DummyAtom:
                def __init__(self, coord):
                    self.coord = coord

                def get_coord(self):
                    return self.coord

            ligand_atoms.append(DummyAtom((x, y, z)))

    # Use NeighborSearch to find protein atoms near ligand atoms
    ns = NeighborSearch(protein_atoms)
    contact_count = 0
    for latom in ligand_atoms:
        close_atoms = ns.search(latom.get_coord(), distance_cutoff)
        contact_count += len(close_atoms)

    return contact_count


def compute_box_from_ligand_coords(coords):
    if not coords:
        return None, None

    x_vals, y_vals, z_vals = zip(*coords)
    center = (
        sum(x_vals) / len(x_vals),
        sum(y_vals) / len(y_vals),
        sum(z_vals) / len(z_vals),
    )

    buffer = 5.0
    x_range = max(x_vals) - min(x_vals) + buffer
    y_range = max(y_vals) - min(y_vals) + buffer
    z_range = max(z_vals) - min(z_vals) + buffer

    box_size = (x_range, y_range, z_range)
    return center, box_size


def parse_hetnam_map(pdb_path: str) -> dict[str, str]:
    """
    Parse HETNAM records into {HET_ID: concatenated_name}, uppercasing names for matching.
    """
    het_parts: dict[str, List[str]] = defaultdict(list)
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("HETNAM"):
                    continue
                het_id = line[11:14].strip().upper()
                if not het_id:
                    continue
                name_part = line[15:].strip()
                if name_part:
                    het_parts[het_id].append(name_part)
    except FileNotFoundError:
        return {}
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("[hetnam.parse] failed for %s err=%s", pdb_path, exc)
        return {}

    hetnam_map: dict[str, str] = {}
    for het_id, parts in het_parts.items():
        joined = " ".join(p.strip() for p in parts if p.strip())
        hetnam_map[het_id] = joined.upper() if joined else ""
    return hetnam_map


def _coords_from_lines(lines: List[str]) -> list[Tuple[float, float, float]]:
    coords: list[Tuple[float, float, float]] = []
    for line in lines:
        if not line.startswith(("HETATM", "ATOM  ")):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append((x, y, z))
        except ValueError:
            logger.warning("[pockets.json] invalid coord line skipped: %s", line.strip())
    return coords


def build_ligand_pockets_manifest(
    pdb_path: str, ligands_dict: Dict[tuple, List[str]], logger=logging
) -> dict:
    exclusion = pdb_fixer.get_ligand_exclusion_spec()
    hetnam_map = parse_hetnam_map(pdb_path)

    raw_exclude_ids = exclusion.get("exclude_het_ids", set())
    raw_exclude_kw = exclusion.get("exclude_het_name_keywords", [])
    raw_solvent_res = exclusion.get("remove_as_solvent_resnames", set())
    raw_glycan_res = exclusion.get("remove_as_glycan_resnames", set())
    exclude_ids = set(raw_exclude_ids) if isinstance(raw_exclude_ids, (set, list, tuple)) else set()
    exclude_kw = list(raw_exclude_kw) if isinstance(raw_exclude_kw, (list, tuple, set)) else []
    solvent_res = set(raw_solvent_res) if isinstance(raw_solvent_res, (set, list, tuple)) else set()
    glycan_res = set(raw_glycan_res) if isinstance(raw_glycan_res, (set, list, tuple)) else set()

    pockets: list[dict] = []
    excluded: list[dict] = []

    stem = Path(pdb_path).stem
    pdb_id = re.sub(r"(?i)(_nolig(_cleaned)?|_cleaned|_fixed)$", "", stem).upper()

    for key, lines in ligands_dict.items():
        chain, resname, resnum = key
        resname_u = (resname or "").strip().upper()
        chain_u = (chain or "").strip() or "-"
        resnum_s = str(resnum).strip() if resnum is not None else ""
        resnum_u = resnum_s or "?"
        het_name_upper = hetnam_map.get(resname_u, "") or ""

        reason = ""
        if resname_u in exclude_ids:
            reason = "resname_in_exclude_het_ids"
        elif resname_u in solvent_res:
            reason = "remove_as_solvent_resname"
        elif resname_u in glycan_res:
            reason = "remove_as_glycan_resname"
        else:
            for kw in sorted(exclude_kw, key=lambda x: len(str(x)), reverse=True):
                kw_up = (kw or "").upper()
                if kw_up and kw_up in het_name_upper:
                    reason = f"het_name_keyword:{kw_up}"
                    break

        if reason:
            excluded.append(
                {
                    "ligand_resname": resname_u,
                    "ligand_chain": chain_u,
                    "ligand_resnum": resnum_u,
                    "reason": reason,
                }
            )
            continue

        coords = _coords_from_lines(lines)
        if not coords:
            logger.warning(
                "[pockets.json] skipping ligand with no coords resname=%s chain=%s resnum=%s",
                resname_u,
                chain_u,
                resnum_u,
            )
            continue

        xs, ys, zs = zip(*coords)
        atom_count = len(coords)
        center = [sum(xs) / atom_count, sum(ys) / atom_count, sum(zs) / atom_count]
        bounds = {
            "min": [min(xs), min(ys), min(zs)],
            "max": [max(xs), max(ys), max(zs)],
        }

        pockets.append(
            {
                "ligand_resname": resname_u,
                "ligand_chain": chain_u,
                "ligand_resnum": resnum_u,
                "protein_chain": chain_u,
                "het_name": het_name_upper,
                "atom_count": atom_count,
                "center": center,
                "bounds": bounds,
                "coords": [[x, y, z] for x, y, z in coords],
            }
        )

    pocket_chains = sorted(
        {p["protein_chain"] for p in pockets if p.get("protein_chain")}
    )
    return {
        "pdb_id": pdb_id,
        "pocket_chains": pocket_chains,
        "pockets": pockets,
        "excluded": excluded,
    }


def write_pockets_json(
    pdb_path: str, ligands_dict: Dict[tuple, List[str]], out_path: Path, logger=logging
) -> Path:
    manifest = build_ligand_pockets_manifest(pdb_path, ligands_dict, logger=logger)

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    logger.info(
        "[pockets.json] wrote=%s included=%d excluded=%d path=%s",
        Path(pdb_path).name,
        len(manifest.get("pockets", [])),
        len(manifest.get("excluded", [])),
        out,
    )
    return out
