import logging
import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import ligand_pocket
import p2rank_pocket


def select_pocket(
    pdb_cleaned: str,
    pdb_file: str,
    pdb_id: str,
    ligands: Dict[Any, list[str]],
    ligands_dir,
    logger,
    override_box: Optional[object] = None,
    override_meta: Optional[Dict[str, Any]] = None,
) -> Tuple[object, object, str, dict]:
    """
    Returns (center_xyz, box_xyz, source_label, extra_metadata_dict)
    """
    source = "p2rank"
    if override_box is not None:
        if isinstance(override_box, dict):
            center = override_box.get("center")
            box_size = override_box.get("box_size")
        else:
            try:
                center, box_size = override_box  # type: ignore[misc]
            except Exception:
                center, box_size = None, None
        extra = dict(override_meta or {})
        return center, box_size, "pocket_eval_override", extra

    if ligands:
        logger.debug(
            "[activesite.main] ligand_keys=%s",
            list(ligands.keys()),
        )
    try:
        lig_files = sorted(Path(ligands_dir).glob("*.pdb"))
        logger.debug(
            "[activesite.main] ligands_dir_pdb_files=%s",
            [p.name for p in lig_files],
        )
    except Exception as e:
        logger.warning(
            "[activesite.main] unable to list ligands_dir=%s err=%s",
            ligands_dir,
            e,
        )

    if ligands:
        logging.info(f"Ligands removed for {pdb_file}. Ranking ligands by contacts...")
        for key, lines in ligands.items():
            logging.info(f"Ligand {key} has {len(lines)} atoms.")
        ranked_ligands = ligand_pocket.rank_ligands_by_atom_count(ligands)
        if ranked_ligands:
            # Select top ligand by contact count
            top_ligand_key, top_lines = ranked_ligands[0]
            logging.info(
                f"Top ligand by size: {top_ligand_key} with {len(top_lines)} atoms."
            )
            top_ligand_lines = top_lines

            ligand_output_dir = str(ligands_dir)
            os.makedirs(ligand_output_dir, exist_ok=True)
            ligand_path = os.path.join(ligand_output_dir, f"{pdb_id}.pdb")
            with open(ligand_path, "w") as f:
                f.writelines(ligands[top_ligand_key])
            logging.info(f"Top ligand saved to {ligand_path}")

            # Use coordinates of the top ligand only to compute box
            top_ligand_lines = ligands[top_ligand_key]
            top_ligand_coords = []
            for line in top_ligand_lines:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    top_ligand_coords.append((x, y, z))
                except ValueError:
                    logging.warning(
                        f"Invalid coordinates in ligand line: {line.strip()}"
                    )

            if not top_ligand_coords:
                logging.warning(
                    "Top ligand has no valid coordinates. Falling back to P2Rank."
                )
                center, box_size = p2rank_pocket.get_box_from_p2rank_csv(pdb_cleaned)
            else:
                center, box_size = ligand_pocket.compute_box_from_ligand_coords(
                    top_ligand_coords
                )
                source = "ligand_top"

        else:
            logging.warning("No ligands ranked, fallback to P2Rank.")
            center, box_size = p2rank_pocket.get_box_from_p2rank_csv(pdb_cleaned)
    else:
        logger.info(
            "[activesite.main] no_ligands_after_extract pdb_cleaned=%s ligands_dir=%s",
            pdb_cleaned,
            ligands_dir,
        )
        center, box_size = p2rank_pocket.get_box_from_p2rank_csv(pdb_cleaned)

    return center, box_size, source, {}
