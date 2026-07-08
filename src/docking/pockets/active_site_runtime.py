import os
import shutil
import logging
import re
from config.runtime_config import load_config
from config.normalize import _to_bool
from pathlib import Path

from cli.cli_utils import _norm_pdb_id
from path_router import make_paths, expand_variants
from cli.run_manifest_runtime import PocketDetectionEvent, emit_pocket_detection_event
from docking.docking_utils import norm
from protein_prep import pdb_fixer_runtime as pdb_fixer
from docking.pockets import ligand_pocket_runtime as ligand_pocket
from docking.pockets import pocket_policy_runtime as pocket_policy


from protein_prep.pdb_fixer_runtime import (  # noqa: F401  # Back-compat re-exports; prefer importing from the new module.
    assert_no_helium_in_hydrogen_names,
    assert_no_helium_in_pdbqt,
    derive_element,
    ensure_model_records,
    fix_element_columns_in_file,
    fix_ligand_element_columns,
    fix_pdb_elements,
    get_atom_rules,
    load_canonical_cofactors,
    load_canonical_metals,
    load_canonical_waters,
    log_pre_variant_policy_breadcrumb,
    remove_unparsable_hetatms,
    rules_version,
    scan_helium_counts,
    scan_helium_counts_with_hits,
    summarize_ions,
)  # noqa: F401
from docking.pockets.ligand_pocket_runtime import (  # noqa: F401  # Back-compat re-exports; prefer importing from the new module.
    calculate_ligand_protein_contacts,
    compute_box_from_ligand_coords,
    extract_and_remove_ligands,
    rank_ligands_by_atom_count,
)  # noqa: F401
from docking.pockets.p2rank_runtime import detect_pocket, get_box_from_p2rank_csv  # noqa: F401
from docking.pockets.pocket_policy_runtime import select_pocket  # noqa: F401

logger = logging.getLogger(__name__)


def _load_runtime_config() -> dict:
    try:
        cfg = load_config() or {}
        return cfg if isinstance(cfg, dict) else {}
    except Exception:
        return {}


def _default_variant(cfg):
    """Resolve a single variant preference from config/environment."""
    mode = os.environ.get("APO_HOLO_MODE") or (
        cfg.get("APO_HOLO_MODE") if isinstance(cfg, dict) else None
    )
    for v in expand_variants(mode):
        return v
    return None


def prepare_receptor(cfg, paths, logger):
    """
    Run or reuse protein preparation to produce:
      - cleaned PDB without ligands
      - receptor PDBQT
    Returns (cleaned_pdb_path_str, receptor_pdbqt_path_str) or (None, None) on failure.
    """
    variant = _default_variant(cfg)
    ph_token = None

    receptor_cleaned = paths.receptor_cleaned_pdb(variant)
    receptor_pdbqt_path = paths.receptor_pdbqt(variant, ph_token)

    ligands_raw_dir = paths.ligand_output_dir
    p2rank_work_dir = paths.work_dir / "p2rank"

    ligands_raw_dir.mkdir(parents=True, exist_ok=True)
    p2rank_work_dir.mkdir(parents=True, exist_ok=True)

    force_reprocess = bool(_to_bool(cfg.get("FORCE_REPROCESS", False), False))
    logger.info(
        f"FORCE_REPROCESS={force_reprocess} | "
        f"exists(cleaned)={receptor_cleaned.exists()} "
        f"exists(receptor)={receptor_pdbqt_path.exists()}"
    )

    if (
        receptor_cleaned.exists()
        and receptor_pdbqt_path.exists()
        and not force_reprocess
    ):
        logger.info("Reusing existing cleaned PDB and receptor PDBQT.")
        return norm(receptor_cleaned), norm(receptor_pdbqt_path)
    from protein_prep import automate_protein_prep

    result = automate_protein_prep.main(str(paths.nolig_pdb_path))
    if not result or not isinstance(result, tuple) or len(result) != 2:
        logger.warning("Protein prep failed.")
        return None, None

    cleaned_pdb, receptor_pdbqt = result

    # Ensure receptor lives in canonical PDBQT_DIR
    try:
        if Path(receptor_pdbqt).resolve() != receptor_pdbqt_path.resolve():
            from shutil import copy2

            receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
            copy2(receptor_pdbqt, receptor_pdbqt_path)
            receptor_pdbqt = str(receptor_pdbqt_path)
    except Exception as e:
        logger.warning(f"Could not relocate receptor PDBQT: {e}")

    return norm(cleaned_pdb), norm(receptor_pdbqt)


# ---------- high-level pipeline steps ----------


def extract_ligands(cfg, paths, logger):
    """
    Extract and strip ligands from input PDB into clean PDB without ligands.
    Returns the ligand count.
    """
    malformed_log = paths.ligands_mol2_dir / "malformed_ligands.txt"
    if malformed_log.exists():
        malformed_log.unlink()

    ligands_dict, _ = extract_and_remove_ligands(
        str(paths.input_pdb_path),
        str(paths.nolig_pdb_path),
        str(paths.ligand_output_dir),
    )
    logger.info(f"Extracted {len(ligands_dict)} ligands → {paths.ligand_output_dir}")
    return len(ligands_dict)


# ---------- small utils ----------


def get_recenter_params(cfg):
    """Read early/fallback recentering knobs from config with safe defaults."""
    return {
        "EARLY_RECENTER_RATIO": float(cfg.get("EARLY_RECENTER_RATIO", 0.70)),
        "EARLY_RECENTER_MIN_EVAL": int(cfg.get("EARLY_RECENTER_MIN_EVAL", 10)),
        "EARLY_RECENTER_FAR_A": float(cfg.get("EARLY_RECENTER_FAR_A", 15.0)),
        "EARLY_RECENTER_MEDIAN_A": float(cfg.get("EARLY_RECENTER_MEDIAN_A", 10.0)),
        "ALLOW_BOX_EXPAND": bool(cfg.get("ALLOW_BOX_EXPAND", True)),
        "MAX_RECENTER_ATTEMPTS": int(cfg.get("MAX_RECENTER_ATTEMPTS", 3)),
    }


def load_aliases():
    # back-compat shim
    return pdb_fixer.get_atom_rules()


def _sanitize_pdb_id(stem: str) -> str:
    """
    Strip legacy trailing tags from a filename stem, repeatedly, case-insensitively.
    Examples:
      4GRL_CLEANED_cleaned -> 4GRL
      2SRC_nolig_cleaned   -> 2SRC
      3ERT_withH_fixed     -> 3ERT
    """
    import re

    s = stem
    # remove multiple trailing tags if present
    while True:
        s2 = re.sub(r"(?i)(?:_(?:cleaned|nolig|withh|fixed))$", "", s)
        if s2 == s:
            return s
        s = s2


def _canon_pdb_id_from_path(pdb_path: str) -> str:
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    # Drop legacy suffixes like _nolig, _nolig_cleaned, _cleaned
    cleaned = re.sub(r"(?i)(_nolig(_cleaned)?|_cleaned)$", "", stem)
    normalized = _norm_pdb_id(cleaned)
    if normalized:
        return normalized
    return cleaned.upper()


def main(pdb_file):
    exists = os.path.exists(pdb_file)
    size = os.path.getsize(pdb_file) if exists else -1
    stem = os.path.splitext(os.path.basename(pdb_file))[0]
    logger.info(
        "[activesite.main] entry pdb_file=%s exists=%s size=%s stem=%s",
        pdb_file,
        exists,
        size,
        stem,
    )
    pdb_id = _canon_pdb_id_from_path(pdb_file)

    logger.info(
        "[activesite.main] canon_pdb_id=%s from=%s",
        pdb_id,
        pdb_file,
    )

    config = _load_runtime_config()
    paths = make_paths(config, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant = _default_variant(config)
    ph_token = None

    receptor_cleaned = paths.receptor_cleaned_pdb(variant)
    receptor_pdbqt_path = paths.receptor_pdbqt(variant, ph_token)

    ligands_raw_dir = paths.ligand_output_dir
    p2rank_work_dir = paths.work_dir / "p2rank"
    pockets_json_out = p2rank_work_dir / "pockets.json"

    logger.info(
        "[activesite.main] paths base_id=%s receptor_cleaned=%s ligands_raw_dir=%s work_dir=%s",
        pdb_id,
        receptor_cleaned,
        ligands_raw_dir,
        paths.work_dir,
    )

    receptor_cleaned.parent.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    ligands_raw_dir.mkdir(parents=True, exist_ok=True)
    p2rank_work_dir.mkdir(parents=True, exist_ok=True)

    pdb_cleaned = str(receptor_cleaned)
    temp_fixed_pdb_path = paths.work_dir / f"{pdb_id}_fixed.pdb"
    ligands_dir = ligands_raw_dir
    try:
        temp_fixed_pdb_path.parent.mkdir(parents=True, exist_ok=True)

        # Prefer the canonical input PDB (HOLO) for ligand-based box detection.
        # Fallback to the passed-in pdb_file if input_pdb_path is missing.
        src_from_paths = getattr(paths, "input_pdb_path", None)
        if src_from_paths and Path(src_from_paths).exists():
            src_pdb_for_box = Path(src_from_paths)
        else:
            src_pdb_for_box = Path(pdb_file)

        logging.info(
            "[activesite.main] using src_pdb_for_box=%s (exists=%s) to build %s",
            src_pdb_for_box,
            src_pdb_for_box.exists(),
            temp_fixed_pdb_path,
        )

        shutil.copyfile(str(src_pdb_for_box), temp_fixed_pdb_path)
        logging.info(f"Copied PDB for fixing: {temp_fixed_pdb_path}")
        pdb_fixer.fix_pdb_elements(str(temp_fixed_pdb_path))

        ligands, ligand_coords = ligand_pocket.extract_and_remove_ligands(
            str(temp_fixed_pdb_path),
            pdb_cleaned,
            str(ligands_dir),
        )

        try:
            processed_root = getattr(paths, "root_pdb_dir", None)
            if processed_root is None:
                processed_root = Path(paths.work_dir).parent
            pockets_dir = Path(processed_root) / "pockets"
            pockets_json_out = pockets_dir / "pockets.json"
            ligand_pocket.write_pockets_json(
                str(src_pdb_for_box), ligands, pockets_json_out, logger
            )
        except Exception as exc:  # pragma: no cover - avoid blocking pocket selection
            logger.warning("[activesite.main] pockets.json write failed: %s", exc)

        n_lig = len(ligands) if ligands else 0
        logger.info(
            "[activesite.main] extract_and_remove_ligands in=%s out=%s ligands_dir=%s n_ligands=%s",
            temp_fixed_pdb_path,
            pdb_cleaned,
            ligands_dir,
            n_lig,
        )

        center, box_size, source, _ = pocket_policy.select_pocket(
            pdb_cleaned,
            pdb_file,
            pdb_id,
            ligands,
            ligands_dir,
            logger,
        )

        if center and box_size:
            logging.info(f"{pdb_cleaned}: center={center}, box_size={box_size}")
            logger.info(
                "[activesite.main] success center=%s box=%s pdb_cleaned=%s",
                center,
                box_size,
                pdb_cleaned,
            )
            logging.info(
                "[active-site] source=%s receptor=%s center=%s size=%s",
                source,  # e.g. "ligand_top" or "p2rank"
                receptor_pdbqt_path,
                center,
                box_size,
            )
            try:
                run_id = config.get("RUN_ID")
                logging.debug(
                    "[run-manifest.pocket_detection.call] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%r box=%r",
                    run_id,
                    pdb_id,
                    variant,
                    None,
                    source,
                    center,
                    box_size,
                )
                if run_id:
                    emit_pocket_detection_event(
                        config,
                        PocketDetectionEvent(
                            run_id=str(run_id),
                            pdb_id=pdb_id,
                            variant_label=variant,
                            ph_tag=None,
                            method=source,
                            center=center,
                            box_size=box_size,
                        ),
                    )
                else:
                    logging.debug(
                        "[run-manifest.pocket_detection.skip] no RUN_ID for pdb=%s variant=%s ph=%s",
                        pdb_id,
                        variant,
                        None,
                    )
            except Exception as exc:
                logging.warning(
                    "[run-manifest.skip] active-site pdb=%s err=%s",
                    pdb_id,
                    exc,
                )
            return center, box_size, source
        else:
            logging.warning(f"Box not determined for {pdb_cleaned}")
            logger.warning(
                "[activesite.main] failure center=None box=None pdb_cleaned=%s",
                pdb_cleaned,
            )
            return None, None, source

    finally:
        try:
            temp_fixed_pdb_path.unlink()
            logging.info(f"Temporary file removed: {temp_fixed_pdb_path}")
        except OSError:
            logging.warning(f"Could not delete temp file: {temp_fixed_pdb_path}")
