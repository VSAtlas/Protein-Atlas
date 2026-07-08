import argparse
import logging
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from config.runtime_config import load_config
from protein_prep.pdb_fixer_runtime import fix_element_columns_in_file  # noqa: F401
from prep_ligands.prep_ligands_bulk_meeko import prep_ligands_with_meeko
from prep_ligands.prep_ligands_common_validation import read_config
from prep_ligands.prep_ligands_microstates import (
    _collect_only_from_env_and_cli,
    enumerate_ligands_for_docking,
)
from prep_ligands.prep_ligands_crystal import prep_ligands_from_pdb
from protein_prep.layout import fold_legacy_layout

logger = logging.getLogger(__name__)

prep_ligands_with_mgltools = prep_ligands_with_meeko


def _cfg_env_or_default(key: str, default: Optional[str] = None) -> Optional[str]:
    value = os.environ.get(key)
    if value:
        return value
    try:
        cfg: Dict[str, Any] = load_config()
    except Exception:
        cfg = {}
    resolved = cfg.get(key)
    if resolved not in (None, ""):
        return str(resolved)
    return default


def _refresh_intermediates_link(pdb_id: str, ligands_raw: Path) -> None:
    try:
        prepped_root = _cfg_env_or_default("PREPPED_LIGANDS_DIR", "")
        if not prepped_root:
            return
        link = Path(prepped_root) / pdb_id.upper() / "intermediates"
        link.parent.mkdir(parents=True, exist_ok=True)
        if link.exists() or link.is_symlink():
            if link.is_dir() and not link.is_symlink():
                shutil.rmtree(link)
            else:
                link.unlink()
        link.symlink_to(ligands_raw.resolve(), target_is_directory=True)
        logging.info("Debug symlink: %s -> %s", link, ligands_raw)
    except Exception as exc:
        logging.warning("Could not create debug symlink for %s: %s", pdb_id, exc)


def _fix_elements_in_dir(ligands_raw: Path) -> None:
    try:
        for pdb_file in ligands_raw.glob("*.pdb"):
            try:
                fix_element_columns_in_file(pdb_file, pdb_file)
            except Exception as exc:
                logging.warning("Element-fix skipped for %s: %s", pdb_file.name, exc)
    except Exception as exc:
        logging.warning("Bulk element-fix failed in %s: %s", ligands_raw, exc)


__all__ = [
    "prep_ligands_with_meeko",
    "prep_ligands_with_mgltools",
    "prep_ligands_from_pdb",
    "enumerate_ligands_for_docking",
]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Bulk ligand preparation (SDF/restored-SDF->PDBQT via Meeko)"
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--migrate-legacy", metavar="PROCESSED_ROOT")
    parser.add_argument(
        "--expose-intermediates",
        nargs=2,
        metavar=("PROCESSED_ROOT", "PREPPED_LIGANDS_DIR"),
    )
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--extracted", action="store_true")
    parser.add_argument("--only-extracted", nargs="+", default=None)
    parser.add_argument("--test-extracted", action="store_true")
    parser.add_argument("--in-sdf")
    parser.add_argument("--in-sdf-dir")
    parser.add_argument("--in-pdb-dir")
    parser.add_argument("--mol2-dir")
    parser.add_argument("--out-pdbqt-dir")
    parser.add_argument("--status-log")
    parser.add_argument("--rename-prefix", default="")
    parser.add_argument("--rename-pad", type=int, default=5)
    parser.add_argument("--rename-start", type=int, default=1)
    parser.add_argument("--rename-force", action="store_true")
    return parser


def _set_cli_env(args: argparse.Namespace) -> None:
    for env_name, value in {
        "LIGPREP_IN_SDF": args.in_sdf,
        "LIGPREP_IN_SDF_DIR": args.in_sdf_dir,
        "LIGPREP_IN_PDB_DIR": args.in_pdb_dir,
        "LIGPREP_MOL2_DIR": args.mol2_dir,
        "LIGPREP_OUT_DIR": args.out_pdbqt_dir,
        "LIGPREP_STATUS_LOG": args.status_log,
        "LIGPREP_RENAME_PREFIX": args.rename_prefix,
        "LIGPREP_RENAME_PAD": args.rename_pad,
        "LIGPREP_RENAME_START": args.rename_start,
    }.items():
        if value is not None and str(value).strip() != "":
            os.environ[env_name] = str(value)
    if args.rename_force:
        os.environ["LIGPREP_RENAME_FORCE"] = "1"


def _run_extracted_mode(args: argparse.Namespace) -> int:
    if args.only_extracted:
        os.environ["EXTRACT_ONLY"] = " ".join(args.only_extracted)
    if args.test_extracted:
        os.environ["EXTRACT_TEST"] = "1"
    cfg = read_config()
    prep_ligands_from_pdb(
        Path(cfg["EXTRACTED_LIGANDS_DIR"]).resolve(),
        Path(cfg["LIGANDS_MOL2_DIR"]).resolve(),
        Path(cfg["PREPPED_LIGANDS_DIR"]).resolve(),
    )
    return 0


def _run_legacy_migration(root: Path) -> int:
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        match = re.match(r"^([A-Za-z0-9]{4})(?:_.+)?$", child.name)
        if not match:
            continue
        pdb_id = match.group(1).upper()
        try:
            fold_legacy_layout(pdb_id, root)
        except Exception as exc:
            logging.warning("migrate-legacy skip %s: %s", pdb_id, exc)
    return 0


def _run_expose_intermediates(proc_root: Path, prepped_root: Path) -> int:
    os.environ["PREPPED_LIGANDS_DIR"] = str(prepped_root)
    for pdb_dir in sorted(proc_root.iterdir()):
        if not pdb_dir.is_dir():
            continue
        lig_raw = pdb_dir / "ligands_raw"
        if not lig_raw.is_dir():
            continue
        pdb_id = pdb_dir.name.split("_")[0].upper()
        _fix_elements_in_dir(lig_raw)
        _refresh_intermediates_link(pdb_id, lig_raw)
    print("Exposed intermediates under:", prepped_root)
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    _set_cli_env(args)
    print(
        "[paths.effective]",
        f"in_sdf={os.environ.get('LIGPREP_IN_SDF', 'None')}",
        f"in_sdf_dir={os.environ.get('LIGPREP_IN_SDF_DIR', 'None')}",
        f"in_pdb_dir={os.environ.get('LIGPREP_IN_PDB_DIR', 'None')}",
        f"mol2_dir={os.environ.get('LIGPREP_MOL2_DIR', 'None')}",
        f"out_pdbqt_dir={os.environ.get('LIGPREP_OUT_DIR', 'None')}",
        f"status_log={os.environ.get('LIGPREP_STATUS_LOG', 'None')}",
    )
    only_set = _collect_only_from_env_and_cli(args.only if hasattr(args, "only") else None)
    if args.extracted:
        return _run_extracted_mode(args)
    if args.migrate_legacy:
        return _run_legacy_migration(Path(args.migrate_legacy).resolve())
    if args.expose_intermediates:
        proc_root, prepped_root = (Path(part).resolve() for part in args.expose_intermediates)
        return _run_expose_intermediates(proc_root, prepped_root)
    prep_ligands_with_meeko(
        force=args.force,
        only=only_set if only_set else None,
        in_sdf=args.in_sdf,
        in_sdf_dir=args.in_sdf_dir,
        in_pdb_dir=args.in_pdb_dir,
        mol2_dir=args.mol2_dir,
        out_pdbqt_dir=args.out_pdbqt_dir,
        status_log=args.status_log,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
