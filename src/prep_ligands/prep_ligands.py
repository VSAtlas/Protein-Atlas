import logging
import os
import re
import shutil
import sys
import types
from pathlib import Path
from typing import Optional

# Ensure repo root is on sys.path so top-level helpers (activesite, etc.) resolve when run as a script.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
SRC_DIR = ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
if "prep_ligands" not in sys.modules:
    pkg = types.ModuleType("prep_ligands")
    pkg.__path__ = [str(SRC_DIR / "prep_ligands")]
    sys.modules["prep_ligands"] = pkg

from activesite import fix_element_columns_in_file  # noqa: E402,F401
try:
    from prep_ligands.prep_ligands_bulk import (  # noqa: E402,F401
        _audit_protonation_metrics,
        _collect_only_from_env_and_cli,
        _log_malformed,
        _pdbqt_from_mol2_via_obabel,
        _prepare_one,
        _valid_pdbqt,
        enumerate_ligands_for_docking,
        is_valid_ligand,
        prep_ligands_with_mgltools,
        read_config,
    )
except ModuleNotFoundError as exc:
    if exc.name not in {"prep_ligands", "prep_ligands.prep_ligands_bulk"}:
        raise
    from prep_ligands_bulk import (  # type: ignore  # noqa: E402,F401
        _audit_protonation_metrics,
        _collect_only_from_env_and_cli,
        _log_malformed,
        _pdbqt_from_mol2_via_obabel,
        _prepare_one,
        _valid_pdbqt,
        enumerate_ligands_for_docking,
        is_valid_ligand,
        prep_ligands_with_mgltools,
        read_config,
    )

try:
    from prep_ligands.prep_ligands_crystal import prep_ligands_from_pdb  # noqa: E402
except ModuleNotFoundError as exc:
    if exc.name not in {"prep_ligands", "prep_ligands.prep_ligands_crystal"}:
        raise
    from prep_ligands_crystal import prep_ligands_from_pdb  # type: ignore  # noqa: E402


logger = logging.getLogger(__name__)


def _as_path(p) -> Path:
    return p if isinstance(p, Path) else Path(p)


def _cfg_env_or_default(key: str, default: Optional[str] = None) -> Optional[str]:
    """Lightweight config reader that prefers env, then config.txt next to this file, else default."""
    v = os.environ.get(key)
    if v:
        return v
    try:
        # Move up from src/prep_ligands/ to repo root
        root = Path(__file__).resolve().parents[2]
        cfg = root / "config.txt"
        if cfg.is_file():
            for line in cfg.read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, val = line.split("=", 1)
                if k.strip() == key:
                    return val.strip()
    except Exception:
        pass
    return default


def _canon_base(output_root: Path, pdb_id: str) -> Path:
    """Canonical per-protein base dir: processed_pdbs/<PDB>"""
    output_root = _as_path(output_root)
    return (output_root / pdb_id.upper()).resolve()


def _merge_dir(src: Path, dst: Path) -> None:
    """Merge src directory into dst (mkdirs as needed); removes src if emptied."""
    src, dst = src.resolve(), dst.resolve()
    if not src.exists():
        return
    dst.mkdir(parents=True, exist_ok=True)
    for root, dirs, files in os.walk(src):
        r = Path(root)
        rel = r.relative_to(src)
        (dst / rel).mkdir(parents=True, exist_ok=True)
        for d in dirs:
            (dst / rel / d).mkdir(parents=True, exist_ok=True)
        for f in files:
            s = r / f
            t = dst / rel / f
            if t.exists():
                try:
                    if s.stat().st_size == t.stat().st_size:
                        continue
                except Exception:
                    pass
            shutil.move(str(s), str(t))
    try:
        shutil.rmtree(src)
    except Exception:
        pass


def fold_legacy_layout(pdb_id: str, output_root) -> None:
    """
    Override with extended handling to pull uppercase legacy dirs into the canonical tree:
      <PDB>_NOLIG            -> processed_pdbs/<PDB>/nolig
      <PDB>_CLEANED_LIGANDS  -> processed_pdbs/<PDB>/ligands_raw
      <PDB>_nolig(.pdb)      -> processed_pdbs/<PDB>/nolig/<PDB>_nolig_phenix_clean.pdb
    """
    try:
        root = _as_path(output_root).resolve()
        pdb_idU = pdb_id.upper()
        base = _canon_base(root, pdb_idU)
        (base / "nolig").mkdir(parents=True, exist_ok=True)
        (base / "ligands_raw").mkdir(parents=True, exist_ok=True)

        legacy_dirs = [
            (root / f"{pdb_id}_nolig", base / "nolig"),
            (root / f"{pdb_id.lower()}_nolig", base / "nolig"),
            (root / f"{pdb_idU}_NOLIG", base / "nolig"),
            (root / f"{pdb_id}_cleaned_ligands", base / "ligands_raw"),
            (root / f"{pdb_id.lower()}_cleaned_ligands", base / "ligands_raw"),
            (root / f"{pdb_idU}_CLEANED_LIGANDS", base / "ligands_raw"),
        ]
        for src, dst in legacy_dirs:
            if src.exists():
                logging.info("Migrating legacy directory %s -> %s", src, dst)
                _merge_dir(src, dst)

        candidates_files = [
            (root / f"{pdb_id}_nolig.pdb", base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
            (root / f"{pdb_id.lower()}_nolig.pdb", base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
        ]
        for src, dst in candidates_files:
            if src.exists() and not dst.exists():
                logging.info("Moving legacy file %s -> %s", src, dst)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))
    except Exception as e:
        logging.warning("fold_legacy_layout (extended) failed for %s: %s", pdb_id, e)


def _expose_ligand_intermediates_for_debug(pdb_id: str, ligands_raw: Path) -> None:
    """
    If PREPPED_LIGANDS_ROOT is configured, create/update:
      <PREPPED_LIGANDS_ROOT>/<PDB>/intermediates -> <processed_pdbs>/<PDB>/ligands_raw
    so intermediates (sanitized PDB, MOL2) are visible next to final PDBQTs.
    """
    try:
        prepped_root = _cfg_env_or_default("PREPPED_LIGANDS_ROOT", "")
        if not prepped_root:
            return
        dst_dir = _as_path(prepped_root) / pdb_id.upper()
        dst_dir.mkdir(parents=True, exist_ok=True)
        link = dst_dir / "intermediates"
        if link.is_symlink() or link.exists():
            try:
                if link.is_dir() and not link.is_symlink():
                    shutil.rmtree(link)
                else:
                    link.unlink()
            except Exception:
                pass
        link.symlink_to(ligands_raw.resolve(), target_is_directory=True)
        logging.info("Debug symlink: %s -> %s", link, ligands_raw)
    except Exception as e:
        logging.warning("Could not create debug symlink for %s: %s", pdb_id, e)


def _element_fix_all_in_dir(ligands_raw: Path) -> None:
    """
    Run your element-column fixer on every PDB in ligands_raw if the function exists.
    This reduces 'Unknown atom name F29 -> C' noise and improves template matches downstream.
    """
    try:
        fixer = globals().get("fix_element_columns_in_file", None)
        if fixer is None:
            return
        for p in ligands_raw.glob("*.pdb"):
            try:
                fixer(p, p)
            except Exception as e:
                logging.warning("Element-fix skipped for %s: %s", p.name, e)
    except Exception as e:
        logging.warning("Bulk element-fix failed in %s: %s", ligands_raw, e)


__all__ = [
    "prep_ligands_with_mgltools",
    "prep_ligands_from_pdb",
    "enumerate_ligands_for_docking",
]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Bulk ligand preparation (SDF→MOL2→PDBQT)")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Do not skip existing PDBQTs; overwrite outputs even if a valid PDBQT exists.",
    )
    parser.add_argument(
        "--migrate-legacy",
        metavar="PROCESSED_ROOT",
        help="Scan processed_pdbs and migrate *_NOLIG/*_CLEANED_LIGANDS into canonical layout.",
    )
    parser.add_argument(
        "--expose-intermediates",
        nargs=2,
        metavar=("PROCESSED_ROOT", "PREPPED_LIGANDS_ROOT"),
        help="Create/refresh prepped_ligands/<PDB>/intermediates symlinks for all PDBs.",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        help="Limit run to specific ligands (accepts rdk_0004931, rdk_4931, 0004931, 4931; comma/space OK)",
    )
    parser.add_argument(
        "--extracted",
        action="store_true",
        help="Run ligand prep for extracted-from-PDB ligands (ligands_raw/*.pdb).",
    )
    parser.add_argument(
        "--only-extracted",
        nargs="+",
        default=None,
        help="Limit extracted prep to specific ligand basenames (e.g., RXT_A1204). Mirrors EXTRACT_ONLY.",
    )
    parser.add_argument(
        "--test-extracted",
        action="store_true",
        help="Verbose test mode for extracted ligands (mirrors EXTRACT_TEST=1).",
    )
    parser.add_argument("--in-sdf", help="Path to a single SDF file (process only this file)")
    parser.add_argument("--in-sdf-dir", help="Directory of bulk SDFs (replaces default extracted-SDF root)")
    parser.add_argument("--in-pdb-dir", help="Directory of extracted ligand PDBs (crystal-safe path)")
    parser.add_argument("--mol2-dir", help="Directory for MOL2 intermediates")
    parser.add_argument("--out-pdbqt-dir", help="Destination directory for final PDBQTs")
    parser.add_argument("--status-log", help="Basename or full path for the TSV status log")
    parser.add_argument("--rename-prefix", default="", help="Force output ligand stem prefix (e.g., decoys_).")
    parser.add_argument("--rename-pad", type=int, default=5, help="Zero-pad width for rename index.")
    parser.add_argument("--rename-start", type=int, default=1, help="Starting index for renamed ligands.")
    parser.add_argument("--rename-force", action="store_true", help="Always rename ligands, ignoring input names.")

    args = parser.parse_args()

    def _set_env(k, v):
        if v is not None and str(v).strip() != "":
            os.environ[k] = str(v)

    _set_env("LIGPREP_IN_SDF", args.in_sdf)
    _set_env("LIGPREP_IN_SDF_DIR", args.in_sdf_dir)
    _set_env("LIGPREP_IN_PDB_DIR", args.in_pdb_dir)
    _set_env("LIGPREP_MOL2_DIR", args.mol2_dir)
    _set_env("LIGPREP_OUT_DIR", args.out_pdbqt_dir)
    _set_env("LIGPREP_STATUS_LOG", args.status_log)
    _set_env("LIGPREP_RENAME_PREFIX", args.rename_prefix)
    _set_env("LIGPREP_RENAME_PAD", args.rename_pad)
    _set_env("LIGPREP_RENAME_START", args.rename_start)
    if args.rename_force:
        _set_env("LIGPREP_RENAME_FORCE", "1")

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

    if getattr(args, "extracted", False):
        if getattr(args, "only_extracted", None):
            os.environ["EXTRACT_ONLY"] = " ".join(args.only_extracted)
        if getattr(args, "test_extracted", False):
            os.environ["EXTRACT_TEST"] = "1"

        cfg = read_config()
        ligand_extracted_dir = Path(cfg["LIGAND_EXTRACTED_DIR"]).resolve()
        ligands_mol2_dir = Path(cfg["LIGANDS_MOL2_DIR"]).resolve()
        output_ligands_dir = Path(cfg["OUTPUT_LIGANDS_DIR"]).resolve()

        prep_ligands_from_pdb(ligand_extracted_dir, ligands_mol2_dir, output_ligands_dir)
        sys.exit(0)

    if args.migrate_legacy:
        root = Path(args.migrate_legacy).resolve()
        for d in sorted(root.iterdir()):
            if not d.is_dir():
                continue
            m = re.match(r"^([A-Za-z0-9]{4})(?:_.+)?$", d.name)
            if not m:
                continue
            pdb_id = m.group(1).upper()
            try:
                fold_legacy_layout(pdb_id, root)
            except Exception as e:
                logging.warning("migrate-legacy skip %s: %s", pdb_id, e)
        sys.exit(0)

    if args.expose_intermediates:
        proc_root = Path(args.expose_intermediates[0]).resolve()
        prepped_root = Path(args.expose_intermediates[1]).resolve()
        os.environ["PREPPED_LIGANDS_ROOT"] = str(prepped_root)
        for pdb_dir in sorted(proc_root.iterdir()):
            if not pdb_dir.is_dir():
                continue
            pdb_id = pdb_dir.name.split("_")[0].upper()
            lig_raw = pdb_dir / "ligands_raw"
            if lig_raw.is_dir():
                _element_fix_all_in_dir(lig_raw)
                _expose_ligand_intermediates_for_debug(pdb_id, lig_raw)
        print("Exposed intermediates under:", prepped_root)
        sys.exit(0)

    prep_ligands_with_mgltools(
        force=args.force,
        only=(only_set if only_set else None),
        in_sdf=args.in_sdf,
        in_sdf_dir=args.in_sdf_dir,
        in_pdb_dir=args.in_pdb_dir,
        mol2_dir=args.mol2_dir,
        out_pdbqt_dir=args.out_pdbqt_dir,
        status_log=args.status_log,
    )

# Acceptance tests (manual, keep behavior unchanged):
# 1) Bulk subset:
#    export LIGPREP_LIBRARY=cah2
#    python prep_ligands.py --force --in-sdf-dir /home/michael/atlas/code/protein_automation/extracted_ligands/cah2 --only 1 2 3 4 5 6 7 8 9 10
# 2) Crystal direct-path:
#    python prep_ligands.py --extracted --only-extracted /home/michael/atlas/code/protein_automation/processed_pdbs/1BN1/ligands_raw/AL5_A555.pdb --test-extracted
# 3) Optional import/compile:
#    python -m py_compile prep_ligands.py prep_ligands_bulk.py prep_ligands_bulk_sdf.py prep_ligands_microstates.py prep_ligands_crystal.py
#    python - << "PY"
#    import prep_ligands
#    from prep_ligands import prep_ligands_with_mgltools, prep_ligands_from_pdb, enumerate_ligands_for_docking
#    print("imports_ok")
#    PY
